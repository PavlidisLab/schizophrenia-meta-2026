"""Step 2 of the Merikangas tracker: turn retrieved files into one DE gene
list per study contrast, with genes normalized to current HGNC symbols.

    python3 extract.py

Reads manifest.json (one entry per list: which file, which rows, which
columns, what threshold) and writes

    lists/<pmid>__<list_id>.tsv   the list, with #-prefixed provenance lines
    summary.tsv                   one row per list, incl. the count the paper
                                  states and whether the extraction matches it

Gene normalization happens in one place (normalize_gene) and runs in order:
a manifest-supplied override (symbol_overrides), probe -> platform annotation
(Gemma's annotation of the GPL where Gemma has it, else the GEO .annot file;
gemma_first), Ensembl ID -> HGNC, symbol -> HGNC approved / previous /
alias symbol. A printed name that is only a bare clone-library prefix (RP11,
CTD, ...) is read as no name (fix X08).

When those steps leave a row without a symbol (and no override applies), the
fallbacks of fixes X09-X11 (analyses/list_extraction_audit/
proposed_fixes.json) are tried, first hit wins; they reproduce the recovery
routes of analyses/list_extraction_audit/a1_unmapped.py:
  X11  a probe ID (reported ID or printed name) on the list's platform ->
       Gemma's annotation of the same GPL (extraction_audit_2026-09-28/
       platform_annots/<GPL>.gemma.annot.gz), when its gene symbols all name
       one HGNC gene, or else its one NCBI Gene ID carries an HGNC ID. Used
       only when the manifest names a GEO platform keyed by probe ID.
  X10  identifiers: an Ensembl gene ID HGNC lacks -> NCBI gene_info dbXrefs
       (one HGNC gene); an NCBI Gene ID (lists whose spec is a Gemma
       reanalysis or names an Entrez/NCBI/gene-id column) -> HGNC entrez_id,
       else gene_info's HGNC cross-reference; a RefSeq/GenBank accession ->
       HGNC refseq_accession / ena, else the 'GenBank Accession' column of
       every local GEO .annot (unanimous); a UniGene cluster -> the 'UniGene
       ID' column of every local GEO .annot (unanimous).
  X09  names: LOC<GeneID> -> gene_info's HGNC cross-reference; a spreadsheet
       date (1-Sep, Mar-02, 2019-03-05 00:00:00, serial 43530) -> SEPT/
       SEPTIN/MARCH/MARCHF/DEC<day> candidates (and SEP<day>), taken when
       they resolve to exactly one gene; otherwise the printed name -> gene_names.Resolver
       (analyses/gene_names/), tiers exact_approved, unique_alias and
       curated only (an ambiguous name stays unmapped).
An X10 or X11 result is dropped when the printed name itself means other
genes (the resolver's symbol or candidates, or a LOC number's gene), as the
audit's 'conflict' status does. Files a list's fallbacks read are named, with
sha256, in its header.

Fixes 1-5 of the 2026-09-30 Step 2 symbol audit (analyses/step2_symbol_audit_2026-09-30/):
  1  a full date or serial on the 1st is also read with its year as the number (SEPT8 -> 2008-09-01), and MARC<n>
     (MTARC1/2) is tried for March; a date that can mean several genes stays unmapped ('excel-date-ambiguous:').
  2  UniProt accessions in the reported-ID column -> HGNC (X10, uniprot_route).
  3  a probe ID printed in another case matches the platform's probe (platform_probes).
  4  a string that is one gene's previous symbol and another's alias is ambiguous unless HGNC renamed the
     previous-symbol gene after the paper's year (symbol_to_hgnc); the row's own ID may settle it, a name route may not.
  5  where the platform names another gene than the printed name, the route records the paper's gene
     ('platform->approved;paper:MT1M'). The platform reading wins (Paul, 2026-09-30); a probe the platform annotates
     to several genes, none the paper's, stays unmapped instead of falling back to the printed name.
Gemma-first probe mapping (Paul, 2026-09-30; sized by analyses/not_expressed_audit_2026-09-30/gemma_vs_geo_first.tsv):
for a reported probe ID on a manifest platform Gemma annotates (extraction_audit_2026-09-28/platform_annots/, or a
production export data/platforms/<GPL>.gemma-prod.annot.gz where gemma2's is absent), Gemma's alignment-based annotation
decides. Gemma's 'no gene' is unmapped with no fallback to the GEO .annot or the printed name; a probe Gemma annotates to
several genes stays unmapped unless they are one HGNC gene (the printed name does not pick among them); a gene HGNC lacks
is unmapped. map_route records the GEO .annot's and the paper's readings where they differ (';geo:X', ';geo:none',
';geo-multi:A,B', ';paper:X'). The GEO .annot reads only platforms Gemma lacks and probes absent from Gemma's table
(';gemma:absent'). Lists with no manifest platform are unchanged (their probe IDs are not looked up).
Rows are tagged in map_route (';tag:control-probe', ';tag:reference-gene'; row_tags), not dropped.
';tag:printed-name-fallback' (2026-09-30, forking paths): the manifest's platform holds the row's probe ID but its
annotation gives the probe no HGNC gene (no symbol, or a symbol HGNC lacks), and the symbol came from the paper's printed
name (printed_name_fallback). recurrence.py can drop these rows (PRINTED_NAME_FALLBACK=drop); the default keeps them.

    python3 extract.py --lists <name>,<name>    re-extract only these lists (summary.tsv rows of the others are kept)
"""
import csv, datetime, glob, gzip, hashlib, json, operator, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")  # repo data/: lists, manifests, reference, raw, universe sets
OUT = os.path.join(DATA, "lists")
HGNC_FILE = os.path.join(DATA, "reference", "hgnc_complete_set.txt")
OPS = {"<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq}


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


# ---------------------------------------------------------------- HGNC
# approved symbol -> year of HGNC's date_symbol_changed (filled by load_hgnc; symbol_to_hgnc's date rule)
SYMBOL_CHANGED_YEAR = {}


def load_hgnc():
    approved, previous, alias, ensembl = {}, {}, {}, {}
    with open(HGNC_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "Approved":
                continue
            sym = r["symbol"]
            approved[sym.upper()] = (sym, r["hgnc_id"])
            for p in filter(None, r["prev_symbol"].strip('"').split("|")):
                previous.setdefault(p.upper(), set()).add((sym, r["hgnc_id"]))
            for a in filter(None, r["alias_symbol"].strip('"').split("|")):
                alias.setdefault(a.upper(), set()).add((sym, r["hgnc_id"]))
            if r["ensembl_gene_id"]:
                ensembl[r["ensembl_gene_id"]] = (sym, r["hgnc_id"])
            if r["date_symbol_changed"][:4].isdigit():
                SYMBOL_CHANGED_YEAR[sym] = int(r["date_symbol_changed"][:4])
    return approved, previous, alias, ensembl


HGNC = None


def symbol_to_hgnc(s, year=None, official=False):
    """A symbol -> (symbol, hgnc_id, route): approved, else a previous symbol of one gene, else an alias of one gene.
    A string that is one gene's previous symbol and another gene's alias (SEPT2: previous of SEPTIN2, alias of SEPTIN6;
    74 list rows in the 2026-09-30 Step 2 audit) is taken as the previous-symbol gene only when HGNC changed that gene's
    symbol after the paper's year (then the string was its approved symbol when the paper was written; route
    'previous-dated'); otherwise it is 'previous-alias-ambiguous'. HGNC's date_symbol_changed is the gene's latest symbol
    change, which for a gene renamed twice may postdate the retirement of this string. official=True is for symbols
    from a source that writes official symbols only (a platform annotation, NCBI Gene, an HGNC cross-reference): an
    alias is never such a symbol, so a previous symbol there is the previous-symbol gene."""
    approved, previous, alias, _ = HGNC
    k = s.strip().upper()
    if k in approved:
        return approved[k] + ("approved",)
    prev, als = previous.get(k, set()), alias.get(k, set())
    if len(prev) == 1 and als - prev and not official:
        hit = next(iter(prev))
        changed = SYMBOL_CHANGED_YEAR.get(hit[0])
        if year and changed and changed > year:
            return hit + ("previous-dated",)
        return ("", "", "previous-alias-ambiguous:" + ",".join(sorted({h[0] for h in prev | als})))
    for hits, route in ((prev, "previous"), (als, "alias")):
        if len(hits) == 1:
            return next(iter(hits)) + (route,)
        if len(hits) > 1:
            return ("", "", f"{route}-ambiguous:" + ",".join(sorted(h[0] for h in hits)))
    return ("", "", "unmapped")


def paper_year(spec):
    """The paper's year from the list label ('Wu 2012'), as the 2026-09-30 audit read it (a4b_dup_junk.year_of)."""
    m = re.search(r"(19|20)\d{2}", spec.get("label") or "")
    return int(m.group(0)) if m else None


# ---------------------------------------------------------------- platforms
PLATFORMS = {}


def platform(path, id_column=None):
    """probe ID -> list of symbols, from a GEO GPL .annot(.gz) file. With
    id_column (e.g. 'GenBank Accession'), keys are that column instead of the
    probe ID, for papers that report the clone's accession rather than the spot."""
    key = (path, id_column)
    if key not in PLATFORMS:
        opener = gzip.open if path.endswith(".gz") else open
        table, header = {}, None
        with opener(os.path.join(DATA, path), "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("#", "!", "^")):
                    continue
                row = line.rstrip("\n").split("\t")
                if header is None:
                    header = row
                    sym_col = header.index("Gene symbol")
                    id_col = header.index(id_column) if id_column else 0
                    continue
                if len(row) > max(sym_col, id_col) and row[id_col]:
                    syms = [s for s in row[sym_col].split("///") if s]
                    table.setdefault(row[id_col], [])
                    table[row[id_col]] += [s for s in syms if s not in table[row[id_col]]]
        PLATFORMS[key] = table
    return PLATFORMS[key]


PLATFORM_INDEX = {}


def platform_probes(path, id_column, rid):
    """Probe IDs matching a reported ID: the probe itself, probes whose part before the
    first '_' is the ID, or probes containing '-ID_'. Indexed, because a linear scan per
    row takes hours on a full-statistics table against a 54,000-probe platform."""
    return probe_matches((path, id_column), platform(path, id_column), rid)


def probe_matches(key, table, rid):
    """platform_probes on any probe-keyed table (a GEO .annot, Gemma's annotation of the same GPL); key names the
    table in the index caches."""
    if key not in PLATFORM_INDEX:
        prefix, dashed = {}, {}
        for p in table:
            prefix.setdefault(p.split("_")[0], []).append(p)
            for m in re.finditer(r"-([^-]*?)_", p):
                dashed.setdefault(m.group(1), []).append(p)
        PLATFORM_INDEX[key] = (prefix, dashed)
    prefix, dashed = PLATFORM_INDEX[key]
    hits = set(prefix.get(rid, [])) | set(dashed.get(rid, []))
    if rid in table:
        hits.add(rid)
    if not hits and ("-" in rid or "_" in rid):
        # IDs containing the separators themselves: fall back to the plain scan
        hits = {p for p in table if p.split("_")[0] == rid or f"-{rid}_" in p}
    if not hits:
        # the whole ID in another case: 31123247 prints '202170_S_AT' for GPL570's 202170_s_at (3,169 rows; 2026-09-30
        # Step 2 audit). Only a unique case-folded match is taken.
        if key not in PLATFORM_CASEFOLD:
            folded = {}
            for p in table:
                folded.setdefault(p.upper(), []).append(p)
            PLATFORM_CASEFOLD[key] = folded
        cand = PLATFORM_CASEFOLD[key].get(rid.upper(), [])
        if len(cand) == 1:
            hits = set(cand)
    return hits


GEO_ID_COLUMN = {}


def geo_probes_by_column(path, id_column, value):
    """The GEO .annot's probe IDs (first column) whose id_column cell (e.g. 'GenBank Accession') is value, version
    stripped, upper case: for spotted / early arrays whose papers print the clone's accession, so that Gemma's
    annotation of the same probes can be read."""
    key = (path, id_column)
    if key not in GEO_ID_COLUMN:
        opener = gzip.open if path.endswith(".gz") else open
        idx, header = {}, None
        with opener(os.path.join(DATA, path), "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("#", "!", "^")):
                    continue
                row = line.rstrip("\n").split("\t")
                if header is None:
                    header, ci = row, row.index(id_column)
                    continue
                if len(row) > ci and row[0] and row[ci]:
                    for v in row[ci].split("///"):
                        idx.setdefault(v.strip().upper().split(".")[0], set()).add(row[0])
        GEO_ID_COLUMN[key] = idx
    return GEO_ID_COLUMN[key].get(value.strip().upper().split(".")[0], set())


PLATFORM_CASEFOLD = {}


# cells that mean "no gene" and must never be looked up: HGNC keeps NA as a previous
# symbol of XK, so an "NA" cell otherwise maps to XK (77 rows in 7 lists, 2026-09-24 review)
MISSING = {"", "NA", "N/A", "---"}
# spreadsheet error values: no information, and a leading '#' makes every reader skip the row as a comment
EXCEL_ERRORS = {"#N/A", "#VALUE!", "#REF!", "#NAME?", "#DIV/0!", "#NUM!", "#NULL!"}
MISSING |= EXCEL_ERRORS
# footnote marks papers attach to gene names (e.g. Pietersen 2014's GWAS asterisk)
FOOTNOTE_MARKS = "*\u2020\u2021"
GREEK = {"\u03b1": "ALPHA", "\u03b2": "BETA", "\u03b3": "GAMMA", "\u03b4": "DELTA", "\u03b5": "EPSILON",
         "\u03b6": "ZETA", "\u03b7": "ETA", "\u03b8": "THETA", "\u03ba": "KAPPA", "\u03bb": "LAMBDA"}
# hsa-miR-380-3p, hsa-miR-199a*, mir-643, hsa-let-7a-5p, hsamiR137: species prefix, mir/let,
# precursor name, optional arm (-3p/-5p), optional star
MIRNA = re.compile(r"^(?:hsa-?)?(mir|let)-?(.+?)(?:-([35]p))?(\*)?$", re.I)
MIR_FAMILY = None


def is_missing(v):
    return (v or "").strip().upper() in MISSING


def mirna_to_hgnc(name):
    """A miRNA name -> its precursor gene (MIR...). The arm (-3p/-5p) and star are
    dropped and the precursor is looked up as HGNC lists it: first the miRBase name,
    which HGNC keeps as an alias (hsa-mir-433 -> MIR433), then the approved symbol
    built from it (hsa-miR-451a -> MIR451A; hsa-miR-181a-2 -> MIR181A2). A name shared
    by several precursors (hsa-miR-199a -> MIR199A1, MIR199A2; hsa-miR-9 -> MIR9-1,
    MIR9-2, MIR9-3) is left ambiguous, like any other ambiguous symbol. Returns None
    when the name is not a miRNA name or nothing matches."""
    global MIR_FAMILY
    m = MIRNA.match(name.strip())
    if not m:
        return None
    kind, core, arm, star = m.group(1).upper(), m.group(2).upper(), m.group(3), m.group(4)
    tag = "mirna-arm" if (arm or star) else "mirna"
    s, h, route = symbol_to_hgnc(f"HSA-{kind}-{core}")
    if s or "ambiguous" in route:
        return s, h, f"{tag}->{route}"
    approved = HGNC[0]
    base = ("MIRLET" if kind == "LET" else "MIR") + core
    for cand in dict.fromkeys((base, re.sub(r"(?<=[A-Z])-(?=\d+$)", "", base))):
        if cand in approved:
            return approved[cand] + (f"{tag}->precursor",)
    if MIR_FAMILY is None:
        # precursor families: hsa-mir-29b -> {MIR29B1, MIR29B2}, from aliases of the
        # form hsa-mir-29b-1 and from the approved symbols themselves
        MIR_FAMILY = {}
        for a, hits in HGNC[2].items():
            f = re.fullmatch(r"HSA-(MIR|LET)-(.+)-\d+", a)
            if f:
                MIR_FAMILY.setdefault(("MIRLET" if f.group(1) == "LET" else "MIR") + f.group(2), set()).update(hits)
        for k, v in approved.items():
            f = re.fullmatch(r"(MIR.*[A-Z])-?\d+|(MIR.*\d)-\d+", k)
            if f:
                MIR_FAMILY.setdefault(f.group(1) or f.group(2), set()).add(v)
    fam = MIR_FAMILY.get(base, set())
    if len(fam) == 1:
        return next(iter(fam)) + (f"{tag}->precursor",)
    if len(fam) > 1:
        return ("", "", f"{tag}->ambiguous:" + ",".join(sorted(x[0] for x in fam)))
    return None


def name_candidates(name):
    """Spellings to try for a printed gene name that does not resolve as written:
    footnote marks and a trailing "(includes others)" removed; "IPLA2(\u03b3)" read as
    IPLA2GAMMA / IPLA2G (the bare stem would name a different family member)."""
    t = re.sub(r"\s*\(includes others\)\s*$", "", name.strip(), flags=re.I)
    t = t.strip(FOOTNOTE_MARKS + " ")
    m = re.fullmatch(r"(\S+?)\s*\((.)\)", t)
    if m and m.group(2) in GREEK:
        return [m.group(1) + GREEK[m.group(2)], m.group(1) + GREEK[m.group(2)][0]]
    return [t] if t and t != name.strip() else []


def name_to_hgnc(name, year=None):
    """A printed gene name -> (symbol, hgnc_id, route): as written, then as a miRNA
    name, then with footnote marks and qualifiers removed (route prefixed 'cleaned->').
    year is the paper's, for symbol_to_hgnc's previous-vs-alias rule. The printed name itself is kept in gene_reported."""
    if is_missing(name):
        return "", "", "missing"
    r = symbol_to_hgnc(name, year)
    if r[0] or "ambiguous" in r[2]:
        return r
    r = mirna_to_hgnc(name)
    if r:
        return r
    for cand in name_candidates(name):
        s, h, route = symbol_to_hgnc(cand, year)
        if s or "ambiguous" in route:
            return s, h, f"cleaned->{route}"
        r = mirna_to_hgnc(cand)
        if r:
            return r[0], r[1], f"cleaned->{r[2]}"
    return "", "", "unmapped"


# X08: bare BAC/cosmid clone-library prefixes, printed as a whole gene name when a clone-based name (RP11-...)
# was cut short. HGNC keeps several as previous symbols of retinitis pigmentosa loci (RP11 -> PRPF31), so they are
# read as no name. 28822801 FCL/LCL print 'RP11' on 3 rows (extraction_audit_2026-09-28/meaning_changed_rows.tsv).
# RP1, RP2 and RP9 are left out: they are current HGNC approved symbols, printed as those genes in 20+ lists.
CLONE_LIBRARY_PREFIXES = ({f"RP{i}" for i in range(1, 14)} - {"RP1", "RP2", "RP9"}) | {"CTD", "CTC", "CTB", "KB", "GS1"}


def is_clone_prefix(v):
    return (v or "").strip().upper() in CLONE_LIBRARY_PREFIXES


def normalize_gene(row, spec):
    """Return (symbol, hgnc_id, route). One chokepoint for every list: the primary steps, then the X09-X11
    fallbacks for rows they leave without a symbol (module docstring)."""
    s, h, route = _normalize_primary(row, spec)
    if s or route.startswith("override->") or route.startswith(GEMMA_NO_GENE):
        # Gemma's 'no gene' (or several genes, or genes HGNC lacks) is final: no fallback to an ID route or the printed
        # name (Paul, 2026-09-30)
        return s, h, route
    rid = (row.get("id_reported") or "").strip()
    sym = (row.get("gene_reported") or "").strip()
    rid = "" if is_missing(rid) else rid
    sym = "" if is_missing(sym) or is_clone_prefix(sym) else sym
    if route.startswith("platform-multi:"):
        # a probe the platform annotates to several genes: a fallback (Gemma's annotation of the probe, the resolver on
        # the printed name) may pick one of them, never a gene outside them. The platform reading wins over the printed
        # name (Paul, 2026-09-30); 9 rows in 6 lists had taken a printed gene the platform does not name
        # (step2_symbol_audit_2026-09-30/probe_vs_platform_rows.tsv, kind multi_excl).
        allowed = platform_multi_genes(route)
        hit = fallback_gene(rid, sym, spec)
        return hit if hit and hit[0] in allowed else (s, h, route)
    if route.startswith("previous-alias-ambiguous:"):
        # the row's own identifier may settle it, when it names one of the candidates; the name routes may not
        cands = set(route.split(":", 1)[1].split(","))
        hit = fallback_gene(rid, "", spec)  # no printed name: ID routes only
        return hit if hit and hit[0] in cands else (s, h, route)
    hit = fallback_gene(rid, sym, spec)
    return hit if hit else (s, h, route)


def platform_multi_genes(route):
    """The approved symbols of the genes a 'platform-multi:A,B;paper:X' route names (the platform's symbols as written,
    and their current HGNC symbols)."""
    syms = route.split(";")[0].split(":", 1)[1].split(",")
    out = set(syms) | {symbol_to_hgnc(x, official=True)[0] for x in syms}
    for x in syms:
        m = re.fullmatch(r"LOC(\d+)", x)
        if m:
            # a LOC placeholder the platform still carries, under its current HGNC symbol (LOC401320 = GARS1-DT)
            FALLBACK_FILES.add(NCBI_GENE_INFO)
            out.add(xrefs()["ncbi"].get(m.group(1), ("",))[0])
    return out - {""}


def _geo_platform(spec, rid, sym, year):
    """The GEO .annot reading of a reported probe ID (the platform step before 2026-09-30's Gemma-first rule; still the
    reading for platforms Gemma lacks): (symbol, hgnc_id, route), or None when the .annot gives the probe no gene."""
    table = platform(spec["platform"], spec.get("platform_id_column"))
    probes = platform_probes(spec["platform"], spec.get("platform_id_column"), rid)
    syms = sorted({s for p in probes for s in table[p]})
    # the printed name's own gene, recorded beside a platform reading that names another gene (fix 5 of the
    # 2026-09-30 Step 2 audit: 162 rows where paper and platform disagree)
    named, named_id, named_route = name_to_hgnc(sym, year) if sym and syms else ("", "", "")
    # a printed name HGNC gives several genes for is recorded with its candidates
    named_cands = set(named_route.split(":", 1)[1].split(",")) if not named and "ambiguous:" in named_route else set()
    note = f";paper:{named}" if named else (";paper-ambiguous:" + ",".join(sorted(named_cands)) if named_cands else "")
    if len(syms) == 1:
        s, h, route = symbol_to_hgnc(syms[0], official=True)
        if not s and sym:
            # the platform's symbol is not in HGNC (e.g. LOC101928973) but the paper
            # prints an approved one (LINC01732): take the printed name
            if named:
                return named, named_id, f"platform-unmapped->paper-name->{named_route}"
        return s, h, f"platform->{route}" + (note if s and named != s and s not in named_cands else "")
    if len(syms) > 1:
        # the paper's own gene name picks among the probe's genes, when it names one of them; a name HGNC gives
        # several genes for (previous vs alias, fix 4) picks the one of its candidates the probe carries
        pgenes = platform_multi_genes("platform-multi:" + ",".join(syms))
        if named and named in pgenes:
            return named, named_id, "platform-multi->paper-name"
        cands = named_cands & pgenes
        if len(cands) == 1:
            return HGNC[0][next(iter(cands)).upper()] + ("platform-multi->paper-name-candidate",)
        return "", "", "platform-multi:" + ",".join(syms) + note
    return None


def reading_notes(gene, geo, sym, year):
    """map_route notes beside a Gemma reading: ';geo:X' (or ';geo:none', ';geo-multi:A,B') where the GEO .annot reads
    the probe otherwise, ';paper:X' (';paper-ambiguous:A,B') where the printed name means another gene."""
    out = ""
    if geo is None:
        geo_gene = ""
    else:
        geo_gene = geo[0]
    if geo_gene != gene:
        if geo_gene:
            out += f";geo:{geo_gene}"
        elif geo is not None and geo[2].startswith("platform-multi:"):
            out += ";geo-multi:" + geo[2].split(";")[0].split(":", 1)[1]
        elif gene:
            out += ";geo:none"
    if sym:
        named, _, named_route = name_to_hgnc(sym, year)
        cands = set(named_route.split(":", 1)[1].split(",")) if not named and "ambiguous:" in named_route else set()
        if named and named != gene:
            out += f";paper:{named}"
        elif cands and gene not in cands:
            out += ";paper-ambiguous:" + ",".join(sorted(cands))
    return out


def _normalize_primary(row, spec):
    rid = (row.get("id_reported") or "").strip()
    sym = (row.get("gene_reported") or "").strip()
    rid = "" if is_missing(rid) else rid
    sym = "" if is_missing(sym) else sym
    year = paper_year(spec)
    override = spec.get("symbol_overrides", {}).get(rid) or spec.get("symbol_overrides", {}).get(sym)
    if override:
        s, h, route = symbol_to_hgnc(override, official=True)  # a curator's choice of symbol
        return s, h, f"override->{route}"
    clone = is_clone_prefix(sym)
    if clone:
        sym = ""  # X08
    if spec.get("platform") and rid:
        geo = _geo_platform(spec, rid, sym, year)
        gem = gemma_first(spec, rid)
        if gem is not None:
            # Gemma's annotation decides (Paul, 2026-09-30); the GEO and printed readings are recorded where they differ
            return gem[0], gem[1], gem[2] + reading_notes(gem[0], geo, sym, year)
        if geo is not None:
            absent = gpl_of(spec) and gemma_annot(gpl_of(spec))[1] is not None
            return geo[0], geo[1], geo[2] + (";gemma:absent" if absent else "")
    if rid.startswith("ENSG"):
        hit = HGNC[3].get(rid.split(".")[0])
        if hit:
            return hit + ("ensembl",)
    if sym:
        if spec.get("provenance") == "gemma":
            # Gemma's gene_symbol column is its platform annotation's official symbol, not a printed name
            s, h, route = symbol_to_hgnc(sym, official=True)
            if s:
                return s, h, route
        return name_to_hgnc(sym, year)
    return "", "", "clone-library-name" if clone else "unmapped"


# ---------------------------------------------------------------- fallbacks X09-X11
# Rows the steps above leave without a symbol. Routes, order and rules follow
# analyses/list_extraction_audit/a1_unmapped.py (recover(), name_candidates()), whose output
# unmapped_rows.tsv is the evidence for X09-X11 in proposed_fixes.json. All inputs are local files, read only.
GENE_NAMES_DIR = os.path.join(HERE, "gene_names")
NCBI_GENE_INFO = os.path.join(GENE_NAMES_DIR, "sources", "Homo_sapiens.gene_info.gz")
GEMMA_ANNOT_DIR = os.path.join(HERE, "extraction_audit_2026-09-28", "platform_annots")
PLATFORM_DIR = os.path.join(DATA, "platforms")
# files the fallbacks read while building the current list (main() resets it per list and names them in the header)
FALLBACK_FILES = set()
_FALLBACK_SHA = {}


def fallback_sha(path):
    """sha256 of a fallback source, computed once per run (the NCBI gene_info file is large)"""
    if path not in _FALLBACK_SHA:
        _FALLBACK_SHA[path] = sha256(path)
    return _FALLBACK_SHA[path]

_RESOLVER = None


def resolver():
    """gene_names.Resolver (HGNC + NCBI Gene + UniProt + curated), loaded on first use (~20 s)."""
    global _RESOLVER
    if _RESOLVER is None:
        if GENE_NAMES_DIR not in sys.path:
            sys.path.insert(0, GENE_NAMES_DIR)
        import gene_names
        _RESOLVER = gene_names.Resolver()
    FALLBACK_FILES.update(os.path.join(GENE_NAMES_DIR, f) for f in ("gene_names_table.tsv", "gene_names_curated.tsv"))
    return _RESOLVER


# X09: the resolver tiers the fix names (gene_names' 'identifier' tier is left out)
RESOLVER_TIERS = {"exact_approved", "unique_alias", "curated"}

_XREF = None


def xrefs():
    """HGNC and NCBI gene_info cross-references, loaded on first use:
    by_id      HGNC ID -> approved symbol
    entrez     NCBI Gene ID -> approved symbol (HGNC entrez_id)
    acc        RefSeq / ENA accession, version stripped, upper case -> set of approved symbols (HGNC refseq_accession, ena)
    ncbi       NCBI Gene ID -> (approved symbol via gene_info's HGNC dbXref or '', NCBI symbol, type_of_gene)
    ncbi_ens   Ensembl gene ID -> set of NCBI Gene IDs (gene_info dbXrefs)"""
    global _XREF
    if _XREF is None:
        by_id, entrez, acc = {}, {}, {}
        with open(HGNC_FILE, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] != "Approved":
                    continue
                by_id[r["hgnc_id"]] = r["symbol"]
                if r["entrez_id"]:
                    entrez[r["entrez_id"]] = r["symbol"]
                for f in ("refseq_accession", "ena"):
                    for a in filter(None, r[f].strip('"').split("|")):
                        acc.setdefault(a.split(".")[0].upper(), set()).add(r["symbol"])
        ncbi, ncbi_ens = {}, {}
        with gzip.open(NCBI_GENE_INFO, "rt", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                hg = re.search(r"HGNC:(HGNC:\d+)", r["dbXrefs"])
                ncbi[r["GeneID"]] = ((by_id.get(hg.group(1)) or "") if hg else "", r["Symbol"], r["type_of_gene"])
                for e in re.findall(r"Ensembl:(ENSG\d+)", r["dbXrefs"]):
                    ncbi_ens.setdefault(e, set()).add(r["GeneID"])
        _XREF = {"by_id": by_id, "entrez": entrez, "acc": acc, "ncbi": ncbi, "ncbi_ens": ncbi_ens}
    return _XREF


UNIPROT_FILE = os.path.join(GENE_NAMES_DIR, "sources", "uniprot_human_reviewed.tsv")
# UniProt accession, optional isoform suffix (UniProt's documented accession format)
UNIPROT_ACC = re.compile(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?")
_UNIPROT = None


def uniprot_genes():
    """X10 (fix 2, 2026-09-30): UniProt accession -> set of approved symbols, from HGNC uniprot_ids and the HGNC
    cross-reference of the local reviewed-human UniProt file (as the 2026-09-30 Step 2 audit read them)."""
    global _UNIPROT
    if _UNIPROT is None:
        by_id, table = xrefs()["by_id"], {}
        with open(HGNC_FILE, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] == "Approved":
                    for u in filter(None, r["uniprot_ids"].strip('"').split("|")):
                        table.setdefault(u, set()).add(r["symbol"])
        with open(UNIPROT_FILE, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for h in re.findall(r"HGNC:\d+", r.get("HGNC") or ""):
                    if h in by_id:
                        table.setdefault(r["Entry"], set()).add(by_id[h])
        _UNIPROT = table
    FALLBACK_FILES.add(UNIPROT_FILE)
    return _UNIPROT


def uniprot_route(x):
    """A cell of UniProt accessions (one, or a protein group 'Q92599-2;A6NMH6;Q92599') -> one HGNC gene, when the
    accessions HGNC or UniProt cross-reference all name the same gene; accessions with no cross-reference (unreviewed
    entries) are passed over."""
    parts = [p for p in re.split(r"[;,|\s]+", x) if p]
    if not parts or not all(UNIPROT_ACC.fullmatch(p) for p in parts):
        return None
    table = uniprot_genes()
    syms = set().union(*(table.get(p.split("-")[0], set()) for p in parts))
    hit = one_hgnc_gene(sorted(syms)) if syms else None
    return hit + ("uniprot->hgnc",) if hit else None


def hgnc_symbol(sym):
    """an approved symbol from a cross-reference -> (symbol, hgnc_id), via symbol_to_hgnc"""
    s, h, _ = symbol_to_hgnc(sym, official=True)
    return (s, h) if s else None


def one_hgnc_gene(symbols):
    """(symbol, hgnc_id) when every symbol maps to HGNC and all name the same gene, else None."""
    hits = [hgnc_symbol(s) for s in symbols]
    if hits and all(hits) and len(set(hits)) == 1:
        return hits[0]
    return None


_GEMMA = {}


def gemma_annot(gpl):
    """probe -> (gene symbols, NCBI Gene IDs) from Gemma's annotation of the GPL, or None when no file is on disk.
    The gemma2 file (GEMMA_ANNOT_DIR, fetched by extraction_audit_2026-09-28/fetch_gemma_annots.py) first; a production
    export in data/platforms only when gemma2's is absent (GPL80's production export, 2026-09-23, gives 6,524 of 7,129
    probes no gene against gemma2's 1,424)."""
    if gpl not in _GEMMA:
        path = os.path.join(GEMMA_ANNOT_DIR, f"{gpl}.gemma.annot.gz")
        if not os.path.exists(path) and os.path.exists(os.path.join(PLATFORM_DIR, f"{gpl}.gemma-prod.annot.gz")):
            path = os.path.join(PLATFORM_DIR, f"{gpl}.gemma-prod.annot.gz")
        table = None
        if os.path.exists(path):
            table = {}
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                    table[r["ElementName"]] = ([s for s in (r["GeneSymbols"] or "").split("|") if s],
                                               [s for s in (r["NCBIids"] or "").split("|") if s])
        _GEMMA[gpl] = (path, table)
    return _GEMMA[gpl]


def gemma_gene_to_hgnc(sym, ncbi_id):
    """One gene of Gemma's annotation -> (symbol, hgnc_id), or None when it is not an HGNC gene: its NCBI Gene ID
    through HGNC's entrez_id, else gene_info's HGNC cross-reference; the symbol (an official one) only when the
    annotation gives no NCBI ID or the ID has no HGNC entry."""
    xr = xrefs()
    if ncbi_id:
        s = xr["entrez"].get(ncbi_id) or xr["ncbi"].get(ncbi_id, ("",))[0]
        if not xr["entrez"].get(ncbi_id):
            FALLBACK_FILES.add(NCBI_GENE_INFO)
        hit = hgnc_symbol(s) if s else None
        if hit:
            return hit
    return hgnc_symbol(sym) if sym else None


def gpl_of(spec):
    m = re.search(r"GPL\d+", os.path.basename(spec.get("platform") or ""))
    return m.group(0) if m else ""


def gemma_first(spec, rid):
    """Gemma-first probe mapping (Paul, 2026-09-30: "we store alignments, and then decide the gene (if any)"). For a
    reported probe ID on a platform Gemma annotates: (symbol, hgnc_id, route) from Gemma's annotation of the probe,
    or None when Gemma has no annotation file for the GPL or no probe matching the ID (then the GEO .annot reads it).
      one HGNC gene         that gene                                  'gemma-platform:GPL570'
      no gene               unmapped, no fallback to GEO or the name   'gemma-platform-nogene:GPL570'
      genes HGNC lacks      unmapped                                   'gemma-platform-nohgnc:GPL570:LOC105369'
      several genes         unmapped unless they are one HGNC gene     'gemma-platform-multi:GPL570:ACTB,ACTG1'
    Several probes matching one ID (GPL80's accession-keyed probes, a reported ID without its _at suffix) pool their
    genes. Lists keyed by another column of the GEO .annot (platform_id_column, e.g. the clone's GenBank accession)
    reach Gemma's annotation through the GEO .annot's probe IDs for that value."""
    gpl = gpl_of(spec)
    if not gpl or not rid:
        return None
    path, table = gemma_annot(gpl)
    if table is None:
        return None
    if spec.get("platform_id_column"):
        probes = {p for p in geo_probes_by_column(spec["platform"], spec["platform_id_column"], rid) if p in table}
    else:
        probes = probe_matches(("gemma", gpl), table, rid)
    if not probes:
        return None
    FALLBACK_FILES.add(path)
    raw, genes, nohgnc = [], set(), []
    for p in sorted(probes):
        syms, ids = table[p]
        for i, s in enumerate(syms):
            if s in raw:
                continue
            raw.append(s)
            hit = gemma_gene_to_hgnc(s, ids[i] if len(ids) == len(syms) else "")
            if hit:
                genes.add(hit)
            else:
                nohgnc.append(s)
    if not raw:
        return "", "", f"gemma-platform-nogene:{gpl}"
    if len(genes) == 1 and not nohgnc:
        return next(iter(genes)) + (f"gemma-platform:{gpl}",)
    if not genes:
        return "", "", f"gemma-platform-nohgnc:{gpl}:" + ",".join(nohgnc)
    return "", "", f"gemma-platform-multi:{gpl}:" + ",".join(sorted({g[0] for g in genes}) + nohgnc)


GEMMA_NO_GENE = ("gemma-platform-nogene:", "gemma-platform-nohgnc:", "gemma-platform-multi:")


_GEO_COLUMN = {}


def geo_column_index(col):
    """X10: value of a GEO .annot column ('GenBank Accession', 'UniGene ID'), version stripped, upper case -> set of
    gene symbols, over every GEO .annot in data/platforms/."""
    files = sorted(glob.glob(os.path.join(PLATFORM_DIR, "*.geo.annot.gz")))
    if col not in _GEO_COLUMN:
        idx = {}
        for p in files:
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                header = None
                for line in fh:
                    if line.startswith(("#", "!", "^")):
                        continue
                    row = line.rstrip("\n").split("\t")
                    if header is None:
                        header = row
                        if col not in header:
                            break
                        ci, si = header.index(col), header.index("Gene symbol")
                        continue
                    if len(row) > max(ci, si) and row[ci] and row[si]:
                        for v in row[ci].split("///"):
                            idx.setdefault(v.strip().upper().split(".")[0], set()).update(s for s in row[si].split("///") if s)
        _GEO_COLUMN[col] = idx
    FALLBACK_FILES.update(files)
    return _GEO_COLUMN[col]


# X11: identifier shapes looked up in a platform annotation (a1_unmapped.PROBE_PATTERNS, probe kinds only)
PROBE_ID = re.compile(r"^(AFFX-\S+|\d+(_PM)?(_[a-z]+)?_(s_|x_|a_|g_|f_|i_|r_|b_)?at|ILMN_\d+|A_\d+_P\d+|GE\d{5,}|\d{7,8})$", re.I)
ENSG_ID = re.compile(r"ENSG\d{11}(\.\d+)?")
ACCESSION = re.compile(r"(GB:)?([A-Z]{1,2}_?\d{5,9})(\.\d+)?", re.I)
CLONE_GENE_NAME = re.compile(r"[A-Z]{1,3}\d{5,6}\.\d+", re.I)
UNIGENE = re.compile(r"Hs\.\d+", re.I)


def entrez_list(spec):
    """X10: NCBI Gene IDs are read only in lists that are a Gemma reanalysis or name an Entrez/NCBI/gene-id column
    (a1_unmapped.main: kinds 'entrez')."""
    return spec.get("provenance") == "gemma" or bool(re.search(r"entrez|ncbi|gene_?id", json.dumps(spec.get("columns", {})), re.I))


def gemma_route(ids, spec):
    """X11: a probe ID -> one HGNC gene through Gemma's annotation of the manifest's GEO platform."""
    m = re.search(r"GPL\d+", os.path.basename(spec.get("platform") or ""))
    if not m or spec.get("platform_id_column"):
        return None
    gpl = m.group(0)
    for x in ids:
        if not PROBE_ID.match(x):
            continue
        path, table = gemma_annot(gpl)
        if table is None:
            return None
        cand = [p for p in (x, x.lower(), x.upper()) if p in table]
        if not cand:
            continue
        FALLBACK_FILES.add(path)
        syms, ncbi_ids = table[cand[0]]
        hit = one_hgnc_gene(syms) if syms else None
        if hit:
            return hit + (f"gemma-platform:{gpl}",)
        xr = xrefs()
        ents = {xr["entrez"].get(n) or xr["ncbi"].get(n, ("",))[0] for n in ncbi_ids} - {""}
        if len(ents) == 1 and len(ncbi_ids) == 1:
            hit = hgnc_symbol(next(iter(ents)))
            if hit:
                return hit + (f"gemma-platform-ncbi:{gpl}",)
    return None


def identifier_route(ids, sym, spec):
    """X10: UniProt / Ensembl / NCBI Gene / RefSeq-GenBank / UniGene identifiers -> one HGNC gene."""
    for x in ids:
        if x != sym:
            # UniProt accessions in the reported-ID column (fix 2, 2026-09-30: 26485546 row 228 prints '2008-09-01' as
            # the gene beside Q92599, SEPTIN8). Tried first because a six-character accession (P12345) also has the
            # shape of a GenBank accession.
            hit = uniprot_route(x)
            if hit:
                return hit
        if ENSG_ID.fullmatch(x):
            e = x.split(".")[0]
            hit = HGNC[3].get(e)
            if hit:
                # an Ensembl ID printed in the gene-name column (the reported-ID column is handled by the primary step)
                return hit + ("ensembl",)
            xr = xrefs()
            FALLBACK_FILES.add(NCBI_GENE_INFO)
            syms = {xr["ncbi"][g][0] for g in xr["ncbi_ens"].get(e, set()) if xr["ncbi"][g][0]}
            if len(syms) == 1:
                hit = hgnc_symbol(next(iter(syms)))
                if hit:
                    return hit + ("ensembl->ncbi-gene-info",)
        elif re.fullmatch(r"\d{1,9}", x) and entrez_list(spec):
            xr = xrefs()
            s = xr["entrez"].get(x)
            route = "entrez->hgnc"
            if not s:
                s, route = xr["ncbi"].get(x, ("",))[0], "entrez->ncbi-gene-info"
                FALLBACK_FILES.add(NCBI_GENE_INFO)
            hit = hgnc_symbol(s) if s else None
            if hit:
                return hit + (route,)
        elif ACCESSION.fullmatch(x) and not (x == sym and CLONE_GENE_NAME.fullmatch(x)):
            # a printed AC123456.3 is an Ensembl clone-based gene name, not an mRNA accession: the clone holds many
            # genes, so HGNC's 'ena' column would name an arbitrary one of them (a1_unmapped.recover)
            a = re.sub(r"^GB:", "", x, flags=re.I).split(".")[0].upper()
            syms, route = set(xrefs()["acc"].get(a, set())), "accession->hgnc"
            if not syms:
                syms, route = set(geo_column_index("GenBank Accession").get(a, set())), "accession->geo-annot"
            hit = one_hgnc_gene(sorted(syms)) if syms else None
            if hit:
                return hit + (route,)
        elif UNIGENE.fullmatch(x):
            syms = geo_column_index("UniGene ID").get(x.upper(), set())
            hit = one_hgnc_gene(sorted(syms)) if syms else None
            if hit:
                return hit + ("unigene->geo-annot",)
    return None


# X09: spreadsheet dates made from gene symbols. The candidate prefixes are the ones the fix names, plus a bare SEP,
# which a1_unmapped.MONTH_GENE also tries: a date that SEP<n> and SEPT<n> read as different genes stays unmapped
# (33805976 st2_p05 '43718' = 10-Sep: SEPTIN10 or MED31, whose NCBI synonym is SEP10; unmapped_rows.tsv).
MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
# MARC<n>: previous symbols of MTARC1/MTARC2, which Excel also turns into March dates (fix 1, 2026-09-30)
MONTH_GENE = {"SEP": ["SEPT", "SEPTIN", "SEP"], "MAR": ["MARCH", "MARCHF", "MARC"], "DEC": ["DEC"]}


def excel_decode(v):
    """(month, day, year or None) when v is a spreadsheet date: 2019-03-05 (00:00:00), 1-Sep, Mar-02, or an Excel
    serial day number (30000-49999, 1982-2036) falling in September, March or December; else None
    (a1_unmapped.excel_decode, plus the year)."""
    t = (v or "").strip()
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})( 00:00:00)?", t)
    if m:
        return MONTHS[int(m.group(2)) - 1], int(m.group(3)), int(m.group(1))
    m = re.fullmatch(r"(\d{1,2})[-. ]([A-Za-z]{3})", t)
    if m and m.group(2).upper() in MONTHS:
        return m.group(2).upper(), int(m.group(1)), None
    m = re.fullmatch(r"([A-Za-z]{3})[-. ](\d{1,2})", t)
    if m and m.group(1).upper() in MONTHS:
        return m.group(1).upper(), int(m.group(2)), None
    m = re.fullmatch(r"([34]\d{4})(\.0)?", t)
    if m:
        d = datetime.date(1899, 12, 30) + datetime.timedelta(days=int(m.group(1)))
        if MONTHS[d.month - 1] in MONTH_GENE:
            return MONTHS[d.month - 1], d.day, d.year
    return None


def excel_numbers(dec):
    """The gene numbers a decoded date can stand for: its day, and, for a full date or serial on the 1st, also the
    last two digits of its year. SEPT8 read as September 2008 is stored as 2008-09-01 (26485546 row 228, beside UniProt
    Q92599 = SEPTIN8; 19255580 serial 38412 = 2005-03-01 on a MARCHF5 probe), while 30053206 and 41639050 carry day-1
    dates whose IDs confirm the day (2026-09-30 Step 2 audit)."""
    mon, day, year = dec
    nums = [day]
    if day == 1 and year is not None and year % 100 not in (0, 1):
        nums.append(year % 100)
    return nums


def resolve_name(name):
    """X09: gene_names.Resolver on a printed name -> (symbol, hgnc_id, route) for the tiers the fix names, else None."""
    r = resolver().resolve(name)
    if r.tier in RESOLVER_TIERS and r.symbol:
        hit = hgnc_symbol(r.symbol)
        if hit:
            return hit + (f"resolver:{r.tier}:{r.route}",)
    return None


def name_route(sym):
    """X09: LOC<GeneID>, spreadsheet dates, then the resolver."""
    m = re.fullmatch(r"LOC(\d+)", sym)
    if m:
        s = xrefs()["ncbi"].get(m.group(1), ("",))[0]
        FALLBACK_FILES.add(NCBI_GENE_INFO)
        hit = hgnc_symbol(s) if s else None
        if hit:
            return hit + ("loc->ncbi-gene-info",)
    if sym.endswith("-"):
        return None  # a truncated name
    dec = excel_decode(sym)
    if dec:
        hits = {}
        for n in excel_numbers(dec):
            for pre in MONTH_GENE.get(dec[0], []):
                r = resolve_name(f"{pre}{n}")
                if r:
                    hits.setdefault(r[:2], f"{pre}{n}")
        if len(hits) == 1:
            (s, h), read_as = next(iter(hits.items()))
            return s, h, f"excel-date:{read_as}"
        if len(hits) > 1:
            # a date that can mean several genes stays unmapped (the row's own ID, tried before this, did not settle it)
            return "", "", "excel-date-ambiguous:" + ",".join(sorted(s for s, _ in hits))
        return None
    return resolve_name(sym)


def printed_name_genes(sym):
    """Genes the printed name itself can mean: the resolver's symbol or its ambiguous candidates, or a LOC number's
    NCBI gene. Empty when the name carries no gene information (a1_unmapped.name_candidates)."""
    if not sym or excel_decode(sym):
        return set()
    m = re.fullmatch(r"LOC(\d+)", sym)
    if m:
        return {xrefs()["ncbi"].get(m.group(1), ("",))[0]} - {""}
    r = resolver().resolve(sym)
    if r.resolved:
        return {r.symbol}
    return {c.get("symbol") for c in r.candidates} - {None, ""}


def fallback_gene(rid, sym, spec):
    """(symbol, hgnc_id, route) from the X11, X10 and X09 fallbacks in that order, or None. rid and sym are the
    reported ID and printed name, with missing values and clone-library prefixes (X08) already blanked."""
    ids = [x for x in (rid, sym) if x]
    for route in (lambda: gemma_route(ids, spec), lambda: identifier_route(ids, sym, spec)):
        hit = route()
        if hit:
            # an ID route that disagrees with the printed name is not taken (a1_unmapped: status 'conflict')
            names = printed_name_genes(sym)
            if names and hit[0] not in names:
                return None
            return hit
    return name_route(sym) if sym else None


# ---------------------------------------------------------------- readers
def dedupe_header(cells):
    """Repeated column names (one block per cohort or disorder, side by side) become
    "name", "name.2", "name.3", ... in every reader, so the first block keeps the plain
    name and a later block can never overwrite it in dict(zip())."""
    seen, header = {}, []
    for c in cells:
        seen[c] = seen.get(c, 0) + 1
        header.append(c if seen[c] == 1 else f"{c}.{seen[c]}")
    return header


def as_dicts(header, grid, **extra):
    """Rows as dicts; short (ragged) rows are padded to the header width, so a missing
    trailing block reads as empty rather than shifting which block a name refers to."""
    w = len(header)
    return [dict(zip(header, list(r) + [None] * (w - len(r))), **extra) for r in grid]


def read_rows(spec):
    path = os.path.join(DATA, spec["file"])
    kind = spec["reader"]
    member = spec.get("member")
    if member:
        # the table is one file inside a supplementary zip ("member" = its name in the archive)
        import io, zipfile
        data = zipfile.ZipFile(path).read(member)
        path = io.BytesIO(gzip.decompress(data) if member.endswith(".gz") else data)
    if kind in ("xlsx", "xls"):
        # "sheets" (a list) stacks several same-layout sheets into one list,
        # e.g. one sheet per cell type when the paper counts DEGs across them
        rows = []
        for name in spec.get("sheets") or [spec["sheet"]]:
            if kind == "xlsx":
                import openpyxl
                ws = openpyxl.load_workbook(path, read_only=True, data_only=True)[name]
                # some writers store <dimension ref="A1"/>, which read-only mode takes as the
                # sheet's extent (Ruzicka 2024 Data S4 read as one cell per sheet)
                ws.reset_dimensions()
                grid = [list(r) for r in ws.iter_rows(values_only=True)]
            else:
                import xlrd
                sh = (xlrd.open_workbook(file_contents=path.getvalue()) if member else xlrd.open_workbook(path)).sheet_by_name(name)
                grid = [sh.row_values(i) for i in range(sh.nrows)]
            h = spec.get("header_row", 1) - 1
            header = dedupe_header([str(c).strip() if c is not None else "" for c in grid[h]])
            # "last_row" (1-based, inclusive) ends the table early, for sheets that hold
            # several tables one below another
            end = spec.get("last_row") or len(grid)
            # "column_offset": {"<header name>": k} for a column printed k rows above the
            # rest of its row, the header cell included (Narayan 2008 Table 1a: its first
            # probe ID sits in the header row, so every probe belongs to the gene below it)
            for col, k in (spec.get("column_offset") or {}).items():
                j = header.index(col)
                grid = [list(r) for r in grid]
                vals = [grid[i - k][j] if i - k >= h else "" for i in range(len(grid))]
                for i in range(h + 1, end):
                    grid[i][j] = vals[i]
            rows += as_dicts(header, [r for r in grid[h + 1:end] if any(v not in (None, "") for v in r)], _sheet=name)
        return rows
    if kind == "docx":
        # a Word table: spec["table"] is its index in the document
        import docx
        t = docx.Document(path).tables[spec["table"]]
        grid = [[c.text.strip() for c in r.cells] for r in t.rows]
        h = spec.get("header_row", 1) - 1
        header = dedupe_header(grid[h])
        return as_dicts(header, [r for r in grid[h + 1:] if any(r)])
    if kind in ("tsv", "csv"):
        # plain or gzip-compressed (full summary statistics often ship as .gz)
        if member:
            fh = io.TextIOWrapper(path, encoding=spec.get("encoding", "utf-8"), errors="replace")
        else:
            opener = gzip.open if path.endswith(".gz") else open
            fh = opener(path, "rt", encoding=spec.get("encoding", "utf-8"), errors="replace")
        with fh:
            lines = [l for l in fh if not l.startswith("#")]
        lines = lines[spec.get("header_row", 1) - 1:]
        grid = list(csv.reader(lines, delimiter="\t" if kind == "tsv" else ","))
        if not grid:
            return []
        # blank lines are skipped, as csv.DictReader did
        return as_dicts(dedupe_header(grid[0]), [r for r in grid[1:] if r])
    if kind == "probe_text":
        # probe IDs in document order; direction from the section headings
        t = open(path, encoding="utf-8").read()
        marks = sorted((t.lower().index(k.lower()), v) for k, v in spec["sections"].items())
        rows = []
        for m in re.finditer(spec["probe_regex"], t):
            d = [v for pos, v in marks if pos < m.start()]
            rows.append({"probe": m.group(0), "direction": d[-1] if d else ""})
        return rows
    if kind == "inline":
        return spec["rows"]
    raise ValueError(kind)


def num(v):
    try:
        return float(str(v).replace("×10−", "e-").replace("−", "-"))
    except (TypeError, ValueError):
        return None


def invert_effect(v, center=0):
    """Case relative to control from control relative to case: a log ratio or
    difference changes sign; a plain ratio (center 1) becomes its reciprocal."""
    x = num(v)
    if x is None:
        return v
    if center == 1:
        return "" if x == 0 else f"{1 / x:.6g}"
    t = str(v).strip().replace("\u2212", "-")
    if x == 0:
        return t
    return t[1:] if t.startswith("-") else "-" + t


def build(spec):
    rows = read_rows(spec)
    cols = spec["columns"]
    for f in spec.get("filters", []):
        if "equals" in f:
            rows = [r for r in rows if str(r.get(f["col"], "")).strip() == f["equals"]]
        elif f.get("nonempty"):
            rows = [r for r in rows if str(r.get(f["col"]) or "").strip()]
        else:
            val = (lambda v: None if v is None else abs(v)) if f.get("abs") else (lambda v: v)
            rows = [r for r in rows if num(r.get(f["col"])) is not None and OPS[f["op"]](val(num(r.get(f["col"]))), f["value"])]
    out = []
    for r in rows:
        # a column spec may be a list: the values are joined with " | " (e.g. "group": ["cell_type", "timepoint"]
        # for tables stacked over cell type and time point)
        rec = {k: " | ".join(str(r.get(x)).strip() for x in (c if isinstance(c, list) else [c]) if r.get(x) not in (None, ""))
               for k, c in cols.items()}
        d = spec.get("direction", {})
        if d.get("from") == "sign":
            # "center": 1 for ratios; "text_values" for non-numeric cells the paper defines (e.g. "large")
            raw = str(r.get(d["col"], "")).strip()
            v = num(raw)
            c = d.get("center", 0)
            rec["direction"] = d.get("text_values", {}).get(raw, "") if v is None else ("up" if v > c else "down" if v < c else "")
        elif d.get("from") == "column":
            v = str(r.get(d["col"], "")).strip()
            rec["direction"] = d.get("map", {}).get(v, v)
        elif d.get("from") == "compare":
            # case vs control expression columns (e.g. FPKM), for tables with no signed effect
            a, b = num(r.get(d["case"])), num(r.get(d["control"]))
            rec["direction"] = "" if a is None or b is None else ("up" if a > b else "down" if a < b else "")
        elif d.get("from") == "constant":
            rec["direction"] = d["value"]
        if d.get("invert"):
            # the table reports control relative to case (e.g. 'CTR vs. SCZ'), whatever the source of the sign
            if rec.get("direction") in ("up", "down"):
                rec["direction"] = {"up": "down", "down": "up"}[rec["direction"]]
            # and the effect is turned round with it, so its sign agrees with direction
            rec["effect"] = invert_effect(rec.get("effect"), d.get("center", 0))
        rec["symbol"], rec["hgnc_id"], rec["map_route"] = normalize_gene(rec, spec)
        rec["map_route"] += row_tags(rec, spec)
        out.append(rec)
    return out


# ---------------------------------------------------------------- tags (2026-09-30)
# Rows are tagged, not dropped: ';tag:control-probe' for array control probe sets (Affymetrix AFFX-, ERCC spike-ins),
# ';tag:reference-gene' for genes on the H4 list of qPCR reference / normalization genes (adopted as is; a trailing *
# there is a prefix match).
CONTROL_PROBE = re.compile(r"^(AFFX-|ERCC-)", re.I)
REFERENCE_GENE_FILE = os.path.join(HERE, "h4", "reference_genes", "reference_gene_list.tsv")
TAG_FILES = set()
_REFERENCE = None


def reference_genes():
    """(exact symbols, prefixes) from the H4 reference-gene list."""
    global _REFERENCE
    if _REFERENCE is None:
        exact, prefixes = set(), []
        with open(REFERENCE_GENE_FILE, encoding="utf-8") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                s = r["symbol"].strip()
                if s.endswith("*"):
                    prefixes.append(s[:-1])
                elif s:
                    exact.add(s)
        _REFERENCE = (exact, tuple(prefixes))
    return _REFERENCE


# routes whose symbol comes from the paper's printed name (symbol_to_hgnc, name_to_hgnc, the X09 name fallbacks), as
# against an identifier route (platform, Gemma platform, Ensembl, Entrez, accession, UniProt, UniGene) or a curator's
# override. LOC<GeneID> printed as the name counts as a printed name.
NAME_ROUTES = {"approved", "previous", "previous-dated", "alias", "loc", "resolver", "excel-date", "cleaned", "mirna", "mirna-arm"}


def printed_name_fallback(rec, spec):
    """True when the list's manifest platform holds the row's probe ID (platform_probes) but annotates it to no gene, or
    to a symbol HGNC lacks ('platform-unmapped->paper-name'), and the row's symbol came from the printed name. The
    2026-09-30 Step 2 audit's 'unmapped_now' rows (step2_symbol_audit_2026-09-30/probe_vs_platform_rows.tsv) whose
    symbol came from the paper, in lists whose manifest names a platform. A probe not on the platform at all is not
    tagged (the audit's 'not_on_platform': the platform cannot speak for it)."""
    route = (rec.get("map_route") or "").split(";")[0]
    if not rec.get("symbol") or not spec or not spec.get("platform"):
        return False
    if route.startswith("platform-unmapped->paper-name"):
        return True
    if re.split(r"->|:", route, maxsplit=1)[0] not in NAME_ROUTES:
        return False
    rid = (rec.get("id_reported") or "").strip()
    if is_missing(rid):
        return False
    probes = platform_probes(spec["platform"], spec.get("platform_id_column"), rid)
    table = platform(spec["platform"], spec.get("platform_id_column"))
    return bool(probes) and not any(table[p] for p in probes)


def row_tags(rec, spec=None):
    tags = []
    if any(CONTROL_PROBE.match(str(rec.get(k) or "").strip()) for k in ("id_reported", "gene_reported")):
        tags.append("control-probe")
    s = rec.get("symbol") or ""
    if s:
        exact, prefixes = reference_genes()
        if s in exact or s.startswith(prefixes):
            tags.append("reference-gene")
            TAG_FILES.add(REFERENCE_GENE_FILE)
    if printed_name_fallback(rec, spec):
        tags.append("printed-name-fallback")
    return ";tag:" + ",".join(tags) if tags else ""


# list-level role (manifest "role"): the study's own new data, or its analysis of data first
# published elsewhere (Paul, 2026-09-24: "keep 'first look at the data' separate"). Only first_look
# lists can represent a study in recurrence.py and h1.py.
ROLES = {"first_look": "the study's own new data", "reanalysis": "the study's analysis of data first published elsewhere"}

FIELDS = ["gene_reported", "id_reported", "symbol", "hgnc_id", "map_route", "direction", "effect", "p", "fdr", "group"]
# "group": cell type or region, for lists stacked from several sheets ("sheets") or a column naming it


def main(only=None):
    """only: list names (<pmid>__<list_id>) to re-extract; the other lists' files are left as they are and their
    summary.tsv rows are kept from the existing file."""
    global HGNC
    HGNC = load_hgnc()
    manifest = json.load(open(os.path.join(DATA, "manifest.json")))
    # batch manifests (manifests/*.json, same schema) let parallel workers avoid editing one file
    for extra in sorted(glob.glob(os.path.join(DATA, "manifests", "*.json"))):
        manifest["lists"] += json.load(open(extra))["lists"]
    os.makedirs(OUT, exist_ok=True)
    today = datetime.date.today().isoformat()
    summary = []
    kept = {}
    if only is not None:
        names = {f"{s['pmid']}__{s['list_id']}" for s in manifest["lists"]}
        if set(only) - names:
            raise SystemExit(f"not in the manifests: {sorted(set(only) - names)}")
        with open(os.path.join(HERE, "summary.tsv"), encoding="utf-8") as fh:
            kept = {r["list"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for spec in manifest["lists"]:
        name = f"{spec['pmid']}__{spec['list_id']}"
        if only is not None and name not in only:
            if name not in kept:
                raise SystemExit(f"{name}: no row in summary.tsv to keep; run extract.py on all lists")
            summary.append({k: v for k, v in kept[name].items() if k not in ("provenance", "role", "bank")})
            continue
        if spec.get("status"):
            summary.append({"list": name, "label": spec["label"], "status": spec["status"],
                            "n_rows": "", "n_genes": "", "n_up": "", "n_down": "", "expected_n": "",
                            "matches_paper": "", "unmapped": "", "universe": spec.get("universe", "")})
            continue
        FALLBACK_FILES.clear()
        TAG_FILES.clear()
        rows = build(spec)
        genes = {r["symbol"] for r in rows if r["symbol"]}
        exp = spec.get("expected_n")
        src = os.path.join(DATA, spec["file"])
        with open(os.path.join(OUT, name + ".tsv"), "w", encoding="utf-8") as fh:
            fh.write(f"# study: {spec['label']} (PMID {spec['pmid']})\n")
            fh.write(f"# list: {spec['list_id']} - {spec['contrast']}\n")
            prov = spec.get("provenance", "published")
            if prov != "published":
                # any non-published list (reproduced by us, reanalysed by Gemma, ...) says so in its header
                what = {"reproduced": "REPRODUCED by us from public data",
                        "gemma": "REANALYSED BY GEMMA (Gemma's own differential expression analysis)"}.get(prov, f"{prov.upper()}")
                fh.write(f"# provenance: {what}, not a list the authors published ({spec.get('reproduction', '')})\n")
            role = spec.get("role", "first_look")
            if role not in ROLES:
                raise SystemExit(f"{name}: role {role!r} is not one of {sorted(ROLES)}")
            fh.write(f"# role: {role} ({ROLES[role]})\n")
            if spec.get("role_note"):
                # dated evidence for the role, where it was corrected or is not obvious from the contrast
                fh.write(f"# role note: {spec['role_note']}\n")
            if spec.get("bank"):
                # one bank's part of a study that spans several banks (lineage node 'banks')
                fh.write(f"# bank: {spec['bank']}\n")
            fh.write(f"# threshold (paper): {spec['threshold']}\n")
            if spec.get("direction", {}).get("invert"):
                how = "reciprocal of" if spec["direction"].get("center", 0) == 1 else "sign reversed from"
                fh.write(f"# effect: {how} the source value; the source reports control relative to case, "
                         f"so effect and direction here are case relative to control\n")
            if spec.get("effect_note"):
                # what the effect column holds, where its name or sign would mislead
                fh.write(f"# effect: {spec['effect_note']}\n")
            fh.write(f"# source: {spec['file']}{' :: ' + spec['member'] if spec.get('member') else ''} {spec.get('sheet') or ', '.join(spec.get('sheets', []))}\n")
            fh.write(f"# source sha256: {sha256(src)}\n")
            fh.write(f"# script: extract.py via manifest.json, run {today}\n")
            fh.write(f"# hgnc: {os.path.basename(HGNC_FILE)} sha256 {sha256(HGNC_FILE)}\n")
            if FALLBACK_FILES:
                # files the X09-X11 fallbacks of normalize_gene read for this list
                fh.write("# fallback sources: " + "; ".join(
                    f"{os.path.relpath(f, os.path.dirname(HERE))} sha256 {fallback_sha(f)}" for f in sorted(FALLBACK_FILES)) + "\n")
            if TAG_FILES:
                # files the map_route tags (';tag:...') were read from
                fh.write("# tag sources: " + "; ".join(
                    f"{os.path.relpath(f, os.path.dirname(HERE))} sha256 {fallback_sha(f)}" for f in sorted(TAG_FILES)) + "\n")
            w = csv.DictWriter(fh, FIELDS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows:
                out = {k: r.get(k, "") for k in FIELDS}
                if is_missing(out["gene_reported"]):
                    out["gene_reported"] = ""
                if str(out["gene_reported"]).startswith("#"):
                    raise SystemExit(f"{name}: gene_reported {out['gene_reported']!r} starts with '#'; readers would skip the row")
                w.writerow(out)
        summary.append({
            # status_note: a dated history note (e.g. a list re-included after a data fix); readers split on ":"
            "list": name, "label": spec["label"], "status": "extracted" + (f": {spec['status_note']}" if spec.get("status_note") else ""),
            "n_rows": len(rows), "n_genes": len(genes),
            "n_up": sum(r.get("direction") == "up" for r in rows),
            "n_down": sum(r.get("direction") == "down" for r in rows),
            "expected_n": exp if exp is not None else "",
            "matches_paper": "" if exp is None else ("yes" if exp == len(rows) else "NO"),
            "unmapped": sum(1 for r in rows if not r["symbol"]),
            "universe": spec.get("universe", ""),
        })
    for s, spec in zip(summary, manifest["lists"]):
        s["provenance"] = spec.get("provenance", "published")
        s["role"] = spec.get("role", "first_look")
        s["bank"] = spec.get("bank", "")
    # atomic: several workers may run extract.py at once
    tmp = os.path.join(HERE, f".summary.{os.getpid()}.tsv")
    with open(tmp, "w", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(summary[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    os.replace(tmp, os.path.join(HERE, "summary.tsv"))
    for s in summary:
        if only is None or s["list"] in only:
            print("\t".join(str(s[k]) for k in ("list", "status", "n_rows", "n_genes", "n_up", "n_down", "expected_n", "matches_paper", "unmapped")))


if __name__ == "__main__":
    # python3 extract.py [--lists <name>,<name>,... | --lists-file <file, one name per line>]
    args, only = sys.argv[1:], None
    if args[:1] == ["--lists"] and len(args) == 2:
        only = set(filter(None, args[1].split(",")))
    elif args[:1] == ["--lists-file"] and len(args) == 2:
        only = {l.strip() for l in open(args[1], encoding="utf-8") if l.strip() and not l.startswith("#")}
    elif args:
        raise SystemExit("usage: extract.py [--lists a,b | --lists-file f]")
    main(only)
