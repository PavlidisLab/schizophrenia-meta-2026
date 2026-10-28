#!/usr/bin/env python3
"""Tidy per-gene effects table of the 383 recurrent genes across the 36 studies that have a full per-gene effect table.

Read-only on every input; no network. Writes only into the folder this script sits in:
  effects_383_by_study.tsv.gz   one row per study x group x gene (of the 383) present in the chosen table
  effects_383_study_index.tsv   one row per study

Table reading, HGNC symbol mapping and the choice of source per study are copied from
build_inventory.py (read_es, read_stage, read_raw, id_to_symbol;
the inventory script itself cannot be imported because it rewrites its own outputs at import). The chosen source per study is read
from that folder's inventory_80_studies.tsv (source_kind, source_of_full_effects).
Run: python3 build_effects_383.py
"""
import csv, gzip, hashlib, json, math, os, glob, datetime, collections
import pandas as pd

HOME = os.path.expanduser("~")
AM = f"{HOME}/analysis-delivery"
REPO = f"{HOME}/Dev/schizophrenia-transcriptome-replication"
OUT = os.path.dirname(os.path.abspath(__file__))
INV_DIR = f"{AM}/fig2c_effects_inventory_2026-10-19"
TODAY = datetime.date.today().isoformat()

P_INV = f"{INV_DIR}/inventory_80_studies.tsv"
P_CAND = f"{INV_DIR}/candidate_sources.tsv"
P_FACTS = f"{INV_DIR}/facts.json"
P_INV_SCRIPT = f"{INV_DIR}/build_inventory.py"
P_BASE = f"{AM}/recurrence_on_80_studies/base80_studies.tsv"
P_GENES = f"{AM}/recurrence_on_80_studies/out/recurrence_genes.tsv"
P_UNITS = f"{AM}/recurrence_on_80_studies/out/recurrence_units.tsv"
P_PRIMARY = f"{AM}/recurrence_on_80_studies/ref_B80_inputs/primary_lists.json"
ES = f"{REPO}/analyses/effect_sizes/data"
P_SOURCES = f"{ES}/sources.tsv"
STAGE = f"{AM}/derived_all_2026-10-14/work/stage"
P_STAGE_SUM = f"{AM}/derived_all_2026-10-14/derived_all_summary.tsv"
P_HGNC = f"{REPO}/data/reference/hgnc_complete_set.txt"
RAW = f"{REPO}/data/raw"

INPUTS = {}
FILE_SHA = {}


def shortpath(path):
    """Path relative to a named root: starts with schizophrenia-transcriptome-replication/ or sz-transcriptome/."""
    path = os.path.normpath(path)
    for root, name in ((REPO, "schizophrenia-transcriptome-replication"), (AM, "sz-transcriptome")):
        if path.startswith(root + "/"):
            return name + "/" + path[len(root) + 1:]
    raise ValueError(path)


def sha_full(path):
    if path not in FILE_SHA:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        FILE_SHA[path] = h.hexdigest()
    return FILE_SHA[path]


def sha(path):
    INPUTS[shortpath(path)] = sha_full(path)[:12]
    return sha_full(path)


def tsv(path):
    return [l for l in open(path, encoding="utf-8") if not l.startswith("#")]


for p in (P_INV, P_CAND, P_FACTS, P_INV_SCRIPT, P_BASE, P_GENES, P_UNITS, P_PRIMARY, P_SOURCES, P_STAGE_SUM, P_HGNC):
    sha(p)

# ------------------------------------------------------------------ the 383 (as in the inventory)
rows = list(csv.DictReader(tsv(P_GENES), delimiter="\t"))
pb = [r for r in rows if r["scenario"] == "published, brain"]
S383 = {r["gene"] for r in pb if float(r["q_bh"]) < 0.1}
UNIV = {r["gene"] for r in pb}
assert len(S383) == 383, len(S383)

APPROVED, ENSG, ENTREZ = set(), {}, {}
for r in csv.DictReader(open(P_HGNC, encoding="utf-8"), delimiter="\t"):
    if r["status"] != "Approved":
        continue
    APPROVED.add(r["symbol"])
    if r["ensembl_gene_id"]:
        ENSG[r["ensembl_gene_id"]] = r["symbol"]
    if r["entrez_id"]:
        ENTREZ[r["entrez_id"]] = r["symbol"]


def to_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def id_to_symbol(g):
    """Copied from build_inventory.py. Return (symbol_or_None, id_kind). HGNC lookup only for Ensembl and Entrez ids."""
    if g in APPROVED or g in UNIV:
        return g, "symbol"
    base_id = g.split(".")[0]
    if base_id in ENSG:
        return ENSG[base_id], "ensembl"
    if g in ENTREZ:
        return ENTREZ[g], "entrez"
    return None, "other"


# ------------------------------------------------------------------ readers (record = group, gene, effect, se, t, p, q)
def read_es(path):
    recs = []
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            recs.append((r["group"], r["gene"], to_float(r["log2fc"]), to_float(r["se"]), to_float(r["t"]), to_float(r["p"]), None))
    return recs


def read_stage_file(unit, f):
    recs = []
    with gzip.open(f, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        has_fdr = "fdr" in rd.fieldnames
        for r in rd:
            recs.append((unit, r["gene"], to_float(r["logfc"]), None, to_float(r["t"]), to_float(r["p"]),
                         to_float(r["fdr"]) if has_fdr else None))
    return recs


def stage_path(unit):
    f = f"{STAGE}/{unit}.tsv.gz"
    if not os.path.exists(f):  # earlier-job fits kept under prev_src (same columns plus fdr)
        f = f"{STAGE}/../prev_src/stage2/{unit}_L1.tsv.gz"
    return os.path.normpath(f)


def read_raw(pmid):
    if pmid == "41942462":
        f = f"{RAW}/41942462/epmc/41467_2026_71281_MOESM5_ESM.xlsx"
        d = pd.read_excel(f, sheet_name="A", header=2)
        recs = [("dopaminergic neurons (Supplementary Data 3 sheet A)", str(r.iloc[0]), to_float(r["logFC"]), None, to_float(r["t"]),
                 to_float(r["P.Value"]), to_float(r["adj.P.Val"])) for _, r in d.iterrows()]
        return f, f + "::A", recs
    if pmid == "42461238":
        f = f"{RAW}/42461238/figshare/iwbp_a_2696844_sm2740.xlsx"
        d = pd.read_excel(f, sheet_name="Supplementary Table 2", header=3)
        recs = []
        for _, r in d.iterrows():
            fc = to_float(r["Fold Sz / Control"])
            lf = math.log2(fc) if fc and fc > 0 else None
            sym = r["Gene Symbol"]
            recs.append(("all", sym.strip() if isinstance(sym, str) else "", lf, None, None, to_float(r["p"]), None))
        return f, f + "::Supplementary Table 2", recs
    if pmid == "33184259":
        f = f"{RAW}/33184259/figshare/NOIseq_DEG_schizophrenia_vs_control.xls"
        recs = []
        with open(f, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                recs.append(("all", r["GeneID"], to_float(r["log2FoldChange(Group2/Group1)"]), None, None, None, None))
        return f, f, recs
    raise KeyError(pmid)


# ------------------------------------------------------------------ method / scale text
STAGE_META = {r["unit"]: r for r in csv.DictReader(open(f"{STAGE}/meta.tsv"), delimiter="\t")}
sha(f"{STAGE}/meta.tsv")
SRC = {r["study"]: r for r in csv.DictReader(tsv(P_SOURCES), delimiter="\t")}

K_TEXT = {
    "K1_gemma_L1_adopted": ("limma moderated t, project L1 fit (SZ + age + sex where recorded); table adopted from fairness_sample_level/results/b_de (derived_all 01_stage_fits.R)",
                            "limma logFC (SZ vs control)"),
    "K2_libd_limma_trend": ("limma-trend (eBayes trend=TRUE) on logCPM, SZ + age + sex where recorded (derived_all 01_stage_fits.R)",
                            "limma logFC on logCPM (SZ vs control)"),
    "K3_pseudobulk_class": ("pseudobulk per cell class: filterByExpr, TMM, voom, limma with SZ + age + sex (+ bank for PsychAD) + log(nuclei) (derived_all 01_stage_fits.R)",
                            "limma-voom logFC on logCPM (SZ vs control)"),
}
PREV_TEXT = ("limma, SZ + age + sex where recorded; earlier derived-list job table adopted unchanged by derived_all (02_build.py label 'E'); "
             "method detail not restated there",
             "limma logFC (SZ vs control)")

PUB = {  # pmid -> (method, effect_scale), from effect_sizes/data/sources.tsv full_conversion and list headers
    "32203495": ("published authors' table; limma (columns SZ_logFC, SZ_t per sources.tsv)", "limma logFC (SZ_logFC)"),
    "34440415": ("published authors' table; DESeq2 (log2FoldChange; no p or SE in Table S2)", "log2 fold change (DESeq2 log2FoldChange)"),
    "36223459": ("published authors' table; DESeq2 per subtype on per-sample pseudobulk (list header)", "log2 fold change (DESeq2 log2FoldChange)"),
    "36319771": ("published authors' table; limma (logFC, t; Supp Data 8)", "limma logFC"),
    "38781388": ("published authors' table; limma pseudobulk per cohort, fixed-effects meta-analysis (Meta_logFC, Meta_se; list header)", "meta-analysed logFC (Meta_logFC); log2 scale not stated"),
    "39526318": ("published authors' table; DESeq2 (log2FoldChange, lfcSE)", "log2 fold change (DESeq2 log2FoldChange)"),
    "41174894": ("published authors' table; limma (Hb_logFC, Hb_t; ETable 7)", "limma logFC (Hb_logFC)"),
    "42778699": ("published authors' table; dreamlet per subclass, meta-analysed across brain banks (list header)", "dreamlet estimate (log2 CPM scale)"),
}
RAWTXT = {
    "41942462": ("published authors' table; limma (columns logFC, t, P.Value, adj.P.Val of Supplementary Data 3 sheet A); only sheet A (dopaminergic neurons) read, sheet B (non-neurons) not read",
                 "limma logFC; scale not stated"),
    "42461238": ("published authors' table; analysis in JMP Genomics 9.0 (test not stated in the table title); probe-set level",
                 "log2 of 'Fold Sz / Control' (fold ratio converted to log2 by this script)"),
    "33184259": ("published authors' table; NOISeq via BGI Genomics (the table holds a NOISeq probability, not a P value)",
                 "log2FoldChange(Group2/Group1) as in the table; group order not checked against SZ/control"),
}
GEMMA_TEXT = ("Gemma result set (per-probe results collapsed to the first single-gene probe per gene by probe id; sources.tsv full_conversion)",
              "Gemma log2 fold change")

# ------------------------------------------------------------------ list structure
primary = json.load(open(P_PRIMARY))["studies"]
sha(P_PRIMARY)
base = {r["pmid"]: r for r in csv.DictReader(tsv(P_BASE), delimiter="\t")}
inv = [r for r in csv.DictReader(tsv(P_INV), delimiter="\t") if r["source_kind"] != "none"]
assert len(inv) == 36, len(inv)
stage_units = collections.defaultdict(list)
for r in csv.DictReader(tsv(P_STAGE_SUM), delimiter="\t"):
    stage_units[r["pmid"]].append(r["unit"])

out_rows = []   # list of dicts
index = []
stats = {}      # label -> dict(n_rows, n_genes, n_p, dups, ...)
HDR_UNDET = []
DUP_NOTES = []
MISMATCH = []

for r in inv:
    pmid, label, kind = r["pmid"], r["study_label"], r["source_kind"]
    srcs = []   # (group, path_for_hash, shown_path, recs, method, scale)
    if kind in ("published", "gemma"):
        f = f"{ES}/full/{pmid}__{kind}.tsv.gz"
        sha(f)
        recs = read_es(f)
        method, scale = (PUB[pmid] if kind == "published" else GEMMA_TEXT)
        srcs.append((None, f, f, recs, method, scale))
    elif kind == "stage":
        for u in stage_units[pmid]:
            f = stage_path(u)
            sha(f)
            recs = read_stage_file(u, f)
            if u in STAGE_META and "prev_src" not in f:
                method, scale = K_TEXT[STAGE_META[u]["kind"]]
            else:
                method, scale = PREV_TEXT
            srcs.append((u, f, f, recs, method, scale))
    elif kind == "raw_supplement":
        f, shown, recs = read_raw(pmid)
        sha(f)
        method, scale = RAWTXT[pmid]
        srcs.append((None, f, shown, recs, method, scale))
    else:
        raise ValueError(kind)
    # sanity: inventory's chosen file is the one used
    inv_path = r["source_of_full_effects"].replace("~", HOME)
    for part in inv_path.split("; "):
        pass
    if kind != "stage":
        assert srcs[0][2] == inv_path.split("::")[0] + ("::" + inv_path.split("::")[1] if "::" in inv_path else ""), (label, srcs[0][2], inv_path)
    else:
        assert [os.path.normpath(s[2]) for s in srcs] == [os.path.normpath(x) for x in inv_path.split("; ")], label

    study_rows = []
    group_names = []
    for gunit, path, shown, recs, method, scale in srcs:
        by = collections.OrderedDict()   # (group, symbol) -> list of records
        for (grp, gene, eff, se, t, p, q) in recs:
            if gene in ("", "nan", None):
                continue
            sym, _k = id_to_symbol(gene)
            if sym is None or sym not in S383:
                continue
            by.setdefault((grp, sym), []).append((eff, se, t, p, q, gene))
        for (grp, sym), cand in by.items():
            if len(cand) > 1:
                if kind == "stage":      # derived-list build (derived_all 02_build.py) kept the first row per gene
                    pick = cand[0]
                else:                    # smallest P, then largest |effect|, then first
                    pick = sorted(enumerate(cand), key=lambda ic: (ic[1][3] is None, ic[1][3] if ic[1][3] is not None else 0,
                                                                   -abs(ic[1][0]) if ic[1][0] is not None else 0, ic[0]))[0][1]
                DUP_NOTES.append((label, grp, sym, len(cand)))
            else:
                pick = cand[0]
            eff, se, t, p, q, gene = pick
            study_rows.append(dict(group_raw=grp, gene=sym, effect=eff, t=t, se=se, p=p, q=q, shown=shown, path=path, method=method, scale=scale))
    groups = list(dict.fromkeys(x["group_raw"] for x in study_rows))
    # groups present in the source (including those with none of the 383): use srcs
    src_groups = []
    for gunit, path, shown, recs, method, scale in srcs:
        src_groups += list(dict.fromkeys(x[0] for x in recs))
    src_groups = list(dict.fromkeys(src_groups))
    n_grp = len(src_groups)
    if n_grp == 1 and pmid != "41942462":
        gnames = {src_groups[0]: "all"}
    else:
        gnames = {g: g for g in src_groups}
    # fig2c group
    ent = primary.get(pmid, {})
    n_units_in_list = int(base[pmid]["n_units_in_list"])
    pooled = n_grp > 1
    reason = ""
    if pooled:
        reason = (f"study list is the union over its {n_units_in_list} list unit(s) ({', '.join(ent.get('primary', []))}); "
                  f"table has {n_grp} groups; no single group produced the list")
        HDR_UNDET.append((label, reason))
    for x in study_rows:
        x["is_fig2c"] = "TRUE"
    genes = {x["gene"] for x in study_rows}
    stats[label] = dict(pmid=pmid, n_rows=len(study_rows), n_genes=len(genes), n_p=sum(x["p"] is not None for x in study_rows),
                        n_groups=n_grp, kind=kind)
    exp = int(r["n_of_383_with_effect"])
    if len(genes) != exp:
        MISMATCH.append((label, exp, len(genes)))
    for x in study_rows:
        out_rows.append(dict(study_label=label, pmid=pmid, gene=x["gene"], effect=x["effect"], t=x["t"], se=x["se"], p=x["p"], q=x["q"],
                             group_or_unit=gnames[x["group_raw"]], is_fig2c_group=x["is_fig2c"], method=x["method"], effect_scale=x["scale"],
                             source_kind=kind, source_file=shortpath(x["shown"].split("::")[0]) + ("::" + x["shown"].split("::")[1] if "::" in x["shown"] else ""),
                             source_sha256=sha_full(x["path"])))
    # index row
    gl = list(gnames.values())
    meths = list(dict.fromkeys(s[4] for s in srcs))
    scales = list(dict.fromkeys(s[5] for s in srcs))
    files = list(dict.fromkeys(shortpath(s[2].split("::")[0]) + ("::" + s[2].split("::")[1] if "::" in s[2] else "") for s in srcs))
    shas = list(dict.fromkeys(sha_full(s[1]) for s in srcs))
    note = [r["access_class"]]
    if r["has_list"] == "no":
        note.append("study list is empty at FDR 5% (dropped from the 72 lists); effects still held")
    if kind == "stage" and not r["list_origin"].startswith("derived"):
        note.append("table is our uniform reanalysis of the dataset; the study's list is the authors' own")
    if kind == "gemma":
        note.append("authors' list pools layers 3 and 5; Gemma contrast match per fairness_lists/a1_datasets.tsv: see inventory note")
    if kind == "raw_supplement":
        note.append("parsed by the inventory/this script from the supplementary file; no gene-level file of this table in the repository tables")
    if pmid == "42461238":
        note.append("probe-set level table; several probe sets per gene resolved by smallest P")
    index.append(dict(study_label=label, pmid=pmid, n_groups=n_grp, groups="; ".join(gl), fig2c_group=("all" if not pooled else "undetermined (all groups TRUE)"),
                      n_genes_fig2c_group=len(genes), source_kind=kind, method=" | ".join(meths), effect_scale=" | ".join(scales),
                      source_file="; ".join(files), source_sha256="; ".join(shas), access_note="; ".join(note)))

# ------------------------------------------------------------------ checks
keys = collections.Counter((x["study_label"], x["group_or_unit"], x["gene"]) for x in out_rows)
dups = [k for k, v in keys.items() if v > 1]
assert not dups, dups[:5]
assert all(x["study_label"] for x in out_rows)


def fmt(v):
    if v is None:
        return ""
    return repr(float(v)) if not isinstance(v, str) else v


COLS = ["study_label", "pmid", "gene", "effect", "t", "se", "p", "q", "group_or_unit", "is_fig2c_group", "method", "effect_scale",
        "source_kind", "source_file", "source_sha256"]
out_rows.sort(key=lambda x: (x["study_label"], x["group_or_unit"], x["gene"]))

H = []
H.append(f"# effects_383_by_study.tsv.gz: per-gene SZ-vs-control effects of the 383 recurrent genes in the 36 studies with a full per-gene effect table")
H.append(f"# script: build_effects_383.py (this folder), run {TODAY}; read-only on inputs, no network")
H.append("# inputs (sha256 first 12 hex; per-row source_sha256 holds the full 64 hex of each study's source table):")
for k, v in sorted(INPUTS.items()):
    H.append(f"#   {k} {v}")
H.append("# definitions:")
H.append("#   the 383 genes = scenario 'published, brain', study tier, q_bh < 0.1 in recurrence_on_80_studies/out/recurrence_genes.tsv (asserted 383)")
H.append("#   source per study = the inventory's chosen source (inventory_80_studies.tsv): authors' published table > our uniform reanalysis (stage) > Gemma > reproduction > raw supplement; derived-list studies use our own fit")
H.append("#   symbols: exact HGNC-approved or universe symbol; Ensembl (version stripped) and Entrez ids mapped through the HGNC file (Maas 2020 Entrez, Singh 2026 Ensembl); copied from build_inventory.py id_to_symbol; never joined by row position")
H.append("#   gene = HGNC symbol after mapping. effect = the table's stored effect (log fold change or the scale in effect_scale); t, se, p, q = stored values copied, nothing computed from t; blank = the source table has none")
H.append("#   se: for rows from effect_sizes/data/full the se column is as stored there (harmonize.py/ash_fit.R states se given > log2fc/t > |log2fc|/qt(1-p/2,df), so some may be derived upstream); stage tables carry no se (se = logFC/t)")
H.append("#   q: stored only where the source has it (stage tables read from derived_all/work/prev_src/stage2 carry an 'fdr' column; Singh 2026 sheet A adj.P.Val); otherwise blank. Maas 2020 holds a NOISeq probability, not a P: p blank")
H.append("#   effect is stored for every row used here, so no row has effect blank because only t and se were stored")
H.append("#   Dean 2026: effect = log2 of the table's 'Fold Sz / Control' ratio (converted here); the fold ratio itself is not stored in the output")
H.append("#   duplicate rule: one row per study x group x gene. Tables in effect_sizes/data/full are already one row per gene per group. Where several source rows map to one gene "
         "(Dean 2026 probe sets; Ensembl or Entrez ids sharing a symbol), stage tables keep the first row in the file (the rule of derived_all 02_build.py, which built the derived lists with drop_duplicates('gene')); "
         "other sources keep the row with the smallest P, then the largest |effect|, then the first. The inventory had no row-choice rule (it counted distinct genes).")
H.append(f"#   rows resolved by the duplicate rule: {len(DUP_NOTES)}" + ("; " + "; ".join(f"{a} {b} {c} x{d}" for a, b, c, d in DUP_NOTES[:40]) if DUP_NOTES else ""))
H.append("#   group_or_unit: region, layer, cell type or dataset-region as named in the source (stage tables: the derived_all unit name); 'all' where the study's table is one comparison. Singh 2026 keeps its sheet name because the file also has a non-neuron sheet that was not read")
H.append("#   is_fig2c_group: TRUE for the group whose result produced the study's list in the recurrence analysis. In every study with more than one group the study's list is the union over groups "
         "(primary_lists.json / base80_studies.tsv n_units_in_list), so no single group can be named and TRUE is set on every group of that study (listed below). Single-group studies: TRUE, determined")
H.append("#   method / effect_scale: from effect_sizes/data/sources.tsv full_conversion, the list headers in data/lists and derived_all scripts (01_stage_fits.R, 02_build.py); 'not stated' not needed except where noted in the text")
H.append("#   is_fig2c_group undetermined (all groups TRUE): " + "; ".join(f"{a} [{b}]" for a, b in HDR_UNDET))
H.append(f"#   rows: {len(out_rows)}; studies: {len(stats)}")
H.append("#   per study: label | pmid | rows | genes (distinct of 383) | groups | rows with p")
for lab, s in stats.items():
    H.append(f"#   {lab} | {s['pmid']} | {s['n_rows']} | {s['n_genes']} | {s['n_groups']} | {s['n_p']}/{s['n_rows']} = {s['n_p'] / max(s['n_rows'], 1):.3f}")
H.append(f"#   check 1 (distinct genes per study vs inventory n_of_383_with_effect): " +
         ("all 36 equal" if not MISMATCH else "MISMATCH " + "; ".join(f"{a}: inventory {b}, here {c}" for a, b, c in MISMATCH)))
H.append("#   check 2 (no duplicate study_label, group_or_unit, gene): passed (asserted)")
H.append("#   check 3 (fraction of rows with p): per study above")
with gzip.open(f"{OUT}/effects_383_by_study.tsv.gz", "wt", encoding="utf-8", newline="") as fh:
    fh.write("\n".join(H) + "\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(COLS)
    for x in out_rows:
        w.writerow([fmt(x[c]) if c in ("effect", "t", "se", "p", "q") else x[c] for c in COLS])

IH = ["# effects_383_study_index.tsv: one row per study of effects_383_by_study.tsv.gz",
      f"# script: build_effects_383.py (this folder), run {TODAY}",
      "# inputs (sha256 first 12 hex): " + "; ".join(f"{k} {v}" for k, v in sorted(INPUTS.items())),
      "# n_groups = groups/units in the study's chosen table; groups = their names; fig2c_group = 'all' for a one-comparison study, 'undetermined (all groups TRUE)' where the study's list pools groups",
      "# n_genes_fig2c_group = distinct of the 383 genes in the rows with is_fig2c_group TRUE; source_file/source_sha256 list one entry per table (several for multi-unit stage tables)"]
ICOLS = ["study_label", "pmid", "n_groups", "groups", "fig2c_group", "n_genes_fig2c_group", "source_kind", "method", "effect_scale", "source_file", "source_sha256", "access_note"]
with open(f"{OUT}/effects_383_study_index.tsv", "w", encoding="utf-8", newline="") as fh:
    fh.write("\n".join(IH) + "\n")
    w = csv.DictWriter(fh, fieldnames=ICOLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(index)

print("rows", len(out_rows), "studies", len(stats), "dups resolved", len(DUP_NOTES), "mismatch", MISMATCH)
print("undetermined", [a for a, _ in HDR_UNDET])
for lab, s in stats.items():
    print(lab, s["n_rows"], s["n_genes"], s["n_groups"], f"{s['n_p'] / max(s['n_rows'], 1):.2f}")
