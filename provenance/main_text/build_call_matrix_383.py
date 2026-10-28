#!/usr/bin/env python3
"""Gene x study call matrix for the 383 recurrent genes (and the 80-study version), same call logic as fig2_b in
<analysis repo>/analyses/figure_data_exports/export_fig02_fig07.py.

Read-only on every input; no network. Writes only into the folder this script sits in:
  fig02_b_gene_by_study_long_383.tsv   383 genes x 72 studies that have a list
  fig02_b_gene_by_study_long_80.tsv    383 genes x 80 studies (the 8 empty-list studies appended)
and the scratch tree work_call_matrix/ (symlinks and one copy of recurrence_fork.py; see make_tree()).

The repo's own loading path is used: recurrence.scenario_selection('published, brain') runs inside work_call_matrix/repo, a copy of the
analysis tree whose primary_lists.json and universes.json are the 80-study files, whose lists/ also holds the derived lists, and whose
universe_sets/ also holds the rebuilt derived__<pmid>.txt sets. The environment is the one the 80-study run used
(RECURRENCE_NK=membership; RECURRENCE_INCLUDE_FILE = the 80 pmids of pipeline_inputs/primary_pmids.txt).
Validation (mandatory): the first 100 genes x 72 studies must equal the manuscript table fig02_b_gene_by_study_long.tsv in all five columns.
Run: python3 build_call_matrix_383.py
"""
import csv, hashlib, os, re, sys, json, collections, datetime

HOME = os.path.expanduser("~")
OUT = os.path.dirname(os.path.abspath(__file__))
AM = f"{HOME}/analysis-delivery"
SREP = f"{HOME}/Dev/srep-mauney-integrated"
MS = f"{HOME}/Documents/manuscripts/schizophrenia-transcriptome-reproducibility-manuscript"
RUN = f"{AM}/recurrence_on_80_studies"
B80 = f"{RUN}/ref_B80_inputs"
TREE = f"{OUT}/work_call_matrix/repo"
SCENARIO = "published, brain"

P_GENES = f"{RUN}/out/recurrence_genes.tsv"
P_UNITS = f"{RUN}/out/recurrence_units.tsv"
P_PRIM = f"{B80}/primary_lists.json"
P_UNIV = f"{B80}/universes.json"
P_PMIDS = f"{B80}/pipeline_inputs/primary_pmids.txt"
P_FORK = f"{B80}/recurrence_fork.py"
P_SETS = f"{RUN}/out/recurrent_sets.tsv"
P_BASE = f"{RUN}/base80_studies.tsv"
P_LINEAGE = f"{AM}/citation_fates_of_primary_studies/inputs/lineage_4f91.json"
P_EXIST = f"{MS}/figures/data/fig02/fig02_b_gene_by_study_long.tsv"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


def rel(p):
    for base, tag in ((SREP, "<analysis repo>/"), (AM, "~/analysis-delivery/"), (MS, "<manuscript repo>/"), (OUT, "")):
        if p.startswith(base + "/"):
            return tag + p[len(base) + 1:]
    raise SystemExit("path outside the allowed roots: " + p)


def rd(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))


def make_tree():
    """The scratch tree is built by the shell steps recorded in this function's docstring; here we only check it is in place."""
    for p in (f"{TREE}/analyses/recurrence.py", f"{TREE}/analyses/primary_lists.json", f"{TREE}/analyses/universes.json", f"{TREE}/corpus/lineage.json"):
        if not os.path.exists(p):
            raise SystemExit("scratch tree missing: " + p)
    if sha(f"{TREE}/analyses/recurrence.py") != sha(P_FORK):
        raise SystemExit("scratch recurrence.py is not recurrence_fork.py")


def main():
    make_tree()
    os.environ["RECURRENCE_NK"] = "membership"
    os.environ["RECURRENCE_INCLUDE_FILE"] = P_PMIDS
    sys.dont_write_bytecode = True
    sys.path.insert(0, f"{TREE}/analyses")
    import recurrence as R

    sel, M, lists, units = R.scenario_selection(SCENARIO)
    assert len(sel) == 72, len(sel)

    # the 383: study-tier set study_brain of the 80-study run
    sets = rd(P_SETS)
    g383 = {r["gene"] for r in sets if r["set"] == "study_brain"}
    assert len(g383) == 383, len(g383)

    grows = [r for r in rd(P_GENES) if r["scenario"] == SCENARIO]
    grows.sort(key=lambda r: (-int(r["n_lists_observed"]), r["gene"]))
    allrank = {r["gene"]: i for i, r in enumerate(grows, 1)}
    nobs = {r["gene"]: int(r["n_lists_observed"]) for r in grows}
    genes = [r["gene"] for r in grows if r["gene"] in g383]
    assert len(genes) == 383
    rank_in = {g: i for i, g in enumerate(genes, 1)}
    ranks = [allrank[g] for g in genes]
    contiguous = ranks == list(range(1, 384))

    # study order of the 72: the manuscript table's column order (cohort cluster, year, id from fig2_b)
    ex_rows = rd(P_EXIST)
    order72 = []
    for r in ex_rows:
        if r["study_id"] not in order72:
            order72.append(r["study_id"])
    assert len(order72) == 72 and set(order72) == set(sel), "study set differs from the manuscript table"

    # the 8 empty-list studies
    base = rd(P_BASE)
    empty = [r for r in base if r["empty_list_dropped"].strip().lower() in ("yes", "true", "1")]
    assert len(empty) == 8, len(empty)
    pl = json.load(open(P_PRIM))["studies"]
    uni = json.load(open(P_UNIV))["lists"]
    order8, M8, label = [], {}, {}
    for r in empty:
        k = r["pmid"]
        assert k not in sel
        names = pl[k]["primary"]
        u = [uni.get(n) for n in names]
        label[k] = r["label"]
        order8.append(k)
        if any(x is None or not x.get("gene_set") for x in u):
            M8[k] = None
            continue
        sets_k = {frozenset(R.read_set(x["gene_set"])) for x in u}
        assert len(sets_k) == 1, k
        M8[k] = set(next(iter(sets_k)))
        d = {}
        for n in names:
            d.update(R.read_list(n))
        assert not d, f"{k} is not empty"
    order80 = order72 + order8

    # base table labels for the header
    lab = {}
    for r in base:
        lab[r["pmid"]] = r["label"]

    outside = []

    def call_of(g, k):
        if k in sel:
            d = lists.get(k, {}).get(g)
            if d is not None and g not in M.get(k, ()):
                outside.append((g, k))  # on the list but outside its tested universe: n_lists_observed does not count it; since 2026-10-19 the cell is called not_tested (fig2_b called it by its direction)
                return "not_tested"
            if d in ("up", "down"):
                return d
            if d == "both":
                return "mixed"
            if d is not None:
                return "hit_no_direction"
            return "tested_not_hit" if g in M.get(k, ()) else "not_tested"
        m = M8[k]
        return "not_tested" if m is None else ("tested_not_hit" if g in m else "not_tested")

    rows72, rows80, mism = [], [], []
    for g in genes:
        hits = 0
        for k in order72:
            c = call_of(g, k)
            hits += c in ("up", "down", "mixed", "hit_no_direction")
            rows72.append(dict(gene=g, study_id=k, call=c, gene_rank=allrank[g], n_lists_observed=nobs[g], rank_within_383=rank_in[g]))
        if hits != nobs[g]:
            mism.append((g, hits, nobs[g]))
        rows80.extend(rows72[-72:])
        for k in order8:
            rows80.append(dict(gene=g, study_id=k, call=call_of(g, k), gene_rank=allrank[g], n_lists_observed=nobs[g], rank_within_383=rank_in[g]))
    out_set = sorted(set(outside))
    assert not mism, mism
    assert out_set == [("CARTPT", "16236523"), ("NPY", "16236523"), ("NR4A3", "16236523")] or [(g, str(k)) for g, k in out_set] == [("CARTPT", "16236523"), ("NPY", "16236523"), ("NR4A3", "16236523")], out_set
    print("hits != n_lists_observed:", mism, "; list genes outside their tested universe (called as fig2_b does, a hit):", out_set)
    assert all(r["call"] in ("tested_not_hit", "not_tested") for r in rows80 if r["study_id"] in order8)

    # direction counts against recurrence_genes.tsv (up/down/both/other_direction), an independent check
    gr = {r["gene"]: r for r in grows}
    dmis = []
    for g in genes:
        c = collections.Counter(r["call"] for r in rows72[(rank_in[g] - 1) * 72: rank_in[g] * 72])
        exp = (int(gr[g]["up"]), int(gr[g]["down"]), int(gr[g]["both"]), int(gr[g]["other_direction"]))
        if (c["up"], c["down"], c["mixed"], c["hit_no_direction"]) != exp:
            dmis.append((g, (c["up"], c["down"], c["mixed"], c["hit_no_direction"]), exp))

    # VALIDATION against the manuscript table: all 7,200 rows, five columns, in order
    keyed = {(r["gene"], r["study_id"]): r for r in rows72}
    bad = 0
    for r in ex_rows:
        n = keyed.get((r["gene"], r["study_id"]))
        if n is None or any(str(n[c]) != r[c] for c in ("call", "gene_rank", "n_lists_observed")):
            bad += 1
    seq_ok = all((r["gene"], r["study_id"]) == (rows72[i]["gene"], rows72[i]["study_id"]) for i, r in enumerate(ex_rows))
    print(f"validation: existing rows {len(ex_rows)}; mismatching rows {bad}; same row order {seq_ok}")
    print(f"hits == n_lists_observed for {len(genes) - len(mism)} of {len(genes)} genes (exceptions above); direction-count mismatches vs recurrence_genes.tsv: {len(dmis)}")
    print(f"gene_rank over the 383 contiguous 1..383: {contiguous}; max gene_rank {max(ranks)}; rows 383-file {len(rows72)}, 80-file {len(rows80)}")
    bad_cells = [(r["gene"], r["study_id"]) for r in ex_rows if keyed[(r["gene"], r["study_id"])]["call"] != r["call"]]
    assert bad == 1 and [(g, str(k)) for g, k in bad_cells] == [("CARTPT", "16236523")] and seq_ok, (bad, bad_cells)
    assert not dmis, dmis[:5]

    inputs = [P_GENES, P_UNITS, P_PRIM, P_UNIV, P_SETS, P_PMIDS, P_BASE, P_LINEAGE, P_FORK]
    callset = ("call: up/down/mixed = direction in the study's primary list (mixed = both within the study's union); tested_not_hit = gene in the list's tested "
               "universe M but not a hit; not_tested = outside it; hit_no_direction = hit with no stated direction; a gene on a study's list but outside its tested universe M is called not_tested (fig2_b, which made the current manuscript table, called these cells by their direction)")
    gap_txt = ("gene_rank is contiguous 1..383 (the 383 are the 383 top genes by n_lists_observed)" if contiguous else
               f"gene_rank has gaps (max {max(ranks)}): the 383 are the q_bh < 0.1 genes, not the top 383 by count; ranks of the 383 run 1..{max(ranks)}")

    exc = "; ".join(f"{g} in study {k}" for g, k in out_set)

    def write(name, rows, studies, extra):
        today = datetime.date.today().isoformat()
        hdr = [f"# figure: fig02 panel b, all-recurrent-genes version ({name})",
               f"# script: build_call_matrix_383.py (this folder), run {today}; same call logic as fig2_b in <analysis repo>/analyses/figure_data_exports/export_fig02_fig07.py",
               f"# population: {len(genes)} genes x {len(studies)} studies; scenario '{SCENARIO}', study tier; genes = set study_brain of recurrent_sets.tsv (q_bh < 0.1; 383); "
               f"studies with a list = the 72 representatives of recurrence.scenario_selection run on the 80-study inputs (RECURRENCE_NK=membership, RECURRENCE_EMPTY_LISTS=drop)" + extra,
               ]
        hdr += [f"# input: {rel(p)} sha256 {sha(p)}" for p in inputs]
        hdr += [f"# reference (validation only): {P_EXIST.replace(MS, '<manuscript repo>')} sha256 {sha(P_EXIST)}; the first 100 genes x 72 studies reproduce its 7,200 rows (gene, study_id, call, gene_rank, n_lists_observed) with 1 intended difference: CARTPT in study 16236523 is not_tested here ('up' in the reference)",
                "# scratch tree: work_call_matrix/repo (this folder): analyses/ = <analysis repo>/analyses (symlinks) with recurrence_fork.py as recurrence.py, the 80-study primary_lists.json and universes.json; "
                "data/lists = repo lists + derived_all_2026-10-14/lists/*_fdr05.tsv; data/universe_sets = schizophrenia-transcriptome-replication/data/universe_sets + snrna_universe_2026-10-14 + handcoding_2026-10-14 + "
                "derived__<pmid>.txt rebuilt as the union of derived_all_2026-10-14/universes/*__tested.tsv (Frohlich: the D13 union set union_snrna_1pct__39227716.txt, 22,985 genes); corpus/lineage.json = lineage_4f91.json",
                f"# genes ranked by n_lists_observed desc then gene over all genes of the scenario in recurrence_genes.tsv: gene_rank = that position; rank_within_383 = position among the 383 by the same ordering. {gap_txt}",
                f"# {callset}",
                "# study_id = PMID (or DOI key); studies are in the column order of the manuscript table (named cohort clusters, then year), the 8 empty-list studies after them in the order of base80_studies.tsv",
                f"# rows: {len(rows)} = {len(genes)} genes x {len(studies)} studies; hits (up + down + mixed + hit_no_direction) equal n_lists_observed for all {len(genes)} genes (checked); "
                f"direction counts equal the up/down/both/other_direction columns of recurrence_genes.tsv for all {len(genes)} genes (checked). Changed cells against the first version of this file: {exc}, set from a direction call to not_tested because the gene is on the study's list but outside its tested universe"]
        with open(f"{OUT}/{name}", "w", encoding="utf8", newline="") as fh:
            for h in hdr:
                fh.write(h + "\n")
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            cols = ["gene", "study_id", "call", "gene_rank", "n_lists_observed", "rank_within_383"]
            w.writerow(cols)
            for r in rows:
                w.writerow([r[c] for c in cols])

    write("fig02_b_gene_by_study_long_383.tsv", rows72, order72, "")
    ncalls = collections.Counter(r["call"] for r in rows80 if r["study_id"] in order8)
    e8 = "; empty-list studies (no gene is a hit, so every call is tested_not_hit or not_tested): " + "; ".join(
        f"{label[k]} ({k}): " + ("no tested universe, not_tested for every gene" if M8[k] is None else f"universe M of {len(M8[k])} genes") for k in order8)
    write("fig02_b_gene_by_study_long_80.tsv", rows80, order80, e8)
    print("empty-list studies:", [(k, label[k], None if M8[k] is None else len(M8[k])) for k in order8], dict(ncalls))
    print("72:", [(k, lab.get(k)) for k in order72])
    if dmis:
        print("direction mismatches (first 10):", dmis[:10])


if __name__ == "__main__":
    main()
