#!/usr/bin/env python3
"""figure_data/forking_paths_v2/ tables from grid.tsv, grid_genes.tsv, paths_hash.tsv and ph_robustness_v2.json (same shape as
figure_data/fig04: tidy tables with '#' provenance headers).
    python3 export_figure_data.py
"""
import csv, datetime, hashlib, json, os, statistics, collections

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.environ.get("FD_OUT", os.path.expanduser("figure_data/forking_paths_v2"))
os.makedirs(OUT, exist_ok=True)
TODAY = datetime.date.today().isoformat()
SCRIPT = "analyses/multiverse_analysis/export_figure_data.py"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


def rd(fn):
    return list(csv.DictReader((l for l in open(os.path.join(HERE, fn), encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


def write(name, panel, fields, rows, inputs, notes):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        fh.write(f"# figure: forking_paths_v2 table {panel}\n# script: {SCRIPT}, run {TODAY}\n")
        fh.write("# population: published-brain primary studies of recurrence.py (rule A, all lists; 62 studies in the 515-study corpus) and the designed paths of forking_paths_v2.py (blocks 1 to 3; scope labels containing 276, 270, 96 and 180 are the labels of the earlier 276-path grid and are kept so downstream readers still match)\n")
        for i in inputs:
            fh.write(f"# input: analyses/multiverse_analysis/{i} sha256 {sha(os.path.join(HERE, i))}\n")
        for n in notes:
            fh.write(f"# {n}\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(name, len(rows))


def f(x):
    return float(x) if x not in ("", None) else None


def med(v):
    return round(statistics.median(v), 3) if v else ""


grid = rd("grid.tsv")
hashes = rd("paths_hash.tsv")
ph = json.load(open(os.path.join(HERE, "ph_robustness_v2.json")))["paths"]
hid = {}
for r in hashes:
    hid.setdefault(r["hash"], len(hid) + 1)
FORKS = ["rule", "threshold", "empty", "pseudo_noexpr", "noncoding", "single_cell", "cut", "study_filter", "list_level"]

# a: design
cnt = collections.defaultdict(collections.Counter)
for r in hashes:
    cfg = r["path"].split("|")
    for fk, v in zip(FORKS, cfg):
        cnt[fk][v] += 1
desc = {"rule": "list selection rule (R2 dropped: selects the same lists as R1 in 326 of 327 studies)",
        "threshold": "list threshold tier; any criterion = 'all'", "empty": "empty lists kept or dropped; varied under fdr only",
        "pseudo_noexpr": "pseudogenes and genes without brain expression; gtex = median TPM >= 1 in >= 1 of 13 GTEx v8 brain regions",
        "noncoding": "HGNC non-coding RNAs dropped or pooled", "single_cell": "single-nucleus studies without or with",
        "cut": "published list | tight (reported FDR < 0.01 else p < 0.001) | |log2FC| >= 0.26 | >= 0.58 | capped at 200 genes by P",
        "study_filter": "all | N >= 20 | 2016 on | cortex only | no severe-error study | adjusted for >= 1 covariate",
        "list_level": "inclusive (default) | partial DE lists dropped (Scarr, Arion, Altar, Shao) | non-DE-style lists dropped (Ling, Kohen) | mixed-diagnosis dropped (Pai, Frohlich) | Chang from Table S2 only | Pai as the three text-named genes | all three drops"}
rows = [{"fork": fk, "setting": v, "n_designed_paths": n, "description": desc[fk] if i == 0 else ""}
        for fk in FORKS for i, (v, n) in enumerate(cnt[fk].items())]
write("forking_v2_a_design.tsv", "a", ["fork", "setting", "n_designed_paths", "description"], rows, ["paths_hash.tsv"],
      ["fork settings counted over the designed paths (block 1: 96, block 2: 180, block 3: 36; 6 specifications occur in both blocks 1 and 2)"])

# b: spec curve
out = []
for i, r in enumerate(grid, 1):
    nl = r["null"] if r["null"] != "de_prior" else "de_prior_" + r["de_prior_source"]
    out.append({"row_id": i, "block": r["block"], "path_id": r["path"], "hash_id": hid[r["hash"]], **{k: r[k] for k in FORKS},
                "tier": r["unit_tier"], "null": nl, "n_studies": r["n_studies"], "n_clusters": r["n_clusters"],
                "obs_ge3": r["obs_ge3"], "exp_ge3": r["exp_ge3"], "obs_over_exp_ge3": r["ratio_ge3"],
                "obs_ge5": r["obs_ge5"], "exp_ge5": r["exp_ge5"], "obs_over_exp_ge5": r["ratio_ge5"],
                "n_q01": r["genes_q10_expr"], "n_q01_uniform": r["genes_q10_uniform"]})
write("forking_v2_b_spec_curve.tsv", "b", ["row_id", "block", "path_id", "hash_id"] + FORKS + ["tier", "null", "n_studies", "n_clusters", "obs_ge3", "exp_ge3",
      "obs_over_exp_ge3", "obs_ge5", "exp_ge5", "obs_over_exp_ge5", "n_q01", "n_q01_uniform"], out, ["grid.tsv", "paths_hash.tsv"],
      ["one row per designed path x unit tier x null; tier study = genes in >= C lists, cohort = genes in >= C cohort clusters; null = uniform | expr | de_prior_ours | de_prior_global",
       "obs_over_exp_geC = observed / expected genes in >= C units under the row's null (blank where expected = 0); n_q01 = genes at BH q < 0.1 under the expression-matched null, n_q01_uniform under the uniform null",
       "hash_id = index of the distinct analysis (paths with equal post-fork inputs share a hash_id and equal outcome columns); block 2 paths with cut published and study_filter all repeat block 1 paths"])

# c: pH
phrows = [{"block": p["block"], "path_id": p["path"], "hash_id": hid.get(next((h["hash"] for h in hashes if h["path"] == p["path"] and h["block"] == p["block"]), ""), ""),
           "tier": p["unit_tier"], "signed_ph_score": p["signed_score"], "n_direction": p["n_direction"],
           "signed_ph_score_allgenes": p["signed_score_allgenes"], "n_direction_allgenes": p["n_direction_allgenes"],
           "n_set": p["n_set"], "n_measured": p["n_measured"], "ph_obs": p["ph_obs"], "ph_exp": p["ph_exp"], "ph_oe": p["ph_oe"]} for p in ph]
write("forking_v2_c_ph_paths.tsv", "c", ["block", "path_id", "hash_id", "tier", "signed_ph_score", "n_direction", "signed_ph_score_allgenes",
      "n_direction_allgenes", "n_set", "n_measured", "ph_obs", "ph_exp", "ph_oe"], phrows, ["ph_robustness_v2.json", "grid.tsv"],
      ["signed_ph_score = mean of -(SZ sign) x control-pH Z over the path's expression-matched-null gene set, as ph_robustness.py of the 2026-10-09b grid: SZ sign from confounds/recurrent_sets.tsv (study_brain), genes without a sign dropped (n_direction)",
       "signed_ph_score_allgenes = the same with the SZ sign (rule-A majority list direction) from recurrence_genes.tsv for every gene that has one (n_direction_allgenes)",
       "ph_obs / ph_exp / ph_oe = observed and expected genes at combined control pH q < 0.05 and their ratio; blank score = no gene with sign and pH value"])

# d: summary of ratios
def srows(scope, sel):
    res = []
    for tier in ("study", "cohort"):
        for nl in ("uniform", "expr", "de_prior_ours", "de_prior_global"):
            for c in (5, 3):
                v = []
                nblank = 0
                for r in sel:
                    rn = r["null"] if r["null"] != "de_prior" else "de_prior_" + r["de_prior_source"]
                    if r["unit_tier"] == tier and rn == nl:
                        x = f(r[f"ratio_ge{c}"])
                        if x is None:
                            nblank += 1
                        else:
                            v.append(x)
                res.append({"scope": scope, "tier": tier, "null": nl, "cutoff": f"{c}+", "n_paths": len(v) + nblank, "n_with_ratio": len(v),
                            "median": med(v), "min": min(v) if v else "", "max": max(v) if v else "",
                            "n_below_1": sum(1 for x in v if x < 1), "n_above_1": sum(1 for x in v if x > 1)})
    return res


def pooled(scope, sel):
    res = []
    for c in (5, 3, 4):
        v = [f(r[f"ratio_ge{c}"]) for r in sel if f(r[f"ratio_ge{c}"]) is not None]
        res.append({"scope": scope, "tier": "all rows", "null": "all", "cutoff": f"{c}+", "n_paths": len(v), "n_with_ratio": len(v), "median": med(v),
                    "min": min(v), "max": max(v), "n_below_1": sum(1 for x in v if x < 1), "n_above_1": sum(1 for x in v if x > 1)})
    return res


seen_spec, uniq = set(), []
for r in grid:
    k = (r["path"], r["unit_tier"], r["null"], r["de_prior_source"])
    if k not in seen_spec:
        seen_spec.add(k)
        uniq.append(r)
seen_h, dist = set(), []
for r in grid:
    k = (r["hash"], r["unit_tier"], r["null"], r["de_prior_source"])
    if k not in seen_h:
        seen_h.add(k)
        dist.append(r)
summ = []
for scope, sel in (("all_designed_paths_276", grid), ("unique_specifications_270", uniq), ("distinct_analyses", dist),
                   ("block1_96", [r for r in grid if r["block"] == "block1"]), ("block2_180", [r for r in grid if r["block"] == "block2"])):
    summ += pooled(scope + "_pooled_rows", sel) + srows(scope, sel)
write("forking_v2_d_summary.tsv", "d", ["scope", "tier", "null", "cutoff", "n_paths", "n_with_ratio", "median", "min", "max", "n_below_1", "n_above_1"], summ,
      ["grid.tsv"],
      ["observed / expected ratio of genes in >= 5 and >= 3 lists (tier study) or cohort clusters (tier cohort) under each chance model; median, minimum and maximum over the scope's paths",
       "scope: all_designed_paths_276 = every designed path; unique_specifications_270 = block 2's 6 repeats of block 1 paths counted once; distinct_analyses = one path per distinct hash of post-fork inputs",
       "*_pooled_rows = all unit tiers and nulls together (the earlier grid's 'ratio_ge3 min 0.521, median 0.931, max 1.801' pools its 2,304 factorial rows this way); n_paths there counts rows"])

# e: per path
byp = collections.defaultdict(dict)
for r in grid:
    nl = r["null"] if r["null"] != "de_prior" else "de_prior_" + r["de_prior_source"]
    byp[(r["block"], r["path"])][(r["unit_tier"], nl)] = r
phd = {(p["block"], p["path"], p["unit_tier"]): p for p in ph}
rows = []
for h in hashes:
    d = byp[(h["block"], h["path"])]
    s_u, c_u = d[("study", "uniform")], d[("cohort", "uniform")]
    row = {"block": h["block"], "path_id": h["path"], "hash": h["hash"], "hash_id": hid[h["hash"]], **dict(zip(FORKS, h["path"].split("|"))),
           "n_studies": s_u["n_studies"], "n_clusters": s_u["n_clusters"],
           "recurrent_genes_study_q10_uniform": s_u["genes_q10_uniform"], "recurrent_genes_study_q10_expr": s_u["genes_q10_expr"],
           "recurrent_genes_cohort_q10_uniform": c_u["genes_q10_uniform"], "recurrent_genes_cohort_q10_expr": c_u["genes_q10_expr"],
           "genes_in_5plus_lists": s_u["obs_ge5"], "genes_in_3plus_lists": s_u["obs_ge3"],
           "genes_in_5plus_clusters": c_u["obs_ge5"], "genes_in_3plus_clusters": c_u["obs_ge3"]}
    for t in ("study", "cohort"):
        p = phd.get((h["block"], h["path"], t), {})
        row[f"signed_ph_score_{t}"] = p.get("signed_score", "")
        row[f"signed_ph_score_allgenes_{t}"] = p.get("signed_score_allgenes", "")
    rows.append(row)
write("forking_v2_e_paths.tsv", "e", ["block", "path_id", "hash", "hash_id"] + FORKS + ["n_studies", "n_clusters", "recurrent_genes_study_q10_uniform",
      "recurrent_genes_study_q10_expr", "recurrent_genes_cohort_q10_uniform", "recurrent_genes_cohort_q10_expr", "genes_in_5plus_lists", "genes_in_3plus_lists",
      "genes_in_5plus_clusters", "genes_in_3plus_clusters", "signed_ph_score_study", "signed_ph_score_allgenes_study", "signed_ph_score_cohort",
      "signed_ph_score_allgenes_cohort"], rows, ["grid.tsv", "paths_hash.tsv", "ph_robustness_v2.json"],
      ["one row per designed path; recurrent genes = genes at BH q < 0.1 under the uniform and the expression-matched null (the earlier grid's 'headline' set), counted by lists (study tier) and by cohort clusters (cohort tier)",
       "signed_ph_score_*: see table c"])
