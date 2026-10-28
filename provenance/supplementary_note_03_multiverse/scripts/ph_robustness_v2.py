#!/usr/bin/env python3
"""pH enrichment of each forking path's gene list, across the v2 grid (276 designed paths); written to ph_robustness_v2.json.
Reuses ph_robustness.py of forking_paths_rerun_2026-10-09b_sn_bulk_only; additions: block in the key, signed_score_allgenes (direction
from recurrence_genes.tsv for every gene) next to signed_score (recurrent_sets.tsv direction only, as in the earlier grid).

For every factorial path and unit tier, on the path's expression-matched-null gene set (grid_genes.tsv):

- signed pH score: mean of -(SZ sign) x Z_pH over the set's genes with a value, where Z_pH is the per-gene Stouffer Z
  over the five donor-separate datasets' control-only pH association (confounds/results/combined_controls_pH.tsv) and
  SZ sign is the rule-A majority direction (confounds/recurrent_sets.tsv, set study_brain; known for the 133
  study-tier genes, so genes without a direction are dropped). Positive when the genes move with low pH the way they
  move in SZ (confounds_report.md's convention).
- pH-gene O/E: observed over expected membership in the control-pH gene set (genes at combined control pH q < 0.05,
  geneset_enrichment.R's gemma_controls_pH_q05 set), the expectation expression-matched: genes binned in deciles of
  expr_percentile, f_bin the bin's fraction of pH-set genes, expected the sum of f_bin over the set's measured genes.

    python3 analyses/forking_paths_v2_2026-10-10/ph_robustness_v2.py
"""
import collections
import csv
import datetime
import hashlib
import json
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.environ.get("FP_GL", "rerun-full/analyses")   # inputs read from the remote copy of the repo, read-only; outputs are written next to this script (changed from os.path.dirname(HERE))


def read(fn, base=HERE):
    lines = [l for l in open(os.path.join(base, fn), encoding="utf-8") if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:64]


def main():
    grid = read("grid.tsv")
    genes = read("grid_genes.tsv")
    sets = collections.defaultdict(set)
    for g in genes:
        sets[(g["path"], g["block"], g["unit_tier"], g["null"])].add(g["gene"])

    # per-gene control pH results (Stouffer Z over the five donor-separate datasets)
    ph = {}  # gene -> (expr_percentile, q, Z)
    for r in read(os.path.join("results", "combined_controls_pH.tsv"), os.path.join(GL, "confounds")):
        if r["factor"] == "pH" and r["q"] and r["expr_percentile"]:
            ph[r["gene"]] = (float(r["expr_percentile"]), float(r["q"]), float(r["Z"]))
    phset = {g for g, (_, q, _) in ph.items() if q < 0.05}

    # expression-matched expectation: deciles of expr_percentile over the measured genes
    bins = collections.Counter(min(int(p * 10), 9) for p, _, _ in ph.values())
    sig_bins = collections.Counter(min(int(ph[g][0] * 10), 9) for g in phset)
    fbin = {b: c / n for b, (c, n) in ((b, (sig_bins[b], bins[b])) for b in bins) if n}

    # SZ direction, rule-A majority (recurrent_sets.tsv, set study_brain)
    sign = {}
    for r in read("recurrent_sets.tsv", os.path.join(GL, "confounds")):
        if r["set"] == "study_brain" and r["tier"] == "study":
            sign[r["gene"]] = 1 if r["sz_direction"] == "up" else (-1 if r["sz_direction"] == "down" else 0)

    # SZ direction for every gene, rule-A majority of the published-brain list directions (recurrence_genes.tsv; agrees with
    # recurrent_sets.tsv where both exist, asserted below); the grid's own score uses recurrent_sets.tsv only
    rg = {r["gene"]: r for r in read("recurrence_genes.tsv", GL) if r["scenario"] == "published, brain"}
    sign_all = {}
    for g, r in rg.items():
        up, down = int(r["up"]), int(r["down"])
        sign_all[g] = 1 if up > down else (-1 if down > up else 0)
    for g, sg in sign.items():
        if g in sign_all and sg != 0:
            assert sign_all[g] == sg, f"direction mismatch for {g}"

    paths, seen = [], set()
    for r in grid:
        key = (r["block"], r["path"], r["unit_tier"])
        if key in seen:
            continue
        seen.add(key)
        s = sets.get((key[1], key[0], key[2], "expr"), set())
        meas = [g for g in sorted(s) if g in ph]
        exp = sum(fbin[min(int(ph[g][0] * 10), 9)] for g in meas)
        obs = sum(1 for g in meas if g in phset)
        vals = [-sign[g] * ph[g][2] for g in s if g in ph and sign.get(g, 0) != 0]
        vals_all = [-sign_all[g] * ph[g][2] for g in s if g in ph and sign_all.get(g, 0) != 0]
        paths.append({"block": key[0], "path": key[1], "unit_tier": key[2], "n_set": len(s), "n_measured": len(meas),
                      "ph_obs": obs, "ph_exp": round(exp, 1),
                      "ph_oe": round(obs / exp, 2) if exp > 0 else None,
                      "signed_score": round(statistics.mean(vals), 2) if vals else None,
                      "n_direction": len(vals),
                      "signed_score_allgenes": round(statistics.mean(vals_all), 2) if vals_all else None,
                      "n_direction_allgenes": len(vals_all)})

    def summ(ps):
        sc = [p["signed_score"] for p in ps if p["signed_score"] is not None]
        o = [p["ph_oe"] for p in ps if p["ph_oe"] is not None]
        return {"n_paths": len(ps), "n_scored": len(sc), "score_positive": sum(1 for v in sc if v > 0),
                "score_min": round(min(sc), 2), "score_median": round(statistics.median(sc), 2),
                "score_max": round(max(sc), 2), "oe_min": round(min(o), 2),
                "oe_median": round(statistics.median(o), 2), "oe_max": round(max(o), 2)}

    head = [f"# script: analyses/forking_paths_v2_2026-10-10/ph_robustness_v2.py, run {datetime.date.today().isoformat()}",
            f"# input: analyses/forking_paths/grid.tsv sha256 {sha(os.path.join(HERE, 'grid.tsv'))}",
            f"# input: analyses/forking_paths/grid_genes.tsv sha256 {sha(os.path.join(HERE, 'grid_genes.tsv'))}",
            f"# input: analyses/confounds/results/combined_controls_pH.tsv sha256 "
            f"{sha(os.path.join(GL, 'confounds', 'results', 'combined_controls_pH.tsv'))}",
            f"# input: analyses/confounds/recurrent_sets.tsv sha256 "
            f"{sha(os.path.join(GL, 'confounds', 'recurrent_sets.tsv'))}"]
    json.dump({"provenance": head,
               "summary": {"study": summ([p for p in paths if p["unit_tier"] == "study"]),
                           "cohort": summ([p for p in paths if p["unit_tier"] == "cohort"])},
               "paths": paths},
              open(os.path.join(HERE, "ph_robustness_v2.json"), "w"), indent=1)
    s = summ([p for p in paths if p["unit_tier"] == "study"])
    print(f"ph_robustness_v2.json: study tier signed score positive in {s['score_positive']} of {s['n_scored']} scored paths "
          f"(median {s['score_median']}, {s['score_min']} to {s['score_max']}); pH-gene O/E median {s['oe_median']} "
          f"({s['oe_min']} to {s['oe_max']})")


if __name__ == "__main__":
    main()