#!/usr/bin/env python3
"""Per-panel tidy exports for manuscript Figures 3 (agreement within and across cohorts) and 5 (controls).

Request: MAB_TO_RA_2026_10_08_FRESH_FIGURES_PANEL_DATA.md (Figs 3 and 5).
Reads existing analysis outputs only (no model fits, no permutation reruns, no network). The Fig 5 positive
controls are age, sex and donor identification; no T1D data. Writes figure_data/fig03/
and fig05/, each TSV with '#' provenance lines, plus _meta.json and STATUS.md. One function per panel; rerun with

    python3 analyses/figure_data_exports/export_fig03_fig05.py
"""
import csv
import hashlib
import json
import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)            # analyses/
REPO = os.path.dirname(AN)
RUN = "2026-10-08"
SCRIPT = "analyses/figure_data_exports/export_fig03_fig05.py"
OUT = os.path.expanduser("figure_data")
POP = ("Fig 3 and 5 use the sample-level Gemma/LIBD/PsychAD datasets (not the 64-study list population), except "
       "Fig 3b (published different-cohort pairs of the 64-study population, h1_pairs.tsv) and Fig 5c negative control "
       "(copy of fig02 panel c; recurrence.tsv sha256 d11761227abd)")

_sha_cache = {}


def sha(path):
    path = os.path.abspath(path)
    if path not in _sha_cache:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _sha_cache[path] = h.hexdigest()
    return _sha_cache[path]


def rel(path):
    path = os.path.abspath(path)
    if path.startswith(REPO):
        return os.path.relpath(path, REPO)
    return path.replace(os.path.expanduser("~"), "~")


def rd(path):
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


class Figure:
    def __init__(self, fig):
        self.fig = fig
        self.dir = os.path.join(OUT, fig)
        os.makedirs(self.dir, exist_ok=True)
        self.meta = {"script": SCRIPT, "run_id": RUN, "figure": fig, "population": POP, "inputs": {}, "panels": {}}
        self.status = []

    def note(self, panel, text):
        self.status.append(f"- **{panel}**: {text}")

    def write(self, name, panel, df, inputs, comments=()):
        for p in inputs:
            self.meta["inputs"][rel(p)] = sha(p)
        lines = [f"# figure: {self.fig} panel {panel}", f"# script: {SCRIPT}, run {RUN}", f"# population: {POP}"]
        lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs]
        lines += [f"# {c}" for c in comments]
        path = os.path.join(self.dir, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            for l in lines:
                fh.write(l + "\n")
            df.to_csv(fh, sep="\t", index=False, lineterminator="\n")
        self.meta["panels"][name] = {"panel": panel, "rows": int(len(df))}
        print(f"{self.fig}/{name}: {len(df)} rows")
        return len(df)

    def finish(self):
        with open(os.path.join(self.dir, "_meta.json"), "w") as fh:
            json.dump(self.meta, fh, indent=1)
        with open(os.path.join(self.dir, "STATUS.md"), "w") as fh:
            fh.write(f"# {self.fig} status ({RUN})\n\n" + "\n".join(self.status) + "\n")


# ---------------------------------------------------------------- inputs
SH = os.path.join(AN, "split_halves", "results")
P_CTL, P_CROSS, P_WITHIN = (os.path.join(SH, f) for f in ("controls_summary.tsv", "cross_summary.tsv", "within_summary.tsv"))
P_PUBPAIRS = os.path.join(SH, "published_pairs.tsv")
P_PUBSUM = os.path.join(SH, "published_summary.tsv")
P_A4G = os.path.join(AN, "fairness_lists", "a4_pairs_gemma.tsv")
P_A4S = os.path.join(AN, "fairness_lists", "a4_summary.tsv")
CM = os.path.join(AN, "casemix", "results")
P_CMPAIRS, P_CMUNITS, P_CMT1, P_CMT2, P_CMT2C = (os.path.join(CM, f) for f in (
    "casemix_pairs_full.tsv", "casemix_units.tsv", "casemix_test1.tsv", "casemix_test2.tsv", "casemix_test2_counts.tsv"))
P_SUBREL = os.path.join(AN, "subject_relations", "subject_relations.tsv")
PIT = os.path.join(AN, "fairness_sample_level", "same_donors_pitt", "results")
P_G = os.path.join(PIT, "g_individuals_summary.tsv")
P_H = os.path.join(PIT, "h_stanley_summary.tsv")
P_I = os.path.join(PIT, "i_dx_neighbours_summary.tsv")
P_E_PSY = os.path.join(AN, "psychad_cellxgene", "results", "e_individuals_identify_auroc_summary.tsv")
P_E_LIBD = os.path.join(AN, "libd_brainseq", "analysis", "results", "e_individuals_regions_identify_auroc_summary.tsv")
P_CAS_PAIRS = os.path.join(AN, "controls_age_sex", "pairs.tsv")
P_CAS_MATCHED = os.path.join(AN, "controls_age_sex", "pairs_matched.tsv")
P_CAS_SUM = os.path.join(AN, "controls_age_sex", "summary.tsv")
CAE = os.path.join(AN, "controls_age_external", "results")
P_CAE_DS, P_CAE_SUM, P_CAE_COMB = (os.path.join(CAE, f) for f in ("per_dataset.tsv", "summary.tsv", "combined.tsv"))
P_FIG02_C = os.path.join(OUT, "fig02", "fig02_c_pairwise_fold_by_era.tsv")
P_FIG02_CP = os.path.join(OUT, "fig02", "fig02_c_pooled_by_era.tsv")

CROSS_LABEL = {"same_subjects": "LIBD same subjects, two regions",
               "separate_subjects": "LIBD subject-separate halves",
               "separate_cohorts": "different cohorts",
               "separate_cohorts_gemma": "different cohorts"}
WITHIN_LABEL = "within cohort (split halves)"
CONTRAST_NAME = {"sz": "dx", "age": "age", "male": "sex"}


# ================================================================== Figure 3
def fig3_acd(F):
    """a/c/d: agreement (Spearman of t) between two matched-size draws, by pairing kind, contrast and n per group."""
    C, X = rd(P_CTL), rd(P_CROSS)
    w = C[(C.gene_set == "auto") & C.contrast.isin(CONTRAST_NAME)].copy()
    w["comparison"] = WITHIN_LABEL
    w["unit"] = w.dataset
    w["libd_dlpfc"] = w.dataset.isin(["LIBD_Phase2_DLPFC", "LIBD_Phase1_DLPFC"]).map({True: "yes", False: "no"})
    x = X[(X.part == "controls") & (X.gene_set == "auto") & X.contrast.isin(CONTRAST_NAME)].copy()
    x["comparison"] = x.kind.map(CROSS_LABEL)
    x["unit"] = x.pair
    x["libd_dlpfc"] = ""
    keep = ["comparison", "kind", "unit", "libd_dlpfc", "contrast", "level", "rho_t", "rho_t_q25", "rho_t_q75", "pi1_excess", "top100_fold"]
    w["kind"] = "split_halves"
    u = pd.concat([w.reindex(columns=keep), x.reindex(columns=keep)], ignore_index=True)
    u = u[u.level.astype(str).isin(["10", "20", "30", "50"])].copy()
    u["level"] = u.level.astype(int)
    u["trait"] = u.contrast.map(CONTRAST_NAME)
    u = u.rename(columns={"level": "n_per_group"}).drop(columns="contrast")
    u["panels"] = np.where(u.comparison == WITHIN_LABEL,
                           np.where(u.libd_dlpfc == "yes", "a, c, d (blue line); other datasets: context", "context"),
                           np.where(u.comparison == "different cohorts", "a, c, d", "a, d"))
    u = u[["comparison", "kind", "unit", "libd_dlpfc", "trait", "n_per_group", "rho_t", "rho_t_q25", "rho_t_q75",
           "pi1_excess", "top100_fold", "panels"]].sort_values(["comparison", "trait", "unit", "n_per_group"])
    F.write("fig03_a_c_d_agreement_units.tsv", "a, c, d", u, [P_CTL, P_CROSS],
            ["one row = median over the repeats of one dataset (within) or one dataset pair (cross) at one n per group; "
             "Spearman of limma t on autosomal genes (gene_set auto), the d_controls design (50 repeats); trait dx = SZ vs control, "
             "age = continuous, sex = male vs female. Blue line of a/c/d = libd_dlpfc yes (LIBD Phase 1 and 2 DLPFC)",
             "kinds: split_halves; same_subjects (dashed orange); separate_subjects (solid orange); separate_cohorts (red, LIBD vs non-LIBD pairs); "
             "separate_cohorts_gemma has no controls-design rows (SZ only), so it is absent here"])
    # summary (medians over units), restricted to the LIBD DLPFC split halves for the within family
    s = u[(u.comparison != WITHIN_LABEL) | (u.libd_dlpfc == "yes")]
    s = (s.groupby(["comparison", "trait", "n_per_group"])
          .agg(median_rho_t=("rho_t", "median"), q25=("rho_t", lambda v: v.quantile(.25)), q75=("rho_t", lambda v: v.quantile(.75)),
               n_units=("rho_t", "size")).reset_index().round(4))
    F.write("fig03_a_c_d_agreement_summary.tsv", "a, c, d", s, [P_CTL, P_CROSS],
            ["medians over units; the within-cohort family is LIBD Phase 1 and Phase 2 DLPFC only, as in the legend (blue)"])
    g = s.set_index(["comparison", "trait", "n_per_group"]).median_rho_t
    chk = []
    for key, old in (((("within cohort (split halves)", "dx", 30)), 0.34), (("different cohorts", "dx", 30), 0.06),
                     (("different cohorts", "dx", 10), 0.10), (("different cohorts", "age", 30), 0.41)):
        chk.append(f"{key[0]} {key[1]} n={key[2]}: {g.get(key):.2f} (report {old})")
    F.note("a, c, d", "medians reproduce split_halves_report.md ('Different cohorts at matched n'): " + "; ".join(chk)
           + "; pair counts of different cohorts 8, 4, 2 at 10, 20, 30 (n_units) as in the legend of c.")
    F.note("a, c, d", "STALE EXPORT RERUN: the earlier figure_data/fig3b/ files (2026-10-01, controls_summary sha c36a3ad40e9e, "
           "cross_summary sha 6612f8bded49) predate the split-halves rebuild of 2026-10-08 14:51 (now 208048b6c1bb and 9eab31cace31). RA reran the repo's "
           "analyses/split_half_reliability_of_datasets/scripts/export_fig3b_mab.py (pandas only, seconds), which overwrote fig3b/ with the current values; the "
           "fig03 files above read the same two inputs. The 2026-10-01 fig3b files are not kept.")
    F.note("d", "the legend's four pairings are comparison = same subjects / split halves / subject-separate / different cohorts (trait dx).")
    n_sh = u[(u.comparison == WITHIN_LABEL)].unit.nunique()
    F.note("a", f"the within-cohort family has {n_sh} datasets in the units file (context); the legend names LIBD DLPFC only (libd_dlpfc yes).")


def fig3_b(F):
    """b: published different-cohort pairs against LIBD/Gemma split halves of the same n."""
    p = rd(P_PUBPAIRS)
    p["pair_id"] = p.study_a.astype(str) + "_" + p.study_b.astype(str)
    rows = []
    for r in p.itertuples():
        base = dict(pair_id=r.pair_id, label_a=r.label_a, label_b=r.label_b, N_a=r.N_a, N_b=r.N_b, n_matched=r.n_eq,
                    level=r.level, n_datasets=r.n_datasets)
        if pd.notna(r.rho_logfc):
            rows.append({**base, "metric": "spearman_rho_log2fc", "value": r.rho_logfc,
                         "libd_percentile": r.rho_pct_median_over_datasets,
                         "libd_percentile_min": r.rho_pct_min_over_datasets, "libd_percentile_max": r.rho_pct_max_over_datasets,
                         "cross_percentile": r.rho_pct_in_cross})
        if pd.notna(r.top_fold):
            rows.append({**base, "metric": f"top{int(r.top_n)}_fold" if pd.notna(r.top_n) else "top_fold", "value": r.top_fold,
                         "libd_percentile": r.fold_pct_median_over_datasets, "libd_percentile_min": None,
                         "libd_percentile_max": None, "cross_percentile": r.fold_pct_in_cross})
    d = pd.DataFrame(rows)
    d["metric_group"] = d.metric.where(d.metric == "spearman_rho_log2fc", "top_n_fold")
    F.write("fig03_b_published_vs_split_halves.tsv", "b", d, [P_PUBPAIRS, P_PUBSUM],
            ["published different-cohort pairs (h1_pairs.tsv scenario 'published, brain', cohort different; 1,369 pairs), n_matched = min(N_a, N_b) / 2 "
             "at the nearest split-half level; libd_percentile = percentile of the published value among each dataset's split-half repeats at that "
             "level (median over datasets that reach it; the 16 datasets of split_halves include 11 Gemma ones, not LIBD only); cross_percentile = among "
             "pooled different-cohort pair repeats; rows only where the metric exists. top_n_fold: top 20 / 50 / 100 chosen by the pair, metric column names n"])
    ps = rd(P_PUBSUM)
    a = ps[ps.level.astype(str) == "all"].iloc[0]
    rho = d[d.metric == "spearman_rho_log2fc"]
    F.note("b", f"{len(rho)} pairs with a rho (published table on both sides; summary 'all' row n_rho {int(a.n_rho)}), median rho {rho.value.median():.3f} "
                f"(published_summary.tsv {a.pub_rho_logfc_median:.3f}), median libd_percentile {rho.libd_percentile.median():.2f} "
                f"(summary {a.rho_pct_median_over_datasets:.2f}); {int((d.metric != 'spearman_rho_log2fc').sum())} top-fold rows.")
    F.note("b", "MISMATCH with Figures.md Fig. 3b legend: 'the 15 pairs with a published table on both sides sit at the 22nd percentile of LIBD split halves'. "
                f"split_halves_report.md and published_summary.tsv (2026-10-08): {int(a.n_rho)} different-cohort pairs with full tables, median log2 FC rho 0.07, "
                "median percentile 0.31 among the split halves of the datasets at their n (16 datasets, LIBD and Gemma) and 0.47 among different-cohort pairs; 0 pairs "
                "below the 5th percentile of every dataset's split halves. No set of 15 pairs or a 22nd percentile exists in these files.")
    F.note("b", "The legend says 'six agreement metrics'; split_halves/results/published_pairs.tsv carries two (Spearman rho of log2 FC, top-N fold). "
                "Other metrics (pi1, top-100 sign) exist for the split halves only (within_summary.tsv), not per published pair. Not faked.")


def fig3_e(F):
    """e: Spearman of Gemma log2 fold changes by subject relation (fairness_lists a4 part B)."""
    d = rd(P_A4G)
    d["pair_id"] = d.dataset_a + "_" + d.dataset_b
    o = d[["pair_id", "dataset_a", "dataset_b", "series_a", "series_b", "bank_a", "bank_b", "donor_relation", "relation_text",
           "evidence", "k_shared", "region_relation", "n_universe", "rho", "top_overlap", "top_expected", "top_fold"]].rename(
        columns={"donor_relation": "relation", "rho": "rho_log2fc"})
    F.write("fig03_e_gemma_rho_by_relation.tsv", "e", o, [P_A4G, P_A4S],
            ["Spearman rho of the two Gemma result sets' log2 fold changes over the genes both tested (n_universe); 35 brain bulk or laser-capture result sets "
             "from 23 GEO series, 595 pairs; relation = a4_donors.py donor_relation from subject_relations (built 2026-10-07 for this run); "
             "top_fold = top-100 overlap over chance"])
    med = o.groupby("relation").rho_log2fc.agg(["count", "median"])
    F.note("e", "medians by relation: " + "; ".join(f"{k} {v['median']:.2f} ({int(v['count'])})" for k, v in med.iterrows()))
    F.note("e", "legend check: same subjects in one GEO series 0.35 (11 pairs, 'different region, same series'): match; partial subject overlap 0.13 (22 pairs, across series): match; "
                "subject-separate different banks 0.11 (260): match. 'same subjects, same region, same series' (2 pairs, 0.41) is a separate relation in the file; the legend's "
                "'same subjects in one GEO series 0.35' is the 11-pair different-region group only.")
    F.note("e", "a4_pairs_gemma.tsv was written 2026-10-07 21:28, before the 10-08 core run; its meta lists primary_lists.json e282904dcacb and universes.json 919c00303ca0, "
                "equal to the current files. subject_relations.tsv has since been rebuilt (sha now a089c58fe494; the run used 5a33c5135851): relation labels may differ "
                "for some pairs; the labels in this file are those of the run.")


def fig3_f(F):
    """f: case mix. Left: rho(t) adjusted for n and technology against |difference in SZ median age|; right: comparable vs random subsets."""
    pr, un = rd(P_CMPAIRS), rd(P_CMUNITS)
    sep = pr[pr.relation == "subject-separate"].copy()
    age = un.set_index("unit").age_sz
    sep["d_median_age"] = (sep.a.map(age) - sep.b.map(age)).abs()
    X = np.column_stack([np.ones(len(sep)), np.log(sep.neff_a), np.log(sep.neff_b), (sep.tech_a == sep.tech_b).astype(float)])
    beta, *_ = np.linalg.lstsq(X, sep.rho_t.values, rcond=None)
    sep["rho_adj"] = sep.rho_t.values - X @ beta
    sep["mount_sinai_pair"] = np.where((sep.a == "PsychAD_MSSM") | (sep.b == "PsychAD_MSSM"), "yes", "no")
    sep["pair_id"] = sep.a + "_" + sep.b
    left = sep[["pair_id", "a", "b", "d_median_age", "rho_t", "rho_adj", "top100_fold", "neff_a", "neff_b", "tech_a", "tech_b", "mount_sinai_pair"]]
    F.write("fig03_f_casemix_pairs.tsv", "f left", left, [P_CMPAIRS, P_CMUNITS, P_CMT1],
            ["subject-separate pairs only (other relations are not tested); d_median_age = |SZ median age at death a - b| (casemix_units.tsv age_sz); "
             "rho_adj = residual of rho_t from an OLS on log neff_a + log neff_b + same technology, as described in casemix_report.md Methods (recomputed here "
             "by least squares; the R script's own residuals are not saved)"])
    rho = pd.Series(sep.d_median_age).rank().corr(pd.Series(sep.rho_adj).rank())
    t1 = rd(P_CMT1)
    ref = t1[(t1.set.str.contains("PsychAD")) & (t1.measure == "rho_t") & (t1.descriptor == "age_sz")]
    F.note("f left", f"{len(sep)} subject-separate pairs; Spearman of rho_adj with d_median_age recomputed {rho:.3f}; casemix_test1.tsv "
                     f"({ref.set.iloc[0] if len(ref) else 'n/a'}): {ref.rho.iloc[0]:.3f}, one-sided p {ref.p_less.iloc[0]:.3f}, two-sided p {ref.p_two.iloc[0]:.3f}. "
                     f"Legend: 281 pairs, Spearman -0.15, P = 0.12: matches the one-sided p (0.116); the file's two-sided p is {ref.p_two.iloc[0]:.2f}, so the legend should say one-sided.")
    t2 = rd(P_CMT2)
    t2 = t2.rename(columns={"observed": "comparable_subsets", "random_median": "random_median", "full_data": "full_data"})
    F.write("fig03_f_comparable_vs_random.tsv", "f right", t2, [P_CMT2, P_CMT2C],
            ["comparable_subsets = statistic refitted on the condition's subjects (the dot); random_* = 100 random subsets with the same SZ and control counts per unit "
             "(5th, median, 95th percentile); full_data = x; condition age65 = age at death <= 65, suicide = suicides excluded where recorded, both"])
    cnt = rd(P_CMT2C)
    cnt_path = os.path.join(F.dir, "fig03_f_comparable_subset_counts.tsv")
    F.write("fig03_f_comparable_subset_counts.tsv", "f right", cnt, [P_CMT2C], ["subjects kept per unit and condition (n_sz, n_ctl) against full data"])
    # staleness: the relation labels the casemix run used
    sr = rd(P_SUBREL)
    rel_now = {}
    for r in sr.itertuples():
        rel_now[(r.dataset_a, r.dataset_b)] = r.relation_group
        rel_now[(r.dataset_b, r.dataset_a)] = r.relation_group
    diff = [(a, b) for a, b in zip(sep.a, sep.b) if rel_now.get((a, b)) not in (None, "subject-separate")]
    F.note("f", f"casemix results built 2026-10-01 against subject_relations sha 2c84e9aa454e; the table is now a089c58fe494. {len(diff)} of the 281 subject-separate pairs carry a "
                "different relation in the current table" + (f" ({'; '.join(a + ' x ' + b for a, b in diff[:6])})" if diff else "") + ". "
                "Pair 'a' counts of the casemix run are used unchanged.")


def _ident_rows():
    """Per relation: identification AUROC (nodx residuals) in both directions, from the four sources."""
    rows = []
    e1 = rd(P_E_LIBD)
    e1 = e1[e1["mode"] == "nodx"].set_index("direction")
    rows.append(dict(source="LIBD regions", relation="LIBD DLPFC x hippocampus", subjects="all", n_candidates=int(e1.n_candidates.iloc[0]),
                     auroc_ab=e1.mean_auroc["ab"], auroc_ba=e1.mean_auroc["ba"], null_ab=e1.null_mean_auroc["ab"], null_ba=e1.null_mean_auroc["ba"],
                     p_ab=e1.p_perm["ab"], p_ba=e1.p_perm["ba"]))
    e2 = rd(P_E_PSY)
    e2 = e2[(e2.level == "all") & (e2["mode"] == "nodx")].set_index("direction")
    rows.append(dict(source="PsychAD x LIBD", relation="HBCC snRNA-seq pseudobulk x LIBD DLPFC", subjects="all", n_candidates=int(e2.n_candidates.iloc[0]),
                     auroc_ab=e2.mean_auroc["ab"], auroc_ba=e2.mean_auroc["ba"], null_ab=e2.null_mean_auroc["ab"], null_ba=e2.null_mean_auroc["ba"],
                     p_ab=e2.p_perm["ab"], p_ba=e2.p_perm["ba"]))
    g = rd(P_G)
    g = g[g["mode"] == "nodx"]
    for (a, b), t in g.groupby(["unit_a", "unit_b"], sort=False):
        t = t.set_index("direction")
        rows.append(dict(source="Pittsburgh", relation=f"{a} x {b}", subjects="all", n_candidates=int(t.n_candidates.iloc[0]),
                         auroc_ab=t.mean_auroc["ab"], auroc_ba=t.mean_auroc["ba"], null_ab=t.null_mean_auroc["ab"], null_ba=t.null_mean_auroc["ba"],
                         p_ab=t.p_perm["ab"], p_ba=t.p_perm["ba"], region_relation=t.region_relation.iloc[0], platform_relation=t.platform_relation.iloc[0]))
    h = rd(P_H)
    for (coll, subj, a, b), t in h.groupby(["collection", "subjects", "study_a", "study_b"], sort=False):
        t = t.set_index("direction")
        rows.append(dict(source="Stanley " + ("Array Collection" if coll.startswith("Array") else "Consortium"), relation=f"{a} x {b}", subjects=subj,
                         n_candidates=int(t.n_candidates.iloc[0]), auroc_ab=t.mean_auroc["ab"], auroc_ba=t.mean_auroc["ba"],
                         null_ab=t.null_mean_auroc["ab"], null_ba=t.null_mean_auroc["ba"], p_ab=t.p_perm["ab"], p_ba=t.p_perm["ba"],
                         region_relation=t.region_relation.iloc[0], platform_relation=t.platform_relation.iloc[0]))
    d = pd.DataFrame(rows)
    d["auroc_identity"] = (d.auroc_ab + d.auroc_ba) / 2
    d["above_null_both_p05"] = ((d.p_ab < 0.05) & (d.p_ba < 0.05)).map({True: "yes", False: "no"})
    return d


def fig3_g(F):
    """g: identification of a donor across datasets (AUROC) against diagnosis AUROC of the nearest wrong candidate."""
    idr = _ident_rows()
    dx = rd(P_I)
    dx = dx[(dx["mode"] == "nodx") & (dx.donors == "all")]
    dxr = {}
    for (src, relname), t in dx.groupby(["source", "relation"]):
        t = t.set_index("direction")
        if relname == "POOLED":
            continue
        dxr[(src, relname)] = (t.dx_auroc["ab"], t.dx_auroc["ba"], t.null_dx_auroc["ab"], t.p_auroc["ab"], t.p_auroc["ba"])
    out = []
    for r in idr.itertuples():
        key = (r.source, r.relation)
        v = dxr.get(key)
        if r.source.startswith("Stanley") and r.subjects != "sz_ctrl":
            v = None   # the diagnosis AUROC is computed on SZ and controls only
        out.append(dict(source=r.source, relation=r.relation, subjects=r.subjects, n_candidates=r.n_candidates,
                        auroc_identity=round(r.auroc_identity, 4), auroc_identity_ab=r.auroc_ab, auroc_identity_ba=r.auroc_ba,
                        p_identity_ab=r.p_ab, p_identity_ba=r.p_ba, above_null_both_p05=r.above_null_both_p05,
                        auroc_dx=(round((v[0] + v[1]) / 2, 4) if v else None), auroc_dx_ab=v[0] if v else None, auroc_dx_ba=v[1] if v else None,
                        auroc_dx_null=v[2] if v else None, p_dx_ab=v[3] if v else None, p_dx_ba=v[4] if v else None))
    o = pd.DataFrame(out)
    F.write("fig03_g_identification_by_relation.tsv", "g", o, [P_G, P_H, P_I, P_E_PSY, P_E_LIBD],
            ["relation = a pair of datasets with donors in both; auroc_identity = mean over donors of 1 - (rank - 1)/(n - 1) of the donor's own profile among the "
             "candidates in the other dataset, averaged over the two directions (residuals without diagnosis, 'nodx'); auroc_dx = same for 'a nearest wrong candidate "
             "shares the donor's diagnosis' (SZ and controls, true partner removed, 'nodx', all donors; chance 0.5); Stanley rows: subjects all_dx or sz_ctrl, auroc_dx "
             "only for sz_ctrl"])
    pool = dx[dx.relation == "POOLED"] if "POOLED" in set(rd(P_I).relation) else None
    allp = rd(P_I)
    pool = allp[(allp.relation == "POOLED") & (allp["mode"] == "nodx") & (allp.donors == "all")]
    pool = pool[["source", "direction", "n_relations", "dx_auroc", "null_dx_auroc", "p_auroc"]]
    F.write("fig03_g_dx_auroc_pooled_by_source.tsv", "g", pool, [P_I], ["POOLED rows of i_dx_neighbours_summary.tsv (nodx, all donors), one row per source and direction"])
    core = o[(o.source.isin(["LIBD regions", "PsychAD x LIBD", "Pittsburgh"]))]
    F.note("g", f"{len(core)} relations in the 'existing two + Pittsburgh 28' set: mean identification AUROC {core.auroc_identity.min():.3f} to {core.auroc_identity.max():.3f}; "
                f"above the shuffle null in both directions at p < 0.05 in {int((core.above_null_both_p05 == 'yes').sum())} of {len(core)} (report: 30 of 30, 0.736 to 0.995; the report range is over the two directions separately, this one over their mean). "
                f"Stanley adds 30 study pairs x 2 subject sets ({int(o.source.str.startswith('Stanley').sum())} rows).")
    pooled = pool.dx_auroc
    F.note("g", f"diagnosis AUROC: pooled per source {pooled.min():.2f} to {pooled.max():.2f} (legend '0.48 to 0.58': match). Per relation (Pittsburgh) it runs "
                f"{core.auroc_dx.min():.2f} to {core.auroc_dx.max():.2f}; the legend range is the pooled one.")
    F.note("g", "legend count: '30 of 30 relations' = LIBD regions 1 + PsychAD x LIBD 1 + Pittsburgh 28. Fig 5c legend counts '58 relations (Pittsburgh 28, Stanley 30)', which leaves out "
                "the LIBD and PsychAD ones; the report counts 30 + 30 Stanley pairs = 60 relations.")


# ================================================================== Figure 5
def _bulk_sep_map():
    sr = rd(P_SUBREL)
    m = {}
    for r in sr.itertuples():
        m[(r.dataset_a, r.dataset_b)] = r.relation_group
        m[(r.dataset_b, r.dataset_a)] = r.relation_group
    return m


def fig5_a(F):
    """a: 192 subject-separate bulk pairs, dx / age / sex."""
    m = _bulk_sep_map()
    p, pm = rd(P_CAS_PAIRS), rd(P_CAS_MATCHED)
    for d in (p, pm):
        d["relation_now"] = [m.get((a, b)) for a, b in zip(d.dataset_a, d.dataset_b)]
    sel = lambda d: d[(d.relation_now == "subject-separate") & (d.compartment_a == "bulk") & (d.compartment_b == "bulk")]
    p, pm = sel(p), sel(pm)
    spec = {"dx": ("sz_auto", "sz_m_auto"), "age": ("age_auto", "age_m_auto"), "sex": ("sex_auto", "sex_m_auto")}
    rows = []
    for trait, (c, cm) in spec.items():
        a = p[p.contrast == c].set_index(["dataset_a", "dataset_b"])
        b = pm[pm.contrast == cm].set_index(["dataset_a", "dataset_b"])
        for k, r in a.iterrows():
            q = b.loc[k] if k in b.index else None
            f200 = r.fold200
            rows.append(dict(pair_id=f"{k[0]}_{k[1]}", trait=trait, rho_t=r.rho_t,
                             rho_t_matched=None if q is None else q.rho_t_mean,
                             top200_fold=f200, top200_log2fold=(math.log2(f200) if f200 and f200 > 0 else None),
                             top200_fold_matched=None if q is None else q.fold200_mean,
                             n_a=r.n_a, n_b=r.n_b, n_genes_shared=r.n_genes_shared,
                             n_matched_a=None if q is None else q.n_matched_a, n_matched_b=None if q is None else q.n_matched_b))
    d = pd.DataFrame(rows)
    F.write("fig05_a_pairs_dx_age_sex.tsv", "a", d, [P_CAS_PAIRS, P_CAS_MATCHED, P_CAS_SUM, P_SUBREL],
            ["192 subject-separate bulk pairs of 22 dataset-regions; relation read from the current subject_relations.tsv (relation_group), not the stale donor_relation column of pairs.tsv; "
             "autosomal genes; dx = sz_auto, age = age_auto (continuous, all subjects), sex = sex_auto; rho_t_matched = mean over 50 matched-size draws (pairs_matched.tsv rho_t_mean, "
             "contrasts sz_m_auto, age_m_auto, sex_m_auto); top200_log2fold = log2 of the full-data top-200 fold over chance, blank when the overlap is 0 (fold 0); "
             "top200_fold_matched = fold200_mean"])
    s = d.groupby("trait").agg(n=("rho_t", "size"), rho=("rho_t", "median"), rho_m=("rho_t_matched", "median"), f=("top200_fold", "median"))
    F.note("a", "medians over " + "; ".join(f"{k}: n {int(v.n)}, rho(t) {v.rho:.3f}, matched {v.rho_m:.3f}, top-200 fold {v.f:.2f}" for k, v in s.iterrows())
           + ". Legend: n = 192, age 0.143, diagnosis 0.108, sex 0.002: match.")
    nz = int((d.top200_fold == 0).sum())
    F.note("a", f"{nz} pair x trait rows have top-200 fold 0 (no overlap), so log2 fold is blank; the plot needs a rule for them.")
    F.note("a", "analyses/sex_genes holds the XIST-list audit (xist_lists.tsv: which SZ lists contain XIST), not pair agreement; not used for a. Sex agreement is controls_age_sex sex_auto.")
    F.note("a", "controls_age_sex results are from 2026-09-25 (summary.tsv) with the relation table of that date; pair selection here uses the 2026-10-08 table and still gives 192 per contrast.")


def fig5_b(F):
    ds = rd(P_CAE_DS)
    F.write("fig05_b_age_external_per_dataset.tsv", "b", ds, [P_CAE_DS],
            ["per dataset, the median over its subject-separate partner pairs of rho_t and top-200 fold, by class: age_age (SZ-study age vs SZ-study age), age_ext (SZ-study age vs "
             "external age), sz_sz (SZ vs SZ), sz_ext, ext_ext; set SZ-study (22 bulk dataset-regions) or external (9 frontal-cortex datasets)"])
    sm = rd(P_CAE_SUM)
    F.write("fig05_b_age_external_summary.tsv", "b", sm, [P_CAE_SUM], ["medians over pairs by class and subject relation (summary.tsv of controls_age_external)"])
    comb = rd(P_CAE_COMB)
    F.write("fig05_b_age_external_combined.tsv", "b", comb, [P_CAE_COMB], ["combined-signature comparisons (combined.tsv)"])
    ss = sm[sm.relation_group == "subject-separate"].set_index("class")
    F.note("b", f"subject-separate medians: age vs external age {ss.rho_t_median['age_ext']:.3f} ({int(ss.n_pairs['age_ext'])} pairs), age vs age {ss.rho_t_median['age_age']:.3f} "
                f"({int(ss.n_pairs['age_age'])}), SZ vs SZ {ss.rho_t_median['sz_sz']:.3f}. Legend: 0.21, 0.14, 0.11: match. SZ-study rows in per_dataset: "
                f"{int((ds.set == 'SZ-study').sum())}; external: {int((ds.set == 'external').sum())} (legend 22 and 9: match).")
    F.note("b", "controls_age_external was built 2026-09-26 (summary.tsv, per_dataset.tsv); its inputs are not tied to the 10-08 shas; no later rebuild found. It also does not use the 64-study list population.")


def fig5_c(F):
    # negative control: the fig02 panel c pairwise folds (copy; 'Top left' of the legend)
    d = pd.read_csv(P_FIG02_C, sep="\t", comment="#")
    F.write("fig05_c_negative_control_pairwise_fold.tsv", "c", d, [P_FIG02_C],
            ["copy of figure_data/fig02/fig02_c_pairwise_fold_by_era.tsv (all rows); fold of every pair of published brain lists over its own hypergeometric null, old-old, new-new, old-new"])
    d2 = pd.read_csv(P_FIG02_CP, sep="\t", comment="#")
    F.write("fig05_c_negative_control_pooled.tsv", "c", d2, [P_FIG02_CP], ["copy of figure_data/fig02/fig02_c_pooled_by_era.tsv"])
    # per-dataset agreement of the SZ effect, age within SZ, age external: same as b
    ds = rd(P_CAE_DS)
    F.write("fig05_c_per_dataset_agreement.tsv", "c", ds, [P_CAE_DS], ["same table as fig05_b_age_external_per_dataset.tsv (bottom left of c)"])
    # identification
    idr = _ident_rows()
    o = idr[["source", "relation", "subjects", "n_candidates", "auroc_identity", "auroc_ab", "auroc_ba", "null_ab", "null_ba", "p_ab", "p_ba", "above_null_both_p05"]]
    F.write("fig05_c_identification_auroc.tsv", "c", o, [P_G, P_H, P_E_PSY, P_E_LIBD],
            ["mean identification AUROC per relation (a to b, b to a, and their mean), residuals without diagnosis; Stanley rows exist for two subject sets (all_dx, sz_ctrl)"])
    st = o[o.source.str.startswith("Stanley")]
    pit = o[o.source == "Pittsburgh"]
    F.note("c", f"identification: Pittsburgh {len(pit)} relations, Stanley {st.relation.nunique()} study pairs ({len(st)} rows over two subject sets), plus LIBD regions 1 and PsychAD x LIBD 1. "
                "Legend 'mean identification AUROC for 58 cross-dataset donor relations (Pittsburgh 28, Stanley 30)': counts match Pittsburgh and Stanley; the two LIBD/PsychAD relations are extra rows, drop them to get 58.")
    F.note("c", "no T1D positive control (removed from the paper by decision); the repo Dev/t1d-islet-positive-control was not read.")
    F.note("c", "negative control copied from fig02 (see fig02/STATUS.md for its checks); the file notes its own input shas.")


def main():
    F3 = Figure("fig03")
    for fn in (fig3_acd, fig3_b, fig3_e, fig3_f, fig3_g):
        fn(F3)
    F3.finish()
    F5 = Figure("fig05")
    for fn in (fig5_a, fig5_b, fig5_c):
        fn(F5)
    F5.finish()


if __name__ == "__main__":
    main()
