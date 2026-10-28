#!/usr/bin/env python3
"""Summary tables for split_halves/: per-dataset split-half curves, positive-control ceilings, cross-dataset pairs at matched
n, and where the published different-cohort pairs fall.

    opt/anaconda3/bin/python3 scripts/summarize.py        (run from split_halves/)

Adopted, not recomputed (read only):
  ../libd_brainseq/analysis/results/a_reps.tsv.gz   BrainSeq Phase 2 DLPFC and hippocampus split halves (SZ), and DLPFC vs
                                                    hippocampus on the same or subject-separate subjects (a cross pair here)
  ../libd_brainseq/analysis/results/d_reps.tsv.gz   BrainSeq Phase 2 positive controls (sz, age, male)
  ../fairness_sample_level/results/a_reps.tsv.gz    Gemma split halves (SZ; kind 'within') for the Gemma datasets in data/,
                                                    and subject-separate Gemma dataset pairs (kind 'cross', relation
                                                    'subject-separate'), same design (fit_sz = ~ sz + age + sex)
Computed here (scripts/split_halves.R): results/within_sz/, within_controls/, cross_sz/, cross_controls/.
Published pairs: ../h1_pairs.tsv (scenario 'published, brain', cohort 'different'), N from ../h2/h2_studies.tsv, placed at
n = min(N_a, N_b) / 2 and the nearest level, as ../libd_brainseq/analysis/scripts/summarize.py section_a (its helpers
reps_summary, med, q, pct are imported from there).

Outputs (results/):
  datasets.tsv            one row per dataset: family, n, largest n per half, where its SZ split halves come from
  within_summary.tsv      dataset x level: medians (and quartiles) over repeats of the agree() metrics, SZ split halves
  controls_summary.tsv    dataset x contrast x gene set x level: medians over repeats (d_controls design)
  cross_summary.tsv       pair x kind x contrast x gene set x level: medians over repeats
  headline.tsv            dataset: SZ split-half rho (log2 FC and t) and top-100 fold at 10, 25, 50 per group; age and sex
                          (autosomal) rho of t at 10 and 20 per group beside SZ in the same design
  published_pairs.tsv     per published pair: level, percentile of its log2 FC Spearman and top-N fold among each
                          dataset's split-half repeats at that level, and among the pooled different-cohort pair repeats
  published_summary.tsv   per level: medians over pairs of the above
"""
import glob
import hashlib
import importlib.util
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SH = os.path.normpath(os.path.join(HERE, ".."))
GL = os.path.normpath(os.path.join(SH, ".."))
LIBD = os.path.join(GL, "libd_brainseq", "analysis")
FSL = os.path.join(GL, "fairness_sample_level")
R = lambda *p: os.path.join(SH, "results", *p)
spec = importlib.util.spec_from_file_location("libd_summarize", os.path.join(LIBD, "scripts", "summarize.py"))
LS = importlib.util.module_from_spec(spec); spec.loader.exec_module(LS)
med, q, pct, reps_summary, KIND_MAP = LS.med, LS.q, LS.pct, LS.reps_summary, LS.KIND_MAP
LEVELS = [5, 8, 10, 12, 15, 20, 25, 30, 50, 75]
NEW_LIBD = ["LIBD_Phase1_DLPFC", "LIBD_Caudate", "LIBD_Habenula"]
P2 = {"DLPFC": "LIBD_Phase2_DLPFC", "HIPPO": "LIBD_Phase2_HIPPO"}
INPUTS = {}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def rd(p, **kw):
    INPUTS[os.path.relpath(p, SH)] = sha(p)
    return pd.read_csv(p, sep="\t", comment="#", low_memory=False, **kw)


def write(df, name, header):
    with open(R(name), "w") as fh:
        fh.write(f"# built {datetime.now():%Y-%m-%d %H:%M} by split_halves/scripts/summarize.py\n")
        for h in header:
            fh.write(f"# {h}\n")
        df.to_csv(fh, sep="\t", index=False, float_format="%.4g")


def cat_dir(d):
    fs = sorted(glob.glob(R(d, "*.tsv.gz")))
    return pd.concat([rd(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def load_within():
    """SZ split-half repeats of every dataset, one table: dataset, level, n_per_group, rep, agree() columns."""
    parts = []
    W = cat_dir("within_sz")
    if len(W):
        W["source"] = "split_halves/results/within_sz"
        parts.append(W.rename(columns={"dataset_a": "dataset"}))
    A = rd(os.path.join(LIBD, "results", "a_reps.tsv.gz"))
    A["kind"] = A.kind.replace(KIND_MAP)
    a = A[A.kind == "split_halves"].copy()
    a["dataset"] = a.region.map(P2); a["source"] = "libd_brainseq/analysis/results/a_reps.tsv.gz"
    parts.append(a)
    prep = rd(R("prep_gemma.tsv"))
    gem = set(prep[prep.used.astype(str) == "True"].dataset)
    F = rd(os.path.join(FSL, "results", "a_reps.tsv.gz"))
    f = F[(F.kind == "within") & F.dataset_a.isin(gem)].copy()
    f["dataset"] = f.dataset_a; f["source"] = "fairness_sample_level/results/a_reps.tsv.gz"
    parts.append(f)
    out = pd.concat(parts, ignore_index=True)
    out["level"] = out.level.astype(str)
    return out, F, A


def family(d):
    return "LIBD" if d.startswith("LIBD") else "Gemma"


def main():
    W, F, A = load_within()
    S = reps_summary(W, ["dataset", "level"])
    S["family"] = S.dataset.map(family)
    write(S, "within_summary.tsv", ["SZ split halves (two disjoint halves of n SZ + n controls of one dataset), medians over repeats (quartiles for rho_t, top100_fold, top100_sign_in_other)",
                                    "sources: results/within_sz (new LIBD datasets), ../libd_brainseq/analysis/results/a_reps.tsv.gz (Phase 2), ../fairness_sample_level/results/a_reps.tsv.gz (Gemma, kind within)"])

    # datasets
    pl = rd(R("prep_libd.tsv")); pg = rd(R("prep_gemma.tsv"))
    pg = pg[pg.used.astype(str) == "True"]
    rows = []
    for r in pl.itertuples():
        rows.append(dict(dataset=r.dataset, family="LIBD", technology="RNA-seq", n_sz=r.n_sz, n_ctl=r.n_ctl, n_male=r.n_male, n_female=r.n_female,
                         genes=r.genes_unique_symbol, covariates=r.covariates, model="limma-trend", cluster="libd",
                         sz_reps_source=W[W.dataset == r.dataset].source.iloc[0] if (W.dataset == r.dataset).any() else ""))
    for r in pg.itertuples():
        rows.append(dict(dataset=r.dataset, family="Gemma", technology="", n_sz=r.n_sz, n_ctl=r.n_ctl, n_male=r.n_male, n_female=r.n_female,
                         genes=r.n_genes, covariates="age + male", model="limma (eBayes, no trend)", cluster=r.cluster,
                         sz_reps_source=W[W.dataset == r.dataset].source.iloc[0] if (W.dataset == r.dataset).any() else ""))
    DS = pd.DataFrame(rows)
    DS["max_n_per_half"] = (np.minimum(DS.n_sz, DS.n_ctl) // 2).astype(int)
    write(DS, "datasets.tsv", ["datasets with at least 20 SZ and 20 controls after loading (split halves of >= 10 per group); results/prep_libd.tsv, prep_gemma.tsv",
                               "Gemma technology: see ../controls_age_sex/data/datasets.tsv and the Gemma record of each series"])

    # positive controls
    C = cat_dir("within_controls")
    C = C.rename(columns={"dataset_a": "dataset"})
    D = rd(os.path.join(LIBD, "results", "d_reps.tsv.gz"))
    D["kind"] = D.kind.replace(KIND_MAP)
    d = D[D.kind == "split_halves"].copy(); d["dataset"] = d.region.map(P2)
    C = pd.concat([C, d], ignore_index=True)
    C["level"] = C.level.astype(str)
    CS = reps_summary(C, ["dataset", "contrast", "gene_set", "level"])
    CS["family"] = CS.dataset.map(family)
    write(CS, "controls_summary.tsv", ["positive controls, split halves of one dataset (d_controls.R design): sz = n SZ + n controls; male = n males + n females (~ male + sz + covariates); age = 2n subjects at random, age continuous (~ age + sz + covariates)",
                                       "gene_set auto = autosomal genes, all = every gene; sources: results/within_controls, ../libd_brainseq/analysis/results/d_reps.tsv.gz (Phase 2)"])

    # cross-dataset pairs
    xs_, xc_ = cat_dir("cross_sz"), cat_dir("cross_controls")
    xs_["part"] = "sz"; xc_["part"] = "controls"
    X = pd.concat([xs_, xc_], ignore_index=True)
    X["gene_set"] = X.gene_set.fillna("all") if "gene_set" in X else "all"
    a = A[A.kind.isin(["same_subjects_diff_region", "separate_subjects_diff_region"])].copy()
    a["kind"] = a.kind.map({"same_subjects_diff_region": "same_subjects", "separate_subjects_diff_region": "separate_subjects"})
    a["pair"] = "P2_DLPFC_HIPPO"; a["contrast"] = "sz"; a["gene_set"] = "all"; a["part"] = "sz"
    d = D[D.kind.isin(["same_subjects_diff_region", "separate_subjects_diff_region"])].copy()
    d["kind"] = d.kind.map({"same_subjects_diff_region": "same_subjects", "separate_subjects_diff_region": "separate_subjects"}); d["pair"] = "P2_DLPFC_HIPPO"; d["part"] = "controls"
    g = F[(F.kind == "cross") & (F.relation == "subject-separate") & F.dataset_a.isin(set(pg.dataset)) & F.dataset_b.isin(set(pg.dataset))].copy()
    g["pair"] = g.dataset_a + "_x_" + g.dataset_b; g["kind"] = "separate_cohorts_gemma"; g["contrast"] = "sz"; g["gene_set"] = "all"; g["part"] = "sz"
    X = pd.concat([X, a, d, g], ignore_index=True)
    X["level"] = X.level.astype(str)
    XS = reps_summary(X, ["part", "pair", "kind", "contrast", "gene_set", "level"])
    write(XS, "cross_summary.tsv", ["two datasets at the same n per group; kinds: same_subjects (same LIBD subjects by BrNum in both datasets), separate_subjects (subject-separate halves of the LIBD subjects in both), separate_cohorts (independent draws from two datasets of different brain banks, each with its own model), separate_cohorts_gemma (fairness_sample_level Gemma pairs with relation subject-separate)",
                                    "part sz = the SZ design of a_split_half.R (levels 5 to 75, 100 repeats); part controls = the d_controls.R design (sz, age, male at 10, 20, 30, 50; 50 repeats; gene_set auto and all); sources: results/cross_sz, cross_controls; ../libd_brainseq/analysis/results/a_reps.tsv.gz and d_reps.tsv.gz (P2_DLPFC_HIPPO); ../fairness_sample_level/results/a_reps.tsv.gz (Gemma pairs)"])

    # headline per dataset
    rows = []
    for ds, g in S.groupby("dataset"):
        d = dict(dataset=ds, family=family(ds))
        for lv in ["10", "25", "50", "max"]:
            s = g[g.level == lv]
            d[f"n_{lv}"] = float(s.n_per_group.iloc[0]) if len(s) else np.nan
            for m in ["rho_logfc", "rho_t", "top100_fold", "pi1", "pi1_null", "pi1_excess"]:
                d[f"{m}_{lv}"] = float(s[m].iloc[0]) if len(s) and m in s else np.nan
        for lv in ["10", "20", "30"]:
            for ct, gs in [("sz", "auto"), ("age", "auto"), ("male", "auto")]:
                c = CS[(CS.dataset == ds) & (CS.contrast == ct) & (CS.gene_set == gs) & (CS.level == lv)]
                d[f"ctl_{ct}_rho_t_{lv}"] = float(c.rho_t.iloc[0]) if len(c) else np.nan
                d[f"ctl_{ct}_rho_logfc_{lv}"] = float(c.rho_logfc.iloc[0]) if len(c) else np.nan
        rows.append(d)
    H = pd.DataFrame(rows).merge(DS[["dataset", "n_sz", "n_ctl", "cluster"]], on="dataset", how="left")
    write(H, "headline.tsv", ["from within_summary.tsv (SZ split halves; *_max = at each dataset's largest n per half) and controls_summary.tsv (ctl_* = same d_controls design, autosomal genes; male = sex)"])

    # published pairs
    h = rd(os.path.join(GL, "h1_pairs.tsv"), dtype=str)
    h = h[(h.scenario == "published, brain") & (h.cohort == "different")].copy()
    st = rd(os.path.join(GL, "h2", "h2_studies.tsv"), dtype=str)
    N = dict(zip(st.pmid, pd.to_numeric(st.n, errors="coerce")))
    h["N_a"] = h.study_a.map(N); h["N_b"] = h.study_b.map(N)
    for c in ["top_n", "top_fold", "rho_any"]:
        h[c] = pd.to_numeric(h[c], errors="coerce")
    h["n_eq"] = np.minimum(h.N_a, h.N_b) / 2
    n_all = len(h); h = h[h.n_eq.notna()].copy()
    h["level"] = h.n_eq.apply(lambda v: str(min(LEVELS, key=lambda L: abs(L - v))))
    h["topN_match"] = h.top_n.apply(lambda v: np.nan if pd.isna(v) else min([20, 50, 100], key=lambda n: abs(np.log(n) - np.log(max(v, 1)))))
    Wl = {k: g for k, g in W.groupby(["dataset", "level"])}
    Xc = X[(X.part == "sz") & (X.contrast == "sz") & (X.gene_set == "all") & X.kind.isin(["separate_cohorts", "separate_cohorts_gemma"])]
    Xl = {k: g for k, g in Xc.groupby("level")}
    dss = sorted(W.dataset.unique())
    rows = []
    for r in h.itertuples():
        fm = f"top{int(r.topN_match)}_fold" if not pd.isna(r.topN_match) else None
        d = dict(study_a=r.study_a, label_a=r.label_a, N_a=r.N_a, study_b=r.study_b, label_b=r.label_b, N_b=r.N_b, n_eq=r.n_eq, level=r.level,
                 top_n=r.top_n, top_fold=r.top_fold, rho_logfc=r.rho_any)
        prho, pfold = [], []
        for ds in dss:
            ref = Wl.get((ds, r.level))
            if ref is None or len(ref) < 20:
                continue
            prho.append(pct(ref.rho_logfc, r.rho_any)); pfold.append(pct(ref[fm], r.top_fold) if fm else np.nan)
        d.update(n_datasets=len(prho), rho_pct_median_over_datasets=med(pd.Series(prho)), rho_pct_min_over_datasets=np.nanmin(prho) if len(prho) and not all(pd.isna(prho)) else np.nan,
                 rho_pct_max_over_datasets=np.nanmax(prho) if len(prho) and not all(pd.isna(prho)) else np.nan,
                 fold_pct_median_over_datasets=med(pd.Series(pfold)))
        xr = Xl.get(r.level)
        d.update(n_cross_reps=0 if xr is None else len(xr), rho_pct_in_cross=pct(xr.rho_logfc, r.rho_any) if xr is not None else np.nan,
                 cross_rho_logfc_median=med(xr.rho_logfc) if xr is not None else np.nan,
                 fold_pct_in_cross=pct(xr[fm], r.top_fold) if (xr is not None and fm) else np.nan)
        rows.append(d)
    P = pd.DataFrame(rows)
    write(P, "published_pairs.tsv", [f"../h1_pairs.tsv scenario 'published, brain', cohort 'different' ({n_all} pairs, {len(P)} with N for both studies in ../h2/h2_studies.tsv); n_eq = min(N_a, N_b) / 2, nearest level",
                                     "rho_pct_*_over_datasets: over the datasets with >= 20 SZ split-half repeats at that level, the fraction of that dataset's repeats with log2 FC Spearman <= the published pair's (rho_any); median, min and max over datasets",
                                     "*_in_cross: the same among the pooled different-cohort pair repeats at that level (results/cross_sz separate_cohorts and the fairness_sample_level Gemma subject-separate pairs); top fold at the nearest of top 20 / 50 / 100"])
    rows = []
    for lv in [str(x) for x in LEVELS] + ["all"]:
        m = P if lv == "all" else P[P.level == lv]
        if not len(m):
            continue
        rr = m[m.rho_logfc.notna()]
        rows.append(dict(level=lv, n_pairs=len(m), n_rho=len(rr), pub_rho_logfc_median=med(rr.rho_logfc),
                         n_datasets_median=med(m.n_datasets), rho_pct_median_over_datasets=med(rr.rho_pct_median_over_datasets),
                         n_rho_below_q05_every_dataset=int((rr.rho_pct_max_over_datasets < 0.05).sum()),
                         cross_rho_logfc_median=med(m.cross_rho_logfc_median), rho_pct_in_cross_median=med(rr.rho_pct_in_cross),
                         n_rho_below_cross_q05=int((rr.rho_pct_in_cross < 0.05).sum()),
                         pub_top_fold_median=med(m.top_fold), fold_pct_median_over_datasets=med(m.fold_pct_median_over_datasets), fold_pct_in_cross_median=med(m.fold_pct_in_cross)))
    write(pd.DataFrame(rows), "published_summary.tsv", ["per level: medians over published pairs of results/published_pairs.tsv; n_rho_below_q05_every_dataset = pairs below the 5th percentile of every dataset's split halves at that level"])
    pd.Series(INPUTS).to_json(R("summarize_inputs.json"), indent=1)
    print(H.to_string())


if __name__ == "__main__":
    main()
