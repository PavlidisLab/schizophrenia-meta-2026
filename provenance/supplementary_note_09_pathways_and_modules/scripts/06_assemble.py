#!/usr/bin/env python3
"""Assembles the bulk TSVs from the per-dataset outputs of 03_bulk_main.py / 08_neuron_astro_score_corr.py (copied to inputs/bulk_out/), pools across datasets with no shared donors, and fits the schizophrenia-effect models with the
estimated agonal score in bulk.  Run locally (one core, a few minutes).
Pooling: datasets are grouped by the registry's cohort cluster (analyses/age_and_sex_control_comparisons/data/datasets.tsv; datasets of one cluster share donors); the pooled value = mean of the cluster means (equal weight per cluster); interval = bootstrap over clusters (2,000).
Median row: lo95/hi95 hold the minimum and maximum over datasets.  GTEx regions share donors with each other and are not pooled.
Bulk schizophrenia-effect models: OLS of the score/SD on sz + age + sex; + estimated agonal (SD); + pH (SD, where recorded); donors with every variable of the richest model; bootstrap stratified by diagnosis (1,000); change = model minus base, same draws."""
import os, sys, glob, numpy as np, pandas as pd
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(H, "../.."))
sys.path.insert(0, H); from snaplib import write_tsv, pct, pool_clusters
rng = np.random.default_rng(20261012); BO = f"{H}/inputs/bulk_out"
def cat(pat):
    fs = sorted(glob.glob(f"{BO}/{pat}")); return pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True), fs
def with_pools(df, keys, cols, extra=None):
    """append pooled and median rows for registry bulk datasets (not GTEx)"""
    out = [df]; reg = df[~df.dataset.str.startswith("GTEx")]
    for k, g in reg.groupby(keys, dropna=False):
        if len(g) < 2: continue
        row = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
        p, lo, hi, nd, nc = pool_clusters(g[cols[0]].values, g.cluster.values, rng)
        r1 = dict(row, effect=p, lo95=lo, hi95=hi, dataset=f"pooled over {nc} cohort clusters ({nd} datasets)", cluster="pooled")
        for c2 in cols[1:]:
            if c2 in g:
                p2, l2, h2, *_ = pool_clusters(g[c2].values, g.cluster.values, rng); r1[c2] = p2; r1[c2.replace("rel", "rel") + "_lo95"] = l2; r1[c2 + "_hi95"] = h2
        if "n" in g: r1["n"] = int(g.n.sum())
        r2 = dict(row, effect=g[cols[0]].median(), lo95=g[cols[0]].min(), hi95=g[cols[0]].max(), dataset=f"median over {len(g)} datasets (lo95/hi95 = minimum/maximum)", cluster="median")
        if "n" in g: r2["n"] = int(g.n.sum())
        out += [pd.DataFrame([r1, r2])]
    return pd.concat(out, ignore_index=True)
DATASETS = f"{REPO}/analyses/age_and_sex_control_comparisons/data/datasets.tsv"
# ---------- tests
T, fs = cat("tests_*.tsv")
T["tmp"] = T.rel
T = with_pools(T.drop(columns=["tmp"]), ["test", "group", "score"], ["effect", "rel"])
# fix interval columns for rel on pooled rows
T = T.rename(columns={"rel_lo95": "rel_lo95", "rel_hi95": "rel_hi95"})
tcols = ["effect", "lo95", "hi95", "dataset", "cluster", "test", "group", "score", "n", "null_mean", "null_lo95", "null_hi95", "rel", "rel_lo95", "rel_hi95"]
T = T[[c for c in tcols if c in T.columns] + [c for c in T.columns if c not in tcols and c != "p_value_tsv_only"] + ["p_value_tsv_only"]]
T["note"] = "effect: SD of score (sz, age per decade, hardy4_vs_2, ventilator, male), SD per SD (slope_on_est_agonal), partial Spearman rho (partial_spearman_*); rel = effect minus the mean of the decile-matched null scores; null_lo95/null_hi95 = 2.5 and 97.5 percentiles of the null effect"
write_tsv(T[~T.test.str.startswith("slope_on_est_agonal")], f"{H}/snap_bulk_tests.tsv", "SNAP score in 22 bulk registry datasets and GTEx cortex: schizophrenia effect, age slope, recorded covariates, against decile-matched random weights (10,000 draws)", "analyses/snap_signal_analysis/03_bulk_main.py, 06_assemble.py", fs[:2], last="p_value_tsv_only")
G = T[T.dataset.str.startswith("GTEx") & ~T.test.str.startswith("slope_on_est_agonal")]
write_tsv(G, f"{H}/snap_gtex_death_variables.tsv", "SNAP-weighted score in GTEx cortex (Cortex, BA9, BA24) against Hardy scale, ventilator death, ischemic time, RIN, age and sex", "analyses/snap_signal_analysis/03_bulk_main.py", fs[:2])
# ---------- coherence
C, fc = cat("coh_*.tsv"); NA, fn = cat("nacorr_*.tsv"); C = pd.concat([C, NA], ignore_index=True)
C.loc[C.statistic.str.startswith("variance_share") | C.set.str.startswith("score_corr"), "top_n"] = np.nan   # number of weighted genes differs by dataset; not a grouping key
C = with_pools(C, ["set", "top_n", "statistic", "variant", "aligned"], ["effect"])
ccols = ["effect", "lo95", "hi95", "dataset", "cluster", "set", "top_n", "statistic", "variant", "aligned", "observed", "matched_mean", "n_samples"]
C = C[[c for c in ccols if c in C.columns] + [c for c in C.columns if c not in ccols and c != "p_value_tsv_only"] + ["p_value_tsv_only"]]
write_tsv(C, f"{H}/snap_bulk_coherence.tsv", "Coexpression of SNAP's top-loading genes in bulk against expression-decile-matched genes (10,000 draws): mean pairwise correlation (all / neuron-loaded / astrocyte-loaded / neuron-astrocyte cross), variance share of the SNAP direction, neuron-astrocyte score correlation", "analyses/snap_signal_analysis/03_bulk_main.py, 08_neuron_astro_score_corr.py, 06_assemble.py", fc[:2] + fn[:1])
# ---------- PCs
P, fp = cat("pcs_*.tsv")
sm = []
for (sc, pc), g in P[~P.dataset.str.startswith("GTEx")].groupby(["score", "pc"]):
    sm.append(dict(effect=g.abs_rho.median(), lo95=g.abs_rho.min(), hi95=g.abs_rho.max(), dataset=f"median |rho| over {len(g)} registry datasets (lo95/hi95 = minimum/maximum)", cluster="median", score=sc, pc=pc, n=int(g.n.sum()), null_mean_abs_rho=g.null_mean_abs_rho.median() if sc == "snap_total" else np.nan))
P = pd.concat([P, pd.DataFrame(sm)], ignore_index=True); P["note"] = "effect: signed Spearman rho of the score with the PC (centred expression, no covariates removed; PC signs are arbitrary); median rows: median |rho|"
write_tsv(P, f"{H}/snap_vs_pcs.tsv", "Spearman correlation of the bulk SNAP score and of the estimated agonal score with PC1 to PC5 of each dataset", "analyses/snap_signal_analysis/03_bulk_main.py, 06_assemble.py", fp[:2])
# ---------- per-sample scores
fsamp = sorted(glob.glob(f"{BO}/persample_*.tsv")); PS = []
for f in fsamp:
    d = pd.read_csv(f, sep="\t", index_col=0); d.index.name = "sample"; d = d.reset_index(); d.insert(0, "dataset", os.path.basename(f)[10:-4]); PS.append(d[["dataset", "sample"] + [c for c in ["sz", "age", "male", "pH", "PMI", "RIN", "agonal", "hardy", "isch_h", "snap_total", "snap_neuron", "snap_astrocyte", "est_agonal_ex1000", "est_agonal_ex200", "est_agonal_none"] if c in d.columns]])
PSA = pd.concat(PS, ignore_index=True)
write_tsv(PSA, f"{H}/snap_bulk_scores.tsv", "SNAP score per sample in the bulk datasets (z-scored within dataset; loading-weighted mean of standardized expression; total = neuron + astrocyte loadings), with the estimated agonal score and recorded covariates", "analyses/snap_signal_analysis/03_bulk_main.py", fsamp[:2], cols_first=("dataset", "sample"), last="est_agonal_none")
# ---------- estimated agonal: slopes and schizophrenia-effect models
rowsA = []
sl = T[T.test.str.startswith("slope_on_est_agonal") & T.score.isin(["total"])].copy()
for _, r in sl.iterrows():
    rowsA.append(dict(effect=r.effect, lo95=r.lo95, hi95=r.hi95, dataset=r.dataset, cohort_group="bulk" if not r.dataset.startswith("GTEx") else "GTEx cortex", model="SNAP total (SD) on estimated agonal score (SD), " + r.test.replace("slope_on_est_agonal_", "") + " (ex1000 / ex200 = SNAP top 1,000 / 200 genes removed from the axis); given age and sex", analysis="slope", group=r.group, n=r.n, null_mean=r.null_mean, rel=r.rel, rel_lo95=r.rel_lo95, rel_hi95=r.rel_hi95, p_value_tsv_only=r.p_value_tsv_only))
sn_sl = pd.read_pickle(f"{H}/_sl.pkl"); sn_sl["analysis"] = "slope"
for _, r in sn_sl.iterrows(): rowsA.append(dict(r))
reg = pd.read_csv(DATASETS, sep="\t", comment="#").set_index("dataset").cluster
res = []
for d, g in PSA.groupby("dataset"):
    if d.startswith("GTEx") or "sz" not in g or g.sz.nunique() < 2 or min(g.sz.sum(), (1 - g.sz).sum()) < 5: continue
    g = g.dropna(subset=["age", "male", "est_agonal_ex1000", "est_agonal_ex200"]).reset_index(drop=True); haspH = "pH" in g and g.pH.notna().sum() >= 20
    if haspH: g = g.dropna(subset=["pH"]).reset_index(drop=True)
    if g.sz.nunique() < 2 or min(g.sz.sum(), (1 - g.sz).sum()) < 5: continue
    models = [("base: age, sex", []), ("+ estimated agonal (top 1,000 removed)", ["est_agonal_ex1000"]), ("+ estimated agonal (top 200 removed)", ["est_agonal_ex200"])] + ([("+ pH", ["pH"]), ("+ estimated agonal (top 1,000 removed) + pH", ["est_agonal_ex1000", "pH"])] if haspH else [])
    def fit(dd):
        y = dd.snap_total.values / dd.snap_total.std(ddof=1); out = []
        for _, ex in models:
            X = np.column_stack([np.ones(len(dd)), dd.sz.values.astype(float), dd.age.values, dd.male.values.astype(float)] + [(dd[c].values - dd[c].mean()) / dd[c].std(ddof=1) for c in ex]); out.append(np.linalg.lstsq(X, y, rcond=None)[0][1])
        return np.array(out)
    est = fit(g); bs = []
    for _ in range(1000):
        ii = np.concatenate([np.where(g.sz.values == k)[0][rng.integers(0, (g.sz.values == k).sum(), (g.sz.values == k).sum())] for k in (0, 1)]); bs.append(fit(g.iloc[ii]))
    bs = np.array(bs)
    for j, (mn, _) in enumerate(models):
        lo, hi = pct(bs[:, j]); dlo, dhi = pct(bs[:, j] - bs[:, 0]); res.append(dict(effect=est[j], lo95=lo, hi95=hi, dataset=d, cluster=reg[d], cohort_group="bulk", model=mn, analysis="sz_effect", change_from_base=est[j] - est[0], change_lo95=dlo, change_hi95=dhi, n=len(g), n_sz=int(g.sz.sum())))
R = pd.DataFrame(res); pooled = []
for mn, g in R.groupby("model"):
    if g.dataset.nunique() < 2: continue
    p, lo, hi, nd, nc = pool_clusters(g.effect.values, g.cluster.values, rng); p2, l2, h2, *_ = pool_clusters(g.change_from_base.values, g.cluster.values, rng)
    pooled.append(dict(effect=p, lo95=lo, hi95=hi, dataset=f"pooled over {nc} cohort clusters ({nd} datasets)", cohort_group="bulk", model=mn, analysis="sz_effect", change_from_base=p2, change_lo95=l2, change_hi95=h2, n=int(g.n.sum()), n_sz=int(g.n_sz.sum())))
A = pd.DataFrame(rowsA + [dict(r) for _, r in R.drop(columns=["cluster"]).iterrows()] + pooled)
adj = pd.read_pickle(f"{H}/_adj.pkl"); adj["analysis"] = "sz_effect"; A = pd.concat([A, adj], ignore_index=True)
acols = ["effect", "lo95", "hi95", "dataset", "cohort_group", "model", "analysis", "group", "n", "n_sz", "change_from_base", "change_lo95", "change_hi95", "null_mean", "rel", "rel_lo95", "rel_hi95"]
A = A[[c for c in acols if c in A.columns] + [c for c in A.columns if c not in acols and c != "p_value_tsv_only"] + ["p_value_tsv_only"]]
write_tsv(A, f"{H}/snap_estimated_agonal.tsv", "SNAP against the estimated agonal score (slopes within controls and cases; schizophrenia effect on SNAP with and without the estimated agonal score and pH)", "analyses/snap_signal_analysis/05_sn_covariates.py, 06_assemble.py", fs[:1] + [f"{H}/inputs/sn_est_agonal.tsv"])
pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_rows", 400)
print(A[A.analysis == "sz_effect"][["dataset", "model", "effect", "lo95", "hi95", "change_from_base", "change_lo95", "change_hi95", "n"]].round(3).to_string())
