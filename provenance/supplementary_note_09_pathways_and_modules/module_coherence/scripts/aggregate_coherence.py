#!/usr/bin/env python3
"""Pool per-dataset module coherence excess over cohort clusters, per module and over papers (rules of module_retest_2026-10-23/scripts/aggregate.py).
Per module and variant: one bulk dataset per cohort cluster (largest n), DerSimonian-Laird random effects, HKSJ t interval on k - 1 df; needs >= 2 clusters for an interval.
Across modules: unit = paper (Gandal 2018, CoPA one each); mean over a unit's modules, then mean over units with a t interval on units; share of modules with interval above 0 with a paper-cluster bootstrap (10000).
Usage: aggregate_coherence.py DIR   (DIR holds per_dataset_coherence.tsv, modules_index.tsv, dataset_lineage.tsv)
"""
import sys, hashlib
import numpy as np, pandas as pd
from scipy import stats
D = sys.argv[1]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
R = pd.read_csv(f"{D}/per_dataset_coherence.tsv", sep="\t", comment="#")
idx = pd.read_csv(f"{D}/modules_index.tsv", sep="\t", comment="#", dtype=str); idx["n_genes"] = idx.n_genes.astype(int)
rng = np.random.default_rng(20261023)
def meta(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); k = len(y)
    if k == 1: return dict(k=1, est=y[0], lo=np.nan, hi=np.nan, I2=np.nan)
    w = 1 / se ** 2; fe = (w * y).sum() / w.sum(); Q = (w * (y - fe) ** 2).sum(); C = w.sum() - (w ** 2).sum() / w.sum(); tau2 = max(0.0, (Q - (k - 1)) / C)
    ws = 1 / (se ** 2 + tau2); mu = (ws * y).sum() / ws.sum(); q = (ws * (y - mu) ** 2).sum() / (k - 1); se_hk = np.sqrt(max(q, 1e-12) / ws.sum()); tc = stats.t.ppf(0.975, k - 1)
    return dict(k=k, est=mu, lo=mu - tc * se_hk, hi=mu + tc * se_hk, I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)
out = []
for (uid, vn), g in R.groupby(["module_uid", "variant"]):
    g = g.assign(nn=g.n_sz + g.n_ctl); one = g.sort_values("nn", ascending=False).drop_duplicates("cluster")
    m = meta(one.excess, one.se)
    out.append(dict(module_uid=uid, variant=vn, n_datasets_independent=len(g), n_clusters=len(one), est=m["est"], lo=m["lo"], hi=m["hi"], I2=m["I2"], null_mean_avg=one.null_mean.mean(), observed_avg=one.observed.mean()))
P = pd.DataFrame(out).merge(idx[["module_uid", "source", "paper", "module", "direction", "lineage", "n_genes"]], on="module_uid")
P["lineage_named"] = np.where(P.lineage.str.contains("unnamed"), "no", "yes")
P["interval_above_0"] = np.where(P.lo.notna(), (P.lo > 0).astype(float), np.nan); P["interval_below_0"] = np.where(P.lo.notna(), (P.hi < 0).astype(float), np.nan)
P = P[["module_uid", "source", "paper", "module", "direction", "lineage", "lineage_named", "n_genes", "variant", "n_datasets_independent", "n_clusters", "observed_avg", "null_mean_avg", "est", "lo", "hi", "I2", "interval_above_0", "interval_below_0"]].round(4)
src = f"per_dataset_coherence.tsv sha256 {sha(f'{D}/per_dataset_coherence.tsv')}; modules_index.tsv {sha(f'{D}/modules_index.tsv')}"
def write(df, name, note):
    with open(f"{D}/{name}", "w") as fh: fh.write(f"# aggregate_coherence.py (module_coherence); {src}\n# {note}\n"); df.to_csv(fh, sep="\t", index=False)
write(P, "module_coherence_pooled.tsv", "est/lo/hi: coherence excess (mean pairwise r among the module's genes minus expression-matched random genes), one bulk dataset per cohort cluster independent of the module's lineage, DL random effects with HKSJ interval; interval_* blank when only one cluster; variant raw or resid (diagnosis, age, sex, region removed)")
def unit(r): return r.source if r.source != "paper_supplement" else r.paper
S = []; U = []
for vn in ["raw", "resid"]:
    for named in ["yes", "all"]:
        for minc in [1, 2]:
            q = P[(P.variant == vn) & P.est.notna() & (P.n_clusters >= minc)]
            if named == "yes": q = q[q.lineage_named == "yes"]
            if len(q) == 0: continue
            q = q.assign(unit=q.apply(unit, axis=1))
            ub = q.groupby("unit").agg(mean_excess=("est", "mean"), k_mod=("est", "size"), share_above_0=("interval_above_0", "mean")).reset_index()
            k = len(ub); mu = ub.mean_excess.mean(); sd = ub.mean_excess.std(ddof=1); tc = stats.t.ppf(0.975, k - 1)
            gro = [h.interval_above_0.values for _, h in q.groupby("unit")]
            bs = [np.nanmean(np.concatenate([gro[i] for i in rng.integers(0, len(gro), len(gro))])) for _ in range(10000)] if minc >= 2 else [np.nan]
            S.append(dict(variant=vn, lineage_named_only=named, min_clusters=minc, n_modules=len(q), n_units=k, mean_excess_over_units=mu, lo=mu - tc * sd / np.sqrt(k), hi=mu + tc * sd / np.sqrt(k),
                          median_module_excess=q.est.median(), min_module_excess=q.est.min(), modules_est_above_0=int((q.est > 0).sum()),
                          modules_interval_above_0=int((q.lo > 0).sum()) if minc >= 2 else np.nan, modules_interval_below_0=int((q.hi < 0).sum()) if minc >= 2 else np.nan,
                          share_interval_above_0=float(np.nanmean(q.interval_above_0)) if minc >= 2 else np.nan, share_lo=float(np.percentile(bs, 2.5)), share_hi=float(np.percentile(bs, 97.5))))
            if vn == "resid" and named == "all" and minc == 1: U = ub
write(pd.DataFrame(S).round(4), "pooled_over_modules.tsv", "unit = paper; mean_excess_over_units with t interval over units; share_interval_above_0 = share of modules (>= 2 clusters) whose pooled interval is above 0, paper-cluster bootstrap interval")
write(U.round(4), "by_paper.tsv", "variant resid, all modules; mean pooled coherence excess per unit")
print(pd.DataFrame(S).round(3).to_string())
# size and source breakdown, resid
q = P[(P.variant == "resid") & P.est.notna()]
print(q.groupby("source").agg(n=("est", "size"), med=("est", "median"), above=("interval_above_0", "sum")).round(3))
print(np.corrcoef(np.log(q.n_genes), q.est)[0, 1])
