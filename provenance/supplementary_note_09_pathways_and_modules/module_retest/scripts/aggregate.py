#!/usr/bin/env python3
"""Pool the per-dataset module effects (per_dataset_module.tsv) over cohort clusters, per module and over modules/papers.

Independence: a dataset counts for a module only when none of its lineage labels (dataset_lineage.tsv) is among the module's discovery lineage labels (modules_index.tsv, bank_overlap.py rules). Modules whose lineage has an 'unnamed'
part are flagged (lineage_named = no): the unnamed part cannot be excluded.
Per module pool (primary): bulk datasets only, one dataset per cohort cluster (the one with the largest n_sz + n_ctl), DerSimonian-Laird random effects with the Hartung-Knapp-Sidik-Jonkman t interval on k - 1 df (as in snap_sz/combine.py).
Sensitivity pool: all bulk datasets, fixed-effect inside a cluster, then random effects across clusters.
Across modules: modules of one paper share genes and donors, so the unit is the paper (Gandal 2018 and Kang and Oldham 2026 one each): mean of its modules' pooled effects; random-effects mean over papers (HKSJ); sign share counted per module and per paper;
the 95% interval of the module-level share is a paper-cluster bootstrap (resample papers, 10000).
    python3 aggregate.py DIR
"""
import sys, hashlib
import numpy as np, pandas as pd
from scipy import stats
D = sys.argv[1]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
R = pd.read_csv(f"{D}/per_dataset_module.tsv", sep="\t", comment="#")
idx = pd.read_csv(f"{D}/modules_index.tsv", sep="\t", comment="#", dtype=str); idx["n_genes"] = idx.n_genes.astype(int)
lin = pd.read_csv(f"{D}/dataset_lineage.tsv", sep="\t", comment="#", dtype=str).set_index("key")
rng = np.random.default_rng(20261023)
def meta(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, se = y[ok], se[ok]; k = len(y)
    if k == 0: return dict(k=0)
    if k == 1: return dict(k=1, est=y[0], lo=np.nan, hi=np.nan, p=np.nan, I2=np.nan)
    w = 1 / se ** 2; fe = (w * y).sum() / w.sum(); Q = (w * (y - fe) ** 2).sum(); C = w.sum() - (w ** 2).sum() / w.sum(); tau2 = max(0.0, (Q - (k - 1)) / C)
    ws = 1 / (se ** 2 + tau2); mu = (ws * y).sum() / ws.sum(); q = (ws * (y - mu) ** 2).sum() / (k - 1); se_hk = np.sqrt(max(q, 1e-12) / ws.sum()); tc = stats.t.ppf(0.975, k - 1)
    return dict(k=k, est=mu, lo=mu - tc * se_hk, hi=mu + tc * se_hk, p=2 * stats.t.sf(abs(mu / se_hk), k - 1), I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)
R = R.merge(lin[["compartment", "cluster", "lineage"]], left_on="key", right_index=True).rename(columns={"lineage": "ds_lineage"})
R = R.merge(idx[["module_uid", "lineage", "source", "paper", "module", "direction"]], on="module_uid")
def indep(r):
    m = set(r.lineage.split(";")) - {"unnamed"}; d = set(r.ds_lineage.split(";")) - {"own"}
    return not (m & d)
R["independent"] = R.apply(indep, axis=1)
R["n"] = R.n_sz + R.n_ctl
B = R[R.compartment == "bulk"]
out = []
for (uid, model), g in B.groupby(["module_uid", "model"]):
    gi = g[g.independent]
    if len(gi) == 0: out.append(dict(module_uid=uid, model=model, n_datasets_independent=0, n_clusters=0, n_datasets_all=len(g))); continue
    one = gi.sort_values("n", ascending=False).drop_duplicates("cluster")
    m1 = meta(one.beta, one.se)
    cl = []
    for c, h in gi.groupby("cluster"):
        w = 1 / h.se ** 2; cl.append((float((w * h.beta).sum() / w.sum()), float(np.sqrt(1 / w.sum()))))
    m2 = meta([a for a, b in cl], [b for a, b in cl])
    out.append(dict(module_uid=uid, model=model, n_datasets_all=len(g), n_datasets_independent=len(gi), n_clusters=len(one), n_sz=int(one.n_sz.sum()), n_ctl=int(one.n_ctl.sum()),
                    est=m1.get("est"), lo=m1.get("lo"), hi=m1.get("hi"), p=m1.get("p"), I2=m1.get("I2"), est_allbulk=m2.get("est"), lo_allbulk=m2.get("lo"), hi_allbulk=m2.get("hi"), p_allbulk=m2.get("p")))
P = pd.DataFrame(out).merge(idx, on="module_uid")
P["lineage_named"] = np.where(P.lineage.str.contains("unnamed"), "no", "yes")
P["matches_paper"] = np.where(P.est > 0, 1, np.where(P.est < 0, 0, np.nan)); P.loc[P.est.isna(), "matches_paper"] = np.nan
P["interval_excludes_0"] = np.where(P.lo.notna(), ((P.lo > 0) | (P.hi < 0)).astype(float), np.nan)
cols = ["module_uid", "source", "paper", "module", "direction", "lineage", "lineage_named", "n_genes", "model", "n_datasets_all", "n_datasets_independent", "n_clusters", "n_sz", "n_ctl", "est", "lo", "hi", "p", "I2", "est_allbulk", "lo_allbulk", "hi_allbulk", "p_allbulk", "matches_paper", "interval_excludes_0"]
P = P[cols].round(4)
hdr = f"# aggregate.py (module_retest); per_dataset_module.tsv sha256 {sha(f'{D}/per_dataset_module.tsv')}; modules_index.tsv {sha(f'{D}/modules_index.tsv')}; dataset_lineage.tsv {sha(f'{D}/dataset_lineage.tsv')}\n# est/lo/hi/p: one bulk dataset per cohort cluster independent of the module's discovery lineage, DL random effects with HKSJ interval; *_allbulk: all independent bulk datasets, fixed effect inside cluster; effect = SZ minus control in SD of the direction-signed module score\n"
with open(f"{D}/module_pooled.tsv", "w") as fh: fh.write(hdr); P.to_csv(fh, sep="\t", index=False)
# ---- across modules and papers
def unit(r): return r.source if r.source != "paper_supplement" else r.paper
S = []
for model in ["m0", "m1", "m2", "m3"]:
    for named in ["yes", "all"]:
        for minc in [1, 2]:
            q = P[(P.model == model) & P.est.notna() & (P.n_clusters >= minc)]
            if named == "yes": q = q[q.lineage_named == "yes"]
            if len(q) == 0: continue
            q = q.assign(unit=q.apply(unit, axis=1))
            ub = q.groupby("unit").agg(mean_est=("est", "mean"), k_mod=("est", "size"), share=("matches_paper", "mean")).reset_index()
            # unit-level random effects over the spread of module means (between-unit variance only): HKSJ on unit means with unit-level SE = SD of the means / sqrt(k)
            k = len(ub); mu = ub.mean_est.mean(); sd = ub.mean_est.std(ddof=1) if k > 1 else np.nan; tc = stats.t.ppf(0.975, k - 1) if k > 1 else np.nan
            lo, hi = (mu - tc * sd / np.sqrt(k), mu + tc * sd / np.sqrt(k)) if k > 1 else (np.nan, np.nan)
            pv = 2 * stats.t.sf(abs(mu / (sd / np.sqrt(k))), k - 1) if k > 1 and sd > 0 else np.nan
            # bootstrap of the module-level sign share over units
            groups = [h.matches_paper.values for _, h in q.groupby("unit")]
            bs = []
            for _ in range(10000):
                pick = rng.integers(0, len(groups), len(groups)); v = np.concatenate([groups[i] for i in pick]); bs.append(np.nanmean(v))
            ns = int((ub.share > 0.5).sum()); nt = int((ub.share != 0.5).sum())
            S.append(dict(model=model, lineage_named_only=named, min_clusters=minc, n_modules=len(q), n_units=k, mean_effect_over_units=mu, lo=lo, hi=hi, p_units_t=pv,
                          module_share_matching=float(np.nanmean(q.matches_paper)), share_lo=float(np.percentile(bs, 2.5)), share_hi=float(np.percentile(bs, 97.5)),
                          units_majority_matching=ns, units_with_majority=nt, sign_test_p_units=float(stats.binomtest(ns, nt, 0.5).pvalue) if nt else np.nan,
                          modules_interval_excludes_0_matching=int(((q.lo > 0)).sum()), modules_interval_excludes_0_opposite=int((q.hi < 0).sum())))
S = pd.DataFrame(S).round(4)
hdr2 = hdr + "# across modules: unit = paper (Gandal 2018, Kang and Oldham 2026 one unit each); mean_effect_over_units = mean of the units' mean module effects, t interval on units - 1 df; share_* = fraction of modules with a pooled effect in the paper's direction, paper-cluster bootstrap interval\n"
with open(f"{D}/pooled_over_modules.tsv", "w") as fh: fh.write(hdr2); S.to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(S.to_string())
print(P.groupby(["model", "source"]).agg(n=("est", "size"), tested=("est", lambda s: s.notna().sum()), k2=("n_clusters", lambda s: (s >= 2).sum())).to_string())

# ---- per paper table (m1 and m3, lineage named or not, all modules)
Q = P[P.est.notna()].assign(unit=lambda d: d.apply(unit, axis=1))
rows = []
for u, h in Q.groupby("unit"):
    r = dict(unit=u, source=h.source.iloc[0], n_modules=h.module_uid.nunique(), lineage=";".join(sorted(set(";".join(h.lineage).split(";")))), lineage_named="yes" if (h.lineage_named == "yes").all() else "no")
    for m in ("m1", "m2", "m3"):
        x = h[h.model == m]
        r.update({f"{m}_modules": len(x), f"{m}_mean_est": x.est.mean() if len(x) else np.nan, f"{m}_share_matching": x.matches_paper.mean() if len(x) else np.nan, f"{m}_n_ci_matching": int((x.lo > 0).sum()), f"{m}_n_ci_opposite": int((x.hi < 0).sum())})
    rows.append(r)
U = pd.DataFrame(rows).round(4)
with open(f"{D}/by_paper.tsv", "w") as fh: fh.write(hdr + "# one row per paper (Gandal 2018 Table S5 and Kang and Oldham 2026 CoPA one row each); est = pooled effect of the module (see module_pooled.tsv), mean over the paper's modules\n"); U.to_csv(fh, sep="\t", index=False)
print(U[["unit", "n_modules", "lineage_named", "m1_mean_est", "m1_share_matching", "m3_mean_est", "m3_share_matching"]].to_string())
