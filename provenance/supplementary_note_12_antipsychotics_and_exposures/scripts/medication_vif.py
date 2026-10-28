"""STAT-012 addendum 3, medication: variance-inflation correction of the gene-set tests in analyses/antipsychotic_medication_analyses.

Sample-level expression is not in the repo (data/resultsets = Gemma contrast tables; data/samples = sample annotations), and 3 of the 20 units
have 2 or 3 animals per group (20 or fewer distinct label permutations), so the VIF route is used.
Estimator. For a class of units, the matrix X (genes x contrasts) holds the per-contrast z (statistic 'dir': sign_g * z; 'resp': |z|).
rho_set = mean pairwise Pearson correlation among the set's genes across the class's contrasts (pairs with >= MINP shared contrasts);
rho_bg = the same, averaged over NBG random matched sets (the null's own gene correlation). The empirical null already carries rho_bg, so
the set's variance relative to the null is c^2 = VIF_set / VIF_bg, VIF = 1 + (m - 1) rho (common.vif), c clipped at 1 (no narrowing).
corrected SD = null SD x c; interval = (obs - null mean) +- 1.96 x corrected SD. Cluster bootstrap over units gives an interval on rho_set.
Run: opt/anaconda3/bin/python3 medication_vif.py
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
MED = os.path.join(REPO, "analyses", "medication")
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import common  # noqa: E402
sys.path.insert(0, os.path.join(MED, "scripts"))
import analysis as A  # noqa: E402  (load, bins_for; no parallel loader)

NPERM, NBG, NBOOT, MINP, SEED = 5000, 300, 300, 4, 20261012
rng = np.random.default_rng(SEED)
Z, unit_class, sets, prior, expr = A.load()
# gene sets exactly as analysed in medication/results/set_tests.tsv: recurrence_genes.tsv has changed since that run (126 genes at q_expr_bh < 0.1
# now, 145 then), so the sets and signs are read back from results/gene_effects.tsv
GE = pd.read_csv(os.path.join(MED, "results", "gene_effects.tsv"), sep="\t", comment="#")
sets = {k: g.set_index("gene")["sz_direction"] for k, g in GE.groupby("set")}
genes = Z.index.values; gidx = {g: i for i, g in enumerate(genes)}
Zv = Z.values; units = list(Z.columns)
b = A.bins_for(genes, prior, expr); by_bin = {k: np.where(b == k)[0] for k in np.unique(b)}
classes = {c: [i for i, u in enumerate(units) if unit_class[u] == c] for c in ("chronic", "acute", "cells")}
L = pd.read_csv(os.path.join(MED, "results", "contrasts_long.tsv.gz"), sep="\t", comment="#")
L["cid"] = L["unit"] + "|" + L["contrast_fv"].astype(str)
CM = L.pivot_table(index="gene", columns="cid", values="z", aggfunc="mean").reindex(genes)   # genes x contrasts
cinfo = L.drop_duplicates("cid").set_index("cid")[["unit", "class"]].reindex(CM.columns)

def draw(idx, n):
    out = np.empty((n, len(idx)), dtype=int)
    for j, i in enumerate(idx):
        pool = by_bin[b[i]]; out[:, j] = pool[rng.integers(0, len(pool), n)]
    return out

def stats(D, sg):
    z = Zv[D]; resp = np.nanmean(np.abs(z), axis=1); ok = sg != 0
    return resp, np.nanmean(z[:, ok, :] * sg[ok][None, :, None], axis=1)

def rho(gi, sg, cols, stat):
    """mean pairwise correlation among genes gi over contrast columns cols"""
    X = CM.values[np.ix_(gi, cols)].astype(float)
    if stat == "resp": X = np.abs(X)
    else:
        ok = sg != 0; X = X[ok] * sg[ok][:, None]
    m = X.shape[0]; M = (~np.isnan(X)).astype(float); X0 = np.nan_to_num(X)
    n = M @ M.T; sx = X0 @ M.T; sxx = (X0 ** 2) @ M.T
    sxy = X0 @ X0.T
    with np.errstate(all="ignore"):
        mx = sx / n; my = sx.T / n
        cov = sxy / n - mx * my
        vx = sxx / n - mx ** 2; vy = sxx.T / n - my ** 2
        r = cov / np.sqrt(vx * vy)
    r[n < MINP] = np.nan; np.fill_diagonal(r, np.nan)
    return common.mean_offdiag(r)[0]

cvec = lambda cls: [j for j, cid in enumerate(CM.columns) if cinfo.loc[cid, "class"] in cls]
def cfac(m, rs, rb): return float(np.sqrt(max(1.0, common.vif(m, rs) / common.vif(m, rb))))

rows, rr = [], []
meta = {}
for sname, sdir in sets.items():
    idx = np.array([gidx[g] for g in sdir.index if g in gidx])
    sg = np.array([{"up": 1, "down": -1}.get(sdir[genes[i]], 0) for i in idx])
    obs_r, obs_d = stats(idx[None, :], sg)
    Dn = draw(idx, NPERM); null_r, null_d = stats(Dn, sg)
    bgsets = [draw(idx, 1)[0] for _ in range(NBG)]
    for stat, obs, null in (("resp", obs_r[0], null_r), ("dir", obs_d[0], null_d)):
        for rclass, rcols_cls in (("own", None), ("all", ("chronic", "acute", "cells"))):
            pass
        # rho per class: own class, and all contrasts (sensitivity)
        rho_tab = {}
        for key, cls in (("chronic", ("chronic",)), ("acute", ("acute",)), ("cells", ("cells",)), ("all", ("chronic", "acute", "cells"))):
            cols = cvec(cls)
            rs = rho(idx, sg, cols, stat)
            rbs = np.array([rho(f, sg, cols, stat) for f in bgsets])
            # cluster bootstrap over units for rho_set
            us = sorted(set(cinfo.iloc[cols]["unit"])); bs = []
            if len(us) >= 3:
                for _ in range(NBOOT):
                    pick = rng.choice(us, len(us)); cc = [j for u in pick for j in cols if cinfo.iloc[j]["unit"] == u]
                    bs.append(rho(idx, sg, cc, stat))
                lo, hi = np.nanpercentile(bs, [2.5, 97.5])
            else: lo = hi = np.nan
            rho_tab[key] = (rs, np.nanmean(rbs), np.nanstd(rbs), lo, hi, len(cols))
            rr.append({"set": sname, "stat": stat, "contrast_class": key, "n_contrasts": len(cols), "rho_set": rs, "rho_bg_mean": np.nanmean(rbs),
                       "rho_bg_sd": np.nanstd(rbs), "rho_set_boot_lo": lo, "rho_set_boot_hi": hi, "excess": rs - np.nanmean(rbs)})
        def emit(label, cls, ulist, key, k_alt=None):
            rs, rb, _, lo, hi, _ = rho_tab[key]
            o = np.nanmean(obs[ulist]); nl = np.nanmean(null[:, ulist], axis=1); eff = o - nl.mean(); sd = nl.std(ddof=1)
            m = np.mean([np.sum(~np.isnan(Zv[idx, u])) for u in ulist]) if label != "unit" else np.sum(~np.isnan(Zv[idx, ulist[0]]))
            if stat == "dir": m = np.mean([np.sum(~np.isnan(Zv[idx[sg != 0], u])) for u in ulist])
            c = cfac(m, rs, rb); c_hi = cfac(m, hi, rb) if np.isfinite(hi) else np.nan
            return dict(rs=rs, o=o, nullmean=nl.mean(), eff=eff, sd=sd, m=m, c=c, c_hi=c_hi, p_old=None, nl=nl)
        for key, cls in (("chronic", "chronic"), ("acute", "acute"), ("cells", "cells")):
            ul = classes[cls]
            for alt in ("own", "all"):
                rk = key if alt == "own" else "all"
                for u in [None] + ([] if cls == "x" else ul):
                    if u is None: e = emit("class", cls, ul, rk); name = f"ALL {cls} ({len(ul)} units)"
                    else: e = emit("unit", cls, [u], rk); name = units[u]
                    half_old, half_new = 1.96 * e["sd"], 1.96 * e["sd"] * e["c"]
                    z_old = e["eff"] / e["sd"]
                    rows.append({"set": sname, "stat": stat, "rho_basis": alt, "class": cls, "unit": name, "m_genes": round(e["m"], 1),
                                 "obs": e["o"], "null_mean": e["nullmean"], "effect": e["eff"], "null_sd": e["sd"],
                                 "ci_old_lo": e["eff"] - half_old, "ci_old_hi": e["eff"] + half_old,
                                 "c": e["c"], "ci_corr_lo": e["eff"] - half_new, "ci_corr_hi": e["eff"] + half_new,
                                 "c_upper_rho": e["c_hi"], "c_naive": float(np.sqrt(common.vif(e["m"], e["rs"]))), "z_naive": z_old / float(np.sqrt(common.vif(e["m"], e["rs"]))),
                                 "ci_sens_lo": e["eff"] - 1.96 * e["sd"] * (e["c_hi"] if np.isfinite(e["c_hi"]) else e["c"]),
                                 "ci_sens_hi": e["eff"] + 1.96 * e["sd"] * (e["c_hi"] if np.isfinite(e["c_hi"]) else e["c"]),
                                 "z_old": z_old, "z_corr": z_old / e["c"],
                                 "p_old_norm": 2 * norm.sf(abs(z_old)), "p_corr_norm": 2 * norm.sf(abs(z_old / e["c"]))})
        # leave-one-out, chronic, rho recomputed without the unit (own-class basis)
        cu = classes["chronic"]
        for drop in cu:
            k = [u for u in cu if u != drop]
            cols = [j for j in cvec(("chronic",)) if cinfo.iloc[j]["unit"] != units[drop]]
            rs = rho(idx, sg, cols, stat); rb = np.nanmean([rho(f, sg, cols, stat) for f in bgsets[:100]])
            o = np.nanmean(obs[k]); nl = np.nanmean(null[:, k], axis=1); eff = o - nl.mean(); sd = nl.std(ddof=1)
            m = np.mean([np.sum(~np.isnan(Zv[idx[(sg != 0) if stat == "dir" else slice(None)], u])) for u in k])
            c = cfac(m, rs, rb)
            rows.append({"set": sname, "stat": stat, "rho_basis": "loo", "class": "chronic", "unit": f"leave out {units[drop]}", "m_genes": round(m, 1),
                         "obs": o, "null_mean": nl.mean(), "effect": eff, "null_sd": sd, "ci_old_lo": eff - 1.96 * sd, "ci_old_hi": eff + 1.96 * sd,
                         "c": c, "ci_corr_lo": eff - 1.96 * sd * c, "ci_corr_hi": eff + 1.96 * sd * c, "c_upper_rho": np.nan, "ci_sens_lo": np.nan, "ci_sens_hi": np.nan,
                         "z_old": eff / sd, "z_corr": eff / sd / c, "p_old_norm": 2 * norm.sf(abs(eff / sd)), "p_corr_norm": 2 * norm.sf(abs(eff / sd / c))})
R = pd.DataFrame(rows); RR = pd.DataFrame(rr)
old = pd.read_csv(os.path.join(MED, "results", "set_tests.tsv"), sep="\t", comment="#")
srcs = {"set_tests.tsv": common.sha(os.path.join(MED, "results", "set_tests.tsv")), "gene_effects.tsv": common.sha(os.path.join(MED, "results", "gene_effects.tsv")), "contrasts_long.tsv.gz": common.sha(os.path.join(MED, "results", "contrasts_long.tsv.gz"))}
hdr = ["claim: medication gene-set tests (responsiveness, direction) recomputed with a variance-inflation correction for gene correlation",
       "script: analyses/gene_correlation_adjusted_tests/addendum3_2026-10-12/medication/medication_vif.py", f"source sha256: {srcs}",
       f"date: {common.today()}; NPERM {NPERM}, NBG {NBG}, NBOOT {NBOOT}, MINP {MINP}, seed {SEED}",
       "effect = observed - null mean (resp: mean |z|; dir: mean sign x z); ci_old = effect +- 1.96 null SD; ci_corr = effect +- 1.96 null SD x c; c = sqrt(VIF_set/VIF_bg), clipped at 1",
       "c_naive = sqrt(VIF_set) alone (precedent form for an independent-gene null; double-counts the background correlation already in the empirical null): sensitivity only",
       "ci_sens uses the upper 95% bootstrap bound of rho_set (chronic only); P columns are normal-approximation on z"]
common.write_tsv(os.path.join(HERE, "medication_vif_results.tsv"), hdr, R.round(5))
common.write_tsv(os.path.join(HERE, "medication_rho_estimates.tsv"), hdr[:4], RR.round(5))

# summary: old P values from set_tests.tsv beside the corrected intervals (own-class rho basis)
def _key(s): return "study_brain" if "study_brain" in s else "expression-matched"
own = R[R.rho_basis == "own"].copy(); own["setk"] = own["set"].map(_key); old["setk"] = old["set"].map(_key)
rows2 = []
for _, r in own.iterrows():
    o = old[(old.setk == r.setk) & (old.unit == r.unit)]
    if len(o):
        rows2.append({**r.to_dict(), "p_old_perm": float(o["resp_p" if r["stat"] == "resp" else "dir_p"].iloc[0])})
S = pd.DataFrame(rows2).drop(columns=["setk"])
common.write_tsv(os.path.join(HERE, "medication_summary_old_vs_corrected.tsv"), hdr, S.round(5))
