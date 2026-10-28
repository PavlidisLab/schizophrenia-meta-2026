#!/usr/bin/env python3
"""SNAP score, tests, coherence, variance share and PC alignment in one bulk dataset (22 bulk registry datasets, GTEx Cortex / BA9 / BA24).  Run on remote:
  python3 03_bulk_main.py DATASET   (env: OMP/OPENBLAS/MKL threads = 1 or 2, run under nice)   Inputs: out/bulk/<dataset>_{expr.tsv.gz,cv.tsv} (01_export_bulk.R), GTEx files of pc1_ph_gtex_2026-10-12, snap_weights.tsv, agonal_axis.tsv.
Weights: w_total (primary), w_neuron, w_astro of 02_weights.py applied to the genes measured in the dataset (0 where a gene has no loading); score = sum(w z)/sum(|w|), z = gene standardized across the
  dataset's samples; the score is divided by its SD across samples.  Null: the same weights permuted among all measured genes within their expression decile (mean expression), 10,000 draws.
Tests (per score type total / neuron / astrocyte): sz = schizophrenia effect (SD) given age and sex; age_control = slope per decade in controls given sex (GTEx: all samples); recorded covariates = partial Spearman given age, sex
  (and diagnosis where both groups exist); GTEx: Hardy 4 vs 2 and ventilator (Hardy 0) vs other contrasts given age, sex, ischemic time; est_agonal = slope of the score (SD) per SD of the estimated agonal score, given age and sex, per group.
Interval of the observed value: 500 bootstrap resamples of samples.  'rel' = observed minus the mean of the null draws (interval: bootstrap of observed, null mean fixed).  P (TSV only): two-sided, from the null draws."""
import os, sys, json, numpy as np, pandas as pd
from scipy import stats, sparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from snaplib import *
H = os.path.dirname(os.path.abspath(__file__)); HOME = os.path.expanduser("snap_signal_analysis"); OUT = f"{HOME}/out"; GT = os.path.expanduser("pc1_ph_gtex_2026-10-12")
BNULL, NBOOT, NBOOT_COH = 10000, 500, 300
GREG = {"GTEx_Cortex": ("expr_Brain__Cortex.tsv.gz", "Brain - Cortex", "Cortex"), "GTEx_Frontal_Cortex_BA9": ("expr_Brain__Frontal_Cortex_BA9.tsv.gz", "Brain - Frontal Cortex (BA9)", "BA9"),
        "GTEx_Anterior_cingulate_cortex_BA24": ("expr_Brain__Anterior_cingulate_cortex_BA24.tsv.gz", "Brain - Anterior cingulate cortex (BA24)", "BA24")}
W = read_commented(f"{H}/snap_weights.tsv", index_col=0); AX = read_commented(f"{H}/agonal_axis.tsv", index_col=0)
def load(name):
    if name in GREG:
        f, smtsd, reg = GREG[name]; E = pd.read_csv(f"{GT}/{f}", sep="\t", index_col=0); s = pd.read_csv(f"{GT}/samples.tsv", sep="\t").set_index("SAMPID")
        s = s[s.SMTSD == smtsd].loc[E.columns]; cv = pd.DataFrame({"age": s.age.values, "male": s.male.values, "RIN": s.SMRIN.values, "hardy": s.DTHHRDY.values, "isch_h": s.isch_h.values}, index=E.columns)
        return E, cv, "gtex", reg
    E = pd.read_csv(f"{OUT}/bulk/{name}_expr.tsv.gz", sep="\t", index_col=0); cv = pd.read_csv(f"{OUT}/bulk/{name}_cv.tsv", sep="\t").set_index("gsm").loc[E.columns]
    cv = cv.rename(columns={"pH": "pH", "PMI": "PMI"}); return E, cv, str(cv.cluster.iloc[0]), None
def axis_for(exclude):
    cols = [c for c in ["Cortex", "BA9", "BA24"] if c != exclude]; return AX[cols].mean(axis=1, skipna=True).dropna()
def est_agonal(Z, genes, axis, drop):
    """higher = more agonal (slow-death-like, the longer terminal phase): weighted mean of z with the Hardy 4 vs 2 t (positive = higher in slow death) as weights; genes in `drop` removed first"""
    ax = axis[~axis.index.isin(drop)]; g = [i for i, x in enumerate(genes) if x in ax.index]; t = ax.reindex(genes[g]).values
    return (Z[:, g] @ t) / np.abs(t).sum()
def ols_std(Y, X, j):
    Y = Y / Y.std(0, ddof=1); return np.linalg.lstsq(X, Y, rcond=None)[0][j]
def mk(rows, **kw): rows.append(kw)
def process(name):
    rng = np.random.default_rng(abs(hash(name)) % (2**31) if False else int(sha_name(name), 16))
    E, cv, cluster, reg = load(name); genes = E.index.values; n = E.shape[1]
    Xr = E.T.values.astype(float); Z = zcols(Xr); dec = deciles(E)
    w = {k: W[c].reindex(genes).fillna(0).values for k, c in [("total", "w_total"), ("neuron", "w_neuron"), ("astrocyte", "w_astro")]}
    sc = {k: (Z @ v) / np.abs(v).sum() for k, v in w.items()}; scs = {k: s / s.std(ddof=1) for k, s in sc.items()}
    # null scores for the three weight vectors with common permutations
    NS = {k: np.empty((n, BNULL), dtype=np.float32) for k in w}
    for a in range(0, BNULL, 500):
        b = min(BNULL, a + 500); P = perm_within_decile(dec, b - a, rng); Zf = Z.astype(np.float32)
        for k, v in w.items(): NS[k][:, a:b] = (Zf @ v[P].T.astype(np.float32)) / np.abs(v).sum()
    NSs = {k: (NS[k] / NS[k].std(0, ddof=1)).astype(np.float64) for k in NS}
    topw = W.w_total.abs().sort_values(ascending=False).index
    ag = {}
    for lab, drop in [("ex1000", topw[:1000]), ("ex200", topw[:200]), ("none", [])]:
        ag[lab] = est_agonal(Z, genes, axis_for(reg), set(drop))
    # ---- covariates
    age = cv.age.values.astype(float); male = cv.male.values.astype(float) if "male" in cv else None
    sz = cv.sz.values.astype(float) if "sz" in cv.columns else None
    both = sz is not None and min(sz.sum(), (1 - sz).sum()) >= 5
    rows = []
    def run(test, group, rowmask, stat, label_extra=None):
        """stat(Y) with Y = n_sub x m score columns (order: observed total, neuron, astrocyte); null via the same stat on null matrices"""
        idx = np.where(rowmask)[0]
        if len(idx) < 15: return
        for k in ["total", "neuron", "astrocyte"]:
            obs = float(stat(scs[k][idx, None], idx)[0]); nl = np.asarray(stat(NSs[k][idx], idx)); nm = nl.mean()
            bs = []
            for _ in range(NBOOT):
                j = rng.integers(0, len(idx), len(idx)); ii = idx[j]
                try: bs.append(float(stat(scs[k][ii, None], ii)[0]))
                except Exception: pass
            bs = np.array(bs); lo, hi = pct(bs) if len(bs) > 50 else (np.nan, np.nan)
            p = (1 + np.sum(np.abs(nl - nm) >= abs(obs - nm))) / (len(nl) + 1)
            mk(rows, effect=obs, lo95=lo, hi95=hi, dataset=name, cluster=cluster, test=test, group=group, score=k, n=len(idx), null_mean=nm, null_lo95=np.percentile(nl, 2.5), null_hi95=np.percentile(nl, 97.5),
               rel=obs - nm, rel_lo95=lo - nm, rel_hi95=hi - nm, p_value_tsv_only=p)
    def ols_stat(cols_fn, j):
        return lambda Y, ii: ols_std(Y, np.column_stack([np.ones(len(ii))] + cols_fn(ii)), j)
    allm = np.ones(n, bool)
    if both:
        run("sz_effect", "all", allm, ols_stat(lambda ii: [sz[ii], age[ii], male[ii]], 1))
    ctl = (sz == 0) if sz is not None else allm
    run("age_per_decade", "controls" if sz is not None else "all", ctl, ols_stat(lambda ii: [age[ii] / 10, male[ii]], 1))
    for v in ["pH", "PMI", "RIN", "agonal"]:
        if v in cv.columns and cv[v].notna().sum() >= 20:
            x = cv[v].values.astype(float); m = ~np.isnan(x) & ~np.isnan(age)
            extra = (lambda ii: [age[ii], male[ii]] + ([sz[ii]] if both else []))
            run("partial_spearman_" + v, "all", m, lambda Y, ii, x=x, ex=extra: partial_spearman(Y, x[ii], ex(ii)))
    if name in GREG:
        hd = cv.hardy.values.astype(float); isch = cv.isch_h.values.astype(float)
        run("partial_spearman_hardy_0to4", "all", ~np.isnan(hd), lambda Y, ii: partial_spearman(Y, hd[ii], [age[ii], male[ii], isch[ii]]))
        run("partial_spearman_ischemic_time", "all", ~np.isnan(isch), lambda Y, ii: partial_spearman(Y, isch[ii], [age[ii], male[ii]]))
        h42 = np.isin(hd, [2, 4]); ind4 = (hd == 4).astype(float)
        run("hardy4_vs_2_SD", "Hardy 2 or 4", h42, ols_stat(lambda ii: [ind4[ii], age[ii], male[ii], isch[ii]], 1))
        ind0 = (hd == 0).astype(float)
        run("ventilator_vs_other_SD", "all", ~np.isnan(hd), ols_stat(lambda ii: [ind0[ii], age[ii], male[ii], isch[ii]], 1))
        run("male_vs_female_SD", "all", allm, ols_stat(lambda ii: [male[ii], age[ii]], 1))
    # SNAP against estimated agonal score (SD per SD), by group
    for lab in ["ex1000", "ex200"]:
        a = ag[lab]; az = (a - a.mean()) / a.std(ddof=1)
        groups = [("controls", sz == 0), ("cases", sz == 1)] if both else []
        groups.append(("all", allm))
        for gname, gm in groups:
            if gname == "all" and both:
                run(f"slope_on_est_agonal_{lab}", "all (given diagnosis)", allm, ols_stat(lambda ii: [az[ii], age[ii], male[ii], sz[ii]], 1))
            else:
                run(f"slope_on_est_agonal_{lab}", gname, np.asarray(gm), ols_stat(lambda ii: [az[ii], age[ii], male[ii]], 1))
    # ---- PCs
    Xc = Xr - Xr.mean(0); U, S, _ = np.linalg.svd(Xc, full_matrices=False); pve = S**2 / (S**2).sum(); PC = U[:, :5] * S[:5]
    prow = []
    for scname, vec in [("snap_total", scs["total"]), ("snap_neuron", scs["neuron"]), ("snap_astrocyte", scs["astrocyte"]), ("est_agonal_ex1000", ag["ex1000"])]:
        for i in range(5):
            r = stats.spearmanr(vec, PC[:, i])[0]; se = np.sqrt((1 + r**2 / 2) / (n - 3)); ci = np.tanh(np.arctanh(r) + np.array([-1, 1]) * 1.96 * se)
            nullabs = np.abs(partial_spearman(NSs["total"][:, :2000], PC[:, i], [])) if scname == "snap_total" else np.nan
            mk(prow, effect=r, lo95=ci[0], hi95=ci[1], dataset=name, cluster=cluster, score=scname, pc=i + 1, pc_variance_share=pve[i], abs_rho=abs(r), n=n,
               null_mean_abs_rho=np.mean(nullabs), p_value_tsv_only=(1 + np.sum(nullabs >= abs(r))) / (len(np.atleast_1d(nullabs)) + 1) if scname == "snap_total" else np.nan)
    # ---- coherence and variance share
    crow = coherence(name, cluster, E, Z, Xr, cv, genes, dec, n, sz, age, male, both, rng)
    ps = pd.DataFrame({"sample": E.columns, "snap_total": scs["total"], "snap_neuron": scs["neuron"], "snap_astrocyte": scs["astrocyte"], "est_agonal_ex1000": ag["ex1000"], "est_agonal_ex200": ag["ex200"],
                       "est_agonal_none": ag["none"], **{f"PC{i+1}": PC[:, i] for i in range(5)}}); ps = pd.concat([ps.set_index("sample"), cv.drop(columns=[c for c in ["cluster", "dataset"] if c in cv.columns])], axis=1)
    ps.to_csv(f"{OUT}/persample_{name}.tsv", sep="\t"); pd.DataFrame(rows).to_csv(f"{OUT}/tests_{name}.tsv", sep="\t", index=False)
    pd.DataFrame(prow).to_csv(f"{OUT}/pcs_{name}.tsv", sep="\t", index=False); pd.DataFrame(crow).to_csv(f"{OUT}/coh_{name}.tsv", sep="\t", index=False)
    print(name, n, len(genes), len(rows), len(crow), flush=True)
def sha_name(s): return __import__("hashlib").sha256(s.encode()).hexdigest()[:8]
def stat_sets(S_all, S_n, S_a, kN, kA, nobs):
    """mean pairwise correlation for columns of S (n x B or n x 1) built from standardized columns; sets all / neuron / astrocyte / cross"""
    n1 = nobs - 1; k = kN + kA; out = {}
    out["all"] = ((S_all**2).sum(0) / n1 - k) / (k * (k - 1))
    out["neuron"] = ((S_n**2).sum(0) / n1 - kN) / (kN * (kN - 1)) if kN > 1 else np.full(S_all.shape[1], np.nan)
    out["astrocyte"] = ((S_a**2).sum(0) / n1 - kA) / (kA * (kA - 1)) if kA > 1 else np.full(S_all.shape[1], np.nan)
    out["cross"] = (S_n * S_a).sum(0) / (n1 * kN * kA) if kN > 0 and kA > 0 else np.full(S_all.shape[1], np.nan)
    return out
def matched_sets(dec, sel, B, rng):
    """B x k array of gene indices: for each selected gene a random gene of its expression decile, without replacement within the draw"""
    M = np.empty((B, len(sel)), dtype=np.int32)
    for d in range(10):
        pos = np.where(dec[sel] == d)[0]
        if len(pos) == 0: continue
        pool = np.where(dec == d)[0]; keys = rng.random((B, len(pool))); ch = np.argpartition(keys, len(pos) - 1, axis=1)[:, :len(pos)] if len(pos) < len(pool) else np.argsort(keys, axis=1)
        M[:, pos] = pool[ch]
    return M
def indicator(M, signs, G):
    B, k = M.shape
    if k == 0: return sparse.csr_matrix((B, G), dtype=np.float32)
    return sparse.csr_matrix((np.broadcast_to(signs.astype(np.float32), (B, k)).ravel(), M.ravel(), np.arange(0, B * k + 1, k)), shape=(B, G))
def coherence(name, cluster, E, Z, Xr, cv, genes, dec, n, sz, age, male, both, rng):
    G = len(genes); rows = []
    wt = W.w_total.reindex(genes).fillna(0).values; wn = W.w_neuron.reindex(genes).fillna(0).values; wa = W.w_astro.reindex(genes).fillna(0).values
    order = np.argsort(-np.abs(wt)); order = order[np.abs(wt[order]) > 0]
    variants = [("raw", Xr, None)]
    cols = [age, male] + ([sz] if both else [])
    if cols[0] is not None and not np.isnan(age).any(): variants.append(("resid_dx_age_sex" if both else "resid_age_sex", residualize(Xr, cols), None))
    if both and (sz == 0).sum() >= 15: variants.append(("controls_resid_age_sex", residualize(Xr[sz == 0], [age[sz == 0], male[sz == 0]]), "ctl"))
    for top in (200, 1000):
        sel = order[:top]; half = np.abs(wn[sel]) > np.abs(wa[sel]); sgn = np.where(half, np.sign(wn[sel]), np.sign(wa[sel])); sgn[sgn == 0] = 1
        iN, iA = np.where(half)[0], np.where(~half)[0]; kN, kA = len(iN), len(iA)
        M = matched_sets(dec, sel, BNULL, rng)
        for vname, X, flag in variants:
            Zv = zcols(X); nn = X.shape[0]
            for aligned in (True, False):
                s = sgn if aligned else np.ones(len(sel))
                def sums(Mm, ss):
                    out = []
                    for sub in (np.arange(len(sel)), iN, iA):
                        I = indicator(Mm[:, sub], ss[sub], G); out.append(np.asarray(I @ Zv.T).T)       # n x B
                    return out
                obsM = sel[None, :]; So = sums(obsM, s); obs = stat_sets(So[0], So[1], So[2], kN, kA, nn)
                Sn = sums(M, s); nul = stat_sets(Sn[0], Sn[1], Sn[2], kN, kA, nn)
                sub200 = M[:200]; bsd = {k: [] for k in obs}
                for _ in range(NBOOT_COH):
                    j = rng.integers(0, nn, nn); Zb = zcols(X[j]); Mb = sub200; ob = []
                    for sub in (np.arange(len(sel)), iN, iA): ob.append(np.asarray(indicator(obsM[:, sub], s[sub], G) @ Zb.T).T)
                    mb = [np.asarray(indicator(Mb[:, sub], s[sub], G) @ Zb.T).T for sub in (np.arange(len(sel)), iN, iA)]
                    o = stat_sets(*ob, kN, kA, nn); m = stat_sets(*mb, kN, kA, nn)
                    for k in o: bsd[k].append(o[k][0] - np.nanmean(m[k]))
                for k in obs:
                    nl = nul[k]; nm = np.nanmean(nl); ef = float(obs[k][0] - nm); b = np.array(bsd[k]); lo, hi = pct(b[~np.isnan(b)]) if np.isfinite(b).sum() > 50 else (np.nan, np.nan)
                    p = (1 + np.sum(np.abs(nl - nm) >= abs(ef))) / (len(nl) + 1)
                    rows.append(dict(effect=ef, lo95=lo, hi95=hi, dataset=name, cluster=cluster, set=k, top_n=top, statistic="mean_pairwise_r", variant=vname, aligned=aligned, observed=float(obs[k][0]),
                                     matched_mean=nm, n_samples=nn, n_genes_neuron=kN, n_genes_astro=kA, p_value_tsv_only=p))
    # variance share of the SNAP direction (unit-length w_total over measured genes); on mean-centred expression and on standardized expression (genes scaled to unit variance, as in analyses/axis_variance_share)
    v = wt / np.linalg.norm(wt)
    def prep(X, scale):
        Xc = X - X.mean(0)
        if scale == "standardized":
            sd_ = Xc.std(0, ddof=1); sd_[sd_ == 0] = np.inf; Xc = Xc / sd_
        return Xc
    for vname, X, flag in variants:
        if vname == "controls_resid_age_sex": continue
        for scale in ("centred", "standardized"):
            Xc = prep(X, scale); tot = (Xc**2).sum(); sh = lambda V: ((Xc @ V)**2).sum(0) / tot
            obs = float(sh(v[:, None])[0]); nl_rand = np.empty(BNULL); nl_dec = np.empty(BNULL)
            for a_ in range(0, BNULL, 500):
                b_ = min(BNULL, a_ + 500); R = rng.standard_normal((G, b_ - a_)); R /= np.linalg.norm(R, axis=0); nl_rand[a_:b_] = sh(R)
                P = perm_within_decile(dec, b_ - a_, rng); nl_dec[a_:b_] = sh(v[P].T)
            bsv = []
            for _ in range(NBOOT_COH):
                j = rng.integers(0, X.shape[0], X.shape[0]); Xb = prep(X[j], scale); bsv.append(((Xb @ v)**2).sum() / (Xb**2).sum())
            bsv = np.array(bsv)
            for kind, nl in (("random_direction", nl_rand), ("decile_matched_weights", nl_dec)):
                nm = nl.mean(); lo, hi = pct(bsv - nm); p = (1 + np.sum(np.abs(nl - nm) >= abs(obs - nm))) / (BNULL + 1)
                rows.append(dict(effect=obs - nm, lo95=lo, hi95=hi, dataset=name, cluster=cluster, set="variance_share:" + kind, top_n=int((wt != 0).sum()), statistic="variance_share (" + scale + " expression)", variant=vname, aligned=np.nan,
                                 observed=obs, matched_mean=nm, n_samples=X.shape[0], n_genes_neuron=np.nan, n_genes_astro=np.nan, p_value_tsv_only=p))
    return rows
if __name__ == "__main__":
    process(sys.argv[1])
