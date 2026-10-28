#!/usr/bin/env python3
"""SNAP score of single-nucleus donors against recorded covariates, cell-type proportions, nuclei counts, the hypoxia score and the estimated agonal score; and the schizophrenia effect with/without the estimated agonal score and pH.
Run locally (one core, a few minutes):  python3 05_sn_covariates.py   (paths relative to the repository)
Scores: arm_b_scores.tsv (Ling's fixed LF4 loadings on donor pseudobulk; 'SNAP total' and the LF4 glutamatergic, gabaergic and astrocyte scores).  Covariates: arm_b_scores.tsv (age, sex, PMI, RIN, pH, log nuclei), donor nuclei per Ling class from
snap_sz/pb/<key>/cols_class.tsv (inputs/), hypoxia_donor_score.tsv (analyses/snap_ph_axes_2026-10-12; CMC and GSE254569 only; the hypoxia score and the SNAP score come from the same counts), sn_est_agonal.tsv (04_sn_agonal.py).
Statistic: partial Spearman correlation given age and sex (ranks; groups: all donors, controls only, and all donors given diagnosis too); 400 donor-bootstrap resamples (1,000 for the slope and schizophrenia-effect tables).  Pool = the repo's primary pool of independent banks
(HBCC_main, MSSM_main, SZBD_noLing, NSW; CMC is shown per cohort only because it shares Mount Sinai donors with MSSM_main): Fisher-z mean with weights n - 3, bootstrap stratified by bank.
Estimated-agonal slope: OLS slope of the score (SD) on the estimated agonal score (SD), given age and sex, in controls and cases separately.
Schizophrenia effect: OLS of score/SD on sz + age + sex + PMI + log nuclei (Arm B primary model; sz includes the 7 schizoaffective donors of GSE254569, as in arm_b_scores.tsv), then + estimated agonal, + estimated agonal + pH, + pH;
on the donors that have every variable of the richest model for that bank; donor bootstrap stratified by diagnosis (1,000); pooled over HBCC_main, MSSM_main and NSW (the banks with an estimated agonal score and no shared donors) with fixed weights (inverse of the base model's OLS variance)."""
import os, sys, numpy as np, pandas as pd
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(H, "../.."))
sys.path.insert(0, H); from snaplib import write_tsv, pct, partial_spearman
rng = np.random.default_rng(20261012); NB = 400
def rc(p, **k): return pd.read_csv(p, sep="\t", comment="#", **k)
S = rc(f"{REPO}/side-studies/coexpression_claims_survey/snap_sz/results/arm_b_scores.tsv", dtype={"donor_id": str})
KEY = {"CMC": "cmc", "NSW": "gse254569", "HBCC_main": "psychad_hbcc", "MSSM_main": "psychad_mssm", "Ling": "ling", "Batiuk": "batiuk", "SZBD_noLing": "szbd"}
CLS = ["glutamatergic", "gabaergic", "astrocyte", "oligodendrocyte", "polydendrocyte", "microglia", "endothelia"]
frames = []
for st, k in KEY.items():
    d = S[S.set == st].copy(); c = rc(f"{H}/inputs/cols_class_{k}.tsv", dtype={"donor_id": str})
    nn = c.pivot_table(index="donor_id", columns="cell_type", values="n_nuclei", aggfunc="sum").reindex(columns=CLS).fillna(0); tot = nn.sum(1)
    d = d.merge(pd.DataFrame({"donor_id": nn.index, "nuclei_total": tot.values, **{f"prop_{c_}": (nn[c_] / tot).values for c_ in CLS}}), on="donor_id", how="left"); d["log_nuclei_total"] = np.log(d.nuclei_total); frames.append(d)
D = pd.concat(frames, ignore_index=True)
hy = rc(f"{REPO}/analyses/snap_ph_axes_2026-10-12/hypoxia_donor_score.tsv", dtype={"donor_id": str}).replace({"set": {"cmc": "CMC", "gse254569": "NSW"}})[["set", "donor_id", "S_hypoxia"]]
ea = rc(f"{H}/inputs/sn_est_agonal.tsv", dtype={"donor_id": str}); D = D.merge(hy, on=["set", "donor_id"], how="left").merge(ea, on=["set", "donor_id"], how="left")
SC = {"SNAP total": "total", "LF4 glutamatergic": "glutamatergic neuron", "LF4 astrocyte": "astrocyte"}
COV = ["pmi", "rin", "ph", "log_nuclei_total", "S_hypoxia", "est_agonal_ex1000", "est_agonal_ex200"] + [f"prop_{c}" for c in CLS]
POOL = ["HBCC_main", "MSSM_main", "SZBD_noLing", "NSW"]
def groups(d):
    return [("all donors", d, False), ("controls", d[d.sz == 0], False)]
def pstat(d, y, v, dx):
    ok = d[[y, v, "age", "male"]].notna().all(1).values; dd = d[ok]
    if len(dd) < 20 or dd[v].nunique() < 3: return None
    return partial_spearman(dd[y].values, dd[v].values, [dd.age.values, dd.male.values.astype(float)] + ([dd.sz.values.astype(float)] if dx else [])), len(dd)
def boot_rows(sets_data, y, v, gname, dx):
    """sets_data: dict set -> df.  returns per-set rows and a pooled row"""
    out = []; est = {}
    for st, d in sets_data.items():
        r = pstat(d, y, v, dx)
        if r is None: continue
        bs = []
        for _ in range(NB):
            dd = d.iloc[rng.integers(0, len(d), len(d))]; q = pstat(dd, y, v, dx); bs.append(q[0] if q else np.nan)
        lo, hi = pct(np.array(bs)[~np.isnan(bs)]); est[st] = r
        out.append(dict(effect=r[0], lo95=lo, hi95=hi, set=st, group=gname, score=SC[y], covariate=v, statistic="partial Spearman (age, sex" + (", diagnosis" if dx else "") + ")", n=r[1], p_value_tsv_only=2 * stats.t.sf(abs(r[0]) * np.sqrt((r[1] - 3) / (1 - r[0] ** 2)), r[1] - 3)))
    ps = [s for s in POOL if s in est]
    if len(ps) >= 2:
        def pool(dd_by):
            zs, ws = [], []
            for s in ps:
                q = pstat(dd_by[s], y, v, dx)
                if q: zs.append(np.arctanh(np.clip(q[0], -.999, .999))); ws.append(q[1] - 3)
            return np.tanh(np.average(zs, weights=ws))
        pt = pool({s: sets_data[s] for s in ps}); bs = [pool({s: sets_data[s].iloc[rng.integers(0, len(sets_data[s]), len(sets_data[s]))] for s in ps}) for _ in range(NB)]
        lo, hi = pct(bs); n = sum(est[s][1] for s in ps)
        out.append(dict(effect=pt, lo95=lo, hi95=hi, set="pooled: " + ", ".join(ps), group=gname, score=SC[y], covariate=v, statistic="partial Spearman (age, sex" + (", diagnosis" if dx else "") + "), Fisher-z pool",
                        n=n, p_value_tsv_only=np.nan))
    return out
rows = []
for y in SC:
    for v in COV:
        for gname, _, dx in groups(D):
            sd = {}
            for st in KEY:
                d = D[D.set == st]; d = d[d.sz == 0] if gname == "controls" else d
                if len(d) >= 20 and d[v].notna().sum() >= 20: sd[st] = d
            if sd: rows += boot_rows(sd, y, v, gname, dx)
cov = pd.DataFrame(rows)
# slopes on the estimated agonal score, controls and cases
def slope(d, y, x="est_agonal_ex1000"):
    dd = d[[y, x, "age", "male"]].dropna()
    if len(dd) < 15: return None
    yy = dd[y].values / dd[y].std(ddof=1); xx = (dd[x].values - dd[x].mean()) / dd[x].std(ddof=1); X = np.column_stack([np.ones(len(dd)), xx, dd.age.values, dd.male.values.astype(float)])
    b, *_ = np.linalg.lstsq(X, yy, rcond=None); r = yy - X @ b; s2 = r @ r / (len(dd) - 4); se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]); return b[1], se, len(dd)
NB = 1000
sl = []
for est in ("est_agonal_ex1000", "est_agonal_ex200"):
    for gname, gf in (("controls", lambda d: d[d.sz == 0]), ("cases", lambda d: d[d.sz == 1])):
        per = {}
        for st in KEY:
            d = gf(D[D.set == st])
            if d[est].notna().sum() < 15: continue
            r = slope(d, "SNAP total", est); bs = [slope(d.iloc[rng.integers(0, len(d), len(d))], "SNAP total", est) for _ in range(NB)]; bs = [b[0] for b in bs if b]; lo, hi = pct(bs); per[st] = (r, d)
            sl.append(dict(effect=r[0], lo95=lo, hi95=hi, dataset=st, cohort_group="single-nucleus donors", model=f"SNAP total (SD) on {est} (SD), given age and sex", group=gname, n=r[2], p_value_tsv_only=2 * stats.t.sf(abs(r[0] / r[1]), r[2] - 4)))
        ps = [s for s in ["HBCC_main", "MSSM_main", "NSW"] if s in per]
        if len(ps) >= 2:
            def pool(dd):
                es = [(slope(dd[s], "SNAP total", est)) for s in ps]; w = np.array([1 / e[1] ** 2 for e in es]); return float(np.sum(w * np.array([e[0] for e in es])) / w.sum())
            pt = pool({s: per[s][1] for s in ps}); bs = [pool({s: per[s][1].iloc[rng.integers(0, len(per[s][1]), len(per[s][1]))] for s in ps}) for _ in range(NB)]; lo, hi = pct(bs)
            sl.append(dict(effect=pt, lo95=lo, hi95=hi, dataset="pooled: " + ", ".join(ps), cohort_group="single-nucleus donors", model=f"SNAP total (SD) on {est} (SD), given age and sex; inverse-variance pool, bootstrap stratified by bank", group=gname, n=sum(per[s][0][2] for s in ps), p_value_tsv_only=np.nan))
# schizophrenia effect with and without the estimated agonal score
def fit(d, extras):
    y = d["SNAP total"].values / d["SNAP total"].std(ddof=1); X = np.column_stack([np.ones(len(d)), d.sz.values.astype(float)] + [d[c].values.astype(float) for c in ["age", "male", "pmi", "lognuc"] + extras])
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; s2 = r @ r / (len(d) - X.shape[1]); return b[1], np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]), len(d) - X.shape[1]
MODELS = [("base: age, sex, PMI, nuclei", []), ("+ estimated agonal (top 1,000 removed)", ["est_agonal_ex1000"]), ("+ estimated agonal (top 200 removed)", ["est_agonal_ex200"])]
rowsA = []; est_by = {}
for st in ["CMC", "NSW", "HBCC_main", "MSSM_main"]:
    d = D[D.set == st].copy(); need = ["age", "male", "pmi", "lognuc", "est_agonal_ex1000", "est_agonal_ex200"] + (["ph"] if d.ph.notna().any() else [])
    d = d.dropna(subset=need).reset_index(drop=True); models = MODELS + ([("+ pH", ["ph"]), ("+ estimated agonal (top 1,000 removed) + pH", ["est_agonal_ex1000", "ph"])] if "ph" in need else [])
    def allfits(dd): return [fit(dd, ex)[0] for _, ex in models]
    est = [fit(d, ex) for _, ex in models]; est_by[st] = (d, models, est)
    bs = []
    for _ in range(NB):
        ii = np.concatenate([np.where(d.sz.values == g)[0][rng.integers(0, (d.sz.values == g).sum(), (d.sz.values == g).sum())] for g in (0, 1)]); bs.append(allfits(d.iloc[ii]))
    bs = np.array(bs)
    for j, (mn, _) in enumerate(models):
        lo, hi = pct(bs[:, j]); dl = bs[:, j] - bs[:, 0]; dlo, dhi = pct(dl)
        rowsA.append(dict(effect=est[j][0], lo95=lo, hi95=hi, dataset=st, cohort_group="single-nucleus donors", model=mn, change_from_base=est[j][0] - est[0][0], change_lo95=dlo, change_hi95=dhi, n=len(d), n_sz=int(d.sz.sum()),
                          p_value_tsv_only=2 * stats.t.sf(abs(est[j][0] / est[j][1]), est[j][2])))
# pooled over HBCC_main, MSSM_main, NSW for the models they share (first three)
pl = ["HBCC_main", "MSSM_main", "NSW"]
def pooled(dd_by, j):
    es = [fit(dd_by[s], MODELS[j][1]) for s in pl]; w0 = np.array([1 / fit(dd_by[s], [])[1] ** 2 for s in pl]); return float(np.sum(w0 * np.array([e[0] for e in es])) / w0.sum())
for j, (mn, _) in enumerate(MODELS):
    pt = pooled({s: est_by[s][0] for s in pl}, j); pt0 = pooled({s: est_by[s][0] for s in pl}, 0); bs = []
    for _ in range(NB):
        dd = {}
        for s in pl:
            d = est_by[s][0]; ii = np.concatenate([np.where(d.sz.values == g)[0][rng.integers(0, (d.sz.values == g).sum(), (d.sz.values == g).sum())] for g in (0, 1)]); dd[s] = d.iloc[ii]
        bs.append((pooled(dd, j), pooled(dd, 0)))
    bs = np.array(bs); lo, hi = pct(bs[:, 0]); dlo, dhi = pct(bs[:, 0] - bs[:, 1])
    rowsA.append(dict(effect=pt, lo95=lo, hi95=hi, dataset="pooled: HBCC_main, MSSM_main, NSW", cohort_group="single-nucleus donors", model=mn, change_from_base=pt - pt0, change_lo95=dlo, change_hi95=dhi, n=sum(len(est_by[s][0]) for s in pl),
                      n_sz=sum(int(est_by[s][0].sz.sum()) for s in pl), p_value_tsv_only=np.nan))
src = [f"{REPO}/side-studies/coexpression_claims_survey/snap_sz/results/arm_b_scores.tsv", f"{H}/inputs/sn_est_agonal.tsv", f"{REPO}/analyses/snap_ph_axes_2026-10-12/hypoxia_donor_score.tsv"] + [f"{H}/inputs/cols_class_{k}.tsv" for k in KEY.values()]
cov.to_pickle(f"{H}/_cov.pkl"); pd.DataFrame(sl).to_pickle(f"{H}/_sl.pkl"); pd.DataFrame(rowsA).to_pickle(f"{H}/_adj.pkl")
write_tsv(cov, f"{H}/snap_covariates_snRNA.tsv", "SNAP score of single-nucleus donors against PMI, RIN, pH, nuclei counts, cell-type proportions, hypoxia score and estimated agonal score (partial Spearman, bootstrap intervals)", "analyses/snap_signal_analysis/05_sn_covariates.py", src)
print(cov[(cov.score == "total") & (cov.group == "all donors")][["set", "covariate", "effect", "lo95", "hi95", "n"]].round(3).to_string())
print(pd.DataFrame(sl).round(3).to_string()); print(pd.DataFrame(rowsA).round(3).to_string())
