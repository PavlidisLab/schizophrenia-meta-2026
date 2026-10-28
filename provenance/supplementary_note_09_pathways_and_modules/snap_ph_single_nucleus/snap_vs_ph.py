#!/usr/bin/env python3
"""SNAP score of single-nucleus donors against measured pH (GSE254569 NSW and brainSCOPE CMC), and against the estimated agonal score and a pH-axis score of the same donors. Run on frink (snap_sz/envs/py).
SNAP scores: arm_b_scores.tsv (Ling's fixed LF4 loadings on donor pseudobulk; columns 'SNAP total', 'SNAP-a' astrocyte, 'SNAP-n' neuron); covariates age, sex (male), PMI, log nuclei, pH from the same file.
Estimated agonal score: sn_est_agonal.tsv (04_sn_agonal.py; primary variant est_agonal_ex1000 = SNAP's top 1,000 genes removed from the axis; higher = more agonal).
pH-axis score (new here): the same construction as 04_sn_agonal.py with weights = control pH axis v2 Stouffer z (axis_v2_controls.tsv.gz; z > 0 = higher at higher pH), so a higher score = more like higher-pH tissue; SNAP's top 1,000 genes removed (variant ex1000).
Statistics per score and donor group (controls; all donors): Spearman rho with pH; partial Spearman given age and sex (and diagnosis in 'all'); OLS slope of the SNAP score (in SD) per SD of the covariate given age and sex (+ diagnosis in 'all'); 95% interval = percentile interval of 2,000 donor
bootstrap resamples (stratified by diagnosis). Usage: snap_vs_ph.py WORKDIR   (WORKDIR holds arm_b_scores.tsv, sn_est_agonal.tsv, axis_v2_controls.tsv.gz, arm_b_donors.tsv and the snap_what_is_it scripts/ folder; snap_sz/pb is read from ~/snap_sz/pb)
"""
import os, sys, hashlib
import numpy as np, pandas as pd
from scipy import stats, sparse as sp
W_ = sys.argv[1]; sys.path.insert(0, f"{W_}/scripts"); sys.dont_write_bytecode = True
from snap_psychad import z_matrix
PB = os.path.expanduser("~/snap_sz/pb"); KEY = {"NSW": "gse254569", "CMC": "cmc"}
def rc(p, **k): return pd.read_csv(p, sep="\t", comment="#", **k)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
S = rc(f"{W_}/arm_b_scores.tsv", dtype={"donor_id": str}); EA = rc(f"{W_}/sn_est_agonal.tsv", dtype={"donor_id": str})
AX = rc(f"{W_}/axis_v2_controls.tsv.gz", index_col=0)["z_stouffer"].dropna(); Wt = rc(f"{W_}/scripts/snap_weights.tsv", index_col=0); top = Wt.w_total.abs().sort_values(ascending=False).index
ax = AX[~AX.index.isin(set(top[:1000]))]
rows = []
for st, key in KEY.items():
    m = sp.load_npz(f"{PB}/{key}/counts_class.npz").tocsc(); genes = pd.read_csv(f"{PB}/{key}/genes.tsv", sep="\t", dtype=str); cols = pd.read_csv(f"{PB}/{key}/cols_class.tsv", sep="\t", dtype={"donor_id": str, "bioassay_id": str})
    donors = [d for d in S[S.set == st].donor_id if d in set(cols.donor_id)]
    num = np.zeros(len(donors)); den = np.zeros(len(donors))
    for c in ["glutamatergic", "gabaergic", "astrocyte"]:
        Z, _ = z_matrix(m, genes, cols, donors, c, ax.index.values); tt = ax.reindex(Z.columns).values
        ok = ~np.isnan(Z.values).all(1); num += np.where(ok, np.nan_to_num(Z.values) @ tt, 0); den += np.where(ok, np.abs(tt).sum(), 0)
    rows.append(pd.DataFrame({"set": st, "donor_id": donors, "ph_axis_score_ex1000": np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)}))
PA = pd.concat(rows)
D = S[S.set.isin(KEY)].merge(EA[["set", "donor_id", "est_agonal_ex1000"]], on=["set", "donor_id"], how="left").merge(PA, on=["set", "donor_id"], how="left")
D.to_csv("donor_table.tsv", sep="\t", index=False)
rng = np.random.default_rng(20261023); NB = 2000
def ranks(x): return stats.rankdata(x)
def resid(y, X): return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
def stat_set(y, x, cov):
    """returns spearman rho, partial spearman, OLS slope (SD per SD)"""
    rho = stats.spearmanr(y, x)[0]
    X = np.column_stack([np.ones(len(y)), cov]) if cov.shape[1] else np.ones((len(y), 1))
    pr = np.corrcoef(resid(ranks(y), np.column_stack([np.ones(len(y)), np.apply_along_axis(ranks, 0, cov)]) if cov.shape[1] else np.ones((len(y), 1))),
                     resid(ranks(x), np.column_stack([np.ones(len(y)), np.apply_along_axis(ranks, 0, cov)]) if cov.shape[1] else np.ones((len(y), 1))))[0, 1]
    xs = (x - x.mean()) / x.std(ddof=1); ys = (y - y.mean()) / y.std(ddof=1)
    Xd = np.column_stack([np.ones(len(y)), xs, cov]) if cov.shape[1] else np.column_stack([np.ones(len(y)), xs])
    b = np.linalg.lstsq(Xd, ys, rcond=None)[0][1]
    return rho, pr, b
out = []
for st in KEY:
    d = D[D.set == st]
    for grp in ("controls", "all"):
        g = d[d.sz == 0] if grp == "controls" else d
        for sc in ["SNAP total", "SNAP-a", "SNAP-n"]:
            for cv, lab in [("ph", "measured pH"), ("est_agonal_ex1000", "estimated agonal score"), ("ph_axis_score_ex1000", "pH-axis score")]:
                need = [sc, cv, "age", "male"] + (["sz"] if grp == "all" else [])
                h = g.dropna(subset=need).reset_index(drop=True)
                if len(h) < 15: continue
                covm = h[["age", "male"] + (["sz"] if grp == "all" else [])].values.astype(float)
                y = h[sc].values.astype(float); x = h[cv].values.astype(float)
                obs = stat_set(y, x, covm); bs = []
                strat = [np.where(h.sz == v)[0] for v in sorted(h.sz.unique())]
                for _ in range(NB):
                    ii = np.concatenate([rng.choice(s, len(s)) for s in strat])
                    try: bs.append(stat_set(y[ii], x[ii], covm[ii]))
                    except Exception: pass
                bs = np.array(bs)
                for k, nm in enumerate(["spearman_rho", "partial_spearman", "slope_SD_per_SD"]):
                    out.append(dict(set=st, group=grp, snap=sc, against=lab, statistic=nm, estimate=obs[k], lo95=np.percentile(bs[:, k], 2.5), hi95=np.percentile(bs[:, k], 97.5), n=len(h), n_sz=int((h.sz == 1).sum())))
R = pd.DataFrame(out).round(4)
hdr = (f"# snap_vs_ph.py (snap_ph_single_nucleus); arm_b_scores.tsv sha256 {sha(f'{W_}/arm_b_scores.tsv')}; sn_est_agonal.tsv {sha(f'{W_}/sn_est_agonal.tsv')}; axis_v2_controls.tsv.gz {sha(f'{W_}/axis_v2_controls.tsv.gz')}; snap_weights.tsv {sha(f'{W_}/scripts/snap_weights.tsv')}\n"
       "# partial Spearman / slope adjusted for age and sex (and diagnosis in group 'all'); slope = SNAP score (SD) per SD of the covariate; ph_axis_score: higher = higher-pH-like; estimated agonal: higher = more agonal; 2,000 donor bootstraps (stratified by diagnosis)\n")
with open("snap_vs_ph.tsv", "w") as fh: fh.write(hdr); R.to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 220); print(R[R.statistic.isin(["spearman_rho", "slope_SD_per_SD"])].to_string())
print(D.groupby("set").agg(n=("sz", "size"), n_ph=("ph", lambda s: s.notna().sum()), n_ctl_ph=("ph", lambda s: s[D.loc[s.index, "sz"] == 0].notna().sum())))
