#!/usr/bin/env python3
"""SNAP split into SNAP-a (astrocyte) and SNAP-n (neuron) parts: association with the estimated agonal score and the pH-axis score, with measured pH (bulk controls), and the schizophrenia effect under four covariate models. Frink.
DIRECTIONS. Scores: Ling's fixed LF4 loadings (snap_weights.tsv: w_total, w_astro, w_neuron) on within-dataset gene z; score = sum(w z)/sum|w|, then divided by its SD; lower SNAP = the direction Ling reports for schizophrenia and ageing.
 estimated agonal score: weights = GTEx slow-minus-fast death t (agonal_axis.tsv, mean of Cortex, BA9, BA24), SNAP's top 1,000 genes removed; HIGHER = more agonal.
 pH-axis score: weights = control pH axis v2 Stouffer z (axis_v2_controls.tsv.gz; z > 0 = gene HIGHER at higher pH), SNAP's top 1,000 genes removed; HIGHER = more like higher-pH tissue, LOWER = more like low-pH tissue. Measured pH: higher = higher pH.
 So a NEGATIVE correlation of a SNAP part with the pH-axis score or with measured pH means the part is HIGHER at LOW pH; a positive one means LOWER at low pH. A negative correlation with the agonal score means LOWER when more agonal.
 SZ effect = schizophrenia minus control in SD of the score (negative = lower in schizophrenia).
MODES: sn (single-nucleus donors, banks HBCC_main, MSSM_main, NSW pooled; CMC and SZBD shown only where scores exist) and bulk (22 brain-bank datasets of out/bulk, pooled over cohort clusters with snaplib.pool_clusters: equal-weight mean of cluster means, 2,000 bootstraps over clusters).
Usage: snap_components.py sn|bulk WORKDIR   (WORKDIR holds arm_b_scores.tsv, sn_est_agonal.tsv, axis_v2_controls.tsv.gz, scripts -> snap_what_is_it scripts, bulk -> snap_what_is_it out/bulk)
"""
import os, sys, hashlib
import numpy as np, pandas as pd
from scipy import stats, sparse as sp
MODE, W_ = sys.argv[1], sys.argv[2]; sys.path.insert(0, f"{W_}/scripts"); sys.dont_write_bytecode = True
from snaplib import zcols, partial_spearman, pool_clusters, pct, read_commented
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
rng = np.random.default_rng(20261023)
SCORES = {"SNAP total": "total", "SNAP-a": "astro", "SNAP-n": "neuron"}
Wt = read_commented(f"{W_}/scripts/snap_weights.tsv", index_col=0); AGX = read_commented(f"{W_}/scripts/agonal_axis.tsv", index_col=0).mean(axis=1, skipna=True).dropna()
PHX = read_commented(f"{W_}/axis_v2_controls.tsv.gz", index_col=0)["z_stouffer"].dropna()
top1000 = set(Wt.w_total.abs().sort_values(ascending=False).index[:1000])
AGX1 = AGX[~AGX.index.isin(top1000)]; PHX1 = PHX[~PHX.index.isin(top1000)]
def zs(x): x = np.asarray(x, float); return (x - x.mean()) / x.std(ddof=1)
def ols(y, cols, j=1):
    X = np.column_stack([np.ones(len(y))] + cols); b = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ b; df = len(y) - np.linalg.matrix_rank(X)
    s2 = (r ** 2).sum() / df; V = s2 * np.linalg.pinv(X.T @ X); return b[j], np.sqrt(V[j, j]), df
def sz_models(y, sz, cov_base, ag, ph):
    """SZ coefficient (SD of score) under the four models"""
    ys = zs(y); out = {}
    for mn, extra in (("base", []), ("+agonal", [zs(ag)]), ("+pH axis", [zs(ph)]), ("+agonal +pH axis", [zs(ag), zs(ph)])):
        b, se, df = ols(ys, [sz.astype(float)] + cov_base + extra); out[mn] = (b, se, df)
    return out
def fisher_pool(rs, ns):
    z = np.arctanh(np.clip(rs, -.999, .999)); w = np.array(ns) - 3; return np.tanh((w * z).sum() / w.sum())
rows = []
if MODE == "sn":
    from snap_psychad import z_matrix
    PB = os.path.expanduser("~/snap_sz/pb"); KEY = {"HBCC_main": "psychad_hbcc", "MSSM_main": "psychad_mssm", "NSW": "gse254569", "CMC": "cmc"}
    S = read_commented(f"{W_}/arm_b_scores.tsv", dtype={"donor_id": str}); EA = read_commented(f"{W_}/sn_est_agonal.tsv", dtype={"donor_id": str})
    frames = []
    for st, key in KEY.items():
        m = sp.load_npz(f"{PB}/{key}/counts_class.npz").tocsc(); genes = pd.read_csv(f"{PB}/{key}/genes.tsv", sep="\t", dtype=str); cols = pd.read_csv(f"{PB}/{key}/cols_class.tsv", sep="\t", dtype={"donor_id": str, "bioassay_id": str})
        donors = [d for d in S[S.set == st].donor_id if d in set(cols.donor_id)]; num = np.zeros(len(donors)); den = np.zeros(len(donors))
        for c in ["glutamatergic", "gabaergic", "astrocyte"]:
            Z, _ = z_matrix(m, genes, cols, donors, c, PHX1.index.values); tt = PHX1.reindex(Z.columns).values
            ok = ~np.isnan(Z.values).all(1); num += np.where(ok, np.nan_to_num(Z.values) @ tt, 0); den += np.where(ok, np.abs(tt).sum(), 0)
        frames.append(pd.DataFrame({"set": st, "donor_id": donors, "ph_axis_score": np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)}))
    D = S[S.set.isin(KEY)].merge(EA[["set", "donor_id", "est_agonal_ex1000"]], on=["set", "donor_id"]).merge(pd.concat(frames), on=["set", "donor_id"])
    D = D.rename(columns={"est_agonal_ex1000": "ag", "ph_axis_score": "phx", "lognuc": "lognuc"})
    D.to_csv("sn_donor_table.tsv", sep="\t", index=False)
    POOL = ["HBCC_main", "MSSM_main", "NSW"]
    # --- correlations
    for sc in SCORES:
        for tgt, lab in (("ag", "estimated agonal score"), ("phx", "pH-axis score")):
            per = {}
            for st in KEY:
                d = D[D.set == st].dropna(subset=[sc, tgt, "age", "male", "sz"]).reset_index(drop=True)
                if len(d) < 20: continue
                r_plain = stats.spearmanr(d[sc], d[tgt])[0]; r_part = partial_spearman(d[sc].values, d[tgt].values, [d.age.values, d.male.values.astype(float), d.sz.values.astype(float)])
                bs = []
                for _ in range(1000):
                    i = np.concatenate([rng.choice(np.where(d.sz == v)[0], (d.sz == v).sum()) for v in (0, 1)]); bs.append((stats.spearmanr(d[sc].values[i], d[tgt].values[i])[0], partial_spearman(d[sc].values[i], d[tgt].values[i], [d.age.values[i], d.male.values[i].astype(float), d.sz.values[i].astype(float)])))
                bs = np.array(bs); per[st] = (d, r_plain, r_part)
                rows.append(dict(table="correlation", set=st, score=sc, against=lab, statistic="spearman", estimate=r_plain, lo95=np.percentile(bs[:, 0], 2.5), hi95=np.percentile(bs[:, 0], 97.5), n=len(d)))
                rows.append(dict(table="correlation", set=st, score=sc, against=lab, statistic="partial_spearman_age_sex_dx", estimate=r_part, lo95=np.percentile(bs[:, 1], 2.5), hi95=np.percentile(bs[:, 1], 97.5), n=len(d)))
            # pooled over banks (Fisher z, weights n-3; bootstrap within bank x diagnosis)
            pb = [k for k in POOL if k in per]
            for k_i, (stat) in enumerate(("spearman", "partial_spearman_age_sex_dx")):
                est = fisher_pool([per[k][1 + k_i] for k in pb], [len(per[k][0]) for k in pb]); bs = []
                for _ in range(1000):
                    rr = []; nn = []
                    for k in pb:
                        d = per[k][0]; i = np.concatenate([rng.choice(np.where(d.sz == v)[0], (d.sz == v).sum()) for v in (0, 1)])
                        rr.append(stats.spearmanr(d[sc].values[i], d[tgt].values[i])[0] if k_i == 0 else partial_spearman(d[sc].values[i], d[tgt].values[i], [d.age.values[i], d.male.values[i].astype(float), d.sz.values[i].astype(float)])); nn.append(len(i))
                    bs.append(fisher_pool(rr, nn))
                rows.append(dict(table="correlation", set="pooled HBCC_main+MSSM_main+NSW", score=sc, against=lab, statistic=stat, estimate=est, lo95=np.percentile(bs, 2.5), hi95=np.percentile(bs, 97.5), n=sum(len(per[k][0]) for k in pb)))
    # --- SZ effect
    for cs, cvn in (("age+sex", ["age", "male"]), ("age+sex+PMI+log nuclei (Arm B primary)", ["age", "male", "pmi", "lognuc"])):
        for sc in SCORES:
            per = {}
            for st in KEY:
                d = D[D.set == st].dropna(subset=[sc, "ag", "phx", "sz"] + cvn).reset_index(drop=True)
                if len(d) < 30 or d.sz.nunique() < 2: continue
                res = sz_models(d[sc].values, d.sz.values, [zs(d[c].values) if c != "male" else d[c].values.astype(float) for c in cvn], d.ag.values, d.phx.values)
                per[st] = (res, len(d), int(d.sz.sum()))
                for mn, (b, se, df) in res.items(): rows.append(dict(table="sz_effect", set=st, score=sc, covariates=cs, model=mn, estimate=b, lo95=b - stats.t.ppf(.975, df) * se, hi95=b + stats.t.ppf(.975, df) * se, n=len(d), n_sz=int(d.sz.sum())))
            pb = [k for k in POOL if k in per]
            for mn in ("base", "+agonal", "+pH axis", "+agonal +pH axis"):
                b = np.array([per[k][0][mn][0] for k in pb]); se = np.array([per[k][0][mn][1] for k in pb]); w = 1 / se ** 2; mu = (w * b).sum() / w.sum(); sem = np.sqrt(1 / w.sum())
                rows.append(dict(table="sz_effect", set="pooled HBCC_main+MSSM_main+NSW (inverse variance)", score=sc, covariates=cs, model=mn, estimate=mu, lo95=mu - 1.96 * sem, hi95=mu + 1.96 * sem, n=sum(per[k][1] for k in pb), n_sz=sum(per[k][2] for k in pb)))
else:
    BD = f"{W_}/bulk"; names = [l.strip() for l in open(f"{W_}/datasets22.txt") if l.strip()]
    IN_AX = {"GSE17612", "GSE202537_caudate", "GSE208338", "GSE21138", "GSE224683", "GSE35977", "GSE53987_striatum", "GSE80655_AnCg", "GSE92538_GPL570", "GSE92538_GPL96"}
    per = []
    for nm in names:
        E = pd.read_csv(f"{BD}/{nm}_expr.tsv.gz", sep="\t", index_col=0); cv = pd.read_csv(f"{BD}/{nm}_cv.tsv", sep="\t").set_index("gsm").loc[E.columns]
        genes = E.index.values; Z = zcols(E.T.values.astype(float)); cluster = str(cv.cluster.iloc[0])
        sc = {}
        for lab, col in SCORES.items():
            w = Wt[f"w_{col}"].reindex(genes).fillna(0).values; s = Z @ w / np.abs(w).sum(); sc[lab] = s / s.std(ddof=1)
        def axscore(ax):
            g = [i for i, x in enumerate(genes) if x in ax.index]; t = ax.reindex(genes[g]).values; return (Z[:, g] @ t) / np.abs(t).sum()
        ag, phx = axscore(AGX1), axscore(PHX1)
        age = pd.to_numeric(cv.age, errors="coerce").values; male = pd.to_numeric(cv.male, errors="coerce").values; sz = pd.to_numeric(cv.sz, errors="coerce").values; ph = pd.to_numeric(cv.pH, errors="coerce").values
        both = np.nansum(sz == 1) >= 5 and np.nansum(sz == 0) >= 5
        ctl = (sz == 0) & ~np.isnan(age) & ~np.isnan(male)
        for lab in SCORES:
            y = sc[lab]
            # measured pH in controls
            m = ctl & ~np.isnan(ph)
            if m.sum() >= 20:
                rows.append(dict(table="bulk_by_dataset", dataset=nm, cluster=cluster, in_axis_construction=nm in IN_AX, score=lab, analysis="controls: measured pH (higher = higher pH), partial Spearman given age, sex", estimate=partial_spearman(y[m], ph[m], [age[m], male[m]]), n=int(m.sum())))
                rows.append(dict(table="bulk_by_dataset", dataset=nm, cluster=cluster, in_axis_construction=nm in IN_AX, score=lab, analysis="controls: measured pH, Spearman", estimate=stats.spearmanr(y[m], ph[m])[0], n=int(m.sum())))
            for tgt, tl in ((ag, "estimated agonal score"), (phx, "pH-axis score")):
                m2 = ~np.isnan(age) & ~np.isnan(male) & (~np.isnan(sz) if both else True)
                cov = [age[m2], male[m2]] + ([sz[m2]] if both else [])
                rows.append(dict(table="bulk_by_dataset", dataset=nm, cluster=cluster, in_axis_construction=nm in IN_AX, score=lab, analysis=f"all samples: {tl}, partial Spearman given age, sex" + (", dx" if both else ""), estimate=partial_spearman(y[m2], tgt[m2], cov), n=int(m2.sum())))
            if both:
                m3 = ~np.isnan(sz) & ~np.isnan(age) & ~np.isnan(male)
                res = sz_models(y[m3], sz[m3], [zs(age[m3]), male[m3]], ag[m3], phx[m3])
                for mn, (b, se, df) in res.items():
                    rows.append(dict(table="bulk_by_dataset", dataset=nm, cluster=cluster, in_axis_construction=nm in IN_AX, score=lab, analysis=f"SZ effect, model {mn} (age, sex)", estimate=b, lo95=b - stats.t.ppf(.975, df) * se, hi95=b + stats.t.ppf(.975, df) * se, n=int(m3.sum()), n_sz=int((sz[m3] == 1).sum())))
    R0 = pd.DataFrame(rows); pooled = []
    for (sc_, an), g in R0.groupby(["score", "analysis"]):
        for flt, fl in (("all datasets", g), ("datasets not used to build the pH axis", g[~g.in_axis_construction])):
            if len(fl) < 2: continue
            p, lo, hi, nd, nc = pool_clusters(fl.estimate.values, fl.cluster.values, rng)
            pooled.append(dict(table="bulk_pooled", dataset=f"pooled over {nc} cohort clusters ({nd} datasets); {flt}", score=sc_, analysis=an, estimate=p, lo95=lo, hi95=hi, n=int(fl.n.sum())))
    rows = rows + pooled
R = pd.DataFrame(rows)
for c in R.select_dtypes("float"): R[c] = R[c].round(4)
src = f"arm_b_scores.tsv {sha(f'{W_}/arm_b_scores.tsv')}; " if MODE == "sn" else ""
hdr = (f"# snap_components.py {MODE} (snap_up_down); {src}axis_v2_controls.tsv.gz {sha(f'{W_}/axis_v2_controls.tsv.gz')}; snap_weights.tsv {sha(f'{W_}/scripts/snap_weights.tsv')}; agonal_axis.tsv {sha(f'{W_}/scripts/agonal_axis.tsv')}\n"
       "# DIRECTIONS: SNAP-a = astrocyte part, SNAP-n = neuron part of Ling's LF4 loadings; lower SNAP = direction of schizophrenia and ageing. Estimated agonal score: HIGHER = more agonal. pH-axis score: HIGHER = higher-pH-like (gene weights z > 0 = higher at higher pH), so NEGATIVE correlation = part is HIGHER at LOW pH. Measured pH: higher = higher pH. SZ effect in SD of the score (negative = lower in schizophrenia).\n"
       "# pooling: single-nucleus = Fisher z weights n-3 with bootstrap within bank x diagnosis (correlations), inverse-variance fixed effect (SZ effect); bulk = equal-weight mean of cluster means, bootstrap over clusters (2,000)\n")
with open(f"snap_components_{MODE}.tsv", "w") as fh: fh.write(hdr); R.to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
print(R[R.table.isin(["correlation", "sz_effect", "bulk_pooled"])].drop(columns=[c for c in ["cluster", "in_axis_construction"] if c in R]).to_string()[:12000])
