#!/usr/bin/env python3
"""Per-dataset schizophrenia effect on the score of each literature coexpression module (direct test; frink).

Dataset packages: artifact_sample_scores_run_2026-10-18/out/ds/<key>.expr.npz + .cov.tsv (Gemma processed data, HGNC symbols, log2-type scale) and out/sample_scores.tsv.gz (per-sample pH_axis and agonal_axis scores).
Module score of a sample = mean of the within-dataset gene z over the module genes present (>= MINGENES genes and >= 50% of the module's genes; genes with > 20% missing or zero variance dropped, as in 04_score_and_stats.py),
multiplied by the module's stated direction (+1 higher, -1 lower in SZ) so that a positive effect = the paper's direction. Effect = SZ minus control in SD of the score.
SZ versus control rows only (dx control, or case with 'schizo' in the diagnosis text; Gemma outlier samples dropped). Models (OLS; SE and t interval from the residual df):
  m0 diagnosis only | m1 + age, sex (+ region when a dataset has several) | m2 = m1 + PMI + pH (measured, only where >= 80% of the rows have each) | m3 = m1 + pH_axis + agonal_axis sample scores (mean z weighted by the control pH axis and the GTEx agonal axis).
Usage: module_retest.py RUN_DIR MODULES_DIR DATASET_LINEAGE OUT
"""
import sys, os, re, json, hashlib, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
RUN, MODDIR, LIN, OUT = sys.argv[1:5]
sys.path.insert(0, os.path.join(RUN, "scripts")); sys.dont_write_bytecode = True
from lib import zscore, ols_beta
MINGENES = int(os.environ.get("MINGENES", 10))
idx = pd.read_csv(f"{MODDIR}/modules_index.tsv", sep="\t", comment="#", dtype=str)
lng = pd.read_csv(f"{MODDIR}/modules_long.tsv", sep="\t", comment="#", dtype=str)
mg = lng.groupby("module_uid")["gene"].apply(list).to_dict()
lin = pd.read_csv(LIN, sep="\t", comment="#", dtype=str)
SS = pd.read_csv(f"{RUN}/out/sample_scores.tsv.gz", sep="\t", usecols=["key", "sample", "pH_axis", "agonal_axis"], dtype={"sample": str})
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
rows = []
for key in lin.key:
    if not os.path.exists(f"{RUN}/out/ds/{key}.expr.npz"): print("no dataset package:", key); continue
    z = np.load(f"{RUN}/out/ds/{key}.expr.npz", allow_pickle=True); genes = np.array(z["genes"], dtype=object); V = z["values"].astype(np.float64); samples = [str(s) for s in z["samples"]]
    cov = pd.read_csv(f"{RUN}/out/ds/{key}.cov.tsv", sep="\t", dtype={"sample": str}).set_index("sample").reindex(samples).reset_index()
    cov = cov.merge(SS[SS.key == key][["sample", "pH_axis", "agonal_axis"]], on="sample", how="left")
    dxb = cov["dx_group"].map({"control": 0.0, "case": 1.0})
    szm = cov["dx_text"].astype(str).str.contains(r"schizo|\bSZ\b|\bSCZ\b|\bSCHIZ", case=False, regex=True, na=False)
    sz = pd.Series(np.nan, index=cov.index); sel = (dxb == 0) | ((dxb == 1) & szm); sz[sel] = ((dxb == 1) & szm)[sel].astype(float)
    outl = cov["outlier"].astype(str).str.lower().isin(["true", "1", "1.0"]) if "outlier" in cov else pd.Series(False, index=cov.index)
    keep = sz.notna() & ~outl
    if (sz[keep] == 1).sum() < 3 or (sz[keep] == 0).sum() < 3: continue
    ok, ok2, X, mu, sd = zscore(V); g_ok = genes[ok][ok2]; X = X[ok2]; mu = mu[ok2]; sd = sd[ok2]
    X = np.where(np.isnan(X), mu[:, None], X); Z = (X - mu[:, None]) / sd[:, None]
    gi = {g: i for i, g in enumerate(g_ok)}
    c = cov[keep].reset_index(drop=True); Zk = Z[:, keep.values]; y_sz = sz[keep].values
    for col in ("age", "pmi_h", "ph"): c[col] = pd.to_numeric(c[col], errors="coerce")
    c["sex_bin"] = c["sex"].map({"M": 1.0, "F": 0.0})
    n = len(c)
    def design(extra):
        cols = [np.ones(n), y_sz]; use = np.ones(n, bool); parts = []
        for nm in extra:
            if nm == "region":
                r = c["region"].astype(object).where(c["region"].notna(), None)
                if r.nunique() > 1:
                    lv = r.value_counts(); r = r.where(r.isin(lv.index[lv >= 3]), "other_small")
                    if r.nunique() > 1: parts.append(("region", pd.get_dummies(r, drop_first=True).values.astype(float)))
                continue
            v = c[nm].values.astype(float)
            if np.isfinite(v).mean() < 0.8 or np.nanstd(v) == 0: return None
            use &= np.isfinite(v); parts.append((nm, v))
        return use, parts
    models = {"m0": [], "m1": ["age", "sex_bin", "region"], "m2": ["age", "sex_bin", "region", "pmi_h", "ph"], "m3": ["age", "sex_bin", "region", "pH_axis", "agonal_axis"]}
    ds_models = {}
    for mn, extra in models.items():
        d = design(extra)
        if d is None: continue
        use, parts = d
        if mn == "m1" and not all(nm in [p[0] for p in parts] for nm in ("age", "sex_bin")): continue
        if mn == "m2" and not all(nm in [p[0] for p in parts] for nm in ("pmi_h", "ph")): continue
        if mn == "m3" and not all(nm in [p[0] for p in parts] for nm in ("pH_axis", "agonal_axis")): continue
        if (y_sz[use] == 1).sum() < 3 or (y_sz[use] == 0).sum() < 3: continue
        cols = [np.ones(use.sum()), y_sz[use]]
        for nm, v in parts:
            v = v[use]
            if v.ndim == 1: cols.append((v - v.mean()) / (v.std(ddof=1) or 1))
            else: cols.extend(list(v.T))
        Xd = np.column_stack(cols)
        if np.linalg.matrix_rank(Xd) < Xd.shape[1]: continue
        ds_models[mn] = (use, Xd)
    for r in idx.itertuples():
        gl = [g for g in dict.fromkeys(mg[r.module_uid]) if g in gi]
        if len(gl) < MINGENES or len(gl) < 0.5 * len(set(mg[r.module_uid])): continue
        S = Zk[[gi[g] for g in gl]].mean(0) * int(r.direction)
        for mn, (use, Xd) in ds_models.items():
            Y = S[use][:, None]; sdv = Y.std(0, ddof=1)
            if sdv[0] == 0: continue
            Ys = (Y - Y.mean(0)) / sdv
            b, se, df = ols_beta(Ys, Xd, 1)
            rows.append(dict(key=key, module_uid=r.module_uid, model=mn, n_sz=int((y_sz[use] == 1).sum()), n_ctl=int((y_sz[use] == 0).sum()), genes_present=len(gl), beta=float(b[0]), se=float(se[0]), df=int(df)))
R = pd.DataFrame(rows)
R["beta"] = R.beta.round(4); R["se"] = R.se.round(4)
hdr = f"# module_retest.py (module_retest); MINGENES {MINGENES}; modules_index.tsv sha256 {sha(f'{MODDIR}/modules_index.tsv')}; modules_long.tsv {sha(f'{MODDIR}/modules_long.tsv')}; sample_scores.tsv.gz {sha(f'{RUN}/out/sample_scores.tsv.gz')}\n# beta = SZ minus control in SD of the direction-signed module score (positive = the paper's direction); models in the script docstring\n"
with open(OUT, "w") as fh: fh.write(hdr); R.to_csv(fh, sep="\t", index=False)
print(len(R), R.groupby("model").size().to_dict(), R.key.nunique())
