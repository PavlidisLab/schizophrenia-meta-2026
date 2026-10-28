#!/usr/bin/env python3
"""Coexpression coherence of each literature module in independent bulk datasets (frink).

Statistic = mean pairwise Pearson r among the module's genes (present genes, >= MINGENES and >= 50% of the module) minus the mean of the same statistic over NNULL random gene sets of the same size
matched gene-by-gene on expression decile (decile of the gene's mean expression in the dataset; module genes excluded from the pool). Unsigned: modules are coexpression modules, so a positive excess = the genes
move together more than expression-matched genes. Same estimator as snap_sz 03_bulk_main.py coherence() (mean pairwise r from column sums of standardized genes).
Variants: raw (SZ and control rows kept as in module_retest.py) | resid (expression residualized on diagnosis, age, sex and region where available, then re-standardized).
Interval: bootstrap over samples (NBOOT) of the observed statistic with the null mean held at the full-sample value; se = SD of the bootstrap.
Only (module, dataset) pairs independent of the module's discovery lineage are scored (rule of aggregate.py in module_retest_2026-10-23).
Usage: module_coherence.py RUN_DIR RETEST_DIR OUT_DIR
"""
import sys, os, hashlib, warnings
import numpy as np, pandas as pd
from scipy import sparse
from multiprocessing import Pool
warnings.filterwarnings("ignore")
RUN, RT, OUT = sys.argv[1:4]
sys.path.insert(0, os.path.join(RUN, "scripts")); sys.dont_write_bytecode = True
from lib import zscore
MINGENES, NNULL, NBOOT = 10, 200, 200
idx = pd.read_csv(f"{RT}/modules_index.tsv", sep="\t", comment="#", dtype=str)
lng = pd.read_csv(f"{RT}/modules_long.tsv", sep="\t", comment="#", dtype=str)
mg = lng.groupby("module_uid")["gene"].apply(lambda s: list(dict.fromkeys(s))).to_dict()
lin = pd.read_csv(f"{RT}/dataset_lineage.tsv", sep="\t", comment="#", dtype=str)
lin = lin[lin.compartment == "bulk"]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
mlin = {r.module_uid: set(r.lineage.split(";")) - {"unnamed"} for r in idx.itertuples()}
def stdcols(X):
    Xc = X - X.mean(1, keepdims=True); sd = Xc.std(1, ddof=1, keepdims=True); sd[sd == 0] = np.inf
    return Xc / sd
def stat(S, k, n): return ((S ** 2).sum(-1) / (n - 1) - k) / (k * (k - 1))
def run(row):
    key, cluster, dslin = row
    rng = np.random.default_rng(abs(hash(key)) % 2**32)
    z = np.load(f"{RUN}/out/ds/{key}.expr.npz", allow_pickle=True); genes = np.array(z["genes"], dtype=object); V = z["values"].astype(np.float64); samples = [str(s) for s in z["samples"]]
    cov = pd.read_csv(f"{RUN}/out/ds/{key}.cov.tsv", sep="\t", dtype={"sample": str}).set_index("sample").reindex(samples).reset_index()
    dxb = cov["dx_group"].map({"control": 0.0, "case": 1.0})
    szm = cov["dx_text"].astype(str).str.contains(r"schizo|\bSZ\b|\bSCZ\b|\bSCHIZ", case=False, regex=True, na=False)
    sz = pd.Series(np.nan, index=cov.index); sel = (dxb == 0) | ((dxb == 1) & szm); sz[sel] = ((dxb == 1) & szm)[sel].astype(float)
    outl = cov["outlier"].astype(str).str.lower().isin(["true", "1", "1.0"]) if "outlier" in cov else pd.Series(False, index=cov.index)
    keep = (sz.notna() & ~outl).values
    if (sz[keep] == 1).sum() < 3 or (sz[keep] == 0).sum() < 3: return []
    ok, ok2, X, mu, sd = zscore(V); g_ok = genes[ok][ok2]; X = X[ok2]; mu = mu[ok2]; sd = sd[ok2]
    X = np.where(np.isnan(X), mu[:, None], X)[:, keep]; c = cov[keep].reset_index(drop=True); y = sz[keep].values; n = X.shape[1]
    gi = {g: i for i, g in enumerate(g_ok)}
    # expression deciles of the raw mean (mu is the per-gene mean on the dataset scale)
    dec = np.minimum((pd.Series(mu).rank(pct=True).values * 10).astype(int), 9)
    cols = [np.ones(n), y]
    age = pd.to_numeric(c["age"], errors="coerce").values
    if np.isfinite(age).mean() >= 0.8 and np.nanstd(age) > 0: cols.append(np.where(np.isfinite(age), age, np.nanmean(age)))
    sx = c["sex"].map({"M": 1.0, "F": 0.0}).values
    if np.isfinite(sx).mean() >= 0.8 and np.nanstd(sx) > 0: cols.append(np.where(np.isfinite(sx), sx, np.nanmean(sx)))
    r = c["region"].astype(object).where(c["region"].notna(), None) if "region" in c else pd.Series([None] * n)
    if r.nunique() > 1:
        lv = r.value_counts(); r = r.where(r.isin(lv.index[lv >= 3]), "other_small")
        if r.nunique() > 1: cols.extend(list(pd.get_dummies(r, drop_first=True).values.astype(float).T))
    D = np.column_stack(cols); Xres = X - (X @ D @ np.linalg.pinv(D.T @ D)) @ D.T if np.linalg.matrix_rank(D) == D.shape[1] else None
    variants = {"raw": X}
    if Xres is not None: variants["resid"] = Xres
    G = len(g_ok); out = []
    # modules independent of this dataset
    mods = []
    for r_ in idx.itertuples():
        if mlin[r_.module_uid] & (set(dslin.split(";")) - {"own"}): continue
        gl = [g for g in mg[r_.module_uid] if g in gi]
        if len(gl) < MINGENES or len(gl) < 0.5 * len(mg[r_.module_uid]): continue
        mods.append((r_.module_uid, np.array([gi[g] for g in gl])))
    for vn, Xv in variants.items():
        Z = stdcols(Xv); Zb_list = [stdcols(Xv[:, rng.integers(0, n, n)]) for _ in range(NBOOT)]
        for uid, I in mods:
            k = len(I); obs = float(stat(Z[I].sum(0), k, n))
            M = np.empty((NNULL, k), dtype=np.int64); excl = np.zeros(G, bool); excl[I] = True
            for d in range(10):
                pos = np.where(dec[I] == d)[0]
                if len(pos) == 0: continue
                pool = np.where((dec == d) & ~excl)[0]
                if len(pool) < len(pos): pool = np.where(dec == d)[0]
                kk = rng.random((NNULL, len(pool))); M[:, pos] = pool[np.argpartition(kk, len(pos) - 1, axis=1)[:, :len(pos)]] if len(pos) < len(pool) else pool[np.argsort(kk, axis=1)]
            ind = sparse.csr_matrix((np.ones(NNULL * k, np.float32), M.ravel(), np.arange(0, NNULL * k + 1, k)), shape=(NNULL, G))
            nul = stat(np.asarray(ind @ Z.astype(np.float32)), k, n); nm = float(nul.mean())
            bs = np.array([float(stat(Zb[I].sum(0), k, n)) - nm for Zb in Zb_list])
            out.append(dict(key=key, cluster=cluster, variant=vn, module_uid=uid, n=n, n_sz=int((y == 1).sum()), n_ctl=int((y == 0).sum()), genes_present=k, observed=obs, null_mean=nm, null_sd=float(nul.std(ddof=1)),
                            excess=obs - nm, se=float(bs.std(ddof=1)), z_vs_null=(obs - nm) / float(nul.std(ddof=1))))
    print(key, len(out), flush=True)
    return out
if __name__ == "__main__":
    rows = [(r.key, r.cluster, r.lineage) for r in lin.itertuples() if os.path.exists(f"{RUN}/out/ds/{r.key}.expr.npz")]
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(run, rows, chunksize=1)
    R = pd.DataFrame([x for r in res for x in r])
    for c in ("observed", "null_mean", "null_sd", "excess", "se", "z_vs_null"): R[c] = R[c].round(5)
    os.makedirs(OUT, exist_ok=True)
    hdr = (f"# module_coherence.py (module_coherence); MINGENES {MINGENES}; NNULL {NNULL}; NBOOT {NBOOT}; modules_index.tsv {sha(f'{RT}/modules_index.tsv')}; modules_long.tsv {sha(f'{RT}/modules_long.tsv')}; dataset_lineage.tsv {sha(f'{RT}/dataset_lineage.tsv')}\n"
           "# excess = mean pairwise r among module genes minus mean over expression-decile-matched random sets (positive = more coherent than matched genes); only datasets independent of the module's discovery lineage; se = bootstrap SD over samples\n")
    with open(f"{OUT}/per_dataset_coherence.tsv", "w") as fh: fh.write(hdr); R.to_csv(fh, sep="\t", index=False)
    print(len(R), R.key.nunique())
