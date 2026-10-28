#!/usr/bin/env python3
"""STAT-012 addendum 2, item 3b: Extended Fig. 9c. Per GTEx v8 tissue, Spearman rho between the slow-death (Hardy 4) against fast-death (Hardy 2) t statistics of the
recurrent genes and their schizophrenia direction (+1 up, -1 down), with the P value of the Spearman correlation test (the test that treats the genes as independent)
and a sample-permutation null that keeps the correlation among genes: the death-group label is permuted among the donors of the two groups in that tissue
(ischemic time, age and sex stay with the donor), t is recomputed for the recurrent genes, rho recomputed.
Sample selection, gene filter (median TPM >= 1 over all the tissue's samples), design and OLS t follow analyses/gtex_death_speed_all_tissues/gtex_death_speed_all_tissues.py
(re-implemented here for the recurrent-gene rows only; the observed rho is checked against figure_data/gtex_death_speed_all_tissues/rho_by_tissue_and_contrast.tsv).
Gene set: recurrent_sets.tsv set study_brain (209 genes), sz_direction up/down.
    python3 gtex_signdir_perm.py TPM.gct.gz ATTR.txt SUBJ.txt recurrent_sets.tsv OUTDIR [NPERM]"""
import sys, io, hashlib, datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
TPM, ATTR, SUBJ, RSET, OUT = sys.argv[1:6]; NPERM = int(sys.argv[6]) if len(sys.argv) > 6 else 2000
SEED, MINN = 20261012, 10
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
rs = pd.read_csv(RSET, sep="\t", comment="#"); rs = rs[(rs.set == "study_brain") & rs.sz_direction.isin(["up", "down"])]
szdir = dict(zip(rs.gene, np.where(rs.sz_direction == "up", 1.0, -1.0))); genes = set(szdir)
attr = pd.read_csv(ATTR, sep="\t", usecols=["SAMPID", "SMTSD", "SMTS", "SMAFRZE", "SMTSISCH"]); subj = pd.read_csv(SUBJ, sep="\t")
s = attr[(attr.SMAFRZE == "RNASEQ") & attr.SMTSISCH.notna()].copy()
s["SUBJID"] = s.SAMPID.str.split("-").str[:2].str.join("-"); s = s.merge(subj, on="SUBJID")
s["age"] = s.AGE.map(lambda a: np.mean([int(x) for x in a.split("-")])); s["isch_h"] = s.SMTSISCH / 60.0; s["male"] = (s.SEX == 1).astype(float)
hdr = pd.read_csv(TPM, sep="\t", skiprows=2, nrows=0).columns
s = s[s.SAMPID.isin(set(hdr))].reset_index(drop=True)
parts = []
for ch in pd.read_csv(TPM, sep="\t", skiprows=2, usecols=["Name", "Description"] + list(s.SAMPID), dtype={c: np.float32 for c in s.SAMPID}, chunksize=2000):
    parts.append(ch[ch.Description.isin(genes)])
t0 = pd.concat(parts); sym = t0.Description.values; Xraw = t0[list(s.SAMPID)].values.astype(np.float64); X = np.log2(Xraw + 1)
print("samples", len(s), "rows kept", len(sym), "genes", len(set(sym)), flush=True)
rng = np.random.default_rng(SEED)
def tstat(Y, M):
    Q, R = np.linalg.qr(M); beta = np.linalg.solve(R, Q.T @ Y.T); res = Y.T - M @ beta
    s2 = (res ** 2).sum(0) / (M.shape[0] - M.shape[1]); cd = np.diag(np.linalg.inv(R.T @ R))[1]
    return beta[1] / np.sqrt(s2 * cd)
def spear(a, b):
    ra, rb = rankdata(a), rankdata(b); ra -= ra.mean(); rb -= rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))
rows = []
for t, sub in s.groupby("SMTSD"):
    ix = sub.index.values; h = sub.DTHHRDY.values
    if (h == 4).sum() < MINN or (h == 2).sum() < MINN: continue
    keep = np.median(Xraw[:, ix], axis=1) >= 1
    sel = np.isin(h, (2, 4)); lab = (h[sel] == 4).astype(float)
    covs = {"isch": sub.isch_h.values[sel], "age": sub.age.values[sel], "sex": sub.male.values[sel]}
    cv = [v for v in covs.values() if np.ptp(v) != 0]
    Y = X[:, ix[sel]][keep]
    names = sym[keep]
    M = np.column_stack([np.ones(len(lab)), lab] + cv)
    # first of duplicated symbols after the filter, as groupby(level=0).first()
    _, first = np.unique(names, return_index=True); first = np.sort(first)
    gs = names[first]; Yg = Y[first]; sd = np.array([szdir[g] for g in gs])
    tv = tstat(Yg, M); rho = spear(tv, sd); n = len(gs)
    tt = rho * np.sqrt((n - 2) / (1 - rho ** 2)); from scipy.stats import t as T; p_par = 2 * T.sf(abs(tt), n - 2)
    perm = np.empty(NPERM)
    rsd = rankdata(sd); rsd = rsd - rsd.mean(); ns = np.sqrt(rsd @ rsd)
    for i in range(NPERM):
        M[:, 1] = rng.permutation(lab); r = rankdata(tstat(Yg, M)); r = r - r.mean(); perm[i] = r @ rsd / (np.sqrt(r @ r) * ns)
    rows.append(dict(tissue=t, SMTS=sub.SMTS.iloc[0], n_slow=int(lab.sum()), n_fast=int((1 - lab).sum()), n_genes=n, rho_observed=rho, p_spearman_test=p_par,
                     perm_mean=perm.mean(), perm_sd=perm.std(ddof=1), sd_if_genes_independent=1 / np.sqrt(n - 1), sd_inflation=perm.std(ddof=1) * np.sqrt(n - 1),
                     perm_z=(rho - perm.mean()) / perm.std(ddof=1), p_perm_one_sided=(1 + (perm >= rho).sum()) / (NPERM + 1), p_perm_two_sided=(1 + (np.abs(perm) >= abs(rho)).sum()) / (NPERM + 1),
                     perm_q975=np.quantile(perm, .975)))
    print(t, n, round(rho, 4), round(p_par, 4), round(rows[-1]["perm_sd"], 3), rows[-1]["p_perm_two_sided"], flush=True)
R = pd.DataFrame(rows)
with open(f"{OUT}/gtex_signdir_perm.tsv", "w") as f:
    f.write(f"# claim: Extended Fig. 9c: Spearman rho between slow-against-fast death t and the SZ direction of the 209 recurrent genes, P < 0.05 in 28 of 32 non-brain tissues; sample-permutation null of the death-group label\n"
            f"# script: analyses/gene_correlation_adjusted_tests/marker_set_and_ph_permutation_tests/gtex_signdir_perm.py, {NPERM} permutations per tissue, seed {SEED}, run {datetime.date.today().isoformat()}\n"
            f"# source gtex_tpm {sha(TPM)}; attributes {sha(ATTR)}; subjects {sha(SUBJ)}; recurrent_sets.tsv {sha(RSET)}\n")
    R.to_csv(f, sep="\t", index=False)
