#!/usr/bin/env python3
"""Estimated agonal score per single-nucleus donor.  Run on remote with snap_sz/envs/py/bin/python (single thread, nice); needs snap_psychad.py (side-studies/coexpression_claims_survey/snap_psychad/, copied beside this script),
agonal_axis.tsv, snap_weights.tsv and snap_sz/pb/<key>/ (build_pseudobulk.py outputs).
Per donor and SNAP cell class (glutamatergic, gabaergic, astrocyte): log2 CPM (prior 0.5) of summed counts, genes with mean CPM >= 1, z across the donor set (z_matrix of snap_psychad.py); estimated agonal score = sum(t z)/sum(|t|)
over the axis genes (t = GTEx Hardy 4 vs 2 t, mean over Cortex, BA9, BA24; positive = slow death), numerator and denominator summed over the three classes; higher = more agonal (slow-death-like).  Variants: SNAP's top 1,000 / top 200 genes by |w_total| removed from the axis; none removed."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from snap_psychad import z_matrix
H = os.path.dirname(os.path.abspath(__file__)); PB = os.path.expanduser("snap_sz/pb")
SETS = {"CMC": "cmc", "NSW": "gse254569", "HBCC_main": "psychad_hbcc", "MSSM_main": "psychad_mssm"}
def rc(p, **k): return pd.read_csv(p, sep="\t", comment="#", **k)
W = rc(f"{H}/snap_weights.tsv", index_col=0); AX = rc(f"{H}/agonal_axis.tsv", index_col=0); t = AX.mean(axis=1, skipna=True).dropna()
top = W.w_total.abs().sort_values(ascending=False).index
D = rc(f"{H}/arm_b_donors.tsv", dtype={"donor_id": str}); out = []
import scipy.sparse as sp
for sname, key in SETS.items():
    m = sp.load_npz(f"{PB}/{key}/counts_class.npz").tocsc(); genes = pd.read_csv(f"{PB}/{key}/genes.tsv", sep="\t", dtype=str); cols = pd.read_csv(f"{PB}/{key}/cols_class.tsv", sep="\t", dtype={"donor_id": str, "bioassay_id": str})
    donors = list(D[D.set == sname].donor_id); donors = [d for d in donors if d in set(cols.donor_id)]
    res = {}
    for lab, drop in [("ex1000", set(top[:1000])), ("ex200", set(top[:200])), ("none", set())]:
        ax = t[~t.index.isin(drop)]; num = np.zeros(len(donors)); den = np.zeros(len(donors))
        for c in ["glutamatergic", "gabaergic", "astrocyte"]:
            Z, _ = z_matrix(m, genes, cols, donors, c, ax.index.values); tt = ax.reindex(Z.columns).values
            ok = ~np.isnan(Z.values).all(1); num += np.where(ok, np.nan_to_num(Z.values) @ tt, 0); den += np.where(ok, np.abs(tt).sum(), 0)
        res["est_agonal_" + lab] = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    out.append(pd.DataFrame({"set": sname, "donor_id": donors, **res})); print(sname, len(donors), flush=True)
pd.concat(out).to_csv(os.path.expanduser("snap_signal_analysis/out/sn_est_agonal.tsv"), sep="\t", index=False)
