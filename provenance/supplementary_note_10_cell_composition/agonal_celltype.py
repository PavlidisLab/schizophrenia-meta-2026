#!/usr/bin/env python3
"""Why are genes pH-sensitive? Transcript features, RNA-quality (RIN) association and cell-type specificity against the control-only pH Z and the GTEx Hardy t.

    python3 feature_analysis.py [inputs_root]        inputs_root: analyses_laptop/analyses (default, the laptop's confounds results) or analyses (remote environment copies)

Targets (per gene): pH Z = control-only association with brain pH, Stouffer-combined over the independent collections (confounds/results/combined_controls_pH.tsv, factor pH; Z > 0 = RNA lower at low pH);
GTEx Hardy t = mean over BA9, BA24 and Cortex of the slow-versus-fast-death t (gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz; positive = higher in slow death); RIN Z = Stouffer (weights sqrt n) of the control-only RIN t of
GSE53987 BA46 / hippocampus / striatum and GSE208338 (confounds/results/within_controls_*.tsv, factor RIN; the model also holds pH and PMI). Expression: expr_percentile of the pH table (LIBD DLPFC mean log CPM), deciles.
Features: gene_features.tsv (build_features.py: MANE Select transcript lengths, GC, AU-rich elements, TargetScan miRNA sites) and mRNA half-life PC1 (Agarwal and Kelley 2022, analyses/ieg_halflife_2026-10-07/raw).
Cell type: PsychAD control specificity (EWCE: mean CPM in class / sum over classes) HBCC and MSSM, analyses/interneuron_ph_2026-10-01/results/specificity_*.tsv.gz.
Statistics. (1) Decile-adjusted Spearman: ranks of feature and target are centred within expression decile, then correlated (equivalent to the partial rank correlation on decile indicators); p from 10,000 permutations of the target
within decile; 95% interval from 1,000 bootstrap resamples of genes. (2) Marker sets: top 100 genes by specificity per class (class mean CPM >= 1) against 10,000 random sets of the same size matched on expression decile;
effect = mean pH Z of the set minus the matched null mean, with a 95% interval by resampling the set's genes. (3) Nested least-squares models of rank(pH Z) on decile indicators + feature ranks + specificity + RIN Z,
R^2 and increments with 500-bootstrap intervals. Seeds 20261009."""
import gzip, hashlib, os, sys, datetime, csv
import numpy as np, pandas as pd
from scipy.stats import rankdata
HERE = "."; SRC = "axes-handoff/analyses/ph_sensitivity_gene_features_2026-10-09"; ROOT = "axes-handoff"
IN = sys.argv[1] if len(sys.argv) > 1 else "analyses_laptop/analyses"; INR = os.path.join(ROOT, IN); TAG = "laptop" if "laptop" in IN else "remote environment"
AN = os.path.join(ROOT, "analyses"); OUT = os.path.join(HERE, f"results_{TAG}"); os.makedirs(OUT, exist_ok=True)
NPERM, NBOOT, SEED = 10000, 1000, 20261009
rng = np.random.default_rng(SEED)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
def rd(p, **k): return pd.read_csv(p, sep="\t", comment="#", **k)
phf = os.path.join(INR, "confounds/results/combined_controls_pH.tsv"); ph = rd(phf); ph = ph[ph.factor == "pH"].set_index("gene")
D = pd.DataFrame({"pH_Z": ph.Z, "expr_pct": ph.expr_percentile, "n_datasets": ph.n_datasets}).dropna(subset=["pH_Z", "expr_pct"])
D["decile"] = np.minimum((D.expr_pct * 10).astype(int), 9)
hf = os.path.join(AN, "gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"); H = pd.read_csv(hf, sep="\t", index_col=0)
D["gtex_hardy_t"] = H[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).reindex(D.index)
rin = []
for ds in ("GSE53987_BA46", "GSE53987_hippocampus", "GSE53987_striatum", "GSE208338"):
    w = rd(os.path.join(INR, f"confounds/results/within_controls_{ds}.tsv")); w = w[w.factor == "RIN"].set_index("gene"); rin.append((w.t, np.sqrt(w.n)))
num = sum((t * wt).reindex(D.index).fillna(0) for t, wt in rin); den = np.sqrt(sum((wt ** 2).reindex(D.index).fillna(0) for t, wt in rin)); cnt = sum(t.reindex(D.index).notna().astype(int) for t, _ in rin)
D["rin_Z"] = (num / den).where(cnt >= 2)
fe = rd(os.path.join(SRC, "gene_features.tsv")).set_index("gene"); FEATS = ["tx_len", "cds_len", "utr5_len", "utr3_len", "n_exons", "gc_tx", "gc_cds", "gc_utr3", "n_ARE_pentamer", "are_density_per_kb", "n_ARE_nonamer", "mirna_cons_sites", "mirna_site_density_per_kb"]
D = D.join(fe[FEATS], how="left")
hl = rd(os.path.join(AN, "ieg_halflife_2026-10-07/raw/all_HLs_human_PC1.txt")).dropna(subset=["external_gene_name"]).drop_duplicates("external_gene_name").set_index("external_gene_name")["human_PC1"]
D["halflife_PC1"] = hl.reindex(D.index); FEATS = FEATS + ["halflife_PC1"]
spec = {c: pd.read_csv(os.path.join(AN, f"interneuron_ph_2026-10-01/results/specificity_{c}.tsv.gz"), sep="\t", comment="#", index_col=0) for c in ("HBCC", "MSSM")}
S = spec["HBCC"]; D["spec_IN"] = S.spec_coarse_IN.reindex(D.index)
for c in ("Astro", "Oligo", "OPC", "Micro", "Vascular"): D[f"spec_{c}"] = S[f"spec_{c}"].reindex(D.index)
# ---- nested models, two targets on the same genes
ranks = lambda df: df.apply(lambda c: pd.Series(rankdata(c), index=c.index))
cols_feat = [f for f in FEATS if f != "halflife_PC1"]; cols_ct = ["spec_IN", "spec_Astro", "spec_Oligo", "spec_OPC", "spec_Micro", "spec_Vascular"]
C = D[["pH_Z", "gtex_hardy_t", "decile", "halflife_PC1"] + cols_feat + cols_ct].dropna()
dec = pd.get_dummies(C.decile, drop_first=True).astype(float).values
def r2(target, cols, ii=None):
    y = rankdata(C[target].values); X = np.column_stack([np.ones(len(C)), dec] + [rankdata(C[c].values) for c in cols])
    if ii is not None: y, X = y[ii], X[ii]
    b = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ b; return 1 - (res @ res) / ((y - y.mean()) @ (y - y.mean()))
MODELS = {"expression deciles": [], "+ transcript features": cols_feat, "+ cell-type specificity (PsychAD HBCC)": cols_ct, "features + cell type": cols_feat + cols_ct}
boots = [rng.integers(0, len(C), len(C)) for _ in range(500)]
rows = []
for target, label in (("pH_Z", "pH Z (control-only, v1 table)"), ("gtex_hardy_t", "GTEx Hardy t (slow vs fast death)")):
    base = r2(target, [])
    for name, cols in MODELS.items():
        v = r2(target, cols); bs = np.array([r2(target, cols, ii) for ii in boots]); bb = np.array([r2(target, [], ii) for ii in boots]) if name != "expression deciles" else np.zeros(len(boots))
        rows.append({"target": label, "model": name, "n_genes": len(C), "R2": v, "R2_ci95_lo": np.percentile(bs, 2.5), "R2_ci95_hi": np.percentile(bs, 97.5), "increment_over_deciles": v - base, "increment_ci95_lo": np.percentile(bs - bb, 2.5), "increment_ci95_hi": np.percentile(bs - bb, 97.5)})
    v = r2(target, cols_feat + cols_ct) - r2(target, cols_feat); bs = np.array([r2(target, cols_feat + cols_ct, ii) - r2(target, cols_feat, ii) for ii in boots])
    rows.append({"target": label, "model": "increment: cell type beyond features", "n_genes": len(C), "R2": v, "R2_ci95_lo": np.percentile(bs, 2.5), "R2_ci95_hi": np.percentile(bs, 97.5)})
R = pd.DataFrame(rows)
hdr = ["# NEW 2026-10-18: nested least squares of rank(target) on expression-decile indicators and rank-transformed predictors, for two targets on the SAME genes (complete for all predictors); targets: pH Z = control-only pH Z of confounds/results/combined_controls_pH.tsv (v1 table, as in the earlier Note 8 analysis), GTEx Hardy t = mean over BA9, BA24, Cortex of the slow-versus-fast-death t",
       "# script: agonal_celltype.py adapted from analyses/ph_sensitivity_gene_features_2026-10-09/feature_analysis.py (data loading and nested-model definitions unchanged); R2 intervals by 500 bootstrap resamples of genes; seed 20261009",
       f"# input: {os.path.relpath(phf, ROOT)} sha256 {sha(phf)}; hardy_t_by_region.tsv.gz sha256 {sha(hf)}; gene_features.tsv sha256 {sha(os.path.join(SRC, 'gene_features.tsv'))}; run {datetime.date.today()}"]
with open(os.path.join(HERE, "nested_models_two_targets.tsv"), "w") as f:
    for h in hdr: f.write(h + "\n")
    R.to_csv(f, sep="\t", index=False)
print(R.round(4).to_string(index=False))
