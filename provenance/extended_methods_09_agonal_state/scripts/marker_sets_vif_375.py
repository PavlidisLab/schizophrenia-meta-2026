#!/usr/bin/env python3
"""STAT-012 addendum 2, item 3a: Fig. 5g. Marker sets (100 most specific genes of each PsychAD class, class mean CPM >= 1, with a control pH Z) against expression-decile-matched
random sets, with a 95% bootstrap interval over genes (analyses/ph_sensitivity_gene_features_2026-10-09/feature_analysis.py lines 86-98). Selection and statistic re-implemented here
from those lines and checked against results_laptop/marker_sets.tsv.
Correction. The control pH Z of a set gene is a Stouffer combination over five collections; its sample-permutation noise (pH permuted among controls within each collection,
ph_count_perm.R, 2,000 permutations, all genes) gives the variance of the set mean, Var_perm(mean Z over the set). The independent-genes value is mean(Var_perm(Z_g)) / m.
VIF_set = Var_perm(set mean) / that value. The bootstrap interval over genes is widened by sqrt(VIF_set) about the point estimate; the corrected z is effect / (bootstrap SE x sqrt(VIF_set)).
Self-contained check: observed set mean against the permutation distribution of the set mean (null centred on zero)."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common import sha, write_tsv, today
REPO = "sztr-rebuild-infra/repo"
PD = "sztr-rebuild-infra/runs/refreeze_375_2026-10-15/vif_ge10/perm_data"   # copy: permutation arrays of X4 (rebuilt control pH Z)
PH = os.path.join(REPO, "analyses/confounds/results/combined_controls_pH.tsv")
MS = os.path.join(REPO, "analyses/ph_sensitivity_gene_features_2026-10-09/results_laptop/marker_sets.tsv")
SP = {c: os.path.join(REPO, f"analyses/interneuron_ph_2026-10-01/results/specificity_{c}.tsv.gz") for c in ("HBCC", "MSSM")}
rng = np.random.default_rng(20261012); NNULL = 10000
ph = pd.read_csv(PH, sep="\t", comment="#"); ph = ph[ph.factor == "pH"].set_index("gene")
D = pd.DataFrame({"pH_Z": ph.Z, "expr_pct": ph.expr_percentile}).dropna(subset=["pH_Z", "expr_pct"]); D["decile"] = np.minimum((D.expr_pct * 10).astype(int), 9)
genes = [l.strip() for l in open(os.path.join(PD, "perm_Z_allgenes_genes.txt"))]
NP = 2000; G = len(genes)
Zp = np.fromfile(os.path.join(PD, "perm_Z_allgenes_int16x1000.bin"), dtype=np.int16).reshape(NP, G).astype(np.float64) / 1000.0   # permutations x genes
gi = {g: i for i, g in enumerate(genes)}
obsall = pd.read_csv(os.path.join(PD, "observed_allgenes.tsv"), sep="\t").set_index("gene")
assert np.allclose(obsall.Z.reindex(D.index).values, D.pH_Z.values, atol=1e-6), "pH Z differs from combined_controls_pH.tsv"
ms = pd.read_csv(MS, sep="\t", comment="#")
spec = {c: pd.read_csv(p, sep="\t", comment="#", index_col=0) for c, p in SP.items()}
rows = []
for coh, Sx in spec.items():
    for cls in ("SST", "PVALB", "VIP", "LAMP5", "IN_other", "EN", "Astro", "Oligo", "OPC", "Micro", "Vascular"):
        ok = Sx[f"cpm_{cls}"] >= 1; top = Sx.loc[ok, f"spec_{cls}"].sort_values(ascending=False).index[:100]
        T = D.loc[D.index.intersection(top)]; rest = D.drop(index=T.index); n = len(T)
        if n < 20: continue
        old = ms[(ms.cohort == coh) & (ms["class"] == cls)].iloc[0]
        by = {k: rest.pH_Z.values[rest.decile.values == k] for k in range(10)}; need = T.decile.value_counts().to_dict()
        sums = np.zeros(NNULL)
        for k, c in need.items(): sums += np.array([rng.choice(by[k], c, replace=False).sum() for _ in range(NNULL)])
        nullm = (sums / n).mean(); obs = T.pH_Z.mean(); eff = obs - nullm
        # bootstrap over genes (same procedure)
        boot = np.array([T.pH_Z.values[rng.integers(0, n, n)].mean() for _ in range(2000)]) - nullm
        lo, hi = np.percentile(boot, [2.5, 97.5]); se_b = boot.std(ddof=1)
        ix = [gi[g] for g in T.index]; sm = Zp[:, ix].mean(1)
        v_perm = sm.var(ddof=1); v_ind = Zp[:, ix].var(0, ddof=1).mean() / n; vif = v_perm / v_ind
        rho_bar = (vif - 1) / (n - 1); f = np.sqrt(vif)
        zc = eff / (se_b * f)
        rows.append(dict(cohort=coh, cls=cls, n_genes=n, mean_pH_Z=obs, matched_null_mean=nullm, effect=eff, effect_ci_lo=lo, effect_ci_hi=hi,
                         stored_effect=old.effect, stored_ci_lo=old.effect_ci95_lo, stored_ci_hi=old.effect_ci95_hi, boot_se=se_b,
                         perm_sd_of_set_mean=np.sqrt(v_perm), indep_sd_of_set_mean=np.sqrt(v_ind), VIF=vif, mean_pairwise_corr=rho_bar,
                         corrected_ci_lo=eff - (eff - lo) * f, corrected_ci_hi=eff + (hi - eff) * f, corrected_z=zc, corrected_p_two_sided=2 * norm.sf(abs(zc)),
                         ci_excludes_zero_corrected=bool((eff - (eff - lo) * f > 0) or (eff + (hi - eff) * f < 0)),
                         selfcontained_z=obs / np.sqrt(v_perm), selfcontained_p_two_sided=(1 + (np.abs(sm) >= abs(obs)).sum()) / (NP + 1)))
        print(coh, cls, n, round(obs, 3), "stored", round(old.mean_pH_Z, 3), "eff", round(eff, 3), "VIF", round(vif, 1), "corr CI", round(rows[-1]["corrected_ci_lo"], 2), round(rows[-1]["corrected_ci_hi"], 2), flush=True)
R = pd.DataFrame(rows).rename(columns={"cls": "class"})
hdr = ["claim: Fig. 5g, marker-set mean control pH z minus expression-decile-matched random sets, 95% bootstrap interval over genes; interval corrected for correlation among set genes",
       "script: analyses/stat012_vif_2026-10-12/addendum2_2026-10-12/marker_sets_vif.py (perm_Z from ph_count_perm.R with SAVE_ALL=1, remote)",
       f"source: {os.path.relpath(PH, REPO)} sha256 {sha(PH)[:12]}; {os.path.relpath(MS, REPO)} sha256 {sha(MS)[:12]}; specificity_HBCC.tsv.gz {sha(SP['HBCC'])[:12]}; specificity_MSSM.tsv.gz {sha(SP['MSSM'])[:12]}; perm_Z_allgenes_int16x1000.bin {sha(os.path.join(PD, 'perm_Z_allgenes_int16x1000.bin'))[:12]}",
       f"date: {today()}; 2,000 sample permutations; matched null redrawn with {NNULL} draws, seed 20261012 (stored values are from the 2026-10-09 run, seed 20261009)"]
write_tsv(os.path.join(HERE, "fig5g_marker_sets_vif.tsv"), hdr, R)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(R[["cohort", "class", "n_genes", "effect", "stored_effect", "effect_ci_lo", "effect_ci_hi", "VIF", "corrected_ci_lo", "corrected_ci_hi", "corrected_z", "selfcontained_z", "selfcontained_p_two_sided"]].round(3).to_string())
