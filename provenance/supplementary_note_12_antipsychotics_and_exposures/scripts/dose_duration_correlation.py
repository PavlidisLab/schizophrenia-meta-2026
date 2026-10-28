"""F_notes_10-12: correlation of lifetime antipsychotic dose with duration of illness among the Stanley GSE35977 cases (aggregates only).
Input: .cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv (SMRI workbook match, laptop-only; donor rows are never written).
Dose = log10(apd_w clipped at 1) as in c15_final_core80.py; duration = dur_w. Groups as in c15: SZ only (49 donors with dose and duration) and SZ+BP.
Also the variance inflation factor of dose in the c15 SZ-only model with duration (score ~ dose + duration + age + sex + pH + PMI), computed on the covariates only.
Run: python3 -I dose_duration_correlation.py > ../tables/dose_duration_correlation.tsv   date 2026-10-19; no random numbers used (seed: none)."""
import os, hashlib, numpy as np, pandas as pd
from scipy import stats
F = os.path.expanduser('.cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv')
sha = hashlib.sha256(open(F, 'rb').read()).hexdigest()
M = pd.read_csv(F, sep='\t', comment='#'); M = M[M.series == 'GSE35977'].copy()
M['logapd'] = np.log10(M.apd_w.clip(lower=1))
print('# script: F/scripts/dose_duration_correlation.py (run 2026-10-19)')
print('# input: .cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv sha256 %s (laptop-only; aggregates only written)' % sha[:12])
print('group\tn\tpearson_r\tpearson_lo95\tpearson_hi95\tspearman_rho\tvif_dose_given_duration_age_sex_pH_PMI\tvif_duration_given_dose_age_sex_pH_PMI')
for grp, dxs in [('SZ only', ['schizophrenia']), ('SZ+BP', ['schizophrenia', 'BP'])]:
    d = M[M.dx_w.isin(dxs)].dropna(subset=['apd_w', 'dur_w', 'pH_w', 'PMI_w', 'age_w', 'sex_w']).copy()
    r, p = stats.pearsonr(d.logapd, d.dur_w); rho, _ = stats.spearmanr(d.logapd, d.dur_w)
    z = np.arctanh(r); se = 1 / np.sqrt(len(d) - 3); lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    def vif(y, others):
        X = np.column_stack([np.ones(len(d))] + others); b = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ b
        r2 = 1 - (res @ res) / (((y - y.mean()) ** 2).sum()); return 1 / (1 - r2)
    cov = [d.age_w.values, (d.sex_w == 'M').astype(float).values, d.pH_w.values, d.PMI_w.values] + ([(d.dx_w == 'BP').astype(float).values] if grp == 'SZ+BP' else [])
    print('%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f' % (grp, len(d), r, lo, hi, rho, vif(d.logapd.values, [d.dur_w.values] + cov), vif(d.dur_w.values, [d.logapd.values] + cov)))
