"""S5: compile bootstrap, halves and permutation outputs into summary tables (local, pandas only). Reads out/tables pulled from the remote environment."""
import pandas as pd, numpy as np, os, datetime, sys
T = sys.argv[1]; O = sys.argv[2]; os.makedirs(O, exist_ok=True)
VARS = ['v1repro', 'v2_bulk', 'v2_bulk_ofc', 'v2_bulk_ofc_lcm', 'v2_bulk_v1first']
LAB = {'v1repro': 'v1 (five datasets, v1 persons)', 'v2_bulk': 'v2 bulk (primary assignment)', 'v2_bulk_ofc': 'v2 bulk + NSW OFC snRNA', 'v2_bulk_ofc_lcm': 'v2 bulk + OFC + Arion LCM neurons', 'v2_bulk_v1first': 'v2 bulk (v1 datasets keep shared persons)'}
q = lambda x: (np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5))
def ci(point, reps, kind='raw'):
    """point +- 1.96 x bootstrap SE (the percentile interval of a donor bootstrap is shifted away from the point estimate for statistics that combine two noisy estimates, so the bootstrap is used for the SE only). kind 'r': on the Fisher z scale."""
    reps = np.asarray(reps, float); reps = reps[np.isfinite(reps)]
    if not np.isfinite(point) or len(reps) < 20: return (np.nan, np.nan)
    if kind == 'r':
        reps = np.clip(reps, -0.9999, 0.9999); se = np.std(np.arctanh(reps), ddof=1); z = np.arctanh(np.clip(point, -0.9999, 0.9999)); return (float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se)))
    se = np.std(reps, ddof=1); return (point - 1.96 * se, point + 1.96 * se)
hd = lambda t: f'# {t}; compiled {datetime.date.today()} by scripts/s5_compile.py from ph_axis_v2/out/tables (s3b_boot.py); intervals are point +- 1.96 x bootstrap SE (200 donor bootstraps; Fisher z scale for correlations), except halves: 2.5% to 97.5% of 200 random splits (the convention of the stored panel)\n'
def w(df, name, title):
    with open(f'{O}/{name}', 'w') as fh: fh.write(hd(title)); df.to_csv(fh, sep='\t', index=False, float_format='%.4g')
# ---------- reliability
rows = []
rng = np.random.default_rng(1)
B = {v: pd.read_csv(f'{T}/boot_{v}.tsv', sep='\t') for v in VARS if os.path.exists(f'{T}/boot_{v}.tsv')}
H = {v: pd.read_csv(f'{T}/halves_{v}.tsv', sep='\t') for v in VARS if os.path.exists(f'{T}/halves_{v}.tsv')}
HS = {v: pd.read_csv(f'{T}/halves_single_{v}.tsv', sep='\t') for v in VARS if os.path.exists(f'{T}/halves_single_{v}.tsv')}
pt = pd.read_csv(f'{T}/reliability_between_collections_point.tsv', sep='\t', comment='#')
ptmap = {'v1repro': 'v1 (five datasets, v1 persons)', 'v2_bulk': 'v2_bulk', 'v2_bulk_ofc': 'v2_bulk_ofc', 'v2_bulk_ofc_lcm': 'v2_bulk_ofc_lcm', 'v2_bulk_v1first': 'v2_bulk_v1first'}
for v in VARS:
    if v not in B: continue
    b = B[v][(B[v].model == 'ctrl') & (B[v].rep > 0)]; p0 = pt[pt.axis == ptmap[v]].iloc[0]
    r = dict(axis=LAB[v], n_collections=int(p0.n_collections), between_mean_pairwise_rho_point=p0.mean_pairwise_spearman, between_rho_lo=ci(p0.mean_pairwise_spearman, b.between_rbar, 'r')[0], between_rho_hi=ci(p0.mean_pairwise_spearman, b.between_rbar, 'r')[1],
             between_SB_point=p0.SB_all, between_SB_lo=ci(p0.SB_all, b.between_SB, 'r')[0], between_SB_hi=ci(p0.SB_all, b.between_SB, 'r')[1], between_SB4_point=p0.SB_4, between_SB4_lo=ci(p0.SB_4, b.between_SB4, 'r')[0], between_SB4_hi=ci(p0.SB_4, b.between_SB4, 'r')[1], between_SB3_point=p0.SB_3, between_SB3_lo=ci(p0.SB_3, b.between_SB3, 'r')[0], between_SB3_hi=ci(p0.SB_3, b.between_SB3, 'r')[1])
    if v in H: h = H[v]; r.update(halves_n_datasets=int(h.n_datasets.iloc[0]), halves_SB_mean=h.SB.mean(), halves_SB_lo=q(h.SB)[0], halves_SB_hi=q(h.SB)[1], halves_rho_mean=h.rho_halves.mean())
    rows.append(r)
R = pd.DataFrame(rows)
# difference v2 - v1 (independent bootstraps)
if 'v1repro' in B:
    for i, v in enumerate(VARS[1:]):
        if v not in B: continue
        b1 = B['v1repro'][(B['v1repro'].model == 'ctrl') & (B['v1repro'].rep > 0)].between_SB.values; b2 = B[v][(B[v].model == 'ctrl') & (B[v].rep > 0)].between_SB.values
        r1 = R[R.axis == LAB['v1repro']].iloc[0].between_SB_point; r2 = R[R.axis == LAB[v]].iloc[0].between_SB_point; se = np.sqrt(np.nanvar(b1, ddof=1) + np.nanvar(b2, ddof=1)); R.loc[R.axis == LAB[v], 'between_SB_diff_vs_v1'] = r2 - r1; R.loc[R.axis == LAB[v], 'between_SB_diff_lo'] = r2 - r1 - 1.96 * se; R.loc[R.axis == LAB[v], 'between_SB_diff_hi'] = r2 - r1 + 1.96 * se
w(R, 'reliability_summary.tsv', 'reliability of the control pH axis: between collections (Spearman-Brown of the mean pairwise Spearman among collection-level z) and between random donor halves (Spearman-Brown)')
rows = []
for v, h in HS.items():
    for d, g in h.groupby('dataset'): rows.append(dict(axis=LAB[v], dataset=d, n_ctrl=int(g.n_ctrl.iloc[0]), SB_mean=g.SB.mean(), SB_lo=q(g.SB)[0], SB_hi=q(g.SB)[1], raw_rho_mean=g.rho_halves.mean()))
w(pd.DataFrame(rows), 'reliability_single_dataset_halves.tsv', 'split-half reliability of one dataset\'s control pH t (datasets with at least 30 controls), 200 random splits, Spearman-Brown')
# ---------- headline
P = {v: pd.read_csv(f'{T}/perm_{v}.tsv', sep='\t') for v in VARS if os.path.exists(f'{T}/perm_{v}.tsv')}
hp = pd.read_csv(f'{T}/headline_point.tsv', sep='\t', comment='#').set_index('axis')
LM = pd.read_csv(f'{T}/loco_headline_summary.tsv', sep='\t', comment='#').set_index(['set', 'axis'])
LM.index = pd.MultiIndex.from_tuples([('5plus' if a == '5plus' else '209', b) for a, b in LM.index])
rows = []
stat_cols = [('5plus_r', 'r_5plus'), ('209_r', 'r_209'), ('5plus_auroc_absz', 'auroc_5plus'), ('209_auroc_absz', 'auroc_209'), ('agonal_rho', 'rho_agonal'), ('5plus_loco_median_r', 'loco_r_5plus'), ('209_loco_median_r', 'loco_r_209'), ('5plus_frac_sign_consistent', 'frac_consistent_5plus'), ('209_frac_sign_consistent', 'frac_consistent_209')]
ptkey = {('v1repro', 'ctrl'): 'v1 (stored)'}
for v in VARS:
    if v not in B: continue
    for model in (['ctrl', 'both', 'case'] if v != 'v1repro' else ['ctrl']):
        b = B[v][(B[v].model == model) & (B[v].rep > 0)]; b0 = B[v][(B[v].model == model) & (B[v].rep == 0)].iloc[0]
        key = 'v1 (stored)' if v == 'v1repro' else f'{v} / {model}'
        p = hp.loc[key]; r = dict(axis=LAB[v], model=model, n_genes_z=int(p.n_genes_z))
        for c, nm in stat_cols:
            kind = 'raw' if c.endswith('auroc_absz') or 'frac' in c else 'r'
            pnt = p[c] if c in p.index else (b0[c] if c in b0.index else np.nan)
            if c.startswith(('5plus_loco', '209_loco')):
                lk = ('v1, other four collections' if v == 'v1repro' else f'{v}, without the tested collection') if model == 'ctrl' else None
                pnt = LM.loc[(c.split('_')[0], lk), 'median'] if lk is not None and (c.split('_')[0], lk) in LM.index else b0[c]
            lo, hi = ci(pnt, b[c], kind); r[nm] = pnt; r[nm + '_lo'] = lo; r[nm + '_hi'] = hi; r[nm + '_boot_mean'] = np.nanmean(b[c])
        r['n5plus_with_both'] = int(p['5plus_n_with_both']); r['n209_with_both'] = int(p['209_n_with_both'])
        for s in ('5plus', '209'):
            r[f'{s}_obs'] = int(p[f'{s}_obs_v1method']); r[f'{s}_exp'] = p[f'{s}_exp_v1method']; r[f'{s}_obs_over_exp'] = p[f'{s}_obs_v1method'] / p[f'{s}_exp_v1method']
            if model == 'ctrl' and os.path.exists(f'{T}/vif_{v}.tsv'):
                vf = pd.read_csv(f'{T}/vif_{v}.tsv', sep='\t'); vf = vf[vf.set == s]
                for _, x in vf.iterrows(): r[f'{s}_z_corrected_{x.cutoff}'] = x.z_corrected; r[f'{s}_VIF_{x.cutoff}'] = x.VIF
                r[f'{s}_r_perm_z'] = vf.r_z_perm.iloc[0]
        rows.append(r)
Hh = pd.DataFrame(rows)
w(Hh, 'headline_summary.tsv', 'headline statistics, v1 beside v2. r: Pearson between the pH z and the schizophrenia z (Z_sz_base, five v1 datasets) over genes in the set; loco: median over the five v1 collections of r with the axis built without the tested collection; obs/exp: genes at q < 0.05 in the low-pH direction against the v1 expression-decile null (exact expectation); z_corrected_*: pH permuted among controls within each dataset (300 permutations), variance inflation of the count beyond the cutoff (q05_boundary: the observed q < 0.05 boundary; matched_rate: the cutoff at which the permutation exceedance rate equals the matched expectation, the stored convention), corrected z = (observed - expected) / (matched-null SD x sqrt(VIF)); r_perm_z: r against its permutation null; point estimates from limma fits, intervals from OLS donor bootstraps')
# differences v2 - v1 (independent bootstraps)
rows = []
for v in VARS[1:]:
    if v not in B: continue
    for model in ['ctrl', 'both', 'case']:
        b2 = B[v][(B[v].model == model) & (B[v].rep > 0)]; b1 = B['v1repro'][(B['v1repro'].model == 'ctrl') & (B['v1repro'].rep > 0)]
        for c, nm in stat_cols:
            p2 = Hh[(Hh.axis == LAB[v]) & (Hh.model == model)][nm].iloc[0]; p1 = Hh[(Hh.axis == LAB['v1repro']) & (Hh.model == 'ctrl')][nm].iloc[0]
            se = np.sqrt(np.nanvar(b2[c].values, ddof=1) + np.nanvar(b1[c].values, ddof=1)); rows.append(dict(axis=LAB[v], model=model, statistic=nm, v2_point=p2, v1_point=p1, diff_v2_minus_v1=p2 - p1, lo=p2 - p1 - 1.96 * se, hi=p2 - p1 + 1.96 * se, note='difference of points +- 1.96 x sqrt(sum of bootstrap variances)'))
w(pd.DataFrame(rows), 'headline_differences_v2_minus_v1.tsv', 'v2 minus v1 for the headline statistics (independent bootstraps of v1 and v2 persons, 5,000 random pairings)')
# ---------- control vs case
rows = []
for v in VARS:
    f = f'{T}/cc_boot_{v}.tsv'
    if not os.path.exists(f): continue
    c = pd.read_csv(f, sep='\t')
    for (u, gs), g in c.groupby(['unit', 'geneset']):
        b = g[g.rep > 0]; p0 = g[g.rep == 0].iloc[0]
        rows.append(dict(axis=LAB[v], unit=u, geneset=gs, n_genes=int(p0.n_genes), rho_t_point=p0.rho_t, rho_t_lo=ci(p0.rho_t, b.rho_t, 'r')[0], rho_t_hi=ci(p0.rho_t, b.rho_t, 'r')[1], rho_slope_point=p0.rho_slope, rho_slope_lo=ci(p0.rho_slope, b.rho_slope, 'r')[0], rho_slope_hi=ci(p0.rho_slope, b.rho_slope, 'r')[1], ratio_raw_point=p0.ratio_raw, ratio_raw_lo=q(b.ratio_raw)[0], ratio_raw_hi=q(b.ratio_raw)[1],
                         ratio_corrected_point=(p0.ratio_corrected if p0.reliability_ctrl >= 0.2 else np.nan), ratio_corrected_lo=(q(b.ratio_corrected)[0] if p0.reliability_ctrl >= 0.2 else np.nan), ratio_corrected_hi=(q(b.ratio_corrected)[1] if p0.reliability_ctrl >= 0.2 else np.nan), reliability_ctrl_slopes=p0.reliability_ctrl))
CC = pd.DataFrame(rows); w(CC, 'control_vs_case_agreement.tsv', 'agreement of the pH association in cases with that in controls: Spearman of t, Spearman of standardized slopes (beta per SD of expression of the dataset), ratio of case to control slope (regression through the origin; corrected = divided by the share of control-slope variance that is not sampling noise); pooled = inverse-variance mean over datasets; intervals: correlations point +- 1.96 x bootstrap SE (Fisher scale), ratios 2.5% to 97.5% of the bootstrap replicates (donor bootstraps within group); point = OLS on original persons; corrected ratio shown only when the reliability of the control slopes is at least 0.2')
print('compiled')
