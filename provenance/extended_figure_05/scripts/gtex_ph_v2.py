#!/usr/bin/env python3
"""G3 sweep 2026-10-18: rho_pH_Z / n_pH_Z of the GTEx death-speed tables on pH axis v2, from the stored per-tissue t matrices (GTEx v8 OLS not rerun).
Same computation as analyses/gtex_death_speed_all_tissues_2026-10-10/gtex_death_speed_all_tissues.py: gene t (first of duplicate symbols) of one tissue and contrast,
genes present in the axis, Spearman rho (scipy.stats.spearmanr), more than 10 genes. The v1 axis (column 4 = Z of the pH rows of combined_controls_pH.tsv) is recomputed
alongside to compare with the stored values (the stored matrices are rounded to 4 decimals).
usage: gtex_ph_v2.py OUT.tsv"""
import sys, pandas as pd, numpy as np
from scipy.stats import spearmanr
D = 'sztr-rebuild-infra/repo/analyses/gtex_death_speed_all_tissues_2026-10-10'; H3 = 'gtexh3/out'
AN = 'axes-handoff/analyses/scan_on_209_2026-10-11/branch/analyses'
V2F = 'nserc_pilots_2026-10-14/ph_axis_v2/out/tables/axis_v2_controls.tsv'
V1F = f'{AN}/confounds/results/combined_controls_pH.tsv'
v1 = {}
for l in open(V1F):
    if l.startswith('#') or l.startswith('gene\t'): continue
    f = l.rstrip('\n').split('\t')
    if f[1] == 'pH':
        try: v1[f[0]] = float(f[3])
        except ValueError: pass
v2 = pd.read_csv(V2F, sep='\t', comment='#').drop_duplicates('gene').set_index('gene').z_stouffer.dropna().to_dict()
MAT = [('slow_4_vs_fast_2', f'{D}/slow_vs_fast_t_by_tissue.tsv.gz'), ('ventilator_0_vs_fast_2', f'{D}/ventilator_vs_fast_t_by_tissue.tsv.gz'),
       ('hardy3_vs_hardy2', f'{H3}/hardy3_vs_hardy2_t_by_tissue.tsv.gz'), ('hardy4_vs_hardy3', f'{H3}/hardy4_vs_hardy3_t_by_tissue.tsv.gz')]
rows = []
for con, fn in MAT:
    M = pd.read_csv(fn, sep='\t', comment='#', index_col=0, compression='gzip')
    for tissue in M.columns:
        ser = M[tissue].dropna(); r = dict(tissue=tissue, contrast=con)
        for nm, ad in (('v1', v1), ('v2', v2)):
            k = [x for x in ser.index if x in ad]
            rho = spearmanr(ser[k].values, [ad[x] for x in k])[0] if len(k) > 10 else np.nan
            r[f'rho_pH_Z_{nm}'] = round(float(rho), 4); r[f'n_pH_Z_{nm}'] = len(k)
        rows.append(r)
pd.DataFrame(rows).to_csv(sys.argv[1], sep='\t', index=False)
print(len(rows), 'rows')
