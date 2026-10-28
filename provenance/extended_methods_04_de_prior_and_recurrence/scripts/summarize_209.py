#!/usr/bin/env python3
"""Summarize the remote draws (run_null.py draw) into the small TSVs kept in the repo.
python3 summarize.py  ->  recurrence_null_summary.tsv, recurrence_null_draws_primary.tsv, sz_observed.tsv, population_contrasts.tsv"""
import os, json, hashlib, datetime, pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()[:12]
d = pd.read_csv(f'{HERE}/_draws_raw.tsv', sep='\t'); sz = json.load(open(f'{HERE}/_sz_observed.json')); cc = json.load(open(f'{HERE}/_common_counts.json'))
PROV = (f"# script: analyses/redo_on_209_2026-10-11/4_recurrence_null/remote_run/summarize_209.py (draws from run_null.py on remote, 2000 draws per row, seed 20261010, run {datetime.date.today()}); "
        f"draws file _draws_raw.tsv sha256 {sha(HERE+'/_draws_raw.tsv')}; SZ observed = analyses/redo_on_209_2026-10-11/4_recurrence_null/sz_gene_counts.tsv sha256 {sha(HERE+'/sz_gene_counts.tsv')} "
        f"(62 'published, brain' study-tier lists, export_sz_inputs.py); Gemma cache = remote analyses/gemma_176_scan_2026-10-10/cache/de_cache.h5; public experiments only\n")
d['variant'] = d['variant'].fillna('random')
rows = []
for (pop, rule, var), g in d.groupby(['pop', 'rule', 'variant']):
    r = dict(population=pop, rule=rule, n_draws=len(g), median_list_size=g.median_list_size.median(), mean_list_size=g.mean_list_size.mean())
    for c, szc in (('n5', 'n5'), ('n10', 'n10'), ('n5_common', None), ('n10_common', None), ('n5_dir', 'n5_dir'), ('n10_dir', 'n10_dir'), ('n_ge2', 'n_ge2'), ('n_ge3', 'n_ge3'),
                   ('major_frac_5', 'major_frac_5'), ('frac_in_sz1478', None), ('overlap_sz176', None), ('auc_prior_5', 'auc_prior_5'), ('mean_prior_pct_5', 'mean_prior_pct_5'), ('frac_top_decile_prior_5', 'frac_top_decile_prior_5')):
        q = g[c].quantile([.025, .5, .975]).values
        r[f'{c}_lo'], r[f'{c}_median'], r[f'{c}_hi'] = q
    r['sz_n5'] = sz['n5']; r['sz_n10'] = sz['n10']
    r['p_draws_ge_sz_n5'] = (g.n5 >= sz['n5']).mean(); r['p_draws_ge_sz_n10'] = (g.n10 >= sz['n10']).mean()
    r['n_draws_ge_sz_n5'] = int((g.n5 >= sz['n5']).sum()); r['n_draws_ge_sz_n10'] = int((g.n10 >= sz['n10']).sum())
    r['exp_overlap_sz1478_median'] = (g.exp_overlap_sz1478).median()
    r['exp_overlap_sz176_median'] = float((g.n5 * cc['sz176_common'] / cc['n_common']).median())
    rows.append(r)
S = pd.DataFrame(rows)
order = {'A': 0, 'A2': 1, 'B2': 2, 'B': 3, 'C': 4, 'N': 5}
S['o'] = S.rule.map(order); S = S.sort_values(['population', 'o']).drop(columns='o')
with open(f'{HERE}/recurrence_null_summary.tsv', 'w') as f:
    f.write(PROV)
    f.write('# rules: A = top-L by p (ties by |t|), L = the 62 SZ list sizes (list_size_in_universe, sz_lists_62.tsv) permuted over 62 random contrasts; A2 = L = SZ rate (list size / n_k) x genes tested in the contrast; B2 = FDR < 0.05 hits (contrasts with >= 1 hit) capped at the size drawn for A; B = all FDR < 0.05 hits (contrasts with >= 1 hit); N = all nominal p < 0.05; C = 36 slots as B (SZ has 36 FDR lists), 26 as N\n')
    f.write('# columns *_lo / *_median / *_hi = 2.5%, 50%, 97.5% over draws; n5, n10 = genes in >= 5 / >= 10 of 62 lists (sign-blind); *_dir = genes with >= 5 / >= 10 lists in the same direction (max of up, down; SZ: the lists\' recorded direction); frac_in_sz1478 = share of the draw\'s 5+ genes that are among the SZ 1,478 (by symbol; 27,146 symbols common to both universes); overlap_sz176 = draw 5+ genes among the SZ recurrent set (209); frac_in_sz1478 / overlap_sz1478 / exp_overlap_sz1478_median keep the old column names and now refer to the SZ 5+ set (1,641)\n')
    S.to_csv(f, sep='\t', index=False, float_format='%.5g')
szo = pd.DataFrame([dict(statistic=k, value=v) for k, v in sz.items()])
with open(f'{HERE}/sz_observed.tsv', 'w') as f:
    f.write(PROV); f.write('# the 62 schizophrenia lists on the same definitions as the draws (computed in run_null.py::sz_observed from sz_gene_counts.tsv); auc_prior_5 here uses n_de_total (missing -> median) over the 33,188 SZ-universe symbols\n'); szo.to_csv(f, sep='\t', index=False, float_format='%.5g')
prim = d[(d.rule.isin(['A', 'A2'])) & (d['pop'].isin(['P1_disease', 'P2_all_postmortem_brain']))]
cols = ['pop', 'rule', 'draw', 'n5', 'n10', 'n5_dir', 'n10_dir', 'n_ge2', 'n_ge3', 'median_list_size', 'frac_in_sz1478', 'overlap_sz176', 'auc_prior_5', 'mean_prior_pct_5', 'major_frac_5']
with open(f'{HERE}/recurrence_null_draws_primary.tsv', 'w') as f:
    f.write(PROV); f.write('# per-draw values for rules A and A2 in populations P1_disease and P2_all_postmortem_brain; all other rows stay in remote _draws_raw.tsv\n'); prim[cols].to_csv(f, sep='\t', index=False, float_format='%.5g')
m = pd.read_csv(f'{HERE}/_pool_meta.tsv', sep='\t')
keep = ['gemma_id', 'short_name', 'factor', 'contrast_label', 'ctype', 'disorder', 'true_disease', 'case_like', 'n_samples', 'platform', 'sz_status', 'region', 'n_tested', 'n_fdr', 'n_nom']
with open(f'{HERE}/population_contrasts.tsv', 'w') as f:
    f.write(PROV); f.write('# the 277 public postmortem-brain contrasts (one per experiment; analyses/gemma_176_scan_nonhits_2026-10-10/postmortem_brain_contrasts_classified.tsv) with the number of genes tested and genes at BH q < 0.05 / p < 0.05 in the cache; pools: P1 = true_disease & not SZ-related, P2 = not SZ-related (sz_status neither and disorder not schizophrenia)\n'); m[keep].to_csv(f, sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(S[['population', 'rule', 'median_list_size', 'n5_lo', 'n5_median', 'n5_hi', 'n10_lo', 'n10_median', 'n10_hi', 'n_draws_ge_sz_n5', 'n_draws_ge_sz_n10']].round(1).to_string())
print(S[['population', 'rule', 'n5_dir_median', 'n10_dir_median', 'n_ge2_median', 'n_ge3_median', 'major_frac_5_median', 'frac_in_sz1478_median', 'exp_overlap_sz1478_median', 'overlap_sz176_median', 'exp_overlap_sz176_median', 'auc_prior_5_median', 'mean_prior_pct_5_median']].round(3).to_string())
print(sz); print(json.dumps(cc.get('pools'), indent=0))
