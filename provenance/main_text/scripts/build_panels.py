#!/usr/bin/env python3
"""Panels a, b, c of Extended Fig. 11 for one gene set from score375.py shards, and the per-contrast table. Writes only into OUT.
Panels a, c: as calibrated_panels/build_calibrated_panels.py (hit = z of the set > 97.5th percentile of z of its own 200 decile-matched sets, and AUROC > 0.5; groups as that script).
Panel b: as redo_on_209_2026-10-11/1_gemma_scan/scan_intervals_209.py (median |rho| of t with the pH and GTEx axes, 2000 bootstrap resamples of experiments, seed 20261010).
usage: build_panels.py SHARDDIR OUT LABEL:group[,LABEL:group...]   e.g. 375:g375,ge10:ge10,ge5:ge5,old209:old209"""
import sys, os, glob, hashlib, datetime, numpy as np, pandas as pd
SD, OUT, SPEC = sys.argv[1:4]; IN = 'gemma_scans_recurrent_sets/inputs'
P_ALL = f'{IN}/all_contrasts_scored_public.tsv'; P_SCAN = f'{IN}/scan176_contrasts.tsv'; P_COI = f'{IN}/coi_scored.tsv'; P_SETS = f'{IN}/../out/sets_scoring.tsv'
P_CLS = 'sztr-rebuild-infra/repo/analyses/gemma_176_scan_figure_2026-10-10/remote_followup/flagged_candidates_classified.tsv'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12]
rd = lambda p: pd.read_csv(p, sep='\t', comment='#', dtype=str, low_memory=False)
key = lambda d: d.gemma_id + '__' + d.rsid + '__' + d.fvid
Zs = [np.load(f, allow_pickle=True) for f in sorted(glob.glob(f'{SD}/shards/shard_*.npz'))]; cat = lambda k: np.concatenate([z[k] for z in Zs])
gid, rsid, fvid = cat('gid'), cat('rsid'), cat('fvid'); z = cat('z').astype(float); au = cat('auroc').astype(float); npres = cat('npres')
rpH, rgt, ng = cat('rho_pH_all').astype(float), cat('rho_gtex_all').astype(float), cat('n_genes'); sid = np.array(Zs[0]['set_id'])
kz = pd.Series(range(len(gid)), index=[f'{a}__{b}__{c}' for a, b, c in zip(gid, rsid, fvid)])
def wilson(k, n, zc=1.959964):
    if n == 0: return np.nan, np.nan
    p = k / n; d = 1 + zc * zc / n; c = (p + zc * zc / (2 * n)) / d; h = zc * np.sqrt(p * (1 - p) / n + zc * zc / (4 * n * n)) / d; return max(0, c - h), min(1, c + h)
S0 = rd(P_SCAN); S = S0[(S0.is_public == 'True') & (S0.picked == 'True') & (S0.is_single_cell == 'False') & S0.n_genes.notna()].copy()
for c in ['in_sz_corpus', 'held_out', 'flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']: S[c] = S[c] == 'True'
S['any_flag'] = S[['flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']].any(axis=1)
S['row'] = kz.reindex(key(S)).values; n_unscored = int(S.row.isna().sum()); S = S[S.row.notna()].copy(); S['row'] = S.row.astype(int)
GA = [('Postmortem brain', 'brain_postmortem'), ('Other brain', 'brain_other'), ('Other tissues', 'other'), ('Cell lines / in vitro', 'in_vitro_or_cell_line'), ('Blood', 'blood')]
PRIMARY = ['hypoxia', 'ischemia', 'hypoxia_mimetic', 'altitude', 'pH_acidosis', 'PMI_agonal']
Cc = pd.read_csv(P_CLS, sep='\t', dtype={'gemma_id': str, 'fvid': str, 'rsid': str}).fillna(''); Co = rd(P_COI); Co = Co[Co.status == 'scored'].copy()
Co = Co.merge(Cc[['gemma_id', 'tissue_class', 'is_single_cell']].drop_duplicates('gemma_id'), on='gemma_id', how='left'); Co = Co[~Co.is_single_cell.astype(bool)].sort_values(['gemma_id', 'rsid', 'fvid'])
coi_h = Co[Co.cls.isin(PRIMARY)].drop_duplicates('gemma_id'); coi_d = Co[Co.cls == 'disease_ischemia_state'].drop_duplicates('gemma_id')
rows_of = lambda d: kz.reindex(key(d)).dropna().astype(int).values
A = rd(P_ALL); A['key'] = key(A)
meta = A[['key', 'gemma_id', 'rsid', 'fvid', 'is_pick', 'tissue_class', 'cond_class', 'in_sz_corpus', 'held_out', 'is_microarray', 'n_samples']].copy()
S0['key'] = key(S0); meta = meta.merge(S0[['key', 'picked', 'is_public', 'is_single_cell', 'case_like']].drop_duplicates('key'), on='key', how='outer')
sums = f'scored shards {SD}: {len(gid)} contrasts; sets_scoring.tsv {sha(P_SETS)}; scan176_contrasts.tsv {sha(P_SCAN)}; all_contrasts_scored_public.tsv {sha(P_ALL)}; coi_scored.tsv {sha(P_COI)}; flagged_candidates_classified.tsv {sha(P_CLS)}'
today = datetime.date.today().isoformat()
def wr(df, fn, stamp, extra=()):
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn, 'w') as f:
        for s in list(stamp) + list(extra): f.write('# ' + s + '\n')
        df.to_csv(f, sep='\t', index=False, float_format='%.5g')
percol = {}
for spec in SPEC.split(','):
    lab, grp = spec.split(':'); j = int(np.where(sid == grp)[0][0]); cols = np.array([i for i, s in enumerate(sid) if str(s).startswith(f'{grp}__indep_')])
    Zn = z[:, cols]; z0 = z[:, j]; ok = np.isfinite(z0) & np.isfinite(Zn).all(1)
    q = np.full(len(z), np.nan); q[ok] = np.percentile(Zn[ok], 97.5, axis=1); hit = ok & (z0 > q) & (au[:, j] > 0.5)
    loo = np.zeros(Zn.shape, bool)
    for c in range(Zn.shape[1]):
        o = np.delete(Zn[ok], c, axis=1); loo[np.where(ok)[0], c] = Zn[ok, c] > np.percentile(o, 97.5, axis=1)
    loo &= (au[:, cols] > 0.5)
    percol[lab] = dict(z=z0, auroc=au[:, j], npres=npres[:, j], q975=q, hit=hit.astype(float) * np.where(ok, 1, np.nan))
    def stats(ix):
        ix = np.asarray([i for i in ix if ok[i]]); n = len(ix); k = int(hit[ix].sum()) if n else 0; lo, hi = wilson(k, n); sr = loo[ix].mean(0) if n else np.full(len(cols), np.nan)
        return dict(k=k, n=n, share=k / n if n else np.nan, wilson95_lo=lo, wilson95_hi=hi, comparison_set_rate=float(np.mean(sr)), comparison_set_rate_lo=float(np.percentile(sr, 2.5)), comparison_set_rate_hi=float(np.percentile(sr, 97.5)))
    out_a = [dict(panel='a', group=g, **stats(S.row[S.tissue_class == t].values)) for g, t in GA]
    GC = [('SZ corpus (picked)', S.row[S.in_sz_corpus].values), ('SZ-annotated held out', S.row[(~S.in_sz_corpus) & S.held_out].values), ('All other experiments', S.row[(~S.in_sz_corpus) & ~S.held_out].values),
          ('Design flagged by regex (picked)', S.row[S.any_flag].values), ('Hypoxia / ischemia / pH contrast scored', rows_of(coi_h)), ('Ischemia as disease state (stroke, MI)', rows_of(coi_d))]
    out_c = [dict(panel='c', group=g, **stats(r)) for g, r in GC]
    cols_out = ['panel', 'group', 'k', 'n', 'share', 'wilson95_lo', 'wilson95_hi', 'comparison_set_rate', 'comparison_set_rate_lo', 'comparison_set_rate_hi']
    nm = {'375': 'the 375 genes (study tier q < 0.1)', 'ge10': 'the 107 genes in 10 or more lists', 'ge5': 'the 2,370 genes in 5 or more lists (sensitivity)', 'old209': 'the frozen 209 (validation)'}[lab]
    st = [f'claim: Extended Fig. 11 panels a and c under the calibrated per-contrast rule ({nm}; hit = z of the set above the 97.5th percentile of the contrast\'s own null of 200 decile-matched comparison sets, and AUROC above 0.5)',
          f'script: gemma_scans_recurrent_sets/scripts/build_panels.py (build_sets_375.py, score375.py; arithmetic of calibrated_panels/score_grouped.py and build_calibrated_panels.py)', 'source: ' + sums, f'date: {today}']
    cmt = ['comparison_set_rate = mean over the 200 comparison sets of the set\'s own hit rate in the group (each set against the other 199 sets of the contrast); _lo, _hi = 2.5th and 97.5th percentile of the per-set rates', 'k, n = contrasts that are hits, contrasts counted (contrasts with no z for the set or a comparison set are not counted)']
    wr(pd.DataFrame(out_a)[cols_out], f'{OUT}/{lab}/gemma_209_scan_calibrated/panel_a.tsv', st, cmt); wr(pd.DataFrame(out_c)[cols_out], f'{OUT}/{lab}/gemma_209_scan_calibrated/panel_c.tsv', st, cmt)
    # ---- panel b
    R = S0[S0.is_public.notna()].copy(); R = R[R.is_public == 'True'].copy()
    for c in ['picked', 'is_single_cell']: R[c] = R[c] == 'True'
    P = R[R.picked & ~R.is_single_cell].copy(); assert P.gemma_id.is_unique
    rw = kz.reindex(key(P)).values; has = ~np.isnan(rw); rwi = np.where(has, rw, 0).astype(int)
    P['auroc_abs_all'] = np.where(has, au[rwi, j], np.nan); P['z_abs_all_prior'] = np.where(has, z0[rwi], np.nan)
    P['abs_rho_pH'] = np.where(has, np.abs(rpH[rwi]), np.nan); P['abs_rho_gtex'] = np.where(has, np.abs(rgt[rwi]), np.nan)
    for c in ['in_sz_corpus', 'held_out', 'case_like', 'flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia', 'is_microarray']: P[c] = P[c] == 'True'
    P['any_flag'] = P[['flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']].any(axis=1)
    groups = [('all picked', P.index == P.index)] + [(f'tissue_class={t}', (P.tissue_class == t).values) for t in ['brain_postmortem', 'brain_other', 'other', 'in_vitro_or_cell_line', 'blood']]
    groups += [('postmortem brain, case-like', ((P.tissue_class == 'brain_postmortem') & P.case_like).values), ('postmortem brain, not case-like', ((P.tissue_class == 'brain_postmortem') & ~P.case_like).values), ('other brain, case-like', ((P.tissue_class == 'brain_other') & P.case_like).values),
               ('SZ corpus (in_sz_corpus)', P.in_sz_corpus.values), ('SZ-annotated held out, not in corpus', ((~P.in_sz_corpus) & P.held_out).values), ('neither corpus nor held out', ((~P.in_sz_corpus) & ~P.held_out).values),
               ('microarray', P.is_microarray.values), ('RNA-seq', (~P.is_microarray).values), ('any design flag (regex)', P.any_flag.values), ('no design flag', (~P.any_flag).values)] + [(f, P[f].values) for f in ['flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']]
    rng = np.random.default_rng(20261010)
    def boot(x, B=2000):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        if len(x) < 2: return np.nan, np.nan, np.nan, len(x)
        idx = rng.integers(0, len(x), (B, len(x))); med = np.median(x[idx], axis=1); return np.median(x), *np.percentile(med, [2.5, 97.5]), len(x)
    M = []
    for name, m in groups:
        d = P[m]
        for l2, col in (('AUROC', 'auroc_abs_all'), ('pH', 'abs_rho_pH'), ('GTEx', 'abs_rho_gtex'), ('z', 'z_abs_all_prior')):
            med, lo, hi, n = boot(d[col]); M.append(dict(group=name, statistic=l2, median=med, boot95_lo=lo, boot95_hi=hi, n=n))
    M = pd.DataFrame(M); rb = []
    for g, t in GA:
        for l2, ax in (('pH', 'pH axis'), ('GTEx', 'GTEx axis')):
            r = M[(M.group == f'tissue_class={t}') & (M.statistic == l2)].iloc[0]; rb.append(dict(panel='b', group=g, axis=ax, median_abs_rho=r['median'], boot95_lo=r.boot95_lo, boot95_hi=r.boot95_hi, n_experiments=int(r.n)))
    stb = [f'figure: gemma_176_scan_figure (panel b) for {nm}; script: gemma_scans_recurrent_sets/scripts/build_panels.py (arithmetic of redo_on_209_2026-10-11/1_gemma_scan/scan_intervals_209.py and build_figure_209.py)',
           'derived from the per-contrast rho of score375.py (cached result sets), picked non-single-cell public contrasts, one per experiment; rho (all shared genes) does not depend on the gene set; bootstrap = 2,000 resamples of experiments, seed 20261010', 'source: ' + sums, f'date: {today}']
    wr(pd.DataFrame(rb), f'{OUT}/{lab}/gemma_176_scan_figure/panel_b.tsv', stb)
    pd.set_option('display.width', 250); print(lab, grp, 'population rows not scored:', n_unscored); print(pd.DataFrame(out_a)[cols_out].round(4).to_string()); print(pd.DataFrame(out_c)[cols_out].round(4).to_string()); print(pd.DataFrame(rb).round(4).to_string())
T = meta[meta.key.isin(kz.index)].copy(); r_ = kz.reindex(T.key).values.astype(int)
T.insert(0, 'row_check', 1); T = T.drop(columns=['row_check', 'key'])
T['n_genes_scored'] = ng[r_]; T['rho_pH_all'] = rpH[r_]; T['rho_gtex_all'] = rgt[r_]
for lab, d in percol.items():
    for k, v in d.items(): T[f'{k}_{lab}'] = v[r_]
# contrasts scored but absent from both metadata tables
miss = kz.index.difference(meta.key)
if len(miss):
    X = pd.DataFrame(dict(gemma_id=[m.split('__')[0] for m in miss], rsid=[m.split('__')[1] for m in miss], fvid=[m.split('__')[2] for m in miss])); rm = kz.reindex(miss).values.astype(int)
    X['n_genes_scored'] = ng[rm]; X['rho_pH_all'] = rpH[rm]; X['rho_gtex_all'] = rgt[rm]
    for lab, d in percol.items():
        for k, v in d.items(): X[f'{k}_{lab}'] = v[rm]
    T = pd.concat([T, X])
wr(T.sort_values(['gemma_id', 'rsid', 'fvid']), f'{OUT}/per_contrast_table.tsv', ['per-contrast scores: z, AUROC (by |t|), genes present, 97.5th percentile of the 200 comparison-set z (q975), hit (1/0; blank = not counted) for each set; rho_*_all = Spearman of t with the control-pH / GTEx axis over shared genes', 'script: gemma_scans_recurrent_sets/scripts/build_panels.py', 'source: ' + sums, f'date: {today}'])
print('per-contrast rows', len(T))
