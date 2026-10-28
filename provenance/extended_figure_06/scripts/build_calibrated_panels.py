#!/usr/bin/env python3
"""Extended Fig. 11 panels a and c under the calibrated per-contrast rule, plus the caption quantities that depend on the hit rule.
Rule (as ../../gene_correlation_adjusted_tests/addendum3_2026-10-12/fig5h_gemma_notes678/gemma_per_contrast_calibration.py, family 'indep'): a contrast is a hit when z of the 209 genes exceeds the 97.5th percentile of the z of that contrast's own
200 decile-matched independent comparison sets, and the AUROC of the 209 is above 0.5. Contrasts where the 209 or any comparison set has no z (fewer than 5 genes with a prior) are not counted in n.
comparison_set_rate: the same rule applied to each comparison set against the other 199 sets of the contrast (leave one out); the rate is computed per set over the contrasts of the group, and the columns give the mean over the 200 sets
and its 2.5th and 97.5th percentiles. Panel b (median |rho| of t with the pH and GTEx axes) does not use the hit rule and is not rebuilt.
usage: build_calibrated_panels.py SCORED.npz OUTDIR   (scored.npz from score_grouped.py on sets_indep.tsv / contrasts_to_score.tsv; run from anywhere)"""
import sys, os, hashlib, datetime, numpy as np, pandas as pd
NPZ, OUT = sys.argv[1:3]; HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(f'{HERE}/../../..')
FD = os.environ.get('FD', os.path.expanduser('figure_data')); WORK = os.environ.get('WORK', HERE)  # WORK: folder holding sets_indep.tsv and contrasts_to_score.tsv
P_ALL = f'{FD}/gemma_209_all_contrasts/all_contrasts_scored_public.tsv'; P_SCAN = f'{FD}/gemma_176_scan/scan176_contrasts.tsv'; P_COI = f'{FD}/gemma_176_scan_figure/coi_scored.tsv'
P_CLS = os.environ.get('CLS', f'{REPO}/analyses/gemma_176_scan_figure_2026-10-10/remote_followup/flagged_candidates_classified.tsv')
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12]
rd = lambda p: pd.read_csv(p, sep='\t', comment='#', dtype=str, low_memory=False)
key = lambda d: d.gemma_id + '__' + d.rsid + '__' + d.fvid
Zf = np.load(NPZ, allow_pickle=True); sid = Zf['set_id']; j = int(np.where(sid == 'g209')[0][0]); cols = np.array([i for i, s in enumerate(sid) if str(s).startswith('indep_')])
kz = pd.Series(range(len(Zf['gid'])), index=[f'{a}__{b}__{c}' for a, b, c in zip(Zf['gid'], Zf['rsid'], Zf['fvid'])])
z = Zf['z'].astype(float); au = Zf['auroc'].astype(float); Zn = z[:, cols]; z209 = z[:, j]
ok = np.isfinite(z209) & np.isfinite(Zn).all(1)
q = np.full(len(z), np.nan); q[ok] = np.percentile(Zn[ok], 97.5, axis=1)
hit = ok & (z209 > q) & (au[:, j] > 0.5)
loo = np.zeros(Zn.shape, bool)
for c in range(Zn.shape[1]):
    o = np.delete(Zn[ok], c, axis=1); loo[np.where(ok)[0], c] = Zn[ok, c] > np.percentile(o, 97.5, axis=1)
loo &= (au[:, cols] > 0.5)
def wilson(k, n, zc=1.959964):
    if n == 0: return np.nan, np.nan
    p = k / n; d = 1 + zc * zc / n; c = (p + zc * zc / (2 * n)) / d; h = zc * np.sqrt(p * (1 - p) / n + zc * zc / (4 * n * n)) / d; return max(0, c - h), min(1, c + h)
def stats(ix):
    """ix: row indices (into the npz) of the group's contrasts -> dict of k, n, share, wilson, comparison-set rate."""
    ix = np.asarray([i for i in ix if ok[i]]); n = len(ix); k = int(hit[ix].sum()) if n else 0; lo, hi = wilson(k, n)
    sr = loo[ix].mean(0) if n else np.full(len(cols), np.nan)
    return dict(k=k, n=n, share=k / n if n else np.nan, wilson95_lo=lo, wilson95_hi=hi, comparison_set_rate=float(sr.mean()), comparison_set_rate_lo=float(np.percentile(sr, 2.5)), comparison_set_rate_hi=float(np.percentile(sr, 97.5)),
                _sr=sr)
def rows_of(d):  # table with gemma_id/rsid/fvid -> npz row indices (unscored dropped)
    k = key(d); return kz.reindex(k).dropna().astype(int).values
# validation against the published scan of the 209 (original rule), on every scored contrast
A = rd(P_ALL); A['key'] = key(A); A = A.set_index('key'); has = A.index.intersection(kz.index); jj = kz[has].values
zp = A.loc[has, 'z_abs_all_prior'].astype(float).values; ap = A.loc[has, 'auroc_abs_all'].astype(float).values
val = dict(n_compared=len(has), z209_max_abs_diff=float(np.nanmax(np.abs(z209[jj] - zp))), auroc209_max_abs_diff=float(np.nanmax(np.abs(au[jj, j] - ap))))
# ---- population of panels a and c: picked, public, not single-cell rows of scan176_contrasts.tsv (the as-scored population, n = 8,702)
S = rd(P_SCAN); S = S[(S.is_public == 'True') & (S.picked == 'True') & (S.is_single_cell == 'False') & S.n_genes.notna()].copy()
for c in ['in_sz_corpus', 'held_out', 'flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']: S[c] = S[c] == 'True'
S['any_flag'] = S[['flag_pH', 'flag_ischemia', 'flag_postmortem_interval', 'flag_hypoxia']].any(axis=1)
S['row'] = kz.reindex(key(S)).values; n_unscored = int(S.row.isna().sum()); S = S[S.row.notna()].copy(); S['row'] = S.row.astype(int)
S['scored_here'] = ok[S.row.values]
GA = [('Postmortem brain', 'brain_postmortem'), ('Other brain', 'brain_other'), ('Other tissues', 'other'), ('Cell lines / in vitro', 'in_vitro_or_cell_line'), ('Blood', 'blood')]
out_a = [dict(panel='a', group=g, **stats(S.row[S.tissue_class == t].values)) for g, t in GA]
# ---- COI groups, as summarize_coi_209.py
PRIMARY = ['hypoxia', 'ischemia', 'hypoxia_mimetic', 'altitude', 'pH_acidosis', 'PMI_agonal']
Cc = pd.read_csv(P_CLS, sep='\t', dtype={'gemma_id': str, 'fvid': str, 'rsid': str}).fillna(''); Co = rd(P_COI); Co = Co[Co.status == 'scored'].copy()
Co = Co.merge(Cc[['gemma_id', 'tissue_class', 'is_single_cell']].drop_duplicates('gemma_id'), on='gemma_id', how='left'); Co = Co[~Co.is_single_cell.astype(bool)].sort_values(['gemma_id', 'rsid', 'fvid'])
coi_h = Co[Co.cls.isin(PRIMARY)].drop_duplicates('gemma_id'); coi_d = Co[Co.cls == 'disease_ischemia_state'].drop_duplicates('gemma_id')
GC = [('SZ corpus (picked)', S.row[S.in_sz_corpus].values), ('SZ-annotated held out', S.row[(~S.in_sz_corpus) & S.held_out].values), ('All other experiments', S.row[(~S.in_sz_corpus) & ~S.held_out].values),
      ('Design flagged by regex (picked)', S.row[S.any_flag].values), ('Hypoxia / ischemia / pH contrast scored', rows_of(coi_h)), ('Ischemia as disease state (stroke, MI)', rows_of(coi_d))]
out_c = [dict(panel='c', group=g, **stats(r)) for g, r in GC]
# ---- caption quantities
Q = []
def addq(quantity, d, note=''):
    Q.append(dict(quantity=quantity, k=d['k'], n=d['n'], value=d['share'], lo=d['wilson95_lo'], hi=d['wilson95_hi'], note=note))
def addset(quantity, d, note=''):
    sr = d['_sr']; Q.append(dict(quantity=quantity, k='', n=d['n'], value=float(sr.mean()), lo=float(np.percentile(sr, 2.5)), hi=float(np.percentile(sr, 97.5)), note=note + f' min {sr.min():.4f} max {sr.max():.4f}'))
aid = {t: stats(S.row[S.tissue_class == t].values) for _, t in GA}
for g, t in GA: addset(f'comparison sets, per-set hit rate, {g}', aid[t], 'mean over 200 sets; lo, hi = 2.5th and 97.5th percentile of per-set rates;')
brain = aid['brain_postmortem']; Q.append(dict(quantity='comparison sets reaching the 209 postmortem-brain share', k=int((brain['_sr'] >= brain['share']).sum()), n=200, value=float((brain['_sr'] >= brain['share']).mean()), lo='', hi='', note='sets whose per-set brain rate is >= the 209 share'))
# all-contrasts check (all_contrasts_scored_public.tsv: every scored contrast; picks = is_pick)
A2 = rd(P_ALL); A2['row'] = kz.reindex(key(A2)).values; A2 = A2[A2.n_genes.notna() & A2.row.notna()].copy(); A2['row'] = A2.row.astype(int); A2 = A2[ok[A2.row.values]].copy(); A2['hit'] = hit[A2.row.values]; A2['is_pick'] = A2.is_pick == 'True'
exps = set(A2[A2.is_pick].gemma_id); Ae = A2[A2.gemma_id.isin(exps)]
def block(lbl, d):
    pk = d[d.is_pick]; kk = int(d.hit.sum()); n = len(d); lo, hi = wilson(kk, n); em = d.groupby('gemma_id').hit.mean()
    ks, ns = int(pk.hit.sum()), len(pk); lo2, hi2 = wilson(ks, ns)
    Q.append(dict(quantity=f'{lbl}: all scored contrasts, pooled', k=kk, n=n, value=kk / n, lo=lo, hi=hi, note=f'{d.gemma_id.nunique()} experiments'))
    Q.append(dict(quantity=f'{lbl}: mean over experiments of the share of their contrasts that are hits', k='', n=d.gemma_id.nunique(), value=float(em.mean()), lo='', hi='', note='experiments with a scored picked contrast'))
    Q.append(dict(quantity=f'{lbl}: picked contrast per experiment', k=ks, n=ns, value=ks / ns, lo=lo2, hi=hi2, note='is_pick of all_contrasts_scored_public.tsv'))
block('all tissues', Ae); block('postmortem brain', Ae[Ae.tissue_class == 'brain_postmortem'])
for g, t in GA[1:]: block(g, Ae[Ae.tissue_class == t])
sx = Ae[(Ae.tissue_class == 'brain_postmortem') & (Ae.cond_class == 'sex')]; kk, n = int(sx.hit.sum()), len(sx); lo, hi = wilson(kk, n)
Q.append(dict(quantity='postmortem brain, sex contrasts (all candidate contrasts with condition class sex)', k=kk, n=n, value=kk / n, lo=lo, hi=hi, note=f'{sx.gemma_id.nunique()} experiments; mean over experiments {sx.groupby("gemma_id").hit.mean().mean():.4f}'))
sxp = sx[sx.is_pick]; Q.append(dict(quantity='postmortem brain, sex contrasts that are the picked contrast', k=int(sxp.hit.sum()), n=len(sxp), value=sxp.hit.mean() if len(sxp) else np.nan, lo='', hi='', note=''))
# ---- write
stamp = [f'claim: Extended Fig. 11 panels a and c under the calibrated per-contrast rule (209 genes; hit = z of the 209 above the 97.5th percentile of the contrast\'s own null of 200 decile-matched comparison sets, and AUROC above 0.5)',
         f'script: analyses/gemma_scan_null_calibration/calibrated_panels/build_calibrated_panels.py (prepare_inputs.py, score_grouped.py on remote)',
         f'source sha256 (first 12): scored.npz {sha(NPZ)}; sets_indep.tsv {sha(f"{WORK}/sets_indep.tsv")}; contrasts_to_score.tsv {sha(f"{WORK}/contrasts_to_score.tsv")}; scan176_contrasts.tsv {sha(P_SCAN)}; all_contrasts_scored_public.tsv {sha(P_ALL)}; coi_scored.tsv {sha(P_COI)}; flagged_candidates_classified.tsv {sha(P_CLS)}',
         f'date: {datetime.date.today().isoformat()}']
def wr(df, fn, extra=()):
    with open(f'{OUT}/{fn}', 'w') as f:
        for s in stamp + list(extra): f.write('# ' + s + '\n')
        df.to_csv(f, sep='\t', index=False, float_format='%.5g')
cols_out = ['panel', 'group', 'k', 'n', 'share', 'wilson95_lo', 'wilson95_hi', 'comparison_set_rate', 'comparison_set_rate_lo', 'comparison_set_rate_hi']
cmt = ['comparison_set_rate = mean over the 200 comparison sets of the set\'s own hit rate in the group (each set against the other 199 sets of the contrast); _lo, _hi = 2.5th and 97.5th percentile of the per-set rates', 'k, n = contrasts that are hits, contrasts counted (contrasts with no z for the 209 or a comparison set are not counted)']
wr(pd.DataFrame(out_a)[cols_out], 'panel_a.tsv', cmt); wr(pd.DataFrame(out_c)[cols_out], 'panel_c.tsv', cmt)
wr(pd.DataFrame(Q), 'caption_quantities.tsv', ['value = share (or mean per-set rate for comparison-set rows); lo, hi = Wilson 95% (or 2.5th and 97.5th percentile of per-set rates for comparison-set rows)'])
meta = pd.DataFrame([dict(item='contrasts in scored.npz', value=len(z)), dict(item='contrasts with a finite z for the 209 and all 200 comparison sets', value=int(ok.sum())), dict(item='panel a/c population rows not scored (no cached result set, fewer than 2,000 genes)', value=n_unscored),
                     dict(item='panel a/c population rows scored but without a finite z', value=int((~S.scored_here).sum())), dict(item='z209 max abs difference from the published scan', value=val['z209_max_abs_diff']),
                     dict(item='AUROC of the 209 max abs difference from the published scan', value=val['auroc209_max_abs_diff']), dict(item='contrasts compared', value=val['n_compared'])])
wr(meta, 'run_counts.tsv')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 90)
print(meta.to_string()); print(pd.DataFrame(out_a)[cols_out].round(4).to_string()); print(pd.DataFrame(out_c)[cols_out].round(4).to_string()); print(pd.DataFrame(Q).round(4).to_string())
