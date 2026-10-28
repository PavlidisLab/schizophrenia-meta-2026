#!/usr/bin/env python3
"""(Copy for the 209 redo, 2026-10-11: only change is the sz_gene_counts column in_176 -> in_recurrent_set; the arrays named sz176 / overlap_sz176 now hold the 209. Run in a COPY of this folder outside sztr-516.)
Empirical recurrence null from random public postmortem-brain contrasts in the Gemma DE cache (remote).
Run on remote:  python3 run_null.py prepare ; python3 run_null.py draw [ndraws] [nproc]
Inputs (copied alongside): sz_lists_62.tsv, sz_gene_counts.tsv (export_sz_inputs.py), postmortem_brain_contrasts_classified.tsv
(analyses/gemma_176_scan_nonhits_2026-10-10), 1_global_our_prior.tsv (analyses/de_prior_lists/our_prior), and the read-only cache
axes-handoff/analyses/gemma_176_scan_2026-10-10/cache/{de_cache.h5,h5_index.tsv}. No network, no model calls."""
import sys, os, json, hashlib, time, pickle, itertools, datetime
import numpy as np, pandas as pd, h5py
from multiprocessing import Pool
from scipy.stats import rankdata
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = 'axes-handoff/analyses/gemma_176_scan_2026-10-10/cache'
sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()[:12]
MINGENES = 2000
SEED = 20261010

def prepare():
    cl = pd.read_csv(f'{HERE}/postmortem_brain_contrasts_classified.tsv', sep='\t', comment='#')
    ix = pd.read_csv(f'{CACHE}/h5_index.tsv', sep='\t')
    ix = ix[ix.picked == True].drop_duplicates('gemma_id')
    ex = pd.read_csv('sztr-rebuild-infra/repo/analyses/de_prior_lists/experiments.tsv', sep='\t', comment='#')[['gemma_id', 'sz_bank_cues', 'is_public']]
    cl = cl.merge(ex, on='gemma_id', how='left'); assert cl.is_public.all()
    m = cl.merge(ix[['gemma_id', 'name', 'n_genes']].rename(columns={'name': 'key'}), on='gemma_id', how='left')
    print('classified', len(cl), 'with h5 key', m.key.notna().sum())
    m = m[m.key.notna() & (m.n_genes >= MINGENES)].reset_index(drop=True)
    pri = pd.read_csv(f'{HERE}/1_global_our_prior.tsv', sep='\t')
    f = h5py.File(f'{CACHE}/de_cache.h5', 'r')['contrasts']
    ids = {}
    recs = []
    for r in m.itertuples():
        g = f[r.key]
        gid = g['gene_ncbi_id'][:]; p = g['pvalue'][:].astype(float); t = g['tstat'][:].astype(float); q = g['qvalue'][:].astype(float)
        sym = g['symbol'][:]
        ok = np.isfinite(p) & np.isfinite(t) & (gid > 0)
        gid, p, t, q, sym = gid[ok], p[ok], t[ok], q[ok], sym[ok]
        o = np.lexsort((-np.abs(t), p))
        for a, s in zip(gid, sym):
            if a not in ids: ids[a] = s.decode() if isinstance(s, bytes) else str(s)
        recs.append(dict(key=r.key, gid=gid[o].astype(np.int64), up=(t[o] > 0), n=len(gid), n_fdr=int((q < 0.05).sum()), n_nom=int((p < 0.05).sum())))
    gids = np.array(sorted(ids)); pos = {a: i for i, a in enumerate(gids)}
    for rc in recs: rc['idx'] = np.array([pos[a] for a in rc['gid']], dtype=np.int32); del rc['gid']
    sym = np.array([ids[a] for a in gids])
    pmap = pri.set_index('ncbi_gene_id')
    prior_n = np.array([pmap.n_de_total.get(a, np.nan) for a in gids]); prior_norm = np.array([pmap.de_prior_norm.get(a, np.nan) for a in gids])
    meta = m.drop(columns=['n_genes']).copy(); meta['n_tested'] = [r['n'] for r in recs]; meta['n_fdr'] = [r['n_fdr'] for r in recs]; meta['n_nom'] = [r['n_nom'] for r in recs]
    pickle.dump(dict(recs=recs, gids=gids, sym=sym, prior_n=prior_n, prior_norm=prior_norm, meta=meta), open(f'{HERE}/prepared.pkl', 'wb'))
    meta.to_csv(f'{HERE}/_pool_meta.tsv', sep='\t', index=False)
    print('pool', len(meta), 'genes', len(gids), 'in prior', np.isfinite(prior_n).sum())

def pops(meta):
    sz = (meta.sz_status != 'neither') | meta.disorder.str.startswith('schizo')
    P = {}
    P['P1_disease'] = meta.index[(meta.true_disease) & ~sz].values
    P['P1c_disease_case_like'] = meta.index[(meta.true_disease) & (meta.case_like) & ~sz].values
    P['P2_all_postmortem_brain'] = meta.index[~sz].values
    for base in ('P1_disease', 'P2_all_postmortem_brain'):   # power-favourable sensitivity: upper half of the pool by number of nominal (p < 0.05) genes
        ix = P[base]; thr = np.median(meta.n_nom.values[ix]); P[base + '_top_half_signal'] = ix[meta.n_nom.values[ix] >= thr]
    return P

def auc_membership(score, member):
    """AUC of score ranking member genes above others (ties averaged), over genes with a finite score."""
    ok = np.isfinite(score); s = score[ok]; mem = member[ok]
    n1 = mem.sum(); n0 = len(mem) - n1
    if n1 == 0 or n0 == 0: return np.nan
    r = rankdata(s)
    return (r[mem].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def metrics(cnt, upc, dnc, D):
    """cnt/upc/dnc: per-gene counts over D['gids'] universe."""
    out = {}
    for k in (5, 10):
        out[f'n{k}'] = int((cnt >= k).sum())
        out[f'n{k}_common'] = int((cnt[D['common']] >= k).sum())
        mx = np.maximum(upc, dnc)
        out[f'n{k}_dir'] = int((mx >= k).sum())
    r5 = cnt >= 5
    out['n_ge2'] = int((cnt >= 2).sum()); out['n_ge3'] = int((cnt >= 3).sum())
    tot = (upc + dnc)[r5].sum(); out['major_frac_5'] = float(np.maximum(upc, dnc)[r5].sum() / tot) if tot > 0 else np.nan
    out['overlap_sz1478'] = int((r5 & D['sz5']).sum()); out['frac_in_sz1478'] = out['overlap_sz1478'] / max(r5.sum(), 1)
    out['overlap_sz176'] = int((r5 & D['sz176']).sum())
    out['exp_overlap_sz1478'] = r5.sum() * D['sz5'][D['common']].sum() / D['common'].sum()
    pn = D['prior_n_f']; out['auc_prior_5'] = auc_membership(pn, r5)
    pnorm = D['prior_norm']; v = pnorm[r5 & np.isfinite(pnorm)]
    out['mean_prior_pct_5'] = float(v.mean()) if len(v) else np.nan
    out['frac_top_decile_prior_5'] = float((v >= 0.9).mean()) if len(v) else np.nan
    return out

def slot_lists(rule, rng, pool, sz, D, forced=None):
    """Return list of (contrast index in recs, take-n) for 62 slots; None if impossible."""
    recs = D['recs']; sizes = sz.list_size_in_universe.values.copy(); rates = sz.rate.values.copy(); kinds = sz.threshold_kind.values
    n_slots = len(sz)
    perm = rng.permutation(n_slots)
    if rule in ('A', 'A2'):
        pick = list(forced) if forced is not None else []
        rest = [i for i in pool if i not in set(pick)]
        if len(rest) < n_slots - len(pick): return None  # rebuild 2026-10-14: population smaller than the number of SZ lists (74): combination skipped
        pick += list(rng.choice(rest, n_slots - len(pick), replace=False))
        pick = np.array(pick); rng.shuffle(pick)
        out = []
        for j, c in enumerate(pick):
            n = recs[c]['n']
            L = sizes[perm[j]] if rule == 'A' else int(round(rates[perm[j]] * n))
            out.append((c, int(min(max(L, 1), n))))
        return out
    nfdr = np.array([recs[c]['n_fdr'] for c in pool]); poolarr = np.array(pool)
    if rule == 'B':
        cand = poolarr[nfdr >= 1]
        if len(cand) < n_slots: return None
        pick = rng.choice(cand, n_slots, replace=False)
        return [(c, recs[c]['n_fdr']) for c in pick]
    if rule == 'B2':   # FDR < 0.05 hits, capped at the size drawn from the SZ list sizes
        cand = poolarr[nfdr >= 1]
        if len(cand) < n_slots: return None
        pick = rng.choice(cand, n_slots, replace=False)
        return [(c, int(min(recs[c]['n_fdr'], sizes[perm[j]]))) for j, c in enumerate(pick)]
    if rule == 'N':
        pick = rng.choice(poolarr, n_slots, replace=False)
        return [(c, recs[c]['n_nom']) for c in pick]
    if rule == 'C':   # SZ threshold mix: slots of kind fdr take FDR<0.05 lists (from contrasts with >= 1 hit), the rest nominal p<0.05
        nf = int((kinds == 'fdr').sum()); cand = poolarr[nfdr >= 1]
        if len(cand) < nf: return None
        pf = rng.choice(cand, nf, replace=False)
        rest = np.setdiff1d(poolarr, pf); pn = rng.choice(rest, n_slots - nf, replace=False)
        return [(c, recs[c]['n_fdr']) for c in pf] + [(c, recs[c]['n_nom']) for c in pn]

def one_combo(args):
    popname, rule, ndraws, forced_name = args
    D = G; meta = D['meta']
    P = pops(meta); pool = list(P[popname])
    forced = None
    if forced_name == 'bankblock':
        bank = meta.index[(meta.sz_bank_cues.notna())].values if 'sz_bank_cues' in meta else []
        forced = [b for b in bank if b in set(pool)]
    sz = D['sz']
    rng = np.random.default_rng(SEED + int(hashlib.md5(f'{popname}{rule}{forced_name}'.encode()).hexdigest()[:6], 16))
    Ng = len(D['gids']); rows = []
    for d in range(ndraws):
        try: sl = slot_lists(rule, rng, pool, sz, D, forced)
        except ValueError: sl = None  # rebuild 2026-10-14: population smaller than the number of SZ lists (74): combination skipped
        if sl is None: return (popname, rule, forced_name, None)
        cnt = np.zeros(Ng, dtype=np.int32); upc = np.zeros(Ng, dtype=np.int32); dnc = np.zeros(Ng, dtype=np.int32)
        sizes = []
        for c, L in sl:
            r = D['recs'][c]; ii = r['idx'][:L]; up = r['up'][:L]
            cnt[ii] += 1; upc[ii[up]] += 1; dnc[ii[~up]] += 1; sizes.append(L)
        m = metrics(cnt, upc, dnc, D)
        m.update(draw=d, pop=popname, rule=rule, variant=forced_name or 'random', median_list_size=float(np.median(sizes)), mean_list_size=float(np.mean(sizes)), n_empty=int(sum(1 for s in sizes if s == 0)))
        rows.append(m)
        if d < 3 and False: pass
    return (popname, rule, forced_name, pd.DataFrame(rows))

G = None
def init(path):
    global G
    D = pickle.load(open(f'{HERE}/prepared.pkl', 'rb'))
    szc = pd.read_csv(f'{HERE}/sz_gene_counts.tsv', sep='\t', comment='#')
    szs = pd.read_csv(f'{HERE}/sz_lists_62.tsv', sep='\t', comment='#')
    symidx = {}
    for i, s in enumerate(D['sym']): symidx.setdefault(s, i)
    Ng = len(D['gids']); sz5 = np.zeros(Ng, bool); sz176 = np.zeros(Ng, bool); szmask = np.zeros(Ng, bool)
    for r in szc.itertuples():
        i = symidx.get(r.symbol)
        if i is not None:
            szmask[i] = True
            if r.n_lists >= 5: sz5[i] = True
            if r.in_recurrent_set: sz176[i] = True
    D['common'] = szmask; D['sz5'] = sz5; D['sz176'] = sz176; D['sz'] = szs
    pn = D['prior_n'].copy(); med = np.nanmedian(pn); pn[~np.isfinite(pn)] = med; D['prior_n_f'] = pn   # absent from prior -> median count, as in the project's DE-prior weights
    G = D

def sz_observed():
    szc = pd.read_csv(f'{HERE}/sz_gene_counts.tsv', sep='\t', comment='#')
    pri = pd.read_csv(f'{HERE}/1_global_our_prior.tsv', sep='\t').drop_duplicates('gene_symbol').set_index('gene_symbol')
    cnt = szc.n_lists.values; up = szc.n_up.values; dn = szc.n_down.values
    r5 = cnt >= 5; out = {}
    for k in (5, 10): out[f'n{k}'] = int((cnt >= k).sum()); out[f'n{k}_dir'] = int((np.maximum(up, dn) >= k).sum())
    out['n_ge2'] = int((cnt >= 2).sum()); out['n_ge3'] = int((cnt >= 3).sum())
    tot = (up + dn)[r5].sum(); out['major_frac_5'] = float(np.maximum(up, dn)[r5].sum() / tot)
    pn = pd.to_numeric(szc.de_prior_count, errors='coerce').values.astype(float); pn[~np.isfinite(pn)] = np.nanmedian(pn)
    out['auc_prior_5'] = auc_membership(pn, r5)
    pnorm = np.array([pri.de_prior_norm.get(s, np.nan) for s in szc.symbol]); v = pnorm[r5 & np.isfinite(pnorm)]
    out['mean_prior_pct_5'] = float(v.mean()); out['frac_top_decile_prior_5'] = float((v >= 0.9).mean()); out['n_universe'] = len(szc)
    return out

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'prepare': prepare()
    else:
        nd = int(sys.argv[2]) if len(sys.argv) > 2 else 2000; npr = int(sys.argv[3]) if len(sys.argv) > 3 else 12
        init(None)
        json.dump(dict(n_common=int(G['common'].sum()), sz5_common=int(G['sz5'][G['common']].sum()), sz176_common=int(G['sz176'][G['common']].sum()), n_universe_h5=int(len(G['gids']))), open(f'{HERE}/_common_counts.json', 'w'))
        P = pops(G['meta']); print({k: len(v) for k, v in P.items()}, flush=True)
        combos = [(p, r, nd, None) for p in P for r in ('A', 'A2', 'B2', 'B', 'N', 'C')]
        t = time.time()
        with Pool(npr, initializer=init, initargs=(None,)) as pool:
            res = pool.map(one_combo, combos, chunksize=1)
        dfs = []
        for p, r, v, df in res:
            print(p, r, v, 'skipped (pool too small)' if df is None else len(df))
            if df is not None: dfs.append(df)
        allr = pd.concat(dfs); allr.to_csv(f'{HERE}/_draws_raw.tsv', sep='\t', index=False)
        json.dump(sz_observed(), open(f'{HERE}/_sz_observed.json', 'w'))
        print('done', round(time.time() - t), 's')
