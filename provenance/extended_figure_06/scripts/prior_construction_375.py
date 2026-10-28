#!/usr/bin/env python3
"""G11 / R190: calibrated per-contrast hit rule for the 277 picked postmortem-brain contrasts under each construction of the DE prior,
for several gene sets at once (g209 = validation, g375, ge10, ge5).
Copy of remote prior_construction_calibrated_2026-10-12/prior_construction_calibrated.py with ONLY these edits (see EDITS.diff):
  - the single 209-gene set is replaced by SETS (names -> gene lists read from recurrent_sets.tsv files); every set is scored in the same pass
    (stacked rows of M) so the leave-out decile vectors are computed once per contrast;
  - comparison sets are drawn per set with a fresh rng(seed) per set (as build_sets_375.py), 200 per set, matched on that set's decile profile
    under the construction's prior, universe = universe_gtex.txt minus the set;
  - 'original comparison sets' job reads the 200 sets per set from gemma_scans_375 out/sets_all.tsv (family indep, group = set name);
  - M is a scipy.sparse CSR matrix (float64) instead of a dense float32 matrix (same sums, speed only);
  - AUROC of the set comes from score_contrast (row of the set) instead of the 209-specific d['sel'] (checked equal for g209);
  - a set whose decile profile exceeds the pool (universe minus set) in some decile under a construction is reported as 'short' (no hit call), not sampled with replacement;
  - output names carry the set. Constructions, seeds, rule and aggregation are unchanged.
usage: OMP_NUM_THREADS=1 nice python prior_construction_375.py OUTDIR NPROC JOBSEL   (JOBSEL = all | pilot)"""
import sys, os, re, json, time, hashlib, numpy as np, pandas as pd
from multiprocessing import Pool
import scipy.sparse as sp
PC = 'axes-handoff/analyses/scan_on_209_2026-10-11/priorcurve'
sys.path.insert(0, f'{PC}/scripts')
import priorlib as PL
from priorlib import E, NR, POOL_ALL, prior_arrays, rank_norm, decile, loo_dec, strat_z, auroc_abs, wilson, _sym2row
OUT = sys.argv[1]; NP = int(sys.argv[2]); JOBSEL = sys.argv[3]
AN = 'axes-handoff/analyses/scan_on_209_2026-10-11/branch/analyses'
NS = 200
RS_NEW = 'sztr-rebuild-infra/runs/final_2026-10-14/forks/profB__current/recurrent_sets.tsv'
RS_OLD = f'{AN}/confounds/recurrent_sets.tsv'
SETS_ALL = 'gemma_scans_recurrent_sets/out/sets_all.tsv'
UNI_F = 'prior_construction_calibrated_2026-10-12/universe_gtex.txt'
D = {}
for gid in E.index:
    fp = f'{PC}/cache/pre/{gid}.npz'
    if os.path.exists(fp):
        z = np.load(fp); D[gid] = {k: z[k] for k in z.files}
SCAN = pd.read_csv('axes-handoff/analyses/scan_on_209_2026-10-11/results/scan209_contrasts.tsv', sep='\t', comment='#', dtype={'gemma_id': str})
SCAN = SCAN[SCAN.picked & ~SCAN.is_single_cell.astype(bool) & SCAN.is_public.astype(bool)]
scored_all = [g for g in SCAN.gemma_id if g in D and 'ra' in D[g]]
TCa = np.array([E.loc[g, 'tissue_class'] for g in scored_all])
scored = [g for g, t in zip(scored_all, TCa) if t == 'brain_postmortem']
pool_all = sorted(g for g in D if g in POOL_ALL)
BR = ('brain_postmortem', 'brain_other')
pool_b = [g for g in pool_all if E.loc[g, 'tissue_class'] in BR]; pool_c = [g for g in pool_all if E.loc[g, 'tissue_class'] not in BR]; B = len(pool_b)
Gn = pd.read_csv(RS_NEW, sep='\t', comment='#'); Go = pd.read_csv(RS_OLD, sep='\t', comment='#')
GENES = {'g209': list(Go[Go.set == 'study_brain'].drop_duplicates('gene').gene),
         'g375': list(Gn[Gn.set == 'study_brain'].drop_duplicates('gene').gene),
         'ge10': list(Gn[Gn.set == 'study_brain_ge10'].drop_duplicates('gene').gene),
         'ge5': list(Gn[Gn.set == 'study_brain_ge5'].drop_duplicates('gene').gene)}
SN = list(GENES)
ROW = {s: np.array(sorted({_sym2row[g] for g in GENES[s] if g in _sym2row}), dtype=np.int64) for s in SN}
UNI = [g.strip() for g in open(UNI_F)]
URALL = {_sym2row[g] for g in UNI if g in _sym2row}
URow = {s: np.array(sorted(URALL - set(ROW[s].tolist())), dtype=np.int64) for s in SN}
SA = pd.read_csv(SETS_ALL, sep='\t')
def cues(g):
    out = set()
    for c in ('sz_bank_cues', 'other_bank_cues'):
        v = E.loc[g, c]
        if isinstance(v, str): out |= {x for x in re.split(r'[,;|\s]+', v.lower()) if x}
    return out
CUE = {g: cues(g) for g in D}; SER = {g: E.loc[g, 'geo_series'] for g in D}
def family(e, S): return [x for x in S if x == e or (SER[x] == SER[e] and isinstance(SER[e], str)) or (CUE[x] & CUE[e])]
PRF = pd.read_csv(f'{AN}/de_prior_lists/our_prior/1_global_our_prior.tsv', sep='\t', comment='#')
DEC_FILE = np.full(NR + 1, -1, dtype=np.int64); v = PRF.de_prior_norm.values; ok_ = np.isfinite(v); DEC_FILE[:NR][ok_] = np.minimum((v[ok_] * 10).astype(int), 9)
def draw_sets(dec, seed):
    """sparse M (float64) with len(SN)*(NS+1) rows; row i<len(SN): the set; then NS comparison rows per set. valid[s] False if a decile pool is short"""
    rr = []; cc = []; valid = {}
    for i, s in enumerate(SN):
        rr.append(np.full(len(ROW[s]), i)); cc.append(ROW[s])
        prof = np.bincount(dec[ROW[s]][dec[ROW[s]] >= 0], minlength=10)
        pools = [URow[s][dec[URow[s]] == d] for d in range(10)]
        short = [d for d in range(10) if len(pools[d]) < prof[d]]
        valid[s] = not short
        if short: continue
        rng = np.random.default_rng(seed)
        for k in range(NS):
            idx = np.concatenate([rng.choice(pools[d], prof[d], replace=False) for d in range(10)])
            rr.append(np.full(len(idx), len(SN) + i * NS + k)); cc.append(idx)
    r = np.concatenate(rr); c = np.concatenate(cc)
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(SN) * (NS + 1), NR)), valid
def sets_from_file():
    rr = []; cc = []
    for i, s in enumerate(SN):
        rr.append(np.full(len(ROW[s]), i)); cc.append(ROW[s])
        T = SA[(SA.group == s) & (SA.family == 'indep')]
        assert len(T) == NS, (s, len(T))
        for k, gs in enumerate(T.genes):
            idx = np.array([_sym2row[g] for g in gs.split(';') if g in _sym2row], dtype=np.int64)
            rr.append(np.full(len(idx), len(SN) + i * NS + k)); cc.append(idx)
    r = np.concatenate(rr); c = np.concatenate(cc)
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(SN) * (NS + 1), NR)), {s: True for s in SN}
def score_contrast(d, M, dec):
    prow = d['prow']; ra = d['ra'].astype(np.float64); n = len(ra); okp = prow >= 0
    rk = np.zeros(NR); pres = np.zeros(NR, np.float32); rk[prow[okp]] = ra[okp]; pres[prow[okp]] = 1
    sdec = dec[:NR]; okd = (sdec >= 0) & (pres > 0)
    ONE = np.zeros((NR, 10), np.float64); ONE[np.where(okd)[0], sdec[okd]] = 1
    Pm = M.multiply(pres.astype(np.float64)[None, :]).tocsr(); npres = np.asarray(Pm.sum(1)).ravel(); K = np.asarray(Pm @ ONE); n1 = K.sum(1)
    rd_ = rk * okd; obs = np.asarray(Pm @ rd_).ravel() / np.maximum(n1, 1)
    N = ONE.sum(0).astype(float); S1 = ONE.T @ rk; S2 = ONE.T @ (rk * rk)
    with np.errstate(divide='ignore', invalid='ignore'):
        mean = np.where(N > 0, S1 / N, 0.0); var = np.where(N > 0, S2 / N - mean ** 2, 0.0)
        vt = (K * var * np.where(N > 1, (N - K) / np.maximum(N - 1, 1), 0.0)).sum(1)
        exp = (K * mean).sum(1) / np.maximum(n1, 1); sd = np.sqrt(vt) / np.maximum(n1, 1)
        z = np.where(sd > 0, (obs - exp) / sd, np.nan)
        z[n1 < 5] = np.nan
        au = (np.asarray(Pm @ rk).ravel() - npres * (npres + 1) / 2) / (npres * (n - npres))
    return z, au
def run(a):
    name, S, mode, seed, fromfile = a
    Ss = set(S) if S is not None else set()
    if S is None: dec = DEC_FILE; nt = sc = None
    else: nt, sc = prior_arrays(D, S); dec = decile(rank_norm(sc, nt))
    M, valid = sets_from_file() if fromfile else draw_sets(dec, seed)
    R = {s: dict(okc=np.zeros(len(scored), bool), hit=np.zeros(len(scored), bool), loo=np.zeros((len(scored), NS), bool), z=np.full(len(scored), np.nan), q=np.full(len(scored), np.nan), au=np.full(len(scored), np.nan)) for s in SN}
    for i, g in enumerate(scored):
        d = D[g]
        if S is not None and g in Ss and mode != 'none':
            drop = family(g, S) if mode == 'family' else [g]; dd = loo_dec(D, g, nt, sc, drop)
        else: dd = dec
        z, au = score_contrast(d, M, dd)
        for j, s in enumerate(SN):
            if not valid[s]: continue
            zi = z[j]; Zn = z[len(SN) + j * NS: len(SN) + (j + 1) * NS]; aus = au[len(SN) + j * NS: len(SN) + (j + 1) * NS]
            if not (np.isfinite(zi) and np.isfinite(Zn).all()): continue
            r = R[s]; r['okc'][i] = True; r['z'][i] = zi; r['au'][i] = au[j]; r['q'][i] = np.percentile(Zn, 97.5); r['hit'][i] = (zi > r['q'][i]) and (au[j] > 0.5)
            for c in range(NS):
                o = np.delete(Zn, c); r['loo'][i, c] = (Zn[c] > np.percentile(o, 97.5)) and (aus[c] > 0.5)
    return name, R, valid
if __name__ == '__main__':
    t0 = time.time()
    print('scored postmortem brain', len(scored), 'pools', len(pool_all), len(pool_b), len(pool_c), {s: (len(ROW[s]), len(URow[s])) for s in SN}, flush=True)
    jobs = [('project static prior, original comparison sets', None, 'none', None, True),
            ('project static prior, re-drawn comparison sets', None, 'none', [20261012, 0], False),
            ('all-experiment prior, leave-one-out', pool_all, 'loo', [20261012, 1], False),
            ('all-experiment prior, whole family left out', pool_all, 'family', [20261012, 1], False),
            ('brain-only prior (k=351), leave-one-out', pool_b, 'loo', [20261012, 2], False),
            ('brain-only prior, no leave-out (self included)', pool_b, 'none', [20261012, 2], False),
            ('brain-only prior, whole family left out', pool_b, 'family', [20261012, 2], False)]
    nd = 100 if JOBSEL == 'all' else 2
    for d in range(nd):
        jobs.append((f'sm_all_{d}', list(np.random.default_rng([20261011, 1, d]).choice(pool_all, B, replace=False)), 'loo', [20261012, 3, d], False))
        jobs.append((f'sm_nonbrain_{d}', list(np.random.default_rng([20261011, 2, d]).choice(pool_c, B, replace=False)), 'loo', [20261012, 4, d], False))
    res = {}; valid_all = {}
    ck = f'{OUT}/ckpt'; os.makedirs(ck, exist_ok=True)
    todo = []
    for j in jobs:
        fp = f'{ck}/{hashlib.md5(j[0].encode()).hexdigest()}.npz'
        if os.path.exists(fp):
            z = np.load(fp, allow_pickle=True); res[j[0]] = z['R'].item(); valid_all[j[0]] = z['valid'].item(); print('ckpt', j[0], flush=True)
        else: todo.append(j)
    with Pool(NP) as P:
        for nm, R, valid in P.imap_unordered(run, todo, chunksize=1):
            res[nm] = R; valid_all[nm] = valid
            np.savez(f'{ck}/{hashlib.md5(nm.encode()).hexdigest()}.npz', R=np.array(R, dtype=object), valid=np.array(valid, dtype=object))
            print(nm, {s: (int(R[s]['okc'].sum()), int(R[s]['hit'].sum())) for s in SN}, round(time.time() - t0), flush=True)
    for s in SN:
        rows = []
        def agg(label, names, desc):
            names = [n for n in names if valid_all[n][s]]
            if not names:
                rows.append(dict(construction=label, k=np.nan, n=0, share=np.nan, wilson95_lo=np.nan, wilson95_hi=np.nan, comparison_set_rate=np.nan, draws=0, draw_share_p2_5=np.nan, draw_share_p97_5=np.nan, description=desc + ' [decile pool short: not scored]')); return
            n_ok = [int(res[n][s]['okc'].sum()) for n in names]; ks = [int(res[n][s]['hit'].sum()) for n in names]
            sh = np.array([k / n for k, n in zip(ks, n_ok)]); n = int(round(np.mean(n_ok))); share = float(sh.mean()); k = share * n
            lo, hi = wilson(k, n)
            csr = np.mean([res[nm][s]['loo'][res[nm][s]['okc']].mean(0).mean() for nm in names])
            rows.append(dict(construction=label, k=round(k, 1) if len(names) > 1 else int(k), n=n, share=share, wilson95_lo=lo, wilson95_hi=hi, comparison_set_rate=float(csr),
                             draws=len(names), draw_share_p2_5=float(np.percentile(sh, 2.5)), draw_share_p97_5=float(np.percentile(sh, 97.5)), description=desc + (f' [{len(names)} of {len(names)+0} draws]' if False else '')))
        for lab, desc in [('project static prior, original comparison sets', 'project prior deciles, the 200 sets of gemma_scans_375 sets_all.tsv (default calibrated construction)'),
                          ('project static prior, re-drawn comparison sets', 'project prior deciles, 200 sets re-drawn by this script'),
                          ('all-experiment prior, leave-one-out', 'prior from the 8,712-experiment pool; each scored experiment scored with itself left out; sets matched on the full-pool deciles'),
                          ('all-experiment prior, whole family left out', 'as above with the experiment and its GEO-series / bank-cue family left out'),
                          ('brain-only prior (k=351), leave-one-out', 'prior from the 351 brain experiments of the pool; leave-one-out'),
                          ('brain-only prior, no leave-out (self included)', 'brain-only prior, scored experiment included'),
                          ('brain-only prior, whole family left out', 'brain-only prior, family left out')]:
            agg(lab, [lab], desc)
        agg('size-matched all-tissue prior (k=351), leave-one-out', [f'sm_all_{d}' for d in range(nd)], f'{nd} random priors of 351 experiments from the pool; sets re-drawn per prior; mean over draws')
        agg('size-matched non-brain prior (k=351), leave-one-out', [f'sm_nonbrain_{d}' for d in range(nd)], f'{nd} random priors of 351 non-brain experiments; sets re-drawn per prior; mean over draws')
        pd.DataFrame(rows).to_csv(f'{OUT}/prior_construction_raw_{s}.tsv', sep='\t', index=False, float_format='%.6g')
    np.savez_compressed(f'{OUT}/per_contrast_hits_{JOBSEL}.npz', gid=np.array(scored), **{f'hit__{s}__{k}': v[s]['hit'] for k, v in res.items() for s in SN}, **{f'ok__{s}__{k}': v[s]['okc'] for k, v in res.items() for s in SN})
    json.dump({k: v for k, v in valid_all.items()}, open(f'{OUT}/valid_{JOBSEL}.json', 'w'))
    print('done', time.time() - t0)
