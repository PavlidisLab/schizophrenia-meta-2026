"""Core of the dose / duration tests: the GSE35977 within-case method
(analyses/clinical_reruns_209_2026-10-11/c15_stanley_within_case.py) applied to any expression matrix + case table.

Score   : per sample, mean over the directional genes of the recurrent set (set study_brain of confounds/recurrent_sets.tsv,
          sz_direction up/down; 'mixed' dropped) of (SZ sign x within-dataset gene z). Genes: Gemma NCBI id -> HGNC symbol, probes
          to one symbol averaged, bottom tenth of genes by mean expression dropped, z over all samples of the dataset (or subfield).
Model   : within cases, OLS score ~ exposure + age + sex + pH + PMI (those that exist for the dataset), HC3 errors.
Null    : 2,000 random gene sets of the same size, matched gene by gene on mean-expression decile x |case-vs-control Welch t| decile,
          signed by their own t; same exposure slope by Frisch-Waugh-Lovell; P = (1 + #null >= observed)/(1 + 2000).
Added   : per-gene regressions (as medication_dose.tsv for GSE21138): share of recurrent genes at nominal P < 0.05 against the
          genome-wide share, one-sided Fisher; sample permutation of the exposure among cases; bootstrap over cases for the slope and the
          Spearman correlation of score with exposure.
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
import numpy as np, pandas as pd, hashlib, json, re
import statsmodels.api as sm
from scipy import stats

def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def load_hgnc(path):
    hg = pd.read_csv(path, sep='\t', usecols=['symbol', 'entrez_id'], dtype=str).dropna()
    return dict(zip(hg.entrez_id.str.replace(r'\.0$', '', regex=True), hg.symbol))

def load_sets(path, key='study_brain'):
    S = pd.read_csv(path, sep='\t', comment='#')
    d = S[S.set == key]
    sg = d.set_index('gene').sz_direction.map({'up': 1, 'down': -1}).dropna()
    return sg, len(d)

def gemma_to_symbol_matrix(path, e2s, samples_json=None, platform_path=None):
    """Gemma processed matrix -> symbol x GSM (probes to one symbol averaged)."""
    raw = pd.read_csv(path, sep='\t', comment='#', dtype={'NCBIid': str, 'Probe': str}, low_memory=False)
    if platform_path is not None and raw.NCBIid.isna().all():
        an = pd.read_csv(platform_path, sep='\t', comment='#', usecols=['ElementName', 'NCBIids'], dtype=str)
        raw['NCBIid'] = raw.Probe.map(dict(zip(an.ElementName, an.NCBIids)))
    raw = raw[raw.NCBIid.notna() & ~raw.NCBIid.astype(str).str.contains(r'[|,]', regex=True)].copy()
    raw['sym'] = raw.NCBIid.str.replace(r'\.0$', '', regex=True).map(e2s)
    raw = raw.dropna(subset=['sym'])
    samp = [c for c in raw.columns if '___' in c]
    ren = {c: c.split('___')[-1] for c in samp}
    if not samp:  # older export: BioAssayId -> accession
        sj = json.load(open(samples_json)); sj = sj.get('data', sj)
        acc = {str(s['id']): s['accession']['accession'] for s in sj}
        samp = [c for c in raw.columns[6:]]
        ren = {c: acc[re.sub(r'.*BioAssayId=([0-9]+).*', r'\1', c)] for c in samp}
    E = raw.groupby('sym')[samp].mean().rename(columns=ren)
    return E

def prep(E, gsms):
    """Subset samples, drop NaN genes, drop bottom tenth by mean, z-score per gene over the given samples."""
    E = E[gsms]
    ok = E.isna().mean() < 0.5
    E = E.loc[:, ok]
    E = E.dropna(axis=0)
    E = E.loc[E.mean(axis=1).rank(pct=True) > 0.1]
    Z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1), axis=0)
    return E, Z

def design(d, covs):
    cols = [np.ones(len(d))]; names = ['const']
    for c in covs:
        if c == 'sex':
            s = (d['sex'] == 'M').astype(float).values
            if s.std() > 0: cols.append(s); names.append('sex')
        else:
            v = d[c].astype(float).values
            if v.std() > 0: cols.append(v); names.append(c)
    return np.column_stack(cols), names

def fwl_beta(y, x, X):
    H = np.eye(len(y)) - X @ np.linalg.pinv(X)
    rx = H @ x
    d = rx @ rx
    return (rx @ (H @ y)) / d if d > 1e-12 else np.nan

def run_test(Z, E, cases, var, covs, ctl_ids, case_ids, label, rng, n_null=2000, n_perm=2000, n_boot=2000, sg=None, var_label='', control_ids_for_sd=None):
    """Z, E: symbol x GSM (z-scored; raw filtered). cases: DataFrame indexed by GSM with the exposure column `var` and covariates.
    ctl_ids: GSM of controls (for the matched null's case-vs-control t; cases = all case GSMs in Z)."""
    g = [x for x in sg.index if x in Z.index]
    sgv = sg[g]
    score_all = Z.loc[g].mul(sgv, axis=0).mean()
    D = cases.copy()
    D['score'] = score_all.reindex(D.index)
    need = [var] + [c for c in covs]
    D = D.dropna(subset=['score'] + need)
    n = len(D)
    out = dict(label=label, n_cases_value=n, genes_measured=len(g), n_directional_genes=len(sg))
    X, names = design(D, covs)
    Xf = np.column_stack([X, D[var].values.astype(float)])
    k = Xf.shape[1]
    out['covariates_used'] = '+'.join(names[1:]) if len(names) > 1 else 'none'
    out['resid_df'] = n - np.linalg.matrix_rank(Xf)
    m = sm.OLS(D.score.values, Xf).fit(cov_type='HC3')
    m0 = sm.OLS(D.score.values, Xf).fit()
    ci = m.conf_int()[-1]
    out.update(est=m.params[-1], se_hc3=m.bse[-1], lo_hc3=ci[0], hi_hc3=ci[1], p_hc3=m.pvalues[-1], p_ols=m0.pvalues[-1])
    # control SD of the score
    cs = score_all.reindex(ctl_ids).dropna()
    out['sd_score_controls'] = cs.std()
    out['est_in_control_sd'] = out['est'] / cs.std()
    # Spearman raw and partial (rank-based, residualised on covariates)
    rho, prho = stats.spearmanr(D.score, D[var])
    out['spearman_raw'] = rho; out['spearman_p'] = prho
    # bootstrap over cases
    y = D.score.values; x = D[var].values.astype(float)
    Xc = X
    bs = []; br = []
    for _ in range(n_boot):
        ix = rng.integers(0, n, n)
        if np.ptp(x[ix]) == 0: continue
        Xb = np.column_stack([Xc[ix], x[ix]])
        bb = np.linalg.lstsq(Xb, y[ix], rcond=None)[0][-1]
        bs.append(bb)
        r_ = stats.spearmanr(y[ix], x[ix])[0]
        br.append(r_)
    bs = np.array(bs); br = np.array(br[:len(bs)])
    out['boot_lo'], out['boot_hi'] = np.nanpercentile(bs, [2.5, 97.5])
    out['spearman_boot_lo'], out['spearman_boot_hi'] = np.nanpercentile(br, [2.5, 97.5])
    # sample permutation of the exposure among cases
    H = np.eye(n) - X @ np.linalg.pinv(X)
    ry = H @ y
    obs = (H @ x) @ ry / ((H @ x) @ (H @ x))
    pb = np.empty(n_perm)
    for i in range(n_perm):
        xp = rng.permutation(x); rx = H @ xp
        pb[i] = rx @ ry / (rx @ rx)
    out['perm_p_two_sided'] = (1 + (np.abs(pb) >= abs(obs) - 1e-12).sum()) / (n_perm + 1)
    out['perm_p_positive'] = (1 + (pb >= obs - 1e-12).sum()) / (n_perm + 1)
    # matched random gene sets
    ctl_ids = [c for c in ctl_ids if c in Z.columns]
    cas_all = [c for c in case_ids if c in Z.columns]
    Zs = Z[list(cas_all) + list(ctl_ids)]
    ca = Zs[cas_all]; co = Zs[list(ctl_ids)]
    tvec = (ca.mean(axis=1) - co.mean(axis=1)) / np.sqrt(ca.var(axis=1) / ca.shape[1] + co.var(axis=1) / co.shape[1])
    dec = pd.qcut(E.mean(axis=1).rank(), 10, labels=False)
    adec = pd.qcut(tvec.abs().rank(method='first'), 10, labels=False)
    key = pd.DataFrame({'e': dec, 'a': adec}).dropna(); key['k'] = key.e.astype(int) * 10 + key.a.astype(int)
    pool = {kk: np.array(d.index) for kk, d in key.groupby('k')}
    sk = key.loc[[x_ for x_ in g if x_ in key.index], 'k'].values
    Zc = Z.loc[:, D.index]
    def b(s): 
        return fwl_beta(s, x, X)
    obs_b = b(Zc.loc[g].mul(sgv, axis=0).mean().values)
    nul = np.empty(n_null)
    tv = tvec.to_dict()
    for i in range(n_null):
        gg = [rng.choice(pool[kk]) for kk in sk]
        ss = np.sign([tv[x_] for x_ in gg])
        s = Zc.loc[gg].values.T @ ss / len(gg) if False else (Zc.loc[gg].values * ss[:, None]).mean(axis=0)
        nul[i] = b(s)
    out['matched_obs'] = obs_b; out['matched_null_mean'] = nul.mean(); out['matched_null_sd'] = nul.std()
    out['matched_null_lo'], out['matched_null_hi'] = np.percentile(nul, [2.5, 97.5])
    out['matched_P_null_ge_obs'] = (1 + (nul >= obs_b).sum()) / (n_null + 1)
    out['matched_P_two_sided_abs'] = (1 + (np.abs(nul) >= abs(obs_b)).sum()) / (n_null + 1)
    # per-gene regressions for all genes (cases with a value): t for exposure
    Y = Zc.values.T  # n x genes
    Xg = Xf
    q, r_ = np.linalg.qr(Xg)
    coef = np.linalg.solve(r_, q.T @ Y)
    res = Y - Xg @ coef
    dfr = n - np.linalg.matrix_rank(Xg)
    s2 = (res ** 2).sum(axis=0) / dfr
    XtXi = np.linalg.inv(Xg.T @ Xg)
    tg = coef[-1] / np.sqrt(s2 * XtXi[-1, -1])
    pg = 2 * stats.t.sf(np.abs(tg), dfr)
    gen = pd.DataFrame({'t': tg, 'p': pg}, index=Zc.index)
    sig = gen.p < 0.05
    rate = sig.mean()
    gi = [x_ for x_ in g]
    obs_sig = int(sig[gi].sum()); n_set = len(gi)
    other_sig = int(sig.sum()) - obs_sig
    ft = stats.fisher_exact([[obs_sig, n_set - obs_sig], [other_sig, (~sig).sum() - (n_set - obs_sig)]], alternative='greater')
    same = int((sig[gi] & (np.sign(gen.loc[gi, 't']) == sgv.loc[gi])).sum())
    out.update(genome_rate_p05=rate, n_set_p05=obs_sig, expected_p05=rate * n_set, fisher_p=ft[1], n_p05_dose_follows_sz=same,
               n_p05_dose_against_sz=obs_sig - same, n_all_genes=len(gen))
    # binomial for the share (secondary): P(X >= obs) given genome rate
    out['binom_p_ge'] = stats.binom.sf(obs_sig - 1, n_set, rate)
    out['mean_signed_t_gene'] = float((gen.loc[gi, 't'] * sgv.loc[gi]).mean())
    return out, D

def fmt_p(p):
    return '<0.001' if p < 0.001 else '%.3f' % p
