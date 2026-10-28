"""C15/ADV-005 and C28: Stanley (GSE35977, parietal cortex) within-case analysis of the per-donor recurrent-gene score against lifetime
antipsychotic dose, illness duration, suicide and smoking at death, on the current 209 (set study_brain) and the 1,641 (study_brain_5plus).
Method as adversarial/2026-10-10/brain_expert/stanley_donor_covariates.py and stanley_dose_null.py (176-era): per-donor score = mean over genes of
(SZ direction sign x within-dataset gene z of expression), GSE35977 array, genes expressed above the 10th percentile; OLS, HC3; covariates age, sex, pH, PMI,
plus a BP indicator when SZ and BP are pooled. Null = 2,000 random gene sets of the same size matched gene-by-gene on expression decile x |case-vs-control t| decile,
signed by their own t, same regression (FWL). Donor covariates come from the SMRI workbook match file (read in place, not copied; only aggregates are written).
Run: python3 -I c15_stanley_within_case.py   (about 1 min). Outputs: c15_stanley_within_case.tsv, c15_stanley_null.tsv"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
import numpy as np, pandas as pd, hashlib, statsmodels.formula.api as smf
REPO = 'schizophrenia-transcriptome-replication/'
A = REPO + 'analyses/'
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12]
F_HG = REPO + 'data/reference/hgnc_complete_set.txt'; F_EXP = A + 'confounds/data/gemma/5949_processed.tsv.gz'
F_SETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'recurrent_sets_core80.tsv'); F_M = os.path.expanduser('~/.cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv')  # PATCH: sets path changed on the previous line
hg = pd.read_csv(F_HG, sep='\t', usecols=['symbol', 'entrez_id'], dtype=str).dropna()
e2s = dict(zip(hg.entrez_id.str.replace(r'\.0$', '', regex=True), hg.symbol))
raw = pd.read_csv(F_EXP, sep='\t', comment='#', dtype={'NCBIid': str}, low_memory=False)
raw = raw[raw.NCBIid.notna() & ~raw.NCBIid.astype(str).str.contains(r'[|,]', regex=True)].copy()
raw['sym'] = raw.NCBIid.str.replace(r'\.0$', '', regex=True).map(e2s); raw = raw.dropna(subset=['sym'])
samp = [c for c in raw.columns if '___' in c]
E = raw.groupby('sym')[samp].mean().rename(columns={c: c.split('___')[-1] for c in samp})
E = E.loc[E.mean(axis=1).rank(pct=True) > 0.1]
Z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1), axis=0)
S = pd.read_csv(F_SETS, sep='\t', comment='#')
M = pd.read_csv(F_M, sep='\t', comment='#'); M = M[M.series == 'GSE35977'].set_index('gsm')
M['sz'] = (M.dx_w == 'schizophrenia').astype(int); M['bp'] = (M.dx_w == 'BP').astype(int); M['mdd'] = (M.dx_w == 'MDD').astype(int)
M['logapd'] = np.log10(M.apd_w.clip(lower=1))
def sets():
    out = {}
    for nm, key in [('study_brain', 'study_brain'), ('study_brain_ge10', 'study_brain_ge10'), ('study_brain_ge5', 'study_brain_ge5')]:  # PATCH three final sets
        d = S[S.set == key]; sg = d.set_index('gene').sz_direction.map({'up': 1, 'down': -1})
        sg = sg.dropna(); sg = sg[[x for x in sg.index if x in Z.index]]
        out[nm] = (sg, len(d))
    return out
rows = []; nulls = []
rng = np.random.default_rng(20261011)
# case-vs-control Welch t on z (SZ+BP vs control), used only for the matched null
ctl = M.index[M.dx == 'control'] if False else M.index[M.dx_w == 'control']; cas = M.index[M.dx_w.isin(['schizophrenia', 'BP'])]
Zs = Z[M.index]
tvec = (Zs[cas].mean(axis=1) - Zs[ctl].mean(axis=1)) / np.sqrt(Zs[cas].var(axis=1) / len(cas) + Zs[ctl].var(axis=1) / len(ctl))
dec = pd.qcut(E.mean(axis=1).rank(), 10, labels=False); adec = pd.qcut(tvec.abs().rank(method='first'), 10, labels=False)
key = pd.DataFrame({'e': dec, 'a': adec}).dropna(); key['k'] = key.e.astype(int) * 10 + key.a.astype(int)
pool = {k: list(d.index) for k, d in key.groupby('k')}
for nm, (sg, nlist) in sets().items():
    g = list(sg.index)
    M['score'] = (Z.loc[g].mul(sg[g], axis=0)).mean().reindex(M.index)
    D0 = M.dropna(subset=['score'])
    for grp, dxs in [('SZ+BP', ['schizophrenia', 'BP']), ('SZ only', ['schizophrenia']), ('SZ+BP+MDD (176-era script included MDD)', ['schizophrenia', 'BP', 'MDD'])]:
        CS = D0[D0.dx_w.isin(dxs)].copy()
        bp = {'SZ+BP': ' + bp', 'SZ only': '', 'SZ+BP+MDD (176-era script included MDD)': ' + bp + mdd'}[grp]
        cov = 'age_w + C(sex_w) + pH_w + PMI_w'
        for term, var, sub, scale in [('antipsychotic_log10_per_SD', 'apd_z', 'apd_w', None), ('duration_per_10y', 'dur10', 'dur_w', None),
                                       ('suicide_yes', 'suicide_w', 'suicide_w', None), ('smoking_at_death_yes', 'smoking_w', 'smoking_w', None)]:
            d = CS.dropna(subset=[sub, 'pH_w', 'PMI_w', 'age_w', 'sex_w']).copy()
            d['apd_z'] = (d.logapd - d.logapd.mean()) / d.logapd.std(); d['dur10'] = d.dur_w / 10
            m = smf.ols('score ~ %s%s + %s' % (var, bp, cov), d).fit(cov_type='HC3'); ci = m.conf_int().loc[var]
            rows.append(dict(set=nm, genes_scored=len(g), group=grp, term=term, adjusted_for='age, sex, pH, PMI' + (', diagnosis indicators' if bp else ''), n=int(m.nobs), est=m.params[var], se=m.bse[var], lo=ci[0], hi=ci[1], p=m.pvalues[var]))
        d = CS.dropna(subset=['apd_w', 'dur_w', 'pH_w', 'PMI_w', 'age_w', 'sex_w']).copy()
        d['apd_z'] = (d.logapd - d.logapd.mean()) / d.logapd.std(); d['dur10'] = d.dur_w / 10
        m = smf.ols('score ~ apd_z + dur10%s + %s' % (bp, cov), d).fit(cov_type='HC3'); ci = m.conf_int()
        for var, term in [('apd_z', 'antipsychotic_log10_per_SD (mutual with duration)'), ('dur10', 'duration_per_10y (mutual with antipsychotic)')]:
            rows.append(dict(set=nm, genes_scored=len(g), group=grp, term=term, adjusted_for='age, sex, pH, PMI' + (', diagnosis indicators' if bp else '') + ', other term', n=int(m.nobs), est=m.params[var], se=m.bse[var], lo=ci.loc[var, 0], hi=ci.loc[var, 1], p=m.pvalues[var]))
        # SZ/BP vs control for scale
    for ngrp, ndx, ncols in [('SZ+BP', ['schizophrenia', 'BP'], ['bp']), ('SZ+BP+MDD', ['schizophrenia', 'BP', 'MDD'], ['bp', 'mdd'])]:
        # matched null (SZ+BP only), FWL
        D = M.dropna(subset=['score', 'apd_w', 'dur_w', 'pH_w', 'PMI_w', 'age_w']); D = D[D.dx_w.isin(ndx)].copy()
        D['apd_z'] = (D.logapd - D.logapd.mean()) / D.logapd.std(); D['dur10'] = D.dur_w / 10
        X = np.column_stack([np.ones(len(D))] + [D[c] for c in ncols] + [ D.age_w, (D.sex_w == 'M').astype(float), D.pH_w, D.PMI_w]); H = np.eye(len(D)) - X @ np.linalg.pinv(X)
        ra = H @ D.apd_z.values; rd = H @ D.dur10.values
        Zc = Z[D.index]
        def b(score, r): return (r @ (H @ score)) / (r @ r)
        sk = key.loc[[x for x in g if x in key.index], 'k']
        oa = b(Zc.loc[g].mul(sg[g], axis=0).mean().values, ra); od = b(Zc.loc[g].mul(sg[g], axis=0).mean().values, rd)
        na = []; nd = []
        for _ in range(2000):
            gg = [rng.choice(pool[k]) for k in sk.values]; ss = np.sign(tvec[gg].values)
            s = Zc.loc[gg].mul(ss, axis=0).mean().values; na.append(b(s, ra)); nd.append(b(s, rd))
        na = np.array(na); nd = np.array(nd)
        for term, o, n_ in [('antipsychotic_log10_per_SD', oa, na), ('duration_per_10y', od, nd)]:
            nulls.append(dict(set=nm, group=ngrp, term=term, n=len(D), observed=o, null_mean=n_.mean(), null_sd=n_.std(), null_95lo=np.percentile(n_, 2.5), null_95hi=np.percentile(n_, 97.5), P_null_ge_obs=(1 + (n_ >= o).sum()) / 2001, n_random_sets=2000))
    # context
    cc = M.dropna(subset=['score'])
    rows.append(dict(set=nm, genes_scored=len(g), group='all', term='n_donors_scored; SD of control scores', adjusted_for='', n=len(cc), est=cc[cc.dx_w == 'control'].score.std(), se=np.nan, lo=np.nan, hi=np.nan, p=np.nan))
R = pd.DataFrame(rows); N = pd.DataFrame(nulls)
hdr = ['# script: analyses/clinical_reruns_209_2026-10-11/c15_stanley_within_case.py', '# run: 2026-10-11',
       '# input: analyses/confounds/recurrent_sets.tsv sha256 %s (sets study_brain 209 incl. 1 mixed dropped, study_brain_5plus 1641)' % sha(F_SETS),
       '# input: analyses/confounds/data/gemma/5949_processed.tsv.gz sha256 %s (GSE35977 parietal cortex array)' % sha(F_EXP),
       '# input: ~/.cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv sha256 %s (SMRI workbook match; not redistributable, aggregates only)' % sha(F_M),
       '# input: data/reference/hgnc_complete_set.txt sha256 %s' % sha(F_HG),
       '# score = mean over genes of (SZ direction sign x within-dataset z); apd = lifetime antipsychotics, workbook column apd_w, described in analyses/stanley_demographics_2026-10-08/README.md as fluphenazine-equivalent mg; log10 then per SD of the cases analysed']
for df, fn in [(R, 'c15_stanley_within_case.tsv'), (N, 'c15_stanley_null.tsv')]:
    with open(fn, 'w') as f:
        f.write('\n'.join(hdr) + '\n'); df.round(4).to_csv(f, sep='\t', index=False)
print(R.round(3).to_string()); print(N.round(4).to_string())
