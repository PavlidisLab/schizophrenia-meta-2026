"""Stanley (GSE35977 parietal cortex, GSE35974 cerebellum): the recurrent-gene score (209 set) against every per-donor variable the SMRI design file / workbooks hold.
Model of the smoking analysis in Extended Methods 12 (c15_stanley_within_case.py): within cases, OLS HC3, score ~ variable + age + sex + pH + PMI (+ BP indicator when SZ and BP are pooled).
Matched null: 2,000 random gene sets matched gene by gene on expression decile x |SZ+BP-vs-control Welch t| decile, signed by own t; slope by Frisch-Waugh-Lovell.
Donors: GEO array -> SMRI A/C code by missed_candidates_2026-10-11/stanley_donor_level/<GSE>_sample_to_acode.tsv (status unique only) -> expdesign.chen.txt row 'chen_A-n' / 'chen_C-n'.
Duration of illness is not in the design file: read from ArrayCollection.xlsx (Acode) / NeuropathologyConsortium.xlsx (Sample), joined after checking age, sex, PMI, pH agree.
Reads donor tables in place (SMRI terms); writes aggregates only. Run: python3 -I -W ignore stanley_all.py <repo> <out_dir>"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dd_core import *
REPO = os.path.expanduser(sys.argv[1]); OUT = os.path.expanduser(sys.argv[2]); os.makedirs(OUT, exist_ok=True)
MAPSRC = sys.argv[3] if len(sys.argv) > 3 else 'matched'   # matched = .cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv (the donors of c15); acode = missed_candidates sample_to_acode (unique rows only)
A = REPO + '/analyses/'; H = lambda p: p.replace(os.path.expanduser('~'), '~')
F = dict(hgnc=REPO + '/data/reference/hgnc_complete_set.txt', sets=(os.environ.get('F_SETS') or A + 'confounds/recurrent_sets.tsv'), design=REPO + '/data/raw/mistry_server/stanley/expdesign.chen.txt',
         m77=A + 'missed_candidates_2026-10-11/stanley_donor_level/GSE35977_sample_to_acode.tsv', m74=A + 'missed_candidates_2026-10-11/stanley_donor_level/GSE35974_sample_to_acode.tsv',
         x77=A + 'confounds/data/gemma/5949_processed.tsv.gz', x74=A + 'confounds/data/gemma/5939_processed.tsv.gz',
         arr=os.path.expanduser('Downloads/ArrayCollection.xlsx'), neu=os.path.expanduser('Downloads/NeuropathologyConsortium.xlsx'))
e2s = load_hgnc(F['hgnc']); sg, nl = load_sets(F['sets'])
D = pd.read_csv(F['design'], sep='\t', comment='#')
D = D.rename(columns={'Age': 'age', 'Sex': 'sex', 'Brain PH': 'pH', 'PMI (h)': 'PMI', 'Lifetime Anti-psychotics': 'apd', 'Suicide Status': 'suicide', 'Smoking at Time of Death': 'smoking',
                      'Rate of Death': 'rod', 'Lifetime Alcohol': 'alcohol', 'Lifetime Drugs': 'drugs', 'Profile': 'profile'})
D['key'] = D.ID.str.replace('chen_', '', regex=False)  # 'A-51', 'C-12'
dxmap = {'Schizophrenia': 'schizophrenia', 'Bipolar': 'BP', 'Depression': 'MDD', 'Control': 'control'}
D['dx'] = D.profile.map(dxmap)
# duration from the workbooks
Ar = pd.read_excel(F['arr']); Ar = Ar[Ar.Acode.astype(str).str.match(r'^A\d+$')].copy(); Ar['key'] = Ar.Acode.str.replace('A', 'A-', regex=False)
Ne = pd.read_excel(F['neu']); Ne['key'] = Ne.Sample.astype(str)
W = pd.concat([Ar[['key', 'DurationOfIllness', 'Age', 'PMI', 'BrainPH', 'LifetimeAntipsychotics']].rename(columns={'DurationOfIllness': 'dur', 'BrainPH': 'pH', 'LifetimeAntipsychotics': 'apd', 'Age': 'age', 'PMI': 'PMI'}),
               Ne[['key', 'Duration of Illness', 'Age', 'PMI', 'Brain pH', 'Antipsychotic']].rename(columns={'Duration of Illness': 'dur', 'Brain pH': 'pH', 'Antipsychotic': 'apd', 'Age': 'age'})]).set_index('key')
chk = D.set_index('key').join(W, rsuffix='_w', how='left')
print('design rows', len(D), 'with a workbook row', int(chk.age_w.notna().sum()), '; design age equals workbook age in', int((chk.age == chk.age_w).sum()), '; |diff| <= 2 in', int(((chk.age - chk.age_w).abs() <= 2).sum()),
      '; PMI equal', int((chk.PMI == chk.PMI_w).sum()), '; pH equal', int(((chk.pH - chk.pH_w).abs() < 1e-9).sum()), '; lifetime antipsychotics equal', int((pd.to_numeric(chk.apd, errors='coerce') == pd.to_numeric(chk.apd_w, errors='coerce')).sum()))
# age, PMI, pH, lifetime antipsychotics and duration from the workbooks (as c15 used workbook values); sex, rate of death, suicide, smoking, alcohol, drugs from the design file
for c in ['age', 'PMI', 'pH', 'apd']:
    D[c] = D.key.map(W[c]).where(D.key.isin(W.index), D[c]) if c != 'apd' else D.key.map(W[c])
D['dur'] = pd.to_numeric(D.key.map(W['dur']), errors='coerce')
D['apd'] = pd.to_numeric(D.apd, errors='coerce'); D['logapd'] = np.log10(D.apd.clip(lower=1))
D['sex'] = D.sex.astype(str)
D['bp'] = (D.dx == 'BP').astype(int)
D['slowdeath'] = np.where(D.rod.isna(), np.nan, (D.rod != 'Sudden').astype(float))
D['suicide_y'] = D.suicide.map({'Yes': 1.0, 'No': 0.0}); D['smoke_y'] = D.smoking.map({'Yes': 1.0, 'No': 0.0})
for v, c in [('alcohol', 'alc'), ('drugs', 'drug')]:
    D[c + '_heavy'] = np.where(D[v].isna(), np.nan, D[v].fillna('').str.startswith('Heavy').astype(float))
    D[c + '_modheavy'] = np.where(D[v].isna(), np.nan, D[v].fillna('').str.match(r'^(Heavy|Moderate)').astype(float))
D = D.set_index('key')
# inventory over all 165 donors by diagnosis
inv = []
for var in ['rod', 'smoking', 'suicide', 'alcohol', 'drugs']:
    for dx in ['schizophrenia', 'BP', 'MDD', 'control']:
        s = D[D.dx == dx][var]
        for v, n in s.value_counts(dropna=False).items(): inv.append(dict(scope='all 165 donors', variable=var, diagnosis=dx, value='(blank)' if pd.isna(v) else v, n=int(n)))
for dx in ['schizophrenia', 'BP', 'MDD', 'control']:
    s = D[D.dx == dx]
    inv.append(dict(scope='all 165 donors', variable='lifetime antipsychotics (number with a value)', diagnosis=dx, value='numeric', n=int(s.apd.notna().sum())))
    inv.append(dict(scope='all 165 donors', variable='duration of illness (number with a value)', diagnosis=dx, value='numeric', n=int(s.dur.notna().sum())))
rng = np.random.default_rng(20261014)
VARS = [('lifetime antipsychotics (log10, per SD)', 'logapd', 'apd'), ('duration of illness (per 10 years)', 'dur10', 'dur'), ('rate of death: not sudden vs sudden', 'slowdeath', None), ('suicide', 'suicide_y', None),
        ('smoking at time of death', 'smoke_y', None), ('alcohol: heavy (past or present)', 'alc_heavy', None), ('alcohol: moderate or heavy', 'alc_modheavy', None),
        ('drugs: heavy (past or present)', 'drug_heavy', None), ('drugs: moderate or heavy', 'drug_modheavy', None)]
rows = []; prov = {}
for ds, mf, xf in [('GSE35977', 'm77', 'x77'), ('GSE35974', 'm74', 'x74')]:
    Eall = gemma_to_symbol_matrix(F[xf], e2s)
    if MAPSRC == 'acode':
        M = pd.read_csv(F[mf], sep='\t', comment='#'); U = M[(M.status == 'unique') & M.acode.notna()].copy()
        U['key'] = U.acode.astype(str).str.replace(r'^([AC])(\d+)$', r'\1-\2', regex=True)
        U = U[U.gsm.isin(Eall.columns) & U.key.isin(D.index)].drop_duplicates('key')
        mapfile, mapsha = H(F[mf]), sha256(F[mf])
    else:
        MF = os.path.expanduser('.cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv')
        M = pd.read_csv(MF, sep='\t', comment='#'); M = M[M.series == ds].copy()
        Dk = D.reset_index()[['key', 'sex', 'age', 'pH', 'PMI', 'dur', 'apd']].copy()
        for c in ['dur', 'apd']: Dk[c] = Dk[c].fillna(-1)
        Dk['pH'] = Dk.pH.round(3)
        Mm = M.rename(columns={'age_w': 'age', 'sex_w': 'sex', 'pH_w': 'pH', 'PMI_w': 'PMI', 'dur_w': 'dur', 'apd_w': 'apd'}).copy()
        for c in ['dur', 'apd']: Mm[c] = Mm[c].fillna(-1)
        Mm['pH'] = Mm.pH.round(3)
        J = Mm.merge(Dk, on=['sex', 'age', 'pH', 'PMI', 'dur', 'apd'], how='left')
        print(ds, 'matched-file rows', len(M), 'joined to one design/workbook donor', int(J.key.notna().sum()), 'duplicates', int(J.duplicated('gsm', keep=False).sum()))
        J = J[J.key.notna()].drop_duplicates('gsm', keep=False)
        U = J[J.gsm.isin(Eall.columns)].drop_duplicates('key')
        mapfile, mapsha = '.cache/stanley_demographics_local/matched_LOCAL_ONLY.tsv', sha256(MF)
    E, Z = prep(Eall, list(Eall.columns))
    C0 = D.loc[U.key].copy(); C0.index = U.gsm.values
    C0['dur10'] = C0.dur / 10
    ctl = list(C0.index[C0.dx == 'control']); casepm = list(C0.index[C0.dx.isin(['schizophrenia', 'BP'])])
    print(ds, 'arrays', Eall.shape[1], 'mapped (unique) donors', len(U), 'by dx', C0.dx.value_counts().to_dict(), 'genes', len(Z), 'directional measured', sum(x in Z.index for x in sg.index))
    for grp in ['SZ+BP', 'SZ only']:
        dxs = ['schizophrenia', 'BP'] if grp == 'SZ+BP' else ['schizophrenia']
        covs = ['age', 'sex', 'pH', 'PMI'] + (['bp'] if grp == 'SZ+BP' else [])
        for label, col, sub in VARS:
            C = C0[C0.dx.isin(dxs)].dropna(subset=['age', 'sex', 'pH', 'PMI']).copy()
            binary = sub is None
            if col == 'logapd': C = C.dropna(subset=['apd']); C['x'] = (C.logapd - C.logapd.mean()) / C.logapd.std(); unit = 'log10 dose, per SD'
            elif col == 'dur10': C = C.dropna(subset=['dur10']); C['x'] = C.dur10; unit = 'per 10 years'
            else: C = C.dropna(subset=[col]); C['x'] = C[col]; unit = 'yes minus no'
            n = len(C); ny = int(C.x.sum()) if binary else None
            rec = dict(dataset=ds, variable=label, group=grp, n_cases_value=n, n_yes=ny, n_no=(n - ny) if binary else None, unit=unit)
            kpar = len(covs) + 2
            if binary and (min(ny, n - ny) < 5 or n < kpar + 10):
                rec.update(tested='no', reason='smaller category %d (< 5) or n %d < %d (model terms %d + 10 residual df)' % (min(ny, n - ny), n, kpar + 10, kpar)); rows.append(rec); continue
            o, _ = run_test(Z, E, C, 'x', covs, ctl, casepm, ds + ' ' + label + ' ' + grp, rng, sg=sg)
            rec.update(tested='yes', **{k: v for k, v in o.items() if k not in ('label', 'n_cases_value')}); rows.append(rec)
    prov[ds] = dict(n_mapped=len(U), expr=H(F[xf]), expr_sha=sha256(F[xf]), map=mapfile, map_sha=mapsha, mapping=MAPSRC)
    inv += [dict(scope=ds + ' mapped donors', variable=v, diagnosis=dx, value='n', n=int(((C0.dx == dx) & C0[c].notna()).sum())) for v, c in [('rate of death', 'rod'), ('smoking', 'smoking'), ('suicide', 'suicide'), ('alcohol', 'alcohol'), ('drugs', 'drugs'), ('lifetime antipsychotics', 'apd'), ('duration', 'dur')] for dx in ['schizophrenia', 'BP', 'MDD', 'control']]
I = pd.DataFrame(inv); I.to_csv(OUT + '/stanley_all_inventory_%s.tsv' % MAPSRC, sep='\t', index=False)
R = pd.DataFrame(rows)
hdr = ['# script: stanley_all.py (core: dd_core.py run_test; model of c15_stanley_within_case.py)', '# run: 2026-10-19 (F_notes_10-12 factcheck resolution); sets file via env F_SETS, otherwise unchanged from dose_duration_tests_2026-10-14/scripts/stanley_all.py',
       '# input: %s sha256 %s (set study_brain; the file c15 used)' % (H(F['sets']), sha256(F['sets'])), '# input: %s sha256 %s (SMRI design file, 165 donors)' % (H(F['design']), sha256(F['design'])),
       '# input: ArrayCollection.xlsx sha256 %s ; NeuropathologyConsortium.xlsx sha256 %s (SMRI; read in place; duration of illness only)' % (sha256(F['arr']), sha256(F['neu']))] + \
      ['# input: %s' % json.dumps(v) + ' (' + k + ')' for k, v in prov.items()] + \
      ['# model: OLS HC3, score ~ variable + age + sex + pH + PMI (+ BP indicator in SZ+BP); est = change in score (mean z) per unit; est_in_control_sd = est / SD of control scores',
       '# matched null: 2000 random gene sets, matched gene by gene on expression decile x |SZ+BP-vs-control Welch t| decile, signed by own t; same slope by FWL; P = (1 + #null >= observed)/2001 (one-sided upper)',
       '# nominal share: per-gene OLS in the same cases and model; share of the recurrent genes at P < 0.05 against the all-gene share; one-sided Fisher',
       '# tested = no: smaller category < 5 or n < model terms + 10']
with open(OUT + '/stanley_all_tests_raw_%s.tsv' % MAPSRC, 'w') as f:
    f.write('\n'.join(hdr) + '\n'); R.to_csv(f, sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 280)
print(R[['dataset', 'variable', 'group', 'n_cases_value', 'n_yes', 'tested', 'est', 'lo_hc3', 'hi_hc3', 'p_hc3', 'matched_null_mean', 'matched_null_lo', 'matched_null_hi', 'matched_P_null_ge_obs', 'n_set_p05', 'expected_p05']].round(3).to_string())
