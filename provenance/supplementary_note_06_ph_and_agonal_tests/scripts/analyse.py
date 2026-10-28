#!/usr/bin/env python3
"""pH ~ cause of death with study fixed effects, controls first, then SZ vs control with and without COD (2026-10-09).
Input donors.tsv (extraction_notes.md). Studies with pH and COD: Arion 2017, Gilabert-Juan 2015, Ishiguro 2010.
slow = respiratory_infection, cancer, other_chronic; fast = cardiac_acute, suicide, accident_trauma, vascular_other;
ambiguous and undetermined excluded from the fast/slow contrast. Outputs results.tsv."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
d = pd.read_csv('donors.tsv', sep='\t', comment='#')
d['ph'] = pd.to_numeric(d.brain_pH, errors='coerce'); d['age'] = pd.to_numeric(d.age_years, errors='coerce')
d = d[d.ph.notna() & d.cod_class.notna() & (d.cod_class != '')]
fs = {'respiratory_infection': 'slow', 'cancer': 'slow', 'other_chronic': 'slow', 'cardiac_acute': 'fast', 'suicide': 'fast', 'accident_trauma': 'fast', 'vascular_other': 'fast'}
d['speed'] = d.cod_class.map(fs)
d['dx'] = d.diagnosis_class
rows = []
def rec(name, m, term, n):
    for t in [k for k in m.params.index if term in k]:
        ci = m.conf_int().loc[t]; rows.append({'model': name, 'term': t, 'n': n, 'est': round(m.params[t], 3), 'lo': round(ci[0], 3), 'hi': round(ci[1], 3), 'p': round(m.pvalues[t], 4)})
c = d[d.dx == 'control']
x = c.dropna(subset=['speed', 'age'])
m = smf.ols('ph ~ C(study) + C(speed, Treatment("fast")) + age', x).fit(); rec('controls: pH ~ study + slow-vs-fast + age', m, 'speed', len(x))
m = smf.ols('ph ~ C(study) + C(cod_class, Treatment("cardiac_acute")) + age', c.dropna(subset=['age'])).fit(); rec('controls: pH ~ study + cod_class + age', m, 'cod_class', len(c.dropna(subset=['age'])))
for s in ['Arion 2017', 'Gilabert-Juan 2015', 'Ishiguro 2010']:
    y = x[x.study == s]
    if y.speed.nunique() == 2 and len(y) > 6:
        m = smf.ols('ph ~ C(speed, Treatment("fast")) + age', y).fit(); rec(f'controls within {s}: slow-vs-fast + age', m, 'speed', len(y))
z = d[d.dx.isin(['SZ', 'control'])].dropna(subset=['age'])
m0 = smf.ols('ph ~ C(study) + C(dx, Treatment("control")) + age', z).fit(); rec('SZ vs control: study + age', m0, 'dx', len(z))
z2 = z.dropna(subset=['speed'])
m1 = smf.ols('ph ~ C(study) + C(dx, Treatment("control")) + age', z2).fit(); rec('SZ vs control, speed-known donors: study + age', m1, 'dx', len(z2))
m2 = smf.ols('ph ~ C(study) + C(dx, Treatment("control")) + C(speed, Treatment("fast")) + age', z2).fit(); rec('SZ vs control, speed-known donors: + slow-vs-fast', m2, 'dx', len(z2)); rec('  slow-vs-fast term', m2, 'speed', len(z2))
print('SZ share slow vs control share slow (speed-known):'); print(pd.crosstab(z2.dx, z2.speed, normalize='index').round(2)); print(pd.crosstab(z2.dx, z2.speed))
r = pd.DataFrame(rows); r.to_csv('results.tsv', sep='\t', index=False); print(r.to_string(index=False))
