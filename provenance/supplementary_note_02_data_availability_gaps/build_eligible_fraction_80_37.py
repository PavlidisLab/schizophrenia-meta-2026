#!/usr/bin/env python3
"""Eligible-fraction size column on the current Table 4 (80 studies) and Table 5 (37 rows) of the manuscript repository.

Eligible = the studies of Table 4 plus the Table 5 rows with eligibility_status E (confirmed) that are not already in Table 4;
rows marked P (provisional) are left out. Recorded subjects = schizophrenia + control subjects as the authors report them:
Table 4 studies from base81.tsv (n_SZ, n_control); Table 5 studies from study_level.tsv of eligible_fraction_2026-10-14/base71
(cases, controls), or from base81.tsv where the study is there (Shao 2008). A study with no recorded size counts 0 in the
recorded sum; the variant 'missing as 92' gives each such study 92 subjects, the value used for the lower bound of the earlier summary.
Writes study_level_80_37.tsv and summary_80_37.md next to this script. Reads only.
"""
import csv, os, hashlib
H = os.path.expanduser('~')
S = H + '/Documents/manuscripts/schizophrenia-transcriptome-reproducibility-manuscript/supplementary/'
A = H + '/<delivery folder>/'
OUT = os.path.dirname(os.path.abspath(__file__)) + '/'


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12]


def rd(p):
    return list(csv.DictReader((l for l in open(p, encoding='utf-8') if not l.startswith('#')), delimiter='\t'))


def num(x):
    try:
        return int(float(x))
    except Exception:
        return None


P4, P5 = S + 'supp_table3_data_included.tsv', S + 'supp_table4_data_needed.tsv'
PB, PL = A + 'base81_2026-10-14/base81.tsv', A + 'eligible_fraction_2026-10-14/base71/study_level.tsv'
t4, t5 = rd(P4), rd(P5)
b81 = {r['pmid']: r for r in rd(PB)}
sl = {r['pmid']: r for r in rd(PL)}
in4 = {r['pmid'] for r in t4}
rows = []
for r in t4:
    b = b81[r['pmid']]
    c, k = num(b['n_SZ']), num(b['n_control'])
    rows.append(dict(table='4', study=r['study'], pmid=r['pmid'], eligible='yes', cases=c, controls=k, size_source='base81.tsv'))
for r in t5:
    p = r['pmid']
    if p in in4:
        continue            # already counted through Table 4
    elig = 'yes' if r['eligibility_status'] == 'E' else 'no (provisional)'
    c = k = None
    src = ''
    if p in b81:
        c, k, src = num(b81[p]['n_SZ']), num(b81[p]['n_control']), 'base81.tsv'
    elif p in sl:
        c, k, src = num(sl[p]['cases']), num(sl[p]['controls']), 'study_level.tsv (base71)'
    if c is None and k is None:
        src = (src + '; ' if src else '') + 'no recorded size'
    rows.append(dict(table='5', study=r['study'], pmid=p, eligible=elig, cases=c, controls=k, size_source=src))
el = [r for r in rows if r['eligible'] == 'yes']
rec = lambda r: (r['cases'] or 0) + (r['controls'] or 0)
miss = [r for r in el if r['cases'] is None and r['controls'] is None]
t4n = sum(rec(r) for r in el if r['table'] == '4')
t5n = sum(rec(r) for r in el if r['table'] == '5')
with open(OUT + 'study_level_80_37.tsv', 'w') as f:
    f.write('# built by build_eligible_fraction_80_37.py; inputs %s sha256 %s; %s %s; %s %s; %s %s\n' % (
        'supp_table3_data_included.tsv', sha(P4), 'supp_table4_data_needed.tsv', sha(P5), 'base81.tsv', sha(PB), 'study_level.tsv', sha(PL)))
    w = csv.DictWriter(f, ['table', 'study', 'pmid', 'eligible', 'cases', 'controls', 'size_source'], delimiter='\t')
    w.writeheader()
    w.writerows(rows)
n_el = len(el)
n4 = len(t4)
sumall = t4n + t5n
md = f"""# Eligible fraction on the 80-study Table 4 and the 37-row Table 5 (2026-10-21)

Inputs: Table 4 sha {sha(P4)}, Table 5 sha {sha(P5)}, base81.tsv sha {sha(PB)}, study_level.tsv sha {sha(PL)}.

- Eligible studies: {n_el} = {n4} in Table 4 + {n_el - n4} confirmed (E) Table 5 rows not in Table 4. Table 5 has {len(t5)} rows; {sum(1 for r in t5 if r['eligibility_status'] != 'E')} are provisional (P) and {sum(1 for r in t5 if r['pmid'] in in4)} are already in Table 4.
- Studies: 80 of {n_el} = {80 / n_el:.1%}.
- Recorded subjects (schizophrenia + control, before quality control): Table 4 {t4n}; eligible Table 5 rows not in Table 4 {t5n}; total {sumall}. Share of the 80: {t4n} of {sumall} = {t4n / sumall:.1%}.
- {len(miss)} eligible studies have no recorded size and count 0: {', '.join(r['study'] for r in miss)}. Counting each as 92 subjects (the earlier lower bound) the total is {sumall + 92 * len(miss)} and the share {t4n / (sumall + 92 * len(miss)):.1%}.
- Shao 2008 is a Table 5 row that is not in the 80; its size comes from base81.tsv.
"""
open(OUT + 'summary_80_37.md', 'w').write(md)
print(md)
