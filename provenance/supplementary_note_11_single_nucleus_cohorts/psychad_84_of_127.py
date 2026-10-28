#!/usr/bin/env python3
"""Reproduce '84 of 127' Mount Sinai (MSSM) PsychAD SZ cases with a dementia term.
Rule (same as scripts/pb_common.R::load_donors): dementia = regex 'dementia|alzheimer' (case-insensitive)
on the CELLxGENE 'disease' label string(s) in donors.tsv; SZ = dx_gemma == 'SZ'; bank = source_bank M.
Input: donors.tsv (build_donors.py, Gemma 2 API 2026-09-25 + PsychAD Supp. Table 1).
Output: aggregate counts only (no donor rows). Usage: psychad_84_of_127.py donors.tsv out.tsv"""
import sys, re, hashlib, datetime, pandas as pd
src, out = sys.argv[1], sys.argv[2]
sha = hashlib.sha256(open(src, 'rb').read()).hexdigest()
d = pd.read_csv(src, sep='\t', dtype=str).fillna('')
rx = re.compile('dementia|alzheimer', re.I)
d['dem'] = d.disease.map(lambda s: bool(rx.search(s)))
d['bank'] = d.source_bank.map({'H': 'HBCC', 'M': 'MSSM', 'R': 'RADC'})
sz = d[d.dx_gemma == 'SZ']
rows = []
for b, g in sz.groupby('bank'):
    rows.append((b, 'SZ cases (dx_gemma SZ)', len(g)))
    rows.append((b, 'SZ with dementia/alzheimer term', int(g.dem.sum())))
    rows.append((b, 'SZ without term (kept)', int((~g.dem).sum())))
m = sz[(sz.bank == 'MSSM') & sz.dem]
terms = {}
for s in m.disease:
    for t in s.split(' || '):
        if rx.search(t): terms[t.strip().lower()] = terms.get(t.strip().lower(), 0) + 1
for t, n in sorted(terms.items()):
    rows.append(('MSSM', f'SZ donors carrying term "{t}" (donor can carry several)', n))
c = d[(d.dx_gemma == 'CTRL')]
for b, g in c.groupby('bank'):
    rows.append((b, 'CTRL (dx_gemma CTRL) with term', int(g.dem.sum())))
with open(out, 'w') as f:
    f.write(f"# script: psychad_84_of_127.py\n# input: {src} sha256 {sha}\n# date: {datetime.date.today()}\n")
    f.write("# rule: dementia = regex dementia|alzheimer, case-insensitive, on CELLxGENE disease field; SZ = dx_gemma SZ\n")
    f.write("# aggregate counts only; no donor-level rows\nbank\tquantity\tn\n")
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
print(open(out).read())
