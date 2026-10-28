"""Merge the adjudicated chunks into covariate_handling_64.tsv (one row per study x covariate, final call + quote),
agreement summary (reader 1 vs reader 2, before adjudication), and comparison with the earlier h2 first-pass calls."""
import csv, glob, collections, hashlib, os, datetime
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
H2 = os.path.join(REPO, 'analyses/h2/h2_studies.tsv')
def rd(f): return list(csv.DictReader((l for l in open(f) if not l.startswith('#')), delimiter='\t'))
def sha(f): return hashlib.sha256(open(f, 'rb').read()).hexdigest()[:12]
files = sorted(glob.glob(os.path.join(HERE, 'adjudicated/chunk[0-9].tsv')))
import json
EXCL = {k for k in json.load(open(os.path.join(REPO, 'corpus/corpus_exclusions.json'))) if not k.startswith('_')}  # left the corpus (2026-10-10)
rows = [r for f in files for r in rd(f) if r['pmid'] not in EXCL]
by = {(r['pmid'], r['covariate']): r for r in rows}
pm = sorted({r['pmid'] for r in rows})
covs = ['age','sex','pmi','rin','medication','smoking','substance_use','cause_of_death_agonal','batch','cell_composition','ancestry','hidden_factors']
out = os.path.join(HERE, 'covariate_handling_64.tsv')
cols = ['pmid','label','covariate','reader1_call','reader2_call','agreed','final_call','quote','source_file','page_or_line','quote_verified','resolution_note']
with open(out, 'w') as fh:
    fh.write(f"# script: analyses/covariate_handling_by_study/finalize.py, {datetime.date.today()}\n")
    for f in files: fh.write(f"# input: {os.path.relpath(f, REPO)} sha256 {sha(f)}\n")
    fh.write("# population: the primary brain studies (recurrence.py scenario 'published, brain'; studies in corpus/corpus_exclusions.json left out); 12 covariates; pH is in analyses/ph_agonal_handling_by_study\n")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n'); w.writeheader(); w.writerows(rows)
n = len(rows); ag = sum(r['agreed'] == 'yes' for r in rows)
print('studies', len(pm), 'calls', n, 'expected', len(pm) * 12, 'reader agreement', ag, f'{ag/n:.3f}')
print('final calls', dict(collections.Counter(r['final_call'] for r in rows)))
pc = collections.defaultdict(lambda: [0, 0])
for r in rows:
    pc[r['covariate']][0] += 1; pc[r['covariate']][1] += r['agreed'] == 'yes'
for k in covs: print(k, pc[k])
# versus h2 first pass
h2 = {r['pmid']: r for r in rd(H2)}
diff = collections.Counter(); tot = collections.Counter(); dl = []
for r in rows:
    h = h2.get(r['pmid'])
    if not h: continue
    old = h.get('cov_' + r['covariate'], '') or 'blank'
    tot[r['covariate']] += 1
    if old != r['final_call']:
        diff[r['covariate']] += 1; dl.append((r['pmid'], r['label'], r['covariate'], old, r['final_call']))
print('differs from h2 first pass:', sum(diff.values()), 'of', sum(tot.values()))
for k in covs: print(k, diff[k], '/', tot[k])
with open(os.path.join(HERE, 'vs_h2_first_pass.tsv'), 'w') as fh:
    fh.write('pmid\tlabel\tcovariate\th2_call\tfinal_call\n')
    for d in dl: fh.write('\t'.join(d) + '\n')
