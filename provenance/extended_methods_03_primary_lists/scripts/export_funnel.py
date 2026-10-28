#!/usr/bin/env python3
"""Literature funnel table for manuscript Figure 1a, on the 515-study corpus (62 primary brain studies).

Writes figure_data/funnel/funnel_table.tsv and funnel_notes.md (same columns as the 2026-10-09 table).
Rows A1 to T3 (publications examined) are carried over unchanged from the existing funnel_table.tsv (repo copy if present, else
figure_data/funnel/). Rows 1 to 11b are recomputed here: steps 1 and 2 to 8, 11 from figure_data/prisma_flow/prisma_flow.tsv
(boxes PF01 to PF09, written by analyses/prisma_flow_table/build_prisma_flow.py), the relation, group and tissue splits from corpus/lineage.json,
steps 9, 10, 11a, 11b from analyses/recurrence_units.tsv joined to lineage.json.
Run order: build_prisma_flow.py, this script, build_prisma_flow.py (its cross-check reads this table). No model calls, no network.
Run: python3 analyses/figure_data_exports/export_funnel.py"""
import csv, json, os, hashlib, collections as C
R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) + '/'
OUT = R + 'figure_data/funnel/'
OLD = OUT + 'funnel_table.tsv'
if not os.path.exists(OLD):
    OLD = os.path.expanduser('figure_data/funnel/funnel_table.tsv')
V = 'Research Assistant/schizophrenia-transcriptome-replication/'
csv.field_size_limit(10**9)


def sha(p):
    return hashlib.sha256(open(R + p, 'rb').read()).hexdigest()[:12]


def tsv(p):
    return list(csv.DictReader([l for l in open(R + p, newline='') if not l.startswith('#')], delimiter='\t', quoting=csv.QUOTE_NONE))


old_lines = open(OLD, newline='').read().split('\n')
old_hdr = [l for l in old_lines if l.startswith('#')]
old_rows = [r for r in csv.reader([l for l in old_lines if not l.startswith('#') and l], delimiter='\t')]
COLS = old_rows[0]
assert COLS == ['step', 'what', 'n', 'unit', 'source file and field', 'sha256 prefix (12)', 'note'], COLS
STREAM = [r for r in old_rows[1:] if r[0][0] in 'ABCDEFGT']
assert [r[0] for r in STREAM][0] == 'A1' and len(STREAM) == 19, len(STREAM)

LIN, PLJ, RUF, PFT, EXC = 'corpus/lineage.json', 'analyses/primary_lists.json', 'analyses/recurrence_units.tsv', 'figure_data/prisma_flow/prisma_flow.tsv', 'corpus/corpus_exclusions.json'
nodes = json.load(open(R + LIN))['nodes']
corp = [x for x in nodes if x.get('tech_tier') and str(x.get('category')) != '2']
excl = json.load(open(R + EXC)); EX = {k for k in excl if not k.startswith('_')}
N = len(corp)
cat2 = [x for x in nodes if str(x.get('category')) == '2']
pf = {r['box_id']: int(r['n']) for r in tsv(PFT)}
assert pf['PF01'] == N, (pf['PF01'], N)
nid = lambda x: x['pmid'] or 'doi_' + x['id'].split(':', 1)[-1].replace('/', '_')
assert EX <= {x['pmid'] for x in cat2}
now = 'recounted 2026-10-11 on the corpus after the %d moves to the targeted-panel / small-study set (%s); %d lineage category 2 nodes in all (%d one-gene studies)' % (len(EX), EXC, len(cat2), len(cat2) - len(EX))
rows = []


def add(step, what, n, unit, src, hashes, note=''):
    rows.append([step, what, str(n), unit, src, hashes, note])


rel = C.Counter(x['relation'] for x in corp)
RELS = [('a', 'primary'), ('b', 'single-gene follow-up'), ('c', 'data reuse'), ('d', 'meta-analysis / combined'), ('e', 'reanalysis of public data'), ('f', 'single-cell / single-nucleus')]
assert sum(rel[r] for _, r in RELS) == N
add('1', 'Corpus: studies in lineage.json (not the 69 outside-source nodes, not the Merikangas landmark, not lineage category 2)', N, 'publications (studies)',
    V + LIN + ' nodes with tech_tier and category != 2; figure_data/prisma_flow/prisma_flow.tsv PF01', sha(LIN),
    '%d nodes = %d + 69 outside-source + 1 Merikangas landmark + %d category 2; %s' % (len(nodes), N, len(cat2), now))
for k, r in RELS:
    add('1' + k, 'Corpus by lineage relation: ' + r, rel[r], 'studies', V + LIN + ' node.relation', sha(LIN), ('1a-1f sum to %d; ' % N if k == 'f' else '') + now)
grp = C.Counter({'1': 0, '1m': 0, '2': 0, '3': 0, 'p': 0, 'sc': 0})
for x in corp:
    grp[str(x['category']) if x.get('category') else 'sc'] += 1
assert sum(grp.values()) == N and set(grp) <= {'1', '1m', '2', '3', 'p', 'sc'}
GL = [('g', '1', '1 transcriptome-wide first look'), ('h', '1m', '1m targeted multigene panel (>=20 genes)'), ('i', '2', '2 one or few genes'), ('j', '3', '3 no new data (reuse, reanalysis, meta-analysis)'),
      ('k', 'p', 'p proteome-wide'), ('l', 'sc', 'sc single-cell / single-nucleus')]
for k, g, lab in GL:
    add('1' + k, 'Corpus by group (lineage node.category; sc = no category): ' + lab, grp[g], 'studies', V + LIN + ' node.category', sha(LIN), ('1g-1l sum to %d; ' % N if k == 'l' else '') + now)
add('1m', 'Corpus by kind: review', 0, 'studies', 'supplementary_methods 1.1 (reviews excluded); lineage.json relation systematic review = 1 (Merikangas 2022 landmark, outside the corpus)', sha(LIN),
    "5 studies take their lists from a review's compilation (provenance review_compilation); they are primary studies")
tis = C.Counter(x.get('tissue') for x in corp)
TL = [('n', 'brain', 'brain'), ('o', 'blood', 'blood'), ('p', 'iPSC', 'ipsc'), ('q', 'LCL', 'lcl'), ('r', 'other', 'other'), ('s', 'mixed', 'mixed')]
assert sum(tis.values()) == N
# the table keeps the letter ids of the 2026-10-09 table; 1m above is the review row, as before
for k, lab, key in TL:
    add('1' + k, 'Corpus by tissue: ' + lab, tis[key], 'studies', V + LIN + ' node.tissue', sha(LIN), ('1n-1s sum to %d; ' % N if k == 's' else '') + now)
PL = json.load(open(R + PLJ))['studies']
add('2', 'Corpus studies with a primary gene list', pf['PF02'], 'studies', V + PLJ + ' studies[*] with in_corpus and primary; PF02', sha(PLJ), now)
add('3', 'In-corpus studies with a primary DE first-look list (after dropping the studies whose primary lists are all non-DE and those whose only list is a reanalysis role)', pf['PF03'], 'studies',
    V + "analyses/recurrence.py load_units()['all']; PF03", sha('analyses/recurrence.py'), now)
add('4', 'After lineage exclusions (meta-lane, multi-series reanalyses, other reanalyses)', pf['PF04'], 'studies', V + "analyses/recurrence.py load_units()['studies']; PF04", sha('analyses/recurrence.py'),
    "the 2026-10-03 funnel note calls this step 'same-donor collapse'; " + now)
add('5', 'After same-donor collapse (same_data_as groups)', pf['PF05'], 'donor groups', V + "analyses/recurrence.py load_units() group_of over PF04; PF05", sha('analyses/recurrence.py'), now)
add('6', 'Published brain studies (scope of the population)', pf['PF06'], 'studies', V + 'analyses/primary_lists.json provenance published and tissue brain among PF04; PF06', sha(PLJ), now)
add('7', 'Published brain, after same-donor collapse', pf['PF07'], 'donor groups', V + 'analyses/recurrence.py load_units() group_of over PF06; PF07', sha('analyses/recurrence.py'), now)
units = [r for r in tsv(RUF) if r['scenario'] == 'published, brain']
mod = C.Counter(r['modelled'] for r in units)
assert len(units) == pf['PF08'] and mod['yes'] == pf['PF09']
add('8', "Published brain units carrying a row in the recurrence scenario 'published, brain'", pf['PF08'], 'units (donor groups)', V + RUF + " scenario 'published, brain'; PF08", sha(RUF),
    '%d donor groups at step 7 minus %d with no row' % (pf['PF07'], pf['PF07'] - pf['PF08']))
# steps 9 and 10: the class assignment is the 2026-10-03 note's, not a field in the data
NOTE10 = ['Vawter 2001', 'Vawter 2002', 'Arion 2007', 'Iwamoto 2005', 'Weidenhofer 2006', 'Bowden 2008', 'Dean 2007', 'Yanagi 2008', 'Puvogel 2022', 'Alsema 2025']
lab_of = {nid(x): x['label'] for x in nodes}
notmod = [r for r in units if r['modelled'] != 'yes']
in10 = [l for l in NOTE10 if l in {lab_of[r['representative']] for r in notmod}]
gone = [l for l in NOTE10 if l not in in10]
for l in in10:
    assert l in {lab_of[r['representative']] for r in notmod}
add('9', 'Not gene-level or otherwise not modelled, other than step 10 (miRNA, sncRNA, circRNA, lncRNA, transposable-element, HERV lists, and units without a usable list)', len(notmod) - len(in10), 'donor groups',
    V + RUF + ' rows with modelled = no (%s), less the step 10 rows' % ', '.join('%d %s' % (v, k) for k, v in mod.items() if k != 'yes'), sha(RUF),
    'class assignment is the 2026-10-03 note\'s, not a field in the data; the recurrence file records only modelled = no; labels: ' + ', '.join(sorted(lab_of[r['representative']] for r in notmod if lab_of[r['representative']] not in in10)))
add('10', 'Genuine gene-level units with tested-gene information missing or unusable', len(in10), 'donor groups',
    V + 'notes/Population funnel ... table (' + ', '.join(in10) + '); each is no in recurrence_units.tsv', sha(RUF),
    'the note lists 10; %s moved out of the corpus (%s)' % (', '.join(gone), EXC) if gone else '')
add('11', 'Population: published brain studies modelled (one list each)', pf['PF09'], 'studies (= donor groups)', V + RUF + " scenario 'published, brain', modelled = yes; PF09", sha(RUF),
    '%d - %d - %d = %d' % (pf['PF08'], len(notmod) - len(in10), len(in10), pf['PF09']))
tt = {nid(x): x.get('tech_tier') for x in nodes}
yes = [r for r in units if r['modelled'] == 'yes']
sn = sorted(lab_of[r['representative']] for r in yes if tt[r['representative']] == 'single_cell')
cl = [r for r in yes if not r['cluster'].startswith('left out')]
add('11a', 'Population: bulk studies', len(yes) - len(sn), 'studies', V + RUF + ' (modelled = yes) joined to lineage.json tech_tier != single_cell', sha(RUF) + '; ' + sha(LIN), 'verified by join')
add('11b', 'Population: single-nucleus / single-cell studies (%s)' % ', '.join(sn), len(sn), 'studies', V + RUF + ' (modelled = yes) joined to lineage.json tech_tier == single_cell', sha(RUF) + '; ' + sha(LIN),
    'verified by join; %d of %d sit in %d cohort clusters (%d multi-bank studies left out)' % (len(cl), len(yes), len({r['cluster'] for r in cl}), len(yes) - len(cl)))
# keep the id order of the 2026-10-09 table: 1, 1a..1l, 1m (review), 1n..1s, 2 ..
hdr = ['# rebuilt 2026-10-11 by analyses/figure_data_exports/export_funnel.py: rows A1 to T3 carried over from the table built 2026-10-09; rows 1 to 11b recomputed on the %d-study corpus (%d primary brain studies) from the files named in column 5; sha256 prefixes of those files as of 2026-10-11; paths relative to the vault root' % (N, pf['PF09']),
       '# rows 1 to 8 and 11 equal prisma_flow.tsv boxes PF01 to PF09 (cross-checked in prisma_flow_check.tsv)']
os.makedirs(OUT, exist_ok=True)
with open(OUT + 'funnel_table.tsv', 'w', newline='') as f:
    f.write('\n'.join(hdr) + '\n')
    w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')
    w.writerow(COLS)
    for r in STREAM:
        w.writerow(r)
    for r in rows:
        w.writerow([c.replace('\t', ' ').replace('\n', ' ') for c in r])
print('funnel rows', len(STREAM) + len(rows), 'corpus', N, 'steps', [(r[0], r[2]) for r in rows if r[0] in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '11a', '11b')])

# ---- funnel_notes.md: the 'Everything examined' paragraph is carried over; the corpus-to-population paragraphs are regenerated from the numbers above
OLDN = OUT + 'funnel_notes.md'
if not os.path.exists(OLDN):
    OLDN = os.path.expanduser('figure_data/funnel/funnel_notes.md')
old = open(OLDN, encoding='utf8').read().split('\n\n')
head = [p for p in old if p.startswith('**Everything examined')]
assert len(head) == 1
g = {r[0]: int(r[2]) for r in rows}
prov = C.Counter(PL[k]['provenance'] for k, v in PL.items() if v.get('in_corpus') and v.get('primary'))
text = ['# Literature funnel: wording per step', 'Rebuilt 2026-10-11 by analyses/figure_data_exports/export_funnel.py on the %d-study corpus, read-only, from the files and fields named in funnel_table.tsv (step ids below match its first column). Nothing was recounted by a model or fetched from the network. The first paragraph (A1 to T3) is carried over from the 2026-10-09 build.' % N,
        head[0],
        '**Corpus (step 1).** %d studies in corpus/lineage.json, after leaving out the 69 outside-source nodes (cited data sources that are not themselves corpus studies), the Merikangas review node, and %d lineage category 2 nodes (%d moved to the targeted-panel / small-study set per corpus/corpus_exclusions.json, %d one-gene studies that stay in the single-gene set). By lineage relation: %d primary, %d single-gene follow-up, %d data reuse, %d meta-analysis or combined, %d reanalysis of public data, %d single-cell or single-nucleus. By group: %d transcriptome-wide first look (1), %d no new data (3), %d proteome-wide (p), %d single-cell (sc). By tissue: %d brain, %d blood, %d iPSC, %d LCL, %d other, %d mixed. Reviews are excluded by rule.' % (
            N, len(cat2), len(EX), len(cat2) - len(EX), g['1a'], g['1b'], g['1c'], g['1d'], g['1e'], g['1f'], g['1g'], g['1j'], g['1k'], g['1l'], g['1n'], g['1o'], g['1p'], g['1q'], g['1r'], g['1s']),
        '**Studies with a primary list (step 2).** %d corpus studies have a primary gene list.' % g['2'],
        '**First-look DE units (steps 3 to 5).** %d studies have a first-look DE primary list. Removing the meta-lane studies, the multi-series reanalyses and the other reanalyses leaves %d studies (step 4); these sit in %d donor groups (same_data_as, step 5).' % (g['3'], g['4'], g['5']),
        '**Published brain scope (steps 6 to 8).** %d of the %d studies are published brain studies; they form %d donor groups, and %d carry a row in the "published, brain" scenario of the recurrence analysis.' % (g['6'], g['4'], g['7'], g['8']),
        '**Lacking information (steps 9 and 10).** %d of the %d units are not modelled: %d are gene-level units where the tested-gene information is missing or unusable (%s) and %d are other non-modelled units. The split between the two is the 2026-10-03 note\'s classification; the data hold only modelled yes/no with a reason class.' % (g['9'] + g['10'], g['8'], g['10'], ', '.join(in10), g['9']),
        '**Population (step 11).** %d published brain studies are modelled with one list each: %d bulk and %d single-cell or single-nucleus (%s), verified by joining recurrence_units.tsv to the lineage technology tier. %d of the %d fall in %d cohort clusters.' % (g['11'], g['11a'], g['11b'], ', '.join(sn), len(cl), len(yes), len({r['cluster'] for r in cl})),
        'Manuscript wording that the sources support: "%d corpus studies, %d with a primary list, %d first-look DE units, %d after removing meta-analyses and reanalyses (%d donor groups), %d published brain donor groups, of which %d lack usable gene-level information or are not modelled, leaving %d."' % (
            N, g['2'], g['3'], g['4'], g['5'], g['7'], g['7'] - g['11'], g['11']),
        '']
open(OUT + 'funnel_notes.md', 'w', encoding='utf8').write('\n\n'.join(text))
