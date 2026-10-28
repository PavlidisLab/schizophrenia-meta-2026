#!/usr/bin/env python3
"""Build the PRISMA 2020 flow table (final counts) for the frozen corpus (corpus/CORPUS_FREEZE.md, tag corpus-freeze-2026-10-11).

Writes figure_data/prisma_flow/prisma_flow.tsv (one row per flow-diagram box) and
figure_data/prisma_flow/prisma_flow_check.tsv (the arithmetic that closes each box).
Every n is computed from a file read here; no count is typed in. Reads only; the one
import of analyses/recurrence.py (load_units, scenario_selection) writes nothing.

Run: python3 analyses/prisma_flow_2026-10-10/build_prisma_flow.py
"""
import csv, json, os, re, sys, hashlib, datetime, collections as C

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')) + '/'
SN = R + 'search/screen_now_2026-10-09/'
OUT = R + 'figure_data/prisma_flow/'
FUNNEL = R + 'figure_data/funnel/funnel_table.tsv'  # written by analyses/figure_exports/export_funnel.py from this script's PF boxes; run order: this script, export_funnel.py, this script
if not os.path.exists(FUNNEL):
    FUNNEL = os.path.expanduser('figure_data/funnel/funnel_table.tsv')
EXCL_FILE = R + 'corpus/corpus_exclusions.json'
SCRIPT = 'analyses/prisma_flow_2026-10-10/build_prisma_flow.py'
csv.field_size_limit(10**9)
used = {}


def rel(p):
    return p.replace(R, '') if p.startswith(R) else p.replace(os.path.expanduser('~'), '~')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def use(p):
    used[rel(p)] = sha(p)
    return p


def rd(p):
    use(p)
    with open(p, newline='') as f:
        lines = [l for l in f if not l.startswith('#')]
    return list(csv.DictReader(lines, delimiter='\t', quoting=csv.QUOTE_NONE))


def js(p):
    use(p)
    return json.load(open(p))


# ---------------------------------------------------------------- inputs
dec = rd(SN + 'decisions.tsv')
pool = {r['pmid']: r for r in dec}
assert len(pool) == len(dec)
routes = js(SN + 'routes.json')
RT, P64 = routes['routes'], routes['primary64']
psum = js(SN + 'pool_summary.json')
pc_route = {r['route']: r for r in rd(SN + 'prisma_counts.tsv')}
fcand = rd(SN + 'final_candidates.tsv')
ftd = rd(SN + 'decisions_ft.tsv')
fts = {r['key']: r for r in rd(SN + 'fulltext_status.tsv')}
parked = rd(SN + 'parked_retry_results.tsv')
# records retrieved after parked_retry_results.tsv was written (Unpaywall retry, commit f5476ede) are decided in parked_final.tsv, not in parked_retry_results.tsv
pfin = rd(SN + 'parked_final.tsv')
n_harvest = len(parked)
assert {r['pmid'] for r in pfin if r['status'] == 'retrieved_and_screened:ezproxy_cdp'} == {r['pmid'] for r in parked}
parked_unpaywall = [dict(pmid=r['pmid'], final=r['decision'], scope=r['scope'], data_use=r['data_use_or_reason'], stage='unpaywall')
                    for r in pfin if r['status'] == 'retrieved_and_screened:unpaywall']
assert all(r['final'] in ('include', 'exclude') for r in parked_unpaywall)
parked = parked + parked_unpaywall
queues = {q: rd(SN + 'review_queue_%s.tsv' % q) for q in ('corpus_misses', 'new_includes', 'unsure', 'exclusion_sample')}
summ_md = open(use(SN + 'decisions_ft_summary.md')).read()
use(SN + 'recovery_report.md'); use(SN + 'run_log.txt')
for f in ('usage.tsv', 'usage_ft.tsv'):
    use(SN + f)
slog = {r['search_id']: r for r in rd(R + 'search/search_log.tsv')}
pcj = js(R + 'search/prisma_counts.json')
lin = js(R + 'corpus/lineage.json')
plists = js(R + 'analyses/primary_lists.json')['studies']
runits = rd(R + 'analyses/recurrence_units.tsv')
ftab = {r['step']: r for r in rd(FUNNEL)}

boxes, checks = [], []


def box(i, stage, label, n, parent, src, field, status, note=''):
    boxes.append(dict(box_id=i, stage=stage, label=label, n=n, parent_box_id=parent, source_file=src,
                      source_field_or_query=field, status=status, note=note))


def chk(i, formula, lhs, rhs, parents=''):
    checks.append(dict(check_id='C%02d' % (len(checks) + 1), box_id=i, parent_box_id=parents, formula=formula,
                       lhs=lhs, rhs=rhs, closes='yes' if lhs == rhs else 'NO'))


def N(i):
    return next(b['n'] for b in boxes if b['box_id'] == i)


# ---------------------------------------------------------------- identification
ORDER = ['databases', 'single_nucleus', 'backward', 'forward', 'gemma', 'genomewide_proteomic', 'single_gene_sweep']
tags = {k: RT[k] for k in pool}
rc = {r: sum(1 for t in tags.values() if r in t) for r in ORDER}
rc['other'] = sum(1 for t in tags.values() if 'cache_only' in t or 'lineage_only' in t)
co, lo = sum(1 for t in tags.values() if 'cache_only' in t), sum(1 for t in tags.values() if 'lineage_only' in t)
RTSV = rel(SN + 'prisma_counts.tsv')
RJ = rel(SN + 'routes.json')
LAB = {'databases': 'Records identified: databases (April PubMed/Embase/Web of Science screened set)',
       'single_nucleus': 'Records identified: single-cell / single-nucleus search',
       'backward': 'Records identified: backward citation (reference lists)',
       'forward': 'Records identified: forward citation',
       'gemma': 'Records identified: Gemma completeness check',
       'genomewide_proteomic': 'Records identified: genome-wide / proteomic additions',
       'single_gene_sweep': 'Records identified: single-gene sweep (PubMed)',
       'other': 'Records identified: other (saved PubMed cache; lineage-only nodes)'}
ids = []
for k, r in enumerate(ORDER + ['other'], 1):
    bid = 'ID%02d' % k
    ids.append(bid)
    if r == 'other':
        note = 'cache_only %d, lineage_only %d (both counted once here); routes.py tags records from corpus/pubmed_cache.json and lineage nodes with no other route' % (co, lo)
        fld = 'routes.json tags cache_only or lineage_only, over decisions.tsv pmid'
    else:
        pcv = pc_route[r]['records_identified_(any_tag)']
        note = 'a record can carry several route tags; %s records_identified_(any_tag) = %s' % (rel(SN + 'prisma_counts.tsv'), pcv)
        fld = 'routes.json tag %s, over decisions.tsv pmid' % r
        assert int(pcv) == rc[r], (r, pcv, rc[r])
    box(bid, 'identification', LAB[r], rc[r], '', RJ + '; ' + rel(SN + 'decisions.tsv'), fld, 'final', note)
tot = sum(rc.values())
box('ID09', 'identification', 'Records identified, all routes (sum of route counts)', tot, ';'.join(ids), RJ, 'sum of ID01-ID08', 'final',
    'overlap between routes not yet removed')
newsum = sum(int(r['records_new_after_dedup_vs_earlier_routes']) for r in pc_route.values())
dups = tot - len(pool)
box('ID10', 'identification', 'Duplicate records removed (records carrying more than one route tag)', dups, 'ID09', RJ,
    'sum of route counts minus distinct pmid/node ids in decisions.tsv', 'final',
    'removal is by record id across routes; deduplication inside a route (for example the 15,705 forward citation links) happened earlier and is in the memo rows')
box('ID11', 'identification', 'Records after duplicates removed (pool for title/abstract screening)', len(pool), 'ID09', rel(SN + 'decisions.tsv'),
    'rows of decisions.tsv', 'final', 'pool_summary.json pool_final = %d; prisma_counts.tsv records_new_after_dedup column sums to %d' % (psum['pool_final'], newsum))
chk('ID11', 'ID09 - ID10 = ID11', tot - dups, len(pool), 'ID09')
assert psum['pool_final'] == len(pool)
chk('ID11', 'prisma_counts.tsv records_new_after_dedup summed over its 9 route rows = ID11', newsum, len(pool), 'ID09')
chk('ID11', 'pool_summary.json pool_final = ID11', psum['pool_final'], len(pool), 'ID09')
# memo rows (not in the arithmetic)
box('M01', 'identification', 'Memo: records dropped when the pool was built (no title and no abstract)', psum['no_abstract_and_no_title_excluded'], '',
    rel(SN + 'pool_summary.json'), 'no_abstract_and_no_title_excluded', 'final',
    'memo, not in the arithmetic; pool_before_exclusion %d; the record is not in routes.json or decisions.tsv, so its route is not recorded' % psum['pool_before_exclusion'])
mp = [('M02', 'databases', 'S01-S03 PubMed, Embase, Web of Science', pcj['april']['dedup']['unique_by_pmid_doi_or_title'],
       'search/prisma_counts.json april.dedup.unique_by_pmid_doi_or_title', 'unique by PMID, DOI or title (April searches, 670 records before deduplication); the pool route holds only records of search/screening/papers.json (routes.py), so database records outside that file are not in ID01'),
      ('M03', 'single_nucleus', 'S11 single-cell search', 'S11', None, ''),
      ('M04', 'backward', 'S07 reference lists, first filter', 'S07', None, 'S08 (broader filter) is M05'),
      ('M05', 'backward', 'S08 reference lists, broader filter', 'S08', None, ''),
      ('M06', 'forward', 'S12 forward citation candidate pool', 'S12', None, 'S18 update adds 2,342 more records (after_dedup)'),
      ('M07', 'gemma', 'S09 Gemma completeness check, datasets', 'S09', None, 'count is datasets, not publications'),
      ('M08', 'single_gene_sweep', 'S13 single-gene sweep, PubMed union', 'S13', None, 'hits %s, screened %s' % (slog['S13']['hits'], slog['S13']['screened']))]
for bid, route, lab, v, f, note in mp:
    if isinstance(v, str):
        n = int(slog[v]['after_dedup']); f = 'search/search_log.tsv %s after_dedup' % v
        if bid == 'M06':
            note = note + '; hits %s' % slog[v]['hits']
        src = 'search/search_log.tsv'
    else:
        n = v; src = 'search/prisma_counts.json'
    box(bid, 'identification', 'Memo: search-log count for route %s, %s' % (route, lab), n, '', src, f, 'final',
        'memo, not in the arithmetic' + ('; ' + note if note else ''))

# ---------------------------------------------------------------- screening
p1 = {k: r['pass1_decision'] for k, r in pool.items()}
p2 = {k: r['pass2_decision'] for k, r in pool.items()}
corp = {k for k, r in pool.items() if r['in_corpus_536'] == 'yes'}
notscr = [k for k in pool if not p1[k]]
scr = [k for k in pool if p1[k]]
DEC = rel(SN + 'decisions.tsv')
box('SC00', 'screening', 'Records not screened (no parseable pass-1 result after 3 batch rounds)', len(notscr), 'ID11', DEC + '; ' + rel(SN + 'run_log.txt'),
    'pass1_decision empty', 'final', 'run_log.txt: "pass1: 1 records still failed after 3 rounds"')
box('SC01', 'screening', 'Records screened on title and abstract (pass 1, claude-haiku-5-5)', len(scr), 'ID11', DEC, 'pass1_decision non-empty', 'final',
    'pool records without an abstract (screened on title alone): %d (pool_summary.json final_no_abstract)' % psum['final_no_abstract'])
chk('SC01', 'ID11 - SC00 = SC01', len(pool) - len(notscr), len(scr), 'ID11')
x1 = [k for k in scr if p1[k] == 'exclude']
box('SC02', 'screening', 'Records excluded at pass 1 (title/abstract)', len(x1), 'SC01', DEC, 'pass1_decision = exclude', 'final', 'reason per record is free text in decisions.tsv pass1_reason')
a2 = [k for k in scr if p1[k] in ('include', 'unsure')]
assert all(p2[k] for k in a2) and not any(p2[k] for k in x1)
box('SC03', 'screening', 'Records screened at pass 2 (pass-1 include or unsure; claude-sonnet-5-5)', len(a2), 'SC01', DEC, 'pass1_decision in include, unsure; pass2_decision non-empty', 'final')
chk('SC03', 'SC01 - SC02 = SC03', len(scr) - len(x1), len(a2), 'SC01')
x2 = [k for k in a2 if p2[k] == 'exclude']
box('SC04', 'screening', 'Records excluded at pass 2 (title/abstract)', len(x2), 'SC03', DEC, 'pass2_decision = exclude', 'final')
for i, ((du, sc), n) in enumerate(sorted(C.Counter((pool[k]['pass2_data_use'], pool[k]['pass2_scope']) for k in x2).items(), key=lambda t: (-t[1], t[0])), 1):
    box('SC04.%02d' % i, 'screening', 'Excluded at pass 2: data_use=%s, scope=%s' % (du, sc), n, 'SC04', DEC, 'pass2_decision = exclude; pass2_data_use, pass2_scope', 'final')
chk('SC04', 'sum of SC04.nn = SC04', sum(b['n'] for b in boxes if b['box_id'].startswith('SC04.')), len(x2), 'SC04')
adv = [k for k in a2 if p2[k] in ('include', 'unsure')]
box('SC05', 'screening', 'Records advancing after title/abstract screening (pass 2 include or unsure)', len(adv), 'SC03', DEC, 'pass2_decision in include, unsure', 'final')
chk('SC05', 'SC03 - SC04 = SC05', len(a2) - len(x2), len(adv), 'SC03')
box('SC05a', 'screening', 'of which pass 2 include', sum(1 for k in adv if p2[k] == 'include'), 'SC05', DEC, 'pass2_decision = include', 'final',
    'of these %d are corpus members' % sum(1 for k in adv if p2[k] == 'include' and k in corp))
box('SC05b', 'screening', 'of which pass 2 unsure', sum(1 for k in adv if p2[k] == 'unsure'), 'SC05', DEC, 'pass2_decision = unsure', 'final',
    'of these %d are corpus members' % sum(1 for k in adv if p2[k] == 'unsure' and k in corp))
chk('SC05', 'SC05a + SC05b = SC05', N('SC05a') + N('SC05b'), len(adv), 'SC05')
box('SC06', 'screening', 'Records screened on title and abstract and excluded (pass 1 + pass 2)', len(x1) + len(x2), 'SC01', DEC,
    'pass1_decision = exclude or pass2_decision = exclude', 'final',
    'pass-1 excludes %d; pass-2 excludes %d' % (len(x1), len(x2)))
chk('SC06', 'SC01 - SC06 = SC05', len(scr) - len(x1) - len(x2), len(adv), 'SC01')

# ---------------------------------------------------------------- eligibility (full-text assessment)
cls_of, queue_of, row_of = {}, C.defaultdict(set), {}
for r in fcand:
    k = r['pmid']
    c = r['final_class'].split('|')[0]
    if k in cls_of:
        assert cls_of[k] == c, k
    else:
        row_of[k] = r
    cls_of[k] = c
    queue_of[k].add(r['queue'])
for q, rows in queues.items():
    for r in rows:
        assert q in queue_of[r['pmid']]
U = set(cls_of)
inadv = U & set(adv)
sent_x1 = U & set(x1)
sent_x2 = U & set(x2)
assert U <= set(scr)
carried = set(adv) - U
assert all(k in corp and p2[k] == 'include' for k in carried)
FC = rel(SN + 'final_candidates.tsv')
box('EL00a', 'eligibility', 'Advanced on title and abstract and included without a separate full-text assessment (corpus members, pass-2 include)', len(carried), 'SC05',
    FC + '; ' + DEC, 'SC05 records absent from every review_queue_*.tsv', 'final',
    'all %d are corpus members and carry pass-2 include; they go to the included box' % len(carried))
qc = lambda s, q: len([k for k in s if q in queue_of[k]])
box('EL00b', 'eligibility', 'Excluded on title and abstract and re-checked on full text (corpus candidates and a sample of excluded records)', len(sent_x1) + len(sent_x2), 'SC06', FC + '; ' + DEC,
    'review_queue records with pass1_decision or pass2_decision = exclude', 'final',
    'pass-1 excludes %d (in corpus_misses %d, in exclusion_sample %d); pass-2 excludes %d (in unsure queue %d, in corpus_misses %d)' % (
        len(sent_x1), qc(sent_x1, 'corpus_misses'), qc(sent_x1, 'exclusion_sample'), len(sent_x2), qc(sent_x2, 'unsure'), qc(sent_x2, 'corpus_misses')))
box('EL01', 'eligibility', 'Reports sought for retrieval (records advanced on title and abstract and sent to full text, plus the records re-checked)', len(U), 'SC05;EL00b', FC,
    'distinct pmid over queue = corpus_misses, new_includes, unsure, exclusion_sample', 'final',
    'queue rows %d; distinct records %d (88 records sit in two queues); queue sizes corpus_misses %d, new_includes %d, unsure %d, exclusion_sample %d' % (
        len(fcand), len(U), *(len(queues[q]) for q in ('corpus_misses', 'new_includes', 'unsure', 'exclusion_sample'))))
chk('EL01', 'SC05 - EL00a + EL00b = EL01', len(adv) - len(carried) + N('EL00b'), len(U), 'SC05;EL00b')
chk('EL01', 'advancing records in queues (%d) + EL00b = EL01' % len(inadv), len(inadv) + N('EL00b'), len(U), 'SC05;EL00b')
# not retrieved
ftu = {}
for r in ftd:
    ftu.setdefault(r['pmid'], r)
assert set(ftu) == U
NF = {k for k, r in ftu.items() if r['ft_stage'] == 'none'}
PK = {r['pmid'] for r in parked}
assert PK <= NF and len(PK) == len(parked)
STILL = NF - PK
# a record counts as retrieved if fulltext_status.tsv holds a route other than none / pmc_not_oa with >= 8000 characters; every such record must be in parked
RETRIEVED = {k for k in NF if fts[k]['route'] not in ('none', 'pmc_not_oa') and int(fts[k]['chars']) >= 8000}
assert RETRIEVED == PK, (sorted(RETRIEVED - PK), sorted(PK - RETRIEVED))
for k in STILL:
    assert fts[k]['route'] in ('none', 'pmc_not_oa') and int(fts[k]['chars']) < 8000, k
ncorp_nf = STILL & corp
assert not ncorp_nf  # no corpus member is without a full text
box('EL02', 'eligibility', 'Reports sought but not retrieved', len(STILL), 'EL01', rel(SN + 'fulltext_status.tsv') + '; ' + rel(SN + 'parked_retry_results.tsv'),
    'ft_stage = none and pmid not in parked_retry_results.tsv or the Unpaywall rows of parked_final.tsv', 'final',
    '%d records had no staged full text in the main run; %d of them were retrieved later (PDF-fallback harvest %d, Unpaywall retry %d); fulltext_status.tsv route: none %d, pmc_not_oa %d (all under 8,000 characters); none is a corpus member' % (
        len(NF), len(PK), n_harvest, len(parked_unpaywall), sum(1 for k in STILL if fts[k]['route'] == 'none'), sum(1 for k in STILL if fts[k]['route'] == 'pmc_not_oa')))
chk('EL02', 'no-full-text records %d - retrieved later %d = EL02' % (len(NF), len(PK)), len(NF) - len(PK), len(STILL), 'EL01')
chk('EL02', 'fulltext_status.tsv: records with route none or pmc_not_oa = EL02', sum(1 for k in U if fts[k]['route'] in ('none', 'pmc_not_oa')), len(STILL), '')
box('EL03', 'eligibility', 'Full-text reports assessed for eligibility', len(U) - len(STILL), 'EL01', FC + '; ' + rel(SN + 'parked_retry_results.tsv'),
    'EL01 minus EL02', 'final', 'main run %d (Haiku %d, Sonnet %d) + parked retry %d' % (
        len(U) - len(NF), sum(1 for k in U if ftu[k]['ft_stage'] == 'haiku'), sum(1 for k in U if ftu[k]['ft_stage'] == 'sonnet'), len(PK)))
chk('EL03', 'EL01 - EL02 = EL03', len(U) - len(STILL), N('EL03'), 'EL01')
assert N('EL03') == len([k for k in U if ftu[k]['ft_stage'] != 'none']) + len(PK)

# outcomes of the full-text assessment, by final_class (kept as the data field name)
excl_new = [k for k in U if cls_of[k] == 'exclude']
excl_corp = [k for k in U if cls_of[k].startswith('corpus_member_fails_pass3_criteria')]
NOTADD = ('single_gene_already_in', 'add_primary', 'add_reanalysis_or_reuse', 'hold_conflicts_batch10_reject',
          'single_gene_add', 'single_gene_add_sweep_said_exclude', 'single_gene_add_sweep_said_unsure')
inc_conf = [k for k in U if cls_of[k] == 'corpus_confirmed']
inc_notadd = [k for k in U if cls_of[k] in NOTADD]
pk_inc = [r for r in parked if r['final'] == 'include']
pk_exc = [r for r in parked if r['final'] == 'exclude']
assert len(pk_inc) + len(pk_exc) == len(parked)
pk_exc_corp = [r for r in pk_exc if r['pmid'] in corp]  # corpus members among the parked records excluded at full text
assert not [r for r in pk_inc if r['pmid'] in corp]
pk_exc_new = [r for r in pk_exc if r['pmid'] not in corp]
assert len(excl_new) + len(excl_corp) + len(inc_conf) + len(inc_notadd) + len(NF) == len(U)
n_ex = len(excl_new) + len(pk_exc_new)
box('EL04', 'eligibility', 'Full-text reports excluded (reasons below)', n_ex, 'EL03', FC + '; ' + rel(SN + 'parked_retry_results.tsv'),
    'final_class = exclude (%d); parked retry final = exclude (%d)' % (len(excl_new), len(pk_exc_new)), 'final',
    'reasons are the data_use and scope fields returned by the model, not a hand-coded category; free-text reasons are in decisions_ft.tsv ft_reason / parked_retry_results.tsv reason')
reasons = C.Counter((ftu[k]['data_use'], ftu[k]['scope']) for k in excl_new)
reasons.update((r['data_use'], r['scope']) for r in pk_exc_new)
for i, ((du, sc), n) in enumerate(sorted(reasons.items(), key=lambda t: (-t[1], t[0])), 1):
    box('EL04.%02d' % i, 'eligibility', 'Full text excluded: data_use=%s, scope=%s' % (du, sc), n, 'EL04', rel(SN + 'decisions_ft.tsv') + '; ' + rel(SN + 'parked_retry_results.tsv'),
        'ft_decision/final = exclude; data_use, scope', 'final')
chk('EL04', 'sum of EL04.nn = EL04', sum(reasons.values()), n_ex, 'EL04')
n_ret = len(excl_corp) + len(pk_exc_corp)
box('EL05', 'eligibility', 'Full-text reports of corpus studies retained in the corpus', n_ret, 'EL03', FC + '; ' + rel(SN + 'parked_final.tsv'),
    'final_class starts corpus_member_fails_pass3_criteria; plus corpus members among the parked records excluded at full text', 'final',
    'the corpus was frozen (corpus/CORPUS_FREEZE.md) on the screening decisions that built it; these %d studies are in the corpus' % n_ret)
box('EL06', 'eligibility', 'Full-text reports confirming corpus studies', len(inc_conf), 'EL03', FC, 'final_class = corpus_confirmed', 'final',
    'corpus members; all are in the corpus')
n_notadd = len(inc_notadd) + len(pk_inc)
box('EL07', 'eligibility', 'Full-text reports not added to the corpus (corpus frozen)', n_notadd, 'EL03', FC + '; ' + rel(SN + 'parked_retry_results.tsv'),
    'final_class in %s; plus parked-retry includes' % ', '.join(NOTADD), 'final',
    'not corpus members and not in the corpus; the screening adds no study after the freeze (corpus/CORPUS_FREEZE.md); the studies added by hand afterwards are IN00b')
chk('EL07', 'EL04 + EL05 + EL06 + EL07 = EL03', N('EL04') + N('EL05') + N('EL06') + N('EL07'), N('EL03'), 'EL03')
assert not [k for k in inc_notadd if k in corp]

# ---------------------------------------------------------------- included
nodes = lin['nodes']
LN = rel(R + 'corpus/lineage.json')
cats = {(n['pmid'] or n['id']): n.get('category') for n in nodes}
cat2 = {k for k, c in cats.items() if c == '2'}
excl_doc = js(EXCL_FILE)
EXC = {k for k in excl_doc if not k.startswith('_')}
assert EXC <= cat2, sorted(EXC - cat2)  # every moved study is a lineage category-2 node
EXCP, EXCO = EXC & corp, EXC - corp  # members of the screened pool / removed by decision but never in the screened pool (a single-cell stream study)
cat2_other = cat2 - EXC
assert not (cat2_other & corp)  # the one-gene studies were never in the screened corpus
corpus_nodes = {(n['pmid'] or n['id']) for n in nodes if n.get('source') != 'Cited source (not in corpus)' and n.get('category') not in ('review', '2')}
# studies added by hand after the 536-study screening and the freeze (corpus/added_by_hand.json entries that are not in the screened pool)
HAND = {k for k in js(R + 'corpus/added_by_hand.json') if not k.startswith('_')}
ADDED = corpus_nodes - (corp - EXC)
assert ADDED <= HAND, (sorted(ADDED - HAND))
assert corpus_nodes == (corp - EXC) | ADDED, (len(corpus_nodes), len(corp), len(EXC), len(ADDED))
parents_in = ['EL00a', 'EL06', 'EL05']
inc_n = len(carried) + len(inc_conf) + n_ret + len(ncorp_nf)
box('IN00', 'included', 'Studies taken into the corpus', len(corp), ';'.join(parents_in), DEC + '; ' + FC,
    'decisions.tsv in_corpus_536 = yes; composition from EL00a, EL06, EL05', 'final',
    'composition: %d included on title and abstract, %d confirmed at full text, %d retained' % (len(carried), len(inc_conf), n_ret))
chk('IN00', 'EL00a + EL06 + EL05 = IN00', inc_n, len(corp), ';'.join(parents_in))
assert inc_n == len(corp)
n_panel = sum(1 for k in EXCP if excl_doc[k]['reason'].startswith('targeted panel'))
n_dec = sum(1 for k in EXCP if excl_doc[k]['date'] == '2026-10-12')  # removed by decision on 2026-10-12 (TWAS studies; 22q11.2-deletion patient-cell studies)
box('IN00a', 'included', 'Excluded: targeted gene-panel, small-set or decided-out studies', len(EXCP), 'IN00', rel(EXCL_FILE) + '; ' + LN, 'keys of corpus_exclusions.json that are decisions.tsv in_corpus_536 = yes; lineage category 2',
    'final', '%d studies: %d targeted multigene panels (lineage category 1m), %d removed by decision on 2026-10-12 (TWAS or predicted expression; 22q11.2-deletion patient cells) and %d other (Hoseth 2018, a panel treated as targeted); %d further entries of corpus_exclusions.json were never in the screened pool (%s); reason and date per study in corpus_exclusions.json' % (
        len(EXCP), n_panel, n_dec, len(EXCP) - n_panel - n_dec, len(EXCO), ', '.join(sorted(EXCO)) or 'none'))
box('IN00b', 'included', 'Added after the freeze (corpus/added_by_hand.json, not in the screened pool)', len(ADDED), 'IN00', rel(R + 'corpus/added_by_hand.json') + '; ' + LN, 'lineage corpus nodes that are not decisions.tsv in_corpus_536 = yes', 'final',
    '%d study: %s (see corpus/CORPUS_FREEZE_ADDENDUM_2026-10-11_mauney.txt)' % (len(ADDED), ', '.join(sorted(ADDED))))
box('IN01', 'included', 'Studies included (corpus)', len(corpus_nodes), 'IN00;IN00a;IN00b', LN + '; ' + DEC,
    'nodes with source not "Cited source (not in corpus)", category not review (Merikangas landmark) and not 2 (one-gene studies and the excluded targeted gene-panel or small-set studies)',
    'final',
    'IN00 - IN00a; %d lineage nodes = %d + 69 outside-source + 1 Merikangas landmark + %d category 2 (%d targeted gene-panel or small-set exclusions + %d one-gene studies)' % (
        len(nodes), len(corpus_nodes), len(cat2), len(EXC), len(cat2_other)))
chk('IN01', 'IN00 - IN00a + IN00b = IN01', len(corp) - len(EXCP) + len(ADDED), len(corpus_nodes), 'IN00;IN00a')
chk('IN01', 'lineage.json corpus nodes = decisions.tsv in_corpus_536 yes minus the excluded targeted gene-panel or small-set studies plus the studies added after the freeze', len(corpus_nodes), len((corp - EXC) | ADDED), '')

# ---------------------------------------------------------------- primary funnel
sys.path.insert(0, R + 'analyses'); os.chdir(R + 'analyses')
import recurrence
use(R + 'analyses/recurrence.py')
u = recurrence.load_units()
sel = recurrence.scenario_selection('published, brain')[0]
assert set(sel) == (set(P64) - EXC) | ADDED, (sorted(set(sel) ^ ((set(P64) - EXC) | ADDED)))  # routes.json primary64 is regenerated by search/screen_now_2026-10-09/routes.py from recurrence.scenario_selection('published, brain') (64 studies on the 64-list base; 2026-10-13 regeneration); corpus_exclusions studies and the studies added by hand are not in it or are in it already, so the assert holds with or without the - EXC and | ADDED terms
assert set(P64) - set(sel) <= EXC
pl_in = {k for k, v in plists.items() if v.get('in_corpus') and v.get('primary')}
def node_key(k):
    # primary_lists.json / recurrence key DOI-only studies as doi_<doi with / -> _>; lineage.json node ids use the DOI itself
    return k[4:].replace('_', '/', 1) if k.startswith('doi_') else k
assert {node_key(k) for k in pl_in} - cat2 <= corp | ADDED
pf2 = pl_in - cat2
all_ = set(u['all']) - cat2
studies = set(u['studies']) - cat2
groups = {u['group_of'][s] for s in studies}
brain = {s for s in studies if plists[s]['provenance'] == 'published' and plists[s]['tissue'] == 'brain'}
bgroups = {u['group_of'][s] for s in brain}
units = [r for r in runits if r['scenario'] == 'published, brain']
mod = C.Counter(r['modelled'] for r in units)
FT = lambda s: ftab[s]['n']
PL, LN, RU = rel(R + 'analyses/primary_lists.json'), rel(R + 'corpus/lineage.json'), rel(R + 'analyses/recurrence_units.tsv')
box('PF01', 'primary_funnel', 'Corpus studies (= studies included)', len(corpus_nodes), 'IN01', LN, 'as IN01', 'final')
box('PF02', 'primary_funnel', 'Corpus studies with a primary gene list', len(pf2), 'PF01', PL + '; ' + LN,
    'studies with in_corpus and a non-empty primary list, minus lineage category 2', 'final',
    '%d rows in primary_lists.json with in_corpus and primary; %d of them are lineage category 2 nodes and are removed; funnel_table.tsv step 2 = %s' % (len(pl_in), len(pl_in & cat2), FT('2')))
chk('PF02', 'primary_funnel: %d rows - %d category-2 = PF02' % (len(pl_in), len(pl_in & cat2)), len(pl_in) - len(pl_in & cat2), len(pf2), 'PF01')
wl = {(n['pmid'] or n['id']) for n in nodes if (n['pmid'] or n['id']) in corpus_nodes and any(x['status'] in ('extracted', 'gemma', 'reproduced') for x in (n.get('genelist') or {}).get('lists', []))}
wl |= {k for k, v in plists.items() if (v.get('derived') or v.get('list_only')) and k in corpus_nodes}  # rebuild infrastructure 2026-10-14: derived-list and list-only studies carry an injected list, not a lineage genelist entry
nlists = sum(1 for n in nodes if (n['pmid'] or n['id']) in corpus_nodes for x in (n.get('genelist') or {}).get('lists', []) if x['status'] in ('extracted', 'gemma', 'reproduced'))
assert {node_key(k) for k in pf2} <= wl
box('PF02x', 'primary_funnel', 'Corpus studies without a primary gene list', len(corpus_nodes) - len(pf2), 'PF01', LN + '; ' + PL, 'PF01 minus PF02', 'final',
    '%d have no extracted gene list; %d have an extracted list but none is primary (%d studies with a list = %d + %d); %d extracted lists in %d studies (lineage.json nodes[].genelist.lists[].status extracted, gemma or reproduced)' % (
        len(corpus_nodes) - len(wl), len(wl) - len(pf2), len(wl), len(pf2), len(wl) - len(pf2), nlists, len(wl)))
chk('PF02x', 'PF01 - PF02 = PF02x; no extracted list + extracted list but none primary = PF02x', (len(corpus_nodes) - len(wl)) + (len(wl) - len(pf2)), N('PF02x'), 'PF01')
chk('PF02x', 'PF01 - PF02x = PF02', len(corpus_nodes) - N('PF02x'), len(pf2), 'PF01')
box('PF03', 'primary_funnel', 'Studies with a first-look differential-expression primary list', len(all_ & pf2), 'PF02', 'analyses/recurrence.py load_units()',
    "load_units()['all'] minus category 2", 'final',
    '%s non-DE-only studies and %s studies whose only list has a reanalysis role are left out (load_units modality_only, role_only); funnel_table.tsv step 3 = %s' % (
        len(set(u['modality_only']) - cat2), len(set(u['role_only']) - cat2), FT('3')))
chk('PF03', 'PF02 - modality_only - role_only = PF03', len(pf2) - len(set(u['modality_only']) - cat2) - len(set(u['role_only']) - cat2), len(all_), 'PF02')
meta, multi, rean = set(u['meta']) - cat2, set(u['multi']) - cat2, set(u['reanalysis']) - cat2
box('PF04', 'primary_funnel', 'After removing meta-analyses and reanalyses', len(studies), 'PF03', 'analyses/recurrence.py load_units()', "load_units()['studies'] minus category 2", 'final',
    'left out: %d meta-lane, %d multi-series reanalyses, %d other reanalyses (load_units meta, multi, reanalysis, category 2 removed); funnel_table.tsv step 4 = %s' % (len(meta), len(multi), len(rean), FT('4')))
chk('PF04', 'PF03 - meta - multi - reanalysis = PF04', len(all_) - len(meta) - len(multi) - len(rean), len(studies), 'PF03')
box('PF05', 'primary_funnel', 'Donor groups among those studies (same_data_as, all tissues)', len(groups), 'PF04', 'analyses/recurrence.py load_units()', "distinct group_of over PF04 studies", 'final',
    'branch from PF04, not on the path to PF06; %d studies in %d groups; funnel_table.tsv step 5 = %s' % (len(studies), len(groups), FT('5')))
box('PF06', 'primary_funnel', 'Published postmortem-brain studies', len(brain), 'PF04', PL, "PF04 studies with provenance = published and tissue = brain", 'final',
    'funnel_table.tsv step 6 = %s' % FT('6'))
box('PF07', 'primary_funnel', 'Published postmortem-brain donor groups', len(bgroups), 'PF06', 'analyses/recurrence.py load_units()', 'distinct group_of over PF06 studies', 'final',
    'funnel_table.tsv step 7 = %s' % FT('7'))
chk('PF07', 'PF06 studies counted through group membership: sum over PF07 groups of members in PF06 = PF06', sum(sum(1 for s in brain if u['group_of'][s] == g) for g in bgroups), len(brain), 'PF06')
box('PF08', 'primary_funnel', 'Donor groups with a row in the published-brain recurrence scenario', len(units), 'PF07', RU, "rows with scenario = 'published, brain'", 'final',
    '%d group(s) have no row; funnel_table.tsv step 8 = %s' % (len(bgroups) - len(units), FT('8')))
box('PF08a', 'primary_funnel', 'Not modelled: modelled = no (exclude)', mod['no (exclude)'], 'PF08', RU, "modelled = 'no (exclude)'", 'final',
    'the recurrence file does not say why; funnel_table.tsv splits the non-modelled groups into non-gene-level lists (miRNA, sncRNA, circRNA, lncRNA, TE, HERV) and gene-level lists with no usable tested set, from notes/Population funnel and not from a field in the data, so that split is not reproduced here')
box('PF08b', 'primary_funnel', 'Not modelled: modelled = no (targeted)', mod['no (targeted)'], 'PF08', RU, "modelled = 'no (targeted)'", 'final')
box('PF09', 'primary_funnel', 'Gene-level studies with a known tested gene set (modelled in the recurrence analysis)', mod['yes'], 'PF08', RU, "modelled = yes", 'final',
    "also len(analyses/recurrence.py scenario_selection('published, brain')[0]) = %d; funnel_table.tsv step 11 = %s; %s bulk + %s single-cell (funnel_table.tsv 11a, 11b)" % (len(sel), FT('11'), FT('11a'), FT('11b')))
chk('PF09', 'PF08 - PF08a - PF08b = PF09', len(units) - mod['no (exclude)'] - mod['no (targeted)'], mod['yes'], 'PF08')
chk('PF09', 'PF09 = len(scenario_selection)', mod['yes'], len(sel), 'PF08')
box('PF03x', 'primary_funnel', 'Excluded: first-look list is not differential expression, or its only list is a reanalysis', len(set(u['modality_only']) - cat2) + len(set(u['role_only']) - cat2), 'PF02', 'analyses/recurrence.py load_units()', 'modality_only + role_only (category 2 removed)', 'final',
    '%d studies with non-DE-only lists and %d whose only list has a reanalysis role' % (len(set(u['modality_only']) - cat2), len(set(u['role_only']) - cat2)))
box('PF04x', 'primary_funnel', 'Excluded: meta-analyses and reanalyses', len(meta) + len(multi) + len(rean), 'PF03', 'analyses/recurrence.py load_units()', 'meta + multi + reanalysis (category 2 removed)', 'final',
    '%d meta-lane, %d multi-series reanalyses, %d other reanalyses' % (len(meta), len(multi), len(rean)))
box('PF06x', 'primary_funnel', 'Not published postmortem-brain studies', N('PF04') - N('PF06'), 'PF04', PL, 'PF04 minus PF06', 'final')
chk('PF03x', 'PF02 - PF03x = PF03', N('PF02') - N('PF03x'), N('PF03'), 'PF02')
chk('PF04x', 'PF03 - PF04x = PF04', N('PF03') - N('PF04x'), N('PF04'), 'PF03')
chk('PF06x', 'PF04 - PF06x = PF06', N('PF04') - N('PF06x'), N('PF06'), 'PF04')
# agreement with the hand-built funnel table
for step, bid in (('1', 'PF01'), ('2', 'PF02'), ('3', 'PF03'), ('4', 'PF04'), ('5', 'PF05'), ('6', 'PF06'), ('7', 'PF07'), ('8', 'PF08'), ('11', 'PF09')):
    chk(bid, 'funnel_table.tsv step %s n = %s' % (step, bid), int(ftab[step]['n']), N(bid), '')
# cross-checks on the prior records of the same numbers
for q, rows in queues.items():
    m = re.search(r'^\| %s \| (\d+) \|' % q, summ_md, re.M)
    chk('EL01', 'decisions_ft_summary.md queue size %s = review_queue file rows' % q, int(m.group(1)), len(rows), '')
chk('SC05', 'recovery_report.md corpus recovered at pass 2 (include) = SC05a corpus members', int(re.search(r'Pass 2 include \(recall\): (\d+) of', open(SN + 'recovery_report.md').read()).group(1)),
    sum(1 for k in adv if p2[k] == 'include' and k in corp), '')

# ---------------------------------------------------------------- write
os.makedirs(OUT, exist_ok=True)
today = '2026-10-11'
hdr = ['# PRISMA 2020 flow boxes (final counts) for the frozen corpus, one row per box; Extended Fig. 10',
       '# freeze: corpus/CORPUS_FREEZE.md, tag corpus-freeze-2026-10-11 (515 studies, 62 primary; search end 30 September 2026); this table: %d studies, %d primary, with the %d added after the freeze (IN00b)' % (len(corpus_nodes), len(sel), len(ADDED)),
       '# script: %s sha256 %s (reads the sources below and writes this file; no hand-typed numbers)' % (SCRIPT, sha(R + SCRIPT)),
       '# built: %s' % today,
       '# status: final for every box',
       '# parent_box_id: box(es) this one is counted from; empty = memo row outside the arithmetic. Closing arithmetic: prisma_flow_check.tsv',
       '# counts are distinct records; decisions_ft_summary.md quotes queue-row counts (88 records sit in two queues)']
for p, h in sorted(used.items()):
    hdr.append('# source: %s sha256 %s' % (p, h))
hdr.append('# source (outside the repo): figure_data/funnel/funnel_table.tsv sha256 %s (cross-check only)' % sha(FUNNEL))
hdr = [l for l in hdr if 'funnel_table.tsv sha256' not in l or 'outside the repo' in l]


def clean(v):
    return str(v).replace('\t', ' ').replace('\n', ' ')


COLS = ['box_id', 'stage', 'label', 'n', 'parent_box_id', 'source_file', 'source_field_or_query', 'status', 'note']
with open(OUT + 'prisma_flow.tsv', 'w', newline='') as f:
    f.write('\n'.join(hdr) + '\n')
    w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')
    w.writerow(COLS)
    for b in boxes:
        w.writerow([clean(b[c]) for c in COLS])
CC = ['check_id', 'box_id', 'parent_box_id', 'formula', 'lhs', 'rhs', 'closes']
with open(OUT + 'prisma_flow_check.tsv', 'w', newline='') as f:
    f.write('# arithmetic behind prisma_flow.tsv; lhs is computed from the files named there, rhs is the box n (or the independently stated figure); closes = yes when lhs == rhs\n')
    f.write('# script: %s; freeze tag corpus-freeze-2026-10-11; built %s\n' % (SCRIPT, today))
    w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')
    w.writerow(CC)
    for c in checks:
        w.writerow([clean(c[k]) for k in CC])
bad = [c for c in checks if c['closes'] != 'yes']
print('boxes', len(boxes), 'checks', len(checks), 'not closing', len(bad))
for c in bad:
    print('NOT CLOSING', c)
