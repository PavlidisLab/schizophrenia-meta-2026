"""F_notes_10-12: four small recounts from files already on disk (no random numbers, seed: none). Run: python3 -I f_misc_checks.py (writes ../tables/*.tsv). Date 2026-10-19.
1 antipsychotic_counts_rules.tsv   Note 10 'N of 80 studies in K of 18 clusters' under each rule, from figure_data_refresh_80_studies/refresh_proteomics_and_ph_data/new_files/antipsychotic_record_counts_core80.tsv
2 shared_pathways_by_group.tsv     Note 11: share of the 710 shared top pathways (sum over pairs) carried by the oxphos / myelination / microglia term groups; H01 a5b_terms.tsv + M09 term_groups.tsv
3 periph_lane_composition.tsv      Note 11: tissue_short of the 'periph' lane nodes in corpus/lineage.json (lane definition used by the Deans overlap)
4 deans_overlap_383.tsv            Note 11: genes of Deans sheet A1 in the 383-gene study set (recurrent_sets_core80.tsv, set study_brain)
Base: core-80 (80 studies, 72 lists, 71 qualifying, 383 genes)."""
import csv, math, itertools, collections, json, hashlib, os
R = os.path.expanduser(''); OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tables') + '/'
def sha(p): return hashlib.sha256(open(os.path.expanduser(p), 'rb').read()).hexdigest()[:12]
def rd(p): return list(csv.DictReader((l for l in open(p) if not l.startswith('#')), delimiter='\t'))
def w(fn, hdr, rows, rowsfmt=None):
    with open(OUT + fn, 'w') as f:
        for h in hdr: f.write('# ' + h + '\n')
        f.write('\n'.join('\t'.join(str(c) for c in r) for r in rows) + '\n')
# 1
P1 = R + 'figure_data_refresh_80_studies/refresh_proteomics_and_ph_data/new_files/antipsychotic_record_counts_core80.tsv'; A = rd(P1)
def clu(r): return r['label'] if r['cohort_status_table4'] in ('independent', 'multi-bank') else r['cohort_status_table4']
rules = [('strict: class testable (rule of counts.tsv)', {'testable'}), ('strict + Chen 2013 (GSE35977, the dataset the test was run in)', {'testable', 'testable_via_SMRI_workbook'}),
         ('strict + testable_provisional (Table 4 text only)', {'testable', 'testable_provisional'}), ('strict + Chen 2013 + provisional', {'testable', 'testable_via_SMRI_workbook', 'testable_provisional'}),
         ('any per-sample record (testable* + record_per_sample_not_testable + record_in_article_table_expression_public_not_counted)', {'testable', 'testable_via_SMRI_workbook', 'testable_provisional', 'record_per_sample_not_testable', 'record_in_article_table_expression_public_not_counted'})]
rows = []
for nm, cl in rules:
    S = [r for r in A if r['class'] in cl]; rows.append((nm, len(S), len(set(clu(r) for r in S if r['cohort_status_table4']!='multi-bank')), len(A), 18))
w('antipsychotic_counts_rules.tsv', ['script: f_misc_checks.py part 1', 'input: %s sha256 %s' % (P1.replace(R, ''), sha(P1)), 'clusters: 5 named (stanley, cmc+mssm+pitt+libd, aus, harvard, pritzker) + 13 independent studies = 18; the 7 multi-bank studies sit in no cluster and are left out of the cluster count (the studies are still counted)'],
  [('rule', 'n_studies', 'n_clusters', 'of_n_studies', 'of_n_clusters')] + rows)
# 2
H = R + 'rerun_on_80_studies/H01/'; T = rd(H + 'a5b_terms.tsv'); G = rd(R + 'rerun_on_80_studies/M09/results_seedA_noextra_new/term_groups.tsv'); g = {x['term']: x['group'] for x in G}
tot = sum(math.comb(int(t['n_lists_top']), 2) for t in T); by = collections.Counter(); nt = collections.Counter()
for t in T:
    if t['term'] in g: by[g[t['term']]] += math.comb(int(t['n_lists_top']), 2); nt[g[t['term']]] += 1
agr = {r['metric']: r for r in rd(H + 'a5b_agreement.tsv')}['shared_pathways_total']
nullm = float(agr['deprior_null_mean']); obs = float(agr['observed'])
rows = [('all terms in >= 2 lists', len(T), tot, '', '')] + [(k, nt[k], by[k], '%.1f%%' % (100 * by[k] / tot), '') for k in ('oxphos', 'myelin', 'microglia')]
rows += [('oxphos + myelin', nt['oxphos'] + nt['myelin'], by['oxphos'] + by['myelin'], '%.1f%%' % (100 * (by['oxphos'] + by['myelin']) / tot), 'share of excess over DE-prior null mean %.1f: %.1f%% (assumes the null contributes about 0 for these terms: deprior_null_mean 0.01 to 0.03 lists per term)' % (nullm, 100 * (by['oxphos'] + by['myelin']) / (obs - nullm)))]
w('shared_pathways_by_group.tsv', ['script: f_misc_checks.py part 2', 'input: H01/a5b_terms.tsv sha256 %s; H01/a5b_agreement.tsv sha256 %s; M09/results_seedA_noextra_new/term_groups.tsv sha256 %s (all under rerun_on_80_studies)' % (sha(H + 'a5b_terms.tsv'), sha(H + 'a5b_agreement.tsv'), sha(R + 'rerun_on_80_studies/M09/results_seedA_noextra_new/term_groups.tsv')),
  'validation: sum over terms of C(n_lists_top,2) = %d equals shared_pathways_total observed %s in a5b_agreement.tsv' % (tot, agr['observed'])],
  [('group', 'n_terms', 'shared_pathway_pairs', 'share_of_total', 'note')] + rows)
# 3
L = os.path.expanduser('schizophrenia-transcriptome-replication/corpus/lineage.json'); nodes = json.load(open(L))['nodes']; p = [n for n in nodes if n.get('lane') == 'periph']
c = collections.Counter(n.get('tissue_short') or '(none)' for n in p)
w('periph_lane_composition.tsv', ['script: f_misc_checks.py part 3', 'input: schizophrenia-transcriptome-replication/corpus/lineage.json sha256 %s; lane label of all %d periph nodes: %s' % (sha(L), len(p), set(n.get('lane_label') for n in p)),
  'the Deans overlap used lanes ipsc + periph (side-studies/deans-2026/overlap.py); periph = 123 lists from these 130 lineage nodes'], [('tissue_short', 'n_nodes')] + c.most_common())
# 4
S = rd(R + 'rerun_on_80_studies/M04_stanley_laptop/recurrent_sets_core80.tsv'); s383 = {x['gene'] for x in S if x['set'] == 'study_brain'}
A1p = R + 'deans_convergent_genes_2026-10-15/inputs/A1_pooled_metanalyzed_DEGs.tsv'; a1 = {x['MarkerName'] for x in csv.DictReader(open(A1p), delimiter='\t')}
ov = sorted(a1 & s383)
w('deans_overlap_383.tsv', ['script: f_misc_checks.py part 4', 'input: recurrent_sets_core80.tsv sha256 %s (set study_brain, %d genes); A1_pooled_metanalyzed_DEGs.tsv sha256 %s (%d symbols); match on gene symbol as written in both files' % (sha(R + 'rerun_on_80_studies/M04_stanley_laptop/recurrent_sets_core80.tsv'), len(s383), sha(A1p), len(a1))], [('n_overlap', 'genes'), (len(ov), ';'.join(ov))])
for fn in ('antipsychotic_counts_rules', 'shared_pathways_by_group', 'periph_lane_composition', 'deans_overlap_383'): print(open(OUT + fn + '.tsv').read())
