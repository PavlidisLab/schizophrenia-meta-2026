#!/usr/bin/env python3
"""Recomputes the A1 counts (node check, blank N counts, covariate readers, curated edges) into tables/.
Run 2026-10-19; no randomness. Base: the 509-study corpus (node check, blank counts, curated edges) and the core-80 base (covariate readers).
"""
import csv, glob, json, collections, hashlib, os
H='schizophrenia-transcriptome-replication/'
M='schizophrenia-transcriptome-reproducibility-manuscript/figures/data/fig01/'
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sha=lambda p: hashlib.sha256(open(p,'rb').read()).hexdigest()[:12]
def tsv(p): return list(csv.DictReader([l for l in open(p) if not l.startswith('#')],delimiter='\t'))
def out(name,inputs,rows,head):
    with open(f'{A}/tables/{name}','w') as f:
        f.write(f'# script: A/scripts/a1_counts.py; date 2026-10-19; seed none\n')
        for p in inputs: f.write(f'# input: {p} sha256 {sha(p)}\n')
        f.write('# base: 509-study corpus (all but covariate readers); covariate readers: core-80\n')
        f.write('\t'.join(head)+'\n')
        for r in rows: f.write('\t'.join(str(x) for x in r)+'\n')
# node check
fs=sorted(glob.glob(H+'corpus/node_check_2026-09-27/checked_*.tsv')); rows=[]
for f in fs: rows+=tsv(f)
R=[('studies',len({r['id'] for r in rows})),('field_checks',len(rows))]
for k,v in sorted(collections.Counter(r['field'] for r in rows).items()): R.append((f'rows_field_{k}',v))
R.append(('rows_with_quote',sum(1 for r in rows if r['quote'].strip()))); R.append(('rows_without_quote',sum(1 for r in rows if not r['quote'].strip())))
for k,v in sorted(collections.Counter(r['field'] for r in rows if not r['quote'].strip()).items()): R.append((f'no_quote_field_{k}',v))
for k,v in sorted(collections.Counter(r['verdict'] for r in rows).items()): R.append((f'verdict_{k}',v))
out('a1_node_check_counts.tsv',fs,R,['quantity','value'])
# blanks
p=M+'fig01_b_corpus_studies.tsv'; r=tsv(p); b=lambda x,k:x[k].strip()==''
R=[('studies',len(r)),('cases_blank',sum(b(x,'n_cases') for x in r)),('controls_blank',sum(b(x,'n_controls') for x in r)),('total_blank',sum(b(x,'n_total') for x in r)),
   ('cases_or_controls_blank',sum(b(x,'n_cases') or b(x,'n_controls') for x in r)),('cases_and_controls_blank',sum(b(x,'n_cases') and b(x,'n_controls') for x in r))]
out('a1_blank_counts.tsv',[p],R,['quantity','value'])
# covariate readers
p=M+'fig01_f_covariate_calls_80.tsv'; r=tsv(p); st=collections.defaultdict(list)
for x in r: st[x['pmid']].append(x)
one=[k for k,v in st.items() if v[0]['coding_two_readers']!='yes']
no2=[st[k][0]['label'] for k,v in st.items() if k not in one and all(x['reader2_call']=='no_second_reader' for x in v)]
no1=[st[k][0]['label'] for k,v in st.items() if all(x['reader1_call']=='(no reading)' for x in v)]
R=[('studies',len(st)),('one_reader_studies (coding_two_readers != yes)',len(one)),('coded_by_two',len(st)-len(one)),('coded_by_two_but_reader2_absent',';'.join(no2)),('coded_by_two_but_reader1_absent',';'.join(no1)),
   ('studies_with_two_usable_readings',len(st)-len(one)-len(no2)-len(no1)),('rows_reader2_no_text',sum(1 for x in r if x['reader2_call']=='no_text')),('rows_reader2_no_call',sum(1 for x in r if x['reader2_call']=='no_call'))]
out('a1_covariate_readers.tsv',[p],R,['quantity','value'])
# curated edges
ce=json.load(open(H+'corpus/curated_edges.json'))['edges']; sc=json.load(open(H+'corpus/sc_banks.json'))['edges']; p=M+'fig01_c_lineage_edges.tsv'; r=tsv(p)
R=[('curated_edges.json edges',len(ce)),('curated_edges.json with evidence',sum(1 for e in ce if e.get('evidence','').strip())),('sc_banks.json edges',len(sc)),('sc_banks.json with evidence',sum(1 for e in sc if e.get('evidence','').strip())),
   ('register rows',len(r)),('register curated=yes',sum(1 for x in r if x['curated']=='yes')),('register curated=yes with evidence',sum(1 for x in r if x['curated']=='yes' and x['evidence'].strip()))]
out('a1_curated_edges.tsv',[H+'corpus/curated_edges.json',H+'corpus/sc_banks.json',p],R,['quantity','value'])
