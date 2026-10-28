import os,sys,hashlib,io
import numpy as np,pandas as pd
from scipy import stats
sys.dont_write_bytecode=True
H=os.path.dirname(os.path.abspath(__file__)); O=H+'/../'; W=O+'work/wb/'
A=pd.read_csv(W+'tables/81_C.tsv',sep='\t')                # persons from papers among the 71 (primary, Stanley by workbook key)
B=pd.read_csv(W+'tables/status81_all.tsv',sep='\t')        # all persons; in-71 record used where a person has one
B['in_81']=np.where(B.status=='81','yes','no')
B.to_csv(O+'donor_level/donor_level_flags.tsv',sep='\t',index=False)
for n in ('81_C','81_U','81_L'): pd.read_csv(W+f'tables/{n}.tsv',sep='\t').to_csv(O+f'donor_level/persons_{n}.tsv',sep='\t',index=False)
def diff(x):
    a=x[x.group=='SZ']; b=x[x.group=='CTRL']
    if len(a)<2 or len(b)<2: return dict(n_SZ=len(a),n_CTRL=len(b),raw=np.nan,lo=np.nan,hi=np.nan)
    X=np.column_stack([np.ones(len(x)),(x.group=='SZ').astype(float)]+[ (x.study==s).astype(float).values for s in sorted(x.study.unique())[1:]])
    y=x.ph.values; beta,res,rk,sv=np.linalg.lstsq(X,y,rcond=None); r=y-X@beta; dfree=len(y)-rk; s2=(r@r)/dfree; cov=s2*np.linalg.pinv(X.T@X); se=np.sqrt(cov[1,1]); t=stats.t.ppf(.975,dfree)
    return dict(n_SZ=len(a),n_CTRL=len(b),raw=beta[1],lo=beta[1]-t*se,hi=beta[1]+t*se)
fpA=pd.read_csv(O+'results/bounds/fit_primary.tsv',sep='\t'); fpA=fpA[fpA.table=='81_C'].iloc[0]
fpB=pd.read_csv(O+'results/all/fit_primary.tsv',sep='\t').iloc[0]
pcA=pd.read_csv(O+'results/bounds/fit_primary_per_cohort.tsv',sep='\t'); pcA=pcA[pcA.table=='81_C'].set_index('cohort')
pcB=pd.read_csv(O+'results/all/fit_primary_per_cohort.tsv',sep='\t').set_index('cohort')
raw=open(O+'results/m/models.tsv').read().split('\n'); out=[l for l in raw if not l.startswith(' - Rescale')]
m=pd.read_csv(io.StringIO('\n'.join(out)),sep='\t',comment='#')
e14=m[m.label.str.startswith('only donors that are samples of the 14')].iloc[0]; o14=m[m.label.str.startswith('drop donors that are samples of the 14')].iloc[0]
rows=[]
for c in sorted(A.cohort.unique()):
    a=diff(A[A.cohort==c]); b=diff(B[B.cohort==c]); e=diff(A[(A.cohort==c)&(A.in_14_dataset_regions==1)]); o=diff(A[(A.cohort==c)&(A.in_14_dataset_regions==0)])
    rows.append(dict(cohort=c,
      A_n_SZ=a['n_SZ'],A_n_CTRL=a['n_CTRL'],A_est=a['raw'],A_lo=a['lo'],A_hi=a['hi'],A_model_est=pcA.loc[c,'shrunken'],A_model_lo=pcA.loc[c,'shrunk_lo'],A_model_hi=pcA.loc[c,'shrunk_hi'],
      B_n_SZ=b['n_SZ'],B_n_CTRL=b['n_CTRL'],B_est=b['raw'],B_lo=b['lo'],B_hi=b['hi'],B_model_est=pcB.loc[c,'shrunken'],B_model_lo=pcB.loc[c,'shrunk_lo'],B_model_hi=pcB.loc[c,'shrunk_hi'],
      expr14_n_SZ=e['n_SZ'],expr14_n_CTRL=e['n_CTRL'],expr14_est=e['raw'],expr14_lo=e['lo'],expr14_hi=e['hi'],other_n_SZ=o['n_SZ'],other_n_CTRL=o['n_CTRL'],other_est=o['raw'],other_lo=o['lo'],other_hi=o['hi']))
nS=lambda d,g:int((d.group==g).sum())
rows.append(dict(cohort='all cohorts pooled (primary model)',A_n_SZ=nS(A,'SZ'),A_n_CTRL=nS(A,'CTRL'),A_est=fpA.estimate,A_lo=fpA.profile_lo,A_hi=fpA.profile_hi,A_model_est=fpA.estimate,A_model_lo=fpA.profile_lo,A_model_hi=fpA.profile_hi,
  B_n_SZ=nS(B,'SZ'),B_n_CTRL=nS(B,'CTRL'),B_est=fpB.estimate,B_lo=fpB.profile_lo,B_hi=fpB.profile_hi,B_model_est=fpB.estimate,B_model_lo=fpB.profile_lo,B_model_hi=fpB.profile_hi,
  expr14_n_SZ=int(((A.group=='SZ')&(A.in_14_dataset_regions==1)).sum()),expr14_n_CTRL=int(((A.group=='CTRL')&(A.in_14_dataset_regions==1)).sum()),expr14_est=e14.estimate,expr14_lo=e14.lo,expr14_hi=e14.hi,
  other_n_SZ=int(((A.group=='SZ')&(A.in_14_dataset_regions==0)).sum()),other_n_CTRL=int(((A.group=='CTRL')&(A.in_14_dataset_regions==0)).sum()),other_est=o14.estimate,other_lo=o14.lo,other_hi=o14.hi))
T=pd.DataFrame(rows)
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
hd=("# per_cohort_81: SZ minus control brain pH (pH units) per cohort and pooled. Persons first; no donor-level rows.\n"
"# sets: A = persons whose pH comes from a paper among the 81 analysed studies (the 71 primary studies plus 10 derived-list studies; Stanley persons resolved by the SMRI workbook key, other persons by identifiers then fingerprints); B = all persons of the 32-study set plus the added papers (a person with a record in a paper among the 81 is counted once with that record); expr14 = persons of A who are samples of the 14 dataset-regions of Fig. 4a; other = the rest of A\n"
"# *_est/_lo/_hi for a cohort row = within-cohort difference (least squares with study as fixed effect) with 95% t interval; *_model_* = primary-model estimate shrunken toward the pooled effect with its 95% interval; pooled row = primary-model estimate with profile-likelihood 95% interval (expr14 and other: the same model refitted on those persons)\n"
f"# script: scripts/aggregate81.py (pooled_ph_difference_by_study); sources: work/wb/tables/81_C.tsv sha256 {sha(W+'tables/81_C.tsv')}, tables/status81_all.tsv sha256 {sha(W+'tables/status81_all.tsv')}; models in results/ (R 4.4.1, lme4 2.0-6)\n# built 2026-10-14\n")
with open(O+'per_cohort_81.tsv','w') as f: f.write(hd); T.round(4).to_csv(f,sep='\t',index=False)
# outside_81_share
rows=[]
for c in sorted(B.cohort.unique()):
    x=B[B.cohort==c]; out=x[x.in_81=='no']
    base=dict(cohort=c,persons_all=len(x),persons_in_81=int((x.in_81=='yes').sum()),persons_outside_81=len(out),share_outside_81=round(len(out)/len(x),3),
      est_A=T.set_index('cohort').loc[c,'A_est'],est_A_lo=T.set_index('cohort').loc[c,'A_lo'],est_A_hi=T.set_index('cohort').loc[c,'A_hi'],est_B=T.set_index('cohort').loc[c,'B_est'],est_B_lo=T.set_index('cohort').loc[c,'B_lo'],est_B_hi=T.set_index('cohort').loc[c,'B_hi'])
    if len(out)==0: rows.append({**base,'paper_outside_81':'(none)','status_of_paper':'','persons_from_paper_SZ':0,'persons_from_paper_CTRL':0})
    for (s,st),g in out.groupby(['study','status']):
        rows.append({**base,'paper_outside_81':s,'status_of_paper':st,'persons_from_paper_SZ':int((g.group=='SZ').sum()),'persons_from_paper_CTRL':int((g.group=='CTRL').sum())})
S=pd.DataFrame(rows)
hd2=("# outside_81_share: for each cohort, the persons and share of its pH data that come only from papers outside the 81 analysed studies (Table 5 papers, papers outside both), by paper, and the cohort's estimate from persons in the 81 (A) against all persons (B); within-cohort difference with study as fixed effect, 95% t interval. Persons counted once; a person in a paper among the 81 and in an outside paper is counted in the 81.\n"
"# papers: GSE92538 = Hagenauer 2018 (outside both), GSE208338 = Worf 2025 (outside both); Chen 2013 and Ramaker 2017 are inside the 81\n"
f"# script: scripts/aggregate81.py; source tables/status81_all.tsv sha256 {sha(W+'tables/status81_all.tsv')}; built 2026-10-14\n")
with open(O+'outside_81_share.tsv','w') as f: f.write(hd2); S.round(4).to_csv(f,sep='\t',index=False)
pd.set_option('display.width',250); print(T.round(3).T.to_string()); print(S.round(3).to_string())
