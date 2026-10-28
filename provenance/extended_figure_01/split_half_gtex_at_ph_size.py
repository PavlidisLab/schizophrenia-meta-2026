#!/usr/bin/env python3
"""GTEx slow-vs-fast-death axis: split-half reliability at the sample size of the pH collections (50 and 45 donors).
Run (remote environment; data are there): python3 split_half_gtex_at_ph_size.py [ANNOT_DIR]  -> files beside this script.
Axis construction is imported-by-copy from analyses/residual_after_ph_agonal_2026-10-09/_remotereliability/scripts/gtex_reliability.py
(per region OLS of log2(TPM+1) on Hardy-4 indicator, ischemic hours, age, sex; genes with region median TPM >= 1; axis = mean of the
three regions' Hardy t over genes present in all three). Differences from that script: halves are PLAIN random halves
(floor(n/2) and the rest), as in analyses/ph_axis_split_half_2026-10-11/split_half_ph_axis.R: no stratification by Hardy group or sex,
no age matching. The donor universe is the same (Hardy 2 or 4 with an ischemic time in at least one of the three regions)."""
import sys, os, hashlib, datetime, numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE=os.path.dirname(os.path.abspath(__file__))
A=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'..','residual_after_ph_agonal_2026-10-09','_remotereliability','gtex_annotations')
NSPLIT=200; NSUB=20; SEED=20261011; SIZES=[50,45]
TP='/space/scratch/npatadia/GTex_data_V10_dec2024/gene_tpm_v10_%s.gct.gz'
REG={'BA9':'brain_frontal_cortex_ba9','BA24':'brain_anterior_cingulate_cortex_ba24','Cortex':'brain_cortex'}
fa=f'{A}/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt'; fs=f'{A}/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt'
att=pd.read_csv(fa,sep='\t',usecols=['SAMPID','SMTSISCH','SMRIN']); subj=pd.read_csv(fs,sep='\t')
subj['age']=subj.AGE.map(lambda a: np.mean([int(x) for x in a.split('-')])); subj['male']=(subj.SEX==1).astype(float)
data={}
for k,f in REG.items():
    hdr=pd.read_csv(TP%f,sep='\t',skiprows=2,nrows=0).columns
    t=pd.read_csv(TP%f,sep='\t',skiprows=2,dtype={c:np.float32 for c in hdr[2:]})
    s=att[att.SAMPID.isin(hdr[2:])&att.SMTSISCH.notna()].copy(); s['SUBJID']=s.SAMPID.str.split('-').str[:2].str.join('-')
    s=s.merge(subj,on='SUBJID'); s['isch']=s.SMTSISCH/60.0
    keep=np.median(t[list(s.SAMPID)].values,axis=1)>=1
    sym=t.Description.values[keep]; X=np.log2(t[list(s.SAMPID)].values[keep]+1).astype(np.float64)
    ix=pd.Series(sym).groupby(sym).head(1).index
    s=s.reset_index(drop=True); s['in']=s.DTHHRDY.isin([2,4])
    data[k]=dict(s=s,X=X[ix.values],sym=sym[ix.values],subj=s.SUBJID.values,h4=(s.DTHHRDY.values==4).astype(float),
                 M0=np.column_stack([s.isch.values,s.age.values,s.male.values]),ok=s['in'].values)
    print(k,'samples',len(s),'genes',len(ix),flush=True)
def ols_t(Y,M):
    Q,Rr=np.linalg.qr(M); beta=np.linalg.solve(Rr,Q.T@Y.T); res=Y.T-M@beta; dof=M.shape[0]-M.shape[1]
    s2=(res**2).sum(0)/dof; cd=np.diag(np.linalg.inv(Rr.T@Rr))[1]; return beta[1]/np.sqrt(s2*cd)
def axis_of(donors):
    T={}
    for k,d in data.items():
        idx=np.where(np.isin(d['subj'],list(donors))&d['ok'])[0]
        M=np.column_stack([np.ones(len(idx)),d['h4'][idx],d['M0'][idx]])
        if len(idx)<=M.shape[1] or np.linalg.matrix_rank(M)<M.shape[1]: raise ValueError(k)
        T[k]=pd.Series(ols_t(d['X'][:,idx],M),index=d['sym'])
    T=pd.DataFrame(T).dropna(); return T.mean(axis=1)
alld=sorted(set.union(*[set(d['s'].loc[d['s']['in'],'SUBJID']) for d in data.values()]))
dinfo=subj.set_index('SUBJID').loc[alld,['DTHHRDY','male','age']]
slow=list(dinfo.index[dinfo.DTHHRDY==4]); fast=list(dinfo.index[dinfo.DTHHRDY==2]); N=len(alld)
# ---- donor counts
rows=[('all regions (any)','donors',N,len(slow),len(fast))]
inall=set.intersection(*[set(d['s'].loc[d['s']['in'],'SUBJID']) for d in data.values()])
rows.append(('in all three regions','donors',len(inall),int((dinfo.loc[sorted(inall)].DTHHRDY==4).sum()),int((dinfo.loc[sorted(inall)].DTHHRDY==2).sum())))
for k,d in data.items():
    s=d['s'][d['s']['in']]; rows.append((k,'donors with Hardy 2/4 and ischemic time',s.SUBJID.nunique(),s[s.DTHHRDY==4].SUBJID.nunique(),s[s.DTHHRDY==2].SUBJID.nunique()))
rng=np.random.default_rng(SEED)
def sub_ratio(n):
    ns=int(round(n*len(slow)/N)); return ns,n-ns
def draw(n):
    ns,nf=sub_ratio(n); return set(rng.choice(slow,ns,replace=False))|set(rng.choice(fast,nf,replace=False))
REJ=[0]
def split_rel(pool):
    while True:
        try: return split_rel1(pool)
        except ValueError: REJ[0]+=1   # a half with a rank-deficient region design (e.g. no slow-death donor in that region) is redrawn
def split_rel1(pool):
    pool=np.array(sorted(pool)); a=set(rng.choice(pool,len(pool)//2,replace=False)); b=set(pool)-a
    x=axis_of(a); y=axis_of(b); g=sorted(set(x.index)&set(y.index)); r=spearmanr(x[g],y[g])[0]; return r,2*r/(1+r),len(g),len(a),len(b)
import time; t0=time.time()
out={}  # label -> list of (sub, r, rsb)
# full-donor plain random halves
full=[split_rel(alld) for _ in range(NSPLIT)]; print('rejected splits so far',REJ[0],'full donors done',round(time.time()-t0),'s',flush=True)
out['full']=[(0,)+f[:2] for f in full]; ngenes=full[0][2]
perreg=[]
for n in SIZES:
    for u in range(NSUB):
        pool=draw(n); res=[split_rel(pool) for _ in range(NSPLIT)]
        out[(n,u)]=[(u,)+f[:2] for f in res]
        for k,d in data.items():
            s=d['s'][d['s']['in']&d['s'].SUBJID.isin(pool)]; perreg.append((n,u,k,s.SUBJID.nunique()))
        print('rej',REJ[0],n,u,round(np.mean([f[1] for f in res]),3),round(time.time()-t0),'s',flush=True)
# ---- tables
def q(v): return np.quantile(v,[0.025,0.975])
pr=pd.DataFrame(perreg,columns=['n','sub','region','donors'])
for n in SIZES:
    ns,nf=sub_ratio(n)
    for k in REG:
        v=pr[(pr.n==n)&(pr.region==k)].donors; rows.append((f'subsample of {n} ({ns} slow, {nf} fast): {k}','mean donors per subsample (of the 20) with Hardy 2/4 and ischemic time',round(v.mean(),1),'',''))
pd.DataFrame(rows,columns=['group','quantity','donors_total','slow_death_hardy4','fast_death_hardy2']).to_csv(f'{HERE}/donor_counts.tsv',sep='\t',index=False)
sh=lambda p: hashlib.sha256(open(p,'rb').read()).hexdigest()[:16]
inp=[fa,fs]+[TP%f for f in REG.values()]
hdr=['# figure: figures_data_quality Supplementary Fig. 1f, GTEx row at the sample size of the pH collections (panel_c.tsv columns)',
 '# script: analyses/gtex_split_at_ph_size_2026-10-11/split_half_gtex_at_ph_size.py',
 f'# seed: {SEED}; {NSPLIT} plain random splits per donor set; {NSUB} subsample redraws per size; run {datetime.date.today()}',
 '# input sha256 (first 16 hex): '+'; '.join(f'{os.path.basename(p)} {sh(p)}' for p in inp),
 '# lo/hi = 2.5% and 97.5% quantiles of the Spearman-Brown values over splits (pooled over the 20 subsamples where stated)']
R=[]; sm=[]
def what(desc,extra=''):
    return f'GTEx death-speed axis (mean Hardy-4-vs-2 t over BA9, BA24, Cortex; {ngenes} genes); {desc}; Spearman-Brown reliability, plain random halves (no stratification by group or sex, no age matching) (mean; 2.5% to 97.5% of splits){extra}'
v=np.array([x[2] for x in out['full']]); l,h=q(v)
R.append((f'GTEx axis, all {N} donors, plain random halves',v.mean(),l,h,what(f'all {N} donors with Hardy 2 or 4 ({len(slow)} slow, {len(fast)} fast), {NSPLIT} splits into halves of {N//2} and {N-N//2}')))
for n in SIZES:
    ns,nf=sub_ratio(n); allv=np.array([x[2] for u in range(NSUB) for x in out[(n,u)]])
    means=np.array([np.mean([x[2] for x in out[(n,u)]]) for u in range(NSUB)])
    for u in range(NSUB):
        v=np.array([x[2] for x in out[(n,u)]]); sm.append((n,u,v.mean(),*q(v),np.mean([x[1] for x in out[(n,u)]])))
    v0=np.array([x[2] for x in out[(n,0)]]); l,h=q(v0)
    R.append((f'GTEx axis, {n} donors (one subsample), plain random halves',v0.mean(),l,h,what(f'first random subsample of {n} donors ({ns} slow, {nf} fast), {NSPLIT} splits into halves of {n//2} and {n-n//2}')))
    l,h=q(allv)
    R.append((f'GTEx axis, {n} donors (mean over {NSUB} subsamples), plain random halves',means.mean(),l,h,what(f'{NSUB} random subsamples of {n} donors ({ns} slow, {nf} fast), {NSPLIT} splits each; reliability = mean of the {NSUB} subsample means (range {means.min():.3f} to {means.max():.3f}); 2.5% to 97.5% of all {NSUB*NSPLIT} splits pooled')))
hdr.append(f'# splits redrawn because a half left a region design rank-deficient: {REJ[0]} (of {len(SIZES)*NSUB*NSPLIT+NSPLIT+REJ[0]} draws)')
open(f'{HERE}/split_half_gtex_at_ph_size.tsv','w').write('\n'.join(hdr)+'\n')
pd.DataFrame(R,columns=['measure','reliability','lo','hi','what']).to_csv(f'{HERE}/split_half_gtex_at_ph_size.tsv',sep='\t',index=False,mode='a',float_format='%.4g')
pd.DataFrame(sm,columns=['n','subsample','mean_rho_SB','q2.5','q97.5','mean_r_spearman']).to_csv(f'{HERE}/split_half_gtex_at_ph_size_subsamples.tsv',sep='\t',index=False)
pd.DataFrame([(('full',0) if k=='full' else k)+x for k,l_ in out.items() for x in l_ for k in [k]] if False else [ (str(k),)+x for k,l_ in out.items() for x in l_ ],columns=['set','subsample','r_spearman','rho_SB_spearman']).to_csv(f'{HERE}/split_half_gtex_at_ph_size_splits.tsv',sep='\t',index=False)
print(pd.DataFrame(R,columns=['m','rel','lo','hi','w']).drop(columns='w').to_string())
