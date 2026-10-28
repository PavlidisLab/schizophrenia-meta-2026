#!/usr/bin/env python3
"""Fig. 2b: Spearman correlation between dataset Jaccard and log2 fold of the pair overlap over the DE-prior expectation, 15 brain meta-analyses (105 pairs).
Fold definitions: (A) Extended Methods 5: log2((overlap + 0.5)/(expected_deprior + 0.5)), pairs with expected_deprior > 0 (104); (B) same, all 105; (C) log2(fold_deprior), pairs with overlap > 0 and expected > 0 (93).
Permutation: relabel the 15 meta-analyses (rows and columns of the 15 x 15 matrix of log folds permuted jointly, Jaccard fixed), 10,000 permutations, seed 20260925; P = (1 + #(null rho >= observed))/(10,001), one-sided, over pairs defined in the permuted matrix.
Input: figures/data/fig02/fig02_e_meta_analysis_pairs.tsv of the manuscript repo (read-only)."""
import pandas as pd, numpy as np, hashlib, sys
from scipy.stats import spearmanr
P='figures/data/fig02/fig02_e_meta_analysis_pairs.tsv'
d=pd.read_csv(P,sep='\t',comment='#'); sha=hashlib.sha256(open(P,'rb').read()).hexdigest()[:12]
d['lf_pc']=np.log2((d.overlap+0.5)/(d.expected_deprior+0.5)); d['lf_raw']=np.log2(d.fold_deprior.where(d.fold_deprior>0))
labs=sorted(set(d.label_a)|set(d.label_b)); ix={l:i for i,l in enumerate(labs)}; n=len(labs); a=d.label_a.map(ix).values; b=d.label_b.map(ix).values; iu=np.triu_indices(n,1)
def run(vals,nperm=10000,seed=20260925):
    F=np.full((n,n),np.nan); J=np.full((n,n),np.nan)
    for i,j,f,jj in zip(a,b,vals,d.dataset_jaccard.values): F[i,j]=F[j,i]=f; J[i,j]=J[j,i]=jj
    ok=~np.isnan(J[iu])&~np.isnan(F[iu]); obs=spearmanr(J[iu][ok],F[iu][ok])[0]; rng=np.random.default_rng(seed); null=[]
    for _ in range(nperm):
        p=rng.permutation(n); Fp=F[np.ix_(p,p)]; okp=~np.isnan(J[iu])&~np.isnan(Fp[iu]); null.append(spearmanr(J[iu][okp],Fp[iu][okp])[0])
    null=np.array(null); return obs,int(ok.sum()),(1+np.sum(null>=obs))/(nperm+1),float(null.mean()),float(null.std())
rows=[]
for name,vals in (('A: log2((overlap+0.5)/(E+0.5)), E>0',np.where(d.expected_deprior>0,d.lf_pc,np.nan)),('B: same, all pairs',d.lf_pc.values),('C: log2(fold_deprior), overlap>0 and E>0',d.lf_raw.values)):
    o,k,p,m,s=run(vals); rows.append(dict(definition=name,n_pairs=k,n_meta_analyses=n,spearman_rho=round(o,4),perm_P=round(p,4),null_mean=round(m,4),null_sd=round(s,4)))
R=pd.DataFrame(rows)
with open('fig2b_meta_rho.tsv','w') as f:
    f.write(f'# NEW 2026-10-18 recomputation of the Fig. 2b meta-analysis correlation; script fig2b_meta_rho.py; input sha256 {sha}; seed 20260925; 10,000 permutations of meta-analysis labels\n'); R.to_csv(f,sep='\t',index=False)
print(R.to_string(index=False))
