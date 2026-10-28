"""B1 'pH score of a list': option evidence only (not a decision). Mann-Whitney U of signed_ph_score between handling groups, studies with scored_ge10 = yes.
input: manuscript repo figures/data/fig07/fig07_c_ph_score_by_handling_v2.tsv; scipy mannwhitneyu two-sided (exact/asymptotic default); no seed. run 2026-10-19. base core-80 (72 lists), pH axis v2.
Also checks the 223-analysis claim against figures/data/forking_paths_v2/forking_v2_c_ph_paths.tsv."""
import pandas as pd, hashlib
from scipy.stats import mannwhitneyu
M='schizophrenia-transcriptome-reproducibility-manuscript/figures/data/'
f=M+'fig07/fig07_c_ph_score_by_handling_v2.tsv'; print('sha256',hashlib.sha256(open(f,'rb').read()).hexdigest()[:12])
d=pd.read_csv(f,sep='\t',comment='#'); print('rows',len(d),'scored_ge10 yes',(d.scored_ge10=='yes').sum())
s=d[d.scored_ge10=='yes']; print(s.ph_call.value_counts().to_dict())
rows=[]
def cmp(name,a,b):
    u=mannwhitneyu(a,b,alternative='two-sided'); rows.append((name,len(a),len(b),a.median(),b.median(),u.statistic,u.pvalue))
adj=s[s.ph_call=='adjusted covariate'].signed_ph_score
cmp('adjusted covariate vs silent',adj,s[s.ph_call=='silent'].signed_ph_score)
cmp('adjusted covariate vs all other scored',adj,s[s.ph_call!='adjusted covariate'].signed_ph_score)
cmp('adjusted or matched vs silent',s[s.ph_call.isin(['adjusted covariate','matched'])].signed_ph_score,s[s.ph_call=='silent'].signed_ph_score)
o=pd.DataFrame(rows,columns=['comparison','n_a','n_b','median_a','median_b','U','P_two_sided']); print(o.to_string())
o.to_csv('tables/b1_04_handling_mw_options.tsv',sep='\t',index=False)
g=pd.read_csv(M+'forking_paths_v2/forking_v2_c_ph_paths.tsv',sep='\t',comment='#')
print('distinct analyses (hash_id):',g.hash_id.nunique(),'| study-tier scored:',g[g.tier=='study'].groupby('hash_id').signed_ph_score.apply(lambda x:x.notna().any()).sum())
