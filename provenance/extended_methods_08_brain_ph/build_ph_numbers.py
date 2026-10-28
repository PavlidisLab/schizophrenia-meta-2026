#!/usr/bin/env python3
"""(1) tidy per-dataset-region table behind '90% of cases within the control pH range' and (2) Pearson correlations of the schizophrenia effect with and without pH for the genes in 10 or more lists (RA, 2026-10-19)."""
import pandas as pd, numpy as np, hashlib, os, datetime, shutil
from scipy import stats
H=os.path.expanduser("~")
PR=f"{H}/Dev/srep-mauney-integrated/analyses/supplement_number_tasks/taskE_ph_overlap_per_region.tsv"
PP=f"{H}/Dev/srep-mauney-integrated/analyses/supplement_number_tasks/taskE_ph_overlap.py"
GS=f"{H}/fig07/fig07_a2_a4_gene_scatter.tsv"
sha=lambda p:hashlib.sha256(open(p,"rb").read()).hexdigest()
# (1) copy of the per-region table with a header naming the claim
raw=open(PR).read().split("\n"); body=[l for l in raw if not l.startswith("#")]; hdr=[l for l in raw if l.startswith("#")]
d=pd.read_csv(PR,sep="\t",comment="#")
tot=pd.DataFrame([dict(dataset_region="all 14 dataset-regions pooled (counts summed)",study="9 studies",n_case=d.n_case.sum(),n_control=d.n_control.sum(),n_case_in_control_range=d.n_case_in_control_range.sum(),share_case_in_control_range=d.n_case_in_control_range.sum()/d.n_case.sum(),n_case_in_control_iqr=d.n_case_in_control_iqr.sum(),share_case_in_control_iqr=d.n_case_in_control_iqr.sum()/d.n_case.sum())])
out=pd.concat([d,tot],ignore_index=True)
h1=["# claim: Supplementary Note 12, 'pH ranges overlap heavily between cases and controls (90% of cases fall within the pH range of the controls in the same dataset)'; the figure is the per-dataset-region pH distribution of the pH panel (fig07 a1, Fig. 4a in the manuscript numbering), not Fig. 3a",
"# script: note12_ph_overlap_and_ph_adjustment_2026-10-19/build_ph_numbers.py copies the per-region table made by analyses/supplement_number_tasks/taskE_ph_overlap.py (run 2026-10-12) and adds a pooled row",
f"# source: analyses/supplement_number_tasks/taskE_ph_overlap_per_region.tsv sha256 {sha(PR)[:12]}; taskE_ph_overlap.py sha256 {sha(PP)[:12]}; the original header follows",
"# pooled row = counts summed over the 14 dataset-regions (GSE208338 shares donors with GSE21138, so the pooled counts are not donor-independent); a case is 'within' if its pH lies inside [min, max] of the controls of the same dataset-region, inclusive; n_case_in_control_iqr uses the control quartiles"]+hdr
with open("ph_case_vs_control_range_per_dataset_region.tsv","w") as f:
    f.write("\n".join(h1)+"\n"); out.to_csv(f,sep="\t",index=False)
print(out.iloc[-1][["n_case","n_case_in_control_range","share_case_in_control_range"]].to_dict(), "datasets:",len(d))
# (2) correlations with and without pH
g=pd.read_csv(GS,sep="\t",comment="#")
def corr(x,y,df):
    s=df[[x,y]].dropna(); r,p=stats.pearsonr(s[x],s[y]); n=len(s); z=np.arctanh(r); se=1/np.sqrt(n-3)
    rho=stats.spearmanr(s[x],s[y])[0]
    return dict(x=x,y=y,n=n,pearson_r=round(r,4),ci95_lo=round(np.tanh(z-1.96*se),4),ci95_hi=round(np.tanh(z+1.96*se),4),p=p,spearman_rho=round(rho,4))
rows=[]
sets={"genes in >= 10 lists with a combined pH z":g[(g.in_ge10=="yes")&g.z_ph_controls.notna()],
      "the 383 recurrent genes with a combined pH z":g[g.z_ph_controls.notna()&g.z_sz_raw.notna()]}
for name,df in sets.items():
    df=df[df.z_ph_controls.notna()]
    for x,y in (("z_ph_controls","z_sz_raw"),("z_ph_controls","z_sz_ph_adjusted"),("z_sz_raw","z_sz_ph_adjusted")):
        rows.append(dict(set=name,**corr(x,y,df)))
res=pd.DataFrame(rows)
h2=["# claim: Extended Methods 8, 'Schizophrenia effect with and without pH' (same five donor-separate datasets with pH; schizophrenia effect with expression ~ diagnosis + age + sex, then with pH added on the same samples, Stouffer z; the comparison covers the genes in 10 or more lists that have a combined pH z); this table gives the Pearson correlations the paragraph says are reported",
"# script: note12_ph_overlap_and_ph_adjustment_2026-10-19/build_ph_numbers.py (reads the delivered gene table; no model fits)",
f"# source: figure_data_core80_2026-10-17/fig07/fig07_a2_a4_gene_scatter.tsv sha256 {sha(GS)[:12]} (columns z_ph_controls = combined control pH z over the donor-separate datasets, z_sz_raw = Stouffer z of the schizophrenia effect without pH, z_sz_ph_adjusted = with pH added; in_ge10 = gene in >= 10 lists)",
"# rows: pearson_r with a Fisher-z 95% interval treating genes as independent (genes are not independent; the interval is a guide), p two-sided, spearman_rho for reference; x and y name the two columns correlated: (pH z, effect without pH), (pH z, effect with pH), (effect without pH, effect with pH)",
f"# run {datetime.date.today()}"]
with open("sz_effect_with_and_without_ph_correlations.tsv","w") as f:
    f.write("\n".join(h2)+"\n"); res.to_csv(f,sep="\t",index=False)
print(res.to_string())
