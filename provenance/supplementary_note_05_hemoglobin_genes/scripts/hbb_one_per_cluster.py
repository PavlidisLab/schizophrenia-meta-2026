# hbb_one_per_cluster.py  (2026-10-19) -- HBB SZ-minus-control log2FC, DerSimonian-Laird, one dataset per cohort cluster (largest n), as 03_summarize.py does for the blood scores.
# input: analyses/hemoglobin_2026-10-11/samplelevel_hb_sz_effect_by_adjustment.tsv (sha256 96c9c6beec14), analyses/controls_age_sex/data/datasets.tsv (cluster). No randomness. Base: 22 bulk sample-level datasets (no study base or axis enters).
import pandas as pd, numpy as np, math, sys, os
R=os.path.expanduser("schizophrenia-transcriptome-replication/analyses")
ad=pd.read_csv(f"{R}/hemoglobin_2026-10-11/samplelevel_hb_sz_effect_by_adjustment.tsv",sep="\t",comment="#")
reg=pd.read_csv(f"{R}/controls_age_sex/data/datasets.tsv",sep="\t",comment="#")[["dataset","cluster"]]
ad=ad[(ad.compartment=="bulk")&(ad.gene=="HBB")&(ad.model=="base")].merge(reg,on="dataset"); ad["se"]=(ad.logfc/ad.t).abs(); ad["n"]=ad.n_sz+ad.n_ctl
def dl(est,se):
    est,se=np.asarray(est,float),np.asarray(se,float); w=1/se**2; fe=(w*est).sum()/w.sum(); Q=(w*(est-fe)**2).sum(); k=len(est)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())) if k>1 else 0
    ws=1/(se**2+tau2); m=(ws*est).sum()/ws.sum(); s=math.sqrt(1/ws.sum()); return m,m-1.96*s,m+1.96*s,k
print("all 22",dl(ad.logfc,ad.se))
sub=ad.sort_values("n",ascending=False).drop_duplicates("cluster"); print("one per cluster",dl(sub.logfc,sub.se)); print(sub[["dataset","cluster","n","logfc"]].to_string())
