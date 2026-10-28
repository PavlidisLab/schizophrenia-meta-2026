#!/usr/bin/env python3
"""Share of brain-bank samples with PMI under 20 h and under 36 h, from heatmap_matrix.npz of pmi_rank_genes.py (50 datasets, samples with PMI recorded, controls and cases; outliers dropped).
Usage: pmi_shares.py DIR"""
import sys, hashlib, numpy as np, pandas as pd
D = sys.argv[1]; z = np.load(f"{D}/heatmap_matrix.npz", allow_pickle=True); h = z["hours"].astype(float); dx = z["dx"].astype(str); ds = z["dataset"].astype(str)
import re
keys_all = sorted(set(ds)); gse = lambda k: (re.search(r"(GSE\d+)", k) or [None, k])[1]; t1 = {gse(k) for k in keys_all if k.startswith("T1_")}
drop = {k for k in keys_all if not k.startswith("T1_") and gse(k) in t1}   # Gemma package of a series that also has a tier-1 package: same donors
keep = ~np.isin(ds, list(drop))
rows = []
for lab, m in (("all samples", keep), ("controls only", keep & (dx == "control")), ("cases only", keep & (dx == "case")), ("all samples, duplicate series kept (earlier version)", np.ones(len(h), bool))):
    x = h[m]; rows.append(dict(set=lab, datasets=len(set(ds[m])), samples=int(m.sum()), share_lt20=(x < 20).mean(), n_lt20=int((x < 20).sum()), share_lt36=(x < 36).mean(), n_lt36=int((x < 36).sum()), median_h=np.median(x), q25_h=np.percentile(x, 25), q75_h=np.percentile(x, 75), min_h=x.min(), max_h=x.max()))
R = pd.DataFrame(rows).round(4); sha = hashlib.sha256(open(f"{D}/heatmap_matrix.npz", "rb").read()).hexdigest()[:12]
with open(f"{D}/pmi_shares.tsv", "w") as fh: fh.write(f"# pmi_shares.py; heatmap_matrix.npz sha256 {sha}; PMI in hours from the sample covariates (time to freezing as recorded by each dataset)\n"); R.to_csv(fh, sep="\t", index=False)
print(sorted(drop)); print(R.to_string())
