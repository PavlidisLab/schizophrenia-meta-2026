#!/usr/bin/env python3
"""Population A: genes per dataset at each cutoff, for the 28 registry datasets of analyses/sample_level_reanalysis_of_datasets/results/b_de/.

Per dataset and covariate level (L0 SZ only; L1 + age + sex; L2 + PMI; L3 + pH, RIN, agonal; L4 + medication where recorded;
L5a, L5b with surrogate variables): the number of genes (rows of the b_de table with a p value) at
  BH q < 0.05, BH q < 0.1 (BH over the genes of that dataset and level), nominal p < 0.05, p < 0.01, p < 0.001.
The b_de p values are limma's (b_covariate_series.R). BH is analyses/gwas_vs_de/scripts/stats.py::bh. Check: q < 0.05 counts are
compared with nq05_* in fairness_sample_level/results/b_levels.tsv (written by the same R script).
python3 a_registry_counts.py  ->  results/a_registry_per_dataset_level.tsv, results/a_registry_summary.tsv
"""
import gzip, os
import numpy as np, pandas as pd
from common import AN, HERE, bh, write_tsv

FS = os.path.join(AN, "fairness_sample_level", "results")
LEVELS = ["L0", "L1", "L2", "L3", "L4", "L5a", "L5b"]
CUTS = [("q05", "BH q < 0.05"), ("q10", "BH q < 0.1"), ("p05", "p < 0.05"), ("p01", "p < 0.01"), ("p001", "p < 0.001")]
levels = pd.read_csv(os.path.join(FS, "b_levels.tsv"), sep="\t").set_index("dataset")
names = sorted(f[:-7] for f in os.listdir(os.path.join(FS, "b_de")) if f.endswith(".tsv.gz"))
assert len(names) == 28 and set(names) == set(levels.index), (len(names), set(names) ^ set(levels.index))

rows = []
for ds in names:
    t = pd.read_csv(os.path.join(FS, "b_de", f"{ds}.tsv.gz"), sep="\t", usecols=["gene"] + [f"p_{l}" for l in LEVELS])
    for lv in LEVELS:
        p = t[f"p_{lv}"].to_numpy(float)
        ok = ~np.isnan(p)
        if ok.sum() == 0:
            continue
        q = bh(p)
        r = dict(dataset=ds, level=lv, n_sz=int(levels.loc[ds, "n_sz"]), n_ctl=int(levels.loc[ds, "n_ctl"]), n_genes_tested=int(ok.sum()),
                 genes_q05=int((q < 0.05).sum()), genes_q10=int((q < 0.1).sum()), genes_p05=int((p < 0.05).sum()),
                 genes_p01=int((p < 0.01).sum()), genes_p001=int((p < 0.001).sum()))
        col = f"nq05_{lv}"
        r["nq05_b_levels"] = int(levels.loc[ds, col]) if col in levels.columns and not pd.isna(levels.loc[ds, col]) else np.nan
        rows.append(r)
D = pd.DataFrame(rows)
D["n_samples"] = D.n_sz + D.n_ctl
D["fewer_than_20_samples"] = D.n_samples < 20
mism = D[(D.nq05_b_levels.notna()) & (D.nq05_b_levels != D.genes_q05)]
print("datasets x levels:", len(D), "; q05 differs from b_levels nq05 in", len(mism), "rows")
if len(mism):
    print(mism[["dataset", "level", "genes_q05", "nq05_b_levels"]].to_string())
src = [("b_levels.tsv", os.path.join(FS, "b_levels.tsv")), (f"b_de/*.tsv.gz ({len(names)} files; one sha256 over the sorted 'name sha256' lines)", [os.path.join(FS, "b_de", f"{n}.tsv.gz") for n in names])]
write_tsv(D, os.path.join(HERE, "results", "a_registry_per_dataset_level.tsv"),
          "genes per registry dataset at BH q < 0.05, q < 0.1 and nominal p < 0.05, 0.01, 0.001, per covariate level (28 datasets, limma SZ vs control)",
          "analyses/studies_without_gene_lists/a_registry_counts.py", src,
          extra=[f"q05 vs b_levels nq05: {len(D) - len(mism)} of {len(D)} rows equal"])

S = []
for lv in LEVELS:
    d = D[D.level == lv]
    for sub, dd in (("all datasets", d), ("excluding datasets with < 20 samples", d[~d.fewer_than_20_samples])):
        for key, lab in CUTS:
            c = dd[f"genes_{key}"]
            S.append(dict(level=lv, population=sub, cutoff=lab, n_datasets=len(dd), n_datasets_zero_genes=int((c == 0).sum()),
                          median_genes=float(c.median()), min_genes=int(c.min()), max_genes=int(c.max())))
S = pd.DataFrame(S)
write_tsv(S, os.path.join(HERE, "results", "a_registry_summary.tsv"),
          "datasets with zero genes, and median and range of genes per dataset, by cutoff and level (28 registry datasets)",
          "analyses/studies_without_gene_lists/a_registry_counts.py", src)
print(S[S.population == "all datasets"].pivot(index="cutoff", columns="level", values="n_datasets_zero_genes"))
