#!/usr/bin/env python3
"""pH by diagnosis and manner of death in GSE208338 (DLPFC) and GSE202537, from results/geo_sample_characteristics.tsv
(GEO series-matrix characteristics, sha256 in results/geo_files.tsv). One row per subject (GSE208338: first sample per 'sample' id;
GSE202537: first row per patients_id). Writes results/ph_by_manner_of_death.tsv and results/ph_natural_death_only.tsv.
Run: python3 analyses/hagihara_ph_meta_replication/scripts/ph_by_manner_of_death.py"""
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.formula.api as smf
R = Path(__file__).resolve().parents[1] / "results"
d = pd.read_csv(R / "geo_sample_characteristics.tsv", sep="\t", low_memory=False)
def wide(g):
    return d[d.series == g].pivot_table(index="gsm", columns="key", values="value", aggfunc="first")
def num(x, cols):
    for c in cols: x[c] = pd.to_numeric(x[c], errors="coerce")
    return x
a = num(wide("GSE208338"), ["ph", "age", "pmi"])
a = a[a.tissue.str.contains("DLPFC")].drop_duplicates("sample")
a["group"] = a.diagnosis.map({"SCZ": "SZ", "CTL": "control"}).fillna("other")
a["manner"] = a["type of death"]
b = num(wide("GSE202537"), ["ph", "age", "pmi"]).drop_duplicates("patients_id")
b["group"] = np.where(b["disease state"].str.contains("schiz"), "SZ", np.where(b["disease state"].str.contains("control"), "control", "other"))
b["manner"] = b["manner of death"]
rows = []
for g, x in (("GSE208338", a), ("GSE202537", b)):
    t = x.groupby(["group", "manner"]).agg(n=("ph", "size"), mean_ph=("ph", "mean"), sd_ph=("ph", "std"), mean_age=("age", "mean"), mean_pmi=("pmi", "mean")).round(3).reset_index()
    t.insert(0, "dataset", g); rows.append(t)
pd.concat(rows).to_csv(R / "ph_by_manner_of_death.tsv", sep="\t", index=False)
out = []
for g, x, sex in (("GSE208338", a, "sex"), ("GSE202537", b, "gender")):
    n = x[(x.manner == "Natural") & x.group.isin(["SZ", "control"])].copy(); n["sz"] = (n.group == "SZ").astype(int)
    m = smf.ols(f"ph ~ sz + age + C({sex})", n).fit()
    out.append(dict(dataset=g, model="natural deaths only; age + sex", n_sz=int(n.sz.sum()), n_ctl=int((1 - n.sz).sum()), coef_sz=m.params["sz"], se=m.bse["sz"], p=m.pvalues["sz"]))
pd.DataFrame(out).to_csv(R / "ph_natural_death_only.tsv", sep="\t", index=False)
print(pd.DataFrame(out).round(4).to_string())
