#!/usr/bin/env python3
"""GSE92538 (Pritzker; Ramaker 2017 / Sutton 2018): diagnosis effect on tissue pH with and without the agonal factor in the model.

    python3 analyses/hagihara_ph_meta_replication/scripts/gse92538_ph_agonal.py       (run from the tree root; reads results/geo_sample_characteristics.tsv)

Input: results/geo_sample_characteristics.tsv (from geo_characteristics.py; GEO series matrices GSE92538-GPL10526 and -GPL17027, sha256 in results/geo_files.tsv).
One row per subject ('subject id'; the series carry technical replicates), pH = 'tissue ph (cerebellum)', agonal = 'agonal factor' (0-4 score), PMI = 'post-mortem interval' (h), age (y).
OLS (statsmodels), HC3 errors: pH ~ C(dx, Treatment('Control')) + age + sex [+ PMI] [+ agonal]. Writes results/gse92538_ph_models.tsv and prints the table."""
import os, sys, csv, collections
import numpy as np, pandas as pd, statsmodels.formula.api as smf
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "..", "results")
L = pd.read_csv(os.path.join(RES, "geo_sample_characteristics.tsv"), sep="\t", dtype=str)
L = L[L.series == "GSE92538"]
W = L.pivot_table(index=["file", "gsm"], columns="key", values="value", aggfunc="first").reset_index()
def num(s): return pd.to_numeric(s.replace({"NA": np.nan, "N/A": np.nan}), errors="coerce")
D = pd.DataFrame({"file": W.file, "gsm": W.gsm, "subject": W["subject id"], "dx": W["diagnosis"], "age": num(W["age"]), "sex": W["gender"], "pmi": num(W["post-mortem interval"]),
                  "ph": num(W["tissue ph (cerebellum)"]), "agonal": num(W["agonal factor"]), "suicide": num(W["suicide (1=yes)"])})
D["platform"] = D.file.str.extract(r"(GPL\d+)")
rows = []
def fit(label, d, formula):
    m = smf.ols(formula, d).fit(cov_type="HC3")
    for k in m.params.index:
        if "dx" in k:
            rows.append({"subset": label, "formula": formula, "term": k.split("T.")[-1].rstrip("]"), "n": int(m.nobs), "coef": m.params[k], "se": m.bse[k], "p": m.pvalues[k]})
    return m
print("subjects with duplicate samples per platform:", D.groupby("platform").subject.apply(lambda s: s.duplicated().sum()).to_dict())
for name, d in [("GPL10526 (Schiz cohorts 1,2; Dep 5,6)", D[D.platform == "GPL10526"]), ("GPL17027", D[D.platform == "GPL17027"]), ("both, one row per subject (GPL17027 preferred)", D.sort_values("platform", ascending=False).drop_duplicates("subject"))]:
    d = d.drop_duplicates("subject") if name.startswith("GPL") else d
    base = d.dropna(subset=["ph", "age", "pmi", "agonal"])
    print(f"\n{name}: {len(d)} subjects; with pH, age, PMI and agonal: {len(base)}; agonal score counts {base.agonal.value_counts().sort_index().to_dict()}; by dx: {base.dx.value_counts().to_dict()}")
    print("  agonal>0 by dx:", base.assign(a=base.agonal > 0).groupby("dx").a.agg(["sum", "count"]).to_dict("index"))
    fit(name, base, "ph ~ C(dx, Treatment('Control')) + age + sex")
    fit(name, base, "ph ~ C(dx, Treatment('Control')) + age + sex + pmi")
    fit(name, base, "ph ~ C(dx, Treatment('Control')) + age + sex + pmi + agonal")
    fit(name, base, "ph ~ C(dx, Treatment('Control')) + age + sex + pmi + I(agonal > 0)")
    m = smf.ols("ph ~ agonal + age + sex + pmi", base).fit(cov_type="HC3"); rows.append({"subset": name, "formula": "ph ~ agonal + age + sex + pmi", "term": "agonal", "n": int(m.nobs), "coef": m.params["agonal"], "se": m.bse["agonal"], "p": m.pvalues["agonal"]})
R = pd.DataFrame(rows); R.to_csv(os.path.join(RES, "gse92538_ph_models.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 70)
print(R.assign(coef=R.coef.round(3), se=R.se.round(3), p=R.p.round(4))[["subset", "formula", "term", "n", "coef", "se", "p"]].to_string(index=False))
