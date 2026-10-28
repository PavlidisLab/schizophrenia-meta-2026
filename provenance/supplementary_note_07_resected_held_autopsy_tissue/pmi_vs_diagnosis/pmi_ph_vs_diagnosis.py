#!/usr/bin/env python3
"""PMI, brain pH and age against schizophrenia diagnosis, per dataset and pooled.
Data: dataset packages of out/ds (cov.tsv: dx_group, dx_text, age, sex, pmi_h, ph, region). Cases = dx_group case with a schizophrenia term in the diagnosis text; controls = dx_group control; Gemma outlier samples dropped.
A GEO series that has a tier-1 package (T1_*) enters only through it; otherwise through its Gemma package (G*_GSE...), so no series enters twice. Regions of one series are different dataset-regions of the same donors.
Per dataset-region (>= MINN cases and controls with the outcome and covariates): OLS of the outcome on diagnosis (+ covariates); estimate = case minus control on the original scale (hours; pH units; years) with its SE; standardized = estimate / pooled within-group SD (unadjusted).
Models: m0 diagnosis only; m1 + sex; m2 + sex + age (not for age); m3 (pH only) + sex + age + PMI. PMI also as log PMI.
Pooling (inverse variance, DerSimonian-Laird tau2, Hartung-Knapp-Sidik-Jonkman t interval on k - 1 df, prediction interval on k - 2 df, I2): 
  cluster: one dataset-region per cohort cluster (largest n; cluster from dataset_lineage.tsv, a dataset without a listed cluster is its own cluster)
  series: one dataset-region per GEO series (largest n)
  all: every dataset-region (regions of one series share donors: shown only to quantify the inflation).
Sets: ph14 = the 14 dataset-regions of Fig. 3a / Note 6 (tier-1 packages of GSE17612, GSE208338, GSE21138, GSE21935, GSE35974, GSE35977, GSE53987 x3, GSE80655 x3, GSE92538 x2); pmi_all = every dataset-region with PMI recorded for >= 80% of samples; pmi_14 = PMI in the 14.
Usage: pmi_ph_vs_diagnosis.py RUN_DIR DATASET_LINEAGE OUT_DIR"""
import sys, os, re, hashlib, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
RUN, LIN, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True); MINN = 5
PH14 = ["T1_GSE17612", "T1_GSE208338", "T1_GSE21138", "T1_GSE21935", "T1_GSE35974", "T1_GSE35977", "T1_GSE53987_BA46", "T1_GSE53987_hippocampus", "T1_GSE53987_striatum", "T1_GSE80655_AnCg", "T1_GSE80655_DLPFC", "T1_GSE80655_nAcc", "T1_GSE92538_GPL570", "T1_GSE92538_GPL96"]
lin = pd.read_csv(LIN, sep="\t", comment="#", dtype=str); clus = dict(zip(lin.key, lin.cluster))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
allkeys = sorted({f[:-len(".cov.tsv")] for f in os.listdir(f"{RUN}/out/ds") if f.endswith(".cov.tsv")})
gse = lambda k: (re.search(r"(GSE\d+)", k) or [None, k])[1]
t1_gse = {gse(k) for k in allkeys if k.startswith("T1_")}
keys = [k for k in allkeys if k.startswith("T1_") or gse(k) not in t1_gse]
def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; df = len(y) - X.shape[1]
    if df < 3: return None
    s2 = (r ** 2).sum() / df; cov = s2 * np.linalg.inv(X.T @ X); return b[1], np.sqrt(cov[1, 1]), df
rows = []
for key in keys:
    c = pd.read_csv(f"{RUN}/out/ds/{key}.cov.tsv", sep="\t", dtype={"sample": str})
    if "dx_group" not in c: continue
    outl = c["outlier"].astype(str).str.lower().isin(["true", "1", "1.0"]) if "outlier" in c else pd.Series(False, index=c.index)
    szm = c["dx_text"].astype(str).str.contains(r"schizo|\bSZ\b|\bSCZ\b|\bSCHIZ", case=False, regex=True, na=False)
    c = c[~outl & ((c.dx_group == "control") | ((c.dx_group == "case") & szm))].copy(); c["case"] = (c.dx_group == "case").astype(float)
    for col in ("age", "pmi_h", "ph"): c[col] = pd.to_numeric(c[col], errors="coerce") if col in c else np.nan
    c["sexb"] = c["sex"].map({"M": 1.0, "F": 0.0}) if "sex" in c else np.nan; c["logpmi"] = np.log(c.pmi_h.where(c.pmi_h > 0))
    for out, cols_models in (("pH", {"m0": [], "m1": ["sexb"], "m2": ["sexb", "age"], "m3": ["sexb", "age", "pmi_h"]}), ("PMI", {"m0": [], "m1": ["sexb"], "m2": ["sexb", "age"]}), ("logPMI", {"m0": [], "m1": ["sexb"], "m2": ["sexb", "age"]}), ("age", {"m0": [], "m1": ["sexb"]})):
        ycol = {"pH": "ph", "PMI": "pmi_h", "logPMI": "logpmi", "age": "age"}[out]
        if c[ycol].notna().mean() < 0.8: continue
        for mn, cv in cols_models.items():
            d = c.dropna(subset=[ycol] + cv)
            if (d.case == 1).sum() < MINN or (d.case == 0).sum() < MINN: continue
            X = np.column_stack([np.ones(len(d)), d.case.values] + [d[v].values for v in cv])
            if np.linalg.matrix_rank(X) < X.shape[1]: continue
            f = ols(d[ycol].values, X)
            if f is None: continue
            g1, g0 = d[d.case == 1][ycol].values, d[d.case == 0][ycol].values; sdp = np.sqrt(((len(g1) - 1) * g1.var(ddof=1) + (len(g0) - 1) * g0.var(ddof=1)) / (len(g1) + len(g0) - 2))
            rows.append(dict(key=key, series=gse(key), cluster=clus.get(key, key), outcome=out, model=mn, n_case=len(g1), n_ctl=len(g0), mean_case=g1.mean(), mean_ctl=g0.mean(), diff=f[0], se=f[1], df=f[2], sd_pooled=sdp, std_diff=f[0] / sdp, std_se=f[1] / sdp))
R = pd.DataFrame(rows); R["n"] = R.n_case + R.n_ctl
hdr = f"# pmi_ph_vs_diagnosis.py; run {RUN}; dataset_lineage.tsv {sha(LIN)}; dedupe: a GEO series with a T1 package enters only through it ({len(keys)} packages from {len(allkeys)})\n# diff = case minus control on the original scale (hours; pH units; years; logPMI in natural-log hours), se its standard error; std_diff = diff / pooled within-group SD (unadjusted); models: m0 diagnosis; m1 + sex; m2 + sex + age; m3 (pH) + sex + age + PMI\n"
with open(f"{OUT}/per_dataset_diagnosis_differences.tsv", "w") as fh: fh.write(hdr); R.round(5).to_csv(fh, sep="\t", index=False)
def meta(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); k = len(y)
    if k < 2: return dict(k=k, est=y[0] if k else np.nan)
    w = 1 / se ** 2; fe = (w * y).sum() / w.sum(); Q = (w * (y - fe) ** 2).sum(); Cc = w.sum() - (w ** 2).sum() / w.sum(); tau2 = max(0.0, (Q - (k - 1)) / Cc); ws = 1 / (se ** 2 + tau2); mu = (ws * y).sum() / ws.sum()
    q = (ws * (y - mu) ** 2).sum() / (k - 1); s = np.sqrt(max(q, 1e-12) / ws.sum()); tc = stats.t.ppf(0.975, k - 1); out = dict(k=k, est=mu, se_hk=s, lo=mu - tc * s, hi=mu + tc * s, tau=np.sqrt(tau2), I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, Q_p=stats.chi2.sf(Q, k - 1), p=2 * stats.t.sf(abs(mu / s), k - 1))
    if k > 2: tp = stats.t.ppf(0.975, k - 2); sp = np.sqrt(tau2 + s ** 2); out.update(pred_lo=mu - tp * sp, pred_hi=mu + tp * sp)
    return out
SETS = {"ph14": ("pH", PH14), "pmi_all": ("PMI", None), "pmi_14": ("PMI", PH14), "logpmi_all": ("logPMI", None), "logpmi_14": ("logPMI", PH14), "age_14": ("age", PH14), "age_all": ("age", None), "ph_all_with_ph": ("pH", None)}
P = []
for sname, (out, lst) in SETS.items():
    for mn in ("m0", "m1", "m2", "m3"):
        g = R[(R.outcome == out) & (R.model == mn)]
        if lst is not None: g = g[g.key.isin(lst)]
        if len(g) < 2: continue
        for rule in ("cluster", "series", "all"):
            one = g if rule == "all" else g.sort_values("n", ascending=False).drop_duplicates("cluster" if rule == "cluster" else "series")
            for scale, (yc, sc) in (("original", ("diff", "se")), ("standardized", ("std_diff", "std_se"))):
                m = meta(one[yc], one[sc]); m.update(set=sname, outcome=out, model=mn, rule=rule, scale=scale, n_people=int(one.n.sum()), n_case=int(one.n_case.sum()), n_ctl=int(one.n_ctl.sum())); P.append(m)
P = pd.DataFrame(P)[["set", "outcome", "model", "rule", "scale", "k", "n_people", "n_case", "n_ctl", "est", "lo", "hi", "pred_lo", "pred_hi", "tau", "I2", "Q_p", "p"]].round(4)
with open(f"{OUT}/pooled_diagnosis_differences.tsv", "w") as fh: fh.write(hdr + "# rules: cluster = one dataset-region per cohort cluster; series = one per GEO series; all = every dataset-region (shared donors, inflated). HKSJ interval, prediction interval on k - 2 df\n"); P.to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 250)
print(P[(P.rule != "all") & (P.model.isin(["m1", "m2"]))].to_string()); print(R.groupby("outcome").key.nunique())
