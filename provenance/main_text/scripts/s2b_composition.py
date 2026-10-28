#!/usr/bin/env python3
"""Step 2b (deterministic). Case/control balance of non-sudden death per dataset, per collection and over unique donors; covariates by death speed.
Odds ratios: Fisher exact (conditional MLE) with its confidence interval; the proportion difference carries a Wilson-based Newcombe 95% interval. Counts are small (1 to 12 non-sudden per group)."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sd_common import *
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from scipy.stats.contingency import odds_ratio
import pickle
P = pickle.load(open(f"{LOCAL}/prep.pkl", "rb")); don = P["don"]
INPUTS = [f"{LOCAL}/donors_LOCAL_ONLY.tsv"]
def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)); return ((c - h) / d, (c + h) / d)
def newcombe(k1, n1, k0, n0):
    l1, u1 = wilson(k1, n1); l0, u0 = wilson(k0, n0); p1, p0 = k1 / n1, k0 / n0
    return p1 - p0, p1 - p0 - np.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2), p1 - p0 + np.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2)
def block(dd, label, kind):
    rows = []
    ctrl = dd[dd.dx == "CTRL"]
    for grp, sel in (("SZ", dd.dx == "SZ"), ("BD", dd.dx == "BD"), ("MDD", dd.dx == "MDD"), ("all cases", dd.dx != "CTRL")):
        g = dd[sel]
        if len(g) == 0: continue
        k1, n1, k0, n0 = int(g.nonsudden.sum()), len(g), int(ctrl.nonsudden.sum()), len(ctrl)
        row = dict(scope=label, scope_type=kind, group=grp, n_group=n1, nonsudden_group=k1, share_group=k1 / n1, n_control=n0, nonsudden_control=k0, share_control=k0 / n0 if n0 else np.nan)
        if n0:
            d, lo, hi = newcombe(k1, n1, k0, n0); row.update(diff=d, diff_lo=lo, diff_hi=hi)
            tab = [[k1, n1 - k1], [k0, n0 - k0]]; row["fisher_P"] = fisher_exact(tab)[1]
            try:
                r = odds_ratio(tab, kind="conditional"); ci = r.confidence_interval(0.95); row.update(odds_ratio=r.statistic, or_lo=ci.low, or_hi=ci.high)
            except Exception: pass
        rows.append(row)
    return rows
R = []
for ds in DATASETS:
    R += block(don.loc[[c for c in P["Z"][ds].columns]], ds, "dataset")
for c, nm in (("A", "Array Collection (unique donors)"), ("C", "Neuropathology Consortium (unique donors)")):
    R += block(don[don.collection == c], nm, "collection")
R += block(don, "all unique donors", "all")
write_tsv(pd.DataFrame(R), "nonsudden_share_by_diagnosis.tsv", "s2b_composition.py", INPUTS, note="Share of donors with non-sudden death (Slower death, Possible anoxia, On respirator) by diagnosis versus control, per dataset, per collection and over the 165 unique donors. diff = share(group) - share(control) with a Newcombe 95% interval; odds ratio = Fisher conditional MLE with exact interval; Fisher P reported, no threshold applied.")
C = []
for nm, g in (("sudden", don[don.nonsudden == 0]), ("non-sudden", don[don.nonsudden == 1])) + tuple((f"non-sudden: {s}", don[don.speed == s]) for s in ["Slower death", "Possible anoxia", "On respirator"]):
    C.append(dict(group=nm, n=len(g), mean_age=g.age.mean(), male_share=g.male.mean(), mean_pH=g.ph.mean(), median_PMI=g.pmi.median(), mean_PMI=g.pmi.mean(), SZ=int((g.dx == "SZ").sum()), BD=int((g.dx == "BD").sum()), MDD=int((g.dx == "MDD").sum()), CTRL=int((g.dx == "CTRL").sum())))
for col in ("ph", "pmi", "age"):
    a, b = don[col][don.nonsudden == 1], don[col][don.nonsudden == 0]
    C.append(dict(group=f"test {col}: non-sudden vs sudden", n=len(don), mean_age=np.nan, male_share=np.nan, mean_pH=a.mean() - b.mean() if col == "ph" else np.nan, median_PMI=np.nan, mean_PMI=np.nan,
                  SZ=np.nan, BD=np.nan, MDD=np.nan, CTRL=mannwhitneyu(a, b).pvalue))
write_tsv(pd.DataFrame(C), "covariates_by_death_speed.tsv", "s2b_composition.py", INPUTS, note="Unique donors (n=165). Rows 'test ...': Mann-Whitney P in the CTRL column, mean difference in mean_pH for pH. Not a model; context for the adjusted fits.")
print(pd.DataFrame(R)[["scope","group","n_group","nonsudden_group","share_group","n_control","nonsudden_control","share_control","odds_ratio","or_lo","or_hi"]].round(3).to_string())
print(pd.DataFrame(C).round(3).to_string())
