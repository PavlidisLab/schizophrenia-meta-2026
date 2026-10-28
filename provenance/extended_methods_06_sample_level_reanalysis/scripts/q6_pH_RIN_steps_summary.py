"""Check 6: covariate series with pH alone and RIN alone added after age + sex + PMI (L2). Levels: L2p = L2 + pH; L2r = L2 + RIN; L3 = L2 + pH, RIN, agonal factor.
Inputs: b_series_pH_RIN_steps/b_pairs.tsv (cd analyses/sample_level_reanalysis_of_datasets; SERIES_LEVELS=L0,L1,L2,L2p,L2r,L3,L4 SERIES_OUT=../manuscript_number_checks/b_series_pH_RIN_steps NCORES=3 Rscript scripts/b_covariate_series.R),
        fairness_sample_level/results/b_pairs.tsv (identity check of the shared levels).
Pair sets: existing_192 = subject-separate pairs of two bulk datasets (Fig. 4d); a step missing for either dataset falls back to L2 for the medians in this set (the L1p rule of q5)
and is NA in the long export; pH_both / RIN_both = pairs where both datasets record the covariate; pH_and_RIN_both = the intersection (the three steps on identical pairs).
Bootstrap: 2,000 resamples, seed 20261013; pair bootstrap over pairs; dataset bootstrap redraws datasets, weights a pair by count_a*count_b, weighted median.
Run: python3 -I q6_pH_RIN_steps_summary.py"""
import os, hashlib, datetime, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); FS = os.path.join(os.path.dirname(HERE), "fairness_sample_level")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
NEW = os.path.join(HERE, "b_series_pH_RIN_steps/b_pairs.tsv"); OLD = os.path.join(FS, "results/b_pairs.tsv"); NEWL = os.path.join(HERE, "b_series_pH_RIN_steps/b_levels.tsv")
FIGDIR = os.path.expanduser("figure_data/fig07")
rd = lambda p: pd.read_csv(p, sep="\t", comment="#")
new, old, lev = rd(NEW), rd(OLD), rd(NEWL)
LV = ["L0", "L1", "L2", "L2p", "L2r", "L3", "L4"]
DESC = {"L0": "no covariates", "L1": "age + sex", "L2": "age + sex + PMI", "L2p": "age + sex + PMI + pH (datasets that record pH)",
        "L2r": "age + sex + PMI + RIN (datasets that record RIN)", "L3": "age + sex + PMI + pH, RIN, agonal factor (where recorded)", "L4": "L3 + medication (GSE21138, GSE138082 only)"}
m = new[new.level.isin(["L0", "L1", "L2", "L3", "L4"])].merge(old, on=["dataset_a", "dataset_b", "level"], suffixes=("_new", "_old"))
chk = {c: float(np.nanmax(np.abs(m[f"{c}_new"] - m[f"{c}_old"]))) for c in ("rho_t", "rho_logfc", "n_common")}
chk["rows_compared"] = len(m); print("identity:", chk)
# L1p of the earlier run must equal... (not rerun here); L2p/L2r independent of it
W = new.set_index(["dataset_a", "dataset_b", "relation", "both_bulk", "level"]).rho_t.unstack("level").reset_index()
W["both_bulk"] = W.both_bulk.astype(str).str.upper() == "TRUE"
NC = new[new.level == "L0"].set_index(["dataset_a", "dataset_b"]).n_common
E = W[(W.relation == "subject-separate") & W.both_bulk].copy()
E["n_common_L0"] = [NC[(a, b)] for a, b in zip(E.dataset_a, E.dataset_b)]
E["has_pH_both"] = E.L2p.notna(); E["has_RIN_both"] = E.L2r.notna(); E["has_both_covs"] = E.has_pH_both & E.has_RIN_both
E["L2p_f"] = E.L2p.fillna(E.L2); E["L2r_f"] = E.L2r.fillna(E.L2)
print("pairs", len(E), "pH both", int(E.has_pH_both.sum()), "RIN both", int(E.has_RIN_both.sum()), "pH and RIN both", int(E.has_both_covs.sum()))
rng = np.random.default_rng(20261013); NB = 2000
def wmed(v, w):
    o = np.argsort(v); v, w = v[o], w[o]; c = np.cumsum(w); return v[np.searchsorted(c, 0.5 * c[-1])]
def boot(df, col):
    v = df[col].values; n = len(v); out = {"pair": [], "dataset": []}
    ds = sorted(set(df.dataset_a) | set(df.dataset_b)); ix = {d: i for i, d in enumerate(ds)}
    ia = df.dataset_a.map(ix).values; ib = df.dataset_b.map(ix).values
    for _ in range(NB):
        out["pair"].append(np.median(v[rng.integers(0, n, n)]))
        cnt = np.bincount(rng.integers(0, len(ds), len(ds)), minlength=len(ds)); w = (cnt[ia] * cnt[ib]).astype(float)
        if w.sum() > 0: out["dataset"].append(wmed(v, w))
    return {k: (np.percentile(a, 2.5), np.percentile(a, 97.5)) for k, a in out.items()}
SETS = [("existing_192_missing_step_falls_back_to_L2", E, True), ("pH_recorded_in_both", E[E.has_pH_both], False),
        ("RIN_recorded_in_both", E[E.has_RIN_both], False), ("pH_and_RIN_recorded_in_both", E[E.has_both_covs], False)]
col = lambda lv, fb: {"L2p": "L2p_f", "L2r": "L2r_f"}.get(lv, lv) if fb else lv
rows, drows = [], []
for name, df, fb in SETS:
    nd = len(set(df.dataset_a) | set(df.dataset_b))
    for lv in LV:
        if df[col(lv, fb)].isna().any(): continue  # step not recorded for every pair of this set
        b = boot(df, col(lv, fb)); v = df[col(lv, fb)]
        rows.append(dict(pair_set=name, level=lv, description=DESC[lv], n_pairs=len(df), n_datasets=nd, median_rho_t=float(v.median()),
            ci95_pair_bootstrap_lo=b["pair"][0], ci95_pair_bootstrap_hi=b["pair"][1], ci95_dataset_bootstrap_lo=b["dataset"][0], ci95_dataset_bootstrap_hi=b["dataset"][1]))
    for lv in ("L2p", "L2r", "L3", "L4"):
        if df[col(lv, fb)].isna().any(): continue
        d = df.assign(diff=df[col(lv, fb)] - df["L2"]); bb = boot(d, "diff")
        drows.append(dict(pair_set=name, contrast=f"{lv} minus L2", n_pairs=len(d), n_datasets=nd, median_of_differences=float(d["diff"].median()),
            ci95_pair_bootstrap_lo=bb["pair"][0], ci95_pair_bootstrap_hi=bb["pair"][1], ci95_dataset_bootstrap_lo=bb["dataset"][0], ci95_dataset_bootstrap_hi=bb["dataset"][1], share_pairs_lower=float((d["diff"] < 0).mean())))
    if name == "pH_and_RIN_recorded_in_both":
        for a, b_ in (("L2p", "L2r"), ("L2p", "L3"), ("L2r", "L3")):
            d = df.assign(diff=df[b_] - df[a]); bb = boot(d, "diff")
            drows.append(dict(pair_set=name, contrast=f"{b_} minus {a}", n_pairs=len(d), n_datasets=nd, median_of_differences=float(d["diff"].median()),
                ci95_pair_bootstrap_lo=bb["pair"][0], ci95_pair_bootstrap_hi=bb["pair"][1], ci95_dataset_bootstrap_lo=bb["dataset"][0], ci95_dataset_bootstrap_hi=bb["dataset"][1], share_pairs_lower=float((d["diff"] < 0).mean())))
S, Dd = pd.DataFrame(rows), pd.DataFrame(drows)
prov = [f"# script: analyses/manuscript_number_checks/q6_pH_RIN_steps_summary.py, run {datetime.date.today().isoformat()}; seed 20261013; 2,000 bootstrap resamples",
        f"# input: analyses/manuscript_number_checks/b_series_pH_RIN_steps/b_pairs.tsv sha256 {sha(NEW)} (analyses/sample_level_reanalysis_of_datasets/scripts/b_covariate_series.R with SERIES_LEVELS=L0,L1,L2,L2p,L2r,L3,L4); analyses/sample_level_reanalysis_of_datasets/results/b_pairs.tsv sha256 {sha(OLD)}",
        f"# shared levels L0,L1,L2,L3,L4 against the existing results/b_pairs.tsv: {chk['rows_compared']} rows compared; maximum absolute difference rho_t {chk['rho_t']:.2g}, rho_logfc {chk['rho_logfc']:.2g}, n_common {chk['n_common']:.0f}",
        "# 192 pairs = subject-separate pairs of two bulk datasets (Fig. 4d); pair bootstrap resamples pairs, dataset bootstrap redraws datasets with replacement and weights a pair by the product of its two datasets' counts (the wider interval); populations: the 28 datasets with sample-level data, old corpus status as in fig07/STATUS_2026-10-12.md"]
def wr(path, df):
    with open(path, "w") as fh:
        fh.write("\n".join(prov) + "\n"); df.to_csv(fh, sep="\t", index=False, float_format="%.5g", na_rep="NA")
wr(os.path.join(HERE, "q6_pH_RIN_medians.tsv"), S); wr(os.path.join(HERE, "q6_pH_RIN_differences.tsv"), Dd)
wr(os.path.join(FIGDIR, "fig07_a_covariate_series_pH_RIN_steps_summary.tsv"), pd.concat([S.assign(table="medians"), Dd.assign(table="paired_change_from_L2")]).set_index("table").reset_index())
ex = []
for r in E.itertuples():
    for lv in LV:
        ex.append(dict(pair_id=f"{r.dataset_a}_{r.dataset_b}", model_step=lv, rho=getattr(r, lv), model_description=DESC[lv], relation=r.relation, both_bulk="TRUE", n_common=int(r.n_common_L0),
                       pH_recorded_in_both=str(r.has_pH_both).upper(), RIN_recorded_in_both=str(r.has_RIN_both).upper(), pH_and_RIN_recorded_in_both=str(r.has_both_covs).upper()))
exprov = ["# figure: fig07 panel a (Fig. 4d), covariate series with pH alone (L2p) and RIN alone (L2r) added after age + sex + PMI (L2); L3 = L2 + pH, RIN, agonal factor as in the existing series",
          prov[0], prov[1], prov[2],
          "# rho = Spearman correlation of the SZ vs control limma t statistics of the two datasets over shared genes (rho_t); L2p / L2r are NA for a pair in which either dataset does not record pH / RIN (covered for >= 80% of samples); the three *_recorded_in_both columns flag the pairs with a value",
          "# population: 192 subject-separate pairs of two bulk datasets from the 28 datasets with sample-level data (old corpus status); the existing fig07 files are not changed"]
with open(os.path.join(FIGDIR, "fig07_a_covariate_series_pH_RIN_steps.tsv"), "w") as fh:
    fh.write("\n".join(exprov) + "\n"); pd.DataFrame(ex).to_csv(fh, sep="\t", index=False, float_format="%.10g", na_rep="NA")
json.dump(chk, open(os.path.join(HERE, "q6_identity_check.json"), "w"))
pd.set_option("display.width", 250)
print(S.drop(columns="description").round(3).to_string()); print(Dd.round(3).to_string())
