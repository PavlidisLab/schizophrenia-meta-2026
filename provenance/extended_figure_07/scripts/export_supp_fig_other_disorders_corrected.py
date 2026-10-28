#!/usr/bin/env python3
"""Extended Fig. 8 a, b, c with correlation-aware intervals: panel TSVs, README.md, STATUS.md, check_status.py.

Reads result TSVs only (no model fits, no network). Writes to figure_data/supp_fig_other_disorders_corrected/
and copies the TSVs to analyses/figure_data_exports/out_supp_fig_other_disorders_corrected/.

  a = number of the 209 among DE genes against the matched expectation   <- addendum4 item_e_panel_a_corrected.tsv
  b = direction share: all measured genes (addendum3 concordance_vif.tsv) and DE / nominal sets (addendum4 item_e_panel_b_direction.tsv)
  c = pH-axis, slow-death-axis and Hatano pH-axis effect                  <- addendum3 axis_z_vif.tsv (+ permutation check file)
series_set is read from the source files (restricted | all series | reference); the sources take it from other_disorders_shared_controls/series_classification.tsv.
"""
import csv, hashlib, os, shutil, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
AN = os.path.join(REPO, "analyses")
S12 = os.path.join(AN, "gene_correlation_adjusted_tests")
A3 = os.path.join(S12, "addendum3_2026-10-12/other_disorders_and_scan")
A4 = os.path.join(S12, "axes_and_recurrence_dependence_tests/axis_prg_proteomics_rodent/item_e_de_restricted")
OUT = os.path.expanduser("figure_data/supp_fig_other_disorders_corrected")
COPY = os.path.join(HERE, "out_supp_fig_other_disorders_corrected")
SCRIPT = "analyses/figure_data_exports/export_supp_fig_other_disorders_corrected.py"
TODAY = datetime.date.today().isoformat()

SRC = {
    "a": os.path.join(A4, "item_e_panel_a_corrected.tsv"),
    "b_de": os.path.join(A4, "item_e_panel_b_direction.tsv"),
    "b_all": os.path.join(A3, "concordance_vif.tsv"),
    "c": os.path.join(A3, "axis_z_vif.tsv"),
    "c_perm": os.path.join(A3, "axis_z_vif_permutation_check.tsv"),
    "classif": os.path.join(AN, "other_disorders_shared_controls/series_classification.tsv"),
    "old_a": os.path.expanduser("figure_data/supp_fig_other_disorders/panel_a_enrichment_vs_null.tsv"),
    "enr209": os.path.join(AN, "other_disorders_209_2026-10-12/enrichment_209_by_disorder.tsv"),
    "dir209": os.path.join(AN, "other_disorders_209_2026-10-12/direction_209_by_disorder.tsv"),
    "axes209": os.path.join(AN, "other_disorders_209_2026-10-12/axes_within_disorder_209.tsv"),
}


def sha_full(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(p):
    return os.path.relpath(p, REPO) if p.startswith(REPO) else p.replace(os.path.expanduser("~"), "~")


def rows(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


def f(x, nd=4):
    """number as text with nd decimals; empty stays empty"""
    if x is None or x == "":
        return ""
    return f"{float(x):.{nd}f}".rstrip("0").rstrip(".") if nd else str(x)


def num(x):
    return "" if x in (None, "") else x


def header(claim, srcs):
    h = [f"# claim: {claim}", f"# script: {SCRIPT} (reads result TSVs only), run {TODAY}"]
    for k in srcs:
        h.append(f"# input: {rel(SRC[k])} sha256 {sha_full(SRC[k])[:12]}")
    h.append(f"# date: {TODAY}")
    h.append("# columns: effect, then lo95_corrected / hi95_corrected (primary interval), old_lo95 / old_hi95 (independent-genes interval), vif, method, n, source_file; p_value_tsv_only is the last column and is not for the figure or text")
    return h


def write(name, claim, srcs, cols, data, extra_comments=()):
    p = os.path.join(OUT, name)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        for l in header(claim, srcs) + list(extra_comments):
            fh.write(l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in data:
            w.writerow([r.get(c, "") for c in cols])
    return p


def excl(lo, hi):
    if lo == "" or hi == "":
        return ""
    return "yes" if (float(lo) > 0 or float(hi) < 0) else "no"


# ---- panel a
A_COLS = ["table", "disorder", "series_set", "fdr_thr", "observed", "matched_expected", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95",
          "vif", "method", "n", "source_file", "excludes_zero_corrected", "excludes_zero_old", "vif_min_over_variants", "vif_max_over_variants", "ratio", "ratio_lo_corrected", "ratio_hi_corrected",
          "z_corrected", "proxy_series", "jackknife_lo", "jackknife_hi", "p_value_kind", "p_value_tsv_only"]
a_data = []
for r in rows(SRC["a"]):
    a_data.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], fdr_thr=r["fdr_thr"], observed=r["observed"], matched_expected=r["matched_expected"],
                       matched_sd=r["matched_sd"], effect=r["effect"], lo95_corrected=r["corr_ci_lo"], hi95_corrected=r["corr_ci_hi"], old_lo95=r["old_ci_lo"], old_hi95=r["old_ci_hi"],
                       vif=r["VIF"], method="permutation", n=r["n_209_measured"], source_file=rel(SRC["a"]),
                       excludes_zero_corrected="yes" if r["excludes_zero_corrected"] == "True" else "no",
                       excludes_zero_old="yes" if r["excludes_zero_old"] == "True" else "no",
                       vif_min_over_variants=r["VIF_min_over_variants"], vif_max_over_variants=r["VIF_max_over_variants"], ratio=r["ratio"],
                       ratio_lo_corrected=r["ratio_ci_lo"], ratio_hi_corrected=r["ratio_ci_hi"], z_corrected=r["z_corrected"], proxy_series=r["proxy"],
                       jackknife_lo=r["jack_ci_lo"], jackknife_hi=r["jack_ci_hi"], p_value_kind="old matched-set permutation P, uncorrected", p_value_tsv_only=r["perm_p_old"]))

# ---- panel b
B_COLS = ["table", "disorder", "series_set", "set", "concordant", "n_with_direction", "share", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95",
          "vif", "method", "n", "source_file", "excludes_zero_corrected", "share_wilson_lo", "share_wilson_hi", "excess_over_50", "over50_lo_corrected", "over50_hi_corrected",
          "over50_excludes_zero_corrected", "note", "p_value_kind", "p_value_tsv_only"]
b_data = []
for r in rows(SRC["b_all"]):
    b_data.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], set="all measured", concordant=r["concordant"], n_with_direction=r["n_genes"],
                       share=r["frac_concordant"], matched_mean=r["null_frac_mean"], matched_sd=r["null_frac_sd"], effect=r["effect_vs_matched"],
                       lo95_corrected=r["effect_ci95_lo_corrected"], hi95_corrected=r["effect_ci95_hi_corrected"], old_lo95=r["effect_ci95_lo_old"], old_hi95=r["effect_ci95_hi_old"],
                       vif=r["VIF_concordance"], method="VIF-calibrated", n=r["n_genes"], source_file=rel(SRC["b_all"]),
                       excludes_zero_corrected=excl(r["effect_ci95_lo_corrected"], r["effect_ci95_hi_corrected"]),
                       excess_over_50=r["excess_over_50"], over50_lo_corrected=r["ci95_over_50_lo_corrected"], over50_hi_corrected=r["ci95_over_50_hi_corrected"],
                       over50_excludes_zero_corrected=excl(r["ci95_over_50_lo_corrected"], r["ci95_over_50_hi_corrected"]),
                       note="interval of the share minus the matched mean; VIF from the five-collection residual correlation",
                       p_value_kind="corrected one-sided P against the matched mean", p_value_tsv_only=r["p_corrected_matched_one_sided"]))
for r in rows(SRC["b_de"]):
    nd = int(r["n_with_direction"])
    d = dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], set=r["set"], concordant=r["concordant"], n_with_direction=r["n_with_direction"],
             share=r["share"], matched_mean=r["matched_mean"], matched_sd=r["matched_sd"], effect=r["difference"], lo95_corrected=r["diff_B_lo"], hi95_corrected=r["diff_B_hi"],
             old_lo95=r["diff_A_lo"], old_hi95=r["diff_A_hi"], vif="", method="Wilson", n=r["n_with_direction"], source_file=rel(SRC["b_de"]),
             excludes_zero_corrected=excl(r["diff_B_lo"], r["diff_B_hi"]), share_wilson_lo=r["wilson_lo"], share_wilson_hi=r["wilson_hi"],
             excess_over_50=("" if nd == 0 else f"{float(r['share']) - 0.5:.4f}"), over50_lo_corrected=r["share_minus_half_wilson_lo"], over50_hi_corrected=r["share_minus_half_wilson_hi"],
             over50_excludes_zero_corrected=excl(r["share_minus_half_wilson_lo"], r["share_minus_half_wilson_hi"]),
             note=r["note"] or "interval of the difference = Wilson SE combined with the matched SD; no VIF (DE set chosen from the data, small denominator)",
             p_value_kind="old matched-set permutation P, uncorrected", p_value_tsv_only=r["perm_p_old"])
    b_data.append(d)

# ---- panel c
C_COLS = ["table", "disorder", "series_set", "axis", "rho_obs", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95", "vif", "method", "n",
          "source_file", "excludes_zero_corrected", "z_old", "z_corrected", "vif_gtex_ba9", "vif_gtex_cortex", "vif_uniform_mean_pairwise", "p_value_kind", "p_value_tsv_only"]
c_data = []
for r in rows(SRC["c"]):
    c_data.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], axis=r["axis"], rho_obs=r["rho_obs"], matched_mean=r["null_mean"], matched_sd=r["null_sd"],
                       effect=r["effect"], lo95_corrected=r["effect_ci95_lo_corrected"], hi95_corrected=r["effect_ci95_hi_corrected"], old_lo95=r["effect_ci95_lo_old"],
                       old_hi95=r["effect_ci95_hi_old"], vif=r["VIF_axis_weighted"], method="VIF-calibrated", n=r["n_genes"], source_file=rel(SRC["c"]),
                       excludes_zero_corrected=excl(r["effect_ci95_lo_corrected"], r["effect_ci95_hi_corrected"]), z_old=r["z_old"], z_corrected=r["z_corrected"],
                       vif_gtex_ba9=r["VIF_axis_weighted_gtexBA9"], vif_gtex_cortex=r["VIF_axis_weighted_gtexCortex"], vif_uniform_mean_pairwise=r["VIF_uniform_mean_pairwise"],
                       p_value_kind="corrected one-sided P (normal tail)", p_value_tsv_only=r["p_corrected_one_sided"]))
P_COLS = ["axis", "n_genes", "rho_observed_SZ_five_collections", "perm_rho_mean", "perm_rho_sd", "perm_var_times_m", "vif_axis_weighted", "sd_if_independent_genes", "sd_predicted_by_vif",
          "method", "source_file"]
p_data = [dict(axis=r["axis"], n_genes=r["n_genes"], rho_observed_SZ_five_collections=r["rho_observed_SZ_five_collections"], perm_rho_mean=r["perm_rho_mean"], perm_rho_sd=r["perm_rho_sd"],
               perm_var_times_m=r["perm_var_times_m"], vif_axis_weighted=r["VIF_axis_weighted_primary"], sd_if_independent_genes=r["sd_if_independent_genes"],
               sd_predicted_by_vif=r["sd_predicted_by_VIF"], method="permutation", source_file=rel(SRC["c_perm"])) for r in rows(SRC["c_perm"])]

os.makedirs(OUT, exist_ok=True)
files = []
files.append(write("panel_a_enrichment_vs_null.tsv",
                   "Ext. Fig. 8a: number of the 209 among DE genes against the matched expectation; effect = observed minus matched mean; corrected interval from the case-label permutation variance inflation (proxy series)",
                   ["a", "old_a", "enr209"], A_COLS, a_data,
                   ["# note: G18-RNAseq rows are the Gandal 2018 PsychENCODE RNA-seq meta-analysis; 'Gemma ... (lists with <= 18.3% DE)' rows drop lists with more than 18.3% DE genes"]))
files.append(write("panel_b_concordance_vs_null.tsv",
                   "Ext. Fig. 8b: share of the 209 in the schizophrenia direction; set 'all measured' = VIF-calibrated interval of share minus matched mean; DE and nominal sets = Wilson-based interval, no VIF",
                   ["b_all", "b_de", "dir209"], B_COLS, b_data,
                   ["# note: for 'all measured' rows over50_* is the VIF-corrected interval of the share minus 0.5; for DE rows it is the Wilson interval of the share minus 0.5"]))
files.append(write("panel_c_axis_effect.tsv",
                   "Ext. Fig. 8c: Spearman rho of each table's signed z with the pH and slow-death axes (and Hatano pH r) over the 209, minus the matched-set mean, with VIF-calibrated interval",
                   ["c", "axes209"], C_COLS, c_data,
                   ["# note: axis 'control pH Z' = pH axis, 'GTEx agonal t' = slow-death axis; the three Hatano pH r axes were in the earlier panel_e_axis_z.tsv and are kept"]))
files.append(write("panel_c_vif_permutation_check.tsv",
                   "Ext. Fig. 8c support: variance of rho under case/control label permutation within the five reanalysed collections against the VIF (direct check of the VIF)",
                   ["c_perm"], P_COLS, p_data, ["# seed 20261012; 2,000 permutations"]))

# ---- verdict counts for README / STATUS
def cnt(data, key="excludes_zero_corrected", sets=None, pred=lambda r: True):
    sel = [r for r in data if r.get(key) in ("yes", "no") and pred(r)]
    return sum(r[key] == "yes" for r in sel), len(sel)

summary = dict(
    a_main=cnt(a_data, pred=lambda r: r["series_set"] == "restricted"),
    b_all=cnt(b_data, pred=lambda r: r["set"] == "all measured"),
    c=cnt(c_data),
    n_a=len(a_data), n_b=len(b_data), n_c=len(c_data), n_p=len(p_data),
    b_nocorr=[(r["table"], r["set"], r["note"]) for r in b_data if r["lo95_corrected"] == ""],
    a_nocorr=[(r["table"], r["fdr_thr"]) for r in a_data if r["lo95_corrected"] == ""],
    c_pos=sum(float(r["effect"]) > 0 for r in c_data),
)
import json
print(json.dumps(summary, indent=1))

# ---- STATUS hash lines
with open(os.path.join(OUT, "STATUS.md"), "w", encoding="utf-8") as fh:
    fh.write(f"""# supp_fig_other_disorders_corrected: status

Generated {TODAY} by `{SCRIPT}` from the corrected result TSVs of `analyses/gene_correlation_adjusted_tests/` (addendum3 other_disorders_and_scan, addendum4 item_e_de_restricted). The earlier folder `supp_fig_other_disorders/` is unchanged. Run `python3 check_status.py` for the live-hash check; it compares the hashes below (full sha256 at generation) with the files on disk now.

## Panel files
""" + "".join(f"- `{os.path.basename(p)}`: {sum(1 for l in open(p) if not l.startswith('#')) - 1} data rows\n" for p in files) + f"""
## Source hashes at generation (checked by check_status.py)
""" + "".join(f"- sha256 {sha_full(SRC[k])} {rel(SRC[k]).replace(os.path.expanduser('~'), '~') if not rel(SRC[k]).startswith('analyses') else 'REPO/' + rel(SRC[k])}\n" for k in SRC) + """
Status at generation: CURRENT (hashes taken at generation).
""")
shutil.copy(os.path.join(HERE, "check_status_supp_fig_other_disorders_corrected.py"), os.path.join(OUT, "check_status.py")) if os.path.exists(os.path.join(HERE, "check_status_supp_fig_other_disorders_corrected.py")) else None
os.makedirs(COPY, exist_ok=True)
for p in files:
    shutil.copy(p, COPY)
