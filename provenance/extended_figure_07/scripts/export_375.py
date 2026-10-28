#!/usr/bin/env python3
"""Extended Fig. 8 a, b, c with correlation-aware intervals for the re-frozen recurrent sets (375 primary; ge10 = 107 genes; ge5 = 2,370 genes, sensitivity).

Adapted from analyses/figure_exports/export_supp_fig_other_disorders_corrected.py (same row logic, same column order for the primary set). Reads result TSVs only.
Usage: export_375.py <by_set root> <output dir> <label:dir>...   e.g.  by_set out  375:s375 ge10:ge10 ge5:ge5
Each <dir> under <by_set root> holds the files of one run: item_e_panel_a_corrected.tsv, item_e_panel_b_direction.tsv, concordance_vif.tsv, axis_z_vif.tsv, axis_z_vif_permutation_check.tsv.
The first label is the primary set (columns as in the earlier files); each further label adds columns with the suffix _<label>.
"""
import csv, hashlib, os, sys, datetime

BY, OUT = sys.argv[1], sys.argv[2]
SETS = [a.split(":") for a in sys.argv[3:]]           # [(label, dir)]
PRIMARY = SETS[0][0]
TODAY = datetime.date.today().isoformat()
SCRIPT = "export_375.py (adapted from analyses/figure_exports/export_supp_fig_other_disorders_corrected.py)"
NAMES = {"a": "item_e_panel_a_corrected.tsv", "b_de": "item_e_panel_b_direction.tsv", "b_all": "concordance_vif.tsv", "c": "axis_z_vif.tsv", "c_perm": "axis_z_vif_permutation_check.tsv"}
SETDESC = {"375": "375 genes (set study_brain, primary)", "ge10": "107 genes (set study_brain_ge10)", "ge5": "2,370 genes (set study_brain_ge5, sensitivity)"}


def sha_full(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rows(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


def excl(lo, hi):
    if lo == "" or hi == "":
        return ""
    return "yes" if (float(lo) > 0 or float(hi) < 0) else "no"


def src(label, d, k):
    return os.path.join(BY, d, NAMES[k])


def rel(label, d, k):
    return f"by_set/{d}/{NAMES[k]}"


# ---------------------------------------------------------------- per-set row dictionaries
def a_rows(label, d):
    out = []
    for r in rows(src(label, d, "a")):
        out.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], fdr_thr=r["fdr_thr"], observed=r["observed"], matched_expected=r["matched_expected"],
                        matched_sd=r["matched_sd"], effect=r["effect"], lo95_corrected=r["corr_ci_lo"], hi95_corrected=r["corr_ci_hi"], old_lo95=r["old_ci_lo"], old_hi95=r["old_ci_hi"],
                        vif=r["VIF"], method="permutation", n=r["n_209_measured"], source_file=rel(label, d, "a"),
                        excludes_zero_corrected="yes" if r["excludes_zero_corrected"] == "True" else "no",
                        excludes_zero_old="yes" if r["excludes_zero_old"] == "True" else "no",
                        vif_min_over_variants=r["VIF_min_over_variants"], vif_max_over_variants=r["VIF_max_over_variants"], ratio=r["ratio"],
                        ratio_lo_corrected=r["ratio_ci_lo"], ratio_hi_corrected=r["ratio_ci_hi"], z_corrected=r["z_corrected"], proxy_series=r["proxy"],
                        jackknife_lo=r["jack_ci_lo"], jackknife_hi=r["jack_ci_hi"], p_value_kind="old matched-set permutation P, uncorrected", p_value_tsv_only=r["perm_p_old"]))
    return out


def b_rows(label, d):
    out = []
    for r in rows(src(label, d, "b_all")):
        out.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], set="all measured", concordant=r["concordant"], n_with_direction=r["n_genes"],
                        share=r["frac_concordant"], matched_mean=r["null_frac_mean"], matched_sd=r["null_frac_sd"], effect=r["effect_vs_matched"],
                        lo95_corrected=r["effect_ci95_lo_corrected"], hi95_corrected=r["effect_ci95_hi_corrected"], old_lo95=r["effect_ci95_lo_old"], old_hi95=r["effect_ci95_hi_old"],
                        vif=r["VIF_concordance"], method="VIF-calibrated", n=r["n_genes"], source_file=rel(label, d, "b_all"),
                        excludes_zero_corrected=excl(r["effect_ci95_lo_corrected"], r["effect_ci95_hi_corrected"]),
                        excess_over_50=r["excess_over_50"], over50_lo_corrected=r["ci95_over_50_lo_corrected"], over50_hi_corrected=r["ci95_over_50_hi_corrected"],
                        over50_excludes_zero_corrected=excl(r["ci95_over_50_lo_corrected"], r["ci95_over_50_hi_corrected"]),
                        note="interval of the share minus the matched mean; VIF from the five-collection residual correlation",
                        p_value_kind="corrected one-sided P against the matched mean", p_value_tsv_only=r["p_corrected_matched_one_sided"]))
    for r in rows(src(label, d, "b_de")):
        nd = int(r["n_with_direction"])
        out.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], set=r["set"], concordant=r["concordant"], n_with_direction=r["n_with_direction"],
                        share=r["share"], matched_mean=r["matched_mean"], matched_sd=r["matched_sd"], effect=r["difference"], lo95_corrected=r["diff_B_lo"], hi95_corrected=r["diff_B_hi"],
                        old_lo95=r["diff_A_lo"], old_hi95=r["diff_A_hi"], vif="", method="Wilson", n=r["n_with_direction"], source_file=rel(label, d, "b_de"),
                        excludes_zero_corrected=excl(r["diff_B_lo"], r["diff_B_hi"]), share_wilson_lo=r["wilson_lo"], share_wilson_hi=r["wilson_hi"],
                        excess_over_50=("" if nd == 0 else f"{float(r['share']) - 0.5:.4f}"), over50_lo_corrected=r["share_minus_half_wilson_lo"], over50_hi_corrected=r["share_minus_half_wilson_hi"],
                        over50_excludes_zero_corrected=excl(r["share_minus_half_wilson_lo"], r["share_minus_half_wilson_hi"]),
                        note=r.get("note") or "interval of the difference = Wilson SE combined with the matched SD; no VIF (DE set chosen from the data, small denominator)",
                        p_value_kind="old matched-set permutation P, uncorrected", p_value_tsv_only=r["perm_p_old"]))
    return out


def c_rows(label, d):
    out = []
    for r in rows(src(label, d, "c")):
        out.append(dict(table=r["table"], disorder=r["disorder"], series_set=r["series_set"], axis=r["axis"], rho_obs=r["rho_obs"], matched_mean=r["null_mean"], matched_sd=r["null_sd"],
                        effect=r["effect"], lo95_corrected=r["effect_ci95_lo_corrected"], hi95_corrected=r["effect_ci95_hi_corrected"], old_lo95=r["effect_ci95_lo_old"],
                        old_hi95=r["effect_ci95_hi_old"], vif=r["VIF_axis_weighted"], method="VIF-calibrated", n=r["n_genes"], source_file=rel(label, d, "c"),
                        excludes_zero_corrected=excl(r["effect_ci95_lo_corrected"], r["effect_ci95_hi_corrected"]), z_old=r["z_old"], z_corrected=r["z_corrected"],
                        vif_gtex_ba9=r["VIF_axis_weighted_gtexBA9"], vif_gtex_cortex=r["VIF_axis_weighted_gtexCortex"], vif_uniform_mean_pairwise=r["VIF_uniform_mean_pairwise"],
                        p_value_kind="corrected one-sided P (normal tail)", p_value_tsv_only=r["p_corrected_one_sided"]))
    return out


def p_rows(label, d):
    return [dict(axis=r["axis"], n_genes=r["n_genes"], rho_observed_SZ_five_collections=r["rho_observed_SZ_five_collections"], perm_rho_mean=r["perm_rho_mean"], perm_rho_sd=r["perm_rho_sd"],
                 perm_var_times_m=r["perm_var_times_m"], vif_axis_weighted=r["VIF_axis_weighted_primary"], sd_if_independent_genes=r["sd_if_independent_genes"],
                 sd_predicted_by_vif=r["sd_predicted_by_VIF"], method="permutation", source_file=rel(label, d, "c_perm")) for r in rows(src(label, d, "c_perm"))]


A_COLS = ["table", "disorder", "series_set", "fdr_thr", "observed", "matched_expected", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95",
          "vif", "method", "n", "source_file", "excludes_zero_corrected", "excludes_zero_old", "vif_min_over_variants", "vif_max_over_variants", "ratio", "ratio_lo_corrected", "ratio_hi_corrected",
          "z_corrected", "proxy_series", "jackknife_lo", "jackknife_hi", "p_value_kind", "p_value_tsv_only"]
B_COLS = ["table", "disorder", "series_set", "set", "concordant", "n_with_direction", "share", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95",
          "vif", "method", "n", "source_file", "excludes_zero_corrected", "share_wilson_lo", "share_wilson_hi", "excess_over_50", "over50_lo_corrected", "over50_hi_corrected",
          "over50_excludes_zero_corrected", "note", "p_value_kind", "p_value_tsv_only"]
C_COLS = ["table", "disorder", "series_set", "axis", "rho_obs", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95", "vif", "method", "n",
          "source_file", "excludes_zero_corrected", "z_old", "z_corrected", "vif_gtex_ba9", "vif_gtex_cortex", "vif_uniform_mean_pairwise", "p_value_kind", "p_value_tsv_only"]
P_COLS = ["axis", "n_genes", "rho_observed_SZ_five_collections", "perm_rho_mean", "perm_rho_sd", "perm_var_times_m", "vif_axis_weighted", "sd_if_independent_genes", "sd_predicted_by_vif",
          "method", "source_file"]
# columns added per further set (suffix _<label>)
A_X = ["n", "observed", "matched_expected", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "old_lo95", "old_hi95", "vif", "excludes_zero_corrected", "z_corrected", "ratio", "jackknife_lo", "jackknife_hi"]
B_X = ["concordant", "n_with_direction", "share", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "vif", "excludes_zero_corrected", "excess_over_50", "over50_lo_corrected", "over50_hi_corrected", "over50_excludes_zero_corrected"]
C_X = ["n", "rho_obs", "matched_mean", "matched_sd", "effect", "lo95_corrected", "hi95_corrected", "vif", "excludes_zero_corrected", "z_corrected"]
P_X = ["n_genes", "rho_observed_SZ_five_collections", "perm_rho_mean", "perm_rho_sd", "perm_var_times_m", "vif_axis_weighted"]


def merge(fn, keys, base_cols, x_cols):
    per = {lab: fn(lab, d) for lab, d in SETS}
    key = lambda r: tuple(r[k] for k in keys)
    order, rowsd = [], {}
    for lab, _ in SETS:                      # primary rows first, then rows only another set has
        for r in per[lab]:
            if key(r) not in rowsd:
                rowsd[key(r)] = {}
                order.append(key(r))
            rowsd[key(r)][lab] = r
    cols = list(base_cols)
    for lab, _ in SETS[1:]:
        cols += [f"{c}_{lab}" for c in x_cols]
    data, only_other = [], 0
    for k in order:
        d = dict(rowsd[k].get(PRIMARY, {}))
        if PRIMARY not in rowsd[k]:
            only_other += 1
            for kk, v in zip(keys, k):
                d[kk] = v
        for lab, _ in SETS[1:]:
            r = rowsd[k].get(lab)
            for c in x_cols:
                d[f"{c}_{lab}"] = r.get(c, "") if r else ""
        data.append(d)
    return cols, data, {lab: len(per[lab]) for lab, _ in SETS}, only_other


def header(claim, ks):
    h = [f"# claim: {claim}", f"# script: {SCRIPT} (reads result TSVs only), run {TODAY}",
         "# gene sets: " + "; ".join(f"{lab}: {SETDESC.get(lab, lab)}" for lab, _ in SETS) + f". Columns without a suffix are the {PRIMARY} set; columns with a suffix _<set> are the other sets."]
    for lab, d in SETS:
        for k in ks:
            p = src(lab, d, k)
            h.append(f"# input ({lab}): {rel(lab, d, k)} sha256 {sha_full(p)[:12]}")
    h.append(f"# date: {TODAY}")
    h.append("# columns: effect, then lo95_corrected / hi95_corrected (primary interval), old_lo95 / old_hi95 (independent-genes interval), vif, method, n, source_file; p_value_tsv_only is not for the figure or text")
    return h


def write(name, claim, ks, cols, data, extra=()):
    p = os.path.join(OUT, name)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        for l in header(claim, ks) + list(extra):
            fh.write(l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in data:
            w.writerow([r.get(c, "") for c in cols])
    return p


os.makedirs(OUT, exist_ok=True)
info = {}
cols, a_data, n_, oo = merge(a_rows, ["table", "fdr_thr"], A_COLS, A_X); info["a"] = (n_, oo)
write("panel_a_enrichment_vs_null.tsv",
      "Ext. Fig. 8a: number of the recurrent genes among DE genes against the matched expectation; effect = observed minus matched mean; corrected interval from the case-label permutation variance inflation (proxy series)",
      ["a"], cols, a_data,
      ["# note: G18-RNAseq rows are the Gandal 2018 PsychENCODE RNA-seq meta-analysis; 'Gemma ... (lists with <= 18.3% DE)' rows drop lists with more than 18.3% DE genes; n = genes of the set measured in the table"])
cols, b_data, n_, oo = merge(b_rows, ["table", "set"], B_COLS, B_X); info["b"] = (n_, oo)
write("panel_b_concordance_vs_null.tsv",
      "Ext. Fig. 8b: share of the recurrent genes in the schizophrenia direction; set 'all measured' = VIF-calibrated interval of share minus matched mean; DE and nominal sets = Wilson-based interval, no VIF",
      ["b_all", "b_de"], cols, b_data,
      ["# note: for 'all measured' rows over50_* is the VIF-corrected interval of the share minus 0.5; for DE rows it is the Wilson interval of the share minus 0.5"])
cols, c_data, n_, oo = merge(c_rows, ["table", "axis"], C_COLS, C_X); info["c"] = (n_, oo)
write("panel_c_axis_effect.tsv",
      "Ext. Fig. 8c: Spearman rho of each table's signed z with the pH and slow-death axes (and Hatano pH r) over the recurrent genes, minus the matched-set mean, with VIF-calibrated interval",
      ["c"], cols, c_data,
      ["# note: axis 'control pH Z' = pH axis, 'GTEx agonal t' = slow-death axis; the three Hatano pH r axes were in the earlier panel_e_axis_z.tsv and are kept"])
cols, p_data, n_, oo = merge(p_rows, ["axis"], P_COLS, P_X); info["p"] = (n_, oo)
write("panel_c_vif_permutation_check.tsv",
      "Ext. Fig. 8c support: variance of rho under case/control label permutation within the five reanalysed collections against the VIF (direct check of the VIF)",
      ["c_perm"], cols, p_data, ["# seed 20261012; 2,000 permutations"])
for k, v in info.items():
    print(k, "rows per set", v[0], "rows present only in a non-primary set:", v[1])
