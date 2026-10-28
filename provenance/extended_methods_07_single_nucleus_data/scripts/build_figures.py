#!/usr/bin/env python3
"""Tables, statistics and figures for the SZ vs control effect-size comparison.

    .claude/skills/plotting/.venv/bin/python effect_sizes/build_figures.py

Reads only files on disk (written by harmonize.py and ash_fit.R):
    data/sources.tsv, data/by_group.tsv, data/viewA.tsv.gz, data/excluded.tsv
Writes:
    effects_by_study.tsv   one row per source (headline studies, pooled analyses, Gemma datasets,
                           Endresz cohorts), with a '#' provenance header
    summary_by_type.tsv    headline studies per type
    stats.json             every number quoted in effect_sizes_report.md
    fig_1..fig_4 .svg      figures (plotting skill helpers: pavlab_stripchart, pavlab_scatter)
    CAPTIONS.md, index.html
"""
import csv
import datetime
import hashlib
import json
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

sys.path.insert(0, os.path.expanduser(".claude/skills/plotting/python"))
sys.path.insert(0, os.path.expanduser(".claude/skills/architecture-figures/python"))
from pavlab_stripchart import pavlab_stripchart  # noqa: E402
from pavlab_scatter import pavlab_scatter  # noqa: E402
from pavlab_style import apply_pavlab_style, lighten  # noqa: E402
from pavlab_arch.palette import OKABE_ITO  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)
ES = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ES, "data")
TODAY = datetime.date.today().isoformat()
SUBTLE, TEXT = "#6b7280", "#1f2937"
SEED = 20260924
N_BOOT = 10000
SD_FLOOR = 0.003  # where SD = 0 (ash put all weight on the null) is drawn on log axes

ARRAY, RNASEQ = "bulk brain microarray", "bulk brain RNA-seq"
BLOOD, IPSC, OTHER = "blood", "iPSC/organoid (bulk)", "other peripheral / cell line"
SC, IPSC_SC = "single-nucleus brain", "iPSC/organoid (single-cell)"
TYPES = [ARRAY, RNASEQ, BLOOD, IPSC, OTHER, SC, IPSC_SC]
SHORT = {ARRAY: "brain\narray", RNASEQ: "brain\nRNA-seq", BLOOD: "blood", IPSC: "iPSC /\norganoid", OTHER: "other\nperipheral",
         SC: "single-\nnucleus", IPSC_SC: "organoid\nscRNA-seq"}
COLOR = {ARRAY: OKABE_ITO["blue"], RNASEQ: OKABE_ITO["sky_blue"], BLOOD: OKABE_ITO["orange"],
         IPSC: OKABE_ITO["reddish_purple"], OTHER: "#999999", SC: OKABE_ITO["vermillion"], IPSC_SC: OKABE_ITO["bluish_green"]}
PROV_RANK = {"published": 0, "reproduced": 1, "gemma": 2}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def read_tsv(name, **kw):
    return pd.read_csv(os.path.join(DATA, name), sep="\t", comment="#", quoting=csv.QUOTE_NONE, dtype={"study": str}, **kw)


def source_line(fig, text, y=-0.02):
    fig.text(0.01, y, text, fontsize=8, color=SUBTLE, ha="left", va="top")


# ------------------------------------------------------------------ load and aggregate
def load():
    src = read_tsv("sources.tsv")
    grp = read_tsv("by_group.tsv")
    A = pd.read_csv(os.path.join(DATA, "viewA.tsv.gz"), sep="\t", dtype={"study": str})
    A["abs"] = A["log2fc"].abs()
    a = A.groupby("study")["abs"].agg(A_n="size", A_median="median").reset_index()
    g = grp.groupby("source").agg(
        n_groups=("group", "size"), B_median=("median_abs_B", "median"), B_p90=("p90_abs_B", "median"),
        n_genes=("n_genes_B", "median"), C_sd=("sd_g", "median"), C_sd_max=("sd_g", "max"),
        C_mean_abs_pm=("mean_abs_pm", "median"), C_frac_lfsr05=("frac_lfsr05", "median"), C_pi0=("pi0", "median"),
        C_n_fitted=("sd_g", lambda s: int(s.notna().sum())), C_se_method=("se_method", "first"),
        C_median_se=("median_se", "median")).reset_index()
    top = grp.loc[grp["sd_g"].notna()].sort_values("sd_g").groupby("source").tail(1)[["source", "group", "sd_g"]]
    top.columns = ["source", "C_top_group", "C_top_sd"]
    g = g.merge(top, on="source", how="left")
    d = src.merge(g, left_on="full_source", right_on="source", how="left").merge(a, on="study", how="left")
    d["A_n"] = d["A_n"].fillna(0).astype(int)
    return d, grp, A


def pick_headline(d, have):
    """One study per donor cluster: published > reproduced > Gemma, then the largest N with the
    view available, then label."""
    h = d[(d["role"] == "headline") & have].copy()
    h["_n"] = h["N"].fillna(-1)
    h["_p"] = h["full_provenance"].map(PROV_RANK).fillna(3)
    h = h.sort_values(["cluster", "_p", "_n", "label"], ascending=[True, True, False, True])
    return h.groupby("cluster").head(1)


def boot_median_ci(x, rng):
    x = np.asarray(x, float)
    if len(x) < 2:
        return (np.nan, np.nan)
    b = np.median(rng.choice(x, size=(N_BOOT, len(x)), replace=True), axis=1)
    return tuple(np.percentile(b, [2.5, 97.5]))


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y))] + X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * sps.t.sf(np.abs(t), dof)
    r2 = 1 - (resid @ resid) / ((y - y.mean()) @ (y - y.mean()))
    lo, hi = beta - sps.t.ppf(0.975, dof) * se, beta + sps.t.ppf(0.975, dof) * se
    return {"n": int(len(y)), "df": int(dof), "r2": float(r2),
            "coef": {nm: {"b": float(b), "se": float(s), "t": float(tt), "p": float(pp), "ci": [float(l), float(h)],
                          "fold": float(10 ** b), "fold_ci": [float(10 ** l), float(10 ** h)]}
                     for nm, b, s, tt, pp, l, h in zip(["intercept"] + names, beta, se, t, p, lo, hi)}}


def regress(h, types, ref):
    """log10(SD) ~ type (reference `ref`) + log10(N), SD > 0 and N known."""
    x = h[h["type"].isin(types) & (h["C_sd"] > 0) & h["N"].notna()].copy()
    others = [t for t in types if t != ref and (x["type"] == t).any()]
    X = [(x["type"] == t).astype(float).values for t in others] + [np.log10(x["N"].values.astype(float))]
    r = ols(np.log10(x["C_sd"].values), X, [f"type[{t}]" for t in others] + ["log10(N)"])
    r["reference"] = ref
    r["excluded_sd0"] = int(((h["type"].isin(types)) & (h["C_sd"] == 0)).sum())
    r["excluded_noN"] = int(((h["type"].isin(types)) & (h["C_sd"] > 0) & h["N"].isna()).sum())
    return r


# ------------------------------------------------------------------ figures
def strip_by_type(ax, df, ycol, types, size_by_n=True, light=None, ylabel="", floor=None):
    """pavlab_stripchart of one point per row, x = type; recoloured per type; optional lighter
    per-group points (light = DataFrame with 'type' and ycol) drawn first."""
    order = [t for t in types if (df["type"] == t).any()]
    labels = [SHORT[t] for t in order]
    rng = np.random.default_rng(SEED)
    if light is not None and len(light):
        lt = light[light["type"].isin(order)]
        xs = np.array([order.index(t) for t in lt["type"]]) + rng.uniform(-0.3, 0.3, len(lt))
        yv = lt[ycol].values.astype(float)
        if floor:
            yv = np.where(yv > 0, yv, floor)
        ax.scatter(xs, yv, s=9, c=[lighten(COLOR[t], 0.6) for t in lt["type"]], linewidths=0, zorder=1)
    dd = df[df["type"].isin(order)].copy()
    dd["_o"] = dd["type"].map({t: i for i, t in enumerate(order)})
    dd = dd.sort_values("_o")
    yv = dd[ycol].values.astype(float)
    zero = ~(yv > 0)
    if floor:
        yv = np.where(zero, floor, yv)
    np.random.seed(SEED)
    pavlab_stripchart([SHORT[t] for t in dd["type"]], yv, order=labels, jitter=0.18, ax=ax, origin_zero=False,
                      ylabel=ylabel, label_size="small")
    sc = ax.collections[-1]
    sc.set_facecolor([COLOR[t] for t in dd["type"]])
    sc.set_edgecolor(["white" if not z else COLOR[t] for t, z in zip(dd["type"], zero)])
    sc.set_linewidth(0.6)
    if size_by_n:
        n = dd["N"].fillna(20).clip(lower=4).values.astype(float)
        sc.set_sizes(12 + 70 * (np.log10(n) - np.log10(4)) / (np.log10(1500) - np.log10(4)))
    sc.set_zorder(3)
    ax.set_yscale("log")
    if floor and zero.any():
        ax.axhline(floor, color="#e5e7eb", lw=1, ls="--", zorder=0)
        ax.text(len(order) - 0.45, floor, "SD = 0\n(all weight on null)", fontsize=7, color=SUBTLE, va="center", ha="left")
    ax.grid(axis="y", color="#e5e7eb", lw=0.6)
    ax.set_axisbelow(True)
    return order


def size_legend(ax, loc="upper left"):
    hs = []
    for n in (10, 100, 1000):
        s = 12 + 70 * (np.log10(n) - np.log10(4)) / (np.log10(1500) - np.log10(4))
        hs.append(ax.scatter([], [], s=s, c="#9ca3af", linewidths=0, label=f"N = {n}"))
    ax.legend(handles=hs, frameon=False, fontsize=8, loc=loc, title="cases + controls", title_fontsize=8)


def fig1(hC, grp, d, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    light = grp.merge(hC[["full_source", "type"]], left_on="source", right_on="full_source")
    light = light[light["type"].isin([SC]) & light["sd_g"].notna()]
    strip_by_type(axes[0], hC, "C_sd", TYPES, light=light.rename(columns={"sd_g": "C_sd"}), floor=SD_FLOOR,
                  ylabel="SD of estimated true log2 fold changes (ashr)")
    axes[0].set_title("a  One study per donor set", loc="left", fontsize=11, fontweight="normal")
    size_legend(axes[0])
    # b: uniform pipelines
    u = d[d["role"].isin(["gemma", "endresz_cohort"]) & d["C_sd"].notna()].copy()
    u["panel"] = np.where(u["role"] == "endresz_cohort", "Endresz\nsnRNA cohorts",
                          np.where(u["type"] == SC, "Gemma\nsnRNA pseudobulk",
                                   np.where(u["type"] == ARRAY, "Gemma\nbrain array",
                                            np.where(u["type"] == RNASEQ, "Gemma\nbrain RNA-seq", "other"))))
    u = u[u["panel"] != "other"]
    porder = ["Gemma\nbrain array", "Gemma\nbrain RNA-seq", "Gemma\nsnRNA pseudobulk", "Endresz\nsnRNA cohorts"]
    pc = {porder[0]: COLOR[ARRAY], porder[1]: COLOR[RNASEQ], porder[2]: COLOR[SC], porder[3]: COLOR[SC]}
    ul = grp.merge(u[["full_source", "panel"]], left_on="source", right_on="full_source")
    ul = ul[ul["panel"].isin(porder[2:]) & ul["sd_g"].notna()]
    rng = np.random.default_rng(SEED)
    xs = np.array([porder.index(p) for p in ul["panel"]]) + rng.uniform(-0.3, 0.3, len(ul))
    axes[1].scatter(xs, np.where(ul["sd_g"] > 0, ul["sd_g"], SD_FLOOR), s=9, c=lighten(COLOR[SC], 0.6), linewidths=0, zorder=1)
    u["_o"] = u["panel"].map({p: i for i, p in enumerate(porder)})
    u = u.sort_values("_o")
    yv = u["C_sd"].values.astype(float)
    zero = ~(yv > 0)
    np.random.seed(SEED)
    pavlab_stripchart(u["panel"].tolist(), np.where(zero, SD_FLOOR, yv), order=porder, jitter=0.18, ax=axes[1],
                      origin_zero=False, label_size="small")
    sc = axes[1].collections[-1]
    sc.set_facecolor([pc[p] for p in u["panel"]])
    sc.set_edgecolor(["white" if not z else pc[p] for p, z in zip(u["panel"], zero)])
    sc.set_linewidth(0.6)
    n = u["N"].fillna(20).clip(lower=4).values.astype(float)
    sc.set_sizes(12 + 70 * (np.log10(n) - np.log10(4)) / (np.log10(1500) - np.log10(4)))
    sc.set_zorder(3)
    axes[1].set_yscale("log")
    axes[1].set_ylim(axes[0].get_ylim())
    axes[1].axhline(SD_FLOOR, color="#e5e7eb", lw=1, ls="--", zorder=0)
    axes[1].grid(axis="y", color="#e5e7eb", lw=0.6)
    axes[1].set_title("b  Same pipeline across data types", loc="left", fontsize=11, fontweight="normal")
    for ax in axes:
        ax.tick_params(axis="x", labelsize=9)
    fig.tight_layout()
    source_line(fig, "Source: effect_sizes/effects_by_study.tsv, data/by_group.tsv (ashr, normal mixture with point mass at 0). "
                     "Large points: per-study median over cell types or regions; small light points: single-nucleus cell types. "
                     f"Built {TODAY}.")
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)


def fig2(hA, A, out):
    rows = [t for t in TYPES if (hA["type"] == t).any()]
    fig, axes = plt.subplots(len(rows), 1, figsize=(15, 2.7 * len(rows)))
    rng = np.random.default_rng(SEED)
    for ax, t in zip(axes, rows):
        s = hA[hA["type"] == t].sort_values(["N"], na_position="first")
        labs = [f"{l} ({'?' if pd.isna(n) else int(n)})" for l, n in zip(s["label"], s["N"])]
        xs, ys = [], []
        for sid, lab in zip(s["study"], labs):
            v = A.loc[A["study"] == sid, "abs"].values
            v = v[v > 0]
            if len(v) > 300:  # display only; medians use every row
                v = rng.choice(v, 300, replace=False)
            xs += [lab] * len(v)
            ys += list(v)
        np.random.seed(SEED)
        pavlab_stripchart(xs, ys, order=labs, jitter=0.25, ax=ax, origin_zero=False, point_size=5, label_size="small")
        ax.collections[-1].set_color(lighten(COLOR[t], 0.35))
        ax.scatter(range(len(labs)), s["A_median"].values, marker="_", s=160, c=TEXT, linewidths=2, zorder=4)
        ax.set_yscale("log")
        ax.set_ylim(0.005, 40)
        ax.grid(axis="y", color="#e5e7eb", lw=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=6.5, rotation=60)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.set_title(f"{t}  (median of study medians {np.median(s['A_median']):.2f})", loc="left", fontsize=10, fontweight="normal")
        ax.set_ylabel("|log2 FC|", fontsize=9)
    fig.tight_layout()
    source_line(fig, "Source: effect_sizes/data/viewA.tsv.gz (significant rows at each study's own threshold, converted to log2 FC). "
                     "Studies ordered by N (cases + controls, in brackets); black bar: study median over all rows; "
                     f"at most 300 rows drawn per study. Built {TODAY}.", y=-0.005)
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)


def fig3(hB, grp, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    light = grp.merge(hB[["full_source", "type"]], left_on="source", right_on="full_source")
    light = light[light["type"] == SC]
    for ax, col, lcol, ttl in [(axes[0], "B_median", "median_abs_B", "a  Median |log2 FC| over all tested genes"),
                               (axes[1], "B_p90", "p90_abs_B", "b  90th percentile of |log2 FC|")]:
        strip_by_type(ax, hB, col, TYPES, light=light.rename(columns={lcol: col}), ylabel="|log2 FC|, all tested genes")
        ax.set_title(ttl, loc="left", fontsize=11, fontweight="normal")
        ax.tick_params(axis="x", labelsize=9)
    size_legend(axes[0])
    fig.tight_layout()
    source_line(fig, "Source: effect_sizes/effects_by_study.tsv (B_median, B_p90), data/by_group.tsv. Raw estimates, not corrected for "
                     f"noise. Large points: per-study median over cell types or regions; light points: single-nucleus cell types. Built {TODAY}.")
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)


def fig4(hC, d, out):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    x = hC[(hC["N"].notna())].copy()
    x["y"] = np.where(x["C_sd"] > 0, x["C_sd"], SD_FLOOR)
    pavlab_scatter(x["N"].values.astype(float), x["y"].values, color=x["type"].values, ax=ax, origin_zero=False,
                   xlabel="N (cases + controls)", ylabel="SD of estimated true log2 FC (ashr)", label_size="small")
    ax.collections[-1].set_color([COLOR[t] for t in x["type"]])
    e = d[(d["role"] == "endresz_cohort") & d["C_sd"].notna()]
    ax.scatter(e["N"], np.where(e["C_sd"] > 0, e["C_sd"], SD_FLOOR), s=28, facecolors="none", edgecolors=COLOR[SC], linewidths=1.2, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    hs = [ax.scatter([], [], s=24, c=COLOR[t], label=t) for t in TYPES if (x["type"] == t).any()]
    hs.append(ax.scatter([], [], s=28, facecolors="none", edgecolors=COLOR[SC], linewidths=1.2, label="Endresz snRNA cohort (sensitivity)"))
    ax.legend(handles=hs, frameon=False, fontsize=8, loc="upper right")
    ax.axhline(SD_FLOOR, color="#e5e7eb", lw=1, ls="--", zorder=0)
    ax.grid(axis="y", color="#e5e7eb", lw=0.6)
    fig.tight_layout()
    source_line(fig, f"Source: effect_sizes/effects_by_study.tsv (C_sd, N). Points on the dashed line have SD = 0. Built {TODAY}.")
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    apply_pavlab_style()
    plt.rcParams.update({"axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
                         "font.family": ["Helvetica", "Arial", "sans-serif"]})
    d, grp, A = load()
    rng = np.random.default_rng(SEED)
    hC = pick_headline(d, d["C_sd"].notna())
    hB = pick_headline(d, d["B_median"].notna())
    hA = pick_headline(d, d["A_n"] > 0)

    # Endresz cohorts, one per donor cluster (MSSM1 and MSSM2 share Mount Sinai donors): the larger N
    endresz_one = (d[(d["role"] == "endresz_cohort") & d["C_sd"].notna()].sort_values(["cluster", "N"], ascending=[True, False])
                   .groupby("cluster").head(1))

    # ---- summary by type
    rows = []
    for t in TYPES:
        c, b, a = hC[hC["type"] == t], hB[hB["type"] == t], hA[hA["type"] == t]
        lo, hi = boot_median_ci(c["C_sd"].dropna(), rng)
        rows.append({"type": t, "n_studies_C": len(c), "median_C_sd": c["C_sd"].median(), "C_sd_ci_lo": lo, "C_sd_ci_hi": hi,
                     "median_C_sd_top_group": c["C_top_sd"].median() if t == SC else np.nan,
                     "median_C_mean_abs_pm": c["C_mean_abs_pm"].median(), "median_C_frac_lfsr05": c["C_frac_lfsr05"].median(),
                     "n_studies_A": len(a), "median_A_abs": a["A_median"].median(), "n_studies_B": len(b),
                     "median_B_median": b["B_median"].median(), "median_B_p90": b["B_p90"].median(),
                     "median_N_C": c["N"].median()})
    summ = pd.DataFrame(rows)

    # sensitivity: uniform pipelines
    sens = {}
    for name, sel in [("gemma_array", (d["role"] == "gemma") & (d["type"] == ARRAY)),
                      ("gemma_rnaseq", (d["role"] == "gemma") & (d["type"] == RNASEQ)),
                      ("gemma_sc", (d["role"] == "gemma") & (d["type"] == SC)),
                      ("gemma_blood", (d["role"] == "gemma") & (d["type"] == BLOOD)),
                      ("gemma_ipsc", (d["role"] == "gemma") & (d["type"] == IPSC)),
                      ("endresz_sc_cohorts", d["role"] == "endresz_cohort"),
                      ("endresz_sc_donor_sets", d["study"].isin(endresz_one["study"])),
                      ("pooled_bulk", d["role"] == "pooled")]:
        x = d[sel & d["C_sd"].notna()]
        lo, hi = boot_median_ci(x["C_sd"], rng)
        sens[name] = {"n": int(len(x)), "median_C_sd": float(x["C_sd"].median()) if len(x) else None, "ci": [lo, hi],
                      "median_C_top_sd": float(x["C_top_sd"].median()) if len(x) else None,
                      "median_B_median": float(x["B_median"].median()) if len(x) else None,
                      "studies": {l: (None if pd.isna(v) else round(float(v), 4)) for l, v in zip(x["label"], x["C_sd"])}}

    # ---- regressions
    brain = [ARRAY, RNASEQ, SC]
    reg = {"brain_ref_rnaseq": regress(hC, brain, RNASEQ),
           "all_types_ref_rnaseq": regress(hC, [t for t in TYPES if (hC["type"] == t).any()], RNASEQ)}
    # sensitivity: single-nucleus studies replaced by the Endresz cohorts
    hE = pd.concat([hC[hC["type"].isin([ARRAY, RNASEQ])], endresz_one])
    reg["brain_endresz_cohorts_ref_rnaseq"] = regress(hE, brain, RNASEQ)
    # sensitivity: Gemma only (one pipeline for bulk and pseudobulk)
    hG = d[(d["role"] == "gemma") & d["C_sd"].notna()]
    reg["gemma_brain_ref_rnaseq"] = regress(hG, brain, RNASEQ)
    # sensitivity: SD = 0 set to the floor instead of dropped
    hF = hC.copy()
    hF.loc[~(hF["C_sd"] > 0), "C_sd"] = SD_FLOOR
    reg["brain_sd0_at_floor_ref_rnaseq"] = regress(hF, brain, RNASEQ)
    # all types without the iPSC/organoid studies; within iPSC, Spearman of log SD on log N (analysis review E3)
    reg["all_types_without_ipsc_ref_rnaseq"] = regress(hC, [t for t in TYPES if (hC["type"] == t).any() and t != IPSC], RNASEQ)
    ip = hC[(hC["type"] == IPSC) & (hC["C_sd"] > 0) & hC["N"].notna()]
    rho = sps.spearmanr(np.log10(ip["N"].astype(float)), np.log10(ip["C_sd"].astype(float)))
    within_ipsc = {"n": int(len(ip)), "spearman_rho": float(rho.correlation), "p": float(rho.pvalue)}
    # likelihood used per headline source (ash_fit.R column lik)
    lik = grp.groupby("source")["lik"].agg(lambda v: ",".join(sorted(set(map(str, v)))))
    t_lik = {r["label"]: lik.get(r["full_source"], "") for _, r in hC.iterrows() if str(lik.get(r["full_source"], "")).startswith("t(")}
    # Mann-Whitney: single-nucleus vs each bulk brain type
    mw = {}
    for t in (ARRAY, RNASEQ):
        a_, b_ = hC.loc[hC["type"] == SC, "C_sd"], hC.loc[hC["type"] == t, "C_sd"]
        u = sps.mannwhitneyu(a_, b_, alternative="two-sided")
        mw[f"sc_vs_{t}"] = {"U": float(u.statistic), "p": float(u.pvalue), "n_sc": int(len(a_)), "n_other": int(len(b_))}
    # published vs Endresz for the same donors
    pairs = {}
    lab = d.set_index("full_source")
    for coh, pub in [("McLean", "38781388__published"), ("MSSM1", "38781388__published"), ("Batiuk", "36223459__published"),
                     ("Frohlich", "39227716__published"), ("HBCC", "42778699__published"), ("MSSM2", "42778699__published"),
                     ("Multiome", "38781369__published")]:
        if f"endresz__{coh}" in lab.index and pub in lab.index:
            pairs[coh] = {"endresz_C_sd": float(lab.loc[f"endresz__{coh}", "C_sd"]), "published": lab.loc[pub, "label"],
                          "published_C_sd": float(lab.loc[pub, "C_sd"])}
    for gname, pub in [("gemma__CMC", "38781369__published"), ("gemma__SZBDMulti-Seq", "38781369__published"),
                       ("gemma__MultiomeBrain", "38781369__published")]:
        if gname in lab.index:
            pairs[gname] = {"gemma_C_sd": float(lab.loc[gname, "C_sd"]), "published": lab.loc[pub, "label"],
                            "published_C_sd": float(lab.loc[pub, "C_sd"])}

    # ---- effects_by_study.tsv
    keep = d[d["role"].isin(["headline", "pooled", "gemma", "endresz_cohort"])].copy()
    keep["headline_A"] = keep["study"].isin(hA["study"])
    keep["headline_B"] = keep["study"].isin(hB["study"])
    keep["headline_C"] = keep["study"].isin(hC["study"])
    cols = ["study", "label", "type", "role", "cluster", "headline_A", "headline_B", "headline_C", "N", "n_cases", "n_controls",
            "N_evidence", "A_n", "A_median", "full_provenance", "full_source", "n_groups", "n_genes", "B_median", "B_p90",
            "C_sd", "C_top_group", "C_top_sd", "C_mean_abs_pm", "C_frac_lfsr05", "C_pi0", "C_n_fitted", "C_se_method",
            "A_conversion", "full_conversion", "flags", "type_evidence", "role_reason"]
    ins = {os.path.relpath(os.path.join(DATA, f), ES): sha(os.path.join(DATA, f))
           for f in ("sources.tsv", "by_group.tsv", "viewA.tsv.gz", "excluded.tsv")}
    hm = json.load(open(os.path.join(DATA, "harmonize_meta.json")))
    with open(os.path.join(ES, "effects_by_study.tsv"), "w", encoding="utf-8") as fh:
        fh.write(f"# built {TODAY} by effect_sizes/build_figures.py (figures: fig_1 to fig_4; see CAPTIONS.md)\n")
        fh.write("# pipeline: harmonize.py -> ash_fit.R -> build_figures.py\n")
        for k, v in ins.items():
            fh.write(f"# input: {k} sha256 {v}\n")
        for k, v in sorted(hm["inputs_sha256"].items()):
            fh.write(f"# upstream: {k} sha256 {v}\n")
        fh.write("# A_median: median |log2FC| of significant rows (view A). B_median, B_p90: |log2FC| over all tested genes (view B), "
                 "median over groups. C_sd: SD of the ashr prior (view C), median over groups; C_top_sd: largest group.\n")
        keep[cols].to_csv(fh, sep="\t", index=False, float_format="%.5g")
    with open(os.path.join(ES, "summary_by_type.tsv"), "w", encoding="utf-8") as fh:
        fh.write(f"# built {TODAY} by effect_sizes/build_figures.py from effects_by_study.tsv (headline studies, one per donor set)\n")
        summ.to_csv(fh, sep="\t", index=False, float_format="%.4g")

    # ---- figures
    fig1(hC, grp, d, os.path.join(ES, "fig_1_true_effect_sd_by_type.svg"))
    fig2(hA, A, os.path.join(ES, "fig_2_significant_gene_effects.svg"))
    fig3(hB, grp, os.path.join(ES, "fig_3_all_gene_effects.svg"))
    fig4(hC, d, os.path.join(ES, "fig_4_sd_vs_sample_size.svg"))

    # ---- stats.json
    def r4(x):
        return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)
    st = {"built": TODAY, "seed": SEED, "n_boot": N_BOOT, "sd_floor": SD_FLOOR,
          "summary_by_type": [{k: (r4(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()} for r in rows],
          "regressions": reg, "mann_whitney": mw, "uniform_pipelines": sens, "same_donor_pairs": pairs,
          "headline_C": hC[["label", "type", "N", "C_sd", "C_top_group", "C_top_sd", "full_provenance"]].to_dict("records"),
          "n_headline": {"A": int(len(hA)), "B": int(len(hB)), "C": int(len(hC))},
          "within_ipsc_sd_vs_n": within_ipsc, "t_likelihood_headline": t_lik,
          "endresz_donor_sets": endresz_one["label"].tolist(),
          "n_sources": {r: int((d["role"] == r).sum()) for r in d["role"].unique()}}
    json.dump(st, open(os.path.join(ES, "stats.json"), "w"), indent=1, default=lambda o: None if pd.isna(o) else o)

    # ---- captions and preview
    s = summ.set_index("type")
    rb = reg["brain_ref_rnaseq"]["coef"]
    caps = [
        ("fig_1_true_effect_sd_by_type.svg", "Figure 1. Noise-corrected spread of SZ vs control effects, by study type (view C).",
         "Each large point is one study (one per donor set; point size = cases + controls): the SD of the distribution of true "
         "log2 fold changes that ashr estimates from every tested gene's log2 FC and SE, taken as the median over cell types "
         "(single-nucleus) or regions; small light points are single-nucleus cell types; log scale; points on the dashed line have "
         "SD = 0. Panel b repeats the measure within two uniform pipelines (Gemma for bulk and brainSCOPE pseudobulk; Endresz 2026 "
         "for seven snRNA-seq cohorts). SE-derived likelihood: t (lik_t) where the SE comes from a t statistic or p value with "
         "df below 30, else normal (ash_fit.R).",
         f"Median SD: brain array {s.loc[ARRAY,'median_C_sd']:.3f} (n = {int(s.loc[ARRAY,'n_studies_C'])}), brain RNA-seq "
         f"{s.loc[RNASEQ,'median_C_sd']:.3f} (n = {int(s.loc[RNASEQ,'n_studies_C'])}), single-nucleus {s.loc[SC,'median_C_sd']:.3f} "
         f"(n = {int(s.loc[SC,'n_studies_C'])}); Endresz cohorts {sens['endresz_sc_cohorts']['median_C_sd']:.3f} "
         f"(n = {sens['endresz_sc_cohorts']['n']}; {sens['endresz_sc_donor_sets']['median_C_sd']:.3f} over "
         f"{sens['endresz_sc_donor_sets']['n']} donor sets, MSSM1 and MSSM2 counted once)."),
        ("fig_2_significant_gene_effects.svg", "Figure 2. |log2 FC| of the genes each study called significant (view A).",
         "One column per study (label with N), grouped by type and ordered by N; points are significant rows (at most 300 drawn), "
         "black bars the study median over all rows; single-cell studies stack cell types; log scale.",
         f"Median of study medians: brain array {s.loc[ARRAY,'median_A_abs']:.2f}, brain RNA-seq {s.loc[RNASEQ,'median_A_abs']:.2f}, "
         f"blood {s.loc[BLOOD,'median_A_abs']:.2f}, iPSC {s.loc[IPSC,'median_A_abs']:.2f}, single-nucleus {s.loc[SC,'median_A_abs']:.2f}; "
         "these are selection-inflated, most for small N."),
        ("fig_3_all_gene_effects.svg", "Figure 3. Raw |log2 FC| over all tested genes (view B).",
         "One point per study with a full table (size = N): (a) median and (b) 90th percentile of |log2 FC| over every tested gene, "
         "median over cell types or regions; light points are single-nucleus cell types; not corrected for sampling noise.",
         f"Median of study medians: brain array {s.loc[ARRAY,'median_B_median']:.3f}, brain RNA-seq {s.loc[RNASEQ,'median_B_median']:.3f}, "
         f"single-nucleus {s.loc[SC,'median_B_median']:.3f}; raw per-gene estimates are dominated by noise, which is largest for "
         "small N and for lowly expressed genes in single cell types."),
        ("fig_4_sd_vs_sample_size.svg", "Figure 4. Noise-corrected SD against sample size.",
         "Each filled point is a headline study (x: cases + controls, log scale; y: ashr SD, log scale), coloured by type; open circles "
         "are the seven Endresz snRNA-seq cohorts; points on the dashed line have SD = 0.",
         f"Brain studies, log10(SD) ~ type + log10(N): single-nucleus vs RNA-seq x{rb[f'type[{SC}]']['fold']:.2f} "
         f"(95% CI {rb[f'type[{SC}]']['fold_ci'][0]:.2f} to {rb[f'type[{SC}]']['fold_ci'][1]:.2f}, p = {rb[f'type[{SC}]']['p']:.2g}); "
         f"log10(N) slope {rb['log10(N)']['b']:.2f} (p = {rb['log10(N)']['p']:.2g})."),
    ]
    with open(os.path.join(ES, "CAPTIONS.md"), "w", encoding="utf-8") as fh:
        fh.write(f"# Captions: effect-size figures\n\nGenerated {TODAY} by `effect_sizes/build_figures.py` from the files it plots; "
                 "rerun the script to refresh. Three lines per figure: what it is, how to read it, what it shows.\n\n")
        for f, a, b, c in caps:
            fh.write(f"## {f}\n\n{a}\n\n{b}\n\n{c}\n\n")
    with open(os.path.join(ES, "index.html"), "w", encoding="utf-8") as fh:
        fh.write("<!doctype html><html><head><meta charset='utf-8'><title>Effect-size figures</title>"
                 "<style>body{font-family:Helvetica,Arial,sans-serif;font-size:10pt;max-width:1100px;margin:24px auto;color:#1f2937}"
                 "img{max-width:100%;border:1px solid #e5e7eb}p{color:#374151}</style></head><body>"
                 f"<h1 style='font-size:14pt;font-weight:normal'>SZ vs control effect sizes: figures (built {TODAY})</h1>")
        for f, a, b, c in caps:
            fh.write(f"<h2 style='font-size:11pt'>{a}</h2><img src='{f}'><p>{b}</p><p>{c}</p>")
        fh.write("</body></html>")
    print(summ.round(4).to_string())
    print(json.dumps(reg["brain_ref_rnaseq"], indent=1)[:1500])


# ------------------------------------------------------------------ report (from the written files only)
def write_report():
    """effect_sizes_report.md from stats.json, summary_by_type.tsv, effects_by_study.tsv, data/excluded.tsv and
    ../donor_overlap/overlap_summary.tsv. Edit this function, not the report."""
    st = json.load(open(os.path.join(ES, "stats.json")))
    summ = pd.read_csv(os.path.join(ES, "summary_by_type.tsv"), sep="\t", comment="#").set_index("type")
    eb = pd.read_csv(os.path.join(ES, "effects_by_study.tsv"), sep="\t", comment="#", dtype={"study": str})
    exc = pd.read_csv(os.path.join(DATA, "excluded.tsv"), sep="\t", comment="#", dtype={"study": str})
    ov = pd.read_csv(os.path.join(os.path.dirname(ES), "donor_overlap", "overlap_summary.tsv"), sep="\t", comment="#")
    shas = {f: sha(os.path.join(ES, f))[:12] for f in ("stats.json", "summary_by_type.tsv", "effects_by_study.tsv")}
    reg = st["regressions"]
    U = st["uniform_pipelines"]

    def f3(x):
        return "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"

    def ci(lo, hi):
        return "" if lo is None or (isinstance(lo, float) and np.isnan(lo)) else f"{lo:.3f} to {hi:.3f}"

    def cf(r, k):
        c = r["coef"].get(k)
        return "" if not c else f"{c['fold']:.2f} [{c['fold_ci'][0]:.2f}, {c['fold_ci'][1]:.2f}], {c['p']:.2g}"

    def pair(a, b):
        x = ov[((ov["unit_a"] == a) & (ov["unit_b"] == b)) | ((ov["unit_a"] == b) & (ov["unit_b"] == a))]
        return x.iloc[0] if len(x) else None

    s_ = lambda t, c: summ.loc[t, c] if t in summ.index else np.nan  # noqa: E731
    rb = reg["brain_ref_rnaseq"]
    sc_rows = [r for r in st["headline_C"] if r["type"] == SC]
    L = ["---", "type: detail", "bib: My Library", f"modified: {TODAY}", "---", "", "> [!note] Paul", ">", "",
         "# Effect sizes in SZ vs control transcriptomics: single-nucleus vs bulk", "",
         "## Executive summary", "",
         "- **Question.** Are SZ vs control effects smaller in single-nucleus data than in bulk tissue (Paul's hypothesis)?",
         f"- **Data.** One proposed primary list or full table per study and donor set; noise-corrected spread of true log2 fold "
         f"changes estimated per study with ashr (view C: {int(s_(SC, 'n_studies_C'))} single-nucleus, {int(s_(RNASEQ, 'n_studies_C'))} "
         f"bulk brain RNA-seq, {int(s_(ARRAY, 'n_studies_C'))} bulk brain microarray studies).",
         f"- **Result.** Median SD of true effects: single-nucleus {f3(s_(SC, 'median_C_sd'))}, bulk RNA-seq {f3(s_(RNASEQ, 'median_C_sd'))}, "
         f"microarray {f3(s_(ARRAY, 'median_C_sd'))}. Single-nucleus vs RNA-seq fold, adjusted for N: {cf(rb, 'type[' + SC + ']')} (fold [95% CI], p).",
         f"- **Same pipeline.** Endresz 2026's uniform reprocessing of single-nucleus cohorts gives a median of "
         f"{f3(U['endresz_sc_cohorts']['median_C_sd'])} over {U['endresz_sc_cohorts']['n']} cohorts.",
         f"- **Main caveat.** Only {int(s_(SC, 'n_studies_C'))} independent single-nucleus donor sets, so the comparison is inconclusive; "
         "the pipeline changes the single-nucleus estimate on the same donors.", "",
         "Paul, 2026-09-24: \"I hypothesize that single-cell data is generally seeing smaller magnitude effects than bulk tissue.\" "
         f"Generated {TODAY} by `build_figures.py` (`write_report`) from `stats.json` (sha256 {shas['stats.json']}), "
         f"`summary_by_type.tsv` ({shas['summary_by_type.tsv']}) and `effects_by_study.tsv` ({shas['effects_by_study.tsv']}); edit the "
         "script, not this file.", "", "Rerun:", "", "```",
         ".venv/bin/python effect_sizes/harmonize.py                       # data/viewA.tsv.gz, data/full/*, data/sources.tsv, data/excluded.tsv",
         "Rscript effect_sizes/ash_fit.R                                   # data/by_group.tsv (views B and C per study x cell type / region)",
         ".claude/skills/plotting/.venv/bin/python effect_sizes/build_figures.py   # tables, stats.json, figures, CAPTIONS.md, index.html, this report",
         "```", "", "Endresz 2026's per-cohort parquet is first converted with `tools/parquet_to_tsv.py` (command in the header of `harmonize.py`).", "",
         "## Results in brief", "",
         f"View C (noise-corrected SD of true SZ vs control log2 fold changes), median over headline studies, one per donor set: "
         f"single-nucleus {f3(s_(SC, 'median_C_sd'))} (n = {int(s_(SC, 'n_studies_C'))}, bootstrap 95% CI {ci(s_(SC, 'C_sd_ci_lo'), s_(SC, 'C_sd_ci_hi'))}), "
         f"bulk brain RNA-seq {f3(s_(RNASEQ, 'median_C_sd'))} (n = {int(s_(RNASEQ, 'n_studies_C'))}, CI {ci(s_(RNASEQ, 'C_sd_ci_lo'), s_(RNASEQ, 'C_sd_ci_hi'))}), "
         f"bulk brain microarray {f3(s_(ARRAY, 'median_C_sd'))} (n = {int(s_(ARRAY, 'n_studies_C'))}, CI {ci(s_(ARRAY, 'C_sd_ci_lo'), s_(ARRAY, 'C_sd_ci_hi'))}). "
         f"Regression of log10(SD) on type and log10(N) over brain studies (n = {rb['n']}): single-nucleus vs RNA-seq fold "
         f"{cf(rb, 'type[' + SC + ']')} (fold [95% CI], p). Endresz 2026's uniform reprocessing: median {f3(U['endresz_sc_cohorts']['median_C_sd'])} over "
         f"{U['endresz_sc_cohorts']['n']} cohorts (CI {ci(*U['endresz_sc_cohorts']['ci'])}), {f3(U['endresz_sc_donor_sets']['median_C_sd'])} over "
         f"{U['endresz_sc_donor_sets']['n']} donor sets (MSSM1 and MSSM2 share Mount Sinai donors and count once). With "
         f"{int(s_(SC, 'n_studies_C'))} independent single-nucleus donor sets the CI of the fold difference spans "
         f"{rb['coef']['type[' + SC + ']']['fold_ci'][0]:.2f} to {rb['coef']['type[' + SC + ']']['fold_ci'][1]:.2f}: the comparison is inconclusive.", "",
         "## Method", "",
         "**Studies.** One proposed primary list per study from `primary_lists.json`, plus the single-cell and organoid lists in "
         f"`manifests/batch9_sc.json` and `batch9b_sc.json`. Excluded rows are in `data/excluded.tsv` ({len(exc)} rows): meta-analyses and "
         "analyses pooling public datasets, imputed expression, other contrasts, Puvogel 2022 (MAST coefficients on pre-selected genes), "
         "Endresz 2026's meta-analysis (its per-cohort results are used instead), and de Jong 2012 and van Eijk 2015, whose only list is "
         "Gemma GSE38485 (processed data log-transformed twice; the Gemma GSE38485 analyses are also left out of the Gemma sensitivity set). "
         "Four consortium analyses that pool cohorts in one model (Gandal 2018 RNA-seq, Xia 2024, Clifton 2023, Dong 2024) are computed "
         "with role `pooled` and left out of the headline.", "",
         "**Donors, not samples.** Mendizabal 2019 (GSE107638: NeuN+ and OLIG2+ nuclei of 50 donors, 23 SZ and 27 control) and Arion 2017 "
         "(GSE87610: two arrays per donor per layer, 17 SZ and 19 control donors) were analyzed by Gemma with each sample as an "
         "observation and no donor term. Their N is set to donors, their full tables (views B and C) and their Gemma sensitivity "
         "entries are left out, and their lists stay in view A.", "",
         "**Types, harmonization, full tables.** As before: array vs RNA-seq from the universe type, then lineage platform text, then the "
         "manifest contrast (`type_evidence`); per-list conversion to log2 FC, SZ relative to control (`A_conversion`); full tables "
         "published > our reproduction > Gemma (`full_conversion`).", "",
         "**View C likelihood.** ashr per study and per cell type or region on (log2 FC, SE). SE from the table, else log2 FC / t, else "
         "|log2 FC| / qt(1 - p/2, df) (qnorm when df is unknown). Where the SE comes from a t statistic or p value with df below 30, the "
         "fit uses a t likelihood (`lik_t`) with that df and a uniform mixture; df is the table's (Student t tests in our reproductions) or, "
         "when a table has both t and p, the df at which 2 pt(-|t|, df) = p (limma moderated df for Gemma). Elsewhere a normal likelihood "
         "and normal mixture. Ishii 2019 (limma eBayes p) and Narla 2017 (Cuffdiff p) use qnorm. DESeq2 Wald statistics (Stern 2024, "
         "Huynh 2024) are referred to a normal distribution, so they keep the normal likelihood. Headline studies fitted with a t "
         "likelihood: " + ", ".join(f"{k} {v}" for k, v in st["t_likelihood_headline"].items()) + ".", "",
         "**One study per donor set.** Studies linked in `primary_lists.json` `same_data_as` form clusters, with no linking through "
         "meta-lane studies or multi-series reanalyses (the `recurrence.py` rule). Single-nucleus studies with shared donors are merged "
         "(`donor_overlap/`): Ruzicka 2024 and Emani 2024 (Emani's SZBDMulti-Seq is Ruzicka's McLean cohort; 76 of Ruzicka's 92 Mount Sinai "
         "donors match brainSCOPE CMC), and PsychAD and Emani 2024 (45 PsychAD Mount Sinai donors match brainSCOPE CMC, 19.3 by chance). "
         "The headline keeps one per cluster and view: published > reproduced > Gemma, then larger N. "
         f"Headline counts: {st['n_headline']['A']} studies in view A, {st['n_headline']['B']} in B, {st['n_headline']['C']} in C.", "",
         f"**Statistics.** Bootstrap 95% CI of the median over studies within type ({st['n_boot']:,} resamples, seed {st['seed']}); OLS of "
         "log10(SD) on type (reference bulk brain RNA-seq) and log10(N); studies with SD = 0 are dropped from the log regression and, as a "
         f"sensitivity, set to {st['sd_floor']}.", "",
         "## Results", "", "### View C: noise-corrected spread", "", "![[fig_1_true_effect_sd_by_type.svg]]", ""]
    caps = open(os.path.join(ES, "CAPTIONS.md"), encoding="utf-8").read().split("## ")[1:]
    capd = {c.split("\n", 1)[0].strip(): " ".join(x.strip() for x in c.split("\n", 1)[1].strip().split("\n\n")) for c in caps}
    L += [capd["fig_1_true_effect_sd_by_type.svg"], "",
          "| Type | n studies (C) | Median view-C SD [95% CI] | Median view-A \\|log2FC\\| (n) | Median view-B median \\|log2FC\\| (n) | Median N (C) |",
          "|---|---:|---|---:|---:|---:|"]
    for t in TYPES:
        if t not in summ.index:
            continue
        r = summ.loc[t]
        L.append(f"| {t} | {int(r['n_studies_C'])} | " + (f"{r['median_C_sd']:.3f} [{r['C_sd_ci_lo']:.3f}, {r['C_sd_ci_hi']:.3f}]" if r["n_studies_C"] > 1 else
                 (f"{r['median_C_sd']:.3f}" if r["n_studies_C"] == 1 else "")) +
                 " | " + ("" if pd.isna(r["median_A_abs"]) else "%.2f" % r["median_A_abs"]) + f" ({int(r['n_studies_A'])}) | "
                 + ("" if pd.isna(r["median_B_median"]) else "%.3f" % r["median_B_median"]) + f" ({int(r['n_studies_B'])}) | "
                 + ("" if pd.isna(r["median_N_C"]) else "%.0f" % r["median_N_C"]) + " |")
    mw = st["mann_whitney"]
    L += ["", "Single-nucleus headline studies: " + ", ".join(f"{r['label']} (N = {int(r['N'])}) {r['C_sd']:.3f}" for r in sc_rows) +
          f". Mann-Whitney, single-nucleus vs RNA-seq p = {mw['sc_vs_' + RNASEQ]['p']:.2g}, vs microarray p = {mw['sc_vs_' + ARRAY]['p']:.2g}.", "",
          "Uniform pipelines (Figure 1b): " + "; ".join(
              f"{lab} {f3(U[k]['median_C_sd'])} (n = {U[k]['n']}" + (f", CI {ci(*U[k]['ci'])})" if U[k]["n"] > 1 else ")")
              for k, lab in (("gemma_array", "Gemma brain arrays"), ("gemma_rnaseq", "Gemma brain RNA-seq"),
                             ("gemma_sc", "Gemma brainSCOPE pseudobulk"), ("endresz_sc_cohorts", "Endresz snRNA-seq cohorts"),
                             ("endresz_sc_donor_sets", "Endresz donor sets"))) + ".", "",
          "Same donors, different pipeline (`stats.json`, `same_donor_pairs`): " + "; ".join(
              f"{v['published']} published {v['published_C_sd']:.3f} vs " + (f"Endresz {k} {v['endresz_C_sd']:.3f}" if "endresz_C_sd" in v
                                                                             else f"{k.replace('gemma__', 'Gemma ')} {v['gemma_C_sd']:.3f}")
              for k, v in st["same_donor_pairs"].items()) + ".", "",
          "Per-cell-type maxima (`C_top_sd` in `effects_by_study.tsv`) are listed but not compared with per-study medians: the maximum of "
          "19 to 35 noisy per-cell-type SDs has no null here, so it is not evidence that effects are concentrated in a few cell types.", "",
          "### Regression on sample size", "", "![[fig_4_sd_vs_sample_size.svg]]", "", capd["fig_4_sd_vs_sample_size.svg"], "",
          "| Model (log10 SD ~ type + log10 N; reference bulk brain RNA-seq) | n | Single-nucleus vs RNA-seq, fold [95% CI], p | Microarray vs RNA-seq, fold [95% CI], p | Per 10-fold N, fold [95% CI], p |",
          "|---|---:|---|---|---|"]
    for k, lab in (("brain_ref_rnaseq", "Brain headline studies"), ("brain_endresz_cohorts_ref_rnaseq", "Brain, single-nucleus replaced by Endresz donor sets"),
                   ("gemma_brain_ref_rnaseq", "Brain, Gemma analyses only"), ("brain_sd0_at_floor_ref_rnaseq", f"Brain, SD = 0 set to {st['sd_floor']}"),
                   ("all_types_ref_rnaseq", "All types"), ("all_types_without_ipsc_ref_rnaseq", "All types without iPSC/organoid")):
        r = reg[k]
        L.append(f"| {lab} | {r['n']} | {cf(r, 'type[' + SC + ']')} | {cf(r, 'type[' + ARRAY + ']')} | {cf(r, 'log10(N)')} |")
    wi = st["within_ipsc_sd_vs_n"]
    L += ["", f"Within the iPSC/organoid studies, Spearman correlation of log SD with log N: rho = {wi['spearman_rho']:.2f} (n = {wi['n']}, "
          f"p = {wi['p']:.2g}).", "",
          "### View A: significant genes", "", "![[fig_2_significant_gene_effects.svg]]", "", capd["fig_2_significant_gene_effects.svg"], "",
          "### View B: all tested genes, raw", "", "![[fig_3_all_gene_effects.svg]]", "", capd["fig_3_all_gene_effects.svg"], ""]
    hl = pair("hbcc_sgacc:HBCC", "libd_benjamin2024:LIBD")
    jin = eb[(eb["study"] == "34957416")]
    L += ["### Figure 5 (bulk vs single-nucleus on the same donors): not made", "",
          "Emani 2024's CMC donors are CommonMind Mount Sinai donors. Jin 2021 (CMC-SCZ, "
          + (f"N = {int(jin['N'].iloc[0])}" if len(jin) and not pd.isna(jin["N"].iloc[0]) else "N not stated")
          + ", DESeq2 bulk table) is a headline bulk RNA-seq study with a full table on CommonMind donors, which include Mount Sinai; "
          "whether its donors include the brainSCOPE 100 was not checked (`donor_overlap/` has no CommonMind bulk donor table: the "
          "clinical file is under Synapse controlled access). Liu 2021's 254-donor Mount Sinai table is in the same position. PsychAD's "
          "HBCC donors overlap the Akula 2021 / Smith 2026 sgACC donors, a different region. The same-donor comparison that is available "
          "is between pipelines, above.", "",
          "## Caveats", "",
          f"- **Few single-nucleus donor sets.** {int(s_(SC, 'n_studies_C'))} in the headline after merging studies that share donors; "
          f"{U['endresz_sc_donor_sets']['n']} Endresz donor sets, which overlap the published studies' donors.",
          "- **Shared donors across bulk headline studies.** "
          + (f"Smith 2026 (HBCC sgACC) matches {int(hl['exact_shared'])} LIBD BrainSeq donors exactly ({hl['exact_chance_mean']} by chance; "
             f"{int(hl['no_pmi_shared'])} more at the weaker no-PMI level) in `donor_overlap/`, so HBCC and LIBD headline studies are not "
             "fully independent." if hl is not None else ""),
          "- **Composition.** A bulk log2 FC mixes within-cell changes with shifts in cell-type proportions; a pseudobulk log2 FC within "
          "one cell type does not.",
          f"- **Array compression.** Microarray log2 FCs are compressed relative to RNA-seq for the same change; medians here: array "
          f"{f3(s_(ARRAY, 'median_C_sd'))}, RNA-seq {f3(s_(RNASEQ, 'median_C_sd'))}; Gemma arrays {f3(U['gemma_array']['median_C_sd'])} vs Gemma "
          f"RNA-seq {f3(U['gemma_rnaseq']['median_C_sd'])}.",
          "- **Pseudobulk normalization and model.** The pipeline moves the single-nucleus SD on the same donors (table of same-donor "
          "pairs above). Which pipeline is closer to the truth is not something these data decide; the published tables are the "
          "headline because that is what the papers report.",
          "- **ashr SD depends on power.** When no gene clears the noise, ashr puts all weight on the null (SD = 0), so small or noisy "
          "studies can read as zero rather than small. SEs derived from p values assume the stated test and df.",
          "- **Units.** Conversions are per list (`A_conversion`, `full_conversion`). Bharadwaj 2024 (transcript level) and Stankovic 2024 "
          "(Seurat avg_log2FC up to 526) contain extreme values from zero-count groups; they affect view A maxima, not medians.",
          "- **Contrasts.** Fröhlich 2024 compares psychiatric cases (36 of 54 SZ) with controls; Ketchesin 2023's psychosis group includes "
          "8 BD with psychosis; some iPSC studies are co-twin or 22q11.2 designs.",
          "- **Single-nucleus lists without full tables** (Reiner 2020, Gerstner 2025, Wu 2026) enter view A only.", "",
          "## Files", "",
          "- `effects_by_study.tsv`: one row per source (headline, pooled, Gemma dataset, Endresz cohort) with N, views A to C, conversions and provenance header.",
          "- `summary_by_type.tsv`, `stats.json`: the numbers above.",
          "- `data/by_group.tsv`: views B and C per study x cell type or region, with the likelihood used (`lik`, `df_lik`).",
          "- `data/excluded.tsv`: every study, list or Gemma dataset left out, and why.",
          "- `CAPTIONS.md`, `index.html`: captions and a preview page.", ""]
    with open(os.path.join(ES, "effect_sizes_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    if "--report" in sys.argv:
        write_report()
    else:
        main()
        write_report()
