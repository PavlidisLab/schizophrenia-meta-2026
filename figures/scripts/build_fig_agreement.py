#!/usr/bin/env python3
"""Figure 2 -- agreement between published brain studies is near chance.

Data: figures/data/fig02/ (64-study population,
analyses/figure_exports/export_fig02_fig07.py, 2026-10-08).

a  genes in >= k of the primary lists, observed vs three nulls (+-1.96 s.d.); corrected null from analyses/code001_nk_correction_2026-10-11 (CODE-001)
b  meta-analysis pairs: overlap over the DE prior against share of input datasets
c  gene-by-study calls for the 80 most recurrent genes, studies (rows) grouped by cohort cluster, with strips for each
   study's cohort cluster and sample size (cases plus controls, fig01/fig01_b_corpus_studies.tsv); right, each study's genes tested and genes reported DE (log; n_tested from the fig02 annotation)
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.colors import LogNorm

import os
import figkit as fk

D = fk.FD / "fig02"
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
rng = np.random.default_rng(1)


def jitter(n, w=0.28):
    return rng.uniform(-w, w, n)


# Panel a reads the corrected recurrence curve (CODE-001: each list's expected size equals its size, n_k = |M_k|).
CODE001_CURVE = fk.FD / "fig02" / "fig02_a_recurrence_curve.tsv"


def panel_a(ax):
    df = pd.read_csv(CODE001_CURVE, sep="\t", comment="#").rename(columns={
        "n_studies_at_least": "k", "observed_genes": "observed",
        "expected_uniform": "exp_uniform", "expected_de_prior": "exp_deprior", "expected_expression": "exp_exprmatched",
        "sd_de_prior": "sd_deprior", "sd_expression": "sd_exprmatched"})
    k = df["k"].to_numpy()
    for key in ("uniform", "deprior", "exprmatched"):
        m, s = df[f"exp_{key}"].to_numpy(), df[f"sd_{key}"].to_numpy()
        col = fk.NULL_COLOURS[key]
        ax.fill_between(k, np.clip(m - 1.96 * s, 0.1, None), m + 1.96 * s, color=fk.lighten(col, 0.8), lw=0)
        ax.plot(k, np.where(m >= 0.1, m, np.nan), "--", color=col, lw=0.9, label=fk.NULL_LABELS[key])
    ax.plot(k, np.where(df["observed"] > 0, df["observed"], np.nan), "-o", color=fk.TEXT, ms=2.5, lw=1.1, label="observed")  # a zero count has no place on a log axis: the curve ends at the last k with a gene
    ax.set_yscale("log")
    ax.set_ylim(0.3, 3e4)
    ax.set_xlim(0.5, k.max() + 0.5)
    ax.set_xticks([1, 4, 8, 12, 16])
    ax.set_title("Genes by the number of lists naming them", pad=4)
    ax.set_xlabel(f"Lists naming the gene (at least k of {fk.n_primary()})")
    ax.set_ylabel("Genes recurring (count, log scale)")
    ax.legend(loc="upper right", handlelength=1.2, fontsize=fk.SMALL, labelspacing=0.3, borderaxespad=0.2)  # empty corner above the curves


def panel_b(ax):
    df = rd("fig02_c_pairwise_fold_by_era.tsv")
    df = df[df["shares_any_gene"] == "yes"]
    groups = ["old-old", "new-new", "old-new"]
    pooled = rd("fig02_c_pooled_by_era.tsv")
    pooled = pooled[pooled["era_split"] == pooled["era_split"].iloc[0]].set_index("pair_type")
    for i, g in enumerate(groups):
        y = df.loc[df["era_pair"] == g, "fold"].to_numpy()
        ax.scatter(i + jitter(len(y)), y, s=1.5, color=fk.lighten(fk.ACCENT, 0.45), lw=0, rasterized=True)
        ax.plot([i - 0.32, i + 0.32], [np.median(y)] * 2, color=fk.TEXT, lw=1.2)
    # Gemma reanalysis: the same 36 pairs, published lists vs uniform reanalysis
    gp = rd("fig02_c_gemma_reanalysis_pairs.tsv").pivot(index="pair_id", columns="method", values="fold").dropna()
    x0, x1 = 3.3, 4.8
    for _, r in gp.iterrows():
        ax.plot([x0, x1], [r["published"], r["gemma_sized"]], color=fk.lighten(fk.SUBTLE, 0.55), lw=0.4)
    pool = rd("fig02_c_gemma_reanalysis_pooled.tsv").set_index("method")
    ax.plot([x0, x1], [pool.loc["published", "estimate"], pool.loc["gemma_sized", "estimate"]],
            "-o", color=fk.TEXT, ms=2.5, lw=1.2)
    ax.axvline(2.95, color=fk.GRID, lw=0.6)
    ax.set_xticks([0, 1, 2, x0, x1])
    ax.set_xticklabels([f"{lab}\n{int(pooled.loc[g, 'n_pairs'])}" for g, lab in zip(groups, ["older", "newer", "mixed"])]
                       + ["published", "Gemma"], fontsize=fk.SMALL)
    ax.text((x0 + x1) / 2, 1.02, f"same {len(gp)} pairs", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=fk.SMALL, color=fk.SUBTLE)
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_yscale("log")
    ax.set_title("Pairs of studies overlap barely above chance", pad=12)
    ax.set_ylabel("Pair overlap (observed /\nexpected by chance)")
    ax.set_xlabel("Era of the two studies (pairs)")
    ax.text(0.02, 0.97, "one point = one pair of studies", transform=ax.transAxes, ha="left", va="top",
            fontsize=fk.SMALL, color=fk.SUBTLE)
    ax.set_xlim(-0.6, 5.4)


def panel_c(ax):
    df = rd("fig02_d_hoffman82_by_group.tsv")
    labels = {"other banks": "Other\nbanks", "shared bank": "Shared\nbank", "CMC-HBCC reanalyses": "Same\ndonors"}
    x = np.arange(len(df))
    w = 0.26
    for off, col, colour, lab in [(-w, "observed", fk.ACCENT, "observed"),
                                  (0, "exp_uniform", fk.NEUTRAL, "random genes"),
                                  (w, "exp_deprior", fk.ACCENT_3, "DE-prior-weighted\nrandom genes")]:
        ax.bar(x + off, df[col], w, color=colour, label=lab)
    for xi, (_, r) in zip(x, df.iterrows()):
        ax.text(xi - w, r["observed"] + 6, f"{int(r['observed'])}", ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{labels[g]}\n({n})" for g, n in zip(df["group"], df["n_studies"])])
    ax.set_title("Hoffman et al.'s 82 genes recur only in the same donors", pad=4)
    ax.set_ylabel("Hits among the 82 genes\n(count over later studies)")
    ax.set_xlabel("Relation of later studies to the CommonMind cohorts (studies)")
    ax.legend(loc="upper left", handlelength=1)


def panel_d(ax):
    df = rd("fig02_e_meta_analysis_pairs.tsv")
    df = df[df["fold_deprior"] > 0]
    same = df["same_technology"] == "yes"
    x = df["dataset_jaccard"] + jitter(len(df), 0.012)
    ax.scatter(x[same], df.loc[same, "fold_deprior"], s=5, color=fk.ACCENT, lw=0, label="same technology")
    ax.scatter(x[~same], df.loc[~same, "fold_deprior"], s=5, color=fk.ACCENT_3, lw=0, label="different technology")
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_yscale("log")
    ax.set_title("Meta-analyses agree in proportion to shared inputs", pad=4)
    ax.set_xlabel("Shared input datasets (Jaccard index)")
    ax.set_ylabel("List overlap (observed /\nDE-prior expectation)")
    ax.legend(loc="lower right", handlelength=0.8, markerscale=1.2,
              title="one point = one pair of meta-analyses", title_fontsize=fk.SMALL)


META = fk.FD / "meta15"


def panel_meta(ax):
    """AUROC of the 383 recurrent genes, ranked by the number of lists containing each gene: one row per meta-analysis list
    (membership of that list) and, at the top, all 15 lists together."""
    pl = pd.read_csv(META / "per_list_auroc.tsv", sep="\t", comment="#")
    pl = pl[pl["set"] == 383].copy()
    po = pd.read_csv(META / "pooled_auroc.tsv", sep="\t", comment="#")
    po = po[(po["variant"] == "all15") & (po["set"] == 383) & (po["ranking_universe"] == "tested_by_at_least_one_list")].iloc[0]
    ic = pd.read_csv(META / "meta15_input_counts.tsv", sep="\t", comment="#").set_index("quantity")["value"]
    n_distinct = int(ic["distinct pooled datasets"])  # studies re-used across meta-analyses counted once
    ns = pd.read_csv(META / "meta15_input_coverage.tsv", sep="\t", comment="#").set_index("meta_analysis")["n_input_datasets"]
    lg = pd.read_csv(META / "meta15_per_meta_analysis_hits.tsv", sep="\t", comment="#")
    lg = lg[lg["set"] == "study_brain_383"].set_index("meta_analysis")["list_genes"]
    assert set(pl["meta_analysis"]) <= set(ns.index)
    pl["lab"] = [f"{m}   {int(ns[m])}/{int(lg[m]):,}" for m in pl["meta_analysis"]]
    pl = pl.sort_values("auroc_observed")
    rows = [(r.lab, r.auroc_observed, r.null_deprior_lo95, r.null_deprior_hi95, r.null_uniform_mean) for r in pl.itertuples()]
    n_genes = int(pd.read_csv(META / "meta15_genes_by_k.tsv", sep="\t", comment="#").query("variant == 'all15' and k_min == 1")["genes_in_k_or_more"].iloc[0])
    rows.append((f"All 15 pooled   {n_distinct} studies / {n_genes:,} genes", po["auroc_observed"], po["null_deprior_lo95"], po["null_deprior_hi95"], po["null_uniform_mean"]))
    bs = pd.read_csv(META / "pooled_auroc_bootstrap.tsv", sep="\t", comment="#")
    bs = bs[(bs["set"] == 383) & bs["method"].str.startswith("3 dataset bootstrap (list weight")].iloc[0]
    pool_ci = (bs["p2_5"], bs["p97_5"])  # resampling the 42 input datasets; each list weighted by the share of its datasets drawn
    y = np.arange(len(rows))
    for yy, (_, obs, lo, hi, uni) in zip(y, rows):
        ax.add_patch(plt.Rectangle((lo, yy - 0.32), hi - lo, 0.64, color=fk.GRID, lw=0, zorder=1))
        ax.plot([uni, uni], [yy - 0.32, yy + 0.32], color=fk.NULL_COLOURS["uniform"], lw=1.0, zorder=2)
        ax.scatter([obs], [yy], s=16, color=fk.TEXT if yy == y[-1] else fk.ACCENT, zorder=3)
    ax.plot(pool_ci, [y[-1], y[-1]], color=fk.TEXT, lw=1.0, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=fk.SMALL)
    ax.axhline(y[-1] - 0.5, color=fk.GRID, lw=0.6)
    ax.axvline(0.5, color=fk.GRID, lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.45, None)
    ax.set_xlabel("Recovery of the 383 genes (AUROC)")
    ax.set_title("The recurrent genes in\npublished meta-analyses", pad=4)
    ax.text(0.99, 0.02, "dot, observed; grey, 95% range expected from the\nDE prior; line, expected from random genes",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=fk.TINY, color=fk.SUBTLE)
    return len(rows) - 1


CALLS = ["not_tested", "tested_not_hit", "up", "down", "mixed", "hit_no_direction"]
CALL_COLOURS = ["white", "#f3f4f6", fk.UP, fk.DOWN, "#7c3aed", "#4b5563"]
CALL_LABELS = ["not measured", "measured, not listed", "up in SZ (hit list)", "down in SZ (hit list)", "mixed", "no direction stated"]


N_GENES = 80  # was 60; widened to 80 with a narrower per-study bar
# panel e row order: "staircase" (the oncoprint "memo sort": studies listing gene 1 first, then those listing gene 2
# among the rest, and so on; within a step, up before down; ties by list size) or "cohort" (grouped by cohort cluster)
ROW_ORDER = "staircase"
TINT = float(os.environ.get("SZ_TINT", "0.05")) or None  # |log2 fold change| from which unlisted cells are coloured
SOLID = bool(os.environ.get("SZ_SOLID"))  # test: colour unlisted cells at full strength
EFFECT_HOLD = set()  # studies not tinted until checked (Maas 2020 checked 2026-10-19: Group2 is schizophrenia)
TINT_FULL = 0.5  # |log2 fold change| at which the colour is at full strength


def tinted_matrix(m, studies, genes, ann):
    """RGB matrix: the call colours, with 'measured, not listed' cells tinted by the study's effect where one is held."""
    from matplotlib.colors import to_rgb
    base = np.array([to_rgb(c) for c in CALL_COLOURS])
    rgb = base[m.astype(int)]
    ef = pd.read_csv(fk.FD / "fig02" / "fig02_f_effects_by_study.tsv.gz", sep="\t", comment="#",
                     usecols=["study_label", "gene", "effect", "is_fig2c_group"])
    ef = ef[~ef["study_label"].isin(EFFECT_HOLD)]
    ef = ef[ef["is_fig2c_group"].astype(str).str.upper() == "TRUE"].groupby(["study_label", "gene"])["effect"].mean()
    lab = dict(zip(ann["study_id"], ann["label"]))
    up, down = np.array(to_rgb(fk.UP)), np.array(to_rgb(fk.DOWN))
    for i, sid in enumerate(studies):
        L = lab[sid]
        for j, g in enumerate(genes):
            if m[i, j] != 1 or (L, g) not in ef.index:
                continue
            e = ef[(L, g)]
            if not np.isfinite(e) or abs(e) < TINT:
                continue
            a = 0.08 + 0.14 * min(abs(e) / TINT_FULL, 1.0)
            rgb[i, j] = (up if e > 0 else down) if SOLID else (1 - a) * base[1] + a * (up if e > 0 else down)
    return rgb
PH = [("strong", "low-pH direction, q < 0.05", fk.ACCENT_3), ("weak", "low-pH direction", fk.lighten(fk.ACCENT_3, 0.6)),
      ("opposite", "opposite direction", fk.NEUTRAL), ("na", "not measured", "white")]


AG = [("strong", "slow-death direction, q < 0.05", "#a21caf"), ("weak", "slow-death direction", fk.lighten("#a21caf", 0.6)),
      ("opposite", "opposite direction", fk.NEUTRAL), ("na", "not measured", "white")]


def ag_class(r):
    if pd.isna(r["agonal_consistent"]):
        return "na"
    if str(r["agonal_consistent"]).upper() != "TRUE":
        return "opposite"
    return "strong" if r["q_agonal"] < 0.05 else "weak"


def ph_class(r):
    if pd.isna(r["ph_confound_consistent"]):
        return "na"
    if str(r["ph_confound_consistent"]) != "True":
        return "opposite"
    return "strong" if r["q_ph_controls"] < 0.05 else "weak"


CLUSTERS = [("cmc+mssm+pitt+libd", "CMC, Mount Sinai, Pitt, LIBD", fk.ACCENT_2), ("stanley", "Stanley", "#0891b2"),
            ("aus", "NSW / Victoria", "#8b5cf6"), ("harvard", "Harvard / McLean", "#ec4899"),
            ("pritzker", "Pritzker", "#84cc16"), ("single", "one-study cluster", fk.NEUTRAL)]


def cluster_group(r):
    return r["cohort_cluster"] if r["cohort_cluster"] in dict((k, 1) for k, _, _ in CLUSTERS) else "single"


def panel_e(fig, spec):
    long = rd("fig02_b_gene_by_study_long.tsv")
    ann = pd.read_csv(D / "fig02_b_study_annotation.tsv", sep="\t", comment="#", dtype={"study_id": str})
    ann["grp"] = ann.apply(cluster_group, axis=1)
    # sample size: corpus table keyed by PMID, by "sc:"+PMID for single-nucleus studies, else by label (DOI-keyed studies)
    corp = pd.read_csv(fk.FD / "fig01" / "fig01_b_corpus_studies.tsv", sep="\t", comment="#", dtype={"study_id": str})
    corp = corp[corp["in_primary_population"] == "yes"]
    by_id = corp.set_index("study_id")["n_total"]
    by_label = corp.drop_duplicates("label").set_index("label")["n_total"]
    norm_label = lambda l: l.replace(" (PsychAD)", "").replace("Frohlich", "Fröhlich")
    ann["n_total"] = [by_id.get(i, by_id.get("sc:" + i, by_label.get(norm_label(l)))) for i, l in zip(ann["study_id"], ann["label"])]
    gorder = {k: i for i, (k, _, _) in enumerate(CLUSTERS)}
    ann = ann.assign(g=ann["grp"].map(gorder)).sort_values(["g", "column_order"], kind="stable").reset_index(drop=True)
    top = long.drop_duplicates("gene").sort_values("gene_rank").head(N_GENES)
    genes = top["gene"].tolist()
    if ROW_ORDER == "staircase":
        rank = {"up": 3, "down": 2, "mixed": 1, "hit_no_direction": 1}
        calls = long.assign(study_id=long["study_id"].astype(str), r=long["call"].map(rank).fillna(0))
        R = calls.pivot_table(index="study_id", columns="gene", values="r", aggfunc="max").reindex(
            index=ann["study_id"], columns=genes).fillna(0).to_numpy()
        keys = [tuple(-R[i]) + (-ann["n_list_genes"].iloc[i],) for i in range(len(ann))]
        ann = ann.iloc[sorted(range(len(ann)), key=lambda i: keys[i])].reset_index(drop=True)
    studies = ann["study_id"].tolist()
    code = {c: i for i, c in enumerate(CALLS)}
    long = long.assign(study_id=long["study_id"].astype(str))
    m = long.assign(c=long["call"].map(code)).pivot(index="study_id", columns="gene", values="c")
    m = m.reindex(index=studies, columns=genes).fillna(0).to_numpy()
    gs = spec.subgridspec(5, 4, height_ratios=[1.3, 0.16, 0.16, 0.16, 5], width_ratios=[0.022, 0.022, 1, 0.10], hspace=0.06,
                          wspace=0.02)
    axn, axp, axg, axi, ax, axr = (fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[2, 2]),
                                   fig.add_subplot(gs[3, 2]), fig.add_subplot(gs[4, 2]), fig.add_subplot(gs[4, 3]))
    axc, axz = fig.add_subplot(gs[4, 0]), fig.add_subplot(gs[4, 1])
    x = np.arange(len(genes))
    # lists per gene, split by the direction each study reports (up above the line, down below)
    dirc = long[long["gene"].isin(genes)].pivot_table(index="gene", columns="call", values="study_id", aggfunc="count").reindex(genes).fillna(0)
    n_up, n_down = dirc.get("up", 0 * dirc.index.size), dirc.get("down", 0)
    n_other = dirc.get("mixed", 0) + dirc.get("hit_no_direction", 0)
    axn.bar(x, n_up, 0.75, color=fk.UP, lw=0)
    axn.bar(x, n_other, 0.75, bottom=n_up, color=fk.NEUTRAL, lw=0)
    axn.bar(x, -n_down, 0.75, color=fk.DOWN, lw=0)
    axn.axhline(0, color=fk.TEXT, lw=0.6)
    axn.set_xlim(-0.5, len(genes) - 0.5)
    axn.set_xticks([])
    axn.set_ylabel("Studies\n(up above,\ndown below)", fontsize=fk.SMALL)
    axn.set_title(f"The {N_GENES} most recurrent genes and the direction each study reports", pad=62)
    axn.set_yticks([-10, -5, 0, 5, 10, 15])
    axn.set_yticklabels(["10", "5", "0", "5", "10", "15"])
    axn.spines["bottom"].set_visible(False)
    # control-pH association of each gene
    ph = pd.read_csv(fk.FD / "fig07" / "fig07_a2_a4_gene_scatter.tsv", sep="\t", comment="#")
    ph = top[["gene"]].merge(ph, on="gene", how="left")
    cls = ph.apply(ph_class, axis=1)
    pcode = {k: i for i, (k, _, _) in enumerate(PH)}
    axp.imshow(cls.map(pcode).to_numpy()[None, :], aspect="auto", interpolation="nearest",
               cmap=ListedColormap([c for _, _, c in PH]), vmin=-0.5, vmax=len(PH) - 0.5)
    axp.set_xticks([])
    axp.set_yticks([0])
    axp.set_yticklabels(["pH assoc."], fontsize=fk.SMALL)
    axp.tick_params(axis="y", length=0)
    axp.grid(False)
    for sp in axp.spines.values():
        sp.set_visible(False)
    # agonal association of each gene (GTEx slow against fast death)
    ag = pd.read_csv(fk.FD / "fig02" / "fig02_b_gene_agonal_axis.tsv", sep="\t", comment="#")
    ag = top[["gene"]].merge(ag, on="gene", how="left")
    acode = {k: i for i, (k, _, _) in enumerate(AG)}
    axg.imshow(ag.apply(ag_class, axis=1).map(acode).to_numpy()[None, :], aspect="auto", interpolation="nearest",
               cmap=ListedColormap([c for _, _, c in AG]), vmin=-0.5, vmax=len(AG) - 0.5)
    axg.set_xticks([])
    axg.set_yticks([0])
    axg.set_yticklabels(["agonal assoc."], fontsize=fk.SMALL)
    axg.tick_params(axis="y", length=0)
    axg.grid(False)
    for sp in axg.spines.values():
        sp.set_visible(False)
    # relative specificity of each gene to interneurons (share of summed class mean CPM, PsychAD HBCC controls)
    iv = pd.read_csv(fk.FD / "fig02" / "fig02_b_gene_interneuron_specificity.tsv", sep="\t", comment="#")
    iv = top[["gene"]].merge(iv[["gene", "spec_coarse_IN"]], on="gene", how="left")
    axi.imshow(iv["spec_coarse_IN"].to_numpy(dtype=float)[None, :], aspect="auto", interpolation="nearest",
               cmap=LinearSegmentedColormap.from_list("in_spec", ["#ffffff", "#0f766e"]), vmin=0, vmax=0.6)
    axi.set_xticks([])
    axi.set_yticks([0])
    axi.set_yticklabels(["interneuron"], fontsize=fk.SMALL)
    axi.tick_params(axis="y", length=0)
    axi.grid(False)
    for sp in axi.spines.values():
        sp.set_visible(False)
    # matrix
    if TINT is None:
        ax.imshow(m, aspect="auto", interpolation="nearest", cmap=ListedColormap(CALL_COLOURS),
                  vmin=-0.5, vmax=len(CALLS) - 0.5)
    else:
        ax.imshow(tinted_matrix(m, studies, genes, ann), aspect="auto", interpolation="nearest")
    ax.grid(False)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=90, fontsize=fk.TINY)
    # gene symbols in the top 10% of the DE prior (our global prior) in colour
    prior_f = D / "fig02_b_gene_de_prior.tsv"
    if prior_f.exists():
        pr = pd.read_csv(prior_f, sep="\t", comment="#").set_index("gene")["top10pct"].astype(str).eq("yes")
        marked = [i for i, g in enumerate(genes) if pr.get(g, False)]
        ax.scatter(marked, [0.0] * len(marked), s=7, color=fk.TEXT, marker="o", lw=0, clip_on=False,
                   transform=ax.get_xaxis_transform() + ScaledTranslation(0, -4 / 72, fig.dpi_scale_trans))
        n_top = int(sum(pr.get(g, False) for g in genes))
        print(f"green top-10% DE prior: {n_top} of {len(genes)}")
    else:
        n_top = None
    ax.tick_params(axis="x", length=0, pad=7)
    ax.set_yticks([])
    bounds = np.flatnonzero(ann["grp"].to_numpy()[1:] != ann["grp"].to_numpy()[:-1]) + 0.5 if ROW_ORDER == "cohort" else []
    for b in bounds:
        for a_ in (ax, axr, axc, axz):
            a_.axhline(b, color="white", lw=1.2, zorder=3)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # each study's cohort cluster and sample size
    gcode = {k: i for i, (k, _, _) in enumerate(CLUSTERS)}
    axc.imshow(ann["grp"].map(gcode).to_numpy()[:, None], aspect="auto", interpolation="nearest",
               cmap=ListedColormap([c for _, _, c in CLUSTERS]), vmin=-0.5, vmax=len(CLUSTERS) - 0.5)
    axc.set_xticks([])
    axc.set_yticks([])
    for sp in axc.spines.values():
        sp.set_visible(False)
    axc.set_ylabel("Studies, grouped by cohort cluster" if ROW_ORDER == "cohort" else "Studies, ordered by the first gene they list",
                   fontsize=fk.SMALL, labelpad=2)
    nz = ann["n_total"].astype(float).to_numpy()
    zim = axz.imshow(nz[:, None], aspect="auto", interpolation="nearest", cmap="Greys", norm=LogNorm(10, 600))
    axz.set_xticks([])
    axz.set_yticks([])
    for sp in axz.spines.values():
        sp.set_visible(False)
    for a_, t in ((axc, "cohort cluster"), (axz, "N")):
        a_.set_xticks([0])
        a_.set_xticklabels([t], rotation=90, fontsize=fk.TINY)
        a_.tick_params(axis="x", length=0, pad=1)
    n_grp = ann["grp"].value_counts()
    chand = [Patch(facecolor=c, lw=0, label=f"{l} ({n_grp.get(k, 0)})") for k, l, c in CLUSTERS]
    ax.legend(handles=chand, title="Cohort cluster", loc="upper left", bbox_to_anchor=(0, -0.2), ncol=3,
              handlelength=1, columnspacing=1, fontsize=fk.SMALL, title_fontsize=fk.SMALL, alignment="left")
    cax = ax.inset_axes([0.78, -0.3, 0.2, 0.03])
    cb = fig.colorbar(zim, cax=cax, orientation="horizontal")
    cb.set_ticks([10, 30, 100, 300])
    cb.set_ticklabels(["10", "30", "100", "300"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=fk.TINY, length=2)
    cb.outline.set_visible(False)
    cax.set_title("N, cases plus controls (log)", fontsize=fk.SMALL, pad=2, loc="left")
    zim.cmap.set_bad(fk.GRID)
    # per study: genes tested (light) with the genes it reports as DE drawn over them (dark), log scale;
    # DE genes are a subset of the tested genes, so the dark bar sits inside the light one
    y = np.arange(len(studies))
    tested = pd.to_numeric(ann.get("n_tested"), errors="coerce").to_numpy() if "n_tested" in ann else np.full(len(studies), np.nan)
    axr.barh(y, np.nan_to_num(tested, nan=0), 0.8, color=fk.NEUTRAL, lw=0)
    axr.barh(y, ann["n_list_genes"].clip(lower=1).to_numpy(), 0.8, color=fk.TEXT, lw=0)
    for yy in y[np.isnan(tested)]:  # tested set not known
        axr.text(1.3, yy, "?", ha="left", va="center", fontsize=fk.TINY, color=fk.SUBTLE)
    axr.set_xscale("log")
    axr.set_xlim(1, 5e4)
    axr.set_xticks([1, 10, 100, 1000, 10000])
    axr.set_xticklabels(["1", "10", "100", "1k", "10k"])
    axr.xaxis.set_minor_formatter(plt.NullFormatter())
    axr.set_ylim(len(studies) - 0.5, -0.5)
    axr.set_yticks([])
    axr.grid(axis="y", visible=False)
    axr.grid(axis="x", visible=True)
    axr.tick_params(axis="x", labelsize=fk.TINY)
    for t_ in (10, 100, 1000, 10000):
        axr.axvline(t_, color="white", lw=0.8, zorder=5)
    axr.set_xlabel("Genes per\nstudy (log)", fontsize=fk.SMALL)
    axr.legend(handles=[Patch(facecolor=fk.NEUTRAL, lw=0, label="tested"), Patch(facecolor=fk.TEXT, lw=0, label="reported DE")],
               loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=fk.TINY, handlelength=1, borderaxespad=0.2, frameon=False)
    axr.spines["left"].set_visible(False)
    from matplotlib.colors import to_rgb
    P = lambda c, l: Patch(facecolor=c, edgecolor=fk.GRID, lw=0.4, label=l)
    head = lambda l: Patch(facecolor="none", edgecolor="none", label=l)
    cc = dict(zip(CALLS, CALL_COLOURS))
    base1 = np.array(to_rgb(cc["tested_not_hit"]))
    tint = lambda col: 0.86 * base1 + 0.14 * np.array(to_rgb(col))
    handles = [head("Not on the hit list:"), P(cc["not_tested"], "gene not measured"), P(cc["tested_not_hit"], "measured")]
    if TINT is not None:
        handles += [P(tint(fk.UP), f"measured, higher in SZ ($\\log_2$ FC ≥ {TINT:g})"), P(tint(fk.DOWN), f"measured, lower in SZ ($\\log_2$ FC ≥ {TINT:g})")]
    handles += [head("On the hit list:"), P(cc["up"], "up in SZ"), P(cc["down"], "down in SZ"), P(cc["mixed"], "mixed"), P(cc["hit_no_direction"], "no direction stated")]
    handles += [P(c, f"pH: {l.replace('low-pH direction', 'low pH')}") for _, l, c in PH[:3]]
    handles += [P(c, f"agonal: {l.replace('slow-death direction', 'slow death')}") for _, l, c in AG[:2]]
    handles += [P("#0f766e", "interneuron specificity (0 to 0.6)")]
    if n_top is not None:
        handles.append(Line2D([], [], marker="o", color=fk.TEXT, lw=0, markersize=3, label="dot: top 10% of DE prior"))
    while len(handles) % 4:  # four columns of equal length
        handles.append(head(""))
    axn.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.05), ncol=4, handlelength=1,
               columnspacing=0.6, fontsize=fk.SMALL)
    counts = cls.value_counts().to_dict()
    return axn, len(genes), len(studies), counts, int(ann["n_total"].isna().sum())


def main():
    # panel_b() and panel_c() are kept for build_supp_study_pairs.py
    fig = fk.figure(210)
    r1 = GridSpec(1, 2, figure=fig, width_ratios=[1, 0.9], wspace=0.95, left=0.085, right=0.985, top=0.95, bottom=0.73)
    r3 = GridSpec(1, 1, figure=fig, left=0.045, right=0.99, top=0.565, bottom=0.165)
    axa, axb = fig.add_subplot(r1[0]), fig.add_subplot(r1[1])
    panel_a(axa)
    n_meta = panel_meta(axb)
    axe, n_genes, n_studies, counts, n_none = panel_e(fig, r3[0])
    for ax, l in [(axa, "a"), (axb, "b")]:
        fk.letter(ax, l)
    fk.letter(axe, "c")
    inputs = sorted(D.glob("*.tsv")) + [CODE001_CURVE, fk.FD / "fig07" / "fig07_a2_a4_gene_scatter.tsv", fk.FD / "fig01" / "fig01_b_corpus_studies.tsv", fk.FD / "fig02" / "fig02_b_gene_agonal_axis.tsv", fk.FD / "fig02" / "fig02_b_gene_interneuron_specificity.tsv"]
    fk.save(fig, ("fig_agreement" if TINT == 0.05 else f"fig_agreement_tint{TINT:g}") if not SOLID else f"fig_agreement_solid{TINT:g}", "figures/scripts/build_fig_agreement.py", inputs,
            notes=f"panel c: {n_genes} genes x {n_studies} studies, rows ordered {ROW_ORDER} ({n_none} studies without N); pH classes {counts}")


if __name__ == "__main__":
    main()
