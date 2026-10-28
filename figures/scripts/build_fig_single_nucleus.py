#!/usr/bin/env python3
"""Extended Fig. 8 -- single-nucleus data do not agree better than bulk and carry the same pH signature.

Data: figures/data/fig09/ (2026-10-08).

a  agreement (Spearman rho of t) between subject-separate cohorts against donors per group:
   192 bulk dataset pairs and 19 single-nucleus cohort pairs (Endresz 2026, median over 14 classes),
   with bin medians (0-10, 11-30, 31+ donors per group)
b  PsychAD at 20 donors per group, per cell class: HBCC vs Mount Sinai for SZ and age, and split
   halves within each bank for SZ; medians over the 8 classes
c  per Endresz subclass (24): Pearson r between each gene's bulk control-pH association and its
   single-nucleus SZ effect, over the 1,431-gene 5+ lists universe and genome-wide
d  noise-corrected SD of true SZ log2 fold changes (ashr) by data type, brain studies only
"""
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

import figkit as fk

D = fk.FD / "fig09"
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
rng = np.random.default_rng(3)
SN = fk.TEXT  # single-nucleus series (bulk = ACCENT)


def panel_a(ax):
    pairs = rd("fig09_a_agreement_pairs.tsv")
    bins = rd("fig09_a_agreement_by_n_bin.tsv")
    bulk = pairs[pairs["data_type"].str.startswith("bulk")]
    sn = pairs[pairs["data_type"].str.startswith("single")]
    # single nucleus: median over the 14 classes within each cohort pair (as in the bin file)
    snp = sn.groupby(["unit_a", "unit_b"], as_index=False).agg(n_per_group=("n_per_group", "first"),
                                                             rho_t=("rho_t", "median"))
    jit = lambda x: x * np.exp(rng.normal(0, 0.03, len(x)))
    ax.scatter(jit(bulk["n_per_group"]), bulk["rho_t"], s=5, color=fk.ACCENT, alpha=0.55, lw=0,
               label=f"bulk tissue, one pair of datasets ({len(bulk)})")
    ax.scatter(jit(snp["n_per_group"]), snp["rho_t"], s=14, facecolor="white", edgecolor=SN, lw=0.8,
               label=f"single nucleus, one pair of cohorts ({len(snp)}),\nmedian over 14 cell classes")
    spans = {"0-10": (4.5, 10.5), "11-30": (10.5, 30.5), "31+": (30.5, 75)}
    for _, r in bins.iterrows():
        x0, x1 = spans[r["donors_per_group_bin"]]
        col = fk.ACCENT if r["data_type"].startswith("bulk") else SN
        ax.plot([x0, x1], [r["median_rho_t"]] * 2, color=col, lw=1.6, solid_capstyle="butt")
    ax.set_xscale("log")
    ax.set_xticks([5, 10, 20, 30, 50])
    ax.set_xticklabels(["5", "10", "20", "30", "50"])
    ax.set_xlim(4.5, 80)
    ax.set_ylim(-0.42, 0.72)
    ax.axhline(0, color=fk.GRID, lw=0.6)
    ax.set_title("Single-nucleus cohorts agree\nno better than bulk datasets", pad=4)
    ax.set_xlabel("Donors per group (smaller group\nin the pair, log scale)")
    ax.set_ylabel("Agreement between cohorts\n(Spearman ρ of SZ t statistics)")
    ax.legend(loc="upper left", fontsize=fk.SMALL, handletextpad=0.3, borderaxespad=0.2)
    ax.text(0.98, 0.03, "one point = one pair of cohorts; bars = bin medians", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=fk.SMALL, color=fk.SUBTLE)
    return bins


def panel_b(ax):
    df = rd("fig09_b_psychad_20_per_group.tsv")
    med = rd("fig09_b_medians_over_8_classes.tsv")
    df = df[df["cell_type"] != "total"]
    series = [  # (contrast, comparison, label, colour, filled)
        ("sz", "HBCC_vs_MSSM_main", "SZ, HBCC vs Mount Sinai", fk.ACCENT, True),
        ("age", "HBCC_vs_MSSM_main", "age, HBCC vs Mount Sinai", fk.ACCENT_2, True),
        ("sz", "split_HBCC_main", "SZ, two halves of HBCC", fk.ACCENT, False),
        ("sz", "split_MSSM_main", "SZ, two halves of Mount Sinai", fk.lighten(fk.ACCENT, 0.5), False),
    ]
    order = {"EN": "excitatory", "IN": "interneuron", "Astro": "astrocyte", "Oligo": "oligodendrocyte",
             "OPC": "OPC", "Immune": "immune", "Endo": "endothelial", "Mural": "mural"}
    classes = list(order)
    ypos = {c: i for i, c in enumerate(classes)}
    ymed = len(classes) + 0.6
    handles = []
    for k, (con, comp, lab, col, filled) in enumerate(series):
        sub = df[(df["contrast"] == con) & (df["comparison"] == comp)]
        y = [ypos[c] for c in sub["cell_type"]]
        kw = dict(s=16, lw=0.8, zorder=3)
        if filled:
            ax.scatter(sub["rho_t"], y, color=col, edgecolor="white", **kw)
            h = Line2D([], [], marker="o", color=col, ls="", ms=4.5)
        else:
            ax.scatter(sub["rho_t"], y, facecolor="white", edgecolor=col, **kw)
            h = Line2D([], [], marker="o", mfc="white", mec=col, ls="", ms=4.5)
        handles.append((h, lab))
        m = med[(med["contrast"] == con) & (med["comparison"] == comp)]["median_over_classes"].iloc[0]
        ax.plot([m, m], [ymed - 0.35, ymed + 0.35], color=col, lw=1.6, zorder=3)
    ax.set_yticks(list(range(len(classes))) + [ymed])
    ax.set_yticklabels([order[c] for c in classes] + ["median"], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.axhline(len(classes) - 0.3, color=fk.GRID, lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Between PsychAD banks,\nSZ agrees less than age", pad=4)
    ax.set_xlabel("Agreement at 20 donors per group\n(Spearman ρ of t statistics, median of 20 draws)\none point = one PsychAD cell class")
    ax.set_xlim(-0.15, 0.62)
    ax.legend([h for h, _ in handles], [l for _, l in handles], loc="lower right", bbox_to_anchor=(1.0, 0.1), fontsize=fk.TINY,
              handletextpad=0.2, borderaxespad=0.2, labelspacing=0.3)
    return med


def panel_d(ax):
    st = rd("fig09_d_effect_sd_by_study.tsv")
    ct = rd("fig09_d_effect_sd_single_nucleus_cell_types.tsv")
    en = rd("fig09_d_endresz_cohorts.tsv")
    summ = rd("fig09_d_summary_by_type.tsv").set_index("type")
    groups = ["bulk brain microarray", "bulk brain RNA-seq", "single-nucleus brain", "Endresz cohorts"]
    labels = ["bulk\narray", "bulk\nRNA-seq", "sn\nstudies", "sn, one\npipeline"]
    floor = 3e-3
    ax.axhline(floor, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.text(-0.55, floor * 1.12, "SD = 0", fontsize=fk.SMALL, color=fk.SUBTLE, va="bottom", ha="left")
    meds = {}
    ns = {}
    for i, g in enumerate(groups):
        if g == "Endresz cohorts":
            sub = en.rename(columns={"cohort": "label"})
            meds[g] = sub["C_sd"].median()
        else:
            sub = st[st["type"] == g]
            meds[g] = summ.loc[g, "median_C_sd"]
        ns[g] = len(sub)
        x = i + rng.uniform(-0.18, 0.18, len(sub))
        yv = sub["C_sd"].clip(lower=floor)
        size = 4 + 2.2 * np.log10(sub["N"].clip(lower=3))**2 * 2
        col = fk.ACCENT if g.startswith("bulk") else SN
        if g == "Endresz cohorts":
            ax.scatter(x, yv, s=size, facecolor="white", edgecolor=col, lw=0.8, zorder=3)
        else:
            ax.scatter(x, yv, s=size, color=col, alpha=0.75, edgecolor="white", lw=0.4, zorder=3)
        if g == "single-nucleus brain":
            xc = i + rng.uniform(-0.3, 0.3, len(ct))
            ax.scatter(xc, ct["sd_g"].clip(lower=floor), s=2.5, color=SN, alpha=0.35, lw=0, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [meds[g]] * 2, color=fk.TEXT, lw=1.6, zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(floor * 0.8, 1.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{l}\n({ns[g]})" for l, g in zip(labels, groups)], fontsize=fk.SMALL)
    ax.set_xlabel("Data type (number of studies)")
    ax.set_xlim(-0.6, len(groups) - 0.2)
    ax.set_title("True SZ effects are no larger\nin single-nucleus data", pad=4)
    ax.set_ylabel("Spread of true SZ effects\n(s.d. of log2 FC, ashr; log scale)")
    ax.text(0.02, 0.03, "large point = one study;\nsmall points = cell types", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=fk.TINY, color=fk.SUBTLE)
    for s, lab in [(4 + 2.2 * 4 * 1, "N = 10"), (4 + 2.2 * 4 * 2, "N = 100"), (4 + 2.2 * 9 * 2, "N = 1,000")]:
        ax.scatter([], [], s=s, color=fk.NEUTRAL, label=lab)
    ax.legend(loc="upper right", fontsize=fk.SMALL, handletextpad=0.2, borderaxespad=0.2, labelspacing=0.6)
    return meds, ns


def main():
    fig = fk.figure(88)
    gs = GridSpec(1, 3, figure=fig, left=0.085, right=0.995, top=0.86, bottom=0.2, wspace=0.5,
                  width_ratios=[1.1, 1, 1])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axd = fig.add_subplot(gs[0, 2])
    bins = panel_a(axa)
    med = panel_b(axb)
    meds, ns = panel_d(axd)
    for ax, s in zip([axa, axb, axd], "abc"):
        fk.letter(ax, s, dy=0.004)
    notes = (f"a bin medians: {bins[['donors_per_group_bin','data_type','n_pairs','median_rho_t']].to_dict('records')}; "
             f"b medians over 8 classes: {med[med['comparison'].str.endswith('main')][['comparison','contrast','median_over_classes']].to_dict('records')}; "
             f"c (was d) medians {meds}, n {ns}")
    fk.save(fig, "fig_single_nucleus", "figures/scripts/build_fig_single_nucleus.py", sorted(D.glob("*.tsv")), notes)
    print(notes)


if __name__ == "__main__":
    main()
