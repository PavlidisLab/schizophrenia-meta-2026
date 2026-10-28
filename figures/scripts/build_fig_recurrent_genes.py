#!/usr/bin/env python3
"""Supplementary Figure 6 (the former Fig. 4) -- the recurrent genes: direction, sources, co-occurrence, DE prior.

Data: figures/data/fig06/ (export_fig06.py, 2026-10-08).

Extended Fig. 3 as drawn: a per recurrent gene, lists reporting it up / down / without direction;
b share of list hits in each gene's majority direction, observed vs random directions, four gene sets;
c AUC of the DE prior for the recurrent gene sets.
Not drawn: strata, per-list AUROC, co-occurrence.
"""
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

import figkit as fk

D = fk.FD / "fig06"
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
_SIZES = pd.read_csv(D / "fig06_b_direction_summary.tsv", sep="\t", comment="#").set_index("gene_set")["genes"]
NG = {k: f"{int(v):,}" for k, v in _SIZES.items()}  # gene-set sizes for labels, from the export
# the recurrent set is the 375 (study tier, q < 0.1)
REC_B = "study tier q<0.1 (uniform null)"
N_LISTS = len(pd.read_csv(fk.FD / "fig02" / "fig02_b_study_annotation.tsv", sep="\t", comment="#"))
rng = np.random.default_rng(1)


def panel_a(ax_up, ax_down):
    df = rd("fig06_a_gene_direction.tsv")
    df["net"] = df["n_up"] - df["n_down"]
    rows = [(ax_up, df[df["net"] > 0].sort_values(["net", "n_studies"], ascending=[False, False])),
            (ax_down, df[df["net"] <= 0].sort_values(["net", "n_studies"], ascending=[True, False]))]
    width = max(len(d) for _, d in rows)
    for k, (ax, d) in enumerate(rows):
        d = d.reset_index(drop=True)
        x = np.arange(len(d))
        ax.bar(x, d["n_up"], 0.8, color=fk.UP, lw=0, label="up in SZ")
        ax.bar(x, -d["n_down"], 0.8, color=fk.DOWN, lw=0, label="down in SZ")
        other = d["n_mixed"] + d["n_nodir"]
        ax.bar(x, other, 0.8, bottom=d["n_up"], color=fk.NEUTRAL, lw=0, label="mixed or no direction")
        ax.axhline(0, color=fk.TEXT, lw=0.5)
        ax.set_xlim(-0.8, width - 0.2)
        ax.set_ylim(-14, 17)
        ax.set_yticks([-10, -5, 0, 5, 10, 15])
        if width > 150:  # too many names to read at this size: hide them, as the lab heatmap helper does
            ax.set_xticks([])
        else:
            ax.set_xticks(x)
            ax.set_xticklabels(d["gene"], rotation=90, fontsize=fk.TINY)
            ax.tick_params(axis="x", length=0, pad=1)
        ax.set_ylabel("Studies reporting the gene\n(up above, down below)", fontsize=fk.SMALL)
        ax.yaxis.set_major_formatter(lambda v, _: f"{abs(v):.0f}")
        ax.text(0.995, 0.97, f"mostly {'up' if k == 0 else 'down'} ({len(d)} genes)", transform=ax.transAxes,
                ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)
    ax_up.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02), ncol=3, handlelength=0.8, fontsize=fk.SMALL, columnspacing=0.8)
    ax_up.set_title(f"The {len(df)} recurrent genes and the direction each study reports (one bar = one gene)", pad=4)


def panel_b(ax):
    df = rd("fig06_b_direction_summary.tsv")
    labels = [f"recurrent\n({NG[REC_B]})", f"recurrent beyond\nexpression level\n({NG['expression-matched q<0.1']})",
              f"10+ lists\n({NG['10 or more studies']})", f"3+ lists\n({NG['3 or more studies']})"]
    order = [REC_B, "expression-matched q<0.1", "10 or more studies", "3 or more studies"]
    df = df.set_index("gene_set").loc[order]
    y = np.arange(len(df))
    for yi, (_, r) in zip(y, df.iterrows()):
        ax.plot([r["share_majority_if_random"], r["share_majority"]], [yi, yi], color=fk.GRID, lw=2, zorder=1)
    ax.scatter(df["share_majority_if_random"], y, s=14, color=fk.NEUTRAL, zorder=2, label="random directions")
    ax.scatter(df["share_majority"], y, s=14, color=fk.ACCENT, zorder=3, label="observed")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fk.TINY)
    ax.set_xlim(0.55, 0.95)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.set_xlabel("Hits in the gene's majority direction\n(share of all list hits in the set)")
    ax.set_title("Hits agree on direction\nfar beyond chance\n(row = gene set)", pad=4)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)
    r0 = df.iloc[0]
    ax.text(r0["share_majority_if_random"], -0.3, "random\ndirections", ha="left", va="bottom", fontsize=fk.SMALL, color=fk.SUBTLE)  # left-anchored: clear of the axis line
    ax.text(r0["share_majority"], -0.3, "observed", ha="center", va="bottom", fontsize=fk.SMALL, color=fk.ACCENT)
    ax.set_ylim(len(df) - 0.5, -1.1)


STRATA = [("size", "large (N >= 150)", "large studies (N ≥ 150)"),
          ("era", "era 2001-2012", "2001–2012"), ("era", "era 2013-2018", "2013–2018"), ("era", "era 2019-2026", "2019–2026"),
          ("tech", "tech array_early", "early arrays"), ("tech", "tech array_whole", "whole-genome arrays"),
          ("tech", "tech rnaseq", "RNA-seq"), ("tech", "tech single_cell", "single-cell"),
          ("theme", "immune / inflammation", "immune / inflammation"),
          ("theme", "synapse / neurotransmission", "synapse"),
          ("theme", "neurodevelopment", "neurodevelopment"),
          ("theme", "mitochondria / energy metabolism", "mitochondria / energy"),
          ("theme", "GABA interneurons", "GABA interneurons"),
          ("theme", "epigenetics / chromatin", "epigenetics"),
          ("theme", "RNA processing / non-coding RNA", "RNA processing"),
          ("theme", "cell death / apoptosis", "cell death"),
          ("theme", "oligodendrocytes / myelin", "oligodendrocytes / myelin"),
          ("theme", "protein homeostasis / ubiquitin", "protein homeostasis"),
          ("theme", "astrocytes", "astrocytes")]


def panel_c(ax):
    df = rd("fig06_c_strata.tsv")
    df = df[df["gene_set"] == "study tier q<0.1"].set_index(["stratum_type", "stratum"])
    rows = [(t, s, l) for t, s, l in STRATA if (t, s) in df.index]
    y, ypos, prev = [], 0, None
    for t, _, _ in rows:
        if prev is not None and t != prev:
            ypos += 0.6
        y.append(ypos)
        ypos += 1
        prev = t
    for yi, (t, s, _) in zip(y, rows):
        r = df.loc[(t, s)]
        sig = pd.notna(r["q_bh"]) and r["q_bh"] < 0.05
        ax.scatter(r["ratio"], yi, s=14, color=fk.ACCENT_3 if t == "theme" else fk.ACCENT, zorder=3,
                   edgecolor=fk.TEXT if sig else "none", lw=0.6)
        ax.text(1.03, yi, f"{int(r['n_studies'])}", fontsize=fk.SMALL, color=fk.SUBTLE, va="center", ha="left",
                transform=ax.get_yaxis_transform())
    ax.axvline(1, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xscale("log", base=2)
    ax.set_xlim(0.5, 4.2)
    ax.set_xticks([0.5, 1, 2, 4])
    ax.set_xticklabels(["½", "1", "2", "4"])
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, _, l in rows], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.set_xlabel("Over- or under-representation of the recurrent genes\n(share of their hits / share of all list entries; 1 = proportional)")
    ax.set_title("No kind of study carries the recurrence\n(one point = one group of studies)", pad=4)
    ax.text(1.03, -1.0, "studies", fontsize=fk.SMALL, color=fk.SUBTLE, ha="left", va="center", transform=ax.get_yaxis_transform())
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)


def panel_e(ax):
    df = rd("fig06_e_de_prior_auc_by_set.tsv").set_index(["gene_set", "reference"])
    sets = [("study_brain (209-equivalent, q<0.1)", f"recurrent\n({NG[REC_B]})"),
            ("expression-matched q<0.1", f"beyond\nexpression\nlevel\n({NG['expression-matched q<0.1']})"),
            ("study_brain_ge10 (primary, 10+ lists)", f"10+ lists\n({NG['10 or more studies']})")]
    refs = [("global prior, n_de_total", "prior from\nall tissues", fk.ACCENT),
            ("brain postmortem, no SZ banks prior, n_de_total", "prior from brain,\nSZ banks left out", fk.ACCENT_3)]
    x = np.arange(len(sets))
    for k, (ref, lab, col) in enumerate(refs):
        v = [df.loc[(s, ref), "auc"] for s, _ in sets]
        ax.bar(x + (k - 0.5) * 0.36, v, 0.34, color=col, label=lab)
        for xi, vi in zip(x, v):
            ax.text(xi + (k - 0.5) * 0.36, vi + 0.01, f"{vi:.2f}", ha="center", va="bottom", fontsize=fk.TINY, color=fk.TEXT, rotation=90)
    ax.axhline(0.5, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_ylim(0.4, 1.3)
    ax.set_yticks(np.arange(0.4, 1.01, 0.1))
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in sets], fontsize=fk.TINY)
    ax.set_ylabel("DE prior picks out the set\n(AUC; 0.5 = chance)")
    ax.set_title("The DE prior picks out\nthe recurrent genes\n(bar = gene set)", pad=4)
    ax.legend(loc="upper left", handlelength=0.8, fontsize=fk.TINY, ncol=1)


def panel_f(ax):
    df = rd("fig06_e_de_prior_per_list_auroc.tsv")
    df = df[(df["reference"] == "global prior, n_de_total") & (df["subset"] == "all lists")]
    groups = [("SZ published brain", "SZ"), ("BD Gemma top 166", "BD"), ("MDD Gemma top 166", "MDD")]
    df = df.set_index("group")
    for i, (g, _) in enumerate(groups):
        r = df.loc[g]
        ax.plot([i, i], [r["q1"], r["q3"]], color=fk.ACCENT if i == 0 else fk.NEUTRAL, lw=2.5)
        ax.scatter(i, r["median"], s=16, color=fk.TEXT, zorder=3)
    ax.axhline(0.5, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{l}\n{int(df.loc[g, 'n_lists'])}" for g, l in groups], fontsize=fk.SMALL)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylabel("Fit of the DE prior to each list\n(per-list AUROC; 0.5 = chance)")
    ax.set_title("DE prior fits\nSZ lists only", pad=4)
    ax.text(0.0, -0.2, "point = median\nbar = IQR\nnumbers = lists", transform=ax.transAxes,
            ha="left", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)


def main():
    # Extended Fig. 3: a: recurrent genes by direction; b: majority-direction share; c: AUC of the DE prior.
    # The strata panel (panel_c, hits by study group) was dropped; panel_c and panel_f are kept unused.
    fig = fk.figure(185)
    g1 = GridSpec(2, 1, figure=fig, hspace=0.72, left=0.075, right=0.99, top=0.96, bottom=0.58)
    g2 = GridSpec(1, 2, figure=fig, width_ratios=[1, 1], wspace=0.45, left=0.12, right=0.97, top=0.4, bottom=0.11)
    axa, axa2 = fig.add_subplot(g1[0]), fig.add_subplot(g1[1])
    axb, axe = (fig.add_subplot(g2[i]) for i in range(2))
    panel_a(axa, axa2)
    panel_b(axb)
    panel_e(axe)
    for ax, l, dx in [(axa, "a", -0.02), (axb, "b", -0.02), (axe, "c", -0.02)]:
        fk.letter(ax, l, dx=dx)
    fk.save(fig, "supp_recurrent_genes", "figures/scripts/build_fig_recurrent_genes.py", sorted(f for f in D.glob("*.tsv") if "cooccurrence" not in f.name),
            notes="Extended Fig. 3: a direction per gene, b direction share, c DE-prior AUC (strata panel dropped 2026-10-11)")


if __name__ == "__main__":
    main()
