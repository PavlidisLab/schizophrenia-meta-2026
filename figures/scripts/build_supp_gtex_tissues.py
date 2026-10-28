#!/usr/bin/env python3
"""Supplementary Fig. 9 -- a slow death in every GTEx tissue: like brain's slow-death effect, like the schizophrenia
genes, but like low brain pH only in brain.

Data: figures/data/gtex_death_speed_all_tissues/ (analyses/gtex_death_speed_all_tissues_2026-10-10/gtex_death_speed_all_tissues.py).
Per tissue, the slow (Hardy 4) against fast natural (Hardy 2) death t of every gene (OLS with ischemic time, age, sex),
then the Spearman rho of those t with
a  the brain slow-death axis (mean t of Cortex, BA9, BA24; circular for those three regions)
b  the control pH axis, sign flipped so that positive means genes move the way they do at low pH
c  the schizophrenia direction of the recurrent genes
One row is one tissue with >= 10 donors in each group (rho_by_tissue_and_contrast.tsv, contrast slow_4_vs_fast_2).
d  per tissue, rho with the schizophrenia direction of the recurrent genes for Hardy 3 vs 2 (length of dying) against Hardy 4 vs 3
   (long illness) (rho_hardy3_contrasts_by_tissue.tsv; 36 tissues with >= 10 Hardy 3 donors)
e  share of the agreement between the 192 bulk schizophrenia dataset pairs removed by each brain axis
   (pair_agreement_hardy3_axis_removal.tsv, set ALL, own genes)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figkit as fk

D = fk.FD / "gtex_death_speed_all_tissues"
F = D / "rho_by_tissue_and_contrast.tsv"
F3 = D / "rho_hardy3_contrasts_by_tissue.tsv"
FP = D / "pair_agreement_hardy3_axis_removal.tsv"
AXIS_REGIONS = {"Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)"}


def load():
    d = pd.read_csv(F, sep="\t", comment="#")
    d = d[d["contrast"] == "slow_4_vs_fast_2"].copy()
    d["brain"] = d["SMTS"] == "Brain"
    d["low_ph"] = -d["rho_pH_Z"]
    d = d.sort_values(["brain", "rho_brain_slow_death_axis"], ascending=[True, True])
    return d.reset_index(drop=True)


def panel(ax, d, col, title, xlabel, show_labels):
    y = np.arange(len(d))
    colours = np.where(d["brain"], fk.ACCENT, fk.NEUTRAL)
    ax.hlines(y, 0, d[col], color=colours, lw=0.8)
    ax.scatter(d[col], y, s=10, c=colours, zorder=3, lw=0)
    if col == "rho_brain_slow_death_axis":
        circ = d["tissue"].isin(AXIS_REGIONS)
        ax.scatter(d.loc[circ, col], y[circ], s=26, facecolors="none", edgecolors=fk.TEXT, lw=0.6, zorder=4)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6)
    nb = int((~d["brain"]).sum())
    ax.axhline(nb - 0.5, color=fk.SUBTLE, lw=0.4, ls="--")
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(d["tissue"].str.replace("Brain - ", "", regex=False) if show_labels else [], fontsize=fk.SMALL)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)


def panel_d(ax):
    d = pd.read_csv(F3, sep="\t", comment="#")
    w = d.pivot_table(index=["tissue", "SMTS"], columns="contrast", values="rho_SZ_direction_176").reset_index()
    brain = w["SMTS"] == "Brain"
    for m, col, lab in [(~brain, fk.NEUTRAL, f"other tissues ({int((~brain).sum())})"), (brain, fk.ACCENT, f"brain ({int(brain.sum())})")]:
        ax.scatter(w.loc[m, "hardy4_vs_hardy3"], w.loc[m, "hardy3_vs_hardy2"], s=12, color=col, lw=0, label=lab)
    lim = (-0.4, 0.9)
    ax.plot(lim, lim, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.axhline(0, color=fk.GRID, lw=0.6); ax.axvline(0, color=fk.GRID, lw=0.6)
    ax.set_xlabel("Long illness: Hardy 4 vs 3\n(rho with schizophrenia direction)")
    ax.set_ylabel("Length of dying: Hardy 3 vs 2\n(rho with schizophrenia direction)")
    ax.set_title(f"The {fk.n_recurrent()} genes correlate more with the length\nof dying than with long illness", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="lower right")
    ax.text(0.02, 0.98, "one point = one tissue", transform=ax.transAxes, fontsize=fk.SMALL, color=fk.SUBTLE, va="top")


def panel_e(ax):
    p = pd.read_csv(FP, sep="\t", comment="#")
    p = p[(p["set"] == "ALL") & (p["geneset"] == "own")]
    axes = [("4v2", "Hardy 4 vs 2\n(agonal axis)"), ("3v2", "Hardy 3 vs 2\n(length of dying)"), ("4v3", "Hardy 4 vs 3\n(long illness)")]
    x = np.arange(len(axes)); w = 0.36
    for j, (lev, lab, col) in enumerate([("L0", "age and sex", fk.ACCENT), ("L4", "all measured covariates", fk.NEUTRAL)]):
        ys = [p[(p["axes_removed"] == a) & (p["level"] == lev)]["share_removed"].iloc[0] for a, _ in axes]
        ax.bar(x + (j - 0.5) * w, ys, w, color=col, label=lab)
        for xi, yi in zip(x + (j - 0.5) * w, ys):
            ax.text(xi, yi + 0.01, f"{yi:.0%}", ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT)
    ax.set_xticks(x); ax.set_xticklabels([a[1] for a in axes], fontsize=fk.SMALL)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_ylim(0, 0.75)
    ax.set_ylabel("Agreement between schizophrenia\ndatasets removed (192 pairs)")
    ax.set_title("With age and sex only, the length-of-dying\naxis takes out the most agreement", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="upper right", ncol=2, title="models adjusted for", title_fontsize=fk.SMALL)


def main():
    d = load()
    fig = fk.figure(235)
    axs = [fig.add_axes(r) for r in ([0.30, 0.37, 0.20, 0.58], [0.54, 0.37, 0.20, 0.58], [0.78, 0.37, 0.20, 0.58])]
    panel(axs[0], d, "rho_brain_slow_death_axis", "Like a slow death\nin brain",
          "Agreement with the brain\nslow-death effect (Spearman rho)", True)
    panel(axs[1], d, "low_ph", "Like low brain pH\nonly in brain",
          "Agreement with the low-pH\npattern (Spearman rho)", False)
    panel(axs[2], d, "rho_SZ_direction_176", f"Like the schizophrenia\ndirection, everywhere",
          f"Agreement with schizophrenia\ndirection, {fk.n_recurrent()} genes (rho)", False)
    axs[1].text(0.98, 0.005, "blue, 13 brain regions\ngrey, 32 other tissues", transform=axs[1].transAxes, ha="right",
                va="bottom", fontsize=fk.SMALL, color=fk.SUBTLE)
    axs[2].text(0.98, 0.005, "one point = one tissue", transform=axs[2].transAxes, ha="right", va="bottom",
                fontsize=fk.SMALL, color=fk.SUBTLE)
    ad = fig.add_axes([0.12, 0.05, 0.33, 0.20]); panel_d(ad)
    ae = fig.add_axes([0.62, 0.05, 0.34, 0.20]); panel_e(ae)
    for ax, l in zip(axs + [ad, ae], "abcde"):
        fk.letter(ax, l)
    fk.save(fig, "supp_gtex_tissues", "figures/scripts/build_supp_gtex_tissues.py", [F, F3, FP],
            notes="Open circles in a: the "
                  "three regions the brain axis is built from (circular).")


if __name__ == "__main__":
    main()
