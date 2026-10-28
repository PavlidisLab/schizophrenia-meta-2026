"""Supplementary Fig. 3: positive control (age) beside diagnosis, and age against external brain data.

Two panels in one row. Data: figures/data/fig05/.
a) 192 subject-separate bulk pairs: rho(t) for diagnosis, age, autosomal sex (full and matched group sizes),
   plus top-200 overlap as log2 fold over chance (full data). Pairs with zero overlap (fold 0, blank log2) are
   drawn as hollow points at the axis floor and counted in n.
b) Per-dataset median rho(t): SZ vs SZ and age vs age within the 22 SZ-study dataset-regions, SZ-study age vs
   9 external frontal-cortex datasets (one point per dataset on either side, 31), and external vs external age (9).
   Black bars are medians over subject-separate pairs (summary.tsv), the numbers quoted in the text.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figkit as fk  # noqa: E402

D = fk.FD / "fig05"
F_A = D / "fig05_a_pairs_dx_age_sex.tsv"
F_B = D / "fig05_b_age_external_per_dataset.tsv"
F_BS = D / "fig05_b_age_external_summary.tsv"

TRAIT_COL = {"dx": fk.ACCENT, "age": fk.ACCENT_2, "sex": fk.NEUTRAL}
TRAIT_LAB = {"dx": "diagnosis", "age": "age", "sex": "sex"}


def jitter(n, w=0.22, seed=0):
    return np.random.default_rng(int(round(seed * 10))).uniform(-w, w, n)


def strip(ax, x, vals, colour, median, floor=None, seed=0):
    vals = np.asarray(vals, float)
    ok = ~np.isnan(vals)
    ax.scatter(x + jitter(ok.sum(), seed=seed), vals[ok], s=4, color=fk.lighten(colour, 0.25), lw=0, zorder=2)
    if floor is not None and (~ok).sum():
        ax.scatter(x + jitter((~ok).sum(), seed=seed + 1), np.full((~ok).sum(), floor), s=6,
                   facecolors="white", edgecolors=colour, lw=0.5, zorder=2)
    ax.hlines(median, x - 0.32, x + 0.32, color="black", lw=1.3, zorder=3)


def panel_a(ax1, ax2, a):
    # left: rho(t), full (filled) and matched group sizes (hollow), by trait
    centres = []
    x = 0
    for t in ["dx", "age"]:
        for col, hollow in [("rho_t", False), ("rho_t_matched", True)]:
            v = a.loc[a.trait == t, col].values.astype(float)
            c = TRAIT_COL[t]
            if hollow:
                ax1.scatter(x + jitter(len(v), seed=x), v, s=5, facecolors="white", edgecolors=fk.lighten(c, 0.15),
                            lw=0.4, zorder=2)
                ax1.hlines(np.median(v), x - 0.32, x + 0.32, color="black", lw=1.3, zorder=3)
            else:
                strip(ax1, x, v, c, np.median(v), seed=x)
            x += 1
        centres.append(x - 1.5)
        x += 0.6
    ax1.axhline(0, color=fk.SUBTLE, lw=0.5, ls=":")
    ax1.set_ylim(-0.48, 0.8)  # room above the points for the legend
    ax1.set_xticks(centres); ax1.set_xticklabels([TRAIT_LAB[t] for t in ["dx", "age"]])
    ax1.set_ylabel("Agreement between two datasets\n(Spearman ρ of moderated t)")
    ax1.set_title("Datasets agree on age a little\nmore than on diagnosis")
    ax1.text(0.02, 0.015, "one point = one pair of datasets (192)", transform=ax1.transAxes, ha="left", va="bottom",
             fontsize=fk.SMALL, color=fk.SUBTLE)
    h = [plt.Line2D([], [], marker="o", ls="", ms=3.5, color=fk.TEXT, label="full data"),
         plt.Line2D([], [], marker="o", ls="", ms=3.5, mfc="white", mec=fk.TEXT, label="matched group sizes")]
    ax1.legend(handles=h, loc="upper right", fontsize=fk.SMALL, handletextpad=0.2, borderaxespad=0.2)
    # right: top-200 overlap, log2 fold over chance, full data
    floor = -1.0
    low = float(np.nanmin(a["top200_log2fold"].values.astype(float)))  # "none" counts go below the lowest point
    for i, t in enumerate(["dx", "age"]):
        v = a.loc[a.trait == t, "top200_log2fold"].values.astype(float)
        strip(ax2, i, v, TRAIT_COL[t], np.nanmedian(v), floor=floor, seed=10 + i)
        ax2.text(i, min(low, floor) - 0.25, f"{int(np.isnan(v).sum())}\nnone", ha="center", va="top", fontsize=fk.SMALL,
                 color=fk.SUBTLE)
    ax2.axhline(0, color=fk.SUBTLE, lw=0.5, ls=":")
    ax2.set_xticks(range(2)); ax2.set_xticklabels(["diag-\nnosis", "age"])
    ax2.set_ylabel("Shared top-200 genes\n(log2 fold over chance)")
    ax2.set_title("Same, by top-gene\noverlap")
    ax2.set_ylim(min(low, floor) - 1.2, None)
    ax2.set_xlim(-0.6, 1.6)


def panel_b(ax, p, s):
    sz = p[p["set"] == "SZ-study"]
    ext = p[p["set"] == "external"]
    med = s[s.relation_group == "subject-separate"].set_index("class").rho_t_median
    groups = [
        ("SZ\nvs SZ", sz.sz_sz_rho_t_median.values, fk.ACCENT, med["sz_sz"]),
        ("age\nvs age", sz.age_age_rho_t_median.values, fk.ACCENT_2, med["age_age"]),
        ("age vs\nexternal",
         np.concatenate([sz.age_ext_rho_t_median.values, ext.age_ext_rho_t_median.values]), fk.ACCENT_2, med["age_ext"]),
        ("external\nvs external", ext.ext_ext_rho_t_median.values, fk.ACCENT_2, med["ext_ext"]),
    ]
    for i, (lab, v, c, m) in enumerate(groups):
        v = v[~np.isnan(v)]
        strip(ax, i, v, c, m, seed=20 + i)
        ax.text(i, -0.22, f"n = {len(v)}", ha="center", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)
        ax.text(i + 0.36, m, f"{m:.2f}", ha="left", va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axhline(0, color=fk.SUBTLE, lw=0.5, ls=":")
    ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups], fontsize=fk.SMALL)
    ax.set_ylabel("Agreement of a dataset with its partners\n(median Spearman ρ of t over pairs)")
    ax.set_title("Age in schizophrenia-study data\nmatches age in external datasets")
    ax.text(0.02, 0.97, "one point = one dataset", transform=ax.transAxes, ha="left", va="top",
            fontsize=fk.SMALL, color=fk.SUBTLE)
    ax.set_ylim(-0.3, 0.66)
    ax.set_xlim(-0.6, 3.75)


def main():
    a = pd.read_csv(F_A, sep="\t", comment="#")
    p = pd.read_csv(F_B, sep="\t", comment="#")
    s = pd.read_csv(F_BS, sep="\t", comment="#")
    fig = fk.figure(74)
    gs = fig.add_gridspec(1, 3, width_ratios=[2.1, 1.05, 1.95], left=0.085, right=0.985, bottom=0.16, top=0.86, wspace=0.45)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    panel_a(ax1, ax2, a)
    panel_b(ax3, p, s)
    fk.letter(ax1, "a"); fk.letter(ax3, "b")
    fk.save(fig, "supp_positive_controls", "figures/scripts/build_supp_positive_controls.py", [F_A, F_B, F_BS],
            notes="Supplementary Fig. 3. a: 192 subject-separate bulk pairs; medians over pairs. Zero-overlap pairs "
                  "drawn hollow at the floor. b: points are per-dataset medians (per_dataset.tsv); black bars are "
                  "medians over subject-separate pairs (summary.tsv): sz_sz, age_age, age_ext, ext_ext.")


if __name__ == "__main__":
    main()
