#!/usr/bin/env python3
"""Supplementary Note 5 Figure 2 -- SNAP in bulk tissue.

Data: data/snap_what_is_it/ (2026-10-12).
a  coexpression of SNAP's top-loading genes (top 1,000 by |weight|, aligned by loading sign) in each bulk dataset:
   mean pairwise correlation minus that of expression-decile-matched genes (snap_bulk_coherence.tsv, variant raw)
b  schizophrenia effect on the loading-weighted SNAP score in each bulk dataset (SD, given age and sex), and with the
   estimated agonal score in the model (SNAP's top 1,000 genes removed from the axis) (snap_estimated_agonal.tsv)
"""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import figkit as fk

ND = fk.note_data(__file__)  # this note's own data; shared tables are read from fk.FD

D = ND / "snap_what_is_it"
FILES = [D / "snap_bulk_coherence.tsv", D / "snap_estimated_agonal.tsv"]
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
BASE, ADJ = fk.ACCENT, fk.ACCENT_2


def label(ds):
    return (ds.replace("GTEx_Anterior_cingulate_cortex_BA24", "GTEx BA24").replace("GTEx_Frontal_Cortex_BA9", "GTEx BA9")
            .replace("GTEx_Cortex", "GTEx cortex").replace("_", " "))


def panel_a(ax):
    c = rd("snap_bulk_coherence.tsv")
    c = c[(c["set"] == "all") & (c["variant"] == "raw") & (c["top_n"] == 1000) & (c["aligned"].astype(str) == "True")]
    per = c[~c["dataset"].str.contains("pooled|median")].sort_values("effect", ascending=False).reset_index(drop=True)
    pooled = c[c["dataset"].str.startswith("pooled")].iloc[0]
    y = np.arange(len(per))
    ax.errorbar(per["effect"], y, xerr=[per["effect"] - per["lo95"], per["hi95"] - per["effect"]], fmt="o", ms=3,
                color=BASE, elinewidth=0.8, capsize=0)
    yp = len(per) + 0.6
    ax.errorbar([pooled["effect"]], [yp], xerr=[[pooled["effect"] - pooled["lo95"]], [pooled["hi95"] - pooled["effect"]]],
                fmt="D", ms=4, color=fk.TEXT, elinewidth=1.0, capsize=0)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels([label(d) for d in per["dataset"]] + ["pooled, 22 brain-bank datasets"], fontsize=fk.TINY)
    ax.set_ylim(yp + 0.8, -0.8)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Mean correlation among SNAP genes, minus\nexpression-matched genes (95% interval)")
    ax.set_title("Coexpression of SNAP genes in bulk tissue", loc="left")


def panel_b(ax):
    e = rd("snap_estimated_agonal.tsv")
    e = e[(e["analysis"] == "sz_effect") & (e["cohort_group"] == "bulk")]
    base = e[e["model"] == "base: age, sex"].set_index("dataset")
    adj = e[e["model"] == "+ estimated agonal (top 1,000 removed)"].set_index("dataset")
    per = [d for d in base.index if not d.startswith("pooled")]
    per = sorted(per, key=lambda d: base.loc[d, "effect"])
    pooled = [d for d in base.index if d.startswith("pooled over 8")]
    rows = per + pooled
    y = np.arange(len(per)).tolist() + [len(per) + 0.6]
    for yi, d in zip(y, rows):
        for off, tab, col, m in [(-0.16, base, BASE, "o"), (0.16, adj, ADJ, "s")]:
            if d not in tab.index:
                continue
            r = tab.loc[d]
            big = d.startswith("pooled")
            ax.errorbar([r["effect"]], [yi + off], xerr=[[r["effect"] - r["lo95"]], [r["hi95"] - r["effect"]]],
                        fmt="D" if big else m, ms=4 if big else 3, color=col, elinewidth=0.8, capsize=0)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([label(d) for d in per] + ["pooled, 22 brain-bank datasets"], fontsize=fk.TINY)
    ax.set_ylim(y[-1] + 0.8, -0.8)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Schizophrenia effect on SNAP score\n(SD, 95% interval)")
    ax.set_title("Schizophrenia effect on SNAP in bulk tissue", loc="left")
    ax.legend(handles=[Line2D([], [], marker="o", color=BASE, lw=0.8, markersize=3.5, label="age and sex"),
                       Line2D([], [], marker="s", color=ADJ, lw=0.8, markersize=3.5,
                              label="also estimated agonal score")],
              frameon=False, fontsize=fk.SMALL, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2,
              title="model", title_fontsize=fk.SMALL, handletextpad=0.3, columnspacing=1.0)


def main():
    fk.apply_style()
    fig = fk.figure(120)
    axa = fig.add_axes([0.20, 0.19, 0.30, 0.71])
    axb = fig.add_axes([0.68, 0.19, 0.30, 0.71])
    panel_a(axa)
    panel_b(axb)
    for ax, l in zip((axa, axb), "ab"):
        fk.letter(ax, l)
    fk.save(fig, "supp_snap_bulk", "supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_snap_bulk.py", FILES,
            notes="Coherence top 1,000 aligned, raw; SZ effect given age and sex, and with est. agonal (top 1,000 removed).")


if __name__ == "__main__":
    main()
