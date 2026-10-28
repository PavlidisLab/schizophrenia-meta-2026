#!/usr/bin/env python3
"""Supplementary Note 8 Figure 1: pH association in single-nucleus data against the bulk axes at four pseudobulk levels, the composition-only simulation, and cell-class composition against pH and schizophrenia.

Data: data/ (summary rows only; no donor-level values): levels_all_rows.tsv, sim_rho_summary.tsv, comp_ph_effect_controls.tsv, comp_dx_effect.tsv.
"""
import hashlib, os, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import figkit as fk

D = fk.note_data(__file__)
ACC = {"all": "#2563eb", "rec383": "#10b981", "rec107": "#f59e0b"}
LAB = {"all": "all shared genes", "rec383": "383 recurrent genes", "rec107": "107 genes in 10+ lists"}
LEVELS = ["total", "total_composition", "class", "subclass"]
XLAB = {"total": "total\n(all nuclei)", "total_composition": "total, composition\nadjusted", "class": "class", "subclass": "subclass"}
GREY, TEXT, SUBTLE, GRID = "#9ca3af", "#1f2937", "#6b7280", "#e5e7eb"
LEVELS_RC = ({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"],
    "font.size": 9, "axes.edgecolor": SUBTLE, "axes.labelcolor": TEXT, "xtick.color": SUBTLE, "ytick.color": SUBTLE,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.titleweight": "normal", "axes.titlesize": 9, "figure.facecolor": "white", "axes.facecolor": "white"})



def key_figure():
    """Supplementary Note 8 Figure 1: the main points in one figure (NSW OFC pH association against the two bulk axes at four levels, the composition-only
    simulation, and the composition effects of pH and schizophrenia)."""
    def rd(f): return pd.read_csv(D / f, sep="\t", comment="#")
    allr = rd("levels_all_rows.tsv"); R = rd("sim_rho_summary.tsv"); ph = rd("comp_ph_effect_controls.tsv"); dx = rd("comp_dx_effect.tsv")
    BL, GR, AM, TX, SU, GD = "#2563eb", "#10b981", "#f59e0b", "#1f2937", "#6b7280", "#e5e7eb"
    ph = ph[ph.model == "crumblr ~ scale(age) + sex01 + scale(ph)"]; dx = dx[dx.term == "dxSchizophrenia"]
    rows = []
    with mpl.rc_context({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"], "svg.fonttype": "none", "font.size": 8.5,
                         "axes.edgecolor": SU, "axes.labelcolor": TX, "xtick.color": SU, "ytick.color": SU, "axes.spines.top": False, "axes.spines.right": False}):
        fig = plt.figure(figsize=(13, 4.6)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.42, left=0.06, right=0.985, top=0.88, bottom=0.27)
        for k, (axis, ttl) in enumerate((("pH", "a  Correlation with the bulk pH axis"), ("agonal", "b  Correlation with the bulk agonal axis"))):
            ax = fig.add_subplot(gs[0, k]); ax.set_title(ttl, loc="left", color=TX, fontsize=9)
            d = allr[(allr.dataset == "NSW_OFC") & (allr.axis == axis) & allr.population.str.startswith("controls")]
            for gset, off, col, lab in (("all", -0.1, BL, "all shared genes"), ("rec383", 0.1, GR, "383 recurrent genes")):
                xs, ys, lo, hi = [], [], [], []
                for i, lev in enumerate(LEVELS):
                    q = d[(d.gene_set == gset) & (d.level == lev)]
                    q = q[q.cell_type == "all"] if lev in ("total", "total_composition") else q[q.cell_type == "MEDIAN_over_cell_types_weighted_nuclei"]
                    r = q.iloc[0]; xs.append(i + off); ys.append(r.rho); lo.append(r.lo); hi.append(r.hi)
                    rows.append(dict(panel=("a", "b")[k], series=lab, level=lev, est=r.rho, lo=r.lo, hi=r.hi))
                ys = np.array(ys, float)
                ax.errorbar(xs, ys, yerr=[ys - np.array(lo, float), np.array(hi, float) - ys], fmt="none", ecolor=col, elinewidth=1.2, capsize=2.5)
                ax.plot(xs, ys, "-o", color=col, ms=5, lw=1.1, label=lab)
            sim = [R[(R.tag == "main_grid") & (R.N == 200) & (R.rna == 1) & (R.scale == 1.0) & (R.profile == pr) & (R.axis == axis) & (R.gene_set == "all")].iloc[0] for pr in ("HBCC", "MSSM")]
            lo_s, hi_s = min(r.rho_mean for r in sim), max(r.rho_mean for r in sim)
            ax.axhspan(lo_s, hi_s, color=AM, alpha=0.25, lw=0, label="composition-only simulation (all genes)")
            for pr, r in zip(("HBCC", "MSSM"), sim): rows.append(dict(panel=("a", "b")[k], series="simulation, all genes, %s profiles" % pr, level="simulated", est=r.rho_mean, lo=r.rho_lo, hi=r.rho_hi))
            ax.axhline(0, color=SU, lw=0.8); ax.set_xticks(range(4)); ax.set_xticklabels([XLAB[l] for l in LEVELS], fontsize=7.5); ax.set_xlim(-0.5, 3.5)
            ax.set_ylabel("Spearman correlation across genes"); ax.grid(axis="y", color=GD, lw=0.6); ax.set_axisbelow(True)
            if k == 0: ax.legend(frameon=False, fontsize=7.5, loc="lower right")
        ax = fig.add_subplot(gs[0, 2]); ax.set_title("c  Cell-class composition in NSW OFC", loc="left", color=TX, fontsize=9)
        p = ph[(ph.dataset == "NSW_OFC") & (ph.level == "class")].set_index("category"); cats = list(p.index)
        a0 = dx[(dx.dataset == "NSW_OFC") & (dx.level == "class") & (dx.model_id == "M0") & (dx.donor_set == "all")].set_index("category")
        a1 = dx[(dx.dataset == "NSW_OFC") & (dx.level == "class") & (dx.model_id == "M1") & (dx.donor_set == "all")].set_index("category")
        ys_ = np.arange(len(cats))[::-1]
        for off, (lab, tab, col, per_sd, mk) in zip([0.25, 0, -0.25], [("pH effect in controls (per SD of pH)", p, AM, True, "o"), ("schizophrenia effect", a0, BL, False, "o"), ("schizophrenia effect, pH in the model", a1, GR, False, "s")]):
            for y, c in zip(ys_, cats):
                r = tab.loc[c]; e, lo, hi = (r.est_per_sd, r.lo_per_sd, r.hi_per_sd) if per_sd else (r.est, r.lo, r.hi)
                ax.plot([lo, hi], [y + off] * 2, color=col, lw=1.2); ax.plot(e, y + off, mk, color=col, ms=4, label=lab if c == cats[0] else None)
                rows.append(dict(panel="c", series=lab, level=c, est=e, lo=lo, hi=hi))
        ax.axvline(0, color=GD, lw=0.8, zorder=0); ax.set_yticks(ys_); ax.set_yticklabels([c.replace("_", " ") for c in cats]); ax.grid(axis="x", color=GD, lw=0.5)
        ax.set_xlabel("difference in centred log-ratio composition"); ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.35, -0.2), ncol=1)
        for a in fig.axes:
            for art in list(a.patches) + list(a.lines) + list(a.texts) + list(a.collections): art.set_clip_on(False)
        fk.save(fig, "fig_composition_key", "supplementary_notes/composition_and_axes/build_supp_composition.py",
                [D / f for f in ("levels_all_rows.tsv", "sim_rho_summary.tsv", "comp_ph_effect_controls.tsv", "comp_dx_effect.tsv")])
    pd.DataFrame(rows).to_csv(D / "fig_composition_key_data.tsv", sep="\t", index=False)


if __name__ == "__main__":
    key_figure()
