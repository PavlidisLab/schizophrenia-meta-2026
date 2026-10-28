#!/usr/bin/env python3
"""Regenerates fig_livpm_vs_axes_383.svg/.png and fig_livpm_data.tsv from the TSVs written by s1_main.py (no numbers typed in). Run: python3 -I s2_figure.py"""
import os, sys
sys.dont_write_bytecode = True
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rd = lambda f: pd.read_csv(os.path.join(OUT, f), sep="\t", comment="#", dtype={"set": str})
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"], "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#6b7280", "axes.labelcolor": "#1f2937", "xtick.color": "#6b7280", "ytick.color": "#6b7280",
                     "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#e5e7eb", "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.titleweight": "normal", "axes.titlelocation": "left"})
BLUE, GREEN, AMBER, RED, GREY, TEXT = "#2563eb", "#10b981", "#f59e0b", "#ef4444", "#6b7280", "#1f2937"
SIG = [("LIV-PM Liharska (ST1)", "LIV-PM\nLiharska"), ("LIV-PM Collado-Torres, no qSVA", "LIV-PM\nno qSVA"), ("LIV-PM Collado-Torres, with qSVA", "LIV-PM\nwith qSVA"),
       ("SCZ1", "SCZ1"), ("SCZ2", "SCZ2"), ("ALZ1", "ALZ1"), ("ALZ2", "ALZ2"), ("PD1", "PD1"), ("PD2", "PD2")]
gw, dr, ct = rd("genomewide_axis_correlations.tsv"), rd("sets_direction.tsv"), rd("celltype_view.tsv")
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={"width_ratios": [1.15, 1.15, 1.2]})
data = []
# A
a = ax[0]
for k, (ax_key, col, off) in enumerate((("pH", BLUE, -0.12), ("agonal", AMBER, 0.12))):
    for i, (full, short) in enumerate(SIG):
        r = gw[(gw.signature == full) & gw.axis.str.startswith(ax_key)].iloc[0]
        a.errorbar(i + off, r.spearman, yerr=[[r.spearman - r.spearman_boot95_low], [r.spearman_boot95_high - r.spearman]], fmt="o", color=col, ms=5, lw=1.2, label=("pH axis v2" if ax_key == "pH" and i == 0 else ("agonal axis" if ax_key == "agonal" and i == 0 else None)))
        data.append(dict(panel="A", signature=full, series=ax_key, value=r.spearman, lo=r.spearman_boot95_low, hi=r.spearman_boot95_high, n=r.n_genes))
a.axhline(0, color=GREY, lw=0.8); a.set_xticks(range(len(SIG))); a.set_xticklabels([s for _, s in SIG], fontsize=8)
a.set_ylabel("Spearman across genes with the axis"); a.set_title("A  Signatures against the pH and agonal axes"); a.legend(frameon=False, fontsize=8, loc="lower right")
# B
b = ax[1]
SB = [("LIV-PM Liharska (ST1)", "LIV-PM\nLiharska"), ("LIV-PM Collado-Torres, no qSVA", "LIV-PM\nno qSVA"), ("LIV-PM Collado-Torres, with qSVA", "LIV-PM\nwith qSVA"),
      ("LIV-PM Liharska (ST1) | pH and agonal axes removed", "LIV-PM\naxes out"), ("SCZ1", "SCZ1"), ("SCZ2", "SCZ2"), ("ALZ1", "ALZ1"), ("ALZ2", "ALZ2")]
for i, (full, short) in enumerate(SB):
    r = dr[(dr.signature == full) & (dr.set == "383")].iloc[0]
    b.add_patch(plt.Rectangle((i - 0.3, r.null_frac_q025), 0.6, r.null_frac_q975 - r.null_frac_q025, color="#e5e7eb", zorder=1))
    b.plot([i - 0.3, i + 0.3], [r.null_frac_mean] * 2, color=GREY, lw=1, zorder=2)
    b.errorbar(i, r.frac_concordant, yerr=[[r.frac_concordant - r.frac_boot95_low], [r.frac_boot95_high - r.frac_concordant]], fmt="o", color=BLUE, ms=5, lw=1.2, zorder=3)
    data.append(dict(panel="B", signature=full, series="383 concordance", value=r.frac_concordant, lo=r.frac_boot95_low, hi=r.frac_boot95_high, n=r.n_genes_with_SZ_sign_measured, null_mean=r.null_frac_mean, null_lo=r.null_frac_q025, null_hi=r.null_frac_q975))
b.axhline(0.5, color=GREY, lw=0.8, ls=":"); b.set_xticks(range(len(SB))); b.set_xticklabels([s for _, s in SB], fontsize=8); b.set_xlim(-0.6, len(SB) - 0.4)
b.set_ylabel("Share of the 383 genes with the SZ direction"); b.set_title("B  Same direction as the recurrent SZ pattern")
b.text(-0.5, 0.915, "grey box: 2.5-97.5% of 10,000 matched gene sets (bar = mean)", ha="left", fontsize=7.5, color=GREY); b.set_ylim(0.43, 0.93)
# C
c = ax[2]
CL = ["Astro", "EN", "IN_other", "LAMP5", "Micro", "OPC", "Oligo", "PVALB", "SST", "VIP", "Vascular"]
for sg, col, off in (("LIV-PM Liharska (ST1)", BLUE, -0.13), ("LIV-PM Collado-Torres, with qSVA", AMBER, 0.13)):
    for i, cl in enumerate(CL):
        r = ct[(ct.signature == sg) & (ct.cell_class == cl)].iloc[0]
        lo, hi = r.mean_boot95_low - r.matched_null_mean, r.mean_boot95_high - r.matched_null_mean
        c.errorbar(i + off, r.effect_over_null, yerr=[[r.effect_over_null - lo], [hi - r.effect_over_null]], fmt="o", color=col, ms=5, lw=1.2,
                   label=(("LIV-PM Liharska" if col == BLUE else "LIV-PM with qSVA") if i == 0 else None))
        data.append(dict(panel="C", signature=sg, series=cl, value=r.effect_over_null, lo=lo, hi=hi, n=r.n_markers, null_mean=r.matched_null_mean))
c.axhline(0, color=GREY, lw=0.8); c.set_xticks(range(len(CL))); c.set_xticklabels(CL, fontsize=8, rotation=40, ha="right")
c.set_ylabel("Marker mean t minus matched-gene mean\n(positive = higher in postmortem)"); c.set_title("C  Cell-class markers (top 100 by specificity)"); c.legend(frameon=False, fontsize=8, loc="lower left")
fig.text(0.01, 0.005, "Sources: Liharska et al. 2025 Supplementary Tables (MOESM8, MOESM12); Collado-Torres et al. 2023 Table S4; pH axis v2 (2026-10-14); GTEx Hardy t; base-80 recurrence 383. Bars: 95% bootstrap over genes. Script: scripts/s2_figure.py", fontsize=7.5, color="#6b7280")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(OUT, "fig_livpm_vs_axes_383.svg")); fig.savefig(os.path.join(OUT, "fig_livpm_vs_axes_383.png"), dpi=200)
pd.DataFrame(data).to_csv(os.path.join(OUT, "fig_livpm_data.tsv"), sep="\t", index=False)
