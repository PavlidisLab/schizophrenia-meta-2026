#!/usr/bin/env python3
"""Supplementary Note 7 Figures 1 and 2 (death signatures): the Yaqubi contrasts and the Liharska resected-against-postmortem signature.

Data: data/ (yaqubi_*.tsv, liharska_*.tsv); the pH axis v2, the agonal axis and the 383 recurrent genes were applied when those tables were made.
Figure 1 (fig_yaqubi_axes), a-c: effect of surgical holding and autopsy (Yaqubi) and of the Dachet ex vivo time course against the pH and agonal axes,
  over all genes (a) and the 383 recurrent genes (b), and the share of the 383 whose schizophrenia sign matches the effect (c).
Figure 3 (fig_celltype_structure), a-d: cell-type structure of each contrast (R2 added by cell-type specificity over all genes and within the 383, the
  composition-only fit, and the implied change in the neuronal share).
Figure 2 (fig_livpm_vs_axes_383), A-C: Liharska and Collado-Torres living-minus-postmortem signatures against the axes (A), their direction against
  the recurrent schizophrenia pattern (B), and cell-class markers (C).
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import figkit as fk

ND = fk.note_data(__file__)
D = ND
FILES_Y = [D / f for f in ("yaqubi_axis_correlation_all_genes.tsv", "yaqubi_axis_correlation_recurrent_sets.tsv",
                           "yaqubi_direction_concordance.tsv", "yaqubi_dachet_axis_correlation.tsv")]
FILES_L = [D / f for f in ("liharska_genomewide_axis_correlations.tsv", "liharska_sets_direction.tsv", "liharska_celltype_view.tsv")]
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#", dtype={"set": str})
BASE = {"font.size": 9, "axes.titlesize": "large", "axes.labelsize": "medium", "xtick.labelsize": "medium", "ytick.labelsize": "medium", "legend.fontsize": "medium", "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#6b7280", "axes.labelcolor": "#1f2937",
        "xtick.color": "#6b7280", "ytick.color": "#6b7280", "axes.titlelocation": "left", "axes.titleweight": "normal"}
STYLE_L = {**BASE, "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#e5e7eb", "grid.linewidth": 0.6, "axes.axisbelow": True}


def figure_yaqubi():
    AL, AS, CO, DA = (rd("yaqubi_axis_correlation_all_genes.tsv"), rd("yaqubi_axis_correlation_recurrent_sets.tsv"),
                      rd("yaqubi_direction_concordance.tsv"), rd("yaqubi_dachet_axis_correlation.tsv"))
    ROWS = [("6 h 4 C", "S_6h4C_vs_0h"), ("24 h 4 C", "S_24h4C_vs_0h"), ("36 h 4 C", "S_36h4C_vs_0h"), ("6 h 20 C", "S_6h20C_vs_0h"), ("24 h 20 C", "S_24h20C_vs_0h"),
            ("36 h 20 C", "S_36h20C_vs_0h"), ("Autopsy, short PMI", "A_SPMI_vs_CTRL"), ("Autopsy, long PMI", "A_LPMI_vs_CTRL")]
    ROWS_D = [(f"Dachet {h} h", h) for h in (1, 2, 4, 8, 12, 24)]
    AX = {"pH axis v2": "#2563eb", "agonal axis (GTEx slow-minus-fast)": "#f59e0b"}
    SET = "383 recurrent (q_bh<0.1)"
    fig, axs = plt.subplots(1, 3, figsize=(13, 5.4), sharey=True)
    labels = [r[0] for r in ROWS] + [r[0] for r in ROWS_D]
    y = np.arange(len(labels))[::-1]
    for k, (an, col) in enumerate(AX.items()):
        off = (k - 0.5) * 0.3
        for i, (lab, c) in enumerate(ROWS):
            r = AL[(AL.effect == "t:" + c) & (AL.axis == an)].iloc[0]
            axs[0].errorbar(r.spearman, y[i] + off, xerr=[[r.spearman - r.spearman_boot95_low], [r.spearman_boot95_high - r.spearman]], fmt="o", color=col, ms=4, lw=1, label=an if i == 0 else None)
            r = AS[(AS.gene_set == SET) & (AS.effect == "t:" + c) & (AS.axis == an) & (AS.statistic == "spearman")].iloc[0]
            axs[1].plot([r.null_q025, r.null_q975], [y[i] + off] * 2, color=col, alpha=0.3, lw=6, solid_capstyle="butt")
            axs[1].errorbar(r.value, y[i] + off, xerr=[[r.value - r.boot95_low], [r.boot95_high - r.value]], fmt="o", color=col, ms=4, lw=1)
        for j, (lab, h) in enumerate(ROWS_D):
            i = len(ROWS) + j
            r = DA[(DA.scope == "all genes") & (DA.hours == h) & (DA.axis == an)].iloc[0]
            axs[0].errorbar(r.spearman, y[i] + off, xerr=[[r.spearman - r.boot95_low], [r.boot95_high - r.spearman]], fmt="s", color=col, ms=4, lw=1, mfc="white")
            r = DA[(DA.scope == SET) & (DA.hours == h) & (DA.axis == an)].iloc[0]
            axs[1].plot([r.null_q025, r.null_q975], [y[i] + off] * 2, color=col, alpha=0.3, lw=6, solid_capstyle="butt")
            axs[1].errorbar(r.spearman, y[i] + off, xerr=[[r.spearman - r.boot95_low], [r.boot95_high - r.spearman]], fmt="s", color=col, ms=4, lw=1, mfc="white")
    for i, (lab, c) in enumerate(ROWS):
        r = CO[(CO.gene_set == SET) & (CO.effect == "t:" + c)].iloc[0]
        axs[2].plot([r.null_q025, r.null_q975], [y[i]] * 2, color="#6b7280", alpha=0.3, lw=6, solid_capstyle="butt")
        axs[2].errorbar(r.frac_concordant, y[i], xerr=[[r.frac_concordant - r.wilson95_low], [r.wilson95_high - r.frac_concordant]], fmt="o", color="#1f2937", ms=4, lw=1)
    for ref, nm in (("ref: low-pH direction (minus pH axis z)", "#2563eb"), ("ref: slow-death direction (GTEx t)", "#f59e0b")):
        r = CO[(CO.gene_set == SET) & (CO.effect == ref)].iloc[0]
        axs[2].axvline(r.frac_concordant, color=nm, lw=1, ls="--")
    axs[0].set_title("a  All genes: effect vs axis"); axs[1].set_title("b  383 recurrent genes: effect vs axis"); axs[2].set_title("c  383 recurrent genes: SZ sign = effect sign")
    for a in axs[:2]:
        a.axvline(0, color="#6b7280", lw=0.6); a.set_xlabel("Spearman rho (bootstrap 95% over genes)")
    axs[2].axvline(0.5, color="#6b7280", lw=0.6); axs[2].set_xlabel("share concordant (Wilson 95%)")
    axs[0].set_yticks(y); axs[0].set_yticklabels(labels)
    for a in axs:
        a.xaxis.grid(True, color="#e5e7eb", lw=0.6); a.set_axisbelow(True)
    axs[2].legend(*axs[0].get_legend_handles_labels(), frameon=False, fontsize=7, loc="lower right")
    fig.text(0.01, 0.005, "Source: Yaqubi 2026 Source Data re-fit; Dachet GSE153669. Open squares = Dachet; grey/colour bars in b, c = matched-null 95% range; dashed lines in c = same test for the pH (blue) and agonal (amber) axes.", fontsize=7, color="#6b7280")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig


def figure_livpm():
    BLUE, GREEN, AMBER, RED, GREY, TEXT = "#2563eb", "#10b981", "#f59e0b", "#ef4444", "#6b7280", "#1f2937"
    SIG = [("LIV-PM Liharska (ST1)", "LIV-PM\nLiharska"), ("LIV-PM Collado-Torres, no qSVA", "LIV-PM\nno qSVA"), ("LIV-PM Collado-Torres, with qSVA", "LIV-PM\nwith qSVA"),
           ("SCZ1", "SCZ1"), ("SCZ2", "SCZ2"), ("ALZ1", "ALZ1"), ("ALZ2", "ALZ2"), ("PD1", "PD1"), ("PD2", "PD2")]
    gw, dr, ct = rd("liharska_genomewide_axis_correlations.tsv"), rd("liharska_sets_direction.tsv"), rd("liharska_celltype_view.tsv")
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
    fig.text(0.01, 0.005, "Sources: Liharska et al. 2025 Supplementary Tables (MOESM8, MOESM12); Collado-Torres et al. 2023 Table S4; pH axis v2 (2026-10-14); GTEx Hardy t; base-80 recurrence 383. Bars: 95% bootstrap over genes.", fontsize=7.5, color="#6b7280")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig


def figure_celltype():
    ROWS = [("Freshly resected vs postmortem, original", "L_ST1"), ("Freshly resected vs postmortem, no qSVA", "L_CT_noqSVA"), ("Freshly resected vs postmortem, with qSVA", "L_CT_qSVA"),
            ("Yaqubi holds vs 0 h, 36 h, 4 \u00b0C", "S_36h4C_vs_0h"), ("Yaqubi holds vs 0 h, 36 h, 20 \u00b0C", "S_36h20C_vs_0h"),
            ("Yaqubi autopsy, short interval", "A_SPMI_vs_CTRL"), ("Yaqubi autopsy, long interval", "A_LPMI_vs_CTRL"),
            ("Dachet 1 h", "D_1h"), ("Dachet 2 h", "D_2h"), ("Dachet 4 h", "D_4h"), ("Dachet 8 h", "D_8h"), ("Dachet 24 h", "D_24h"),
            ("pH axis", "REF_pH_axis_v2"), ("Control pH z", "REF_pH_Z"), ("GTEx slow minus fast death t", "REF_hardy_t")]
    H = rd("headline_by_contrast.tsv").set_index("contrast")
    BLUE, AMBER, GREY, TEXT = "#2563eb", "#f59e0b", "#6b7280", "#1f2937"
    fig, axs = plt.subplots(1, 4, figsize=(15, 5.6), sharey=True)
    y = np.arange(len(ROWS))[::-1]
    for i, (lab, c) in enumerate(ROWS):
        r = H.loc[c]
        col = GREY if c.startswith("REF") else BLUE
        if not c.startswith("REF"):
            axs[0].errorbar(r.R2_added_cell_type, y[i], xerr=[[r.R2_added_cell_type - r.ci95_lo], [r.ci95_hi - r.R2_added_cell_type]], fmt="o", color=col, ms=4, lw=1.2)
        else:
            axs[0].plot(r.R2_added_cell_type, y[i], "o", color=col, ms=4)
        axs[1].plot([r.R2_added_383_matched_null_mean, r.R2_added_383_matched_null_q975], [y[i]] * 2, color=GREY, alpha=0.3, lw=7, solid_capstyle="butt")
        axs[1].plot(r.R2_added_383, y[i], "o", color=col, ms=4)
        axs[2].plot(r.comp_shape_R2, y[i], "o", color=col, ms=4)
        if not c.startswith("REF"):
            axs[3].plot(r.implied_neuron_rel_change_pct, y[i], "o", color=col, ms=4)
    for a, t, x in zip(axs, ("a  All genes", "b  383 recurrent genes", "c  Composition-only fit", "d  Implied neuronal change"),
                       ("R\u00b2 added by cell-type specificity\n(95% interval over genes)", "R\u00b2 added by cell-type specificity\n(grey: matched genes, mean to 95th percentile)", "R\u00b2 of the composition-only fit", "Change in the neuronal share (%)")):
        a.set_title(t); a.set_xlabel(x); a.xaxis.grid(True, color="#e5e7eb", lw=0.6); a.set_axisbelow(True); a.grid(axis="y", visible=False)
    axs[3].axvline(0, color=GREY, lw=0.6)
    axs[0].set_yticks(y); axs[0].set_yticklabels([r[0] for r in ROWS])
    for a in axs: a.axhline(y[len(ROWS) - 3] + 0.5, color="#e5e7eb", lw=0.8)
    fig.text(0.01, 0.005, "Source: death-signature contrasts; PsychAD (HBCC) control cell-class specificity. Grey points, reference axes.", fontsize=7.5, color="#6b7280")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig


def main():
    fk.apply_style()
    with mpl.rc_context(BASE):
        fig = figure_yaqubi()
        fk.save(fig, "fig_yaqubi_axes", "supplementary_notes/death_signatures/build_supp_death_signatures.py", FILES_Y)
    with mpl.rc_context(STYLE_L):
        fig = figure_livpm()
        fk.save(fig, "fig_livpm_vs_axes_383", "supplementary_notes/death_signatures/build_supp_death_signatures.py", FILES_L)
    with mpl.rc_context(BASE):
        fig = figure_celltype()
        fk.save(fig, "fig_celltype_structure", "supplementary_notes/death_signatures/build_supp_death_signatures.py", [D / "headline_by_contrast.tsv"])


if __name__ == "__main__":
    main()
