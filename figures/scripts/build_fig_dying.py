#!/usr/bin/env python3
"""Figure 5 -- the pH signature resembles the dying process; which genes are pH-sensitive.

Data: figures/data/fig08/ and figures/data/fig_interneuron/.

a  GTEx v8 slow vs fast natural death (Hardy 4 vs 2) against the control-pH association, per brain region
c  agreement of 192 subject-separate bulk pairs before and after removing the pH or GTEx agonal axis
d  GSE92538: cerebellar pH against agonal factor, by diagnosis (one point per subject)
e  per-study agonal score of the listed genes
f  other disorders' postmortem brain: correlation of each study's effects on the recurrent genes with the pH and agonal axes,
   minus matched random genes; one dot per study (figures/data/supp_fig_other_disorders_corrected/panel_c_per_study.tsv)
g  pH association of cell-class marker genes (top 100 specific; PsychAD controls) over expression-matched genes

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

import figkit as fk

D8 = fk.FD / "fig08"
DI = fk.FD / "fig_interneuron"
rd = lambda d, f: pd.read_csv(d / f, sep="\t", comment="#")
rng = np.random.default_rng(1)


def panel_a(ax):
    df = rd(D8, "fig08_a_dachet_gene_scatter.tsv").dropna(subset=["ph_assoc", "change_24h"])
    ax.hexbin(df["ph_assoc"], df["change_24h"], gridsize=45, cmap="Greys", mincnt=1, bins="log", linewidths=0)
    rho = df[["ph_assoc", "change_24h"]].corr(method="spearman").iloc[0, 1]
    ax.text(0.97, 0.96, f"ρ = {rho:.2f}\none point = one gene\n({len(df):,}, shaded by density)", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL)
    ax.axhline(0, color=fk.GRID, lw=0.6)
    ax.axvline(0, color=fk.GRID, lw=0.6)
    ax.set_xlabel("Gene's association with pH in controls\n(combined z)")
    ax.set_ylabel("Gene's change 24 h after resection\n(log2, room temperature)")
    # one title for both plots of panel a; it runs on over the right-hand plot
    ax.set_title("Resected cortex left at room temperature changes\nlike low pH, from the first hour", pad=4)
    return rho, len(df)


def panel_b(ax):
    df = rd(D8, "fig08_a_alignment_by_time.tsv")
    for axis, col, lab in [("pH Z", fk.ACCENT_3, "pH in controls"), ("PMI Z", fk.NEUTRAL, "PMI in controls"),
                           ("Gandal SCZ t", fk.ACCENT, "SZ (Gandal)")]:
        s = df[df["axis"] == axis].sort_values("hours")
        ax.plot(s["hours"], s["rho"], "-o", color=col, ms=2.5, lw=1, label=lab)
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8, 12, 24])
    ax.set_xticklabels(["1", "2", "4", "8", "12", "24"])
    ax.set_xlabel("Time after resection (hours at room temperature)")
    ax.set_ylabel("Correlation of the change with\neach reference (Spearman ρ, genes)")
    lo, hi = ax.get_ylim()  # empty band above the lines holds the legend, clear of the PMI line
    ax.set_ylim(lo, hi + 0.30)
    ax.legend(loc="upper right", ncol=1, fontsize=fk.SMALL, handlelength=1, borderaxespad=0.2, labelspacing=0.2,
              title="reference: each gene's\nassociation with", title_fontsize=fk.SMALL, alignment="left")


def panel_c(ax):
    df = rd(D8, "fig08_e_gtex_by_region.tsv").sort_values("rho_hardy_vs_ph")
    y = np.arange(len(df))
    ax.barh(y, df["rho_hardy_vs_ph"], 0.7, color=fk.ACCENT_3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["region"].str.replace("Brain - ", "").str.replace(r" \(.*\)", "", regex=True), fontsize=fk.SMALL)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Slow-death effect against pH\nassociation (Spearman ρ, genes)")
    ax.set_title("Slow death and low pH\nin GTEx brain regions", pad=4)
    lo, hi = ax.get_ylim()  # one empty row above the top bar holds the note
    ax.set_ylim((lo, hi - 1.0) if ax.yaxis_inverted() else (lo, hi + 1.0))
    ax.text(0.03, 0.985, "one bar = one GTEx brain region", transform=ax.transAxes, fontsize=fk.SMALL, color=fk.SUBTLE, va="top")
    ax.tick_params(axis="y", length=0)
    return float(df["rho_hardy_vs_ph"].median()), int((df["rho_hardy_vs_ph"] < 0).sum()), len(df)


def panel_d(ax):
    df = rd(D8, "fig08_b_pair_agreement_gtex_axis.tsv")
    df = df[df["pair_class"] == "both-bulk donor-separate"]
    order = [("none (raw)", "none"), ("pH axis removed", "pH axis"), ("GTEx agonal axis removed", "agonal axis"),
             ("pH and GTEx agonal removed", "both")]
    x = np.arange(len(order))
    # L0 is unadjusted; the dark bars were mislabelled "age, sex"
    for off, lev, col, lab in [(-0.18, "L0", fk.ACCENT, "no covariates"), (0.18, "L4", fk.lighten(fk.ACCENT, 0.5), "all covariates")]:
        s = df[df["level"] == lev].set_index("axis_removed")
        v = [s.loc[k, "median_rho"] for k, _ in order]
        ax.bar(x + off, v, 0.34, color=col, label=lab)
        for xi, vi in zip(x, v):
            ax.text(xi + off, vi + 0.002, f"{vi:.2f}", ha="center", va="bottom", fontsize=5.5, color=fk.TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(["none", "pH", "agonal", "both"], fontsize=fk.SMALL)
    ax.set_xlabel("Axis regressed out of both datasets")
    ax.set_ylabel("Agreement between SZ datasets\n(median Spearman ρ of t, 192 pairs)")
    ax.set_title("Agreement between datasets after\nremoving the pH or agonal axis", pad=4)
    ax.set_ylim(0, 0.14)
    ax.legend(loc="upper right", fontsize=fk.SMALL, handlelength=0.8)


DX = [("Control", "control", fk.NEUTRAL), ("Major Depressive Disorder", "MDD", fk.ACCENT_2),
      ("Bipolar Disorder", "BD", fk.ACCENT), ("Schizophrenia", "SZ", fk.ACCENT_4)]


def panel_e(ax):
    df = rd(D8, "fig08_f_gse92538_subjects.tsv")
    df = df[df["complete_case_for_models"] == "yes"]
    for k, (dx, lab, col) in enumerate(DX):
        s = df[df["diagnosis"] == dx]
        x = s["agonal_factor"] + (k - 1.5) * 0.12 + rng.uniform(-0.04, 0.04, len(s))
        ax.scatter(x, s["ph_cerebellum"], s=5, color=col, lw=0, label=f"{lab} ({len(s)})")
    ax.set_xticks(range(5))
    ax.set_xlabel("Agonal factor score\n(0 to 4; Pritzker)")
    ax.set_ylabel("Brain pH (cerebellum)")
    ax.set_title("Brain pH by agonal\nfactor, Pritzker", pad=4)
    ax.legend(loc="upper right", fontsize=fk.SMALL, handletextpad=0.1, markerscale=1.5,
              title="one point = one subject", title_fontsize=fk.SMALL)
    return len(df)


def panel_f(ax):
    df = rd(D8, "fig08_g_hitlist_agonal_score_by_study.tsv").dropna(subset=["agonal_score"])
    y = df["agonal_score"].to_numpy()
    ax.scatter(rng.uniform(-0.2, 0.2, len(y)), y, s=6, color=fk.ACCENT_3, lw=0)
    ax.plot([-0.3, 0.3], [np.median(y)] * 2, color=fk.TEXT, lw=1.2)
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([0])
    ax.set_xticklabels([f"one point = one study\n{len(y)} studies,\n{(y > 0).sum()} above 0"], fontsize=fk.SMALL)
    ax.set_ylabel("Agonal signature in the study's list\n(mean GTEx slow-death t, signed)")
    ax.set_title("Slow-death\nsignature of each\nstudy's list", pad=4)
    return float(np.median(y)), int((y > 0).sum()), len(y)


CLASSES = [("SST", "SST"), ("PVALB", "PVALB"), ("LAMP5", "LAMP5"), ("VIP", "VIP"), ("IN_other", "other interneurons"),
           ("EN", "excitatory"), ("Micro", "microglia"), ("Vascular", "vascular"), ("OPC", "OPC"),
           ("Oligo", "oligodendrocyte"), ("Astro", "astrocyte")]


def panel_g(ax):
    # intervals widened by sqrt(VIF) for correlation among each class's marker genes
    df = rd(DI, "fig5g_marker_sets_vif.tsv").rename(columns={"corrected_ci_lo": "effect_ci95_lo",
                                                              "corrected_ci_hi": "effect_ci95_hi"})
    y = np.arange(len(CLASSES))
    for off, cohort, col in [(-0.17, "HBCC", fk.ACCENT), (0.17, "MSSM", fk.lighten(fk.ACCENT, 0.5))]:
        s = df[df["cohort"] == cohort].set_index("class")
        v = [s.loc[c] for c, _ in CLASSES]
        ax.errorbar([r["effect"] for r in v], y + off,
                    xerr=[[r["effect"] - r["effect_ci95_lo"] for r in v], [r["effect_ci95_hi"] - r["effect"] for r in v]],
                    fmt="o", ms=2.5, color=col, elinewidth=0.8, capsize=0, label="PsychAD " + ("Mount Sinai" if cohort == "MSSM" else cohort))
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l in CLASSES], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("pH sensitivity of marker genes\n(control-pH z minus matched)")
    ax.set_title("Cell-type-specific differences\nin pH association", pad=4)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", fontsize=fk.TINY, handletextpad=0.2, handlelength=1.2, frameon=False, borderaxespad=0.1)
    ax.set_xlim(right=ax.get_xlim()[1] + 2.0)

# panel d shows schizophrenia as one dot per donor group of the 26 bulk
# datasets beside the other disease classes; data figures/data/fig4d_sz_dots/ (copied from RA's folder, sha256 in SHA256SUMS)
D4D = fk.FD / "fig4d_sz_dots"
SURVEY_ROWS = [("SZ", "Schizophrenia", False),
               ("SUD", "Substance use", False), ("MS", "Multiple sclerosis", False), ("BD", "Bipolar disorder", False),
               ("AD", "Alzheimer disease", False), ("ASD", "Autism", False),
               ("MDD", "Major depression", False), ("PD", "Parkinson disease", False),
               ("OTHER", "Other", False)]


def panel_survey(ax):
    d = pd.read_csv(D4D / "figure_data_fig4d_v4.tsv", sep="\t", comment="#")
    sm = pd.read_csv(D4D / "figure_data_fig4d_v4_summary.tsv", sep="\t", comment="#").set_index("class_shown")
    lanes = [("effect_pH", "pH", "pH axis", fk.ACCENT_3, "o", -0.17), ("effect_agonal", "agonal", "agonal axis", fk.ACCENT_2, "s", 0.17)]
    ys = np.arange(len(SURVEY_ROWS))
    for y, (cls, lab, hollow) in zip(ys, SURVEY_ROWS):
        m = d[d["class_shown"] == cls]
        for col, key, _, c, mk, off in lanes:
            x = m[col].to_numpy(float)
            jit = rng.uniform(-0.08, 0.08, len(x))
            ax.scatter(x, y + off + jit, s=5, marker=mk, color=c, alpha=0.6, lw=0, zorder=2)
            ax.plot([sm.loc[cls, f"boot_lo_{key}"], sm.loc[cls, f"boot_hi_{key}"]], [y + off] * 2, color=c, lw=0.9, zorder=2)
            ax.scatter([sm.loc[cls, f"median_{key}"]], [y + off], s=22, marker="|", color=fk.TEXT, linewidths=1.6, zorder=3)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.axhline(0.5, color=fk.GRID, lw=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{lab} ({int(sm.loc[cls, 'n'])})" for cls, lab, _ in SURVEY_ROWS], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Correlation of the disease effects with the axis,\nminus matched random genes (Spearman \u03c1)")
    ax.set_title("pH and agonal axes in other brain conditions", pad=12)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker=mk, color=c, lw=0, markersize=3.5, label=l) for _, _, l, c, mk, _ in lanes],
              loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=fk.TINY, handletextpad=0.2, borderaxespad=0.1)


DC = fk.REPO / "supplementary_notes" / "composition_and_axes" / "data"
LEVELS = [("Whole sample", "total"), ("Whole sample, cell-class\nfractions added", "total_composition"),
          ("Within cell class", "class"), ("Within cell subclass", "subclass")]
SIM_ROWS = [("Simulation (HBCC)", "HBCC"), ("Simulation (Mount Sinai)", "MSSM")]


def panel_levels(ax):
    d = pd.read_csv(DC / "levels_summary_main_rows.tsv", sep="\t", comment="#")
    d = d[(d["dataset"] == "NSW_OFC") & (d["population"] == "controls") & (d["axis"] == "pH") & (d["gene_set"] == "all")].set_index("level")
    sim = pd.read_csv(DC / "sim_fraction_of_ceiling.tsv", sep="\t", comment="#")
    sim = sim[sim["rna"] == 1].set_index("profile")
    rows = [(l, d.loc[k, "rho"], d.loc[k, "lo"], d.loc[k, "hi"], fk.ACCENT) for l, k in LEVELS]
    rows += [(l, sim.loc[k, "rho_mean"], sim.loc[k, "rho_lo"], sim.loc[k, "rho_hi"], fk.SUBTLE) for l, k in SIM_ROWS]
    y = np.arange(len(rows)) + (np.arange(len(rows)) >= len(LEVELS)) * 0.5
    for yy, (_, m, lo, hi, c) in zip(y, rows):
        ax.plot([lo, hi], [yy, yy], color=c, lw=0.9)
        ax.scatter([m], [yy], s=16, color=c, zorder=3)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Agreement with the bulk pH axis\n(Spearman \u03c1 over genes)")
    ax.set_title("pH association in single-\nnucleus data at each level\nof pooling (32 controls)", pad=4)
    return float(d.loc["total", "rho"])


# label offsets in data units, per panel (pH axis, agonal axis)
MARKERS_G = [{"SST": (5.6, -14), "PVALB": (2.4, -19), "GAD1": (6.3, 9), "AQP4": (-4.2, -15), "GFAP": (-5.2, 22)},
             {"SST": (-6.0, -14), "PVALB": (-4.2, -19), "GAD1": (-6.2, -3.5), "AQP4": (8.0, -4), "GFAP": (2.5, 22)}]


def panel_resected(ax):
    """Freshly resected against postmortem cortex (Liharska et al.) against the pH axis and the agonal axis, over the recurrent genes."""
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    from scipy.stats import spearmanr
    d = rd(D8, "fig08_h_liharska_per_gene_383.tsv")
    fig = ax.figure
    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=ax.get_subplotspec(), wspace=0.06)
    ax.remove()
    out = []
    for k, (col, xl, lab) in enumerate([("ph_axis_v2_z", "pH axis (z)", "pH axis"), ("agonal_axis_t", "Agonal axis (t)", "agonal axis")]):
        a = fig.add_subplot(sub[k])
        v = d[d["in_374_with_direction"].astype(str).str.upper() == "TRUE"].dropna(subset=["liv_pm_t", col])
        up = v["sz_direction"] == "up"
        dn = v["sz_direction"] == "down"
        a.scatter(v.loc[dn, col], v.loc[dn, "liv_pm_t"], s=6, color=fk.DOWN, lw=0, alpha=0.75)
        a.scatter(v.loc[up, col], v.loc[up, "liv_pm_t"], s=6, color=fk.UP, lw=0, alpha=0.75)
        r = spearmanr(v[col], v["liv_pm_t"])[0]
        # cell-type marker genes named, as in Fig. 3b: interneuron markers (lower in schizophrenia) and astrocyte markers (higher)
        pts = v.set_index("gene")
        for g, off in MARKERS_G[k].items():
            if g in pts.index:
                gx, gy = float(pts.loc[g, col]), float(pts.loc[g, "liv_pm_t"])
                a.scatter([gx], [gy], s=12, facecolor="none", edgecolor=fk.TEXT, lw=0.7, zorder=5)
                a.annotate(g, (gx, gy), xytext=off, textcoords="data", fontsize=fk.TINY, color=fk.TEXT,
                           ha="center", va="center", zorder=6, bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.9),
                           arrowprops=dict(arrowstyle="-", color=fk.TEXT, lw=0.4, shrinkA=6, shrinkB=2))
        a.axhline(0, color=fk.TEXT, lw=0.6)
        a.axvline(0, color=fk.TEXT, lw=0.6)
        a.set_xlabel(xl, fontsize=fk.SMALL)
        a.grid(False)
        a.text(0.5, 1.0, f"ρ = {r:.2f}", transform=a.transAxes, ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT)
        if k == 0:
            a.set_ylabel("Postmortem minus freshly\nresected cortex (t)")
            a.set_title("pH and agonal axes compared to the\npattern of Liharska et al.\n", loc="left", pad=4)
        else:
            a.set_yticklabels([])
        out.append((a, r))
    a = out[1][0]
    from matplotlib.lines import Line2D
    a.legend(handles=[Line2D([], [], marker="o", color=fk.UP, lw=0, markersize=3, label="higher"),
                      Line2D([], [], marker="o", color=fk.DOWN, lw=0, markersize=3, label="lower")],
             frameon=False, fontsize=fk.TINY, loc="lower right", handletextpad=0.1, borderaxespad=0.2, title="In schizophrenia", title_fontsize=fk.TINY)
    return out[0][1], out[1][1], out[0][0]


def main():
    fig = fk.figure(232)
    r1 = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.15, 0.6], wspace=1.0, left=0.185, right=0.975, top=0.945, bottom=0.76)
    r2 = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.0], wspace=0.6, left=0.20, right=0.975, top=0.635, bottom=0.44)
    r3 = GridSpec(1, 2, figure=fig, width_ratios=[0.55, 1.45], wspace=0.3, left=0.17, right=0.975, top=0.315, bottom=0.065)
    axa, axb, axd = (fig.add_subplot(r1[i]) for i in range(3))
    axe, axf = (fig.add_subplot(r2[i]) for i in range(2))
    axg, axh = (fig.add_subplot(r3[i]) for i in range(2))
    gtex = panel_c(axa)
    panel_d(axb)
    agon = panel_f(axd)
    panel_survey(axe)
    panel_g(axf)
    rho_tot = panel_levels(axg)
    rho_ph, rho_ag, axh = panel_resected(axh)
    pf, ph_ = axg.get_position(), axh.get_position()  # f and g share the same baseline
    axg.set_position([pf.x0, ph_.y0, pf.width, pf.y1 - ph_.y0])
    for ax, l in zip([axa, axb, axd, axe, axf, axg, axh], "abcdefg"):
        fk.letter(ax, l)
    fig.texts[-1].set_y(fig.texts[-2].get_position()[1])  # f and g letters on one line
    fk.save(fig, "fig_dying", "figures/scripts/build_fig_dying.py",
            sorted(D8.glob("*.tsv")) + [DI / "fig5g_marker_sets_vif.tsv", D4D / "figure_data_fig4d_v4.tsv", D4D / "figure_data_fig4d_v4_summary.tsv",
                                         DC / "levels_summary_main_rows.tsv", DC / "sim_fraction_of_ceiling.tsv",
                                         D8 / "fig08_h_liharska_per_gene_383.tsv"],
            notes=f"a median {gtex[0]:.3f}, {gtex[1]}/{gtex[2]} negative; "
                  f"c median {agon[0]:.3f}, {agon[1]}/{agon[2]} positive; f whole-sample rho {rho_tot:.3f}; g Liharska rho pH {rho_ph:.2f}, agonal {rho_ag:.2f}")


if __name__ == "__main__":
    main()
