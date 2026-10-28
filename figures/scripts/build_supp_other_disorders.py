#!/usr/bin/env python3
"""Extended Fig. 7 -- the 209 recurrent schizophrenia genes in other disorders' postmortem brain.

Data: figures/data/supp_fig_other_disorders/
(export_supp_fig_other_disorders.py, 2026-10-09, from analyses/other_disorders_176_2026-10-09/, seed 20261009).
a  number of the 176 among DE genes against the matched null, per table (panel_a_enrichment_vs_null.tsv)
b  share of the 176 in the schizophrenia direction, among DE genes and among all measured genes, against the matched null
   (panel_b_concordance_vs_null.tsv)
c  z of the Spearman rho between the disorder's signed z and the pH and agonal axes, against the matched null (panel_e_axis_z.tsv)
d  Hagihara 2018 Table S1 pH difference from controls, SZ and BD (panel_d_hagihara_ph_by_diagnosis.tsv)
Not drawn: per-list Gemma counts (panel_c_gemma_per_list.tsv); the numbers are quoted in the legend draft.
Rows of a to c: one table each. DE is FDR < 0.05, except Gemma autism (FDR < 0.10; at 0.05 one pool gene in 15 is one of the 176).
Gemma bipolar, depression and Alzheimer rows are the pooled lists with at most 18.3% of genes DE.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import figkit as fk

D = fk.FD / "supp_fig_other_disorders"
# a to c: effects with intervals that allow for correlation among genes; d is unchanged
DC = fk.FD / "supp_fig_other_disorders_corrected"
FILES = {"a": DC / "panel_a_enrichment_vs_null.tsv", "b": DC / "panel_b_concordance_vs_null.tsv",
         "d": D / "panel_d_hagihara_ph_by_diagnosis.tsv"}
rd = lambda k: pd.read_csv(FILES[k], sep="\t", comment="#")
PS = {k: pd.read_csv(DC / f"panel_{k}_per_study.tsv", sep="\t", comment="#") for k in "abc"}
PS_DOT = 5  # marker area of one study

GEMMA_BD = "Gemma BD pooled [restricted, 7 lists] (all lists)"
GEMMA_BD_ALL = "Gemma BD pooled [every list, 28] (lists with <= 18.3% DE)"
GEMMA_MDD = "Gemma MDD pooled (lists with <= 18.3% DE)"
GEMMA_ASD = "Gemma ASD pooled [restricted, 6 lists] (all lists)"
GEMMA_AD = "Gemma AD pooled, independent contrasts (n=28) (lists with <= 18.3% DE)"

# (group header, [(table, label, fdr)]) -- main rows use only series sharing no donors with schizophrenia data;
# tables dropped by that rule are grouped as a comparison
GROUPS = [
    ("Schizophrenia (for comparison)", [("G18-array SCZ", "Gandal 2018, microarray", 0.05),
                                               ("G18-RNAseq SCZ", "PsychENCODE, RNA-seq", 0.05)]),
    ("Bipolar disorder, no shared donors", [(GEMMA_BD, "Gemma, 7 lists", 0.05)]),
    ("Autism, no shared donors", [("G18-array ASD", "Gandal 2018, microarray", 0.05),
                                  (GEMMA_ASD, "Gemma, 6 lists (FDR 0.10)", 0.1)]),
    ("Shared donors (comparison)", [("G18-array BD", "bipolar, Gandal 2018 microarray", 0.05),
                                    ("G18-RNAseq BD", "bipolar, PsychENCODE RNA-seq", 0.05),
                                    (GEMMA_BD_ALL, "bipolar, Gemma, 28 lists", 0.05),
                                    ("G18-RNAseq ASD", "autism, PsychENCODE RNA-seq", 0.05)]),
    ("Alzheimer disease", [(GEMMA_AD, "Gemma, 28 contrasts", 0.05)]),
    ("Major depression", [("G18-array MDD", "Gandal 2018, microarray", 0.05),
                          (GEMMA_MDD, "Gemma, pooled lists", 0.05)]),
]
REF = fk.TEXT
OBS = fk.ACCENT
NULLC = fk.SUBTLE


def layout_rows():
    """y position per table, header y per group (rows top to bottom, a gap and a header line per group)."""
    rows, headers, y = [], [], 0.0
    for g, items in GROUPS:
        headers.append((y, g))
        y += 1.0
        for t, lab, fdr in items:
            rows.append(dict(table=t, label=lab, fdr=fdr, y=y, ref=g.startswith("Schizophrenia")))
            y += 1.0
        y += 0.35
    return pd.DataFrame(rows), headers, y - 0.35


def fmt_p(p):
    return f"{p:.2g}" if p >= 0.001 else f"{p:.0e}".replace("e-0", "e-")


def _rows_axis(ax, rows, headers, ymax, labels):
    ax.set_ylim(ymax + 0.1, -0.7)
    ax.set_yticks(rows["y"])
    ax.set_yticklabels(rows["label"] if labels else [], fontsize=fk.SMALL)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    for yh, g in headers:
        if labels:
            ax.text(-0.02, yh, g, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=fk.SMALL,
                    color=fk.TEXT, fontweight="bold")
        ax.axhline(yh - 0.55, color=fk.GRID, lw=0.5) if yh > 0 else None


def _effect_row(ax, e, y, ref, lo_col="lo95_corrected", hi_col="hi95_corrected", clear_col="excludes_zero_corrected",
                marker="o", s=20, dy=0.0, colour=None):
    """effect with its correlation-aware 95% interval; blue where the interval excludes zero, grey where not"""
    clear = str(e[clear_col]).lower() in ("yes", "true")
    col = colour if colour is not None else (REF if ref else OBS)  # trends only, not coloured by interval
    ax.plot([e[lo_col], e[hi_col]], [y + dy] * 2, color=col, lw=1.0, solid_capstyle="butt", zorder=2,
            alpha=1.0 if (clear or ref or colour is not None) else 0.8)
    ax.scatter(e["effect"], y + dy, s=s, marker=marker, color=col, zorder=3, lw=0)


def members(table, which):
    """per-study rows (383 set) of one figure row: the table itself when it is a single table, else the lists pooled into it"""
    d = PS[which]
    d = d[d["set"] == 383]
    pooled = d["pooled_into"].fillna("").str.split(" \\| ")
    hit = d[(d["table"] == table) | pooled.apply(lambda l: table in l)]
    return hit


def _offsets(vals, width):
    """spread dots over the row height by rank, so overlapping studies stay visible"""
    order = np.argsort(np.argsort(vals))
    levels = np.linspace(-width, width, 5)
    return levels[order % 5]


def _dots(ax, x, y, ref, colour=None, width=0.3):
    x = np.asarray(x, float)
    if not len(x):
        return
    col = colour if colour is not None else (REF if ref else OBS)
    ax.scatter(x, y + _offsets(x, width), s=PS_DOT, color=col, alpha=0.55, lw=0, zorder=2)


def _pooled(ax, x, y, ref, colour=None, marker="D", s=22):
    col = colour if colour is not None else (REF if ref else OBS)
    ax.scatter([x], [y], s=s, marker=marker, color=col, edgecolor="white", linewidth=0.5, zorder=4)


def panel_a(ax, rows, headers, ymax):
    a = rd("a").set_index(["table", "fdr_thr"])
    xs = [0.0]
    for _, r in rows.iterrows():
        e = a.loc[(r["table"], r["fdr"])]
        m = members(r["table"], "a")
        m = m[m["fdr"] == r["fdr"]]
        single = len(m) <= 1
        if not single:
            _dots(ax, m["effect"], r["y"], r["ref"])
        _pooled(ax, e["effect"], r["y"], r["ref"], marker="o" if single else "D", s=18 if single else 22)
        xs += list(m["effect"]) + [e["effect"]]
        ax.text(0.99, r["y"], f"{int(e['observed'])} vs {e['matched_expected']:.1f}", transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=fk.TINY, color=fk.TEXT)
    _rows_axis(ax, rows, headers, ymax, labels=True)
    lo, hi = min(xs), max(xs)
    span = hi - lo
    ax.set_xlim(lo - 0.04 * span, hi + 0.30 * span)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xlabel("Schizophrenia-recurrent genes among the DE genes,\nobserved minus expected by chance (count)")
    ax.set_title(f"Are the {fk.n_recurrent()} schizophrenia-recurrent genes\namong each disorder's DE genes?", loc="left", pad=4)


def panel_b(ax_de, ax_all, rows, headers, ymax):
    b = rd("b").set_index(["table", "set"])
    for ax, which in ((ax_de, "de"), (ax_all, "all")):
        for _, r in rows.iterrows():
            st = ("DE q<0.10" if r["fdr"] == 0.1 else "DE q<0.05") if which == "de" else "all measured"
            e = b.loc[(r["table"], st)]
            if pd.isna(e["lo95_corrected"]) or not e["n_with_direction"]:
                ax.text(0.03, r["y"], "no DE genes", transform=ax.get_yaxis_transform(), va="center", ha="left",
                        fontsize=fk.TINY, color=fk.SUBTLE)
                continue
            m = members(r["table"], "b")
            m = m[(m["subset"] == st) & (m["n_with_direction"] > 0)]
            single = len(m) <= 1
            if not single:
                _dots(ax, m["effect"], r["y"], r["ref"])
            _pooled(ax, e["effect"], r["y"], r["ref"], marker="o" if single else "D", s=18 if single else 22)
            ax.text(0.02, r["y"] - 0.4, f"{int(e['concordant'])} of {int(e['n_with_direction'])}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=fk.TINY, color=fk.SUBTLE)
        _rows_axis(ax, rows, headers, ymax, labels=False)
        ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax_de.set_xlabel("Share in SZ direction,\nminus matched\n(DE genes)")
    ax_all.set_xlabel("Share in SZ direction,\nminus matched\n(all measured)")
    ax_de.set_title(f"Do they change in the same direction\nas in schizophrenia?", loc="left", pad=4)


def panel_d(ax):
    d = rd("d")
    per = d[d["model"].str.startswith("per dataset")].copy()
    per["dataset"] = per["model"].str.replace("per dataset: ", "", regex=False).str.replace("Dean et al, 2016", "Dean 2016")
    order = sorted(per["dataset"].unique(), key=lambda s: (s.startswith("SMRI"), s.startswith("NSW"), s))
    col = {"SZ": REF, "BD": OBS}
    off = {"SZ": -0.18, "BD": 0.18}
    labels, ys = [], []
    for i, ds in enumerate(order):
        sub = per[per["dataset"] == ds]
        for _, r in sub.iterrows():
            ax.scatter(r["coef_pH_units"], i + off[r["dx"]], s=14, color=col[r["dx"]], zorder=3, lw=0)
        labels.append(ds)
        ys.append(i)
    y0 = len(order) + 0.6
    pooled = [("all rows, dataset fixed effect", "all 10 datasets, pooled"),
              ("datasets with SZ, BD and controls (GSE35977, GSE53987, SMRI Collection A, SMRI Collection C)",
               "4 datasets with both, pooled")]
    for k, (mk, lab) in enumerate(pooled):
        yk = y0 + k * 1.1
        for dx in ("SZ", "BD"):
            r = d[(d["model"] == mk) & (d["dx"] == dx)].iloc[0]
            ax.plot([r["ci_low"], r["ci_high"]], [yk + off[dx]] * 2, color=col[dx], lw=1.0, zorder=2)
            ax.scatter(r["coef_pH_units"], yk + off[dx], s=24, marker="D", color=col[dx], zorder=3, lw=0)
        labels.append(lab)
        ys.append(yk)
    yk = y0 + 2 * 1.1
    r = d[d["dx"] == "SZ-BD"].iloc[0]
    ax.plot([r["ci_low"], r["ci_high"]], [yk] * 2, color=fk.ACCENT_3, lw=1.0, zorder=2)
    ax.scatter(r["coef_pH_units"], yk, s=24, marker="D", color=fk.ACCENT_3, zorder=3, lw=0)
    ax.text(r["ci_low"] - 0.01, yk, f"{r['coef_pH_units']:.3f}, P {r['p']:.2f}", va="center", ha="right",
            fontsize=fk.TINY, color=fk.TEXT)
    labels.append("schizophrenia minus bipolar,\nsame 4 datasets")
    ys.append(yk)
    ax.axhline(len(order) - 0.2, color=fk.GRID, lw=0.5)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=fk.SMALL)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(yk + 0.7, -0.7)
    ax.set_xlim(-0.40, 0.16)
    ax.set_xticks([-0.4, -0.3, -0.2, -0.1, 0, 0.1])
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Brain pH, case minus control\n(pH units; Hagihara 2018 Table S1)")
    ax.set_title("Brain pH, schizophrenia and bipolar disorder\nagainst controls, by dataset", loc="left", pad=4)
    ax.legend(handles=[Line2D([], [], marker="o", color=REF, lw=0, markersize=3.6, label="schizophrenia"),
                       Line2D([], [], marker="o", color=OBS, lw=0, markersize=3.6, label="bipolar disorder"),
                       Line2D([], [], marker="D", color=fk.SUBTLE, lw=1.0, markersize=3.6, label="pooled, 95% CI")],
              frameon=False, fontsize=fk.SMALL, loc="upper left", bbox_to_anchor=(-0.02, -0.17), handletextpad=0.2,
              borderaxespad=0.1, labelspacing=0.25, ncol=2, columnspacing=1.0)
    ax.text(1.0, -0.245, "one point = one dataset", transform=ax.transAxes, ha="right", va="center",
            fontsize=fk.SMALL, color=fk.SUBTLE)


def legend_ab(fig):
    h = [Line2D([], [], marker="o", color=OBS, lw=1.0, markersize=4, label="other disorders (observed minus matched random genes)"),
         Line2D([], [], marker="o", color=REF, lw=1.0, markersize=4, label="schizophrenia, for comparison"),
         Line2D([], [], marker="o", color=fk.SUBTLE, lw=0, markersize=2.5, label="one study"),
         Line2D([], [], marker="D", color=fk.SUBTLE, lw=0, markersize=4, label="pooled")]
    fig.legend(handles=h, frameon=False, fontsize=fk.SMALL, loc="upper left", bbox_to_anchor=(0.23, 0.535), ncol=4,
               handletextpad=0.3, columnspacing=1.2)


def main():
    fk.apply_style()
    fig = fk.figure(215)
    rows, headers, ymax = layout_rows()
    top, h = 0.61, 0.33
    ax_a = fig.add_axes([0.235, top, 0.33, h])
    ax_b1 = fig.add_axes([0.615, top, 0.155, h])
    ax_b2 = fig.add_axes([0.82, top, 0.155, h])
    ax_d = fig.add_axes([0.235, 0.135, 0.40, 0.29])
    panel_a(ax_a, rows, headers, ymax)
    panel_b(ax_b1, ax_b2, rows, headers, ymax)
    panel_d(ax_d)
    legend_ab(fig)
    for ax, l in zip((ax_a, ax_b1), "ab"):
        fk.letter(ax, l)
    fk.letter(ax_d, "c", dx=-0.06)
    fk.save(fig, "supp_other_disorders", "figures/scripts/build_supp_other_disorders.py", list(FILES.values()) + [DC / f"panel_{k}_per_study.tsv" for k in "ab"],
            notes="a-c corrected intervals; d re-exported 2026-10-12. Main rows use series sharing no donors with schizophrenia data; dropped tables shown as a comparison. The 209 are not independent of Gandal 2018 (Stanley and SMRI array collections, "
                  "BrainGVEX) or of the Gemma BD and MDD lists (Stanley); SZ rows are a reference. Gemma BD and MDD rows are the "
                  "pooled lists with at most 18.3% DE; Gemma autism row is FDR 0.10. ")


if __name__ == "__main__":
    main()
