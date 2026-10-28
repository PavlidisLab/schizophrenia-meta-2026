#!/usr/bin/env python3
"""Supplementary Note 5 Figure 3 -- SNAP's astrocyte and neuron components against pH and agonal state.

Data: data/snap_ph_parts/ (snap_components_sn.tsv, snap_components_bulk.tsv, snap_vs_ph.tsv).
a  correlation of SNAP-a, SNAP-n and SNAP total with the pH-axis score, measured pH and the estimated agonal score,
   in 610 single-nucleus donors (3 banks), the Frohlich controls with pH, and bulk datasets
b  schizophrenia effect on each component (SD), in the base model and with the agonal score or the pH-axis score added
"""
import numpy as np
import pandas as pd

import figkit as fk

ND = fk.note_data(__file__)
D = ND / "snap_ph_parts"
FILES = [D / "snap_components_sn.tsv", D / "snap_components_bulk.tsv", D / "snap_vs_ph.tsv"]
SN = pd.read_csv(FILES[0], sep="\t", comment="#")
BK = pd.read_csv(FILES[1], sep="\t", comment="#")
FR = pd.read_csv(FILES[2], sep="\t", comment="#")
PARTS = [("SNAP-a", fk.ACCENT_3, "SNAP-a (astrocyte)"), ("SNAP-n", fk.ACCENT, "SNAP-n (neuron)"), ("SNAP total", fk.NEUTRAL, "SNAP total")]


def sn_corr(score, against):
    r = SN[(SN["table"] == "correlation") & SN["set"].str.startswith("pooled") & (SN["score"] == score)
           & (SN["against"] == against) & (SN["statistic"] == "spearman")].iloc[0]
    return r["estimate"], r["lo95"], r["hi95"]


def bk(score, analysis_start, dataset_start):
    r = BK[(BK["table"] == "bulk_pooled") & BK["dataset"].str.startswith(dataset_start) & (BK["score"] == score)
           & BK["analysis"].str.startswith(analysis_start)]
    r = r[r["dataset"].str.contains("all datasets")].iloc[0]
    return r["estimate"], r["lo95"], r["hi95"]


def fr(score, against):
    r = FR[(FR["set"] == "NSW") & (FR["group"] == "controls") & (FR["snap"] == score) & (FR["against"] == against)
           & (FR["statistic"] == "spearman_rho")].iloc[0]
    return r["estimate"], r["lo95"], r["hi95"]


def sz_sn(score, model):
    r = SN[(SN["table"] == "sz_effect") & SN["set"].str.startswith("pooled") & (SN["score"] == score) & (SN["model"] == model)].iloc[0]
    return r["estimate"], r["lo95"], r["hi95"]


def sz_bk(score, model):
    return bk(score, f"SZ effect, model {model} (age, sex)", "pooled over 7")


def rows_dots(ax, rows, parts, offsets):
    for yi, row in enumerate(rows):
        for (key, col, _), off in zip(parts, offsets):
            v, lo, hi = row[1][key]
            ax.errorbar(v, yi + off, xerr=[[v - lo], [hi - v]], fmt="o", ms=3, color=col, elinewidth=0.8, capsize=0)


def panel_a(ax):
    rows = []
    for lab, f in [
        ("pH-axis score\n(single-nucleus, 610)", lambda s: sn_corr(s, "pH-axis score")),
        ("pH-axis score\n(bulk, 22 datasets)", lambda s: bk(s, "all samples: pH-axis score", "pooled over 7")),
        ("measured pH, controls\n(Fröhlich, 33)", lambda s: fr(s, "measured pH")),
        ("measured pH, controls\n(bulk, 13 datasets)", lambda s: bk(s, "controls: measured pH (higher", "pooled over 5")),
        ("agonal score\n(single-nucleus, 610)", lambda s: sn_corr(s, "estimated agonal score")),
        ("agonal score\n(bulk, 22 datasets)", lambda s: bk(s, "all samples: estimated agonal score", "pooled over 7")),
    ]:
        rows.append((lab, {p[0]: f(p[0]) for p in PARTS}))
    rows_dots(ax, rows, PARTS, [-0.22, 0, 0.22])
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=fk.SMALL)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    for y in (1.5, 3.5):
        ax.axhline(y, color=fk.GRID, lw=0.6)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel("Correlation (Spearman ρ)")
    ax.set_title("SNAP components against pH and agonal scores", loc="left")
    ax.grid(axis="y", visible=False)


def panel_b(ax):
    models = [("base", "age, sex", fk.ACCENT), ("+agonal", "+ agonal score", fk.ACCENT_2), ("+pH axis", "+ pH-axis score", fk.ACCENT_3)]
    rows = []
    for dt, fn in (("single-nucleus", sz_sn), ("bulk", sz_bk)):
        for key, _, lab in [("SNAP-n", 0, "SNAP-n"), ("SNAP-a", 0, "SNAP-a"), ("SNAP total", 0, "SNAP total")]:
            rows.append((f"{lab}\n({dt})", {m[0]: fn(key, m[0]) for m in models}))
    parts = [(m[0], m[2], m[1]) for m in models]
    rows_dots(ax, rows, parts, [-0.22, 0, 0.22])
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=fk.SMALL)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.axhline(2.5, color=fk.GRID, lw=0.6)
    ax.set_xlabel("Schizophrenia effect on the score (SD)")
    ax.set_title("Schizophrenia effect on each component", loc="left")
    ax.grid(axis="y", visible=False)


def main():
    from matplotlib.lines import Line2D
    fk.apply_style()
    fig = fk.figure(105)
    A = fig.add_axes([0.20, 0.25, 0.27, 0.62])
    B = fig.add_axes([0.70, 0.25, 0.27, 0.62])
    panel_a(A)
    panel_b(B)
    mod = [("age, sex", fk.ACCENT), ("+ agonal score", fk.ACCENT_2), ("+ pH-axis score", fk.ACCENT_3)]
    for ax, items, ttl in ((A, [(l, c) for _, c, l in PARTS], None), (B, mod, "model:")):
        h = [Line2D([], [], marker="o", ls="", color=c, label=l, ms=3.5) for l, c in items]
        ax.legend(handles=h, frameon=False, fontsize=fk.SMALL, loc="upper center", bbox_to_anchor=(0.42, -0.26), ncol=3,
                  columnspacing=1.0, handletextpad=0.2)
    for ax, l in zip((A, B), "ab"):
        fk.letter(ax, l, dx=-0.02)
    fk.save(fig, "supp_snap_ph", "supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_snap_ph.py", FILES,
            notes="Panels: a correlations with pH and agonal scores, b schizophrenia effect by model.")


if __name__ == "__main__":
    main()
