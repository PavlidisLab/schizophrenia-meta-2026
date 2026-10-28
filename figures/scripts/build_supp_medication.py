#!/usr/bin/env python3
"""Supplementary Fig. 7 -- antipsychotic treatment and the 176 recurrent genes.

Data: figures/data/supp_fig_medication/ (export_supp_fig_medication.py, 2026-10-09;
final for rodent and cell results). Only the 176-gene set (study_brain, uniform q < 0.1) is drawn.
a  direction score per chronic rodent unit against its matched-null mean, and pooled (a_unit_tests.tsv)
b  pooled chronic direction p with one unit left out, rows in the order of a (b_leave_one_out.tsv)
c  share of drug-changed genes (p < 0.05) that change in the same direction as in schizophrenia: rodent and human cell classes (a_unit_tests.tsv),
   dose within schizophrenia cases, GSE21138 (h_dose_within_sz.tsv)
d  per-gene mean drug z over the chronic units, by the gene's schizophrenia direction (c_gene_effects.tsv)
Unit drugs and species for the row labels come from e_units.tsv.
"""
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

import figkit as fk

D = fk.FD / "supp_fig_medication"
FILES = {k: D / f for k, f in {"a": "a_unit_tests.tsv", "b": "b_leave_one_out.tsv", "c": "c_gene_effects.tsv",
                                 "e": "e_units.tsv", "h": "h_dose_within_sz.tsv"}.items()}
rd = lambda k: pd.read_csv(FILES[k], sep="\t", comment="#", dtype={"enrichment_p": str, "concordance_p_depletion": str})
STRICT = "recurrent (study_brain, uniform q < 0.1)"
POOLED = "ALL chronic (15 units)"
SIG = fk.ACCENT
SAME, OPP = fk.ACCENT, fk.NEUTRAL


def fmt_p(p):
    return f"{p:.2g}" if p >= 0.001 else f"{p:.1e}"


def unit_rows():
    """Chronic units of the 176-gene set, sorted by observed minus null; labels with drug and species."""
    a = rd("a")
    a = a[(a["set"] == STRICT) & (a["class"] == "chronic") & (a["unit"] != POOLED)].copy()
    e = rd("e").set_index("unit")

    def lab(u):
        r = e.loc[u]
        drugs = []
        for d in str(r["drugs"]).split(";"):
            d = d.strip().split(" ")[0]
            if d not in drugs:
                drugs.append(d)
        geo, region = u.split("_", 1)
        region = region.replace("_", " ").replace("3regions", "3 regions").replace("PFC", "prefrontal")
        return f"{geo} {region}, {r['species']}\n{', '.join(drugs)}"

    a["label"] = a["unit"].map(lab)
    return a.sort_values("dir_minus_null").reset_index(drop=True)


def panel_a(ax, rows):
    a = rd("a")
    pooled = a[(a["set"] == STRICT) & (a["unit"] == POOLED)].iloc[0]
    y = np.arange(len(rows))
    for yi, (_, r) in zip(y, rows.iterrows()):
        ax.plot([r["dir_null_mean"], r["dir_obs"]], [yi, yi], color=fk.GRID, lw=1.2, zorder=1)
        ax.scatter(r["dir_null_mean"], yi, s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, zorder=2)
        ax.scatter(r["dir_obs"], yi, s=18, color=SIG, zorder=3)  # trends, not coloured by P
    yp = len(rows) + 0.6
    ax.plot([pooled["dir_null_mean"], pooled["dir_obs"]], [yp, yp], color=fk.GRID, lw=1.2, zorder=1)
    ax.scatter(pooled["dir_null_mean"], yp, s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, zorder=2)
    ax.scatter(pooled["dir_obs"], yp, s=26, marker="D", color=fk.TEXT, zorder=3)
    ax.text(pooled["dir_obs"] + 0.06, yp, f"{pooled['dir_obs']:.2f} vs {pooled['dir_null_mean']:.2f}",
            va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels(list(rows["label"]) + ["all 15 pooled"], fontsize=fk.TINY)
    ax.set_ylim(yp + 0.8, -0.8)
    ax.set_xlim(-0.45, 1.0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Direction score\n(mean SZ sign x drug z; > 0, same direction as SZ)")
    ax.set_title("Direction score of the recurrent genes in\nchronic rodent experiments", loc="left")
    ax.scatter([], [], s=18, color=SIG, label="observed")
    ax.scatter([], [], s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, label="matched random genes (mean)")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="upper right", bbox_to_anchor=(1.0, 0.95), handletextpad=0.2, borderaxespad=0.1)
    ax.text(0.99, 0.995, "one row = one chronic rodent experiment", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_b(ax, rows):
    b = rd("b")
    b = b[b["set"] == STRICT].set_index("left_out")
    pooled_p = rd("a").query("set == @STRICT and unit == @POOLED")["dir_p"].iloc[0]
    y = np.arange(len(rows))
    for yi, u in zip(y, rows["unit"]):
        p = b.loc[u, "dir_p"]
        above = p > 0.05
        ax.barh(yi, -np.log10(p), 0.6, color=fk.NEUTRAL)
        ax.text(-np.log10(p) + 0.05, yi, fmt_p(p), va="center", fontsize=fk.SMALL, color=fk.TEXT)
    yp = len(rows) + 0.6
    ax.barh(yp, -np.log10(pooled_p), 0.6, color=fk.TEXT)
    ax.text(-np.log10(pooled_p) + 0.05, yp, f"{fmt_p(pooled_p)} (none left out)", va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(-np.log10(0.05), color=fk.SUBTLE, lw=0.6, ls="--")
    ax.text(-np.log10(0.05) - 0.03, len(rows) - 0.2, "P 0.05", fontsize=fk.SMALL, color=fk.SUBTLE, ha="right", va="center")
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels([])
    ax.set_ylim(yp + 0.8, -0.8)
    ax.set_xlim(0, 2.9)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Pooled P with that experiment left out\n(−log10)")
    ax.set_title("Pooled P with each\nexperiment left out", loc="left")


def panel_c(ax):
    a = rd("a").query("set == @STRICT").set_index("unit")
    h = rd("h").set_index("set").loc["study_brain"]
    rows = []
    for u, l in [("ALL chronic (15 units)", "rodent brain, chronic\n(15 experiments)"), ("ALL acute (2 units)", "rodent brain, acute\n(2 experiments)"),
                 ("ALL cells (3 units)", "human neural cells\n(3 experiments)")]:
        r = a.loc[u]
        rows.append((l, int(r["n_sig_same_as_sz"]), int(r["n_sig_opposite_sz"]),
                     f"score {r['dir_obs']:.2f}\nvs {r['dir_null_mean']:.2f} matched"))
    same = int(h["n_p05_same_sign_as_sz"])
    rows.append(("dose within SZ cases\n(GSE21138, 26 cases)", same, int(h["n_set_p05"]) - same,
                 f"{int(h['n_set_p05'])} genes vs\n{h['expected_p05']:.1f} expected"))
    y = np.arange(len(rows))
    for yi, (_, r) in zip(y, rows.iterrows()):
        ax.plot([r["dir_null_mean"], r["dir_obs"]], [yi, yi], color=fk.GRID, lw=1.2, zorder=1)
        ax.scatter(r["dir_null_mean"], yi, s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, zorder=2)
        ax.scatter(r["dir_obs"], yi, s=18, color=SIG, zorder=3)  # trends, not coloured by P
    yp = len(rows) + 0.6
    ax.plot([pooled["dir_null_mean"], pooled["dir_obs"]], [yp, yp], color=fk.GRID, lw=1.2, zorder=1)
    ax.scatter(pooled["dir_null_mean"], yp, s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, zorder=2)
    ax.scatter(pooled["dir_obs"], yp, s=26, marker="D", color=fk.TEXT, zorder=3)
    ax.text(pooled["dir_obs"] + 0.06, yp, f"{pooled['dir_obs']:.2f} vs {pooled['dir_null_mean']:.2f}",
            va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels(list(rows["label"]) + ["all 15 pooled"], fontsize=fk.TINY)
    ax.set_ylim(yp + 0.8, -0.8)
    ax.set_xlim(-0.45, 1.0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Direction score\n(mean SZ sign x drug z; > 0, same direction as SZ)")
    ax.set_title("Direction score of the recurrent genes in\nchronic rodent experiments", loc="left")
    ax.scatter([], [], s=18, color=SIG, label="observed")
    ax.scatter([], [], s=14, facecolor="white", edgecolor=fk.SUBTLE, lw=0.7, label="matched random genes (mean)")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="upper right", bbox_to_anchor=(1.0, 0.95), handletextpad=0.2, borderaxespad=0.1)
    ax.text(0.99, 0.995, "one row = one chronic rodent experiment", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_b(ax, rows):
    b = rd("b")
    b = b[b["set"] == STRICT].set_index("left_out")
    pooled_p = rd("a").query("set == @STRICT and unit == @POOLED")["dir_p"].iloc[0]
    y = np.arange(len(rows))
    for yi, u in zip(y, rows["unit"]):
        p = b.loc[u, "dir_p"]
        above = p > 0.05
        ax.barh(yi, -np.log10(p), 0.6, color=fk.NEUTRAL)
        ax.text(-np.log10(p) + 0.05, yi, fmt_p(p), va="center", fontsize=fk.SMALL, color=fk.TEXT)
    yp = len(rows) + 0.6
    ax.barh(yp, -np.log10(pooled_p), 0.6, color=fk.TEXT)
    ax.text(-np.log10(pooled_p) + 0.05, yp, f"{fmt_p(pooled_p)} (none left out)", va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(-np.log10(0.05), color=fk.SUBTLE, lw=0.6, ls="--")
    ax.text(-np.log10(0.05) - 0.03, len(rows) - 0.2, "P 0.05", fontsize=fk.SMALL, color=fk.SUBTLE, ha="right", va="center")
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels([])
    ax.set_ylim(yp + 0.8, -0.8)
    ax.set_xlim(0, 2.9)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Pooled P with that experiment left out\n(−log10)")
    ax.set_title("Pooled P with each\nexperiment left out", loc="left")


def panel_c(ax):
    a = rd("a").query("set == @STRICT").set_index("unit")
    h = rd("h").set_index("set").loc["study_brain"]
    rows = []
    for u, l in [("ALL chronic (15 units)", "rodent brain, chronic\n(15 experiments)"), ("ALL acute (2 units)", "rodent brain, acute\n(2 experiments)"),
                 ("ALL cells (3 units)", "human neural cells\n(3 experiments)")]:
        r = a.loc[u]
        rows.append((l, int(r["n_sig_same_as_sz"]), int(r["n_sig_opposite_sz"]),
                     f"score {r['dir_obs']:.2f}\nvs {r['dir_null_mean']:.2f} matched"))
    same = int(h["n_p05_same_sign_as_sz"])
    rows.append(("dose within SZ cases\n(GSE21138, 26 cases)", same, int(h["n_set_p05"]) - same,
                 f"{int(h['n_set_p05'])} genes vs\n{h['expected_p05']:.1f} expected"))
    y = np.arange(len(rows)) + np.array([0, 0, 0, 0.3])
    for yi, (l, s, o, note) in zip(y, rows):
        t = s + o
        ax.barh(yi, s / t, 0.62, color=SAME)
        ax.barh(yi, o / t, 0.62, left=s / t, color=OPP)
        ax.text(0.02, yi, f"{s}", va="center", ha="left", fontsize=fk.SMALL, color="white")
        ax.text(0.98, yi, f"{o}", va="center", ha="right", fontsize=fk.SMALL, color=fk.TEXT)
        ax.text(1.04, yi, note, va="center", ha="left", fontsize=fk.TINY, color=fk.TEXT, transform=ax.get_yaxis_transform())
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=fk.TINY)
    ax.set_ylim(y[-1] + 0.6, -0.6)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Genes changed by drug or dose at P < 0.05\n(share in the same direction as SZ; numbers, gene counts)")
    ax.set_title("Direction of drug-changed genes in other\nsettings", loc="left", pad=16)
    ax.legend(handles=[Patch(color=SAME, label="same direction as SZ"), Patch(color=OPP, label="opposite direction")], frameon=False,
              fontsize=fk.SMALL, loc="lower left", bbox_to_anchor=(-0.02, 0.99), ncol=2, handlelength=1.0, borderaxespad=0.1)


def panel_d(ax):
    c = rd("c")
    c = c[(c["set"] == STRICT) & (c["sz_direction"].isin(["down", "up"])) & c["chronic_mean_z"].notna()]
    rng = np.random.default_rng(1)
    label = {"SERPINA3", "HBA2", "NNAT", "ATP2B2", "SLCO1C1", "WSB1"}
    for i, d in enumerate(["down", "up"]):
        g = c[c["sz_direction"] == d]
        z = g["chronic_mean_z"].to_numpy()
        szway = (z < 0) if d == "down" else (z > 0)
        x = i + rng.uniform(-0.18, 0.18, len(z))
        ax.scatter(x[szway], z[szway], s=7, color=SAME, alpha=0.85, lw=0)
        ax.scatter(x[~szway], z[~szway], s=7, color=OPP, alpha=0.85, lw=0)
        ax.text(i, 1.85, f"{szway.sum()} of {len(z)}\nsame as SZ", ha="center", fontsize=fk.SMALL, color=fk.TEXT)
        for xi, zi, gn in zip(x, z, g["gene"]):
            if gn in label:
                ax.text(xi + (0.06 if xi >= i else -0.06), zi, gn, fontsize=fk.TINY, color=fk.TEXT, va="center",
                        ha="left" if xi >= i else "right")
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["lower in SZ", "higher in SZ"], fontsize=fk.SMALL)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-1.6, 2.3)
    ax.set_ylabel("Drug effect on the gene\n(mean z over chronic rodent experiments)")
    ax.set_title("Drug effect on each recurrent gene", loc="left")
    ax.text(0.99, 0.02, "one point = one gene", transform=ax.transAxes, ha="right", va="bottom", fontsize=fk.SMALL, color=fk.SUBTLE)


def main():
    fk.apply_style()
    fig = fk.figure(175)
    rows = unit_rows()
    # leave-one-out panel (old b) removed: the text calls the direction result preliminary
    ax_a = fig.add_axes([0.20, 0.50, 0.40, 0.44])
    ax_d = fig.add_axes([0.73, 0.50, 0.25, 0.44])
    ax_c = fig.add_axes([0.20, 0.07, 0.25, 0.27])
    panel_a(ax_a, rows)
    panel_d(ax_d)
    panel_c(ax_c)
    for ax, l in zip((ax_a, ax_d, ax_c), "abc"):
        fk.letter(ax, l)
    fk.save(fig, "supp_medication", "figures/scripts/build_supp_medication.py", list(FILES.values()),
            notes="Final for rodent and cell results. No multiplicity adjustment across sets or classes.")


if __name__ == "__main__":
    main()
