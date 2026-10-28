#!/usr/bin/env python3
"""Regenerates fig_ieg_across_contrasts.{svg,png} and fig_ieg_across_contrasts_caption.txt from ieg_table.tsv and per_gene_7.tsv (no numbers typed here).
One row per contrast: coloured marks = the 7 rapid primary response genes' effects; black diamond = their mean; grey bar = 95% range of the means of
10,000 expression-decile-matched random sets (variant A). Dachet rows (log2FC) have their own axis.   python3 make_figure.py
"""
import hashlib, os, sys
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = pd.read_csv(f"{OUT}/ieg_table.tsv", sep="\t", comment="#")
G = pd.read_csv(f"{OUT}/per_gene_7.tsv", sep="\t", comment="#")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"], "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#6b7280", "axes.labelcolor": "#1f2937", "xtick.color": "#6b7280", "ytick.color": "#1f2937", "font.size": 7.5})
GENES = ["EGR1", "NR4A1", "DUSP1", "FOS", "MAFF", "PPP1R15A", "GADD45G"]
COL = {"EGR1": "#0072B2", "NR4A1": "#56B4E9", "DUSP1": "#009E73", "FOS": "#CC79A7", "MAFF": "#E69F00", "PPP1R15A": "#D55E00", "GADD45G": "#7A3E9D"}
szdir = G.drop_duplicates("gene").set_index("gene").sz_direction.to_dict()
S = T[T.set.str.startswith("S7") & T.null_variant.str.startswith("A")].set_index("contrast")
GROUPS = [("Yaqubi surgical holds", "Surgical holds (Yaqubi)"), ("Yaqubi autopsy", "Autopsy vs frozen (Yaqubi)"), ("Dachet time course (1 patient)", "Ex vivo time course (Dachet)"),
          ("Living vs postmortem", "Living vs postmortem (Liharska)"), ("GTEx brain regions", "GTEx slow vs fast death, brain"), ("GTEx other tissues", "GTEx slow vs fast death, other tissues"),
          ("Pritzker agonal score", "Pritzker agonal factor"), ("SZ case-control comparator", "SZ case-control, 22 bulk datasets")]
rows = []
for gk, gl in GROUPS:
    names = [c for c in G[G.contrast_group == gk].contrast.drop_duplicates() if not c.endswith("cortical mean (Cortex, BA9, BA24)") and "reproduction (original)" not in c]
    for c in names:
        rows.append((gk, gl, c))
n = len(rows)
fig = plt.figure(figsize=(8.6, 0.155 * n + 2.6))
gs = fig.add_gridspec(1, 2, width_ratios=[4.2, 1.15], wspace=0.04, left=0.30, right=0.985, top=1 - 0.5 / (0.155 * n + 2.6), bottom=0.9 / (0.155 * n + 2.6))
ax = fig.add_subplot(gs[0]); ad = fig.add_subplot(gs[1], sharey=ax)
ypos = {}; y = 0; prev = None; seps = []
for gk, gl, c in rows:
    if prev is not None and gk != prev:
        y += 1.7; seps.append((y - 0.9, gl))
    ypos[c] = y; y += 1; prev = gk
labels = {}
for gk, gl, c in rows:
    lab = c.replace("GTEx slow vs fast death: ", "").replace("Yaqubi: ", "").replace("Dachet: ", "").replace("Liharska: living vs postmortem ", "").replace("Pritzker agonal factor ", "agonal factor ")
    lab = lab.replace("SZ case-control: mean t over the 22 bulk datasets ", "mean t, 22 datasets ")
    labels[c] = lab
for gk, gl, c in rows:
    a = ad if gk.startswith("Dachet") else ax
    yy = ypos[c]; r = S.loc[c]
    a.barh(yy, r.null_hi95 - r.null_lo95, left=r.null_lo95, height=0.7, color="#e5e7eb", zorder=1, linewidth=0)
    for g in GENES:
        e = G[(G.contrast == c) & (G.gene == g)].effect.values[0]
        if np.isnan(e):
            continue
        a.scatter(e, yy, s=13, color=COL[g], marker="o" if szdir[g] == "down" else "^", zorder=3, linewidths=0)
    a.scatter(r.mean_stat, yy, s=26, color="black", marker="D", zorder=4, linewidths=0)
for a in (ax, ad):
    a.axvline(0, color="#6b7280", lw=0.6, zorder=0); a.grid(axis="x", color="#e5e7eb", lw=0.5); a.set_axisbelow(True)
ax.set_xscale("symlog", linthresh=3, linscale=1.0)
ax.set_xticks([-20, -10, -3, 0, 3, 10, 20]); ax.set_xticklabels(["-20", "-10", "-3", "0", "3", "10", "20"])
ax.set_yticks([ypos[c] for _, _, c in rows]); ax.set_yticklabels([labels[c] for _, _, c in rows], fontsize=6.5)
ax.set_ylim(y - 0.2, -1.2)
ax.set_xlabel("effect on the gene (t; symmetric-log scale)", fontsize=7.5)
ad.set_xlabel("log2 change from 0 h", fontsize=7.5); plt.setp(ad.get_yticklabels(), visible=False); ad.tick_params(axis="y", length=0); ad.spines["left"].set_visible(False)
# hide unused rows on each axis visually by lightly marking groups
first = {}
for gk, gl, c in rows:
    first.setdefault(gl, ypos[c])
for gl, yy in first.items():
    ax.text(-0.015, yy - 1.15, gl, transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=7, fontweight="bold", color="#1f2937")
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o" if szdir[g] == "down" else "^", color=COL[g], ls="", ms=4, label=f"{g} (SZ {szdir[g]})") for g in GENES]
h += [Line2D([], [], marker="D", color="black", ls="", ms=4, label="mean of the 7"), plt.Rectangle((0, 0), 1, 1, color="#e5e7eb", label="95% range, matched random sets")]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=6.5, bbox_to_anchor=(0.6, 0.012))
fig.text(0.02, 0.003, f"Source: ieg_table.tsv (sha256 {sha(f'{OUT}/ieg_table.tsv')}), per_gene_7.tsv (sha256 {sha(f'{OUT}/per_gene_7.tsv')}); scripts/make_figure.py", fontsize=6, color="#6b7280")
for a in (ax, ad):
    for art in list(a.patches) + list(a.lines) + list(a.texts) + list(a.collections):
        art.set_clip_on(False)
fig.savefig(f"{OUT}/fig_ieg_across_contrasts.svg"); fig.savefig(f"{OUT}/fig_ieg_across_contrasts.png", dpi=200)
cap = ("Figure. The seven rapid primary response genes (EGR1, NR4A1, DUSP1, FOS, MAFF, PPP1R15A, GADD45G; the Tyssowski rapid primary response genes among the 383 recurrent genes) in every contrast.\n"
       "What it is: one row per contrast (surgical holds and autopsy against immediately frozen tissue, the Dachet ex vivo time course, living against postmortem cortex, slow against fast death in GTEx brain regions and other tissues, the Pritzker agonal factor, and the mean schizophrenia case-control t over 22 bulk datasets). Each coloured mark is one gene's effect (t statistic; Dachet rows, log2 change from 0 h, own axis at right); circles are genes lower in schizophrenia lists, triangles genes higher; the black diamond is the mean of the seven; the grey bar is the 95% range of the means of 10,000 random gene sets matched on expression decile.\n"
       "How to read it: positive means higher in the later, dead, slow-death, postmortem, warmer or case condition. A diamond outside the grey bar means the seven as a group sit away from expression-matched genes; marks of one colour that fall on the same side across rows go the same way in those contrasts. The horizontal axis is linear within -3 to 3 and logarithmic beyond.\n"
       "What it shows: see SUMMARY_for_writeup.md; the direction of each gene in each contrast, and whether it matches the gene's schizophrenia direction, is in per_gene_7.tsv. Contrasts mix handling, cohort and tissue; seven genes is a small set.\n")
open(f"{OUT}/fig_ieg_across_contrasts_caption.txt", "w").write(cap)
print("figure rows:", n)
