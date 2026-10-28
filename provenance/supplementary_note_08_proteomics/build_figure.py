#!/usr/bin/env python3
"""Step 6 (2026-10-14): figure and figure data. Reads the results TSVs in the output folder. Run: PYTHONDONTWRITEBYTECODE=1 python3 -I build_figure.py"""
import os, re, hashlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.environ["OUTDIR"]
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial"], "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#6b7280", "axes.labelcolor": "#1f2937", "xtick.color": "#6b7280", "ytick.color": "#6b7280", "axes.titlelocation": "left", "axes.titleweight": "normal",
                     "axes.titlesize": 8.5, "axes.grid": False, "figure.facecolor": "white"})
BLUE, AMBER, GREY, TEXT, SUB, EM = "#2563eb", "#f59e0b", "#e5e7eb", "#1f2937", "#6b7280", "#10b981"
def rd(f): return pd.read_csv(f"{OUT}/{f}", sep="\t", comment="#")
U = rd("overlap_209_by_unit.tsv"); P = rd("overlap_209_pooled_and_sensitivity.tsv"); A = rd("axis_effect_correlation_by_unit.tsv"); AP = rd("axis_effect_correlation_pooled.tsv")
L = rd("axis_direction_lists_pooled.tsv"); CMP = rd("axis_comparison_with_rna_pooled.tsv")
data = []
fig = plt.figure(figsize=(7.6, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1.15, 1], hspace=0.45, wspace=0.95)
# ---------- a
ax = fig.add_subplot(gs[0, 0])
order = ["Testa 2026", "Koopmans 2026", "MacDonald 2020", "Focking 2015 (PSD)", "Pinacho 2016", "Vidal-Domenech 2020", "Nucifora 2019", "Alnafisah 2022", "Krishnamurthy 2013", "Vera-Montecinos 2021"]
rows = []
for u in order:
    r = U[U.unit == u].iloc[0]; lo, hi = [int(x) for x in r.null95_matched.split(" to ")]; E = r.expected_matched
    rows.append((u, r.hits_in_209, E, lo, hi, int(r.list_proteins_in_universe), r.p_matched))
for lab, key in [("9 units pooled", "pooled, 9 units"), ("without Testa 2026", "9 units without Testa 2026"), ("without Pittsburgh group", "9 units without the Pittsburgh group")]:
    r = P[P.unit.str.startswith(key)].iloc[0]; lo, hi = [int(x) for x in r.null95_matched.split(" to ")]
    rows.append((lab, r.hits_in_209, r.expected_matched, lo, hi, np.nan, r.p_matched))
ys = np.arange(len(rows))[::-1].astype(float); ys[-3:] -= 0.6
for y, (u, o, E, lo, hi, n, p) in zip(ys, rows):
    ax.axhline(y, color=GREY, lw=0.5, zorder=0)
    ax.plot([lo / E, hi / E], [y, y], color="#9ca3af", lw=4.5, solid_capstyle="butt", zorder=1)
    ax.plot(o / E, y, "o", color=BLUE, ms=4.5, zorder=3)
    ax.text(4.0, y, f"{int(o)} (exp. {E:.1f})", va="center", ha="left", fontsize=6.5, color=TEXT)
    data.append({"panel": "a", "row": u, "series": "hits among the 209 / matched-null expected", "value": round(o / E, 3), "lo": round(lo / E, 3), "hi": round(hi / E, 3), "n": int(o), "note": f"expected {E:.2f}; null 95% {lo} to {hi}; P matched {p:.2g}"})
ax.axvline(1, color=SUB, lw=0.7, ls=":")
ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=7); ax.set_xlim(0, 5.4); ax.set_xticks([0, 1, 2, 3]); ax.set_xlabel("Hits among the 209 / expected\nfrom matched random proteins", fontsize=7.5)
ax.set_title("a  Hits among the recurrent genes", fontsize=8.5, loc="left", pad=10)
ax.text(4.0, ys[0] + 0.6, "hits (expected)", fontsize=6.3, color=SUB, ha="left")
ax.tick_params(axis="y", length=0)
# ---------- b
ax = fig.add_subplot(gs[0, 1])
units = ["Testa 2026", "Koopmans 2026", "MacDonald 2020 homogenate", "MacDonald 2020 synaptosome", "Nucifora 2019", "Focking 2015 PSD", "Pinacho 2016", "Vidal-Domenech 2020", "Alnafisah 2022", "Vera-Montecinos 2021"]
yb = np.arange(len(units))[::-1].astype(float) + 1.6
for y, u in zip(yb, units):
    ax.axhline(y, color=GREY, lw=0.5, zorder=0)
    for dy, axn, col in [(0.16, "pH", BLUE), (-0.16, "agonal", AMBER)]:
        r = A[(A.unit == u) & (A.axis == axn) & (A.set == "all tested proteins")].iloc[0]
        ax.plot([r.ci_lo, r.ci_hi], [y + dy, y + dy], color=col, lw=1.2, zorder=2); ax.plot(r.rho, y + dy, "o", color=col, ms=3.5, zorder=3)
        data.append({"panel": "b", "row": u, "series": axn + " axis", "value": round(r.rho, 3), "lo": round(r.ci_lo, 3), "hi": round(r.ci_hi, 3), "n": int(r.n), "note": "Spearman, effect vs axis, all tested proteins"})
y0 = 0.4
ax.axhline(y0 + 0.3, color=SUB, lw=0.5, zorder=0)
for dy, axn, col in [(0.16, "pH", BLUE), (-0.16, "agonal", AMBER)]:
    r = AP[(AP.axis == axn) & (AP.units == "all 10 full-table units")].iloc[0]
    ax.plot([r.ci_lo, r.ci_hi], [y0 + dy, y0 + dy], color=col, lw=1.6); ax.plot(r.pooled_rho_mean_over_donor_groups, y0 + dy, "D", color=col, ms=4)
    data.append({"panel": "b", "row": "pooled over 6 donor groups", "series": axn + " axis", "value": round(r.pooled_rho_mean_over_donor_groups, 3), "lo": round(r.ci_lo, 3), "hi": round(r.ci_hi, 3), "n": 10, "note": "mean over donor groups; interval bootstraps proteins only"})
    rr = CMP[(CMP.axis == axn) & CMP.quantity.str.startswith("CommonMind")].iloc[0].pooled_rho_mean_over_donor_groups
    ax.axvline(rr, color=col, lw=0.8, ls="--", alpha=0.8)
    data.append({"panel": "b", "row": "RNA reference (CommonMind t, same proteins), pooled", "series": axn + " axis", "value": rr, "lo": np.nan, "hi": np.nan, "n": 10, "note": "dashed line"})
ax.axvline(0, color=SUB, lw=0.7, ls=":")
ax.set_yticks(list(yb) + [y0]); ax.set_yticklabels([u.replace(" 2020", " 2020").replace("Focking 2015 PSD", "Focking 2015 (PSD)") for u in units] + ["Pooled"], fontsize=7); ax.tick_params(axis="y", length=0)
ax.set_xlim(-0.35, 0.5); ax.set_xlabel("Spearman rho, protein effect vs axis", fontsize=7.5, labelpad=2)
ax.set_title("b  Effect sizes along the two axes", fontsize=8.5, loc="left", pad=10)
ax.plot([], [], "o", color=BLUE, ms=3.5, label="pH axis"); ax.plot([], [], "o", color=AMBER, ms=3.5, label="agonal axis"); ax.plot([], [], ls="--", color=SUB, lw=0.8, label="RNA, same proteins")
ax.legend(frameon=False, fontsize=6.5, loc="upper left", ncol=3, handlelength=1.2, bbox_to_anchor=(-0.2, -0.2), columnspacing=1.0)
# ---------- c
ax = fig.add_subplot(gs[1, :])
sets = [("the 31 units of the side study (13 donor groups, as Extended Fig. 4)", "31 units, 13 donor groups"), ("all units (31 + 5 added)", "36 units, 17 donor groups"), ("all units without the Pittsburgh group", "36 units without Pittsburgh group")]
xs = []; x = 0
for sname, lab in sets:
    for axn, col, key in [("pH", BLUE, "pH_matched"), ("agonal", AMBER, "agonal")]:
        r = L[(L.set == sname) & (L.axis == key) & (L.top_n == 200)].iloc[0]
        ax.bar(x, r.share, width=0.8, color=col, zorder=2)
        ax.plot([x - 0.4, x + 0.4], [r.expected_share] * 2, color=TEXT, lw=1.3, zorder=3)
        ax.text(x, max(r.share, r.expected_share) + 0.03, f"{int(r.concordant)}/{int(r.hits)}", ha="center", va="bottom", fontsize=6.8, color=TEXT)
        data.append({"panel": "c", "row": lab, "series": axn + " axis, top 200", "value": r.share, "lo": r.expected_share, "hi": np.nan, "n": int(r.hits), "note": f"concordant {int(r.concordant)}; value = share, lo column = share expected from label shuffle; excess {r.excess}; shuffle P {r.p_one_sided_shuffle:.2g}; jackknife over groups {r.excess_ci_lo_jackknife_groups} to {r.excess_ci_hi_jackknife_groups}"})
        x += 1
    xs.append(x - 1.5); x += 0.8
ax.set_xticks(xs); ax.set_xticklabels([s[1] for s in sets], fontsize=7)
ax.set_ylim(0, 1.08); ax.set_xlim(-0.8, x - 0.5); ax.set_ylabel("Share of hits in the predicted direction", fontsize=7.5)
ax.yaxis.grid(True, color=GREY, lw=0.5, zorder=0); ax.set_axisbelow(True)
ax.set_title("c  Direction of protein changes among the 200 genes at the top of each axis", fontsize=8.5, loc="left")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[Patch(color=BLUE, label="pH axis (collection sharing the bank left out)"), Patch(color=AMBER, label="agonal axis"), Line2D([], [], color=TEXT, lw=1.3, label="expected from shuffling labels")], frameon=False, fontsize=6.8, loc="upper right", ncol=3, bbox_to_anchor=(1.0, 1.0))
fig.text(0.01, 0.005, "Source: protein tables of 10 brain studies, lists of 36 units; axes from control-brain pH and GTEx slow-minus-fast death; 209 recurrent genes.", fontsize=6.3, color=SUB)
fig.subplots_adjust(left=0.17, right=0.98, top=0.94, bottom=0.09)
for a_ in fig.axes:
    for art in list(a_.patches) + list(a_.lines) + list(a_.texts) + list(a_.collections): art.set_clip_on(False)
svg = f"{OUT}/fig_proteomics_axes.svg"; fig.savefig(svg, format="svg")
t_ = open(svg).read(); t_ = re.sub(r'<clipPath.*?</clipPath>', '', t_, flags=re.S); t_ = re.sub(r' clip-path="url\(#[^)]*\)"', '', t_); open(svg, "w").write(t_); fig.savefig(f"{OUT}/fig_proteomics_axes.png", dpi=200)
D = pd.DataFrame(data)
with open(f"{OUT}/fig_proteomics_axes_data.tsv", "w") as f:
    f.write("# figure: proteomics side-study figure (panels a, b, c); script: scripts/build_figure.py; built 2026-10-14\n")
    f.write("# sources: overlap_209_by_unit.tsv, overlap_209_pooled_and_sensitivity.tsv, axis_effect_correlation_by_unit.tsv, axis_effect_correlation_pooled.tsv, axis_comparison_with_rna_pooled.tsv, axis_direction_lists_pooled.tsv (sha256 of each in README.md)\n")
    f.write("# columns: panel; row; series; value (a: hits/expected; b: Spearman rho; c: share concordant); lo, hi (a: matched-null 95% range divided by expected; b: 95% bootstrap interval; c: lo = expected share); n (a: hits; b: proteins; c: hits in the top 200); note\n")
    D.to_csv(f, sep="\t", index=False)
t = open(svg).read(); fams = sorted(set(re.findall(r"font-family:\s*([^;\"]+)", t))); print("font families in svg:", fams, "| DejaVu present:", "DejaVu" in t, "| clipPath:", "clipPath" in t)
