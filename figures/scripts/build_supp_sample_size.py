"""Supplementary Fig. 1: agreement against sample size, and differences between patients.

One matplotlib figure (figkit), five panels, from figures/data/fig03/.
a  agreement vs subjects per group: LIBD DLPFC split halves vs different cohorts, diagnosis and age
b  the 76 published different-cohort pairs with a results table on both sides: percentile of each
   pair's rho among split halves of its own size
c  Spearman rho of log2 FC of Gemma brain result sets by subject relation
d  case mix: agreement vs difference in cases' median age at death; comparable vs random subsets
e  identification of a donor across datasets vs diagnosis, per relation
f  reliability of the GTEx slow-vs-fast-death axis and the control-pH axis, within one dataset
 and between pH collections (the pH-slope rows of the
   export are a different quantity and are left out), from figures/data/figures_data_quality/panel_c.tsv
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figkit as fk  # noqa: E402

D = fk.FD / "fig03"
SCRIPT = "figures/scripts/build_supp_sample_size.py"
rng = np.random.default_rng(3)


def tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(D / name, sep="\t", comment="#")


def jitter(n: int, w: float = 0.18) -> np.ndarray:
    return rng.uniform(-w, w, n)


# ---------------------------------------------------------------- data
summ = tsv("fig03_a_c_d_agreement_summary.tsv")
units = tsv("fig03_a_c_d_agreement_units.tsv")
pub = tsv("fig03_b_published_vs_split_halves.tsv")
rel = tsv("fig03_e_gemma_rho_by_relation.tsv")
cm_pairs = tsv("fig03_f_casemix_pairs.tsv")
cm_rand = tsv("fig03_f_comparable_vs_random.tsv")
ident = tsv("fig03_g_identification_by_relation.tsv")
REL_TSV = fk.FD / "figures_data_quality" / "panel_c.tsv"
SPLIT_TSV = fk.FD / "figures_data_quality" / "split_half_ph_axis.tsv"
relia = pd.read_csv(REL_TSV, sep="\t", comment="#")
split_ph = pd.read_csv(SPLIT_TSV, sep="\t", comment="#")
GSIZE_TSV = fk.FD / "figures_data_quality" / "split_half_gtex_at_ph_size.tsv"
gtex_size = pd.read_csv(GSIZE_TSV, sep="\t", comment="#")
INPUTS = [D / n for n in [
    "fig03_a_c_d_agreement_summary.tsv", "fig03_a_c_d_agreement_units.tsv",
    "fig03_b_published_vs_split_halves.tsv", "fig03_e_gemma_rho_by_relation.tsv",
    "fig03_f_casemix_pairs.tsv", "fig03_f_comparable_vs_random.tsv",
    "fig03_g_identification_by_relation.tsv"]] + [REL_TSV, SPLIT_TSV, GSIZE_TSV]

fig = fk.figure(222)
gs = GridSpec(3, 6, figure=fig, left=0.085, right=0.985, top=0.935, bottom=0.075,
              hspace=0.8, wspace=2.0, width_ratios=[1, 1, 1, 1.15, 1, 1])
ax_a = fig.add_subplot(gs[0, 0:2])
ax_b = fig.add_subplot(gs[0, 2:4])
ax_c = fig.add_subplot(gs[0, 4:6])
# d's two halves get a wider gap so the right half's row labels clear the left half's points
gsd = gs[1, 0:4].subgridspec(1, 2, wspace=0.95, width_ratios=[1, 0.85])
ax_d1 = fig.add_subplot(gsd[0])
ax_d2 = fig.add_subplot(gsd[1])
ax_e = fig.add_subplot(gs[1, 4:6])
ax_f = fig.add_subplot(gs[2, 2:5])

# ---------------------------------------------------------------- a
ax = ax_a
COL = {"dx": fk.ACCENT, "age": fk.ACCENT_2}
STYLE = {"within cohort (split halves)": "-", "different cohorts": "--"}
for comp, ls in STYLE.items():
    for trait, col in COL.items():
        s = summ[(summ.comparison == comp) & (summ.trait == trait)].sort_values("n_per_group")
        u = units[(units.comparison == comp) & (units.trait == trait)]
        if comp.startswith("within"):
            u = u[u.libd_dlpfc == "yes"]
        ax.scatter(u.n_per_group + jitter(len(u), 1.0), u.rho_t, s=6, color=fk.lighten(col, 0.55),
                   lw=0, zorder=2)
        ax.plot(s.n_per_group, s.median_rho_t, ls, color=col, lw=1.2, zorder=3,
                marker="o", ms=3.5, mec="white", mew=0.4)
        if comp == "different cohorts" and trait == "dx":
            npairs = ", ".join(str(int(v)) for v in s.n_units)
            nper = ", ".join(str(int(v)) for v in s.n_per_group)
ax.set_xticks([10, 20, 30, 50])
ax.set_xlim(6, 54)
ax.set_ylim(-0.05, 0.98)  # room above the lines for the legend
# pair counts as one note in the empty lower right, not on the lines
ax.text(0.98, 0.03, f"different cohorts:\n{npairs} pairs at\n{nper} per group", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=fk.TINY, color=fk.SUBTLE)
ax.set_xlabel("Subjects per group (cases or controls)")
ax.set_ylabel("Agreement between two halves\n(Spearman ρ of limma t)")
ax.set_title("Agreement grows with n within a cohort,\nnot between cohorts", fontsize=8)
handles = [Line2D([], [], color=fk.ACCENT, lw=1.2, label="schizophrenia"),
           Line2D([], [], color=fk.ACCENT_2, lw=1.2, label="age"),
           Line2D([], [], color=fk.TEXT, lw=1.2, ls="-", label="within cohort (LIBD DLPFC)"),
           Line2D([], [], color=fk.TEXT, lw=1.2, ls="--", label="different cohorts")]
ax.legend(handles=handles, loc="upper left", fontsize=fk.SMALL, handlelength=2.2,
          borderaxespad=0.2, labelspacing=0.3, title="point = one dataset or dataset pair", title_fontsize=fk.TINY)

# ---------------------------------------------------------------- b
ax = ax_b
rho = pub[pub.metric == "spearman_rho_log2fc"]
n_rho = len(rho)
pct = rho.libd_percentile.to_numpy()
bins = np.linspace(0, 1, 11)
ax.hist(pct, bins=bins, color=fk.ACCENT, edgecolor="white", lw=0.5, zorder=2)
ax.axhline(n_rho / 10, color=fk.SUBTLE, ls=":", lw=0.8, zorder=3)
ax.text(0.98, n_rho / 10 + 0.4, "expected if like\nsplit halves", ha="right", va="bottom",
        fontsize=fk.TINY, color=fk.SUBTLE, transform=ax.get_yaxis_transform())
med_pct = float(np.median(pct))
med_rho = float(rho.value.median())
ax.axvline(med_pct, color=fk.TEXT, lw=1.0, zorder=4)
ax.text(med_pct + 0.02, ax.get_ylim()[1] * 0.97, f"median {med_pct:.2f}", ha="left", va="top",
        fontsize=fk.SMALL, color=fk.TEXT)
ax.text(0.98, 0.70, f"{n_rho} published pairs,\ndifferent cohorts,\nmedian ρ = {med_rho:.2f}",
        ha="right", va="top", fontsize=fk.SMALL, color=fk.TEXT, transform=ax.transAxes)
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xlabel("Rank among split halves of the pair's\nown size (percentile)")
ax.set_title("Published pairs agree less than split\nhalves of one dataset at their size", fontsize=8)
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
ax.set_ylabel("Published study pairs (count)")

# ---------------------------------------------------------------- c
ax = ax_c
ORDER = [
    ("same subjects, same region, same series", "same subjects,\nsame region"),
    ("same subjects, different region, same series", "same subjects,\ndifferent region"),
    ("partial subject overlap, same series", "partial overlap,\nsame series"),
    ("partial subject overlap, across series", "partial overlap,\nacross series"),
    ("same collection, subject overlap unknown (Stanley, same sub-collection)", "same Stanley\nsub-collection"),
    ("same collection, subject overlap unknown", "same collection,\noverlap unknown"),
    ("subject overlap unknown, bank not stated", "bank not stated,\noverlap unknown"),
    ("subject-separate (bank not stated; fingerprint found no excess)", "no shared subjects,\nbank not stated"),
    ("subject-separate (different banks)", "no shared subjects,\ndifferent banks"),
]
ax.grid(False)
ax.grid(True, axis="x", color=fk.GRID, lw=0.4)
ax.axvline(0, color=fk.SUBTLE, lw=0.6, ls=":")
ylabels = []
for i, (key, lab) in enumerate(ORDER):
    y = len(ORDER) - 1 - i
    v = rel.loc[rel.relation == key, "rho_log2fc"].to_numpy()
    shared = i < 4
    ax.scatter(v, y + jitter(len(v), 0.22), s=5, lw=0,
               color=fk.ACCENT if shared else fk.NEUTRAL, alpha=0.8, zorder=2)
    ax.plot([np.median(v)] * 2, [y - 0.3, y + 0.3], color="black", lw=1.2, zorder=3)
    ylabels.append((y, f"{lab} ({len(v)})"))
ax.set_yticks([y for y, _ in ylabels])
ax.set_yticklabels([l for _, l in ylabels], fontsize=fk.TINY)
ax.set_ylim(-0.6, len(ORDER) - 0.4)
ax.set_xlabel("Agreement between two result\nsets (Spearman ρ of log2 FC)")
ax.set_title("Agreement rises with shared\ndonors (point = one pair)", fontsize=8)
ax.tick_params(axis="y", length=0)

# ---------------------------------------------------------------- d, left
ax = ax_d1
ms = cm_pairs.mount_sinai_pair == "yes"
ax.scatter(cm_pairs.d_median_age[~ms], cm_pairs.rho_adj[~ms], s=7, lw=0, color=fk.ACCENT, alpha=0.7,
           label=f"other pairs ({(~ms).sum()})")
ax.scatter(cm_pairs.d_median_age[ms], cm_pairs.rho_adj[ms], s=9, lw=0, color=fk.ACCENT_3,
           label=f"pairs with PsychAD Mount Sinai ({ms.sum()})")
ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls=":")
ax.set_xlabel("Difference in cases' median age at death\nbetween the two datasets (years)")
ax.set_ylabel("Agreement between the two datasets\n(ρ of t, adjusted for n and technology)")
ax.set_title("Agreement does not fall as case mix\ndiffers (point = one pair of datasets)", fontsize=8)
lo, hi = ax.get_ylim()  # room above and below the points for the note and the legend
ax.set_ylim(lo - 0.22 * (hi - lo), hi + 0.25 * (hi - lo))
ax.legend(loc="lower left", fontsize=fk.TINY, borderaxespad=0.2, handletextpad=0.3)
ax.text(0.02, 0.97, f"{len(cm_pairs)} pairs with no subjects in common, Spearman −0.15, one-sided P = 0.12",
        ha="left", va="top", fontsize=fk.TINY, transform=ax.transAxes)

# ---------------------------------------------------------------- d, right
ax = ax_d2
COND = {"age65": "age ≤ 65", "no_suicide": "no suicides", "age65_no_suicide": "both"}
rows = []
for stat, grp in [("median_rho_t", "all pairs"), ("median_rho_t_affected", "affected pairs")]:
    for cond, clab in COND.items():
        r = cm_rand[(cm_rand.condition == cond) & (cm_rand.statistic == stat)].iloc[0]
        n = int(cm_rand[(cm_rand.condition == cond) &
                        (cm_rand.statistic == stat.replace("median_rho_t", "n_pairs"))].iloc[0].comparable_subsets)
        rows.append((f"{grp}: {clab} ({n})", r))
ax.grid(False)
ax.grid(True, axis="x", color=fk.GRID, lw=0.4)
for i, (lab, r) in enumerate(rows):
    y = len(rows) - 1 - i
    ax.plot([r.random_q05, r.random_q95], [y, y], color=fk.GRID, lw=5, solid_capstyle="butt", zorder=1)
    ax.plot([r.random_median] * 2, [y - 0.22, y + 0.22], color=fk.SUBTLE, lw=1.0, zorder=2)
    ax.scatter([r.full_data], [y], marker="x", s=22, color=fk.TEXT, lw=0.9, zorder=3)
    ax.scatter([r.comparable_subsets], [y], s=22, color=fk.ACCENT, zorder=4, lw=0)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([lab for lab, _ in rows][::-1], fontsize=fk.TINY)
ax.tick_params(axis="y", length=0)
ax.set_ylim(-0.7, len(rows) + 2.6)  # empty band above the rows holds the legend
ax.set_xlabel("Median agreement over dataset pairs\n(Spearman ρ of t)")
ax.set_title("Comparable patient subsets agree\nno better than random subsets", fontsize=8)
handles = [Line2D([], [], marker="o", color=fk.ACCENT, ls="", ms=4.5, label="comparable subset"),
           Line2D([], [], marker="x", color=fk.TEXT, ls="", ms=4.5, mew=0.9, label="full data"),
           Line2D([], [], color=fk.GRID, lw=5, label="random subsets, 5th–95th pct."),
           Line2D([], [], color=fk.SUBTLE, lw=1.0, label="random subsets, median")]
ax.set_xlim(-0.13, 0.145)  # room left of the bars for the legend
ax.legend(handles=handles, loc="upper left", fontsize=fk.TINY, borderaxespad=0.2,
          handlelength=1.6, labelspacing=0.3)

# ---------------------------------------------------------------- e
ax = ax_e
core = ident[(ident.subjects == "all")]
stan = ident[(ident.subjects == "sz_ctrl")]
groups = [("Pittsburgh", core[core.source == "Pittsburgh"]),
          ("LIBD,\nPsychAD", core[core.source.isin(["LIBD regions", "PsychAD x LIBD"])]),
          ("Stanley", stan)]
xt, xl = [], []
for i, (lab, g) in enumerate(groups):
    x0 = i * 1.0
    idv = g.auroc_identity.to_numpy()
    dxv = g.auroc_dx.dropna().to_numpy()
    ax.scatter(x0 - 0.18 + jitter(len(idv), 0.1), idv, s=9, lw=0, color=fk.ACCENT, zorder=2)
    ax.scatter(x0 + 0.18 + jitter(len(dxv), 0.1), dxv, s=9, lw=0, color=fk.NEUTRAL, zorder=2)
    for xx, v in [(x0 - 0.18, idv), (x0 + 0.18, dxv)]:
        if len(v):
            ax.plot([xx - 0.14, xx + 0.14], [np.median(v)] * 2, color="black", lw=1.2, zorder=3)
    xt.append(x0)
    xl.append(f"{lab}\n({len(g)})")
ax.axhline(0.5, color=fk.SUBTLE, lw=0.6, ls=":")
ax.text(-0.45, 0.505, "chance", ha="left", va="bottom", fontsize=fk.TINY, color=fk.SUBTLE)
ax.set_xticks(xt)
ax.set_xticklabels(xl, fontsize=fk.SMALL)
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(0.35, 1.02)
ax.set_ylabel("Recovery in the other dataset\n(mean AUROC; 0.5 = chance)")
ax.set_xlabel("Pairs sharing donors (source)")
ax.set_title("A donor is found across\ndatasets, a diagnosis is\nnot (point = one pair)", fontsize=8)
handles = [Line2D([], [], marker="o", color=fk.ACCENT, ls="", ms=4, label="same donor found"),
           Line2D([], [], marker="o", color=fk.NEUTRAL, ls="", ms=4, label="diagnosis found")]
ax.legend(handles=handles, loc="lower left", fontsize=fk.TINY, borderaxespad=0.2, handletextpad=0.3)

# ---------------------------------------------------------------- f
ax = ax_f
# group each row by its "what" text; the measure column mixes group headings into some labels
HEADS = ["GTEx slow-vs-fast-death axis (one dataset, donors split)",
         "Control-pH axis (one collection, donors split)",
         "Control-pH axis (between collections)"]


def rel_group(what: str) -> int:
    if "GTEx" in what or "three-region" in what:
        return 0
    return 2 if what.startswith("pH axis") else -1


GROUPS = [[h, []] for h in HEADS]
for _, r in relia.iterrows():
    parts = [x.strip() for x in str(r.measure).split("   ") if x.strip()]
    lab = " ".join(parts) if parts[-1].startswith("(") else parts[-1]
    lab = lab.split(", ", 1)[-1] if lab.startswith("Control-pH axis") else lab
    g = rel_group(str(r.what))
    if g >= 0:
        GROUPS[g][1].append((lab, r))
# GTEx subsampled to the pH collections' sizes (mean over 20 subsamples), after the two-halves row
for i, (_, r) in enumerate(gtex_size[gtex_size.measure.str.contains("mean over 20 subsamples")].iterrows()):
    n = re.search(r"(\d+) donors", r.measure).group(1)
    GROUPS[0][1].insert(1 + i, (f"subsampled to {n} donors", r))
GROUPS[0][1][0] = ("all 285 donors", GROUPS[0][1][0][1])
for _, r in split_ph.iterrows():
    m = re.search(r"halves of the (\S+) controls \(n=(\d+)\)", r.measure)
    GROUPS[1][1].append((f"{m.group(1).replace('_', ' ')}, {m.group(2)} controls", r))
GCOL = [fk.ACCENT_3, fk.ACCENT, fk.ACCENT_2]
ax.grid(False)
ax.grid(True, axis="x", color=fk.GRID, lw=0.4)
y = 0
yt, yl = [], []
for gi, (head, rows_) in enumerate(GROUPS):
    col = GCOL[gi % len(GCOL)]
    ax.text(-0.02, y, head, ha="right", va="center", fontsize=fk.SMALL, color=fk.TEXT,
            fontweight="bold", transform=ax.get_yaxis_transform())
    y -= 1
    for lab, r in rows_:
        hollow = "overstated" in lab or "alone" in lab
        if pd.notna(r.lo):
            ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.0, zorder=2)
        ax.scatter([r.reliability], [y], s=20, zorder=3, color="white" if hollow else col,
                   edgecolor=col, lw=0.9)
        if pd.notna(r.lo):  # label right of the interval so it touches neither line nor point
            ax.text(r.hi + 0.02, y, f"{r.reliability:.2f}", ha="left", va="center",
                    fontsize=fk.TINY, color=fk.TEXT)
        else:
            ax.text(r.reliability + 0.02, y, f"{r.reliability:.2f}", ha="left", va="center",
                    fontsize=fk.TINY, color=fk.TEXT)
        yt.append(y)
        yl.append(lab.replace("(overstated)", "").strip())
        y -= 1
    y -= 0.3
ax.set_yticks(yt)
ax.set_yticklabels(yl, fontsize=fk.TINY)
ax.tick_params(axis="y", length=0)
ax.set_ylim(y + 0.6, 0.6)
ax.set_xlim(-0.1, 1.08)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xlabel("Reliability (1 = fully repeatable)")
ax.set_title("Reliability of the slow-death and pH axes\n"
             "(point = one estimate; hollow = overstated or a single collection)", fontsize=8)

# ---------------------------------------------------------------- letters, save
for a, s in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d1, "d"), (ax_e, "e")]:
    fk.letter(a, s)
fk.letter(ax_f, "f", dx=-0.3)

fk.save(fig, "supp_sample_size", SCRIPT, INPUTS,
        notes=f"b: {n_rho} pairs, median percentile {med_pct:.2f}, median rho {med_rho:.3f}; "
              "e: Pittsburgh/LIBD/PsychAD rows use subjects='all', Stanley rows subjects='sz_ctrl'.")
print("supp_sample_size built", n_rho, med_pct, med_rho)
