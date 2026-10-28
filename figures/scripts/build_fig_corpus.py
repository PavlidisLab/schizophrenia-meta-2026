#!/usr/bin/env python3
"""Figure 1 -- the corpus.

Data: figures/data/fig01/ (export_fig01_fig11.py, 2026-10-08). Panel a is the pipeline schematic built by build_fig1_pipeline.py.

Panel letters. The panel_* function
names below are the old letters.
Figure 1 keeps a and b; the old panels c to g are Extended Fig. 2 a to e and the old h is Extended Fig. 2f (supp_corpus).
a  pipeline schematic, with counts from figures/data/funnel/funnel_table.tsv (screened = PRISMA SC01, corpus 1, the primary brain studies row 11)
b  lineage: studies by year and source lane, with reuse links (panel_h)
c  corpus studies per year, by tissue (panel_b)
d  technology by era (panel_c)
e  case + control N per study by technology (panel_d)
f  reported vs distinct subjects in the primary brain studies (panel_e)
g  covariates adjusted for, matched or tested, the primary studies, by era (panel_f)
h  the primary studies: pooled by a meta-analysis, reanalysed only, never reused, by bank (panel_g)
"""
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import figkit as fk

D = fk.FD / "fig01"
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
# the analysed base: 71 studies with the authors' list + 10 with a list derived by us (80); panels e to h are drawn on it
N_PRIMARY = len(pd.read_csv(D / "fig01_d_neglected_datasets_80.tsv", sep="\t", comment="#"))
rng = np.random.default_rng(1)

TISSUE = [("brain", "brain", fk.ACCENT), ("blood", "blood", fk.ACCENT_4),
          ("cells", "iPSC / LCL", fk.ACCENT_2), ("other", "other / mixed", fk.NEUTRAL)]

TISSUE_OF = {"brain": "brain", "blood": "blood", "ipsc": "cells", "lcl": "cells", "mixed": "other", "other": "other"}

# technology has its own hues, apart from the tissue colours (panels c, e) and the study-kind colours (panel b)
TECH_CYAN, TECH_PINK, TECH_BROWN = "#0891b2", "#db2777", "#92400e"
TECH = [("array", "microarray", fk.lighten(TECH_CYAN, 0.55)), ("rnaseq", "RNA-seq", TECH_CYAN),
        ("single_cell", "single-cell", TECH_PINK), ("proteomics", "proteomics", TECH_BROWN),
        ("other", "other", fk.NEUTRAL)]


def tech_group(t):
    if str(t).startswith("array"):
        return "array"
    return t if t in ("rnaseq", "single_cell", "proteomics") else "other"


BANK_SHORT = {
    "Stanley Medical Research Institute (SMRI): Array Collection and Neuropathology Consortium": "Stanley",
    "CommonMind Consortium (CMC), PsychENCODE Consortium (PEC), NIMH Human Brain Collection Core (HBCC)": "CMC / HBCC",
    "New South Wales Brain Tissue Resource Centre (NSW TRC), Victorian Brain Bank Network (VBBN)": "NSW / Victoria",
    "University of Pittsburgh Brain Tissue Donation Program (Conte Center, CCNMD)": "Pittsburgh",
    "Lieber Institute for Brain Development (LIBD), BrainSeq Consortium, NIMH Clinical Brain Disorders Branch (CBDB)": "LIBD",
    "Mount Sinai / JJ Peters VA Medical Center NIH Brain and Tissue Repository (MSSM)": "Mount Sinai",
    "Harvard Brain Tissue Resource Center (HBTRC), McLean Hospital": "Harvard / McLean",
    "Pritzker Neuropsychiatric Disorders Research Consortium, UC Irvine Brain Donor Program": "Pritzker",
    "Other brain samples": "other brain",
    "iPSC-derived cells and organoids": "iPSC / organoids",
    "Blood, LCL, fibroblast, other": "blood, LCL, other",
    "Meta-analyses and combined analyses": "meta-analyses",
}


# plain labels
PHASES = [
    ("Search", ["PubMed, EMBASE, WoS", "citing and cited papers",
                           "single-gene and single-cell studies", "screening (PRISMA)"]),
    ("Annotation", ["technology and tissue", "brain bank and cohort", "kind of result reported",
                    "data reuse between studies"]),
    ("Data retrieval", ["full text and supplements", "gene lists and directions",
                        "raw data reanalysed", "validation and provenance"]),
    ("Analyses", ["recurrence, pairwise agreement", "sample size, donors, case mix",
                  "pH, agonal state, covariates", "positive and negative controls"]),
]
PHASE_COLOURS = [fk.ACCENT, fk.ACCENT_3, "#7c3aed", fk.ACCENT_2]


FUNNEL = fk.FD / "funnel" / "funnel_table.tsv"
PRISMA = fk.FD / "prisma_flow" / "prisma_flow.tsv"


def phase_counts():
    f = pd.read_csv(FUNNEL, sep="\t", comment="#", dtype=str).set_index("step")
    n = lambda k: int(str(f.loc[k, "n"]).replace(",", ""))
    # records screened on title and abstract: the PRISMA box SC01, so Fig. 1a and the PRISMA figure agree
    # (funnel row T2 counts distinct PMIDs held in candidate files plus lineage-only PMIDs, a different quantity)
    pf = pd.read_csv(PRISMA, sep="\t", comment="#", dtype=str).set_index("box_id")
    sc01 = int(str(pf.loc["SC01", "n"]).replace(",", ""))
    # 2026-10-12: funnel_table.tsv not regenerated for the integrated corpus; corpus and primary counts come from
    # the regenerated fig01_b corpus table and fig02 study annotation
    n_corpus = len(pd.read_csv(fk.FD / "fig01" / "fig01_b_corpus_studies.tsv", sep="\t", comment="#"))
    return [f"{sc01:,} screened", f"{n_corpus:,} studies", None, f"{N_PRIMARY:,} studies"]


def panel_a(ax):
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    counts = phase_counts()
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    w, gap = 0.86, 0.14
    for i, ((title, items), col) in enumerate(zip(PHASES, PHASE_COLOURS)):
        x = i + gap / 2
        ax.add_patch(FancyBboxPatch((x, 0.02), w, 0.96, boxstyle="round,pad=0,rounding_size=0.04",
                                    facecolor=fk.lighten(col, 0.88), edgecolor=col, lw=0.8))
        head = f"{i + 1}. {title}" + (f" ({counts[i]})" if counts[i] else "")
        ax.text(x + 0.04, 0.87, head, fontsize=7.4, fontweight="bold", color=fk.TEXT, va="center")
        for j, it in enumerate(items):
            ax.text(x + 0.04, 0.67 - j * 0.17, "\u2022 " + it, fontsize=6.2, color=fk.TEXT, va="center")
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + w + 0.005, 0.5), (x + w + gap - 0.005, 0.5), arrowstyle="-|>",
                                         mutation_scale=6, color=fk.SUBTLE, lw=0.8))


def panel_b(ax, c):
    years = np.arange(c["year"].min(), c["year"].max() + 1)
    bottom = np.zeros(len(years))
    for key, lab, col in TISSUE:
        n = c[c["tgroup"] == key].groupby("year").size().reindex(years, fill_value=0).to_numpy()
        ax.bar(years, n, 0.85, bottom=bottom, color=col, label=lab, lw=0)
        bottom += n
    ax.axvline(2015.5, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_title("Studies per year", pad=4)
    ax.set_xlabel("Year of publication")
    ax.set_ylabel("Studies (count; one bar = one year)")
    ax.set_xticks([1990, 2000, 2010, 2020])
    ax.legend(loc="upper left", handlelength=0.8)
    ax.set_xlim(1983.5, 2027)


def panel_c(ax, c):
    eras = [("old", "<2016"), ("new", "2016+")]
    for i, (e, _) in enumerate(eras):
        sub = c[c["era_2016"] == e]
        share = sub["tg"].value_counts(normalize=True)
        b = 0
        for key, lab, col in TECH:
            v = share.get(key, 0)
            ax.bar(i, v, 0.7, bottom=b, color=col, lw=0)
            b += v
    ax.set_xticks(range(2))
    ax.set_xticklabels([f"{l}\n({(c['era_2016'] == e).sum()})" for e, l in eras])
    ax.set_title("Technology by era", pad=4)
    ax.set_ylabel("Studies (share of each era)")
    ax.set_xlabel("Era of publication (studies)")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    share = c[c["era_2016"] == "new"]["tg"].value_counts(normalize=True)
    b, ys = 0, []
    for key, lab, col in TECH:
        v = share.get(key, 0)
        ys.append([b + v / 2, lab, col])
        b += v
    for i in range(1, len(ys)):  # keep labels at least 9% apart
        ys[i][0] = max(ys[i][0], ys[i - 1][0] + 0.09)
    for y, lab, col in ys:
        ax.text(1.38, y, lab, fontsize=fk.SMALL, color=fk.TEXT, va="center")
    ax.set_xlim(-0.45, 2.6)
    ax.spines["bottom"].set_bounds(-0.45, 1.45)


def panel_d(ax, c):
    sub = c[c["relation"].isin(["primary", "single-cell / single-nucleus"]) & c["n_total"].notna()]
    groups = [("array", "array"), ("rnaseq", "RNA-\nseq"), ("single_cell", "single\ncell")]
    strips = [("brain", fk.ACCENT), ("other tissues", fk.NEUTRAL)]  # brain in the colour of panel c, everything else grey
    step = 0.3
    for i, (g, _) in enumerate(groups):
        for k, (lab, col) in enumerate(strips):
            is_brain = (sub["tgroup"] == "brain") if lab == "brain" else (sub["tgroup"] != "brain")
            yt = sub.loc[(sub["tg"] == g) & is_brain, "n_total"].to_numpy(float)
            xc = i + (k - 0.5) * step
            ax.scatter(xc + rng.uniform(-0.11, 0.11, len(yt)), yt, s=2.5, color=col, lw=0)
            if len(yt) >= 3:
                ax.plot([xc - 0.13, xc + 0.13], [np.median(yt)] * 2, color=fk.TEXT, lw=1.2)
    ax.set_yscale("log")
    ax.set_ylim(top=8000)  # room for the legend above the largest studies
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{l}\n({(sub['tg'] == g).sum()})" for g, l in groups])
    ax.set_title("Sample size by technology and tissue", pad=4)
    ax.set_ylabel("Cases + controls per study (log scale)")
    ax.set_xlabel("Technology (studies)")
    ax.text(0.03, 0.03, "one point = one study", transform=ax.transAxes, fontsize=fk.SMALL, color=fk.SUBTLE)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", color=col, lw=0, markersize=3, label=lab) for lab, col in strips],
              loc="upper right", frameon=False, handletextpad=0.1, borderaxespad=0.2, fontsize=fk.SMALL)


def panel_e(ax):
    d = rd("fig01_b_distinct_subjects_80.tsv")
    d = d[(d["set"] == "published brain") & (d["level"] == "total") & (d["tier"] == "authors+derived")].set_index("measure")
    rows = [("all", "all"), ("SZ", "SZ"), ("CTRL", "controls")]
    x = np.arange(len(rows))
    ax.bar(x - 0.18, [d.loc[m, "naive"] for m, _ in rows], 0.34, color=fk.NEUTRAL, label="sum of reported N")
    lo = np.array([d.loc[m, "lower"] for m, _ in rows], float)
    hi = np.array([d.loc[m, "upper"] for m, _ in rows], float)
    ax.bar(x + 0.18, hi, 0.34, color=fk.lighten(fk.ACCENT, 0.6), label="distinct, upper bound")
    ax.bar(x + 0.18, lo, 0.34, color=fk.ACCENT, label="distinct, lower bound")
    for xi, m in zip(x, [m for m, _ in rows]):
        ax.text(xi - 0.18, d.loc[m, "naive"] + 80, f"{int(d.loc[m, 'naive']):,}", ha="center", va="bottom",
                fontsize=fk.SMALL, color=fk.TEXT, rotation=90)
        ax.text(xi + 0.18, d.loc[m, "upper"] + 80, f"{int(d.loc[m, 'lower']):,}–{int(d.loc[m, 'upper']):,}",
                ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in rows])
    ax.set_title("Reported versus distinct subjects", pad=4)
    ax.set_ylabel(f"Subjects in the {N_PRIMARY} studies (count)")
    ax.set_xlabel("Subject group")
    ax.set_ylim(0, d["naive"].max() * 2.0)  # room for the legend above the rotated labels
    ax.legend(loc="upper right", handlelength=0.8, fontsize=fk.SMALL)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:,.0f}")


COVS = [("age", "age"), ("sex", "sex"), ("pmi", "PMI"), ("ph", "pH"), ("rin", "RIN"),
        ("agonal", "agonal state"), ("cause_of_death", "cause of death"), ("medication", "medication"),
        ("smoking", "smoking"), ("batch", "batch"), ("cell_composition", "cell composition")]
ACTED = {"adjusted", "matched_or_compared", "sensitivity",
         "adjusted covariate", "matched", "sensitivity only", "not covariate in primary; sensitivity with pH"}


def _death_calls(f, col):
    d = pd.read_csv(fk.FD / "fig01" / f, sep="\t", comment="#", dtype=str)
    return d.set_index("pmid")[col]


DEATH_CALLS = {"agonal": _death_calls("calls_agonal_state_80.tsv", "call_a"),
               "cause_of_death": _death_calls("calls_cause_manner_of_death_80.tsv", "call_b")}


def panel_f(ax):
    h = rd("fig01_b_reporting_completeness_80.tsv")  # all 80 analysed studies (the 7 derived-list studies without a list included)
    y = np.arange(len(COVS))
    for off, (e, col, lab) in zip([-0.2, 0.2], [("old", fk.lighten(fk.ACCENT, 0.5), "before 2016"),
                                                 ("new", fk.ACCENT, "2016 on")]):
        sub = h[h["era_2016"] == e]
        share = []
        for c, _ in COVS:
            # all covariates: full-text calls for the 80 (_v2 columns; pH from the 2026-10-09 reading, the other
            # 12 from the two-reader re-call of 2026-10-10, figures/data/fig01/fig01_f_covariate_calls_80.tsv)
            if c in DEATH_CALLS:  # agonal state and cause of death are separate calls (RA, 2026-10-18; two blind readers)
                v = DEATH_CALLS[c].reindex(sub["pmid"].astype(str))
                share.append(v.isin(ACTED | {"exclusion_rule"}).mean())
                continue
            v = sub[f"cov_{c}_v2"]
            v = v[v != "not_applicable"]
            share.append(v.isin(ACTED).mean() if len(v) else np.nan)
        ax.barh(y + off, share, 0.38, color=col, label=f"{lab} ({len(sub)})")
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l in COVS], fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.set_title("Covariates handled, by era", pad=14)
    ax.set_xlabel("Studies that adjusted for, matched on, tested\nor (agonal state, cause of death) excluded on\nthe covariate (share of era)")
    ax.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=2, handlelength=0.8, columnspacing=0.8,
              fontsize=fk.SMALL, borderaxespad=0.2)
    ax.tick_params(axis="y", length=0)


def panel_g(ax):
    """The primary studies by cohort cluster (cohort_cluster column of fig01_d), coloured by earlier reuse (fig01_d).
    A study that pools several banks is counted in the cluster of each of its banks."""
    d = rd("fig01_d_neglected_datasets_80.tsv")
    names = {"cmc+mssm+pitt+libd": "CMC, Mount Sinai, Pitt,\nLIBD (merged)", "stanley": "Stanley",
             "aus": "NSW / Victoria", "harvard": "Harvard / McLean", "pritzker": "Pritzker"}
    d["status"] = np.where(d["ever_pooled"] == "yes", "pooled",
                           np.where(d["ever_reused_any"] == "yes", "reanalysed only", "never reused"))
    rows = []
    for _, r in d.iterrows():
        for c in str(r["cohort_cluster"]).split("; "):
            rows.append((r["label"], c, r["status"]))
    x = pd.DataFrame(rows, columns=["label", "cluster", "status"])
    n_one = int((~x["cluster"].isin(names)).sum())  # each such membership is a cluster of one study
    ONE = f"one study each\n({n_one} clusters)"
    x["grp"] = [names.get(c, ONE) for c in x["cluster"]]
    order = [ONE] + x[x["grp"] != ONE].groupby("grp").size().sort_values().index.tolist()
    st = [("pooled", "pooled", fk.ACCENT), ("reanalysed only", "reanalysed only", fk.ACCENT_3),
          ("never reused", "never reused", fk.NEUTRAL)]
    left = np.zeros(len(order))
    for s_, lab, col in st:
        n = x[x["status"] == s_].groupby("grp").size().reindex(order, fill_value=0).to_numpy()
        n_studies = int((d["status"] == s_).sum())
        ax.barh(range(len(order)), n, 0.7, left=left, color=col, label=f"{lab} ({n_studies})")
        left += n
    for yi, tot in enumerate(left):
        ax.text(tot + 0.3, yi, str(int(tot)), va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=fk.SMALL)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    n_g = N_PRIMARY
    ax.set_title(f"The {n_g} by cohort cluster", pad=4)
    ax.set_xlabel(f"Studies (of the {n_g}; a study that pools banks is counted in each of its clusters)")
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.42), handlelength=0.8, fontsize=fk.SMALL)
    ax.set_xlim(0, float(left.max()) + 2.5)  # room for the count printed after the longest bar
    ax.tick_params(axis="y", length=0)
    assert len(d) == N_PRIMARY

REL = [("primary", "primary", fk.ACCENT), ("single-cell / single-nucleus", "single-cell", fk.ACCENT_3),
       ("reanalysis of public data", "reanalysis", fk.ACCENT_2),
       ("meta-analysis / combined", "meta-analysis", fk.ACCENT_4)]
# meta-analyses lane in the middle, on a light tint
LANES = ["Stanley", "CMC / HBCC", "Pittsburgh", "LIBD", "Mount Sinai", "NSW / Victoria", "meta-analyses",
         "Harvard / McLean", "Pritzker", "other brain", "iPSC / organoids", "blood, LCL, other"]


def panel_h(ax):
    n = rd("fig01_c_lineage_nodes_80.tsv")
    e = rd("fig01_c_lineage_edges.tsv")
    n = n[(n["in_corpus"] == "yes") & n["year"].notna()].copy()
    n["lane"] = n["bank_label"].map(BANK_SHORT)
    n = n[n["lane"].isin(LANES)]
    n["y"] = n["lane"].map({l: i for i, l in enumerate(LANES)}) + rng.uniform(-0.3, 0.3, len(n))
    n["x"] = n["year"] + rng.uniform(-0.35, 0.35, len(n))
    pos = n.set_index("id")[["x", "y"]]
    mi = LANES.index("meta-analyses")
    ax.axhspan(mi - 0.45, mi + 0.45, color=fk.ACCENT_4, alpha=0.07, lw=0, zorder=0)
    k = 0
    for _, r in e.iterrows():
        if r["from"] in pos.index and r["to"] in pos.index:
            a, b = pos.loc[r["from"]], pos.loc[r["to"]]
            ax.plot([a["x"], b["x"]], [a["y"], b["y"]], color=fk.GRID, lw=0.3, zorder=1)
            k += 1
    for rel, lab, col in REL + [(None, "other", fk.NEUTRAL)]:
        sub = n[n["relation"] == rel] if rel else n[~n["relation"].isin([r for r, _, _ in REL])]
        ax.scatter(sub["x"], sub["y"], s=10, color=col, lw=0, zorder=2, label=f"{lab} ({len(sub)})")
    pop = n[n["in_primary_population"] == "yes"]
    ax.scatter(pop["x"], pop["y"], s=20, facecolor="none", edgecolor=fk.TEXT, lw=0.6, zorder=3,
               label=f"in the {N_PRIMARY} ({len(pop)})")
    # this study: one node, linked to the primary studies whose lists it analyses
    # placed at the right edge, vertically in the middle of the lanes
    ty, tx = (len(LANES) - 1) / 2, 2027.6
    THIS_COL = "#7c3aed"  # violet-600, not used elsewhere in the panel
    # curved links (cubic Bezier: leave each study horizontally, arrive at the star level)
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    for _, r in pop.iterrows():
        x0, y0 = r["x"], r["y"]
        dx = tx - x0
        verts = [(x0, y0), (x0 + 0.55 * dx, y0), (tx - 0.35 * dx, ty), (tx, ty)]
        ax.add_patch(PathPatch(Path(verts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                               facecolor="none", edgecolor=THIS_COL, lw=0.45, alpha=0.6, zorder=2.5))
    ax.scatter([tx], [ty], s=60, marker="*", color=THIS_COL, edgecolor=fk.TEXT, lw=0.4, zorder=4,
               label=f"this study (purple lines to the {len(pop)} it analyses)")
    ax.text(tx + 0.5, ty, "this\nstudy", ha="left", va="center", fontsize=fk.SMALL, color=THIS_COL)
    lanes = LANES
    ax.set_yticks(range(len(lanes)))
    ax.set_yticklabels(lanes, fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_title("Every corpus study and its data-reuse links (one point = one study; grey line = shared data or donors)", pad=34)
    ax.set_xlabel("Year of publication")
    ax.set_xlim(1997.5, 2029.6)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, handletextpad=0.1, columnspacing=1, borderaxespad=0.2)
    return k, len(n)


def main():
    c = rd("fig01_b_corpus_studies.tsv")
    c["tgroup"] = c["tissue"].map(TISSUE_OF)
    c["tg"] = c["technology"].map(tech_group)
    inputs = sorted(D.glob("*.tsv")) + [fk.FD / "fig02" / "fig02_b_study_annotation.tsv", FUNNEL, PRISMA]

    # Figure 1: a (schematic), b (lineage)
    fig = fk.figure(112)
    gA = GridSpec(1, 1, figure=fig, left=0.012, right=0.988, top=0.985, bottom=0.745)
    gB = GridSpec(1, 1, figure=fig, left=0.12, right=0.985, top=0.585, bottom=0.135)
    axa, axb = (fig.add_subplot(g[0]) for g in (gA, gB))
    panel_a(axa)
    n_edges, n_nodes = panel_h(axb)
    for ax, l in zip([axa, axb], "ab"):
        fk.letter(ax, l)
    fk.save(fig, "fig_corpus", "figures/scripts/build_fig_corpus.py", inputs,
            notes=f"b (lineage): {n_nodes} corpus studies with year and lane, {n_edges} links drawn")

    # Extended Fig. 2: a to f (the panel_b, panel_c, panel_d, panel_e, panel_f and panel_g functions)
    fig = fk.figure(205)
    gB = GridSpec(1, 3, figure=fig, width_ratios=[1.25, 0.95, 1.0], wspace=0.5, left=0.075, right=0.93, top=0.955, bottom=0.73)
    gC = GridSpec(1, 2, figure=fig, width_ratios=[1.1, 1.1], wspace=0.45, left=0.085, right=0.93, top=0.62, bottom=0.40)
    gD = GridSpec(1, 1, figure=fig, left=0.12, right=0.93, top=0.30, bottom=0.07)
    axa, axb, axc = (fig.add_subplot(gB[i]) for i in range(3))
    axd, axe = (fig.add_subplot(gC[i]) for i in range(2))
    axf = fig.add_subplot(gD[0])
    panel_b(axa, c)
    panel_c(axb, c)
    panel_d(axc, c)
    panel_e(axd)
    panel_f(axe)
    panel_g(axf)
    for ax, l in zip([axa, axb, axc, axd, axe, axf], "abcdef"):
        fk.letter(ax, l)
    fk.save(fig, "supp_corpus", "figures/scripts/build_fig_corpus.py", inputs)


if __name__ == "__main__":
    main()
