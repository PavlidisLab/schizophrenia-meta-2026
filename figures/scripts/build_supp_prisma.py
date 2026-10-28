#!/usr/bin/env python3
"""Supplementary Fig. (PRISMA) -- literature flow, PRISMA 2020 style, from the search to the corpus and on to the primary brain studies.

Data: figures/data/prisma_flow/prisma_flow.tsv
(analyses/prisma_flow_2026-10-10/build_prisma_flow.py; arithmetic checked in prisma_flow_check.tsv).
Every number is read from the TSV by box_id; every box is a final count on the frozen corpus.
The flow's arithmetic is asserted below as in prisma_flow_check.tsv.
"""
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import figkit as fk

F = fk.FD / "prisma_flow" / "prisma_flow.tsv"
T = pd.read_csv(F, sep="\t", comment="#").set_index("box_id")
n = lambda b: int(T.loc[b, "n"])
assert (T["status"] == "final").all(), T.loc[T["status"] != "final", "status"]

ROUTES = [("ID01", "databases (PubMed, Embase, Web of Science)"), ("ID02", "single-cell / single-nucleus search"),
          ("ID03", "reference lists of included studies"), ("ID04", "citing papers"),
          ("ID05", "Gemma datasets"), ("ID06", "genome-wide and proteomic additions"),
          ("ID07", "single-gene literature sweep"), ("ID08", "other (saved records, data lineage)")]


def box(ax, x, y, w, h, text, col=fk.TEXT, fc="white", size=6.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.008", facecolor=fc,
                                edgecolor=col, lw=0.8))
    ax.text(x + 0.012, y + h - 0.012, text, fontsize=size, color=fk.TEXT, va="top", ha="left", linespacing=1.3)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, color=fk.SUBTLE, lw=0.8))


def check():
    n_ = lambda b: int(T.loc[b, "n"])
    assert n_("ID09") - n_("ID10") == n_("ID11")
    assert n_("EL01") == n_("SC05") - n_("EL00a") + n_("EL00b")
    assert n_("EL03") == n_("EL01") - n_("EL02") == n_("EL04") + n_("EL05") + n_("EL06") + n_("EL07")
    assert n_("IN00") == n_("EL00a") + n_("EL06") + n_("EL05")
    assert n_("IN01") == n_("IN00") - n_("IN00a") + n_("IN00b")
    assert n_("EL07") == sum(n_(f"EL07.0{k}") for k in range(1, 7))
    assert n_("SC06") == n_("SC02") + n_("SC04")
    assert n_("EL04") == sum(n_(b) for b in T.index if b.startswith("EL04."))


def main():
    check()
    fig = fk.figure(235)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    L, W, R, RW = 0.10, 0.44, 0.60, 0.38
    L2 = L + 0.15  # main eligibility column, right of the box for studies included on title and abstract
    stages = [(0.985, 0.80, "Identification"), (0.78, 0.60, "Screening"), (0.58, 0.245, "Eligibility"),
              (0.235, 0.005, "Included")]
    for top, bot, lab in stages:
        ax.add_patch(FancyBboxPatch((0.005, bot), 0.06, top - bot, boxstyle="round,pad=0,rounding_size=0.008",
                                    facecolor=fk.lighten(fk.ACCENT, 0.85), edgecolor="none"))
        ax.text(0.035, (top + bot) / 2, lab, rotation=90, ha="center", va="center", fontsize=7, color=fk.TEXT)
    # identification
    routes = "\n".join(f"  {lab}: {n(b):,}" for b, lab in ROUTES)
    box(ax, L, 0.81, W, 0.17, f"Records identified, by route ({n('ID09'):,}; a record can come by several routes)\n{routes}")
    box(ax, R, 0.865, RW, 0.06, f"Duplicate records removed: {n('ID10'):,}\nNo title or abstract: {int(T.loc['M01','n'])}")
    arrow(ax, L + W, 0.895, R, 0.895)
    # screening (models named in Extended Methods 1)
    box(ax, L, 0.70, W, 0.07, f"Records screened on title and abstract: {n('SC01'):,}\n"
                              f"  (pass 1, claude-haiku-5-5; pass 2 on {n('SC03'):,} records, claude-sonnet-5-5;\n   October 2026 rescreen of the saved pool)")
    arrow(ax, L + W / 2, 0.81, L + W / 2, 0.77)
    box(ax, R, 0.70, RW, 0.07, f"Records excluded: {n('SC06'):,}\n  first pass {n('SC02'):,}; second pass {n('SC04'):,}\n"
                               f"Not screened (no usable result): {n('SC00')}")
    arrow(ax, L + W, 0.735, R, 0.735)
    box(ax, L, 0.61, W, 0.07, f"Advanced after title and abstract: {n('SC05'):,}\n"
                              f"  ({n('SC05a'):,} include, {n('SC05b'):,} unsure)")
    arrow(ax, L + W / 2, 0.70, L + W / 2, 0.68)
    # eligibility
    box(ax, L, 0.455, 0.135, 0.125, f"Included on title\nand abstract, no\nseparate full-text\nassessment\n(corpus studies):\n{n('EL00a'):,}",
        col=fk.ACCENT)
    arrow(ax, L + 0.0675, 0.61, L + 0.0675, 0.58)
    box(ax, L2, 0.49, L + W - L2, 0.09, f"Reports sought for full-text review: {n('EL01'):,}\n"
                                         f"  {n('SC05') - n('EL00a'):,} advanced, not yet included\n"
                                         f"  {n('EL00b'):,} excluded on title and abstract,\n  re-checked on full text")
    arrow(ax, L2 + (L + W - L2) / 2, 0.61, L2 + (L + W - L2) / 2, 0.58)
    box(ax, R, 0.535, RW, 0.04, f"Reports not retrieved: {n('EL02'):,}")
    arrow(ax, L + W, 0.555, R, 0.555)
    box(ax, L2, 0.40, L + W - L2, 0.065, f"Full-text reports assessed: {n('EL03'):,}")
    arrow(ax, L2 + (L + W - L2) / 2, 0.49, L2 + (L + W - L2) / 2, 0.465)
    excl = [b for b in T.index if b.startswith("EL04.")]
    nice = {"new_data": "new data", "reanalysis": "reanalysis", "meta_analysis": "meta-analysis", "other": "other",
            "targeted_few_genes": "few genes", "transcriptome_wide": "genome-wide", "unclear": "scope unclear"}
    def lab(b):
        du, sc = T.loc[b, "label"].split(": ", 1)[1].replace("data_use=", "").split(", scope=")
        return f"{nice.get(du, du)}, {nice.get(sc, sc)}"
    ex_txt = "\n".join(f"  {lab(b)}: {n(b)}" for b in excl)
    box(ax, R, 0.325, RW, 0.2, f"Full-text reports excluded: {n('EL04'):,}\n(use of data; scope)\n{ex_txt}")
    arrow(ax, L + W, 0.44, R, 0.44)
    box(ax, R, 0.224, RW, 0.098, f"Reports of studies not yet in the corpus: {n('EL07'):,}\n"
                                 f"  added to the corpus: {n('EL07.01'):,}\n"
                                 f"  met the criteria, then set aside (fixed gene panels): {n('EL07.06'):,}\n"
                                 f"  version of a study already in, or already a node: {n('EL07.02') + n('EL07.03'):,}\n"
                                 f"  no schizophrenia-versus-control comparison: {n('EL07.04'):,}\n"
                                 f"  few-gene (targeted) studies, set aside: {n('EL07.05'):,}")
    arrow(ax, L + W, 0.415, R, 0.29)
    box(ax, L2, 0.265, L + W - L2, 0.105, f"Full-text reports of corpus studies\n"
                                          f"  confirming them: {n('EL06'):,}\n"
                                          f"  not meeting the full-text criteria,\n  retained in the corpus: {n('EL05'):,}",
        col=fk.ACCENT)
    arrow(ax, L2 + (L + W - L2) / 2, 0.40, L2 + (L + W - L2) / 2, 0.37)
    # included
    box(ax, L, 0.165, W, 0.055, f"Studies taken into the corpus: {n('IN00'):,}\n"
                               f"  ({n('EL00a')} + {n('EL06')} + {n('EL05')})", col=fk.ACCENT)
    arrow(ax, L + 0.0675, 0.455, L + 0.0675, 0.22)
    arrow(ax, L2 + (L + W - L2) / 2, 0.265, L2 + (L + W - L2) / 2, 0.22)
    box(ax, R, 0.17, RW, 0.045, f"Excluded: targeted gene-panel, small-set\n  or decided-out studies: {n('IN00a')}")

    arrow(ax, L + W, 0.1925, R, 0.1925)
    box(ax, L, 0.075, W, 0.065, f"Studies included (corpus): {n('IN01'):,}\n  ({n('IN00')} - {n('IN00a')}, + {n('IN00b.02')} added from the rescreen\n   and {n('IN00b.01')} by hand)", col=fk.ACCENT)
    arrow(ax, L + W / 2, 0.165, L + W / 2, 0.14)
    funnel = [("PF02", "with a primary gene list", "PF02x", "without"),
              ("PF03", "differential expression from new data", "PF03x", "left out"),
              ("PF04", "not meta-analyses or reanalyses", "PF04x", "removed"),
              ("PF06", "postmortem brain studies with a list", "PF06x", "other"),
              ("PF07", "donor groups", None, None), ("PF09", "modelled (tested genes known)", None, None),
              ("PF10", f"studies analysed (modelled, plus {n('PF10') - n('PF09')} derived-list)", None, None)]
    ft = "\n".join(f"  {lab}: {n(b):,}" + (f" ({n(x):,} {w})" if x else "") for b, lab, x, w in funnel)
    box(ax, R, 0.005, RW, 0.15, f"From the corpus to the studies analysed\n{ft}", col=fk.ACCENT)
    arrow(ax, L + W, 0.1175, R, 0.1175)
    ax.text(L, 0.045, "Route counts are records in the screening pool,\nafter the earlier searches' own deduplication.",
            fontsize=6, color=fk.SUBTLE, va="top")
    fk.save(fig, "supp_prisma", "figures/scripts/build_supp_prisma.py", [F],
            notes="final counts on the frozen corpus; arithmetic asserted in check().")


if __name__ == "__main__":
    main()
