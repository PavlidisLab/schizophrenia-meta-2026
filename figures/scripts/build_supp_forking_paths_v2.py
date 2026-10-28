#!/usr/bin/env python3
"""Supplementary Fig. 2 (v2): the recurrence result over the distinct ways of building the lists.

Data: figures/data/forking_paths_v2/ (analyses/forking_paths_v2_2026-10-10/). One analysis = one distinct hash of post-fork inputs.
a  the choices and their settings (forking_v2_a_design.tsv)
b  genes in 10 or more lists, observed over expected under each chance model, one point per distinct analysis,
   published lists only (forking_v2_b_spec_curve.tsv; study tier)
c  the same against our DE prior, published lists only, by which studies are included
d  signed control-pH score of the genes each analysis selects (forking_v2_c_ph_paths.tsv; study tier)
"""
import numpy as np
import pandas as pd

import figkit as fk

D = fk.FD / "forking_paths_v2"
FILES = {k: D / f for k, f in {"a": "forking_v2_a_design.tsv", "e": "forking_v2_e_paths.tsv", "b": "forking_v2_b_spec_curve.tsv",
                                 "c": "forking_v2_c_ph_paths.tsv"}.items()}
rd = lambda k: pd.read_csv(FILES[k], sep="\t", comment="#")

NULLS = [("uniform", "random genes", fk.ACCENT), ("expr", "expression-matched random genes", fk.ACCENT_2),
         ("de_prior_ours", "DE prior (ours)", fk.ACCENT_3), ("de_prior_global", "DE prior (external)",
                                                             fk.lighten(fk.ACCENT_3, 0.45))]
CUTS = [("published", "as published"), ("fc026", "1.2-fold"), ("fc058", "1.5-fold"), ("tight", "tighter threshold"),
        ("cap200", "top 200 genes")]
FILTERS = [("all", f"all {fk.n_primary()}"), ("n20", "N ≥ 20"), ("from2016", "2016 on"), ("cortex", "cortex"),
           ("covariate_adjusted", "covariate-\nadjusted")]  # error-study setting removed from the grid
SET_LAB = {"list_level": "which lists are kept", "rule": "list rule", "threshold": "threshold", "empty": "empty lists (FDR only)",
           "pseudo_noexpr": "pseudogenes / not in brain", "noncoding": "non-coding RNA",
           "single_cell": "single-nucleus studies", "cut": "list cut", "study_filter": "studies included"}


def distinct(b, tier="study"):
    return b[b["tier"] == tier].drop_duplicates(["hash_id", "null"])


def panel_a(ax):
    a = rd("a")
    ax.set_axis_off()
    e = rd("e")
    nb = e.groupby("block")["path_id"].nunique()
    ax.set_title(f"Choices crossed ({len(e)} combinations, {e['hash_id'].nunique()} giving distinct lists)", loc="left")
    y = 0.95
    for fork, g in a.groupby("fork", sort=False):
        nice = {"fdr": "computed FDR only", "nominal": "nominal significance", "all": "any author criterion",
                "keep": "keep", "drop": "drop", "gtex": "keep only genes expressed in GTEx brain", "pooled": "pool with coding",
                "without": "leave out", "with": "add", "R1": "R1 (one list per study)", "A": "A (union of a study's lists)",
                "inclusive": "all", "no_partial": "drop partial", "no_nonDE": "drop non-DE-style",
                "no_mixed": "drop mixed-diagnosis", "chang_s2only": "Chang S2 only",
                "pai_text": "Pai text genes", "no_flagged": "all three drops"}
        if fork in ("cut", "study_filter"):
            nice.update({k: v.replace("-\n", "-").replace("\n", " ") for k, v in CUTS + FILTERS})
        sets = "; ".join(nice.get(str(x), str(x)) for x in g["setting"])
        ax.text(0.0, y, SET_LAB.get(fork, fork), fontsize=fk.SMALL, fontweight="bold", color=fk.TEXT, va="top",
                transform=ax.transAxes)
        ax.text(0.30, y, sets, fontsize=fk.SMALL, color=fk.TEXT, va="top", transform=ax.transAxes)
        y -= 0.1
    ax.text(0.0, y - 0.02, f"Block 1: rule x threshold x pseudogenes x non-coding x single-nucleus ({nb.get('block1', 0)}).\n"
                           f"Block 2: rule x threshold x list cut x studies included ({nb.get('block2', 0)}).\n"
                           f"Block 3: rule x threshold x which lists are kept ({nb.get('block3', 0)}).",
            fontsize=fk.SMALL, color=fk.SUBTLE, va="top", transform=ax.transAxes)


def panel_b(ax, b):
    d = distinct(b)
    d = d[d["cut"] == "published"]
    for j, (nk, nl, col) in enumerate(NULLS):
        v = d[d["null"] == nk]["obs_over_exp_ge10"]
        x = j + np.random.default_rng(j).uniform(-0.18, 0.18, len(v))
        ax.scatter(x, v, s=6, color=col, lw=0, alpha=0.7)
        ax.hlines(v.median(), j - 0.28, j + 0.28, color=fk.TEXT, lw=1)
    ax.set_yscale("log")
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xticks(range(len(NULLS)))
    ax.set_xticklabels([nl.replace("expression-matched random genes", "expression-matched\nrandom genes") for _, nl, _ in NULLS], fontsize=fk.SMALL)
    ax.set_xlabel(f"Chance model (published lists; {d['hash_id'].nunique()} distinct analyses)")
    ax.set_ylabel("Genes in 10 or more lists\n(observed over expected, log)")
    ax.set_title("Published lists, by chance model", loc="left")
    ax.text(0.99, 0.97, "one point = one analysis; bar, median", transform=ax.transAxes, ha="right", va="top",
            fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_c(ax, b):
    d = distinct(b)
    d = d[(d["null"] == "de_prior_ours") & (d["cut"] == "published")]
    for i, (fk_, _) in enumerate(FILTERS):
        v = d[d["study_filter"] == fk_]["obs_over_exp_ge10"]
        x = i + np.random.default_rng(i).uniform(-0.15, 0.15, len(v))
        ax.scatter(x, v, s=6, color=fk.ACCENT_3, lw=0)
        ax.hlines(v.median(), i - 0.25, i + 0.25, color=fk.TEXT, lw=1)
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls="--")
    n = d.groupby("study_filter").size()
    ax.set_xticks(range(len(FILTERS)))
    ax.set_xticklabels([f"{l}\n({n.get(k, 0)})" for k, l in FILTERS], fontsize=fk.SMALL)
    ax.set_yscale("log")
    ax.set_ylabel("Genes in 10 or more lists\n(observed over DE prior, log)")
    ax.set_xlabel("Studies included (published lists; analyses)")
    ax.set_title("Published lists, by studies included", loc="left")


def panel_d(ax):
    c = rd("c")
    c = c[c["tier"] == "study"].drop_duplicates("hash_id").merge(
        rd("b")[["hash_id", "cut"]].drop_duplicates("hash_id"), on="hash_id", how="left")
    for i, (ck, _) in enumerate(CUTS):
        v = c[c["cut"] == ck]["signed_ph_score"].dropna()
        x = i + np.random.default_rng(100 + i).uniform(-0.15, 0.15, len(v))
        ax.scatter(x, v, s=6, color=fk.ACCENT, lw=0)
        ax.hlines(v.median(), i - 0.25, i + 0.25, color=fk.TEXT, lw=1)
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xticks(range(len(CUTS)))
    ax.set_xticklabels([l for _, l in CUTS], fontsize=fk.SMALL, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Signed pH score of the\nrecurrent genes")
    ax.set_title("Every analysis carries\nthe pH signature", loc="left")
    n_pos = int((c["signed_ph_score"] > 0).sum())
    ax.text(0.02, 0.97, f"positive in {n_pos} of {int(c['signed_ph_score'].notna().sum())}", transform=ax.transAxes,
            fontsize=fk.SMALL, color=fk.SUBTLE, va="top")


def main():
    fk.apply_style()
    b = rd("b")
    fig = fk.figure(205)
    axa = fig.add_axes([0.07, 0.74, 0.90, 0.21])
    axb = fig.add_axes([0.11, 0.42, 0.86, 0.23])  # raised: its x label cleared panel letter d
    # bottom row raised so panel c's x label and panel d's tick labels stay on the page
    axc = fig.add_axes([0.11, 0.10, 0.48, 0.19])
    axd = fig.add_axes([0.70, 0.10, 0.26, 0.19])
    panel_a(axa)
    panel_b(axb, b)
    panel_c(axc, b)
    panel_d(axd)
    for ax, l in zip([axa, axb, axc, axd], "abcd"):
        fk.letter(ax, l)
    fk.save(fig, "supp_forking_paths_v2", "figures/scripts/build_supp_forking_paths_v2.py", list(FILES.values()),
            notes="Forking paths v2; study tier; distinct analyses by input hash.")


if __name__ == "__main__":
    main()
