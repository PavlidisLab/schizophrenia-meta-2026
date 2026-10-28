#!/usr/bin/env python3
"""Supplementary Note 4 Figure 1 -- overlap of 105 pairs of brain meta-analyses against the share of input datasets they have in common.

Data: figures/data/fig02/fig02_e_meta_analysis_pairs.tsv
"""
import build_fig_agreement as fa
import figkit as fk

def main():
    fig = fk.figure(90)
    axc = fig.add_axes([0.2, 0.20, 0.75, 0.66])
    fa.panel_d(axc)
    axc.set_title("Overlap of meta-analysis pairs\nby shared input datasets", pad=6)
    fk.save(fig, "supp_study_pairs", "figures/scripts/build_supp_study_pairs.py", [fa.D / "fig02_e_meta_analysis_pairs.tsv"])


if __name__ == "__main__":
    main()
