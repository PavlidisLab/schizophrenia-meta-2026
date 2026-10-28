#!/usr/bin/env python3
"""Extended Fig. 6 -- where the recurrent genes stand out across Gemma experiments.

Data: a scan of one contrast per public human Gemma experiment (picked by build_test_prior.pick, 2026-10-10), summarized in
figures/data/gemma_176_scan_figure/ and figures/data/gemma_209_scan_calibrated/.
Panels a and c use the calibrated per-contrast rule: a contrast is a hit when the recurrent set's z is above the 97.5th percentile of that contrast's own 200
DE-prior-decile-matched comparison sets and AUROC > 0.5. The dotted line is the comparison sets' own hit rate (~3%).
a  share of hits by tissue class (Wilson 95%)                         panel_a.tsv (calibrated)
b  median |Spearman rho| of each contrast's t with the control-pH and GTEx slow-death axes (bootstrap 95%)  panel_b.tsv
c  schizophrenia experiments, other experiments, hypoxia/ischemia/acidosis comparisons, ischemic disease  panel_c.tsv (calibrated)

"""
import pandas as pd

import figkit as fk

DA = fk.FD / "gemma_176_scan_figure"          # panel b (does not depend on the hit rule)
DC = fk.FD / "gemma_209_scan_calibrated"      # panels a and c, calibrated rule
SET_LABEL = "recurrent genes"
FILES = [DC / "panel_a.tsv", DA / "panel_b.tsv", DC / "panel_c.tsv"]
rd = lambda p: pd.read_csv(p, sep="\t", comment="#")

SHORT = {"Cell lines / in vitro": "cell lines, in vitro", "SZ corpus (picked)": "schizophrenia, in the corpus",
         "SZ-annotated held out": "schizophrenia, not in corpus", "All other experiments": "all other experiments",
         "Hypoxia / ischemia / pH contrast scored": "hypoxia, ischemia or acidosis",
         "Ischemia as disease state (stroke, MI)": "ischemic disease"}

def share_rows(ax, labels, k, n, lo, hi, cols, reference, background=None):
    y = range(len(labels))
    ax.barh(list(y), [a / b for a, b in zip(k, n)], 0.6, color=cols, lw=0)
    ax.errorbar([a / b for a, b in zip(k, n)], list(y), xerr=[[a / b - l for a, b, l in zip(k, n, lo)],
                [h - a / b for a, b, h in zip(k, n, hi)]], fmt="none", ecolor=fk.TEXT, lw=0.7, capsize=1.5)
    for yi, a, b, h in zip(y, k, n, hi):
        ax.text(h + 0.015, yi, f"{a / b:.0%} (n = {b:,})", va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(reference, color=fk.SUBTLE, lw=0.6, ls=":")  # comparison sets' own hit rate
    if background is not None:
        ax.axvline(background, color=fk.ACCENT_3, lw=0.8, ls="--")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.8)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Contrasts where the set exceeds the DE prior")


def panel_a(ax):
    a = rd(DC / "panel_a.tsv")
    cols = [fk.ACCENT if g == "Postmortem brain" else fk.lighten(fk.ACCENT, 0.5) if g == "Other brain" else fk.NEUTRAL
            for g in a["group"]]
    share_rows(ax, [SHORT.get(g, g.lower() if g != "Postmortem brain" else g) for g in a["group"]], a["k"], a["n"],
               a["wilson95_lo"], a["wilson95_hi"], cols, a["comparison_set_rate"].mean())
    ax.set_title("The recurrent genes stand out in postmortem brain", loc="left")


def panel_b(ax):
    b = rd(DA / "panel_b.tsv")
    groups = list(dict.fromkeys(b["group"]))
    for j, (axis, lab, col) in enumerate([("pH axis", "low pH in control brain", fk.ACCENT_3),
                                          ("GTEx axis", "slow death (GTEx)", fk.ACCENT_4)]):
        s = b[b["axis"] == axis].set_index("group").loc[groups]
        y = [i + (j - 0.5) * 0.36 for i in range(len(groups))]
        ax.barh(y, s["median_abs_rho"], 0.34, color=col, lw=0, label=lab)
        ax.errorbar(s["median_abs_rho"], y, xerr=[s["median_abs_rho"] - s["boot95_lo"], s["boot95_hi"] - s["median_abs_rho"]],
                    fmt="none", ecolor=fk.TEXT, lw=0.6, capsize=1.2)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([SHORT.get(g, g.lower() if g != "Postmortem brain" else g) for g in groups])
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Median |Spearman ρ| with each pattern")
    ax.set_title("Postmortem brain contrasts correlate\nwith low pH and slow death", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="lower right")


def panel_c(ax):
    c = rd(DC / "panel_c.tsv")
    c = c[~c["group"].str.startswith("Design flagged")]
    cols = [fk.ACCENT if g.startswith("SZ") else fk.ACCENT_3 if ("ischemia" in g.lower() or "hypoxia" in g.lower())
            else fk.NEUTRAL for g in c["group"]]
    share_rows(ax, [SHORT.get(g, g) for g in c["group"]], c["k"], c["n"], c["wilson95_lo"], c["wilson95_hi"], cols,
               c["comparison_set_rate"].mean())
    ax.set_title("Schizophrenia and oxygen-loss experiments", loc="left")


def main():
    fk.apply_style()
    fig = fk.figure(150)
    # panels a to c are on the current recurrent set

    axa = fig.add_axes([0.20, 0.60, 0.27, 0.32])
    axb = fig.add_axes([0.66, 0.60, 0.30, 0.32])
    axc = fig.add_axes([0.20, 0.10, 0.27, 0.33])
    panel_a(axa)
    panel_b(axb)
    panel_c(axc)
    for ax, l in zip([axa, axb, axc], "abc"):
        fk.letter(ax, l)
    fk.save(fig, "supp_gemma_scan", "figures/scripts/build_supp_gemma_scan.py", FILES,
            notes="panels a and c: calibrated per-contrast rule; panel b does not depend on the hit rule.")


if __name__ == "__main__":
    main()
