#!/usr/bin/env python3
"""Figure 7 -- the recurrent genes track brain pH in controls.

Data: figures/data/fig07/ (2026-10-08).

a  pH, schizophrenia minus control, per dataset-region (14)
b  per gene: association with pH in controls vs raw schizophrenia effect (genes in 5+ lists)
c  the same against the pH-adjusted schizophrenia effect

d  agreement of subject-separate bulk dataset pairs as covariates are added (28 datasets with sample-level data)
e  signed pH score of each study's list by how the study handled pH (lists with >= 10 scored genes)
(the former panel f, recurrent genes' combined z with and without covariates, is Supplementary Fig. 5: build_supp_covariate_genes.py)
"""
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

import figkit as fk

D = fk.FD / "fig07"
rd = lambda f: pd.read_csv(D / f, sep="\t", comment="#")
rng = np.random.default_rng(1)


def panel_a(ax):
    df = rd("fig07_a1_ph_difference.tsv")
    df = df[df["in_14_dataset_regions"] == "yes"].sort_values("delta_pH")
    cols = [fk.ACCENT if p < 0.05 else fk.lighten(fk.ACCENT, 0.6) for p in df["p"]]
    y = np.arange(len(df))
    ax.barh(y, df["delta_pH"], color=cols, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(df["dataset_region"].str.replace("_", " "), fontsize=fk.SMALL)
    ax.invert_yaxis()
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlabel("Brain pH difference (mean SZ − mean control)")
    ax.tick_params(axis="y", length=0)
    ax.set_title("pH difference in the 14\ndatasets with expression data", pad=4)
    ax.text(0.03, 0.03, "one bar = one dataset\n(brain region)\ndark: P < 0.05", transform=ax.transAxes, fontsize=fk.SMALL, color=fk.SUBTLE, va="bottom")


# label positions in data coordinates (x = pH association, y = SZ effect), each in an empty part of the panel and
# on a white background so the dots do not run through the letters
MARKERS = [("SST", 6.0, -7.4), ("PVALB", 3.0, -7.4), ("GAD1", 5.8, 1.7), ("AQP4", -2.2, 6.6), ("GFAP", -5.6, -2.2)]


def scatter(ax, ycol, ylabel, labels=()):
    # the 5-or-more population (ge5 file) plus any recurrent gene outside it (main file); in_176 marks the 375
    df = pd.concat([rd("fig07_a2_a4_gene_scatter.tsv"), rd("fig07_a2_a4_gene_scatter_ge5.tsv")]).drop_duplicates("gene")
    df = df.dropna(subset=["z_ph_controls", "z_sz_raw", "z_sz_ph_adjusted"])
    rec = df["in_176"] == "yes"
    ax.scatter(df.loc[~rec, "z_ph_controls"], df.loc[~rec, ycol], s=1.5, color=fk.NEUTRAL, lw=0,
               rasterized=True, label=f"other genes in\n5+ lists ({(~rec).sum():,})")
    ax.scatter(df.loc[rec, "z_ph_controls"], df.loc[rec, ycol], s=3, color=fk.ACCENT, lw=0,
               label=f"recurrent ({rec.sum()})")
    # line, r and n are for the recurrent genes (as in the text); the grey genes are shown for context only
    x, y = df.loc[rec, "z_ph_controls"].to_numpy(), df.loc[rec, ycol].to_numpy()
    b, a = np.polyfit(x, y, 1)
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, a + b * xs, color=fk.TEXT, lw=0.8, ls="--")
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.97, 0.96, f"r = {r:.2f}\nn = {rec.sum():,} recurrent", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL,
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.9), zorder=7)
    # a few cell-type marker genes named: interneuron markers high-pH / down in SZ,
    # astrocyte markers low-pH / up in SZ
    for g, dx, dy in labels:
        row = df[df["gene"] == g]
        if len(row):
            gx, gy = float(row["z_ph_controls"].iloc[0]), float(row[ycol].iloc[0])
            ax.scatter([gx], [gy], s=12, facecolor="none", edgecolor=fk.TEXT, lw=0.7, zorder=5)
            ax.annotate(g, (gx, gy), xytext=(dx, dy), textcoords="data", fontsize=fk.SMALL,
                        color=fk.TEXT, ha="center", va="center", zorder=6,
                        bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.9),
                        arrowprops=dict(arrowstyle="-", color=fk.TEXT, lw=0.4, shrinkA=6, shrinkB=2))
    ax.axhline(0, color=fk.GRID, lw=0.6)
    ax.axvline(0, color=fk.GRID, lw=0.6)
    ax.set_xlabel("Gene's pH association\nin controls (combined z)")
    ax.set_ylabel(ylabel)
    return r, int(rec.sum())


def panel_d(ax):
    df = rd("fig07_a3_observed_vs_expected.tsv").set_index("set")
    # cohort-cluster bar dropped: the cohort-cluster set depends on a few studies
    sets = [("study_brain", "recurrent"), ("study_brain_5plus", "5+ lists")]
    x = np.arange(len(sets))
    fold = [df.loc[s, "observed"] / df.loc[s, "expected"] for s, _ in sets]
    ax.bar(x, fold, 0.6, color=fk.ACCENT)
    for xi, (s, _) in zip(x, sets):
        r = df.loc[s]
        ax.text(xi, r["observed"] / r["expected"] + 0.1, f"{int(r['observed'])}\n/ {r['expected']:.3g}",
                ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{l}\n{int(df.loc[s, 'n_measured']):,}" for s, l in sets], fontsize=fk.SMALL)
    ax.set_ylabel("pH-associated genes in the set\n(observed / expression-matched expectation)")
    ax.set_xlabel("Gene set (genes measured)")
    ax.set_title("The recurrent genes are enriched\nfor pH-associated genes", pad=4)
    ax.set_ylim(0, max(fold) * 1.3)


# numbers: of the 22 bulk datasets, how many record the covariate
STEPS = [("L0", "none"), ("L1", "age,\nsex"), ("L2", "+PMI\n(22)"), ("L2p", "+pH\n(18)"), ("L3", "+RIN (10),\nagonal\n(2)"), ("L5b", "+medica-\ntion (4),\nSVs")]


def panel_e(ax):
    df = rd("fig07_b_covariate_series_pairs.tsv")
    df = df[(df["relation"] == "subject-separate") & (df["both_bulk"].astype(str).str.upper() == "TRUE")]
    w = df.pivot(index="pair_id", columns="model_step", values="rho")
    # pH alone after PMI; fig07_a_covariate_series_pH_RIN_steps.tsv. Pairs where either dataset lacks pH keep their PMI-level value, the same rule as the other steps.
    st = rd("fig07_a_covariate_series_pH_RIN_steps.tsv")
    st = st[st["model_step"] == "L2p"].set_index("pair_id")["rho"]
    w["L2p"] = st.reindex(w.index)
    w["L2p"] = w["L2p"].fillna(w["L2"])
    w = w[[s for s, _ in STEPS]].dropna()
    x = np.arange(len(STEPS))
    for _, r in w.iterrows():
        ax.plot(x, r.to_numpy(), color=fk.lighten(fk.SUBTLE, 0.75), lw=0.3)
    ax.plot(x, w.median().to_numpy(), "-o", color=fk.ACCENT, ms=2.5, lw=1.2)
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in STEPS], fontsize=fk.SMALL)
    ax.set_ylabel("Agreement between two datasets\n(Spearman ρ of SZ t statistics)")
    ax.set_xlabel("Covariates in the model, cumulative\n(datasets of 22 that record the covariate)")
    ax.set_title("Adjusting for pH lowers study agreement", pad=4)
    ax.set_ylim(w.min().min() - 0.02, w.max().max() + 0.08)
    ax.text(0.98, 0.97, f"one line = one of {len(w)} pairs of datasets\nwith no subjects in common; blue: median", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)
    return len(w), w.median().round(3).to_dict()


def panel_covariates(ax):
    df = rd("fig07_b_covariate_series_genes.tsv")
    for d, col in [("up", fk.UP), ("down", fk.DOWN), ("mixed", fk.NEUTRAL)]:
        s = df[df["sz_direction"] == d]
        ax.scatter(s["z_L1"], s["z_L3"], s=3, color=col, lw=0, label=f"{d} in SZ")
    lim = [min(df["z_L1"].min(), df["z_L3"].min()) - 0.3, max(df["z_L1"].max(), df["z_L3"].max()) + 0.3]
    ax.plot(lim, lim, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Recurrent gene's SZ effect\n(combined z, age and sex adjusted)")
    ax.set_ylabel("SZ effect with all covariates\n(combined z, 7 datasets)")
    ax.legend(loc="upper left", handletextpad=0.1, markerscale=1.5, title="one point = one recurrent gene", title_fontsize=fk.SMALL)
    ax.set_title("Covariates shrink the recurrent genes' effects", pad=4)


HANDLING = [("adjusted", "adjusted", ["adjusted covariate"]),
            ("sens", "sensitivity\nonly", ["sensitivity only", "not covariate in primary; sensitivity with pH"]),
            ("matched", "matched", ["matched"]),
            ("reported", "reported\nonly", ["reported only"]),
            ("discussed", "discussed,\nnot used", ["discussed as confound without adjustment", "considered as candidate covariate"]),
            ("silent", "not\nmentioned", ["silent", "unread"])]


def panel_g(ax):
    # per-study pH handling re-read from full texts and supplements
    df = rd("fig07_c_ph_score_by_handling_v2.tsv")
    df = df[df["scored_ge10"] == "yes"]
    for i, (_, _, calls) in enumerate(HANDLING):
        y = df.loc[df["ph_call"].isin(calls), "signed_ph_score"].to_numpy()
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(y)), y, s=16, color=fk.ACCENT_3, lw=0)
        if len(y):
            ax.plot([i - 0.3, i + 0.3], [np.median(y)] * 2, color=fk.TEXT, lw=1.2)
    ax.axhline(0, color=fk.SUBTLE, lw=0.6, ls=":")
    ax.set_xticks(range(len(HANDLING)))
    ax.set_xticklabels([f"{l}\n({df['ph_call'].isin(c).sum()})" for _, l, c in HANDLING], fontsize=fk.SMALL)
    ax.set_ylabel("pH signature in the study's list\n(mean control-pH z, signed by reported direction)")
    ax.set_xlabel("How the study handled pH (studies)")
    ax.set_title("pH score of each study's list,\nby how the study handled pH", pad=4)
    ax.text(0.98, 0.97, "one point = one study; bar: median", transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_f(ax):
    """pH of cases minus controls by cohort: all persons, persons in the primary studies, persons in the 14 dataset-regions."""
    d = pd.read_csv(fk.FD / "ph_by_cohort" / "per_cohort_80.tsv", sep="\t", comment="#")
    pool = d[d["cohort"].str.startswith("all cohorts")].iloc[0]
    rest = d[~d["cohort"].str.startswith("all cohorts")].copy()
    order = ["CMC/MSSM/Pitt/LIBD", "Australian", "Stanley", "Pritzker", "other single-study", "Harvard"]
    rest = rest.set_index("cohort").loc[order].reset_index()
    rows = [("All cohorts (mixed model)", pool, True)] + [(r["cohort"].replace("MSSM", "Mount Sinai"), r, False) for _, r in rest.iterrows()]
    y = np.arange(len(rows))[::-1].astype(float)
    y[1:] -= 0.6
    def draw(yi, r, pre, marker, colour, mfc, dy, size):
        n_sz, n_c = r[f"{pre}_n_SZ"], r[f"{pre}_n_CTRL"]
        if pd.isna(n_sz) or n_sz == 0:
            return
        est = r[f"{pre}_model_est"] if (pre in ("A", "B") and pooled_row(r)) else r[f"{pre}_est"]
        lo = r[f"{pre}_model_lo"] if (pre in ("A", "B") and pooled_row(r)) else r[f"{pre}_lo"]
        hi = r[f"{pre}_model_hi"] if (pre in ("A", "B") and pooled_row(r)) else r[f"{pre}_hi"]
        if pd.isna(est):
            return
        ax.plot([lo, hi], [yi + dy] * 2, color=colour, lw=1.1)
        ax.plot(est, yi + dy, marker, color=colour, mfc=mfc, ms=size)
    def pooled_row(r):
        return str(r["cohort"]).startswith("all cohorts")
    for yi, (lab, r, is_pool) in zip(y, rows):
        draw(yi, r, "B", "o", fk.ACCENT, fk.ACCENT, 0.2, 4.2)
        draw(yi, r, "A", "o", fk.ACCENT_3, "white", 0.0, 4.2)
        draw(yi, r, "expr14", "D", fk.TEXT, fk.TEXT, -0.2, 3.4)
    ax.axvline(0, color=fk.SUBTLE, lw=0.7, ls=(0, (3, 2)))
    def n_(r, pre):
        return "" if pd.isna(r[f"{pre}_n_SZ"]) or r[f"{pre}_n_SZ"] == 0 else f"{int(r[f'{pre}_n_SZ'])}/{int(r[f'{pre}_n_CTRL'])}"
    labels = ["All cohorts" if is_pool else lab for lab, _, is_pool in rows]  # people per group: Supplementary Note 6
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fk.SMALL)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Brain pH, schizophrenia minus control (95% interval)")
    ax.set_title("Brain pH difference by cohort", pad=18)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=fk.ACCENT, marker="o", lw=1.1, label="all with a pH"),
                       Line2D([], [], color=fk.ACCENT_3, marker="o", mfc="white", lw=1.1, label="analysed studies"),
                       Line2D([], [], color=fk.TEXT, marker="D", ms=3.4, lw=1.1, label="datasets in panel a")],
              loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, handlelength=1.2, columnspacing=0.8, fontsize=fk.SMALL, borderaxespad=0.1)


def panel_subclass(ax):
    df = pd.read_csv(fk.FD / "fig09" / "fig09_c_ph_vs_sn_effect_by_subclass.tsv", sep="\t", comment="#").sort_values("r_ph_vs_sn_effect_ge10_universe")
    y = np.arange(len(df))
    ax.scatter(df["r_ph_vs_sn_effect_genomewide"], y, s=12, facecolor="white", edgecolor=fk.ACCENT_3, lw=0.8,
               zorder=3, label="all genes")
    ax.scatter(df["r_ph_vs_sn_effect_ge10_universe"], y, s=14, color=fk.ACCENT_3, edgecolor="white",
               lw=0.5, zorder=4, label="genes in 10+ lists")
    ax.set_yticks(y)
    ax.set_yticklabels(df["subclass"], fontsize=fk.TINY)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.tick_params(axis="y", length=0)
    n_neg = int((df["r_ph_vs_sn_effect_ge10_universe"] < 0).sum())
    ax.set_title(f"pH signature in {n_neg} of {len(df)}\nsingle-nucleus subclasses", pad=4)
    ax.set_xlabel("Correlation of pH association\nwith single-nucleus SZ effect (r)")
    ax.legend(loc="upper left", fontsize=fk.TINY, handletextpad=0.2, borderaxespad=0.2)
    ax.set_ylim(-0.7, len(df) - 0.3)
    lo = min(df["r_ph_vs_sn_effect_ge10_universe"].min(), df["r_ph_vs_sn_effect_genomewide"].min())
    hi = max(df["r_ph_vs_sn_effect_ge10_universe"].max(), df["r_ph_vs_sn_effect_genomewide"].max())
    ax.set_xlim(min(-0.29, lo - 0.01), max(0.05, hi + 0.01))
    return df


def main():
    # former panel d (two observed/expected ratios) dropped, numbers stay in the text
    inputs = sorted(D.glob("*.tsv"))
    # Figure 3: a, b and c (the panel_f function: pH by cohort)
    fig = fk.figure(150)
    r1 = GridSpec(1, 2, figure=fig, width_ratios=[1, 1], wspace=0.7, left=0.16, right=0.78, top=0.925, bottom=0.62)
    r3 = GridSpec(1, 1, figure=fig, left=0.20, right=0.78, top=0.43, bottom=0.075)
    axa, axb = (fig.add_subplot(r1[i]) for i in range(2))
    axc = fig.add_subplot(r3[0])
    panel_a(axa)
    rb, n = scatter(axb, "z_sz_raw", "Gene's SZ effect\n(combined z, age and sex adjusted)", labels=MARKERS)
    axb.set_ylim(-9, axb.get_ylim()[1])  # room below the points for the labels
    axb.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), markerscale=2, handletextpad=0.1, frameon=False, fontsize=fk.SMALL,
               title="one point = one gene", title_fontsize=fk.SMALL, alignment="left")
    axb.set_title("Schizophrenia effect against\npH association", pad=4)
    panel_f(axc)
    for ax, l in zip([axa, axb, axc], "abc"):
        fk.letter(ax, l)
    fk.save(fig, "fig_ph", "figures/scripts/build_fig_ph.py", inputs, notes=f"b r={rb:.3f} n={n}")

    # Extended Fig. 4: a (the panel_e function: agreement as covariates are added), b (the panel_g function: pH score by handling)
    fig = fk.figure(75)
    r2 = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.15], wspace=0.35, left=0.09, right=0.985, top=0.88, bottom=0.14)
    axa, axb = fig.add_subplot(r2[0]), fig.add_subplot(r2[1])
    n_pairs, medians = panel_e(axa)
    panel_g(axb)
    for ax, l in zip([axa, axb], "ab"):
        fk.letter(ax, l)
    fk.save(fig, "supp_ph_covariates", "figures/scripts/build_fig_ph.py", inputs,
            notes=f"a {n_pairs} subject-separate bulk pairs, medians {medians}")


if __name__ == "__main__":
    main()
