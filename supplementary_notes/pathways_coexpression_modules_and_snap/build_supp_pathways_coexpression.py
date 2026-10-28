#!/usr/bin/env python3
"""Supplementary Note 6 figure -- pathways and coexpression modules: no stronger agreement than genes; what agrees tracks pH.

Data: data/suppnote6_pathways_coexpression/ (export_suppnote1.py, 2026-10-09;
PROVISIONAL).
a  GO pathway and gene agreement of the 62 lists over two nulls (suppnote6_a_pathway_vs_gene_agreement.tsv)
b  control-pH association of the pathway term groups (suppnote6_d_pH_test_term_groups.tsv)
c  coexpression claims by whether they re-test the schizophrenia association (suppnote6_e_coexpression_claim_classes.tsv,
   suppnote6_f_retested_claims.tsv)
d  SNAP rediscovery: best match with Ling's factor in schizophrenia datasets
    (panel_k_snap_best_match_chu_vs_sz.tsv)
e  SNAP schizophrenia effect per bank and pooled (suppnote6_g_snap_other_snrnaseq.tsv) beside the specificity
   difference, SNAP cell types minus other cell types (panel_l_snap_specificity.tsv)
e  modules from 20 papers on the signed pH and agonal axes; Gandal 2018 modules (suppnote6_h_..., suppnote6_i_...)
f  the 196 modules by result in other cohorts (same direction as reported, no clear difference, opposite), without and
   with the pH and agonal axis scores in the model (module_retest/pooled_over_modules.tsv)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figkit as fk

ND = fk.note_data(__file__)  # this note's own data; shared tables are read from fk.FD

D = ND / "suppnote6_pathways_coexpression"
if not D.exists():
    D = fk.FD / "suppnote1_pathways_coexpression"
FILES = {k: D / f for k, f in {"a": "suppnote6_a_pathway_vs_gene_agreement.tsv",
                                 "d": "suppnote6_d_pH_test_term_groups_ten_regions.tsv",
                                 "e": "suppnote6_e_coexpression_claim_classes.tsv",
                                 "f": "suppnote6_f_retested_claims.tsv",
                                 "h": "suppnote6_h_modules_signed_axes_by_paper.tsv",
                                 "i": "suppnote6_i_gandal_copa_axes.tsv",
                                 "g": "suppnote6_g_snap_other_snrnaseq.tsv",
                                 "k": "panel_k_snap_best_match_chu_vs_sz.tsv",
                                 "l": "panel_l_snap_specificity.tsv"}.items()}
FILES["r2"] = ND / "module_retest" / "pooled_over_modules.tsv"
rd = lambda k: pd.read_csv(FILES[k], sep="\t", comment="#")


def panel_a(ax):
    # suppnote5_coexpr_null (2026-10-13): 62 lists of the 64, 1,891 pairs, 1,000 draws per null.
    # Third null keeps coexpression: each list's hits shifted along the GTEx cortex coexpression ordering within
    # 40 DE-prior bins, lists of one cohort cluster shifted together (R2). Whiskers: delete-one jackknife over cohort clusters.
    a = pd.read_csv(ND / "suppnote5_coexpr_null" / "suppnote5_a_coexpr_null_current.tsv", sep="\t", comment="#")
    levels = [("mean_pathway_jaccard", "GO pathways"), ("mean_gene_jaccard", "genes")]
    nulls = [("uniform_independent", "random genes", fk.ACCENT),
             ("de_prior_independent", "DE-prior-weighted random genes", fk.ACCENT_3),
             ("rot_R2_DEprior40_clusters_together", "DE-prior-matched, controls for gene correlations", fk.ACCENT_2)]
    w = 0.26
    for j, (nk, nl, col) in enumerate(nulls):
        for i, (mk, _) in enumerate(levels):
            r = a[(a["metric"] == mk) & (a["null"] == nk)].iloc[0]
            x = i + (j - 1) * w
            v, lo, hi = r["obs_over_null"], r["jack_cluster_lo"], r["jack_cluster_hi"]
            ax.bar(x, v, w * 0.95, color=col, label=nl if i == 0 else None)
            ax.errorbar(x, v, yerr=[[v - lo], [hi - v]], color=fk.TEXT, lw=0.6, capsize=1.5)
            ax.text(x, hi + 0.6, f"{v:.2f}", ha="center", va="bottom", fontsize=fk.SMALL, color=fk.TEXT, rotation=90)
    ax.axhline(1, color=fk.SUBTLE, lw=0.6, ls="--")
    ax.set_xticks(range(len(levels)))
    n, nl_ = int(a["n_pairs"].iloc[0]), int(a["n_lists"].iloc[0])
    ax.set_xticklabels([f"{l}\n({n:,} pairs of {nl_} lists)" for _, l in levels], fontsize=fk.SMALL)
    ax.set_ylim(0, 30)
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_ylabel("Agreement of two lists\n(mean Jaccard, observed over null)")
    ax.set_title("Pathways and genes shared between lists", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="upper right", title="null: lists replaced by",
              title_fontsize=fk.SMALL)
    ax.text(0.5, 1.3, "chance", fontsize=fk.SMALL, color=fk.SUBTLE, ha="center", va="bottom")
    ax.text(0.98, 0.42, "one bar = observed over one null\nwhisker = 95% interval,\ncohort clusters left out in turn",
            transform=ax.transAxes, ha="right", va="top", fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_b(ax):
    d = pd.read_csv(D / "suppnote6_d_pH_test_term_groups_ten_regions.tsv", sep="\t", comment="#").set_index("set")
    # gene counts in the labels are the genes entering each combined z (a pH z in at least 8 of the 10 dataset-regions)
    ng = d["n_genes_with_combined_z"]
    rows = [("oxphos_terms", f"oxidative phosphorylation,\nall term genes ({int(ng['oxphos_terms'])})"),
            ("oxphos_hits", f"oxidative phosphorylation,\nlist genes ({int(ng['oxphos_hits'])})"),
            ("myelin_terms", f"myelination,\nall term genes ({int(ng['myelin_terms'])})"),
            ("myelin_hits", f"myelination,\nlist genes ({int(ng['myelin_hits'])})"),
            ("microglia_terms", f"microglial activation,\nall term genes ({int(ng['microglia_terms'])})"),
            ("microglia_hits", f"microglial activation,\nlist genes ({int(ng['microglia_hits'])})"),
            ("recurrent_q10", f"recurrent genes ({int(ng['recurrent_q10'])})")]
    col = {"oxphos": fk.ACCENT_3, "myelin": fk.ACCENT, "microglia": fk.ACCENT_2, "recurrent": fk.NEUTRAL}
    y = np.arange(len(rows))
    for yi, (k, _) in zip(y, rows):
        p = d.loc[k]
        c = col[k.split("_")[0]]
        hatch = "////" if k.endswith("_hits") else None
        ax.barh(yi, p["combined_perm_z_ten"], 0.6, color=c if hatch is None else "white", edgecolor=c, hatch=hatch, lw=0.6)
        ax.plot([p["drop_one_min"], p["drop_one_max"]], [yi, yi], color=fk.TEXT, lw=0.9, solid_capstyle="butt")
        xt = max(p["combined_perm_z_ten"], p["drop_one_max"], 0) + 0.25
        ax.text(xt, yi, f"{p['combined_perm_z_ten']:.2f}", va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    for v in (-1.96, 1.96):
        ax.axvline(v, color=fk.SUBTLE, lw=0.5, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([r[1] for r in rows], fontsize=fk.SMALL, linespacing=0.85)
    ax.set_ylim(len(rows) + 2.2, -0.6)
    ax.set_xlim(-4.2, 8.0)
    ax.set_xlabel("Control-pH association (combined z)")
    ax.set_title("Shared terms and pH in control brains", loc="left")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=fk.NEUTRAL, label="bar: one gene set, ten dataset-regions"),
                       Line2D([], [], color=fk.TEXT, lw=0.9, label="range leaving out one region at a time"),
                       Patch(facecolor="white", edgecolor=fk.SUBTLE, hatch="////", label="only the terms' genes on the lists")],
              frameon=True, facecolor="white", edgecolor="none", framealpha=1, fontsize=fk.SMALL, loc="lower right", ncol=1)


def panel_c(ax):
    e = rd("e").set_index("sz_retest_class")
    f = rd("f")
    rows = [("held, independent donors", "held in other donors"),
            ("mixed or not held, independent", "mixed or not held,\nother donors"),
            ("held, partly shared donors", "held, partly\nshared donors"),
            ("same donors or cohort split", "same donors or\ncohort split"),
            ("structure only", "module structure only"),
            ("none offered", "none offered")]
    retest = {"held, independent donors", "mixed or not held, independent", "held, partly shared donors"}
    y = np.arange(len(rows))
    for yi, (k, _) in zip(y, rows):
        n = int(e.loc[k, "claims_all"])
        if k in retest:
            sub = f[f["sz_retest"] == k]
            other = int((sub["same_lineage_in_replication"] == 0).sum())
            ax.barh(yi, n - other, 0.6, color=fk.ACCENT)
            ax.barh(yi, other, 0.6, left=n - other, color=fk.ACCENT_2)
        else:
            ax.barh(yi, n, 0.6, color=fk.NEUTRAL)
        ax.text(n + 1.5, yi, str(n), va="center", fontsize=fk.SMALL, color=fk.TEXT)
    ax.set_yticks(y)
    ax.set_yticklabels([r[1] for r in rows], fontsize=fk.SMALL)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(0, 85)
    ax.set_xlabel(f"Claims (count of {int(e.loc['total', 'claims_all'])} claims, {80} papers)")
    ax.set_title("Coexpression claims, by re-test", loc="left")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=fk.ACCENT, label="re-test, same bank lineage"),
                       Patch(color=fk.ACCENT_2, label="re-test, other bank lineage"),
                       Patch(color=fk.NEUTRAL, label="no re-test in other donors")],
              frameon=False, fontsize=fk.SMALL, loc="upper right", title="one bar = one kind of claim",
              title_fontsize=fk.SMALL)


def panel_snap_match(ax):
    k = rd("k")
    k["best_abs_rho"] = pd.to_numeric(k["best_abs_rho"])
    chu = k[k["source"].str.startswith("Chu")]
    med = float(chu[chu["statistic"].str.startswith("median")]["best_abs_rho"].iloc[0])
    mx = float(chu[chu["statistic"].str.startswith("maximum")]["best_abs_rho"].iloc[0])
    own = k[~k["source"].str.startswith("Chu")].set_index("dataset")
    rows = [("Ling 2024 (Gemma counts through this pipeline)", "Ling 2024, discovery data\n(positive control)", fk.ACCENT_2),
            ("SZBDMulti-Seq minus Ling donors", "SZBDMulti-Seq minus Ling donors", fk.ACCENT),
            ("GSE254569", "GSE254569", fk.ACCENT),
            ("PsychAD Mount Sinai", "PsychAD Mount Sinai", fk.ACCENT),
            ("PsychAD HBCC", "PsychAD HBCC", fk.ACCENT),
            ("brainSCOPE CMC", "brainSCOPE CMC (sensitivity set)", fk.NEUTRAL)]
    y = np.arange(len(rows))
    for yi, (key, _, c) in zip(y, rows):
        r = own.loc[key]
        v = float(r["best_abs_rho"])
        ax.plot([0, v], [yi, yi], color=c, lw=0.8, zorder=2)
        ax.plot(v, yi, "o", color=c, ms=4.5, zorder=3)
        ax.text(v + 0.02, yi, f"{v:.2f} ({int(r['n_donors_in_peer'])} donors)", va="center", fontsize=fk.SMALL,
                color=fk.TEXT)
    ax.set_yticks(y)
    ax.set_yticklabels([r[1] for r in rows], fontsize=fk.SMALL)
    ax.set_ylim(len(rows) + 0.9, -1.0)
    ax.set_xlim(0, 1.3)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Best match with Ling's SNAP factor\n(absolute Spearman rho of gene loadings)")
    ax.set_title("SNAP re-derived in other data", loc="left")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.text(1.0, 1.0, "one point = one dataset", transform=ax.transAxes, ha="right", va="top",
            fontsize=fk.SMALL, color=fk.SUBTLE)


def panel_snap_effect(ax, ax2):
    g = rd("g")
    l = rd("l")
    ds = g[g["row_type"] == "dataset"].set_index("set")
    pl = g[g["row_type"] == "pooled"].set_index("dataset")
    dl = l[l["row_type"] == "dataset"].set_index("set")
    pdl = l[l["row_type"].str.startswith("pooled")].iloc[0]
    rows = [("ds", "Ling", "Ling 2024 (discovery)", fk.ACCENT_2),
            ("ds", "HBCC_main", "PsychAD HBCC", fk.ACCENT),
            ("ds", "MSSM_main", "PsychAD Mount Sinai", fk.ACCENT),
            ("ds", "MSSM_nodementia", "  without dementia-labelled\n  schizophrenia donors", fk.NEUTRAL),
            ("ds", "SZBD_noLing", "SZBDMulti-Seq\nminus Ling donors", fk.ACCENT),
            ("ds", "NSW", "GSE254569", fk.ACCENT),
            ("pool", "primary: SNAP total SZ", "pooled, four banks", fk.TEXT),
            ("pool", "(a) MSSM without dementia-labelled SZ: SNAP total SZ",
             "pooled, Mount Sinai without\ndementia-labelled donors", fk.NEUTRAL)]
    y = np.arange(len(rows))
    for yi, (kind, key, _, c) in zip(y, rows):
        if kind == "ds":
            r = ds.loc[key]
            est, lo, hi = (float(r[x]) for x in ("snap_sz_sd", "snap_sz_lo", "snap_sz_hi"))
            ax.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt="o", color=c, ms=3.5, lw=0.8, capsize=0)
            ax.text(-2.05, yi, f"{int(r['n_sz'])} / {int(r['n_ctrl'])}", va="center", ha="left",
                    fontsize=fk.TINY, color=fk.SUBTLE)
            q = dl.loc[key]
            d, dlo, dhi = (float(q[x]) for x in ("delta_on_minus_off", "delta_lo", "delta_hi"))
            ax2.errorbar(d, yi, xerr=[[d - dlo], [dhi - d]], fmt="o", color=c, ms=3.5, lw=0.8, capsize=0)
        else:
            r = pl.loc[key]
            est, lo, hi = (float(r[x]) for x in ("snap_sz_sd", "snap_sz_lo", "snap_sz_hi"))
            ax.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt="D", color=c, ms=4, lw=1.0, capsize=0,
                        mfc=c if c == fk.TEXT else "white")
            if key.startswith("primary"):
                d, dlo, dhi = (float(pdl[x]) for x in ("delta_on_minus_off", "delta_lo", "delta_hi"))
                ax2.errorbar(d, yi, xerr=[[d - dlo], [dhi - d]], fmt="D", color=c, ms=4, lw=1.0, capsize=0)
    for a in (ax, ax2):
        a.axvline(0, color=fk.TEXT, lw=0.6)
        a.axhline(5.5, color=fk.SUBTLE, lw=0.4)
        a.set_ylim(len(rows) - 0.4, -1.0)
        a.grid(axis="y", visible=False)
        a.grid(axis="x", visible=True)
    ax.text(-2.05, -0.75, "SZ / controls", fontsize=fk.TINY, color=fk.SUBTLE, ha="left", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([r[2] for r in rows], fontsize=fk.SMALL)
    ax.set_xlim(-2.1, 0.6)
    ax.set_xlabel("Schizophrenia effect\non SNAP score (SD)")
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.tick_params(axis="y", length=0)
    ax2.set_xlim(-0.75, 0.45)
    ax2.set_xlabel("SNAP cell types\nminus other\ncell types (SD)")
    ax.set_title("Schizophrenia effect on the SNAP score", loc="left")
    # (annotation removed: it covered the last row; the legend says what a point is)


def panel_modules(ax):
    h = rd("h")
    p = h[h["row_type"] == "paper"].copy()
    for c in ("n_modules", "median_signed_pH", "median_signed_agonal"):
        p[c] = pd.to_numeric(p[c])
    ax.axhline(0, color=fk.TEXT, lw=0.6)
    ax.axvline(0, color=fk.TEXT, lw=0.6)
    on = p["median_signed_pH"] > 0
    ax.scatter(p.loc[on, "median_signed_pH"], p.loc[on, "median_signed_agonal"], s=8 + 4 * p.loc[on, "n_modules"],
               color=fk.ACCENT_3, alpha=0.8, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(p.loc[~on, "median_signed_pH"], p.loc[~on, "median_signed_agonal"], s=8 + 4 * p.loc[~on, "n_modules"],
               color=fk.NEUTRAL, alpha=0.9, edgecolor="white", lw=0.4, zorder=3)
    for _, r in p[~on].iterrows():
        ax.annotate(r["paper"], (r["median_signed_pH"], r["median_signed_agonal"]), xytext=(5, -2),
                    textcoords="offset points", ha="left", fontsize=fk.TINY, color=fk.SUBTLE)
    g = p[p["paper"].str.startswith("Gandal")].iloc[0]
    ax.annotate("Gandal 2018", (g["median_signed_pH"], g["median_signed_agonal"]), xytext=(-6, 8),
                textcoords="offset points", ha="right", fontsize=fk.TINY, color=fk.SUBTLE)
    n_ph, n_ag, n_p = int(on.sum()), int((p["median_signed_agonal"] > 0).sum()), len(p)
    allrow = h[h["row_type"].str.startswith("all")].iloc[0]
    # counts the plot already shows are not repeated in it
    gi = rd("i")
    gz = gi[(gi["set"] == "Gandal modules") & (gi["comparison"] == "SCZ")].set_index("axis")
    ag, ph = gz.loc["mean_agonal_t"], gz.loc["mean_pH_Z"]
    ax.text(0.02, 0.98, "Gandal 2018, 34 modules:\nschizophrenia association against\n"
                        f"agonal axis rho {float(ag['rho']):+.2f} (p {float(ag['p']):.3f})\n"
                        f"pH axis rho {float(ph['rho']):+.2f} (p {float(ph['p']):.2f})\n\none point = one paper\nsize = number of modules",
            transform=ax.transAxes, ha="left", va="top", zorder=4,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5), fontsize=fk.TINY, color=fk.SUBTLE)
    ax.set_xlim(-2.6, 5.0)
    ax.set_ylim(-4.3, 6.4)
    ax.set_xlabel("pH axis (median signed control-pH z;\n> 0 = low-pH side)")
    ax.set_ylabel("Agonal axis (median signed GTEx\nslow-death t; > 0 = agonal side)")
    ax.set_title("Coexpression modules on the pH and agonal axes", loc="left")
    ax.grid(axis="x", visible=True)


def panel_retest(ax):
    pm = pd.read_csv(FILES["r2"], sep="\t", comment="#")
    pm = pm[(pm["lineage_named_only"] == "all") & (pm["min_clusters"] == 1)].set_index("model")
    cols = [("m1", "age, sex,\nregion"), ("m3", "+ pH and\nagonal scores")]
    segs = [("modules_interval_excludes_0_matching", "as reported", fk.ACCENT_2),
            (None, "unclear", fk.NEUTRAL),
            ("modules_interval_excludes_0_opposite", "opposite", fk.ACCENT_4)]
    for i, (m, _) in enumerate(cols):
        r = pm.loc[m]
        n = int(r["n_modules"])
        counts = [int(r[segs[0][0]]), None, int(r[segs[2][0]])]
        counts[1] = n - counts[0] - counts[2]
        bottom = 0
        for (k, lab, col), c in zip(segs, counts):
            ax.bar(i, c, 0.62, bottom=bottom, color=col, label=lab if i == 0 else None)
            if c >= 10:
                ax.text(i, bottom + c / 2, str(c), ha="center", va="center", fontsize=fk.SMALL, color="white" if col != fk.NEUTRAL else fk.TEXT)
            bottom += c
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], fontsize=fk.SMALL)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 196)
    ax.set_ylabel("Modules (of 196)")
    ax.set_title("Modules tested in other cohorts", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, title="direction in\nother cohorts", title_fontsize=fk.SMALL)
    ax.grid(axis="x", visible=False)


def panel_deans(ax):
    """Deans et al. 2026 convergent genes (363, pooled screen) in the project's patient-derived and peripheral lists."""
    o = pd.read_csv(DEANS, sep="\t", comment="#")
    lanes = [("ipsc", "iPSC-derived lists (106)", fk.ACCENT), ("periph", "blood, LCL and fibroblast lists (123)", fk.ACCENT_3)]
    for lane, lab, col in lanes:
        r = o[o["lane"] == lane].sort_values("k_lists")
        ax.plot(r["k_lists"], r["observed"], "-o", color=col, ms=3.5, lw=1.4, label=f"{lab}: observed")
        ax.plot(r["k_lists"], r["expected_deprior"], "--", color=col, lw=1.0, label="expected, DE prior")
        ax.plot(r["k_lists"], r["expected_uniform"], ":", color=col, lw=1.0, label="expected, random draws")
    ax.set_xlabel("Lists naming the gene (at least k)")
    ax.set_ylabel("Convergent genes (count of 350 tested)")
    ax.set_title("Genes changed by risk-gene perturbations, in patient-derived lists", loc="left")
    ax.legend(frameon=False, fontsize=fk.SMALL, ncol=1, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.text(0.01, 0.03, "one point = number of convergent genes in at least k lists", transform=ax.transAxes,
            fontsize=fk.SMALL, color=fk.SUBTLE)


DEANS = D / "panel_o_deans_convergent_genes_recurrence.tsv"


def main():
    fk.apply_style()
    fig = fk.figure(200)
    A = fig.add_axes([0.10, 0.745, 0.33, 0.195])
    B = fig.add_axes([0.68, 0.695, 0.30, 0.245])
    D = fig.add_axes([0.27, 0.405, 0.20, 0.215])
    E = fig.add_axes([0.70, 0.405, 0.12, 0.215])
    E2 = fig.add_axes([0.85, 0.405, 0.10, 0.215])
    F = fig.add_axes([0.20, 0.075, 0.40, 0.215])
    G = fig.add_axes([0.70, 0.075, 0.13, 0.215])
    panel_a(A); panel_b(B); panel_snap_match(D); panel_snap_effect(E, E2); panel_modules(F); panel_retest(G)
    FILES["a"] = ND / "suppnote5_coexpr_null" / "suppnote5_a_coexpr_null_current.tsv"
    for ax, l in zip((A, B, D, E, F, G), "abcdef"):
        fk.letter(ax, l)
    fk.save(fig, "supp_pathways_coexpression", "supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_pathways_coexpression.py",
            list(FILES.values()),
            notes="Panels: a pathways shared, b pH of pathway groups, c SNAP match, d SNAP effect, e modules on the axes, f modules tested in other cohorts.")


if __name__ == "__main__":
    main()
