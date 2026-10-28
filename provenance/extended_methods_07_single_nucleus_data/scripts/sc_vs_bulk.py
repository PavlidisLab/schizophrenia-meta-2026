#!/usr/bin/env python3
"""Single-nucleus agreement beside bulk on the same axes: summary tables and figures for single_cell_report.md.

    .claude/skills/plotting/.venv/bin/python single_cell/sc_vs_bulk.py      (light: reads result tables only)

Inputs (nothing hardcoded; every number comes from these files):
  single_cell/results/sc_split_half_reps.tsv.gz   PsychAD per-class split halves, SZ / age / sex (scripts/sc_split_half.R, remote)
  single_cell/results/sc_cross_pairs.tsv           per class x source pair, full n (scripts/sc_cross_cohort.R, remote)
  single_cell/results/sc_cohort_tier.tsv           cohort-tier recurrence per class (sc_cohort_tier.py)
  libd_brainseq/analysis/results/a_summary.tsv, d_summary.tsv     bulk LIBD split halves (SZ; SZ, age, sex)
  fairness_sample_level/results/a_summary.tsv      bulk Gemma datasets: split halves and subject-separate pairs at matched n
  fairness_sample_level/results/b_pairs.tsv        bulk Gemma dataset pairs at full n (level L2: sz + age + sex + PMI)
  fairness_sample_level/results/datasets_used.tsv  n per group of those datasets
Outputs: single_cell/results/sc_split_half_summary.tsv, sc_vs_bulk.tsv, figures/fig_*.svg/.png, figures/CAPTIONS.md,
         sc_vs_bulk_meta.json
"""
import datetime
import hashlib
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.expanduser(".claude/skills/plotting/python"))
from pavlab_style import apply_pavlab_style  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.path.dirname(HERE)
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
ACCENT, ACCENT_2, ACCENT_3, ACCENT_4 = "#2563eb", "#10b981", "#f59e0b", "#ef4444"
TEXT, SUBTLE, GRID = "#1f2937", "#6b7280", "#e5e7eb"
CLASS_LABEL = {"total": "whole tissue", "EN": "excitatory", "IN": "inhibitory", "Astro": "astrocyte", "OPC": "OPC",
               "Oligo": "oligodendrocyte", "Immune": "immune", "Mural": "mural", "Endo": "endothelial"}
CAPTIONS = {}
INPUTS = []


def rd(rel, **kw):
    p = os.path.join(GL, rel)
    INPUTS.append(rel)
    if p.endswith(".gz"):
        return pd.read_csv(p, sep="\t", comment="#", compression="gzip", **kw)
    return pd.read_csv(p, sep="\t", comment="#", **kw)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def style():
    apply_pavlab_style()
    plt.rcParams.update({"axes.titleweight": "normal", "axes.titlesize": 11, "axes.titlelocation": "left", "font.size": 10,
                         "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "svg.fonttype": "none",
                         "font.family": ["Helvetica", "Arial", "sans-serif"]})


def save(fig, name, source, caption, y=-0.02):
    fig.text(0.01, y, source, fontsize=8, color=SUBTLE, ha="left", va="top")
    for ax in fig.axes:
        ax.set_clip_on(False)
        for a in list(ax.patches) + list(ax.lines) + list(ax.texts) + list(ax.collections):
            a.set_clip_on(False)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".svg"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    CAPTIONS[name] = caption


def q(v, p):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return float(np.quantile(v, p)) if len(v) else np.nan


def main():
    style()
    # ---------------- sc split halves, summarized
    S = rd("single_cell/results/sc_split_half_reps.tsv.gz")
    S = S[S.gene_set == "auto"]
    g = S.groupby(["cell_type", "contrast", "comparison", "n_per_group"])
    SS = g.agg(reps=("rho_t", "size"), rho_t=("rho_t", "median"), rho_t_q25=("rho_t", lambda v: q(v, .25)),
               rho_t_q75=("rho_t", lambda v: q(v, .75)), rho_logfc=("rho_logfc", "median"),
               top100_sign_in_other=("top100_sign_in_other", "median"), top200_fold=("top200_fold", "median"),
               fdr_overlap=("fdr_overlap", "median"), frac_reps_fdr_any=("fdr_n_a", lambda v: float(np.mean(np.asarray(v) > 0)))).reset_index()
    with open(os.path.join(RES, "sc_split_half_summary.tsv"), "w") as f:
        f.write(f"# built {datetime.date.today().isoformat()} by single_cell/sc_vs_bulk.py from results/sc_split_half_reps.tsv.gz; gene_set auto "
                "(not on X or Y); medians and quartiles over replicates (HBCC_vs_*: two pairs per replicate)\n")
        SS.round(4).to_csv(f, sep="\t", index=False)

    # ---------------- bulk references
    LA = rd("libd_brainseq/analysis/results/a_summary.tsv")
    LA = LA[(LA.kind == "split_halves") & (LA.level != "max")].copy()
    LA["n"] = LA.n_per_group.astype(int)
    LD = rd("libd_brainseq/analysis/results/d_summary.tsv")
    LD = LD[(LD.kind == "split_halves") & (LD.gene_set == "auto")].copy()
    FA = rd("fairness_sample_level/results/a_summary.tsv")
    FA = FA[FA.level != "max"].copy()
    FA["n"] = FA.level.astype(int)
    BP = rd("fairness_sample_level/results/b_pairs.tsv")
    DU = rd("fairness_sample_level/results/datasets_used.tsv")
    nmin = dict(zip(DU.dataset, DU[["n_sz", "n_ctl"]].min(axis=1)))
    BP = BP[(BP.level == "L2") & (BP.both_bulk.astype(str) == "True")].copy()
    BP["n_min"] = [min(nmin.get(a, np.nan), nmin.get(b, np.nan)) for a, b in zip(BP.dataset_a, BP.dataset_b)]
    CP = rd("single_cell/results/sc_cross_pairs.tsv")
    CP["n_min"] = CP[["n_sz_a", "n_ctl_a", "n_sz_b", "n_ctl_b"]].min(axis=1)
    CT = rd("single_cell/results/sc_cohort_tier.tsv")

    # ---------------- sc_vs_bulk.tsv: one row per (domain, comparison, relation, cell class, n)
    rows = []
    for _, r in LA.iterrows():
        rows.append(dict(domain="bulk", data=f"LIBD BrainSeq {r.region}", comparison="split halves", contrast="sz", relation="split halves",
                         cell_class="bulk tissue", n_per_group=r.n, units=1, rho_t=r.rho_t, rho_t_q25=r.rho_t_q25, rho_t_q75=r.rho_t_q75,
                         source="libd_brainseq/analysis/results/a_summary.tsv"))
    for _, r in LD.iterrows():
        rows.append(dict(domain="bulk", data=f"LIBD BrainSeq {r.region}", comparison="split halves", contrast=r.contrast, relation="split halves",
                         cell_class="bulk tissue", n_per_group=int(r.n_per_group), units=1, rho_t=r.rho_t, rho_t_q25=r.rho_t_q25,
                         rho_t_q75=r.rho_t_q75, source="libd_brainseq/analysis/results/d_summary.tsv"))
    for _, r in FA[FA.relation.isin(["split_half", "subject-separate"])].iterrows():
        rows.append(dict(domain="bulk", data="Gemma bulk datasets (median over units)", comparison="split halves" if r.relation == "split_half" else "cross-dataset",
                         contrast="sz", relation=r.relation, cell_class="bulk tissue", n_per_group=r.n, units=int(r.n_units), rho_t=r.rho_t,
                         source="fairness_sample_level/results/a_summary.tsv"))
    for (ct, con, cmp, n), r in SS.set_index(["cell_type", "contrast", "comparison", "n_per_group"]).iterrows():
        rows.append(dict(domain="single nucleus", data="PsychAD pseudobulk", comparison="split halves" if cmp.startswith("split") else "cross-bank",
                         contrast=con, relation=cmp, cell_class=CLASS_LABEL.get(ct, ct), n_per_group=n, units=int(r.reps), rho_t=r.rho_t,
                         rho_t_q25=r.rho_t_q25, rho_t_q75=r.rho_t_q75, source="single_cell/results/sc_split_half_summary.tsv"))
    for rel, d in BP.groupby("relation"):
        rows.append(dict(domain="bulk", data="Gemma bulk dataset pairs, full n (L2)", comparison="full-n pairs", contrast="sz", relation=rel,
                         cell_class="bulk tissue", n_per_group=float(d.n_min.median()), units=len(d), rho_t=d.rho_t.median(),
                         rho_t_q25=q(d.rho_t, .25), rho_t_q75=q(d.rho_t, .75), source="fairness_sample_level/results/b_pairs.tsv"))
    for (grp, rel, cls), d in CP.assign(grp=CP.group_a + "/" + CP.group_b).groupby(["grp", "relation", "gemma_class"]):
        rows.append(dict(domain="single nucleus", data=f"source pairs {grp}", comparison="full-n pairs", contrast="sz", relation=rel,
                         cell_class=cls, n_per_group=float(d.n_min.median()) if d.n_min.notna().any() else np.nan, units=len(d),
                         rho_t=d.rho_t.median(), rho_t_q25=q(d.rho_t, .25), rho_t_q75=q(d.rho_t, .75),
                         source="single_cell/results/sc_cross_pairs.tsv"))
    VB = pd.DataFrame(rows)
    with open(os.path.join(RES, "sc_vs_bulk.tsv"), "w") as f:
        f.write(f"# built {datetime.date.today().isoformat()} by single_cell/sc_vs_bulk.py; rho_t = Spearman rho of t (agree(), "
                "fairness_sample_level/scripts/common.R) between two SZ vs control fits; medians (and quartiles where given) over units / replicates\n")
        VB.round(4).to_csv(f, sep="\t", index=False)

    # ---------------- Figure 1: agreement against n, split halves and cross-cohort, SZ beside age and sex
    cls_order = ["total", "EN", "IN", "Astro", "Oligo", "OPC", "Immune", "Mural", "Endo"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    panels = [("sz", "split", "a  SZ, split halves within a bank"), ("sz", "cross", "b  SZ, HBCC vs Mount Sinai at the same n"),
              ("age", "split", "c  Age and sex beside SZ, HBCC vs Mount Sinai")]
    for ax, (con, kind, title) in zip(axes, panels):
        ax.axhline(0, color=GRID, lw=0.8)
        if kind == "split" and con == "sz":
            for cmp, col in (("split_HBCC_main", ACCENT), ("split_MSSM_main", ACCENT_3)):
                for ct in cls_order:
                    d = SS[(SS.contrast == "sz") & (SS.comparison == cmp) & (SS.cell_type == ct)].sort_values("n_per_group")
                    if len(d):
                        ax.plot(d.n_per_group, d.rho_t, "-", color=col, lw=1.6 if ct == "total" else 0.7, alpha=1 if ct == "total" else 0.55)
            d = LA[LA.region == "DLPFC"].sort_values("n")
            ax.plot(d.n, d.rho_t, "-o", color=TEXT, ms=3, lw=1.6)
            d = FA[FA.relation == "split_half"].sort_values("n")
            ax.plot(d.n, d.rho_t, "--o", color=SUBTLE, ms=3, lw=1.2)
            ax.text(0.02, 0.97, "bulk LIBD DLPFC", color=TEXT, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.91, "bulk Gemma datasets (median)", color=SUBTLE, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.85, "PsychAD HBCC, per class", color=ACCENT, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.79, "PsychAD Mount Sinai, per class", color=ACCENT_3, transform=ax.transAxes, va="top", fontsize=8)
        elif kind == "cross":
            for cmp, col in (("HBCC_vs_MSSM_main", ACCENT_4), ("HBCC_vs_MSSM_nodementia", ACCENT_2)):
                for ct in cls_order:
                    d = SS[(SS.contrast == "sz") & (SS.comparison == cmp) & (SS.cell_type == ct)].sort_values("n_per_group")
                    if len(d):
                        ax.plot(d.n_per_group, d.rho_t, "-", color=col, lw=1.6 if ct == "total" else 0.7, alpha=1 if ct == "total" else 0.55)
            d = FA[FA.relation == "subject-separate"].sort_values("n")
            ax.plot(d.n, d.rho_t, "--o", color=SUBTLE, ms=3, lw=1.2)
            ax.text(0.02, 0.97, "bulk subject-separate dataset pairs (median)", color=SUBTLE, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.91, "HBCC vs Mount Sinai, per class", color=ACCENT_4, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.85, "HBCC vs Mount Sinai without dementia-term SZ", color=ACCENT_2, transform=ax.transAxes, va="top", fontsize=8)
        else:
            for con2, col in (("age", ACCENT), ("male", ACCENT_3), ("sz", SUBTLE)):
                for ct in cls_order:
                    d = SS[(SS.contrast == con2) & (SS.comparison == "HBCC_vs_MSSM_main") & (SS.cell_type == ct) & (SS.n_per_group <= 30)].sort_values("n_per_group")
                    if len(d):
                        ax.plot(d.n_per_group, d.rho_t, "-o", color=col, ms=2.5, lw=1.6 if ct == "total" else 0.7, alpha=1 if ct == "total" else 0.55)
            ax.text(0.02, 0.97, "age", color=ACCENT, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.91, "sex (autosomal genes)", color=ACCENT_3, transform=ax.transAxes, va="top", fontsize=8)
            ax.text(0.02, 0.85, "SZ", color=SUBTLE, transform=ax.transAxes, va="top", fontsize=8)
        ax.set_title(title)
        ax.set_xlabel("donors per group in each half")
    axes[0].set_ylabel("Spearman rho of t between the two fits")
    save(fig, "fig_1_agreement_vs_n",
         "Source: single_cell/results/sc_split_half_summary.tsv; libd_brainseq/analysis/results/a_summary.tsv, d_summary.tsv; "
         "fairness_sample_level/results/a_summary.tsv (single_cell/sc_vs_bulk.py)",
         "Agreement between two SZ vs control fits (Spearman rho of the t statistics, median over replicates) against the number of "
         "donors per group in each fit. Thick line: whole-tissue pseudobulk; thin lines: the eight PsychAD cell classes. (a) Split halves "
         "within one PsychAD bank, beside bulk split halves (LIBD DLPFC; median over Gemma bulk datasets). (b) One half from HBCC against "
         "one from Mount Sinai at the same n, beside bulk dataset pairs with separate subjects. (c) HBCC against Mount Sinai at the same n for age "
         "and sex (autosomal genes), with SZ in grey.")

    # ---------------- Figure 2: full-n cohort pairs, against n and per Gemma class, beside bulk dataset pairs
    E = CP[(CP.group_a == "endresz") & (CP.group_b == "endresz")].copy()
    Epair = E.groupby(["source_a", "source_b", "relation"]).agg(rho_t=("rho_t", "median"), n_min=("n_min", "first")).reset_index()
    bulk = BP[BP.relation == "subject-separate"].copy()
    bins = [(0, 10), (11, 30), (31, 1000)]
    for lo, hi in bins:
        b_ = bulk[(bulk.n_min >= lo) & (bulk.n_min <= hi)]
        s_ = Epair[(Epair.relation == "subject-separate") & (Epair.n_min >= lo) & (Epair.n_min <= hi)]
        for dom, d in (("bulk", b_), ("single nucleus", s_)):
            VB.loc[len(VB)] = dict(domain=dom, data="Gemma bulk dataset pairs (L2)" if dom == "bulk" else "Endresz cohort pairs (median over classes)",
                                   comparison=f"full-n pairs, smaller group {lo}-{hi if hi < 1000 else 'max'}", contrast="sz",
                                   relation="subject-separate", cell_class="all", n_per_group=float(d.n_min.median()) if len(d) else np.nan,
                                   units=len(d), rho_t=d.rho_t.median() if len(d) else np.nan, rho_t_q25=q(d.rho_t, .25), rho_t_q75=q(d.rho_t, .75),
                                   source="fairness_sample_level/results/b_pairs.tsv" if dom == "bulk" else "single_cell/results/sc_cross_pairs.tsv")
    with open(os.path.join(RES, "sc_vs_bulk.tsv"), "w") as f:
        f.write(f"# built {datetime.date.today().isoformat()} by single_cell/sc_vs_bulk.py; rho_t = Spearman rho of t (agree(), "
                "fairness_sample_level/scripts/common.R) between two SZ vs control fits; medians (and quartiles where given) over units / replicates\n")
        VB.round(4).to_csv(f, sep="\t", index=False)
    order = (E[E.relation == "subject-separate"].groupby("gemma_class").rho_t.median().sort_values(ascending=False).index.tolist())
    relcol = {"subject-separate": ACCENT, "partial subject overlap": ACCENT_3, "same bank, subject overlap unknown": ACCENT_2}
    fig, (ax0, ax) = plt.subplots(1, 2, figsize=(14, 4.4), gridspec_kw=dict(width_ratios=[1, 2.2]))
    ax0.axhline(0, color=GRID, lw=0.8)
    ax0.scatter(bulk.n_min, bulk.rho_t, s=10, color=SUBTLE, alpha=0.6, lw=0)
    for rel, col in relcol.items():
        d = Epair[Epair.relation == rel]
        ax0.scatter(d.n_min, d.rho_t, s=28, color=col, lw=0)
    ax0.set_xscale("log")
    ax0.set_xticks([5, 10, 20, 50])
    ax0.set_xticklabels(["5", "10", "20", "50"])
    ax0.minorticks_off()
    ax0.set_xlabel("smaller group size in the pair (log scale)")
    ax0.set_ylabel("Spearman rho of t")
    ax0.set_title("a  Full-data pairs against n")
    ax0.text(0.02, 0.97, f"grey: bulk dataset pairs, separate subjects (n={len(bulk)})", color=SUBTLE, transform=ax0.transAxes, va="top", fontsize=8)
    ax0.text(0.02, 0.91, f"colour: Endresz cohort pairs, median over classes (n={len(Epair)})", color=ACCENT, transform=ax0.transAxes, va="top", fontsize=8)
    ax.axhline(0, color=GRID, lw=0.8)
    rng = np.random.default_rng(20260927)
    ax.scatter(rng.uniform(-0.25, 0.25, len(bulk)), bulk.rho_t, s=6, color=SUBTLE, alpha=0.5, lw=0)
    ax.plot([-0.3, 0.3], [bulk.rho_t.median()] * 2, color=TEXT, lw=1.5)
    for i, c in enumerate(order, start=1):
        d = E[E.gemma_class == c]
        for rel, col in relcol.items():
            dd = d[d.relation == rel]
            ax.scatter(i + rng.uniform(-0.25, 0.25, len(dd)), dd.rho_t, s=14, color=col, lw=0, alpha=0.9)
        dsep = d[d.relation == "subject-separate"]
        if len(dsep):
            ax.plot([i - 0.3, i + 0.3], [dsep.rho_t.median()] * 2, color=TEXT, lw=1.5)
    ax.set_xticks(range(len(order) + 1))
    ax.set_xticklabels([f"bulk pairs"] + order, rotation=45, ha="right")
    ax.set_title("b  Endresz cohort pairs per Gemma cell class")
    for k, (rel, col) in enumerate(relcol.items()):
        ax.text(0.99, 0.97 - 0.06 * k, rel, color=col, transform=ax.transAxes, ha="right", va="top", fontsize=8)
    save(fig, "fig_2_cohort_pairs",
         "Source: single_cell/results/sc_cross_pairs.tsv; fairness_sample_level/results/b_pairs.tsv (level L2), datasets_used.tsv (single_cell/sc_vs_bulk.py)",
         "Agreement between full-data SZ vs control results of two donor sets (Spearman rho of t). (a) Against the smaller group size of "
         "the pair: grey, pairs of bulk Gemma datasets with separate subjects (model sz + age + sex + PMI); coloured, pairs of Endresz 2026 "
         "cohorts (one annotation and pipeline), median over the Gemma cell classes, coloured by donor relation. (b) The Endresz pairs per "
         "Gemma cell class (each point one cohort pair, median over the subclasses of the class); black bars, medians of the subject-separate "
         "pairs; left, the bulk pairs of (a).", y=-0.12)

    # ---------------- Figure 3: cohort-tier recurrence per class (Endresz cohorts), observed vs expected
    C = CT[(CT.analysis == "endresz") & (CT.tier == "cohort") & CT.obs_ge2.notna()].copy()
    C = C.sort_values("exp_ge2", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 3.8))
    x = np.arange(len(C))
    w = 0.27
    for k, (col, lab, colr) in enumerate((("obs_ge2", "observed", TEXT), ("exp_ge2", "expected, uniform", ACCENT), ("exp_prior_ge2", "expected, DE prior", ACCENT_3))):
        ax.bar(x + (k - 1) * w, C[col].astype(float), w, color=colr, label=lab)
    for xi, (o, nc) in enumerate(zip(C.obs_ge2, C.n_clusters)):
        ax.text(xi - w, float(o) + 1, f"{int(o)}\n({int(nc)} sets)", ha="center", va="bottom", fontsize=7, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(C.gemma_class, rotation=45, ha="right")
    ax.set_ylabel("genes in 2 or more donor sets")
    ax.set_title("Recurrence across independent donor sets (Endresz cohorts, lists at q < 0.1)")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "fig_3_cohort_tier",
         "Source: single_cell/results/sc_cohort_tier.tsv (single_cell/sc_cohort_tier.py; sc_vs_bulk.py)",
         "Genes at q < 0.1 in two or more independent donor sets (Mount Sinai cohorts merged, HBCC and MultiomeBrain merged), per Gemma "
         "cell class with at least two donor sets carrying DE genes; observed beside the expectation under the uniform null and under the "
         "DE prior. Numbers above the bars: observed count and donor sets with DE genes in that class.", y=-0.16)

    with open(os.path.join(FIG, "CAPTIONS.md"), "w") as f:
        f.write("# Captions: single_cell figures\n\nBuilt by `single_cell/sc_vs_bulk.py`; every number is read from the result tables named in each figure's source line.\n\n")
        for k, v in CAPTIONS.items():
            f.write(f"## {k}\n\n{v}\n\n")
    idx = "".join(f"<h3>{k}</h3><img src='{k}.svg' style='max-width:100%'><p>{v}</p>" for k, v in CAPTIONS.items())
    open(os.path.join(FIG, "index.html"), "w").write(f"<!doctype html><meta charset=utf-8><title>single_cell figures</title><body style='font-family:Arial;max-width:1100px'>{idx}</body>")
    json.dump(dict(script="single_cell/sc_vs_bulk.py", built=datetime.datetime.now().isoformat(timespec="seconds"),
                   inputs={p: sha(os.path.join(GL, p)) for p in sorted(set(INPUTS))},
                   outputs={p: sha(os.path.join(HERE, p)) for p in ["results/sc_split_half_summary.tsv", "results/sc_vs_bulk.tsv"] +
                            [f"figures/{k}.svg" for k in CAPTIONS]}),
              open(os.path.join(HERE, "sc_vs_bulk_meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
