#!/usr/bin/env python3
"""Supplementary Note 6 (coexpression), extra panels j to o for a figure that carries SNAP.

Request: RA task of 2026-10-09 (handoff RA_TO_MAB_2026_10_09_NOTE6_COEXPRESSION_COMPLETE.md).
Companion to export_suppnote1.py (panels a to i); that script and its files are not touched. Reads existing outputs only
(no model fits, no network). Never overwrites: a panel file that already exists is skipped.

    python3 analyses/figure_data_exports/export_suppnote6_coexpression_extra.py

Panels:
  j  SNAP rediscovery (Arm A): best-matching PEER factor per dataset, all runs of arm_a_summary.tsv
  k  best-match correlation with Ling's SNAP loadings: Chu et al. seven non-SZ datasets (summary statistics as printed in
     the preprint text; per-dataset values are only in its Supplemental Figure 3C) and this project's SZ datasets
  l  specificity: Ling's SNAP weights on SNAP cell types against other cell types (Arm B), per dataset and pooled
  m  PsychAD projection of Ling's SNAP, SNAP-a and SNAP-n: SZ and age effects, dementia-labelled sensitivity, cross-cell-type coordination
  n  lab's earlier SNAP-in-SZ analysis (Rui, 2025) beside this project's best |Pearson| with Ling's feature set
  o  Deans et al. 2026 convergent genes (363 pooled-screen meta-DEGs) recurrence in the project's iPSC-type and peripheral lists
"""
import csv
import hashlib
import os

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
REPO = os.path.dirname(AN)
RUN = "2026-10-09"
SCRIPT = "analyses/figure_data_exports/export_suppnote6_coexpression_extra.py"
FIG = "suppnote6_pathways_coexpression"
OUT = os.path.join(os.path.expanduser("figure_data"), FIG)
CO = os.path.join(REPO, "side-studies", "coexpression")
SNAP = os.path.join(CO, "snap_sz", "results")
PSY = os.path.join(CO, "snap_psychad", "snap_psychad_results.tsv")
DEANS = os.path.join(REPO, "side-studies", "deans-2026", "results")
CHU_PDF = os.path.join(REPO, "data", "raw", "doi_10.64898_2026.06.17.731965",
                       "Chu et al. - 2026 - Evaluation of analysis modes for RNA coexpression in single-cell and bulk tissue.pdf")
RUI = os.path.expanduser("Documents/manuscripts/eric/xCell vs xBulk/ling-case-study-files-rui/correlation_summary_scz.csv")
STATUS = os.path.join(OUT, "STATUS.md")

_sha = {}


def sha(p):
    p = os.path.abspath(p)
    if p not in _sha:
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for c in iter(lambda: fh.read(1 << 20), b""):
                h.update(c)
        _sha[p] = h.hexdigest()
    return _sha[p]


def rel(p):
    p = os.path.abspath(p)
    return os.path.relpath(p, REPO) if p.startswith(REPO) else p.replace(os.path.expanduser("~"), "~")


def rd(p, delim="\t"):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter=delim, quoting=csv.QUOTE_NONE))


def num(x):
    return float(x) if x not in ("", None) else None


def write(name, panel, fields, rows, inputs, comments=()):
    path = os.path.join(OUT, name)
    if os.path.exists(path):
        print("exists, skipped:", name)
        return None
    os.makedirs(OUT, exist_ok=True)
    lines = [f"# figure: {FIG} panel {panel}", f"# script: {SCRIPT}, run {RUN}"]
    lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs] + [f"# {c}" for c in comments]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for l in lines:
            fh.write(l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for r in rows:
            w.writerow(["" if v is None else v for v in r])
    print("wrote", name, len(rows), "rows")
    return (name, panel, len(rows))


def panel_j():
    p = os.path.join(SNAP, "arm_a_summary.tsv")
    cols = ["run", "set", "route", "k", "classes", "features", "donors", "best_factor", "best_rho", "rho_astrocyte",
            "rho_glutamatergic", "best_abs_pearson", "call", "n_sz", "n_ctrl", "factor_sz", "factor_sz_t", "factor_sz_p",
            "factor_age_ctrl_per_decade"]
    rows = [[r[c] for c in cols] for r in rd(p)]
    return write("panel_j_snap_arm_a_rediscovery.tsv", "j", cols, rows, [p],
                 ["all 57 completed runs of Arm A (PEER on donor pseudobulk; routes resid_top50 = primary, resid_lingset, logcpm_top50, nemo); k = number of PEER factors",
                  "best_rho = Spearman of the best-matching factor's loadings with Ling's LF4 (sign arbitrary); call = pre-set call from PLAN.md (found / partial / not found)",
                  "primary rows for the figure: route resid_top50, k 10; Ling 'found' rows are the positive control"])


def panel_k():
    a = os.path.join(SNAP, "arm_a_summary.tsv")
    fields = ["source", "dataset", "n_donors_in_peer", "sz_donors_present", "best_abs_rho", "statistic", "note"]
    rows = [["Chu et al. 2026 (preprint)", "7 other single-nucleus datasets (Lau, Lim, Nagy/Pineda, Pineda, Ramos, ROSMAP, Velmeshev)", None, "no",
             0.15, "median best match over datasets", "Results, 'Case study': median best match 0.15 over 210 comparisons"],
            ["Chu et al. 2026 (preprint)", "7 other single-nucleus datasets", None, "no",
             0.33, "maximum best match", "highest correlation 0.33, LF3 in the Lim data set; per-dataset values are in Supplemental Figure 3C only (figure, no table on disk)"]]
    keep = {"HBCC_main": "PsychAD HBCC", "MSSM_main": "PsychAD Mount Sinai", "SZBD_noLing": "SZBDMulti-Seq minus Ling donors",
            "NSW": "GSE254569", "CMC": "brainSCOPE CMC", "Ling": "Ling 2024 (Gemma counts through this pipeline)"}
    for r in rd(a):
        if r["route"] == "resid_top50" and r["k"] == "10" and r["set"] in keep:
            rows.append(["this project, Arm A (residuals, top 50% genes, K = 10)", keep[r["set"]], r["donors"], "yes" if r["set"] != "Ling" else "yes (positive control)",
                         abs(float(r["best_rho"])), "best |Spearman| with Ling LF4 loadings", "pre-set call: " + r["call"]])
    return write("panel_k_snap_best_match_chu_vs_sz.tsv", "k", fields, rows, [a, CHU_PDF],
                 ["Chu rows are typed from the preprint text (PDF page text, sha above); they are not read from a table",
                  "Chu et al. is a bioRxiv preprint, not peer reviewed; its seven datasets have no schizophrenia donors"])


def panel_l():
    d = os.path.join(SNAP, "arm_b_delta.tsv")
    m = os.path.join(SNAP, "snap_sz_meta.tsv")
    fields = ["row_type", "set", "n_sz", "n_ctrl", "sz_effect_on_snap_cell_types_sd", "sz_effect_on_other_cell_types_sd", "delta_on_minus_off",
              "delta_lo", "delta_hi", "off_classes"]
    rows = [["dataset", r["set"], r["n_sz"], r["n_ctrl"], r["B_on"], r["B_off"], r["delta"], r["delta_lo"], r["delta_hi"], r["off_classes"]]
            for r in rd(d)]
    for r in rd(m):
        if r["pool"] == "primary" and r["quantity"].startswith("Delta"):
            rows.append(["pooled, primary pool of four banks (random effects, Hartung-Knapp interval)", r["sets"], None, None, None, None,
                         r["estimate"], r["lo_hksj"], r["hi_hksj"], ""])
    return write("panel_l_snap_specificity.tsv", "l", fields, rows, [d, m],
                 ["same Ling LF4 astrocyte and glutamatergic weights applied gene by gene to microglia, endothelia and oligodendrocyte expression; interval = donor bootstrap (1,000) for datasets"])


def panel_m():
    cols = None
    rows = []
    with open(PSY, encoding="utf-8") as fh:
        rr = csv.reader((l for l in fh if not l.startswith("#")), delimiter="\t", quoting=csv.QUOTE_NONE)
        cols = next(rr)
        rows = [r for r in rr]
    return write("panel_m_snap_psychad_projection.tsv", "m", cols, rows, [PSY],
                 ["Ling's LF4, SNAP-a and SNAP-n loadings projected onto PsychAD pseudobulk; null = loadings permuted within cell type (1,000); null_median is not zero by construction",
                  "set Mount-Sinai-without-dementia leaves out the SZ donors with a dementia term; 'concerted rho' rows are the astrocyte-neuron correlation across donors"])


def panel_n():
    a = os.path.join(SNAP, "arm_a_summary.tsv")
    fields = ["dataset_rui", "rui_best_abs_corr_k10", "rui_best_variable", "this_project_set", "this_project_donors", "this_project_best_abs_pearson_ling_feature_set_k10"]
    mp = {"CMC": ["CMC"], "SZBDMulti-Seq": ["SZBD_all", "SZBD_noLing"], "Frolich": ["NSW"], "Batiuk": ["Batiuk"], "MultiomeBrain": ["Multiome"]}
    ours = {}
    for r in rd(a):
        if r["route"] == "resid_lingset" and r["k"] == "10":
            ours[r["set"]] = r
    rows = []
    with open(RUI, encoding="utf-8") as fh:
        rui = list(csv.DictReader(fh))
    for r in rui:
        for s in mp.get(r["dataset"], []):
            o = ours.get(s)
            rows.append([r["dataset"], r["max_abs_corr"], r["max_corr_var"], s, o["donors"] if o else None, o["best_abs_pearson"] if o else None])
    return write("panel_n_snap_lab_earlier_analysis.tsv", "n", fields, rows, [RUI, a],
                 ["Rui Xiang Yu, 2025: SCTransform residuals, all donors, Ling's gene_class features, Pearson correlation of loadings (correlation_summary_scz.csv, max_abs_corr = best of Var1-Var10)",
                  "this project: analytic Pearson residuals, donors with at least 20 nuclei in all six classes, Ling's feature set, K = 10, route resid_lingset",
                  "Rui's PsychAD run does not exist; her slide deck rui-sept5th-2025.pdf (slides 15, 24, 37) is not tabulated here"])


def panel_o():
    r = os.path.join(DEANS, "recurrence_by_lane.tsv")
    d = os.path.join(DEANS, "direction_agreement.tsv")
    dd = rd(d)[0]
    rr = rd(r)
    fields = list(rr[0].keys())
    rows = [[x[c] for c in fields] for x in rr]
    return write("panel_o_deans_convergent_genes_recurrence.tsv", "o", fields, rows, [r, d],
                 ["Deans et al. 2026 (Nat Commun) pooled-screen meta-analysed DEGs, 363 genes (sheet A1 of the authors' gene-list workbook), in lists from the project's iPSC-type (ipsc) and peripheral (periph) studies",
                  "observed = number of the tested convergent genes in at least k lists; expected under the uniform and DE-prior nulls of analyses/recurrence.py; side-studies/deans-2026/README.md, files dated 2026-10-03",
                  f"direction agreement with Deans' direction over genes with a direction in these lists: {dd['same']} same, {dd['opposite']} opposite of {dd['total']} ({dd['pct_same']}%)",
                  "the authors sent corrected meta-analysis supplementary files in January 2026 (side_projects_inventory_2026-10-02.md); whether the workbook used here predates the correction is not checked"])


def main():
    done = [x for x in (panel_j(), panel_k(), panel_l(), panel_m(), panel_n(), panel_o()) if x]
    if done and not any("Appended 2026-10-09: panels j to o" in l for l in open(STATUS, encoding="utf-8")):
        with open(STATUS, "a", encoding="utf-8") as fh:
            fh.write(f"\n\n## Appended {RUN}: panels j to o (extra coexpression panels for a figure with SNAP)\n\n"
                     f"Script: {SCRIPT} (new; export_suppnote1.py untouched). Reads existing outputs only; each TSV has '#' provenance lines. "
                     "Files are not overwritten. Added by the RA task of 2026-10-09.\n\n")
            names = {"j": "SNAP rediscovery (Arm A), all runs; source snap_sz/results/arm_a_summary.tsv (built 2026-09-28)",
                     "k": "best-match correlation with Ling's SNAP: Chu et al. summary (typed from the preprint text) and this project's SZ datasets",
                     "l": "specificity of SNAP weights, on-target against other cell types, per dataset; source arm_b_delta.tsv (2026-09-28)",
                     "m": "PsychAD projection of SNAP, SNAP-a, SNAP-n, dementia sensitivity, cross-cell-type coordination; source snap_psychad_results.tsv (2026-09-27)",
                     "n": "lab's earlier SNAP-in-SZ analysis (Rui 2025) beside this project's best |Pearson|; source correlation_summary_scz.csv (2026-01-13) and arm_a_summary.tsv",
                     "o": "Deans et al. 2026 convergent genes recurrence in iPSC-type and peripheral study lists; source side-studies/deans-2026/results/ (2026-10-03); optional, not part of the coexpression claim chain"}
            for name, panel, n in done:
                fh.write(f"- **{panel}** `{name}` ({n} rows): {names[panel]}.\n")
            fh.write("\nNOT exported (no table on disk): Chu et al. per-dataset best-match correlations (Supplemental Figure 3C of the preprint, figure only); "
                     "Ling's own discovery statistics beyond those recomputed in snap_psychad (Ling's published donor scores are not in the repo as a table); "
                     "SNAP lists in the DE recurrence with and without Ling (forking_paths_rerun_2026-10-09, not tabulated for this note).\n")


if __name__ == "__main__":
    main()
