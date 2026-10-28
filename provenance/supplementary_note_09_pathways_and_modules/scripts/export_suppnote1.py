#!/usr/bin/env python3
"""Supplementary Note 6 (pathways and coexpression): tidy panel exports.

Request: MAB_TO_RA_2026_10_09_DISCUSSION_SUPPLEMENTARY_NOTES.md, item 1.
Population: the 64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrence.tsv sha256
d11761227abd; recurrent set = the 176 study-tier genes at q < 0.1 (uniform null).

Reads existing analysis outputs only (no model fits, no network). Exception: panel b reads the theme-agreement null from
figure_data/suppnote6_pathways_coexpression/recomputed/, produced on 2026-10-09 by running the unchanged
analyses/uniform_vs_published_list_agreement/a5a_themes.py (1,000 permutations, seed 20260925) in a scratch copy of the tree, because the
repository copy of that output predates the 2026-10-08 additions to corpus/themes.json.

    python3 analyses/figure_data_exports/export_suppnote1.py
"""
import csv
import hashlib
import json
import os
import statistics
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
REPO = os.path.dirname(AN)
RUN = "2026-10-09"
SCRIPT = "analyses/figure_data_exports/export_suppnote1.py"
FIG = "suppnote6_pathways_coexpression"
OUT = os.path.join(os.path.expanduser("figure_data"), FIG)
POP = ("64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrence.tsv sha256 d11761227abd")
CO = os.path.join(REPO, "side-studies", "coexpression")
SNAP = os.path.join(CO, "snap_sz", "results")
FL = os.path.join(AN, "fairness_lists")
PHP = os.path.join(AN, "confounds", "ph_pathways_direction_2026-09-30", "results")
MOD = os.path.join(AN, "coexpression_modules_2026-10-08")
AX = os.path.join(AN, "coexpression_axes_2026-10-08")
RECOMP = os.path.join(OUT, "recomputed")

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


def rd(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t", quoting=csv.QUOTE_NONE))


def f(x):
    return float(x) if x not in ("", None) else None


meta = {"script": SCRIPT, "run_id": RUN, "figure": FIG, "population": POP, "inputs": {}, "panels": {}}


def write(name, panel, fields, rows, inputs, comments=()):
    os.makedirs(OUT, exist_ok=True)
    for p in inputs:
        meta["inputs"][rel(p)] = sha(p)
    lines = [f"# figure: {FIG} panel {panel}", f"# script: {SCRIPT}, run {RUN}", f"# population: {POP}"]
    lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs] + [f"# {c}" for c in comments]
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="") as fh:
        for l in lines:
            fh.write(l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for r in rows:
            w.writerow(["" if r.get(k) is None else r.get(k) for k in fields])
    meta["panels"][name] = {"panel": panel, "rows": len(rows)}
    print(name, len(rows))


# a ---------------------------------------------------------------- GO pathway agreement against gene agreement
def panel_a():
    p = os.path.join(FL, "a5b_agreement.tsv")
    rows = []
    for r in rd(p):
        for null, pre in (("uniform", "uniform_"), ("DE prior", "deprior_")):
            rows.append({"metric": r["metric"], "level": "pathway" if "pathway" in r["metric"] else "gene", "n_lists": r["n_lists"],
                         "n_pairs": r["n_pairs"], "observed": r["observed"], "null": null, "null_mean": r[pre + "null_mean"],
                         "null_lo95": r[pre + "null_lo95"], "null_hi95": r[pre + "null_hi95"], "ratio": r[pre + "ratio"],
                         "p_upper": r[pre + "p_upper"]})
    write("suppnote6_a_pathway_vs_gene_agreement.tsv", "a",
          ["metric", "level", "n_lists", "n_pairs", "observed", "null", "null_mean", "null_lo95", "null_hi95", "ratio", "p_upper"],
          rows, [p, os.path.join(FL, "a5b_go_meta.json")],
          ["100 null replicates per null, seed 20260925; p_upper = 0.0099 is the smallest value 100 replicates give",
           "a5b_agreement.tsv was run 2026-10-07; its input stamps predate the 2026-10-08 header-only restamp of recurrence.tsv; the 64 study ids, primary_lists.json, universes.json and the DE prior are byte-identical to the stamped ones (see STATUS.md)"])


# b ---------------------------------------------------------------- reported themes against nulls
def panel_b():
    p = os.path.join(RECOMP, "a5a_themes_null.tsv")
    p2 = os.path.join(RECOMP, "a5a_themes_vs_genes.tsv")
    rows = [dict(r) for r in rd(p)]
    write("suppnote6_b_themes_vs_null.tsv", "b",
          ["study_set", "n_studies", "n_excluded_no_theme", "null", "metric", "observed", "null_mean", "null_lo95", "null_hi95", "ratio_obs_null", "p_upper"],
          rows, [p, p2, os.path.join(RECOMP, "a5a_themes_meta.json")],
          ["recomputed 2026-10-09 from the unchanged analyses/uniform_vs_published_list_agreement/a5a_themes.py (1,000 permutations, seed 20260925; themes.json sha256 fc982383bb44, generated 2026-10-08)",
           "null 'marginal' = themes drawn at their marginal rates; null 'fixed' = each study keeps its own number of themes"])


# c ---------------------------------------------------------------- GO terms that carry the agreement
def panel_c():
    p = os.path.join(FL, "a5b_terms.tsv")
    tg = os.path.join(PHP, "term_groups.tsv")
    grp = {r["term"]: r["group"] for r in rd(tg)}
    rows = []
    for r in rd(p):
        rows.append({"term": r["term"], "name": r["name"], "group_in_pH_test": grp.get(r["term"], ""), "n_lists_top": r["n_lists_top"],
                     "uniform_null_mean": r["uniform_null_mean"], "deprior_null_mean": r["deprior_null_mean"],
                     "deprior_p_upper": r["deprior_p_upper"], "lists": r["lists"]})
    rows.sort(key=lambda r: -int(r["n_lists_top"]))
    write("suppnote6_c_go_terms_in_top_pathways.tsv", "c",
          ["term", "name", "group_in_pH_test", "n_lists_top", "uniform_null_mean", "deprior_null_mean", "deprior_p_upper", "lists"],
          rows, [p, tg],
          ["terms in the top pathways of >= 2 of the 62 lists; null means = expected number of lists; group_in_pH_test from the 2026-10-04 term groups (a5b_terms.tsv at commit 841baa2); list membership per term has since changed for some terms"])


# d ---------------------------------------------------------------- pH test of the term groups in controls
def panel_d():
    p = os.path.join(PHP, "controls_sets_combined.tsv")
    p2 = os.path.join(PHP, "controls_sets_by_dataset.tsv")
    keep = {"oxphos_hits", "oxphos_terms", "myelin_hits", "myelin_terms", "microglia_hits", "microglia_terms", "recurrent_q10"}
    rows = [r for r in rd(p) if r["set"] in keep and r["combination"] in ("PRIMARY", "ALT")]
    rows2 = [r for r in rd(p) if r["set"] in keep and r["combination"] == "PRIMARY, per gene"]
    perg = {r["set"]: r for r in rows2}
    out = []
    for r in rows:
        g = perg.get(r["set"], {}) if r["combination"] == "PRIMARY" else {}
        out.append({"combination": r["combination"], "set": r["set"], "n_datasets": r["n_datasets"], "n_positive": r["n_positive"],
                    "combined_perm_z": r["combined_perm_z"], "combined_p_one_sided": r["combined_p_one_sided"],
                    "n_genes_combined": g.get("n_genes_combined", ""), "n_q05_confound_sign": g.get("n_q05_confound_sign", "")})
    write("suppnote6_d_pH_test_term_groups.tsv", "d",
          ["combination", "set", "n_datasets", "n_positive", "combined_perm_z", "combined_p_one_sided", "n_genes_combined", "n_q05_confound_sign"],
          out, [p, p2],
          ["analysis run 2026-10-03/04; limma in controls, expression ~ pH + age + sex; z signed positive when the genes move with low pH the way the lists say they move in schizophrenia",
           "PRIMARY = GSE21138, GSE35977, GSE53987_BA46, GSE17612, GSE92538_GPL96; ALT = GSE208338, GSE35977, GSE202537_caudate, GSE17612, GSE80655_AnCg",
           "'*_hits' = list genes (membership from primary_lists.json sha 516475518ece, before the 2026-10-07 restamp); '*_terms' = all annotated genes of the term group (current library)"])
    # per dataset-region
    rows3 = [r for r in rd(p2) if r["factor"] == "pH" and r["set"] in ("oxphos_terms", "myelin_terms", "microglia_terms", "recurrent_q10")]
    write("suppnote6_d_pH_test_by_dataset.tsv", "d",
          ["dataset", "set", "n_ctl", "n_set_measured", "signed_score", "matched_mean", "matched_p", "perm_z", "perm_p"],
          rows3, [p2], ["per dataset-region signed score and sample-permutation z (pH among controls); matched = 10,000 expression-decile matched gene sets"])


# e ---------------------------------------------------------------- claim classes
def panel_e():
    p = os.path.join(CO, "coexpression_inventory_meta.json")
    m = json.load(open(p))
    order = ["held, independent donors", "mixed or not held, independent", "held, partly shared donors", "structure only", "same donors or cohort split", "none offered"]
    rows = []
    for c in order:
        rows.append({"sz_retest_class": c, "claims_all": m["claims_by_sz_retest"].get(c, 0), "claims_part_AB": m["claims_by_sz_retest_A_B"].get(c, 0),
                     "claims_part_C": m["claims_by_sz_retest_C"].get(c, 0)})
    rows.append({"sz_retest_class": "total", "claims_all": m["n_claims"], "claims_part_AB": m["n_claims_by_part"]["A/B (seed)"], "claims_part_C": m["n_claims_by_part"]["C (scope widening)"]})
    write("suppnote6_e_coexpression_claim_classes.tsv", "e", ["sz_retest_class", "claims_all", "claims_part_AB", "claims_part_C"], rows,
          [p, os.path.join(CO, "coexpression_inventory.tsv")],
          [f"headline claims = postmortem brain; {m['n_claims']} claims from {m['n_papers']} papers; inventory built {m['built']}; part A/B = seed search, part C = scope widening",
           "claims_by_sz_retest_A_B, _C and n_claims_by_part read from the meta file"])


# f ---------------------------------------------------------------- the claims re-tested in other donors
def panel_f():
    p = os.path.join(CO, "bank_overlap_claims.tsv")
    inv = os.path.join(CO, "coexpression_inventory.tsv")
    I = {r["claim_id"]: r for r in rd(inv)}
    keep = ("held, independent donors", "mixed or not held, independent", "held, partly shared donors")
    rows = []
    for r in rd(p):
        if r["sz_retest"] in keep:
            i = I[r["claim"]]
            rows.append({"claim": r["claim"], "paper": r["paper"], "pmid": i["pmid"], "discovery_banks": r["discovery"], "replication_banks": r["replication"],
                         "sz_retest": r["sz_retest"], "same_lineage_in_replication": r["same_lineage_in_replication"],
                         "discovery_data": i["discovery_data"], "replication_data": i["replication_data"], "independence": i["independence"]})
    write("suppnote6_f_retested_claims.tsv", "f",
          ["claim", "paper", "pmid", "discovery_banks", "replication_banks", "sz_retest", "same_lineage_in_replication", "discovery_data", "replication_data", "independence"],
          rows, [p, inv, os.path.join(CO, "bank_overlap_summary.tsv")],
          ["bank lineages are assigned by keyword rules in side-studies/coexpression_claims_survey/bank_overlap.py (CMC lineage = CommonMind, Mount Sinai, Penn, Pittsburgh, HBCC, LIBD, PsychENCODE; Stanley/SMRI incl. BrainGVEX)"])


# g ---------------------------------------------------------------- SNAP in other single-nucleus data
def panel_g():
    p = os.path.join(SNAP, "snap_sz_per_dataset.tsv")
    pm = os.path.join(SNAP, "snap_sz_meta.tsv")
    pa = os.path.join(SNAP, "arm_a_summary.tsv")
    pb = os.path.join(SNAP, "arm_b_delta.tsv")
    A = {}
    for r in rd(pa):
        if r["route"] == "resid_top50" and r["k"] == "10":
            A[r["set"]] = r
    D = {r["set"]: r for r in rd(pb)}
    role = {"Ling": "positive control", "HBCC_main": "primary", "MSSM_main": "primary", "SZBD_noLing": "primary", "NSW": "primary",
            "MSSM_nodementia": "sensitivity", "SZBD_all": "sensitivity", "CMC": "sensitivity", "Multiome": "sensitivity", "Batiuk": "sensitivity (neurons only)"}
    rows = []
    for r in rd(p):
        s = r["set"]
        a = A.get(s) or A.get({"HBCC_main": "HBCC_main"}.get(s, s)) or {}
        d = D.get(s, {})
        se = f(r["snap_sz_se"])
        rows.append({"row_type": "dataset", "set": s, "dataset": r["dataset"], "role": role.get(s, ""), "n_sz": r["n_sz"], "n_ctrl": r["n_ctrl"],
                     "snap_sz_sd": r["snap_sz"], "snap_sz_lo": round(f(r["snap_sz"]) - 1.96 * se, 4), "snap_sz_hi": round(f(r["snap_sz"]) + 1.96 * se, 4),
                     "snap_sz_null_median": r["snap_sz_null_median"], "snap_sz_p_perm": r["snap_sz_p_perm"],
                     "snap_age_ctrl_per_decade": r["snap_age_ctrl"], "snap_a_sz": r["snap_a_sz"], "snap_n_sz": r["snap_n_sz"],
                     "delta_on_minus_off": d.get("delta", ""), "delta_lo": d.get("delta_lo", ""), "delta_hi": d.get("delta_hi", ""),
                     "armA_best_rho_resid_top50_k10": a.get("best_rho", ""), "armA_call": a.get("call", "")})
    for r in rd(pm):
        if r["pool"] == "primary" or r["quantity"] == "SNAP total SZ":
            rows.append({"row_type": "pooled", "set": r["sets"], "dataset": f"{r['pool']}: {r['quantity']}", "role": "random effects (DerSimonian-Laird; Hartung-Knapp interval)",
                         "n_sz": "", "n_ctrl": "", "snap_sz_sd": r["estimate"], "snap_sz_lo": r["lo_hksj"], "snap_sz_hi": r["hi_hksj"], "snap_sz_p_perm": r["p_hksj"],
                         "snap_age_ctrl_per_decade": f"I2 {r['I2']}"})
    write("suppnote6_g_snap_other_snrnaseq.tsv", "g",
          ["row_type", "set", "dataset", "role", "n_sz", "n_ctrl", "snap_sz_sd", "snap_sz_lo", "snap_sz_hi", "snap_sz_null_median", "snap_sz_p_perm",
           "snap_age_ctrl_per_decade", "snap_a_sz", "snap_n_sz", "delta_on_minus_off", "delta_lo", "delta_hi", "armA_best_rho_resid_top50_k10", "armA_call"],
          rows, [p, pm, pa, pb],
          ["dataset rows: SZ effect in SD of Ling's fixed SNAP score; interval = estimate +/- 1.96 se (as combine.py); snap_sz_p_perm = permutation of the loadings, not SZ versus control",
           "pooled rows: snap_sz_p_perm holds the Hartung-Knapp p; interval is the Hartung-Knapp interval; one donor set per bank",
           "built 2026-09-28; the 2026-09-27 PsychAD file (snap_psychad_results.tsv) is reproduced exactly by this arm (HBCC -0.604)"])


# h ---------------------------------------------------------------- other coexpression modules on the pH and agonal axes
def panel_h():
    p = os.path.join(MOD, "modules_axes_all.tsv")
    R = rd(p)
    D = [r for r in R if r["sign"] not in ("0", "")]
    rows = []
    papers = sorted(set(r["paper"] for r in D))
    for pp in papers:
        m = [r for r in D if r["paper"] == pp]
        rows.append({"row_type": "paper", "paper": pp, "n_modules": len(m), "median_signed_pH": round(statistics.median(f(r["sig_pH"]) for r in m), 3),
                     "median_signed_agonal": round(statistics.median(f(r["sig_ag"]) for r in m), 3),
                     "modules_pH_side": sum(f(r["sig_pH"]) > 0 for r in m), "modules_agonal_side": sum(f(r["sig_ag"]) > 0 for r in m)})
    n = len(papers)
    ph = sum(r["median_signed_pH"] > 0 for r in rows)
    ag = sum(r["median_signed_agonal"] > 0 for r in rows)

    def sign_p(k):  # two-sided exact binomial, p = 0.5
        tail = sum(comb(n, i) for i in range(max(k, n - k), n + 1)) / 2 ** n
        return min(1.0, 2 * tail)
    rows.append({"row_type": "all direction-stated modules", "paper": f"{n} papers", "n_modules": len(D),
                 "median_signed_pH": round(statistics.median(f(r["sig_pH"]) for r in D), 3),
                 "median_signed_agonal": round(statistics.median(f(r["sig_ag"]) for r in D), 3),
                 "modules_pH_side": sum(f(r["sig_pH"]) > 0 for r in D), "modules_agonal_side": sum(f(r["sig_ag"]) > 0 for r in D)})
    rows.append({"row_type": "papers on the pH / agonal side (sign test, two-sided)", "paper": f"pH {ph} of {n}, p {sign_p(ph):.4f}; agonal {ag} of {n}, p {sign_p(ag):.4f}", "n_modules": ""})
    write("suppnote6_h_modules_signed_axes_by_paper.tsv", "h",
          ["row_type", "paper", "n_modules", "median_signed_pH", "median_signed_agonal", "modules_pH_side", "modules_agonal_side"], rows, [p],
          ["signed score = for a module stated higher (+1) or lower (-1) in schizophrenia: sig_pH = -sign x mean control pH Z, sig_agonal = sign x mean GTEx agonal t; positive = module moves the way schizophrenia would if it were low pH / agonal",
           "modules with a stated direction and linked to schizophrenia by their paper; modules from one paper share genes and donors, so counts are indicative"])


# i ---------------------------------------------------------------- Gandal and CoPA modules on the axes
def panel_i():
    p = os.path.join(AX, "module_axes_summary.tsv")
    rows = [{"set": r["set"], "comparison": r["comparison"], "axis": r["axis"], "n_or_medians": r["n_or_medians"], "rho": r["rho"], "p": r["p"]} for r in rd(p)]
    write("suppnote6_i_gandal_copa_axes.tsv", "i", ["set", "comparison", "axis", "n_or_medians", "rho", "p"], rows,
          [p, os.path.join(AX, "gandal_modules_axes.tsv"), os.path.join(AX, "copa_modules_axes.tsv")],
          ["Gandal 2018 Table S5: 34 modules, Spearman across modules of the module's disease t against its mean axis value; CoPA (Kang and Oldham 2026 Table S14): Mann-Whitney of module mean axis, lower versus higher in schizophrenia"])


if __name__ == "__main__":
    for fn in (panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, panel_g, panel_h, panel_i):
        fn()
    with open(os.path.join(OUT, "_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
