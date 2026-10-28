#!/usr/bin/env python3
"""Per-panel tidy exports for manuscript Figures 2 (agreement) and 7 (pH).

Request: MAB_TO_RA_2026_10_08_FRESH_FIGURES_PANEL_DATA.md, Priority 1.
Population: the 64 published brain studies (recurrence.py scenario 'published, brain', study tier);
recurrent set = the 176 study-tier genes at q < 0.1 (uniform null).

Reads existing analysis outputs only (no model fits, no permutation reruns, no network). Writes new files to
figure_data/fig02/ and fig07/, each TSV with '#' provenance lines, plus
_meta.json and STATUS.md per folder. One function per panel; rerun with

    python3 analyses/figure_data_exports/export_fig02_fig07.py
"""
import csv
import collections
import hashlib
import json
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)            # analyses/
REPO = os.path.dirname(AN)
sys.path.insert(0, AN)
import recurrence as R  # noqa: E402  (the repo's own loading path; light, ~2 s)

RUN = "2026-10-08"
SCRIPT = "analyses/figure_data_exports/export_fig02_fig07.py"
OUT = os.path.expanduser("figure_data")
SCENARIO = "published, brain"
POP = ("64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrent set = the 176 "
       "study-tier genes at q < 0.1 (uniform null); recurrence.tsv sha256 d11761227abd")

# ---------------------------------------------------------------- helpers
_sha_cache = {}


def sha(path):
    path = os.path.abspath(path)
    if path not in _sha_cache:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _sha_cache[path] = h.hexdigest()
    return _sha_cache[path]


def rel(path):
    return os.path.relpath(os.path.abspath(path), REPO)


def rd(path):
    """Rows of a TSV, '#' lines skipped, values kept as the strings in the file."""
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


class Figure:
    def __init__(self, fig):
        self.fig = fig
        self.dir = os.path.join(OUT, fig)
        os.makedirs(self.dir, exist_ok=True)
        self.meta = {"script": SCRIPT, "run_id": RUN, "figure": fig, "population": POP, "inputs": {}, "panels": {}}
        self.status = []

    def note(self, panel, text):
        self.status.append(f"- **{panel}**: {text}")

    def write(self, name, panel, fields, rows, inputs, extra_comments=()):
        """Write one TSV with provenance comment lines; record row count and input hashes."""
        for p in inputs:
            self.meta["inputs"][rel(p)] = sha(p)
        lines = [f"# figure: {self.fig} panel {panel}",
                 f"# script: {SCRIPT}, run {RUN}",
                 f"# population: {POP}"]
        lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs]
        lines += [f"# {c}" for c in extra_comments]
        path = os.path.join(self.dir, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            for l in lines:
                fh.write(l + "\n")
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(fields)
            n = 0
            for r in rows:
                w.writerow(["" if r.get(f) is None else r.get(f) for f in fields])
                n += 1
        self.meta["panels"][name] = {"panel": panel, "rows": n}
        print(f"{self.fig}/{name}: {n} rows")
        return n

    def finish(self):
        with open(os.path.join(self.dir, "_meta.json"), "w") as fh:
            json.dump(self.meta, fh, indent=1)
        with open(os.path.join(self.dir, "STATUS.md"), "w") as fh:
            fh.write(f"# {self.fig} status ({RUN})\n\n" + "\n".join(self.status) + "\n")


def scenario_row(path, scenario=SCENARIO, tier="study"):
    for r in rd(path):
        if r["scenario"] == scenario and r["tier"] == tier:
            return r
    raise KeyError(scenario)


P_REC = os.path.join(AN, "recurrence.tsv")
P_GENES = os.path.join(AN, "recurrence_genes.tsv")
P_UNITS = os.path.join(AN, "recurrence_units.tsv")
P_LINEAGE = os.path.join(REPO, "corpus", "lineage.json")
P_TECH = os.path.join(AN, "technology.json")
P_UNIV = os.path.join(AN, "universes.json")
P_PRIM = os.path.join(AN, "primary_lists.json")


def lineage_years():
    nodes = json.load(open(P_LINEAGE))["nodes"]
    return {str(n.get("id")): n.get("year") for n in nodes}, {str(n.get("pmid")): n.get("year") for n in nodes}


# ================================================================== Figure 2
def fig2_a(F):
    """Recurrence curve: observed vs three nulls at each k (gene_coverage/recurrence_curve.tsv)."""
    p = os.path.join(AN, "gene_coverage", "recurrence_curve.tsv")
    rows = rd(p)
    out = [{"k": r["n_studies_at_least"], "observed": r["observed_genes"],
            "exp_uniform": r["expected_uniform"], "sd_uniform": r["sd_uniform"],
            "exp_deprior": r["expected_de_prior"], "sd_deprior": r["sd_de_prior"],
            "exp_exprmatched": r["expected_expression"], "sd_exprmatched": r["sd_expression"]} for r in rows]
    F.write("fig02_a_recurrence_curve.tsv", "a",
            ["k", "observed", "exp_uniform", "sd_uniform", "exp_deprior", "sd_deprior", "exp_exprmatched", "sd_exprmatched"],
            out, [p, P_REC],
            ["k = genes found DE in at least k of the 64 studies; sd = sqrt(sum P_g (1 - P_g)) over genes, from recurrence_curve.py"])
    # cross-check against recurrence.tsv (study tier, 'published, brain')
    s = scenario_row(P_REC)
    bad = []
    for r in rows:
        k = r["n_studies_at_least"]
        if k in ("2", "3", "5", "10"):
            for a, b in (("observed_genes", "obs_ge" + k), ("expected_uniform", "exp_ge" + k),
                         ("expected_de_prior", "exp_prior_ge" + k), ("expected_expression", "exp_expr_ge" + k)):
                if abs(float(r[a]) - float(s[b])) > 0.051:
                    bad.append(f"k={k} {a} {r[a]} vs recurrence.tsv {b} {s[b]}")
    F.note("a", "recurrence_curve.tsv agrees with recurrence.tsv at k = 2, 3, 5, 10 (all three nulls)." if not bad
           else "MISMATCH curve vs recurrence.tsv: " + "; ".join(bad))
    F.note("a", f"recurrence.tsv alone has only k >= 2, 3, 5, 10 and no SDs; the full curve and SDs come from gene_coverage/recurrence_curve.tsv "
                f"(run {RUN}, same input shas). Genes q < 0.1: uniform {s['genes_q10']}, DE prior {s['genes_q10_prior']}, "
                f"expression-matched {s['genes_q10_expr']}; report/report.md states 99 for the expression-matched set "
                f"(lines 95, 168, 173), recurrence.tsv and recurrence_report.md state {s['genes_q10_expr']}.")


def fig2_b(F):
    """Gene x study matrix for the top 100 genes + study annotation (same loading path as recurrence_gene_list_matrix.py)."""
    sel, M, lists, units = R.scenario_selection(SCENARIO)
    grows = [r for r in rd(P_GENES) if r["scenario"] == SCENARIO]
    grows.sort(key=lambda r: (-int(r["n_lists_observed"]), r["gene"]))
    top = grows[:100]
    ynode, ypmid = lineage_years()
    named = {"cmc+mssm+pitt+libd", "stanley", "aus", "harvard", "pritzker"}

    def year(k):
        y = ynode.get(k) or ypmid.get(k)
        if y:
            return int(y)
        m = re.search(r"(\d{4})", units["labels"].get(k, ""))
        return int(m.group(1)) if m else None

    cluster = {k: R.cluster_of(units, k) for k in sel}
    order = sorted(sel, key=lambda k: (cluster[k] if cluster[k] in named else "~other", year(k) or 0, k))
    out, mism, nodir = [], [], 0
    for rank, tr in enumerate(top, 1):
        g = tr["gene"]
        hits = 0
        for k in order:
            d = lists.get(k, {}).get(g)
            if d in ("up", "down"):
                call = d
            elif d == "both":
                call = "mixed"
            elif d is not None:
                call = "hit_no_direction"
                nodir += 1
            elif g in M.get(k, ()):
                call = "tested_not_hit"
            else:
                call = "not_tested"
            if call in ("up", "down", "mixed", "hit_no_direction"):
                hits += 1
            out.append({"gene": g, "study_id": k, "call": call, "gene_rank": rank, "n_lists_observed": tr["n_lists_observed"]})
        if hits != int(tr["n_lists_observed"]):
            mism.append(f"{g}: {hits} vs {tr['n_lists_observed']}")
    F.write("fig02_b_gene_by_study_long.tsv", "b", ["gene", "study_id", "call", "gene_rank", "n_lists_observed"], out,
            [P_GENES, P_UNITS, P_PRIM, P_UNIV],
            ["top 100 genes of scenario 'published, brain' ranked by n_lists_observed desc then gene (recurrence_genes.tsv); "
             "call: up/down/mixed = direction in the study's primary list (mixed = both within the study's union); "
             "tested_not_hit = gene in the list's tested universe but not a hit; not_tested = outside it; "
             "hit_no_direction = hit with no stated direction",
             "built with recurrence.scenario_selection (same path as recurrence_gene_list_matrix.py)"])
    ann = []
    _uni = json.load(open(P_UNIV))["lists"]
    nk = R.load_lists(units, _uni)[2]
    tclass = {r["pmid"]: r["tested_class"] for r in csv.DictReader((l for l in open(os.path.join(AN, "corpus_516_prep_2026-10-10", "tested_genes_per_study_516.tsv"), encoding="utf8") if not l.startswith("#")), delimiter="\t")}
    no_tested = []
    for i, k in enumerate(order, 1):
        st = units["studies"].get(k, {})
        L = lists.get(k, {})
        Mk = M.get(k)
        if Mk is None or k not in nk:
            no_tested.append(k)
        ann.append({"study_id": k, "label": units["labels"].get(k, ""), "cohort_cluster": cluster[k],
                    "year": year(k), "technology": st.get("technology", ""),
                    "n_list_genes": len(set(L) & Mk) if Mk is not None else len(L),
                    "column_order": i,
                    "n_tested": ("" if k not in nk else round(nk[k])),
                    "n_tested_class": tclass.get(k, ""),
                    "n_list_genes_all": len(L)})
    F.write("fig02_b_study_annotation.tsv", "b",
            ["study_id", "label", "cohort_cluster", "year", "technology", "n_list_genes", "column_order", "n_tested", "n_tested_class", "n_list_genes_all"], ann,
            [P_UNITS, P_LINEAGE, P_TECH, P_UNIV],
            ["year: corpus/lineage.json node year (label year if the node is missing); technology: technology.json via recurrence_units; "
             "column_order: cohort cluster (named clusters first) then year, as in recurrence_gene_list_matrix.py",
             "n_tested = n_k of the study's universe as recurrence.py uses it (universes.json, recurrence.load_lists; same value as tested_genes_per_study_516.tsv n_tested); "
             "n_tested_class from that file (measured_gene_set = gene set held; assumed_rnaseq = the 25,137-gene RNA-seq default assumed, not measured); "
             "n_list_genes = |list genes in the study's tested set| (the DE count on the same footing as n_tested); n_list_genes_all = all genes in the study's primary list(s), the earlier n_list_genes definition",
             f"studies without a tested set: {len(no_tested)}"])
    _diff = sum(1 for a in ann if a["n_list_genes"] != a["n_list_genes_all"])
    F.note("b", f"n_tested added (n_k per study, recurrence.py basis); n_list_genes is now |L intersect M| and differs from the earlier all-list-genes count (kept as n_list_genes_all) in {_diff} of {len(ann)} studies; "
                f"{len(no_tested)} studies have no tested set. n_tested class counts: {dict(collections.Counter(a['n_tested_class'] for a in ann))}.")
    F.note("b", ("cell counts (up + down + mixed + no-direction hits) equal n_lists_observed for all 100 genes."
                 if not mism else "MISMATCH hits vs n_lists_observed: " + "; ".join(mism[:10])))
    if nodir:
        F.note("b", f"{nodir} gene x study hits have no stated direction and are labelled hit_no_direction "
                    f"(recurrence_gene_list_matrix.py draws them as tested_not_hit).")
    multi = [k for k in order if cluster[k] not in named and not str(cluster[k]).startswith("brain_other")]
    F.note("b", f"{len(order)} studies in the matrix (expected 64). Cohort cluster = recurrence.cluster_of; "
                f"corpus/lineage.json has changed since the recurrence run (sha now {sha(P_LINEAGE)[:12]}, run used 6f4ba89306a7; "
                f"commit 2026-10-08 18:22 adds Muralidharan 2025); it is read here only for publication years.")


def fig2_c(F):
    p = os.path.join(AN, "h1_pairs.tsv")
    rows = [r for r in rd(p) if r["scenario"] == SCENARIO]
    cmap = {"same": "yes", "different": "no", "multi-bank": "multi-bank"}
    out = [{"pair_id": f"{r['study_a']}_{r['study_b']}", "era_pair": r["era_2015"], "fold": r["fold"],
            "shares_any_gene": "yes" if int(r["overlap"]) > 0 else "no", "same_cohort": cmap[r["cohort"]],
            "overlap": r["overlap"], "expected": r["expected"], "era_pair_2020": r["era_2020"],
            "label_a": r["label_a"], "label_b": r["label_b"], "year_a": r["year_a"], "year_b": r["year_b"]} for r in rows]
    F.write("fig02_c_pairwise_fold_by_era.tsv", "c",
            ["pair_id", "era_pair", "fold", "shares_any_gene", "same_cohort", "overlap", "expected", "era_pair_2020",
             "label_a", "label_b", "year_a", "year_b"], out, [p],
            ["pair_id = study_a_study_b (PMID keys; A is the earlier study); era_pair: old <= 2015, new >= 2016 (era_2015), era_pair_2020 = sensitivity split old <= 2020; "
             "fold = overlap / expected under the pair's shared tested universe; shares_any_gene = overlap > 0; "
             "same_cohort yes/no from h1.py 'cohort' (recurrence.cluster_of), 'multi-bank' = a multi-bank study without per-bank lists"])
    # pooled folds by era pair (h1_summary.tsv)
    ps = os.path.join(AN, "h1_summary.tsv")
    keep = [r for r in rd(ps) if r["scenario"] == SCENARIO and r["cohort"] == "all"]
    F.write("fig02_c_pooled_by_era.tsv", "c",
            ["era_split", "pair_type", "n_pairs", "n_studies", "median_fold", "null_median_fold", "null_median_fold_lo",
             "null_median_fold_hi", "pooled_fold", "n_pairs_overlap0", "pooled_top_fold", "median_concordance", "n_concordance"],
            keep, [ps], ["era_split 2015 = old <= 2015 / new >= 2016; null_median_fold_lo/hi = interval of the era-label permutation null; all cohorts"])
    n_overlap0 = sum(1 for r in out if r["shares_any_gene"] == "no")
    F.note("c", f"{len(out)} pairs (64 choose 2 = 2016); {n_overlap0} share no gene. "
                f"Pooled fold all pairs {[r['pooled_fold'] for r in keep if r['era_split']=='2015' and r['pair_type']=='all'][0]} (fairness_lists report: 1.35).")
    # Gemma reanalysis (fairness_lists/a1)
    a1s, a1b, a1p = (os.path.join(AN, "fairness_lists", f) for f in ("a1_summary.tsv", "a1_boot.tsv", "a1_pairs.tsv"))
    srows = [r for r in rd(a1s) if r["pair_class"] == "brain-brain, different groups" and r["method"] in ("published", "gemma_sized")]
    brows = [r for r in rd(a1b) if r["pair_class"] == "brain-brain, different groups" and r["comparison"] == "gemma_sized minus published"]
    res = []
    for r in srows:
        res.append({"kind": "pooled", "method": r["method"], "statistic": "pooled_fold", "estimate": r["pooled_fold"],
                    "boot_lo95": "", "boot_hi95": "", "n_pairs": r["n_pairs"], "n_studies": r["n_studies"]})
    for r in brows:
        if r["statistic"] in ("pooled_fold", "median_log2_fold"):
            res.append({"kind": "difference", "method": "gemma_sized minus published", "statistic": r["statistic"],
                        "estimate": r["estimate"], "boot_lo95": r["boot_lo95"], "boot_hi95": r["boot_hi95"],
                        "n_pairs": r["n_pairs"], "n_studies": ""})
    F.write("fig02_c_gemma_reanalysis_pooled.tsv", "c",
            ["kind", "method", "statistic", "estimate", "boot_lo95", "boot_hi95", "n_pairs", "n_studies"], res, [a1s, a1b],
            ["brain-brain pairs from different groups, 25 datasets with Gemma reanalyses; method published = the study's own list, gemma_sized = Gemma list cut to the same size; "
             "interval = dataset-level bootstrap (1000, seed 20260925) of the difference gemma minus published, not of each pooled fold"])
    prow = [r for r in rd(a1p) if r["tissue_pair"] == "brain-brain" and r["same_group"] == "no" and r["method"] in ("published", "gemma_sized")]
    pout = [{"pair_id": f"{r['study_a']}_{r['study_b']}", "method": r["method"], "fold": r["fold"], "overlap": r["overlap"],
             "expected": r["expected"], "label_a": r["label_a"], "label_b": r["label_b"]} for r in prow]
    F.write("fig02_c_gemma_reanalysis_pairs.tsv", "c", ["pair_id", "method", "fold", "overlap", "expected", "label_a", "label_b"],
            pout, [a1p], ["per-pair folds behind the pooled rows; brain-brain, different groups"])
    npair = sorted({r["pair_id"] for r in pout})
    F.note("c", f"Gemma reanalysis comes from fairness_lists/a1_* (run 2026-10-07 21:25, one day before the 10-08 core run). Its primary_lists.json and universes.json "
                f"shas match the current files; its meta lists h1_pairs.tsv/recurrence.tsv shas from the earlier run, but a1_uniform_vs_published.py does not read those "
                f"two files. {len(npair)} distinct brain-brain pairs. No interval exists for each pooled fold separately.")


def fig2_d(F):
    ps = os.path.join(AN, "..", "side-studies", "hoffman-2021", "followup_summary.tsv")
    pst = os.path.join(AN, "..", "side-studies", "hoffman-2021", "followup_summary_studies.tsv")
    pdt = os.path.join(AN, "..", "side-studies", "hoffman-2021", "followup_data.tsv")
    label = {"no Hoffman bank": "other banks", "shares a Hoffman bank or donors": "shared bank",
             "CommonMind / HBCC data (same donors)": "CMC-HBCC reanalyses"}
    rows = rd(ps)
    out = [{"group": label[r["group"]], "group_in_source": r["group"], "n_studies": r["studies"], "observed": r["overlap"],
            "exp_uniform": r["expected"], "exp_deprior": r["expected_de_prior"], "n_same_dir": r["same_direction"],
            "n_opposite": r["opposite_direction"], "studies_p_lt_0.05": r["studies_p_lt_0.05"]} for r in rows]
    order = {"other banks": 0, "shared bank": 1, "CMC-HBCC reanalyses": 2}
    out.sort(key=lambda r: order[r["group"]])
    F.write("fig02_d_hoffman82_by_group.tsv", "d",
            ["group", "group_in_source", "n_studies", "observed", "exp_uniform", "exp_deprior", "n_same_dir", "n_opposite", "studies_p_lt_0.05"],
            out, [ps, pdt],
            ["gene set both_82: the 82 genes significant (FDR 5%) in both Hoffman et al. 2021 cohorts; later studies 2021 on, meta-analyses and preprints excluded; "
             "observed = summed overlap of the 82 with each study's list; expectations = uniform hypergeometric within each study's tested universe / DE prior, summed; "
             "n_same_dir/n_opposite = overlap genes in / against Hoffman's direction"])
    srows = rd(pst)
    sout = [{"group": label[r["group"]], "study": r["study"], "year": r["year"], "testable": r["testable"], "observed": r["overlap"],
             "exp_uniform": r["expected"], "p_hyper": r["p_hyper"], "exp_deprior": r["expected_de_prior"],
             "n_same_dir": r["same_direction"], "n_opposite": r["opposite_direction"]} for r in srows]
    F.write("fig02_d_hoffman82_by_study.tsv", "d",
            ["group", "study", "year", "testable", "observed", "exp_uniform", "p_hyper", "exp_deprior", "n_same_dir", "n_opposite"],
            sout, [pst], ["one row per later study"])
    got = {r["group"]: (r["n_studies"], r["observed"], r["exp_uniform"], r["exp_deprior"], r["n_same_dir"], r["n_opposite"]) for r in out}
    F.note("d", "MISMATCH with report/report.md line 183 (Hoffman paragraph): the report text gives other banks 29 studies, 108 listed, 72.5 / 105.3 expected, 81 same / 13 opposite; "
                "shared bank 13 studies, 97, 53.1 / 77.6; CMC-HBCC reanalyses 9, 393, 142.5 / 199.0, 390 same. followup_summary.tsv (run 2026-10-08) gives "
                "other banks " + "/".join(got['other banks']) + ", shared bank " + "/".join(got['shared bank']) + ", CMC-HBCC " + "/".join(got['CMC-HBCC reanalyses']) + " "
                "(order: n_studies, observed, exp_uniform, exp_deprior, same, opposite). Only the CMC-HBCC row matches (opposite 3 in the file; the report gives only the 390 same). "
                "The exported numbers are the file's.")


def fig2_e(F):
    p = os.path.join(AN, "h3", "h3_pairs.tsv")
    rows = [r for r in rd(p) if r["tissue_set"] == "brain" and r["both_cross_cohort"] == "yes"]
    out = []
    for r in rows:
        out.append({"pair_id": f"{r['unit_a']}_{r['unit_b']}", "share_inputs": "yes" if int(r["n_shared_datasets"]) > 0 else "no",
                    "fold_deprior": r["fold_deprior"], "same_technology": "yes" if r["tech_a"] == r["tech_b"] else "no",
                    "shares_cohort_cluster": "yes" if int(r["n_shared_cohorts"]) > 0 else "no",
                    "n_shared_datasets": r["n_shared_datasets"], "dataset_jaccard": r["dataset_jaccard"],
                    "n_shared_cohorts": r["n_shared_cohorts"], "overlap": r["overlap"], "expected_deprior": r["expected_deprior"],
                    "fold_uniform": r["fold_uniform"], "label_a": r["label_a"], "label_b": r["label_b"],
                    "tech_a": r["tech_a"], "tech_b": r["tech_b"]})
    F.write("fig02_e_meta_analysis_pairs.tsv", "e",
            ["pair_id", "share_inputs", "fold_deprior", "same_technology", "shares_cohort_cluster", "n_shared_datasets", "dataset_jaccard",
             "n_shared_cohorts", "overlap", "expected_deprior", "fold_uniform", "label_a", "label_b", "tech_a", "tech_b"], out, [p],
            ["pairs of the 16 brain cross-cohort meta-analyses with lists (tissue_set brain, both_cross_cohort yes); pair_id = unit_a_unit_b; "
             "share_inputs = at least one input dataset in common; same_technology = tech_a == tech_b; shares_cohort_cluster = n_shared_cohorts > 0"])
    z = [r for r in rows if r["n_shared_datasets"] == "0"]
    c0 = [r for r in rows if r["n_shared_cohorts"] == "0"]
    pool = lambda rs: sum(float(r["overlap"]) for r in rs) / sum(float(r["expected_deprior"]) for r in rs)
    units = {r["unit_a"] for r in rows} | {r["unit_b"] for r in rows}
    F.note("e", f"{len(rows)} pairs, {len(units)} meta-analyses; {len(z)} pairs share no input dataset (pooled fold over the DE prior {pool(z):.2f}) and "
                f"{len(c0)} share no cohort cluster ({pool(c0):.2f}); report/report.md gives 120 / 62 / 1.22 / 7 / 1.06: match.")


# ================================================================== Figure 7
P_CONF = os.path.join(AN, "confounds")
P_RGPH = os.path.join(P_CONF, "results", "recurrent_genes_pH.tsv")
P_CTLPH = os.path.join(P_CONF, "results", "combined_controls_pH.tsv")
P_CSUM = os.path.join(P_CONF, "results", "combined_summary.tsv")
P_GD = os.path.join(P_CONF, "results", "group_differences.tsv")


def fig7_a1(F):
    rows = [r for r in rd(P_GD) if r["variable"] == "pH"]
    out = [{"dataset_region": r["dataset"], "delta_pH": r["diff"], "p": r["t_p"], "wilcox_p": r["wilcox_p"], "n_sz": r["n_sz"],
            "n_ctl": r["n_ctl"], "mean_sz": r["mean_sz"], "mean_ctl": r["mean_ctl"],
            "in_14_dataset_regions": "no" if r["dataset"] == "GSE208338_noshared" else "yes"} for r in rows]
    F.write("fig07_a1_ph_difference.tsv", "a1",
            ["dataset_region", "delta_pH", "p", "wilcox_p", "n_sz", "n_ctl", "mean_sz", "mean_ctl", "in_14_dataset_regions"], out, [P_GD],
            ["delta_pH = mean pH in SZ minus mean pH in controls; p = two-sided t-test (t_p of group_differences.tsv), wilcox_p = Wilcoxon; "
             "GSE208338_noshared is the same dataset with donors shared with GSE21138 removed (a sensitivity variant, flagged in_14_dataset_regions = no)"])
    main = [r for r in out if r["in_14_dataset_regions"] == "yes"]
    neg = sum(float(r["delta_pH"]) < 0 for r in main)
    sig = [r["dataset_region"] for r in main if float(r["p"]) < 0.05]
    F.note("a1", f"{len(main)} dataset-regions (plus the GSE208338_noshared variant); pH lower in SZ in {neg} of {len(main)}; t-test p < 0.05 in {len(sig)} ({', '.join(sig)}). "
                 "confounds_report.md: 14 of 14 lower, 4 significant (GSE35977, GSE17612, GSE21935, GSE80655_AnCg).")


def fig7_a2(F):
    rows = rd(P_RGPH)
    out = []
    for r in rows:
        sets = r["sets"].split(",")
        out.append({"gene": r["gene"], "z_ph_controls": r["Z_pH_controls"], "z_sz_raw": r["Z_sz_base"], "z_sz_ph_adjusted": r["Z_sz_pH"],
                    "in_176": "yes" if "study_brain" in sets else "no", "in_5plus": "yes" if "study_brain_5plus" in sets else "no",
                    "sz_direction": r["sz_direction"], "q_ph_controls": r["q_pH_controls"], "n_datasets_ph": r["n_datasets_pH"],
                    "n_datasets_sz": r["n_datasets_sz"], "ph_confound_consistent": r["pH_confound_consistent"]})
    F.write("fig07_a2_a4_gene_scatter.tsv", "a2/a4",
            ["gene", "z_ph_controls", "z_sz_raw", "z_sz_ph_adjusted", "in_176", "in_5plus", "sz_direction", "q_ph_controls",
             "n_datasets_ph", "n_datasets_sz", "ph_confound_consistent"], out, [P_RGPH],
            ["z_ph_controls = Stouffer Z of control-only pH association over the 5 donor-separate datasets (Z > 0: higher expression at higher pH); "
             "z_sz_raw / z_sz_ph_adjusted = Stouffer Z (weights sqrt n) of the SZ effect, limma expr ~ dx + age + sex without / with pH, over the datasets with pH; "
             "in_176 = study_brain set (176 genes), in_5plus = study_brain_5plus (1,478 genes); blank = not measured in enough datasets",
             "source scripts: confounds/scripts/combine_datasets.R (run 2026-10-08)"])
    n176 = sum(r["in_176"] == "yes" for r in out)
    n5 = sum(r["in_5plus"] == "yes" for r in out)
    miss = sum(1 for r in out if r["in_176"] == "yes" and not r["z_ph_controls"])
    F.note("a2/a4", f"{len(out)} genes: the 5+ set ({n5} genes), of which {n176} are the 176 recurrent genes. {miss} of the 176 have no control pH Z (measured in fewer than 4 of 5 datasets). "
                    "z_sz_ph_adjusted comes from confounds/results (run 2026-10-08), not the older ph_adjustment_2026-10-01 (2026-10-03).")


def fig7_a3(F):
    ctl = {r["gene"]: r for r in rd(P_CTLPH) if r["factor"] == "pH"}
    sets = {}
    for r in rd(P_RGPH):
        for s in r["sets"].split(","):
            sets.setdefault(s, []).append(r)
    summ = [r for r in rd(P_CSUM) if r["variant"] == "primary_without_GSE208338"
            and r["test"].startswith("combined control pH q < 0.05, expression-matched null")]
    names = {"cohort_brain": "cohort tier, 10 genes", "study_brain": "176 study-tier genes", "cohort_brainFDR": "cohort tier, FDR lists",
             "study_brainFDR": "study tier, FDR lists", "study_brain_5plus": "5+ studies, 1,478 genes"}
    out = []
    check = []
    for r in summ:
        s = r["set"]
        members = sets[s]
        meas = [g["gene"] for g in members if g["gene"] in ctl]
        any_sign = sum(float(ctl[g]["q"]) < 0.05 for g in meas)
        cons = sum(1 for g in members if g["gene"] in ctl and float(ctl[g["gene"]]["q"]) < 0.05 and g["pH_confound_consistent"] == "TRUE")
        check.append(f"{s}: {cons} vs summary {r['observed']}")
        out.append({"set": s, "set_label": names.get(s, s), "n_set": len(members), "n_measured": len(meas), "observed": r["observed"],
                    "expected": r["expected"], "p": r["combined_p_one_sided"], "observed_any_sign": any_sign})
    F.write("fig07_a3_observed_vs_expected.tsv", "a3",
            ["set", "set_label", "n_set", "n_measured", "observed", "expected", "p", "observed_any_sign"], out, [P_CSUM, P_CTLPH, P_RGPH],
            ["observed = set genes at combined control pH q < 0.05 with the confound-consistent sign (the SZ-up genes higher at low pH, SZ-down lower); "
             "expected = mean of 10,000 expression-matched random sets; p = their one-sided p (combined_summary.tsv, variant primary_without_GSE208338: 5 donor-separate datasets)",
             "n_measured = set genes with a combined control pH row (measured in >= 4 of 5 datasets), counted here from combined_controls_pH.tsv and recurrent_genes_pH.tsv; "
             "observed_any_sign = set genes at q < 0.05 whatever the sign, counted here"])
    ok = all(c.split(": ")[1].split(" vs summary ")[0] == c.split("summary ")[1] for c in check)
    F.note("a3", ("observed (confound-consistent) recounted from gene files equals combined_summary.tsv for every set (" + "; ".join(check) + "). "
                  if ok else "MISMATCH recount vs summary: " + "; ".join(check) + ". ")
           + "Report text: cohort tier 5 vs 0.99 (9 measured), 5+ set 541 vs 245.82 (1,431 measured), 647 any sign: match. "
             "The request names 'set | n_measured | observed | expected | p'; no n_measured column exists in a file, so it is counted from the gene tables.")


P_FSL = os.path.join(AN, "fairness_sample_level", "results")
LEVEL_DESC = {"L0": "no covariates", "L1": "age + sex", "L2": "age + sex + PMI", "L3": "age + sex + PMI + pH, RIN, agonal factor (where recorded)",
              "L4": "L3 + medication (GSE21138, GSE138082 only)", "L5a": "age + sex + surrogate variables",
              "L5b": "L4 + surrogate variables"}


def fig7_b(F):
    pp, pg, pl, pi = (os.path.join(P_FSL, f) for f in ("b_pairs.tsv", "b_recurrent_genes.tsv", "b_levels.tsv", "b_independent_set.tsv"))
    rows = rd(pp)
    out = [{"pair_id": f"{r['dataset_a']}_{r['dataset_b']}", "model_step": r["level"], "rho": r["rho_t"], "model_description": LEVEL_DESC[r["level"]],
            "relation": r["relation"], "both_bulk": r["both_bulk"], "n_common": r["n_common"]} for r in rows]
    F.write("fig07_b_covariate_series_pairs.tsv", "b",
            ["pair_id", "model_step", "rho", "model_description", "relation", "both_bulk", "n_common"], out, [pp],
            ["rho = Spearman correlation of the SZ vs control limma t statistics of the two datasets (rho_t); levels L0-L5b = covariate adjustment series (scripts/b_covariate_series.R); "
             "relation = subject relation of the two datasets (subject-separate pairs are the ones plotted in the report); both_bulk TRUE = neither is laser-captured",
             "POPULATION: the 28 datasets with sample-level data (Gemma / GEO matrices), not the 64 published studies; the recurrent genes of panel b2 are the current 176"])
    g = rd(pg)
    wide = {}
    for r in g:
        w = wide.setdefault(r["gene"], {"gene": r["gene"], "sz_direction": r["sz_direction"]})
        w["z_" + r["level"]] = r["signed_stouffer_z"]
        w["n_datasets_" + r["level"]] = r["n_datasets"]
    gout = []
    for gene, w in wide.items():
        w["z_age_sex"] = w.get("z_L1")
        w["z_all_covariates"] = w.get("z_L3")
        gout.append(w)
    fields = ["gene", "sz_direction", "z_age_sex", "z_all_covariates"] + [f"z_{l}" for l in LEVEL_DESC] + ["n_datasets_L1", "n_datasets_L3"]
    F.write("fig07_b_covariate_series_genes.tsv", "b",
            fields, gout, [pg, pi],
            ["z = Stouffer z of the SZ effect over the independent dataset set (b_independent_set.tsv), signed by the gene's recurrent SZ direction (positive = SZ direction); "
             "z_age_sex = L1; z_all_covariates = L3 (PMI, pH, RIN and agonal factor added, as in the report legend 'PMI, pH, RIN and agonal factor added'); "
             "L4 (adds medication where recorded) and L5a/L5b (surrogate variables) are in the z_L* columns",
             "the gene set is the 176 recurrent genes"])
    lv = rd(pl)
    F.write("fig07_b_levels_by_dataset.tsv", "b",
            ["dataset", "compartment", "n_sz", "n_ctl", "n_genes", "L2_adds", "L3_adds", "L4_adds", "nsv_L5a", "nsv_L5b"] + [f"nq05_{l}" for l in LEVEL_DESC],
            lv, [pl], ["datasets entering the series and the covariates each level adds"])
    ind = [r["dataset"] for r in rd(pi) if r["in_independent_set"] == "TRUE"]
    F.note("b", f"Covers the {len(lv)} datasets with sample-level data in fairness_sample_level (378 dataset pairs x 7 levels = {len(rows)} rows), NOT the 64 published studies. "
                f"The gene z are over the independent set of {len(ind)} donor-separate datasets ({', '.join(ind)}). "
                f"b_pairs.tsv, b_recurrent_genes.tsv written {RUN} 16:05 and carry the current 176 genes ({len(gout)} genes here); "
                "the 47-study run's results (report/report.md item 19, figure 'generated 2026-10-01') are superseded for the gene panel. "
                "Rerunning on the 64 studies is not possible (the series needs sample-level expression and covariates, which exist for these 28 datasets only). "
                "BLAS_SENSITIVITY.md in that folder: medians and adjusted-vs-unadjusted counts are stable across linear-algebra libraries; L5a/L5b of single pairs are not.")
    F.note("b", "The request's 'z_all_covariates' is taken as L3 to match the report legend (PMI, pH, RIN, agonal factor); if L4 or L5b is meant, use the z_L4 / z_L5b columns.")


def fig7_c(F):
    units = [x for x in rd(P_UNITS) if x["scenario"] == SCENARIO and x["modelled"].lower().startswith("y")]
    p_h2 = os.path.join(AN, "h2", "h2_studies.tsv")
    p_old = os.path.join(AN, "ph_by_adjustment_2026-10-06", "ph_by_adjustment.tsv")
    h2 = {x["pmid"]: x for x in rd(p_h2)}
    z = {r["gene"]: float(r["Z"]) for r in rd(P_CTLPH) if r["factor"] == "pH" and r["Z"]}
    pl = json.load(open(P_PRIM))["studies"]
    handling = {"adjusted": "adjusted", "sensitivity": "tested", "matched_or_compared": "matched", "reported_only": "reported",
                "not_mentioned": "not_mentioned"}
    old = {r["pmid"]: r for r in rd(p_old)}
    out, diffs = [], 0
    for u in units:
        k = u["representative"]
        raw = h2[k]["cov_ph"] or "not_mentioned"
        genes = {}
        for ln in pl[k]["primary"]:
            for g, d in R.read_list(ln).items():
                genes[g] = d
        vals = [-{"up": 1, "down": -1}.get(d, 0) * z[g] for g, d in genes.items() if {"up": 1, "down": -1}.get(d, 0) and g in z]
        score = round(statistics.mean(vals), 3) if vals else ""
        o = old.get(k)
        if o and (o["score"] != str(score) and o["score"] != ""):
            diffs += 1
        out.append({"study_id": k, "label": u["label"], "ph_handling": handling.get(raw, "not_applicable"), "signed_ph_score": score,
                    "n_genes": len(genes), "n_scored": len(vals), "scored_ge10": "yes" if len(vals) >= 10 else "no",
                    "ph_status_in_h2": raw})
    F.write("fig07_c_ph_score_by_adjustment.tsv", "c",
            ["study_id", "label", "ph_handling", "signed_ph_score", "n_genes", "n_scored", "scored_ge10", "ph_status_in_h2"], out,
            [p_h2, P_CTLPH, P_UNITS, P_PRIM, os.path.join(AN, "h2", "cov_corrections.tsv")],
            ["signed_ph_score = mean over the study's primary-list genes with a direction and a control-pH Z of -(SZ sign in the study) x Z_pH; positive = the list moves with low pH "
             "the way it moves in SZ (convention of ph_by_adjustment_2026-10-06/ph_by_adjustment.py, reimplemented here against the 2026-10-08 combined_controls_pH.tsv)",
             "ph_handling from h2/h2_studies.tsv cov_ph (model reading of the main text with the hand corrections of h2/cov_corrections.tsv): "
             "adjusted; tested = sensitivity analysis; matched = matched/compared groups; reported = reported only; not_mentioned. Studies with < 10 scored genes are listed (scored_ge10 = no) but not tested"])
    cnt = {}
    for r in out:
        cnt[r["ph_handling"]] = cnt.get(r["ph_handling"], 0) + 1
    F.note("c", f"{len(out)} studies; handling counts {cnt} (report/report.md line 255: 13 adjusted, 6 tested, 11 matched, 7 reported, 27 not mentioned). "
                f"The existing ph_by_adjustment_2026-10-06/ph_by_adjustment.tsv was built before the 2026-10-08 combined_controls_pH.tsv; scores recomputed here, "
                f"{diffs} of {len(out)} scores differ from that older file (0 = identical to the third decimal). The existing file was not edited. "
                "Median-by-group figures in report.md (0.45 to 0.54 vs 0.86; p 0.87) come from the older scores; not rechecked here.")


def main():
    F2 = Figure("fig02")
    for fn in (fig2_a, fig2_b, fig2_c, fig2_d, fig2_e):
        fn(F2)
    F2.finish()
    F7 = Figure("fig07")
    for fn in (fig7_a1, fig7_a2, fig7_a3, fig7_b, fig7_c):
        fn(F7)
    F7.finish()


if __name__ == "__main__":
    main()
