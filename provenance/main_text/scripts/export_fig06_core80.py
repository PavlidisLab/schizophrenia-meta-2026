#!/usr/bin/env python3
"""Per-panel tidy exports for manuscript Figure 6 (recurrent genes: direction, sources, co-occurrence, DE prior).

Request: MAB_TO_RA_2026_10_08_FRESH_FIGURES_PANEL_DATA.md (Fig 6).
Population: the 64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrent set = the 176
study-tier genes at q < 0.1 (uniform null).

Reads existing analysis outputs only (no model fits, no permutation reruns, no network). Writes new files to
figure_data/fig06/, each TSV with '#' provenance lines, plus _meta.json and STATUS.md.
One function per panel; rerun with

    python3 analyses/figure_exports/export_fig06.py
"""
import csv
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
AN = "core80_work_2026-10-16/P11/b/T2/analyses"      # core-80: private tree (was os.path.dirname(HERE))
REPO = os.path.dirname(AN)
import sys  # noqa: E402
sys.path.insert(0, AN)
sys.path.insert(0, HERE)
import _run_stamp  # noqa: E402
import recurrence as R  # noqa: E402

RUN = _run_stamp.run_date()
SCRIPT = "core80_work_2026-10-16/P11/b/f6/export_fig06_core80.py (copy of analyses/figure_exports/export_fig06.py; paths and population stamp changed, panel d dropped)"
OUT = os.environ["FIGURE_DATA_OUT"]   # core-80: no default
POP = _run_stamp.population(len(R.scenario_selection("published, brain")[0]))
# Panel a and the 'frozen analysis set' row of panel b use the frozen 209 (FIXED_STUDY_SET=analyses/confounds/recurrent_sets.tsv when the
# satellites are run); every 'q<0.1' row is the in-place recurrence_genes.tsv set on the current base.
POP += "; panel a = the study-tier q<0.1 set of the in-place recurrence_genes.tsv (383), rows named 'q<0.1' = the in-place q_bh / q_expr_bh sets on this base"

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


P_DIRG = os.path.join(AN, "recurrence_direction", "direction_genes.tsv")
P_DIRS = os.path.join(AN, "recurrence_direction", "direction_summary.tsv")
P_STRATA = os.path.join(AN, "recurrence_sources", "strata.tsv")
P_THEME = os.path.join(AN, "recurrence_sources", "theme_strata.tsv")
COOC_DIR = os.environ.get("COOC_DIR") or os.path.join(AN, "recurrence_sources")   # panel d inputs; COOC_DIR = analyses/cooccurrence_64_2026-10-13/out/cooccurrence_209 for the 64-list rerun
P_COO = os.path.join(COOC_DIR, "cooccurrence_pairs.tsv")
P_COOC = os.path.join(COOC_DIR, "cooccurrence_pairs_cohort.tsv")
P_COOM = os.path.join(COOC_DIR, "cooccurrence_meta.json")
L7 = "rerun_on_80_studies/L07/recurrence_priors"
P_AUC = os.path.join(L7, "set_auc.tsv")
P_TOP = os.path.join(L7, "set_tests.tsv")
P_LIST = os.path.join(L7, "list_auroc_summary.tsv")
P_GENES = os.path.join(AN, "recurrence_genes.tsv")


def fig6_a(F):
    rows = rd(P_DIRG)
    out = [{"gene": r["gene"], "n_studies": r["n_studies"], "n_up": r["up"], "n_down": r["down"], "n_mixed": r["both"],
            "n_nodir": r["no_direction"], "majority": r["majority"], "share_majority": r["share_majority"],
            "sign_test_p": r["sign_test_p"]} for r in rows]
    bad = [r["gene"] for r in rows if int(r["up"]) + int(r["down"]) + int(r["both"]) + int(r["no_direction"]) != int(r["n_studies"])]
    F.write("fig06_a_gene_direction.tsv", "a", ["gene", "n_studies", "n_up", "n_down", "n_mixed", "n_nodir", "majority", "share_majority", "sign_test_p"],
            out, [P_DIRG, P_GENES],
            ["one row per gene of the frozen analysis set (analyses/confounds/recurrent_sets.tsv, study_brain), sorted as in the source; n_mixed = studies reporting the gene in both directions ('both' in the source); "
             "n_nodir = hit with no stated direction"])
    F.note("a", f"{len(out)} genes (frozen analysis set; expected 209); up + down + mixed + no-direction equals n_studies for every gene" if not bad
           else f"MISMATCH up + down + mixed + nodir != n_studies for {bad}")
    tot = {k: sum(int(r[k]) for r in rows) for k in ("up", "down", "both", "no_direction", "n_studies")}
    F.note("a", f"column sums: up {tot['up']}, down {tot['down']}, mixed {tot['both']}, no direction {tot['no_direction']}, n_studies {tot['n_studies']}.")
    return rows


def fig6_b(F):
    rows = rd(P_DIRS)
    out = [dict(r) for r in rows]
    F.write("fig06_b_direction_summary.tsv", "b",
            ["gene_set", "genes", "list_hits_with_direction", "in_majority_direction", "share_majority", "share_majority_if_random",
             "genes_unanimous", "genes_split_40pct_or_more", "genes_sign_test_p05"], out, [P_DIRS],
            ["share_majority_if_random = expected share in the majority direction if each list's direction were random (source column); "
             "rows: frozen analysis set (209), study tier q<0.1 uniform null, expression-matched q<0.1, 3 or more studies, 5+ studies"])
    F.note("b", "report/report.md lines 172 and 187: 1,315 of 1,540 hits (85%) vs 63% random, 48 unanimous, 10 split; 99 expression-matched 89% vs 63%; 5,477 genes in 3 or more lists 80% vs 70%: all match direction_summary.tsv. "
                "Note the report text and legends say 'the 99 expression-matched' and '5,477 genes in 3 or more lists'; the 5+ studies row (1,478 genes, 81% vs 67%) is not quoted in the report's direction bullet.")
    return rows


def fig6_c(F):
    out = []
    for r in rd(P_STRATA):
        out.append({"gene_set": r["gene_set"], "n_genes": r["n_genes"], "stratum_type": r["stratum"].split(" ")[0] if r["stratum"].startswith(("era", "tech")) else "size",
                    "stratum": r["stratum"], "n_studies": r["studies"], "list_slots": r["list_slots"], "share_slots": r["share_of_slots"],
                    "hits": r["hits"], "share_hits": r["share_of_hits"], "ratio": r["ratio"], "perm_p": "", "q_bh": ""})
    for r in rd(P_THEME):
        out.append({"gene_set": r["gene_set"], "n_genes": "", "stratum_type": "theme", "stratum": r["theme"], "n_studies": r["n_studies"],
                    "list_slots": "", "share_slots": r["share_of_slots"], "hits": "", "share_hits": r["share_of_hits"], "ratio": r["ratio"],
                    "perm_p": r["p_perm"], "q_bh": r["q_bh"]})
    F.write("fig06_c_strata.tsv", "c",
            ["gene_set", "n_genes", "stratum_type", "stratum", "n_studies", "list_slots", "share_slots", "hits", "share_hits", "ratio", "perm_p", "q_bh"],
            out, [P_STRATA, P_THEME],
            ["stratum_type size / era / tech from strata.tsv (no permutation p is stored for these: perm_p empty); stratum_type theme from theme_strata.tsv "
             "(perm_p = 5000 permutations of theme membership across studies, seed 20260929; q_bh = BH over themes within a gene set)",
             "share_slots = share of list slots, share_hits = share of the set's list hits; ratio = share_hits / share_slots"])
    n_theme = sum(1 for r in out if r["stratum_type"] == "theme")
    F.note("c", f"{len(out)} rows: size/era/tech strata from strata.tsv for 5 gene sets, {n_theme} theme rows with perm_p. strata.tsv stores no permutation p for the size, "
                "era or technology strata; only the themes have perm_p (recurrence_sources.py, NPERM_THEME = 5000). A permutation p for size or era would need a new "
                "run of recurrence_sources.py with a permutation step added; not run.")
    F.note("c", "report/report.md lines 172-174 vs the files. MATCH: large studies 17% of slots / 16% of hits (ratio 0.91), expression-matched 16% (0.95); era 2001-2012 20% of hits, 18% of slots; "
                "37 of 176 and 18 of 99 genes in no large study's list (recurrence_sources_report.md); immune 26 studies 45% / 49%, 1.09, p 0.24; immune+stress+microglia+astrocytes 28 studies 46% / 51%, 1.11, p 0.18; "
                "expression-matched immune 1.24, p 0.087; mitochondria 9 studies 24% / 22%, 0.91, p 0.71 (0.96, p 0.60 expression-matched); oligodendrocytes 5 studies 2.22, p 0.029, q 0.37; microglia 2.2, p 0.17, q 0.45.")
    F.note("c", "MISMATCH 1: the report names 8 large studies; strata.tsv has 9 (large = N >= 150; recurrence_sources_report.md lists Fromer 2016, Jaffe 2018, Collado-Torres 2019, Jaffe 2020, Benjamin 2022, Schulmann 2023, Ling 2024, Dean 2026 and Lee 2026 (PsychAD)). The report's list omits Lee 2026.")
    F.note("c", "MISMATCH 2: the report's 'robust' set is 8 genes (7% of hits, ratio 0.43, 2 of 8 in no large list; immune 1.47, p 0.040; mitochondria 22%, 0.92, p 0.61). The files' 'robust set (impact)' is 35 genes: large studies 7.2% of hits, ratio 0.42; immune 1.27, p 0.079; mitochondria 22%, 0.92, p 0.6129. Only the mitochondria numbers agree. The robust set is exported with the other sets; the figure need not use it.")


def _coo(F, path, unit, name, panel, tier_label):
    rows = rd(path)
    out = []
    zero = 0
    for r in rows:
        a, b, t, ob = int(r[f"{unit}_a"]), int(r[f"{unit}_b"]), int(r[f"{unit}_both_tested"]), int(r[f"{unit}_both"])
        zero += ob == 0
        out.append({"gene_a": r["gene_a"], "gene_b": r["gene_b"], f"{unit}_both_tested": t, f"{unit}_a": a, f"{unit}_b": b, f"{unit}_both": ob,
                    "expected_uniform": round(a * b / t, 4), "fold_uniform": round(ob / (a * b / t), 4) if a * b else "",
                    **({"expected_swap": r["expected_swap"]} if "expected_swap" in r else {}), "p_swap": r["p_swap"], "q_bh": r["q_bh"]})
    m = json.load(open(P_COOM))
    has_swap = "expected_swap" in rows[0]
    stale = m["inputs"]["primary_lists.json"] != sha(os.path.join(AN, "primary_lists.json"))[:12]
    F.write(name, panel, ["gene_a", "gene_b", f"{unit}_both_tested", f"{unit}_a", f"{unit}_b", f"{unit}_both", "expected_uniform", "fold_uniform"] + (["expected_swap"] if has_swap else []) + ["p_swap", "q_bh"],
            out, [path, P_COOM],
            [f"tier: {tier_label}; n_randomizations {m['n_randomizations']} (swaps per randomization {m['swaps_per_randomization']}, seed {m['seed']})",
             "expected_uniform = a * b / both_tested (computed here from the columns, the same formula cooccurrence_figure.py uses); fold_uniform = both / expected_uniform",
             (f"refreshed on the 64 base: analyses/cooccurrence_64_2026-10-13/cooccurrence_64.py, {m['n_lists']} lists, primary_lists.json {m['inputs']['primary_lists.json']}, recurrence_genes.tsv {m['inputs']['recurrence_genes.tsv']}, universes.json {m['inputs']['universes.json']}, {m['n_genes_44']} genes ({m.get('gene_set', 'q<0.1 set')}); {len(rows)} pairs" if not stale else
              "DATA NOT REFRESHED ON THE 64 BASE: cooccurrence.py (remote-only) was last run on recurrence_genes.tsv sha256 " + m["inputs"]["recurrence_genes.tsv"] + " (" + str(m["n_lists"]) + " lists, primary_lists.json " + m["inputs"]["primary_lists.json"] + ", " + str(m["n_genes_44"]) + " genes; " + str(len(rows)) + " pairs); the current base is 64 lists, primary_lists.json 5bdce9b94e15, recurrence_genes.tsv " + sha(P_GENES)[:12]),
             ("expected_swap = mean over the randomizations of the swap-null count of lists (clusters) holding both genes; expected_uniform and fold_uniform use the a * b / both_tested formula" if has_swap else "expected_swap (swap-null mean co-occurrence per pair) is NOT stored: cooccurrence.py keeps only the exceedance count per pair")])
    return rows, zero


def fig6_d(F):
    r1, z1 = _coo(F, P_COO, "lists", "fig06_d_cooccurrence_study_tier.tsv", "d", "study tier")
    r2, z2 = _coo(F, P_COOC, "clusters", "fig06_d_cooccurrence_cohort_tier.tsv", "d", "cohort tier")
    m = json.load(open(P_COOM))
    F.note("d", f"study tier {len(r1)} pairs ({z1} with zero co-occurrence), cohort tier {len(r2)} pairs ({z2} with zero); n_randomizations {m['n_randomizations']} per tier "
                f"(swaps per randomization {m['swaps_per_randomization']}, seed {m['seed']}); pairs at q < 0.1: {m['pairs_q_lt_0.1']} (study), {m['pairs_q_lt_0.1_cohort']} (cohort).")
    if "expected_swap" in r1[0]:
        F.note("d", f"refreshed on the 64 base ({m['n_lists']} lists, {m['n_genes_44']} genes, {len(r1)} pairs): analyses/cooccurrence_64_2026-10-13/cooccurrence_64.py; expected_swap stored (mean over {m['n_randomizations']} randomizations of the swap-null count of lists holding both genes). p_swap = (exceedances + 1) / (n_randomizations + 1), smallest attainable {1/(m['n_randomizations']+1):.2g}; smallest attainable q over {len(r1)} pairs = {len(r1)/(m['n_randomizations']+1):.2f}.")
    else:
        F.note("d", "expected_swap NOT stored anywhere. cooccurrence.py (run_chain) accumulates only ge[a][b], the number of randomizations with a count >= observed; the "
                    "randomized tables are discarded. p_swap = (ge + 1) / (n_randomizations + 1). To get expected_swap, run_chain would have to also accumulate the sum of "
                    "the randomized pair counts (Yn @ Yn.T) and divide by n_randomizations; that is a full rerun of cooccurrence.py (40,000 randomizations x 20,000 swaps, "
                    "minutes with the matrix-product version, two processes) and was not run. expected_uniform and fold_uniform are exported instead and are labelled as uniform.")
        # the input hash recorded in the co-occurrence TSV header versus the current recurrence_genes.tsv
        hdr = open(P_COO).readline()
        F.note("d", f"provenance: the co-occurrence TSVs were written 2026-10-08 09:53 from recurrence_genes.tsv sha256 {m['inputs']['recurrence_genes.tsv']}; the current "
                    f"recurrence_genes.tsv is {sha(P_GENES)[:12]} (rewritten 14:31). The gene set of the pairs equals the 176 genes of direction_genes.tsv (same symbols); "
                    "whether cooccurrence.py was rerun against the new file is not recorded.")
    F.note("d (report figure)", "analyses/recurrence_sources/cooccurrence_figure.py edited and analyses/recurrence_sources/fig_cooccurrence.svg regenerated (2026-10-08, report build not run). Before: log2 x-axis from 1/8 to 64 "
                "(data maximum fold 6.56 study, 8.5 cohort), the 273 zero-co-occurrence study-tier pairs clamped into the 1/8 bin, about 60% of the canvas empty at the right (the footer text set the bounding box). "
                "After: linear fold axis limited to the data range (0 to 6.75 study, 0 to 8.5 cohort, 0.25 bins), the 273 zero-co-occurrence pairs as their own labelled bar ('none', cohort tier 0), "
                "fold = observed / (a * b / tested) with no floor on the expectation (the old code floored it at 0.5, which changed the fold of 419 cohort-tier pairs; cohort median 2.33 before, 2.36 after; study median 2.14 both), "
                "Helvetica/Arial, svg.fonttype none, y grid only, normal-weight left-aligned titles, footer wrapped. report/report.md line 178 still says 'histogram, log2 x-axis' and 'generated 2026-10-03'; "
                "figures/fig_cooccurrence.pdf is built from the SVG by the report build and was not regenerated.")
    return r1, r2, z1, z2


def fig6_e(F):
    rows = rd(P_AUC)
    out = [{k: r[k] for k in ("gene_set", "reference", "n_universe", "n_set", "n_set_in_universe", "median_set", "median_universe", "auc", "p_mannwhitney_two_sided", "provisional")} for r in rows]
    F.write("fig06_e_de_prior_auc_by_set.tsv", "e",
            ["gene_set", "reference", "n_universe", "n_set", "n_set_in_universe", "median_set", "median_universe", "auc", "p_mannwhitney_two_sided", "provisional"], out, [P_AUC],
            ["AUC of the reference (DE prior or expression) ranking the set's genes above the other genes of the universe; 0.5 = chance; p_mannwhitney_two_sided of 0 is the file's underflow"])
    top = rd(P_TOP)
    F.write("fig06_e_de_prior_top_fraction.tsv", "e", list(top[0].keys()), top, [P_TOP], ["observed vs expected set genes in the top 10% / 25% of the reference"])
    lst = rd(P_LIST)
    F.write("fig06_e_de_prior_per_list_auroc.tsv", "e", list(lst[0].keys()), lst, [P_LIST],
            ["per-list AUROC of the prior over each list's genes, summarised by group (SZ published brain, BD, depression Gemma lists)"])
    F.note("e", "DE-prior AUC by set from recurrence_priors/set_auc.tsv (3 sets x 7 reference/measure combinations), plus top-fraction counts and per-list AUROC summaries. "
                "report/report.md line 182 quotes the global n_de_total prior only.")
    return rows


def main():
    F = Figure("fig06")
    for fn in (fig6_a, fig6_b, fig6_c, fig6_e):   # core-80: panel d (co-occurrence) is not in the current paper
        fn(F)
    return F


if __name__ == "__main__":
    F = main()
    F.finish()
