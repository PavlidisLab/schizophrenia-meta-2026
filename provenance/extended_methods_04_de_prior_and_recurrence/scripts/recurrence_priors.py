"""Do the recurrent genes, and the individual published lists, sit high in the DE priors and in expression level?

Hypergeometric enrichment and AUC, set-level and per list as in Crow et al. 2019. Scenario 'published, brain' of recurrence.py (the modelled lists, rule A primary lists,
each restricted to the genes it tested).

References (each tested on its own):
  - DE priors in PRIOR_DIR (PRIOR_FILES): the global prior and the postmortem-brain prior with and without the SZ banks;
    scores n_de_total and de_prior_norm, joined to HGNC symbols by recurrence.de_prior_table() (NCBI gene ID, symbol
    fallback), the join recurrence.py uses;
  - expression level: LIBD BrainSeq DLPFC mean log CPM, hippocampus where DLPFC has none (recurrence.expression_levels(),
    recurrence.EXPR_FILE).
Higher score = more often DE / higher expression.

1. Set tests. Sets from recurrence_genes.tsv: q_bh < 0.1 (uniform null), q_expr_bh < 0.1 (expression-matched null) and
   genes in 5 or more of the lists (n_lists_observed >= 5, "5+ studies"). Universe: genes measured by at least one of the
   lists that the reference covers. Hypergeometric test of the set in the reference's top 10% and top 25% (top = genes
   whose score is at least the score at rank ceil(f x N), so ties can make the top larger than f x N; its size is
   reported), both tails; Mann-Whitney of the set against the rest of the universe, with AUC = U / (n_set x n_rest).
2. Per list, as Crow et al. 2019: AUROC of the reference ranking for membership of the list among the genes the list
   tested (and the reference covers). The same for the non-SZ comparison lists of de_prior_lists/checks_2026-09-29/
   (Gemma bipolar disorder and major depression postmortem brain lists at FDR 0.05, and their top-166 sensitivity lists).

The prior files follow recurrence.py's DE_PRIOR_SOURCE: our own 2-fold Gemma build
(de_prior_lists/our_prior/); DE_PRIOR_SOURCE=global gives the external files of de_prior_lists/de_prior_ranks/, built at
|log2FC| >= 2 (4-fold) where 2-fold was meant (Neera, 2026-09-29), whose rows are marked provisional.

    python3 analyses/recurrent_genes_vs_de_prior/recurrence_priors.py            # or PRIOR_DIR=<dir> python3 ...
    -> set_tests.tsv, set_auc.tsv, list_auroc.tsv, list_auroc_summary.tsv, recurrence_priors_meta.json, recurrence_priors_report.md
"""
import csv, datetime, hashlib, json, math, os, sys
import numpy as np
from scipy.stats import hypergeom, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
sys.path.insert(0, AN)
import recurrence as R  # noqa: E402

SCEN = "published, brain"
# The DE-prior files: one directory, three files, following recurrence.py's DE_PRIOR_SOURCE (default "ours" since
# 2026-09-30: our own 2-fold Gemma build in de_prior_lists/our_prior/; "global": the external 4-fold build). PRIOR_DIR=<dir>
# points at another copy of the same three files.
if R.DE_PRIOR_SOURCE == "ours":
    PRIOR_DIR = os.environ.get("PRIOR_DIR", os.path.join(AN, "de_prior_lists", "our_prior"))
    PRIOR_FILES = {"global": "1_global_our_prior.tsv", "brain postmortem": "2_brain_postmortem_our_prior.tsv",
                   "brain postmortem, no SZ banks": "2b_brain_postmortem_no_sz_banks_our_prior.tsv"}
    PRIOR_PROVISIONAL = ""
else:
    PRIOR_DIR = os.environ.get("PRIOR_DIR", os.path.join(AN, "de_prior_lists", "de_prior_ranks"))
    PRIOR_FILES = {"global": "1_global_de_prior.tsv", "brain postmortem": "2_brain_postmortem_de_prior.tsv",
                   "brain postmortem, no SZ banks": "2b_brain_postmortem_no_sz_banks_de_prior.tsv"}
    PRIOR_PROVISIONAL = ("provisional: the prior files were built at |log2FC| >= 2 (4-fold) where 2-fold was meant "
                         "(Neera, 2026-09-29)")
TOP = (0.10, 0.25)
CHECKS = os.path.join(AN, "de_prior_lists", "checks_2026-09-29")
# Crow et al. 2019 [@crowPredictabilityHumanDifferential2019], PNAS 116(13), verbatim, from the PDF in Zotero (62VCXIKF)
CROW = [("p. 6493", "We determine how well we recover DE genes for each dataset using the ranking from the compendium-wide list "
                    "(leaving out the study to be predicted) and report performance as the AUROC."),
        ("p. 6498", "For each experiment, DE results were thresholded to only include genes with an absolute log2 fold change >2 "
                    "and those with FDR of <0.05."),
        ("p. 6493", "Thus, 0.5 is random, 1.0 is perfect, and 0.7–0.8 is generally considered high performance. We find that the "
                    "global DE prior has remarkably high performance (mean AUROC, 0.83 ± 0.1; Fig. 2B; all scores, Dataset S1)."),
        ("p. 6498", "Reducing the threshold to include all genes with log2 fold change of >|1| at FDR < 0.05 would have increased "
                    "the total number of datasets in the compendium from 635 to 719 while slightly decreasing the performance of "
                    "the prior (AUROC, 0.79 ± 0.1).")]
sha = lambda p, n=12: hashlib.sha256(open(p, "rb").read()).hexdigest()[:n]


def references():
    """[(ref_id, label, {gene: score}, source path, provisional)]"""
    out = []
    keep = R.DE_PRIOR
    try:
        for name, fn in PRIOR_FILES.items():
            path = os.path.join(PRIOR_DIR, fn)
            R.DE_PRIOR = path  # de_prior_table() reads the file R.DE_PRIOR names (same columns in every group file)
            t = R.de_prior_table()
            prov = bool(PRIOR_PROVISIONAL)
            out.append((f"{name} prior, n_de_total", f"{name} DE prior, n_de_total", {g: v[0] for g, v in t.items()}, path, prov))
            out.append((f"{name} prior, de_prior_norm", f"{name} DE prior, de_prior_norm", {g: v[1] for g, v in t.items()}, path, prov))
    finally:
        R.DE_PRIOR = keep
    out.append(("expression", "expression level (LIBD DLPFC mean log CPM, else hippocampus)", R.expression_levels(), R.EXPR_FILE, False))
    return out


def auc(pos, neg):
    """AUC = U / (n_pos n_neg) (ties count one half) and the two-sided Mann-Whitney p."""
    if not len(pos) or not len(neg):
        return float("nan"), float("nan")
    # asymptotic p: with an untied score (our prior's de_prior_norm) and a list of 8 genes or fewer, method="auto" picks the
    # exact distribution against ~30,000 other genes, which ran for over 40 minutes and was killed (2026-09-30)
    r = mannwhitneyu(pos, neg, alternative="two-sided", method="asymptotic")
    return r.statistic / (len(pos) * len(neg)), r.pvalue


def read_symbols(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path, encoding="utf-8") if l.strip() and not l.startswith("#")]
    return {r[0] for r in rows[1:]} if rows and rows[0][0] == "symbol" else {r[0] for r in rows}


def summarize(vals):
    v = np.array([x for x in vals if not math.isnan(x)])
    if not len(v):
        return {"n_lists": 0}
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return {"n_lists": len(v), "median": round(float(med), 3), "q1": round(float(q1), 3), "q3": round(float(q3), 3),
            "mean": round(float(v.mean()), 3), "sd": round(float(v.std(ddof=1)), 3) if len(v) > 1 else "",
            "n_above_0.5": int((v > 0.5).sum()), "min": round(float(v.min()), 3), "max": round(float(v.max()), 3)}


def main():
    sel, M, lists, units = R.scenario_selection(SCEN)
    Lr = {k: set(lists[k]) & M[k] for k in sel}
    measured = set().union(*(M[k] for k in sel))
    rg_path = os.path.join(AN, "recurrence_genes.tsv")
    rg = [r for r in csv.DictReader((l for l in open(rg_path) if not l.startswith("#")), delimiter="\t") if r["scenario"] == SCEN]
    sets = {"study tier q<0.1": {r["gene"] for r in rg if float(r["q_bh"]) < 0.1},
            "expression-matched q<0.1": {r["gene"] for r in rg if r.get("q_expr_bh") and float(r["q_expr_bh"]) < 0.1},
            "5+ studies": {r["gene"] for r in rg if int(r["n_lists_observed"] or 0) >= 5}}
    # the lists here and recurrence_genes.tsv must describe the same scenario
    obs = {g: sum(1 for k in sel if g in Lr[k]) for g in measured}
    if {g for g, n in obs.items() if n >= 5} != sets["5+ studies"]:
        raise SystemExit("recurrence_genes.tsv n_lists_observed does not match the lists of recurrence.scenario_selection; rerun recurrence.py")
    refs = references()
    set_rows, auc_rows, list_rows = [], [], []
    for ref_id, label, score, path, prov in refs:
        U = sorted(g for g in measured if g in score)
        x = np.array([score[g] for g in U], dtype=float)
        N = len(U)
        order = np.sort(x)[::-1]
        for sname, S in sets.items():
            inS = np.array([g in S for g in U])
            n = int(inS.sum())
            a, p = auc(x[inS], x[~inS])
            auc_rows.append({"gene_set": sname, "reference": ref_id, "n_universe": N, "n_set": len(S), "n_set_in_universe": n,
                             "median_set": round(float(np.median(x[inS])), 4) if n else "", "median_universe": round(float(np.median(x)), 4),
                             "auc": round(a, 4), "p_mannwhitney_two_sided": f"{p:.3g}", "provisional": PRIOR_PROVISIONAL if prov else ""})
            for f in TOP:
                thr = order[math.ceil(f * N) - 1]
                top = x >= thr
                K, k = int(top.sum()), int((top & inS).sum())
                e = n * K / N
                set_rows.append({"gene_set": sname, "reference": ref_id, "top": f"{f:.0%}", "threshold": round(float(thr), 4),
                                 "n_universe": N, "n_top": K, "n_set_in_universe": n, "observed": k, "expected": round(e, 2),
                                 "fold": round(k / e, 2) if e else "", "p_enrichment": f"{hypergeom.sf(k - 1, N, K, n):.3g}",
                                 "p_depletion": f"{hypergeom.cdf(k, N, K, n):.3g}", "provisional": PRIOR_PROVISIONAL if prov else ""})
        # per list (Crow et al. 2019): SZ lists of the scenario
        for k in sel:
            Uk = [g for g in M[k] if g in score]
            pos = np.array([score[g] for g in Uk if g in Lr[k]], dtype=float)
            neg = np.array([score[g] for g in Uk if g not in Lr[k]], dtype=float)
            a, p = auc(pos, neg)
            list_rows.append({"group": "SZ published brain", "list": "; ".join(units["studies"][k]["primary"]),
                              "study": units["studies"][k]["label"], "reference": ref_id, "n_tested": len(Uk), "n_in_list": len(pos),
                              "auroc": round(a, 4), "p_mannwhitney_two_sided": f"{p:.3g}", "in_global_prior_matrix": "",
                              "provisional": PRIOR_PROVISIONAL if prov else ""})
        # non-SZ comparison: Gemma BD and MDD postmortem brain lists (checks_2026-09-29)
        meta = [r for r in csv.DictReader((l for l in open(os.path.join(CHECKS, "check1_lists.tsv")) if not l.startswith("#")), delimiter="\t")]
        for r in meta:
            uni = read_symbols(os.path.join(CHECKS, "universes", r["list_id"] + ".txt"))
            for kind, d in (("FDR 0.05", "lists"), ("top 166", "lists_top")):
                fp = os.path.join(CHECKS, d, r["list_id"] + ".tsv")
                if not os.path.exists(fp):
                    continue
                L = read_symbols(fp) & uni
                if not L:
                    continue
                Uk = [g for g in uni if g in score]
                pos = np.array([score[g] for g in Uk if g in L], dtype=float)
                neg = np.array([score[g] for g in Uk if g not in L], dtype=float)
                a, p = auc(pos, neg)
                list_rows.append({"group": f"{r['disorder']} Gemma {kind}", "list": r["list_id"], "study": r["short_name"], "reference": ref_id,
                                  "n_tested": len(Uk), "n_in_list": len(pos), "auroc": round(a, 4), "p_mannwhitney_two_sided": f"{p:.3g}",
                                  "in_global_prior_matrix": r["in_prior_matrix"], "provisional": PRIOR_PROVISIONAL if prov else ""})
    summ = []
    groups = sorted({r["group"] for r in list_rows}, key=lambda g: (not g.startswith("SZ"), g))
    for ref_id, *_ in refs:
        for g in groups:
            v = [r["auroc"] for r in list_rows if r["group"] == g and r["reference"] == ref_id]
            summ.append({"group": g, "reference": ref_id, "subset": "all lists"} | summarize(v))
            if not g.startswith("SZ"):
                v = [r["auroc"] for r in list_rows if r["group"] == g and r["reference"] == ref_id and r["in_global_prior_matrix"] == "False"]
                summ.append({"group": g, "reference": ref_id, "subset": "dataset not in the global prior's input matrix"} | summarize(v))
    stamp = datetime.date.today().isoformat()
    ins = {os.path.relpath(p, AN): sha(p) for p in [rg_path, os.path.join(AN, "primary_lists.json"), os.path.join(AN, "universes.json"),
                                                    R.LINEAGE, R.EXPR_FILE, os.path.join(CHECKS, "check1_lists.tsv")]
           + [os.path.join(PRIOR_DIR, f) for f in PRIOR_FILES.values()]}
    head = [f"# script: analyses/recurrent_genes_vs_de_prior/recurrence_priors.py, run {stamp}; scenario '{SCEN}', {len(sel)} lists, {len(measured)} genes measured",
            "# inputs sha256: " + "; ".join(f"{k} {v}" for k, v in ins.items()),
            f"# DE prior: {R.DE_PRIOR_SOURCE} ({os.path.relpath(PRIOR_DIR, AN)}/)" + (f"; DE-prior rows are {PRIOR_PROVISIONAL}" if PRIOR_PROVISIONAL else "")]

    def write(name, rows):
        cols = list(dict.fromkeys(c for r in rows for c in r))
        with open(os.path.join(HERE, name), "w", newline="") as fh:
            fh.write("\n".join(head) + "\n")
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
            w.writeheader(); w.writerows(rows)
    write("set_tests.tsv", set_rows)
    write("set_auc.tsv", auc_rows)
    write("list_auroc.tsv", list_rows)
    write("list_auroc_summary.tsv", summ)
    json.dump({"generated": stamp, "script": "analyses/recurrent_genes_vs_de_prior/recurrence_priors.py", "scenario": SCEN, "n_lists": len(sel),
               "genes_measured": len(measured), "gene_sets": {k: len(v) for k, v in sets.items()}, "prior_dir": os.path.relpath(PRIOR_DIR, AN),
               "prior_files": PRIOR_FILES, "prior_provisional": PRIOR_PROVISIONAL, "top_fractions": list(TOP), "inputs_sha256": ins,
               "crow_2019_quotes": [{"page": p, "quote": q} for p, q in CROW]},
              open(os.path.join(HERE, "recurrence_priors_meta.json"), "w"), indent=1)
    report(stamp, sel, measured, sets, refs, set_rows, auc_rows, summ, ins)


def report(stamp, sel, measured, sets, refs, set_rows, auc_rows, summ, ins):
    ref_ids = [r[0] for r in refs]
    short = {"study tier q<0.1": "44-type set (study tier q < 0.1)", "expression-matched q<0.1": "expression-matched q < 0.1", "5+ studies": "5+ studies"}
    L = ["---", "type: detail", f"modified: {stamp}", "---", "", "> [!note] Paul", ">", "",
         "# Recurrent genes and published lists against the DE priors and expression level", "",
         f"Built by `recurrence_priors.py` ({stamp}); edit the script, not this file. Scenario '{SCEN}' of `recurrence.py`: {len(sel)} lists, "
         f"{len(measured)} genes measured by at least one. Sets: " + "; ".join(f"{k} ({len(v)} genes)" for k, v in sets.items()) + ". "
         "Universe for each test: the measured genes the reference covers (its size is in each row). References: "
         + "; ".join(f"`{r[0]}`" for r in refs) + f". Prior files: `{os.path.relpath(PRIOR_DIR, AN)}/` ("
         + ", ".join(PRIOR_FILES.values()) + ").", "",
         (f"**Every DE-prior number here is {PRIOR_PROVISIONAL}.** " if PRIOR_PROVISIONAL else
          f"DE prior: `{R.DE_PRIOR_SOURCE}` (recurrence.py DE_PRIOR_SOURCE). ")
         + "Rerun: `python3 analyses/recurrent_genes_vs_de_prior/recurrence_priors.py` "
         "(or with `DE_PRIOR_SOURCE=global`, or `PRIOR_DIR=<dir>` for another copy of the three files).", "",
         "## Set AUC (Mann-Whitney of the set against the rest of the universe)", "",
         "| reference | " + " | ".join(f"{s} ({len(v)})" for s, v in sets.items()) + " |", "|---|" + "---:|" * len(sets)]
    for ref in ref_ids:
        cells = []
        for s in sets:
            r = next(x for x in auc_rows if x["gene_set"] == s and x["reference"] == ref)
            cells.append(f"{r['auc']:.3f} (p {r['p_mannwhitney_two_sided']}; {r['n_set_in_universe']} of {r['n_universe']})")
        L.append(f"| {ref} | " + " | ".join(cells) + " |")
    L += ["", "## Hypergeometric: the set in the reference's top 10% and top 25%", "",
          "| reference | top | set | universe | top genes | set genes in universe | observed | expected | fold | p enrichment | p depletion |",
          "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in set_rows:
        L.append(f"| {r['reference']} | {r['top']} | {r['gene_set']} | {r['n_universe']} | {r['n_top']} | {r['n_set_in_universe']} | {r['observed']} | "
                 f"{r['expected']} | {r['fold']} | {r['p_enrichment']} | {r['p_depletion']} |")
    L += ["", "Top = genes whose score is at least the score at rank ceil(f × N); ties (n_de_total has many) make the top larger than f × N.", "",
          "## Per-list AUROC (as Crow et al. 2019)", "",
          "For each list, the AUROC of the reference ranking for list membership among the genes the list tested and the reference covers. "
          "Non-SZ comparison: Gemma bipolar disorder (BD) and major depression (MDD) postmortem brain lists of "
          "`../de_prior_lists/checks_2026-09-29/` (FDR 0.05 lists with at least one gene, and top-166 sensitivity lists). The SZ series are "
          "held out of the priors (`../de_prior_lists/README.md`, `held_out.tsv`); some BD and MDD datasets are in the prior's input, so their "
          "rows are also shown without them (`in_global_prior_matrix`, from `check1_lists.tsv`, refers to the global prior).", "",
          "| group | reference | subset | lists | median (IQR) | mean ± SD | lists with AUROC > 0.5 |", "|---|---|---|---:|---|---|---:|"]
    for r in summ:
        if not r.get("n_lists"):
            continue
        L.append(f"| {r['group']} | {r['reference']} | {r['subset']} | {r['n_lists']} | {r['median']:.3f} ({r['q1']:.3f}-{r['q3']:.3f}) | "
                 f"{r['mean']:.3f} ± {r['sd'] if r['sd'] == '' else format(r['sd'], '.3f')} | {r['n_above_0.5']} |")
    L += ["", "Crow et al. 2019 [@crowPredictabilityHumanDifferential2019], verbatim, for comparison:", ""]
    L += [f"> {q} ({p})" for p, q in CROW]
    L += ["", "Files: `set_auc.tsv`, `set_tests.tsv`, `list_auroc.tsv` (one row per list and reference), `list_auroc_summary.tsv`, "
          "`recurrence_priors_meta.json` (inputs sha256: " + ", ".join(f"{k} {v}" for k, v in ins.items()) + ")."]
    open(os.path.join(HERE, "recurrence_priors_report.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
