#!/usr/bin/env python3
"""Per-panel tidy exports for manuscript Figure 1 (b, c, d) and Figure 11 (a, plus the periphery tables).

Request: MAB_TO_RA_2026_10_08_FRESH_FIGURES_PANEL_DATA.md, Priority 2.
Population: the 64 published brain studies (recurrence.py scenario 'published, brain', study tier).

Reads existing analysis outputs only (no model fits, no permutation reruns, no network, no data/raw). Writes new files to
figure_data/fig01/ and fig11/, each TSV with '#' provenance lines, plus _meta.json and
STATUS.md per folder. One function per panel; rerun with

    python3 analyses/figure_exports/export_fig01_fig11.py
"""
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
AN = 'core80_work_2026-10-16/P11/b/T4/analyses'   # core-80: private tree (was dirname(HERE))
REPO = os.path.dirname(AN)
sys.path.insert(0, AN)
import _run_stamp  # noqa: E402
import recurrence as R  # noqa: E402  (the repo's own loading path; light)

RUN = _run_stamp.run_date()
SCRIPT = "core80_work_2026-10-16/P11/b/f1/export_fig01_core80.py (copy of analyses/figure_exports/export_fig01_fig11.py; panels b and c only; paths and population stamp changed)"
OUT = os.environ["FIGURE_DATA_OUT"]   # core-80: no default
SCENARIO = "published, brain"
POP = _run_stamp.population(len(R.scenario_selection("published, brain")[0]))

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


def oneline(s):
    return re.sub(r"\s+", " ", str(s)).strip() if s is not None else ""


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


P_LINEAGE = os.path.join(REPO, "corpus", "lineage.json")
P_TECH = os.path.join(AN, "technology.json")
P_STATS = os.path.join(AN, "corpus_stats.json")
P_H2 = os.path.join(AN, "h2", "h2_studies.tsv")
P_DS = os.path.join(AN, "subject_relations", "distinct_subjects.tsv")
P_UNITS = os.path.join(AN, "recurrence_units.tsv")
P_PRIM = os.path.join(AN, "primary_lists.json")
P_UNIV = os.path.join(AN, "universes.json")
P_NEG = os.path.join(AN, "neglected_data_2026-10-03", "results.tsv")
SMK = os.path.join(AN, "smoking_signature_2026-10-05")
PER = os.path.join(AN, "brain_genes_periphery_2026-10-01", "results")
LINEAGE_RUN_SHA = "6f4ba89306a7"

# ------------------------------------------------------------------ shared loads
_L = json.load(open(P_LINEAGE))
NODES = _L["nodes"]
EDGES = _L["edges"]
SEL = R.scenario_selection(SCENARIO)[0]          # 64 study keys (PMIDs)
SELSET = set(SEL)
CORPUS = [n for n in NODES if n.get("tech_tier")]  # study nodes of the corpus (sources / non-study nodes carry no tech_tier)
if os.environ.get("EXCLUDE_GROUP2"):
    # one-gene (candidate-gene) studies, lineage category '2', stay in the single-gene set and are not corpus studies
    CORPUS = [n for n in CORPUS if str(n.get("category")) != "2"]
TECH = json.load(open(P_TECH))["studies"]


def in_pop(n):
    return n["id"] in SELSET or n.get("pmid") in SELSET


def era(y, cut):
    return "" if y in (None, "") else ("old" if int(y) <= cut else "new")


def banks(n):
    b = n.get("banks") or []
    return ";".join(b)


# ================================================================== Figure 1
def fig1_b(F):
    h2 = {r["pmid"]: r for r in rd(P_H2)}
    out = []
    for n in CORPUS:
        t = TECH.get(n["id"]) or TECH.get(n.get("pmid")) or {}
        out.append({
            "study_id": n["id"], "pmid": n.get("pmid", ""), "label": n["label"], "year": n.get("year"),
            "era_2016": era(n.get("year"), 2015), "era_2020": era(n.get("year"), 2019),
            "relation": n.get("relation"), "modality": n.get("modality"),
            "tissue": n.get("tissue"), "tissue_short": n.get("tissue_short"), "brain_class": n.get("brain_class"),
            "bank": n.get("lane"), "bank_label": n.get("lane_label"), "bank_multi": banks(n),
            "technology": n.get("tech_tier"),
            "technology_json": t.get("technology", ""), "array_generation": t.get("array_generation", ""),
            "n_total": n.get("n"), "n_cases": n.get("cases"), "n_controls": n.get("controls"),
            "in_primary_population": "yes" if in_pop(n) else "no",
            "in_h2_reporting_table": "yes" if (n.get("pmid") in h2 or n["id"] in h2) else "no",
            "times_cited": n.get("times_cited"),
        })
    fields = list(out[0].keys())
    F.write("fig01_b_corpus_studies.tsv", "b", fields, out, [P_LINEAGE, P_TECH, P_H2, P_UNITS, P_PRIM, P_UNIV],
            ["one row per study node of corpus/lineage.json that carries a tech_tier (the corpus; 641 nodes in all, the rest are cited sources "
             "with no tech_tier); era_2016 old = year <= 2015, era_2020 old = year <= 2019 (the h2 splits)",
             "bank = lineage 'lane' (cohort/bank row), bank_label its description, bank_multi = lineage 'banks' list for studies spanning banks; "
             "technology = lineage tech_tier; technology_json / array_generation = analyses/technology.json where the study has a primary list",
             "n_total / n_cases / n_controls as in the lineage node (blank = not recorded there); in_primary_population = the study key "
             f"is one of the {len(SEL)} units of recurrence.scenario_selection('{SCENARIO}')",
             f"lineage.json sha256 now {sha(P_LINEAGE)[:12]}; the 2026-10-08 recurrence run used {LINEAGE_RUN_SHA}"])
    npop = sum(r["in_primary_population"] == "yes" for r in out)
    byera = Counter(r["era_2016"] for r in out)
    F.note("b", f"{len(out)} corpus studies (nodes with a tech_tier); {npop} flagged as primary population (expected {len(SEL)}). "
                f"corpus_stats.json also says 571 studies; report/report.md says 571 in most places but 570 in the Figure 1 caption "
                f"('full corpus, 570 studies: 202 old era, 368 new'): from the lineage nodes the 2016 split gives old {byera['old']}, new {byera['new']} "
                f"over {sum(byera.values())} studies with a year.")
    pk = [k for k in SEL if not any(n["id"] == k for n in CORPUS)]
    F.note("b", f"{len(pk)} of the 64 population keys are PMIDs whose lineage node id is a DOI or other id (study_id is the node id, pmid the PMID; match on pmid for these): "
                + ", ".join(pk))
    noN = sum(1 for r in out if r["n_cases"] in (None, ""))
    F.note("b", f"n_cases / n_controls are blank for {noN} of {len(out)} studies (not recorded in the lineage node); n_total is blank for "
                f"{sum(1 for r in out if r['n_total'] in (None, ''))}.")
    F.note("b", f"lineage.json was edited after the 2026-10-08 recurrence run: sha now {sha(P_LINEAGE)[:12]}, the run used {LINEAGE_RUN_SHA} "
                "(commit 2026-10-08 18:22 adds Muralidharan 2025). corpus_stats.json already records the new sha (df1d54611414); technology.json "
                f"and subject_relations/distinct_subjects.tsv record {LINEAGE_RUN_SHA}. "
                "The population flag does not depend on lineage.json (it comes from recurrence.scenario_selection).")
    # reporting completeness
    rows = rd(P_H2)
    if os.environ.get("EXCLUDE_GROUP2"):
        rows = [r for r in rows if str(r.get("category")) != "2"]
    popn = {k for k in SEL}
    for r in rows:
        r["in_primary_population"] = "yes" if r["pmid"] in popn else "no"
    cols = ["pmid", "label", "year", "era_2016", "era_2020", "set", "category", "tissue", "lane", "technology", "cases", "controls", "n",
            "list_status", "primary_provenance", "threshold_kind", "stratum", "in_primary_population"]
    covs = [c for c in rows[0].keys() if c.startswith("cov_") and c not in ("cov_src",)]
    # cov_ph_v2: pH handling from the full-text re-read of the 64 primary studies (2026-10-09); other studies keep the h2 cov_ph value
    P_PH64 = os.path.join(AN, "ph_agonal_handling_64_2026-10-09", "ph_agonal_handling_64.tsv")
    ph64 = {r["pmid"]: r["ph_call"] for r in csv.DictReader((l for l in open(P_PH64, encoding="utf8") if not l.startswith("#")), delimiter="\t")}
    for r in rows:
        if r["pmid"] in ph64:
            r["cov_ph_v2"], r["cov_ph_v2_source"] = ph64[r["pmid"]], "full-text re-read 2026-10-09 (64 primary studies)"
        else:
            r["cov_ph_v2"], r["cov_ph_v2_source"] = r.get("cov_ph", ""), ""
    # cov_<name>_v2: the 12 other covariates, re-called for the 64 primary studies from full text by two independent readers plus adjudication (2026-10-10)
    P_COV64 = os.path.join(AN, "covariate_handling_64_2026-10-10", "covariate_handling_64.tsv")
    cov64 = {(r["pmid"], r["covariate"]): r["final_call"] for r in csv.DictReader((l for l in open(P_COV64, encoding="utf8") if not l.startswith("#")), delimiter="\t")}
    COV12 = ["age", "sex", "pmi", "rin", "medication", "smoking", "substance_use", "cause_of_death_agonal", "batch", "cell_composition", "ancestry", "hidden_factors"]
    for r in rows:
        for c in COV12:
            r[f"cov_{c}_v2"] = cov64.get((r["pmid"], c), r.get(f"cov_{c}", ""))
        r["cov_v2_source"] = "full-text re-call 2026-10-10 (64 primary studies)" if (r["pmid"], "age") in cov64 else "h2 first pass"
    cols += covs + ["mt_final", "list_published", "deposited", "cov_ph_v2", "cov_ph_v2_source"] + [f"cov_{c}_v2" for c in COV12] + ["cov_v2_source"]
    F.write("fig01_b_reporting_completeness.tsv", "b", cols, rows, [P_H2, P_UNITS, P_COV64],
            ["one row per study of h2/h2_studies.tsv (the DE-reporting subset, not the whole corpus); cov_* = model status of the covariate in the main SZ vs "
             "control DE analysis: adjusted | sensitivity | matched_or_compared | reported_only | not_applicable | not_mentioned (blank = no text on disk); "
             "cov_n_applicable / cov_n_adjusted / cov_n_considered count over the stratum's applicable covariates",
             "cov_* values are the model reading with the hand corrections of h2/cov_corrections.tsv (see the h2_studies.tsv header)",
             f"cov_ph_v2 = pH handling call from analyses/ph_agonal_handling_64_2026-10-09/ph_agonal_handling_64.tsv (sha256 {sha(P_PH64)[:12]}) for the 64 primary studies; for other studies the h2 cov_ph value; cov_<name>_v2 = final call from analyses/covariate_handling_64_2026-10-10/covariate_handling_64.tsv (two independent readers, adjudicated) for the 64 primary studies, the h2 value for the others; the original cov_* columns are the h2 first pass"])
    F.note("b", f"reporting completeness comes from h2/h2_studies.tsv ({len(rows)} studies, sha {sha(P_H2)[:12]}); {sum(r['in_primary_population']=='yes' for r in rows)} "
                f"of the 64 population studies are in it. No separate report/figures.tsv table of reporting exists (report/figures.tsv lists figure files only).")
    # distinct subjects
    ds = rd(P_DS)
    F.write("fig01_b_distinct_subjects.tsv", "b", list(ds[0].keys()), ds, [P_DS],
            ["naive = sum of each study's own N; lower = maximal overlap within each bank or collection; upper = only measured overlap and fixed collection "
             "sizes collapsed; subjects, not samples. Rows: set (published brain = the 64; first-look = all tissues), level (total or cluster), measure (all / SZ / CTRL)",
             "same file as analyses/subject_relations/distinct_subjects.tsv, rows unchanged"])
    tot = {(r["set"], r["measure"]): r for r in ds if r["level"] == "total"}
    # report.md numbers
    exp = {("published brain", "all"): (5396, 1846, 4192), ("published brain", "SZ"): (2137, 732, 1685),
           ("published brain", "CTRL"): (2606, 857, 2046), ("first-look", "all"): (21078, 12649, 17135)}
    bad = []
    for k, (a, b, c) in exp.items():
        r = tot[k]
        if (int(float(r["naive"])), int(float(r["lower"])), int(float(r["upper"]))) != (a, b, c):
            bad.append(f"{k}: file {r['naive']}/{r['lower']}/{r['upper']} vs report {a}/{b}/{c}")
    F.note("b", "distinct subjects match report/report.md lines 118/128 (5,396; 1,846 to 4,192; SZ 2,137 / 732 to 1,685; controls 2,606 / 857 to 2,046; first-look 21,078, 12,649 to 17,135)."
           if not bad else "MISMATCH distinct subjects vs report/report.md: " + "; ".join(bad))
    F.note("b", f"distinct_subjects.tsv was built {RUN} 14:38 from lineage {LINEAGE_RUN_SHA}; the earlier export figure_data/distinct_subjects/ is from 2026-10-01 (lineage sha 6320781968316b4db78bfaf50b0c8b47146f40509fbb6537951004be92e59c42) and is superseded by this one.")


def fig1_c(F):
    inpop = {n["id"]: ("yes" if in_pop(n) else "no") for n in NODES}
    inc = {n["id"] for n in CORPUS}
    ids = {n["id"]: n for n in NODES}
    out = []
    for e in EDGES:
        out.append({"from": e["from"], "to": e["to"], "link_type": e["kind"],
                    "from_label": ids[e["from"]]["label"], "to_label": ids[e["to"]]["label"],
                    "curated": "yes" if e.get("curated") else "no", "match": e.get("match", ""),
                    "evidence": oneline(e.get("evidence", ""))})
    F.write("fig01_c_lineage_edges.tsv", "c", ["from", "to", "link_type", "from_label", "to_label", "curated", "match", "evidence"], out,
            [P_LINEAGE],
            ["from = the study whose data or subjects are reused, to = the study that reuses them (build_lineage.py writes them in this direction; "
             "'reviewed' goes from a reviewed study to the review, 'parent' from the source node to the child); link_type = lineage edge 'kind'",
             "evidence is the free text of the edge, whitespace collapsed",
             f"lineage.json sha256 now {sha(P_LINEAGE)[:12]}; the 2026-10-08 recurrence run used {LINEAGE_RUN_SHA}"])
    nrows = []
    for n in NODES:
        nrows.append({"id": n["id"], "pmid": n.get("pmid", ""), "label": n["label"], "year": n.get("year"), "relation": n.get("relation"),
                      "bank": n.get("lane"), "bank_label": n.get("lane_label"), "tissue": n.get("tissue"), "technology": n.get("tech_tier"),
                      "n_total": n.get("n"), "n_cases": n.get("cases"), "n_controls": n.get("controls"),
                      "times_cited": n.get("times_cited"), "in_corpus": "yes" if n["id"] in inc else "no",
                      "in_primary_population": inpop[n["id"]]})
    F.write("fig01_c_lineage_nodes.tsv", "c", list(nrows[0].keys()), nrows, [P_LINEAGE, P_UNITS],
            ["every node of lineage.json (in_corpus = has a tech_tier: the study nodes; the other nodes are cited sources outside the corpus); "
             "relation = lineage relation (primary, meta-analysis / combined, reanalysis of public data, source, single-cell / single-nucleus, ...)"])
    kinds = Counter(e["link_type"] for e in out)
    F.note("c", f"{len(out)} edges, {len(nrows)} nodes ({len(CORPUS)} in the corpus); link types: " + "; ".join(f"{k} {v}" for k, v in kinds.most_common()) + ".")
    F.note("c", f"edge direction as stored in lineage.json ('from' / 'to'); edge types read from the data (corpus/build_lineage.py writes 'parent', 'pooled study', 'same dataset', "
                "'reviewed', 'resource paper' and the curated kinds from corpus/curated_edges.json). "
                f"lineage.json is {sha(P_LINEAGE)[:12]} against {LINEAGE_RUN_SHA} used by the 2026-10-08 run: the graph here includes the later edit.")
    stale = Counter(n["relation"] for n in NODES if n["id"] not in inc)
    F.note("c", "the " + str(sum(stale.values())) + " non-corpus nodes are " + ", ".join(f"{k} {v}" for k, v in stale.items()) + ".")


def fig1_d(F):
    rows = rd(P_NEG)
    out = []
    for r in rows:
        nm = int(r["n_prior_meta"])
        nr = int(r["n_prior_reuse_only"])
        out.append({"study_id": r["pmid"], "label": r["label"], "bank": r["lane"], "dataset_identity": r["dataset_identity"],
                    "ever_pooled": "yes" if nm > 0 else "no", "ever_reanalysed": "yes" if nr > 0 else "no",
                    "n_prior_meta": nm, "n_prior_reanalysis_only": nr,
                    "ever_reused_any": "yes" if (nm + nr) > 0 else "no",
                    "first_prior_meta": r["first_prior_meta"], "match_levels": r["match_levels"],
                    "in_primary_population": "yes" if r["pmid"] in SELSET else "no"})
    F.write("fig01_d_neglected_datasets.tsv", "d", list(out[0].keys()), out, [P_NEG],
            ["ever_pooled = n_prior_meta > 0 (pooled by a prior meta-analysis or combined analysis); ever_reanalysed = n_prior_reuse_only > 0 "
             "(a prior single-dataset reanalysis or data-reuse paper that is not a meta-analysis); ever_reused_any = either",
             "counts are floors: the re-use layer is read by a model from the full texts available (analyses/neglected_data_2026-10-03/README.md, caveats)"])
    n = len(out)
    pooled = sum(r["ever_pooled"] == "yes" for r in out)
    neg = [r for r in out if r["ever_pooled"] == "no"]
    neg_no_reuse = sum(r["ever_reused_any"] == "no" for r in neg)
    once = sum(r["n_prior_meta"] == 1 for r in out)
    pop_in = sum(r["in_primary_population"] == "yes" for r in out)
    nsel = len(SELSET)
    F.note("d", f"{n} units; ever pooled {pooled}, never pooled {len(neg)}; of the never pooled, {neg_no_reuse} never re-used at all; pooled exactly once {once}; "
                f"{pop_in} of {n} study ids are in the current {nsel}-study selection.")
    F.note("d", f"earlier states of this table: 28 never pooled of 64 (2026-10-03 README: 36 pooled, 24 never re-used, 7 pooled once); 26 never pooled of 62 "
                f"(2026-10-12 delivery before this one: 36 pooled, 22 never re-used, 7 once). Compare: pooled {pooled}, never pooled {len(neg)}, "
                f"never re-used {neg_no_reuse}, pooled once {once} of {n}.")
    F.note("d", f"results.tsv regenerated {RUN} by analyses/neglected_data_2026-10-03/prior_meta_use.py from the current recurrence_units.tsv "
                f"(units = scenario '{SCENARIO}', modelled = yes); its units equal the current selection."
           if set(r["study_id"] for r in out) == SELSET else
           "the results.tsv units differ from the current selection: only in file: "
           + ", ".join(sorted({r['study_id'] for r in out} - SELSET)) + "; only in selection: " + ", ".join(sorted(SELSET - {r['study_id'] for r in out})))
    # population funnel note
    units = rd(P_UNITS)
    pb = [u for u in units if u["scenario"] == SCENARIO]
    pj = json.load(open(P_PRIM))
    F.note("d", f"'Population funnel' note (notes/, 2026-10-03) quotes 94 collapsed published-brain units, 92 donor groups and 64 in the population (and 327 studies with a primary list, not rechecked: "
                f"no single table holds it). recurrence_units.tsv now has {len(pb)} rows for scenario '{SCENARIO}', of which modelled=yes {sum(u['modelled']=='yes' for u in pb)}; "
                f"the note's 94 is not reproduced by that row count ({len(pb)}), so the 94 / 92 figures are not exported (no file holds them other than the note).")


def fig11_a(F):
    meta = json.load(open(os.path.join(SMK, "overlap_meta.json")))
    sigs = meta["signatures"]
    srows = [{"signature": k, "n_genes": v["raw_n"], "kind": "protein" if v["proteome"] else "gene (expression or methylation)",
              "arm": "tested against 31 brain proteomics units (proteomics.tsv)" if v["proteome"] else "tested against recurrent-gene sets by era (overlap_era.tsv)",
              "note": v["note"]} for k, v in sigs.items()]
    F.write("fig11_a_signatures.tsv", "a", ["signature", "n_genes", "kind", "arm", "note"], srows,
            [os.path.join(SMK, "overlap_meta.json")], ["all six signatures; n_genes = signature size from overlap_meta.json 'raw_n'"])
    era_rows = rd(os.path.join(SMK, "overlap_era.tsv"))
    out = []
    for r in era_rows:
        s = r["signature"]
        out.append({"signature": s, "n_genes": sigs[s]["raw_n"], "analysis": r["analysis"], "scenario": r["scenario"], "split": r["split"],
                    "era": r["era"], "k": r["k"], "n_studies": r["n_studies"], "n_recurrent": r["n_recurrent"],
                    "recurrent_in_signature": r["recurrent_in_signature"], "pct_recurrent": r["frac_recurrent"],
                    "n_background": r["n_background"], "background_in_signature": r["background_in_signature"],
                    "pct_background": r["frac_background"], "OR": r["fisher_or"], "p": r["fisher_p"], "note": r["note"]})
    F.write("fig11_a_smoking_recurrence.tsv", "a",
            ["signature", "n_genes", "analysis", "scenario", "split", "era", "k", "n_studies", "n_recurrent", "recurrent_in_signature", "pct_recurrent",
             "n_background", "background_in_signature", "pct_background", "OR", "p", "note"], out, [os.path.join(SMK, "overlap_era.tsv"), os.path.join(SMK, "overlap_meta.json")],
            ["pct_recurrent / pct_background are the file's frac_recurrent / frac_background as fractions (0.069 = 6.9%), rounded to three places in the source; "
             "k = recurrence threshold (2, 3, or n160 = the K matching the Merikangas 160-gene list); analysis 'old vs new recurrent sets' rows hold 'a vs b' values in the era-pair columns",
             "recurrent sets here are era-based protein-coding recurrence sets (K studies of one representative per same_data_as group), NOT the 176-gene study-tier set",
             "Corlin 2020 and Xiao 2025 are protein signatures: they have no recurrent-gene rows, only fig11_a_smoking_proteomics.tsv"])
    pairs = rd(os.path.join(SMK, "overlap_pairs.tsv"))
    F.write("fig11_a_smoking_cross_signature.tsv", "a", list(pairs[0].keys()), pairs, [os.path.join(SMK, "overlap_pairs.tsv")],
            ["pairwise overlap of the smoking signatures against the HGNC protein-coding background"])
    prot = rd(os.path.join(SMK, "proteomics.tsv"))
    F.write("fig11_a_smoking_proteomics.tsv", "a", list(prot[0].keys()), prot, [os.path.join(SMK, "proteomics.tsv")],
            ["protein signatures (Xiao 2025, Corlin 2020) against the 31 brain proteomics units; universe = union of the units' tested gene sets"])
    # status
    F.note("a", "smoking run: overlap_era.tsv / overlap_pairs.tsv / proteomics.tsv / overlap_meta.json are dated 2026-10-08 14:37 and record recurrence_units.tsv sha c71d007ba4b6 "
                f"(current {sha(P_UNITS)[:12]}), primary_lists.json e282904dcacb (current {sha(P_PRIM)[:12]}) and lineage.json {LINEAGE_RUN_SHA}. The signature stage itself ran 2026-10-05.")
    F.note("a", "the smoking analysis does NOT use the 176-gene study-tier set. It tests era-based recurrent sets (K = 2, 3 and the Merikangas-160-matched K) in three scenarios "
                "('published only, all tissues', 'published, brain', 'published, blood'), one representative per same_data_as group, protein-coding lists. "
                "A smoking-versus-176 table would need a new run: not produced.")
    F.note("a", "all six signatures: four (Charlesworth 2010, Huan 2016, Joehanes 2016, Vink 2017) have recurrent-gene tests; Corlin 2020 and Xiao 2025 are plasma-protein signatures "
                "and exist only in the per-unit proteomics table (no pct_recurrent / pct_background / OR for them in the recurrence sense).")
    # README checks
    readme = open(os.path.join(SMK, "README.md"), encoding="utf-8").read()
    idx = {(r["signature"], r["scenario"], r["split"], r["era"], r["k"]): r for r in era_rows}
    bad, nchk = [], 0
    for m in re.finditer(r"^- (?P<sig>[A-Za-z]+ \d{4}): recurrent-in-signature (?P<a>\d+) of (?P<n>\d+) \((?P<f>[\d.]+)% vs (?P<b>[\d.]+)%\), OR (?P<or>[\d.]+), p=(?P<p>[\d.e-]+)", readme, re.M):
        hit = [r for r in era_rows if r["signature"] == m["sig"] and r["scenario"] == "published, brain" and r["recurrent_in_signature"] == m["a"]
               and r["n_recurrent"] == m["n"] and float(r["fisher_or"]) == float(m["or"])]
        nchk += 1
        if not hit or float(hit[0]["fisher_p"]) != float(m["p"]):
            bad.append(f"brain bullet {m['sig']} {m['a']}/{m['n']}: no matching row")
    for m in re.finditer(r"^\| (?P<sig>[A-Za-z]+ \d{4}) \| (?P<a>\d+) \((?P<f>[\d.]+)%\) \| (?P<b>[\d.]+)% \| OR (?P<or>[\d.]+), p=(?P<p>[\d.e-]+) \|", readme, re.M):
        hit = [r for r in era_rows if r["signature"] == m["sig"] and r["scenario"] == "published, blood" and r["recurrent_in_signature"] == m["a"]
               and float(r["fisher_or"]) == float(m["or"])]
        nchk += 1
        if not hit or float(hit[0]["fisher_p"]) != float(m["p"]):
            bad.append(f"blood table {m['sig']} {m['a']}: no matching row")
    pm = {r["pair"]: r for r in pairs}
    for m in re.finditer(r"^\| (?P<pair>[A-Za-z]+ \d{4} x [A-Za-z]+ \d{4}) \| (?P<o>\d+) \| (?P<e>[\d.]+) \| (?P<f>[\d.]+) \| p=(?P<p>[\d.e-]+) \|", readme, re.M):
        r = pm.get(m["pair"])
        nchk += 1
        if not r or int(r["overlap"]) != int(m["o"]) or float(r["expected"]) != float(m["e"]) or float(r["fold"]) != float(m["f"]) or float(r["fisher_p"]) != float(m["p"]):
            bad.append(f"pair {m['pair']}")
    F.note("a", f"README.md numbers checked against the TSVs: {nchk} rows checked, " + ("all match." if not bad and nchk else ("MISMATCH: " + "; ".join(bad) if bad else "none parsed.")))
    old = [float(r["fisher_p"]) for r in era_rows if r["scenario"] == "published, blood" and r["split"] == "2016" and r["era"].startswith("old") and r["k"] == "3"
           and r["analysis"] == "era recurrence"]
    F.note("a", "README section A heading says every signature is enriched in blood (new era 2016, K = 3); the old-era blood rows at K = 3 (a handful of recurrent genes) have Fisher p "
                + ", ".join(f"{x:g}" for x in old) + " in the TSV.")


def fig11_periphery(F):
    P = os.path.join(PER, "set_overlap.tsv")
    so = rd(P)
    F.write("fig11_periphery_set_overlap.tsv", "periphery", list(so[0].keys()), so, [P],
            ["all scenarios, sets and nulls as in the source: obs_hits = observed (gene, cluster) hits, exp_hits = expected under the null, obs_exp = ratio; "
             "ci95 = 95% interval of obs/exp over 2,000 cluster resamples; null: uniform | de_prior | expr_libd | de_prior_blood",
             "set ids keep their historical names: study_133 is the 176-gene study-tier set, lists5plus_1151 the 1,478-gene 5+ set (README of the periphery folder)"])
    P2 = os.path.join(PER, "direction.tsv")
    d = rd(P2)
    F.write("fig11_periphery_direction.tsv", "periphery", list(d[0].keys()), d, [P2],
            ["hits whose cluster direction and brain direction are both up or down; concordant / discordant out of hits_with_direction; binom_p two-sided vs 0.5"])
    P3 = os.path.join(PER, "ph_contrast.tsv")
    c = rd(P3)
    F.write("fig11_periphery_ph_contrast.tsv", "periphery", list(c[0].keys()), c, [P3],
            ["recurrent genes with control pH |Z| < 1 (lt1) against >= 1 (ge1): ratio_lt1_over_ge1 = obs_exp_lt1 / obs_exp_ge1; ci95 by cluster resampling; "
             "p_conditional_binomial as defined in the source header"])
    # status
    hdr = open(P).read().splitlines()[:3]
    shas = dict(re.findall(r"analyses/(\w+\.json) ([0-9a-f]{12})", hdr[1]))
    shas["lineage.json"] = re.search(r"corpus/lineage.json ([0-9a-f]{12})", hdr[1]).group(1)
    cur = {"primary_lists.json": sha(P_PRIM)[:12], "universes.json": sha(P_UNIV)[:12], "lineage.json": sha(P_LINEAGE)[:12]}
    diff = [f"{k} run {shas.get(k)} current {v}" for k, v in cur.items() if shas.get(k) != v]
    g = os.path.join(AN, "recurrence_genes.tsv")
    F.note("periphery", "STALE: set_overlap / direction / ph_contrast were written 2026-10-03 (population of that date). Inputs recorded in their headers differ from today's files: "
                        + "; ".join(diff) + f"; recurrence_genes.tsv run 5a0991977a32 (build_sets.py) current {sha(g)[:12]}. The tables are exported as they are; a rerun (scripts/overlap.py) is needed for the 2026-10-08 population.")
    sets = rd(os.path.join(PER, "..", "sets.tsv"))
    s133 = {r["gene"] for r in sets if r["set"] == "study_133"}
    cur176 = {r["gene"] for r in rd(g) if r["scenario"] == SCENARIO and float(r["q_bh"]) < 0.1}
    F.note("periphery", f"the 176-gene set in sets.tsv (study_133) has {len(s133)} genes; recurrence_genes.tsv 'published, brain' q_bh < 0.1 now gives {len(cur176)}; "
                        f"in both {len(s133 & cur176)}, only in sets.tsv {len(s133 - cur176)}, only now {len(cur176 - s133)}.")
    # README checks: obs/exp ratios in the main table and the direction table
    readme = open(os.path.join(PER, "..", "README.md"), encoding="utf-8").read()
    so_i = {(r["scenario"], r["set"], r["null"]): r for r in so}
    di = {(r["scenario"], r["set"]): r for r in d}
    bad, n = [], 0
    nullcols = {"blood": ("uniform", "de_prior", "de_prior_blood"), "ipsc_neural": ("uniform", "de_prior", "expr_libd"), "other_cells": ("uniform", "de_prior")}
    setmap = {"176 up": "study_133_up", "176 down": "study_133_down", "careful 25": "careful_25", "interneuron 14": "interneuron_14", "IEG 4": "ieg_4", "5+ lists 1,478": "lists5plus_1151"}
    for m in re.finditer(r"^\| (?P<set>176 up|176 down|careful 25|interneuron 14|IEG 4|5\+ lists 1,478) \| (?P<b>[^|]+)\| (?P<i>[^|]+)\| (?P<o>[^|]+)\|$", readme, re.M):
        for scen, cell in (("blood", m["b"]), ("ipsc_neural", m["i"]), ("other_cells", m["o"])):
            obs = re.match(r"\s*([\d,]+)", cell)
            ratios = re.findall(r"([\d.]+)x", cell)
            if not obs:
                continue
            key = setmap[m["set"]]
            rowsx = [so_i.get((scen, key, nl)) for nl in nullcols[scen]]
            n += 1
            if rowsx[0] is None or int(rowsx[0]["obs_hits"]) != int(obs.group(1).replace(",", "")):
                bad.append(f"{m['set']} {scen} obs {obs.group(1)} vs {rowsx[0] and rowsx[0]['obs_hits']}")
                continue
            for nl, rx in zip(nullcols[scen], ratios):
                r = so_i.get((scen, key, nl))
                n += 1
                if r is None or abs(float(r["obs_exp"]) - float(rx)) > 0.0051:
                    bad.append(f"{m['set']} {scen} {nl}: README {rx} vs file {r and r['obs_exp']}")
    for m in re.finditer(r"^\| (?P<set>176 up|176 down|careful 25|interneuron 14|IEG 4|5\+ lists) \| (?P<br>[\d,]+)/(?P<brn>[\d,]+) \| (?P<b>[^|]+)\| (?P<i>[^|]+)\| (?P<o>[^|]+)\|$", readme, re.M):
        key = setmap.get(m["set"], "lists5plus_1151" if m["set"] == "5+ lists" else None)
        for scen, cell in (("brain", f"{m['br']}/{m['brn']}"), ("blood", m["b"]), ("ipsc_neural", m["i"]), ("other_cells", m["o"])):
            mm = re.match(r"\s*([\d,]+)/([\d,]+)", cell)
            if not mm:
                continue
            r = di.get((scen, key))
            n += 1
            if not r or int(r["concordant"]) != int(mm.group(1).replace(",", "")) or int(r["hits_with_direction"]) != int(mm.group(2).replace(",", "")):
                bad.append(f"direction {m['set']} {scen}: README {mm.group(1)}/{mm.group(2)} vs file {r and (r['concordant'], r['hits_with_direction'])}")
    F.note("periphery", f"periphery README tables (obs/expected ratios, direction counts) checked against the TSVs: {n} cells, "
                        + ("all match." if not bad else f"{len(bad)} mismatch(es): " + "; ".join(bad[:12])))
    F.note("periphery", "ph_top200 rows and the 'In the low-pH direction' table of the README were not parsed; they are in direction.tsv (set ph_top200) as exported.")


def fig11_pending(F):
    F.note("b (iPSC / LCL recurrence within scenario)", "NOT PRODUCED: pending. FRZ is running that analysis.")
    F.note("c (proteomics)", "NOT PRODUCED: no analysis exists. The smoking folder's proteomics.tsv (smoking plasma proteins in 31 brain proteomics units) is exported with panel a "
                             "and is a different question.")


def main():
    F1 = Figure("fig01")
    for fn in (fig1_b, fig1_c):   # core-80: panel b corpus table and panel c lineage; d, fig11 not run here
        fn(F1)
    F1.finish()


if __name__ == "__main__":
    main()
