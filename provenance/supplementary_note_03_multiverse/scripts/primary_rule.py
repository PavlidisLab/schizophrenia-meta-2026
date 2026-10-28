#!/usr/bin/env python3
"""Primary-list rule: candidate written rules, and recurrence and H1 under each (sensitivity analysis).

Nothing here changes the project defaults: primary_lists.json, universes.json, recurrence*.tsv and h1_* in
analyses/ are read, never written. Every output goes to analyses/primary_list_rule_sensitivity/.

    python3 primary_rule.py select [--inputs DIR]                 candidates.tsv, choice_dimensions.tsv,
                                                                   rule_choices.tsv, <rule>/primary_lists.json,
                                                                   <rule>/universes.json, sets/*.txt
    python3 primary_rule.py run [--inputs DIR] [--code DIR]       <rule>/recurrence*.tsv, <rule>/h1_*.tsv
                                                                   (current, R1, R2 through the unmodified
                                                                   recurrence.main and h1.main), R3/*weighted*.tsv
    python3 primary_rule.py check [--inputs DIR] [--code DIR]     R3_code_check.tsv: R3's weighted code with one list per
                                                                   study against current/ (must match exactly)
    python3 primary_rule.py report                                 primary_rule_report.md, from the tsv files

--inputs: a directory holding primary_lists.json, universes.json, technology.json, manifest.json, manifests/,
lists/, universe_sets/, reference/ and corpus/lineage.json (default: the live analyses/, its data items in ../data/, and
../corpus/).
--code: where recurrence.py, h1.py, build_universes.py and extract.py are imported from (default: analyses/).
The run of 2026-09-24/25 used a snapshot of both, taken at 23:34 on 2026-09-24 while another agent was rerunning
the main analyses; inputs_sha256.tsv lists every input file with its sha256.

How the existing scripts are reused without editing them: for each rule a temporary run directory is built with
the rule's primary_lists.json and universes.json, symlinks to the (snapshot) lists, manifests, universe sets and
reference files, and symlinks to the live raw/, gemma/, transcribed/, reproduce/, manual/ and platforms/ (read
only). The module-level HERE and DATA of recurrence, h1, build_universes and extract (and the paths derived from them) are
pointed at that directory, recurrence.main() and h1.main() are called, and their output files are copied here.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import importlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
from functools import reduce

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.path.dirname(HERE)
DATA = os.path.join(os.path.dirname(GL), "data")
SETS = os.path.join(HERE, "sets")
RULES = ("current", "R1", "R2", "R3", "A")
# rule A (Paul, 2026-09-25): list kinds that count as the authors' SZ finding besides DE
A_MODALITIES = ("DE", "gene selection", "other statistic")
RANK = {"published": 0, "reproduced": 1, "gemma": 2}
LIVE_DATA = ("raw", "gemma", "transcribed", "reproduce", "manual", "platforms")
TODAY = datetime.date.today().isoformat()  # was hard-coded "2026-09-25" (code audit CA26, 2026-09-26)

RULE_TEXT = {
    "R1": ("One list per study, chosen by ordered criteria. Candidates are the lists already attached to the study "
           "in primary_lists.json (primary and alternates) whose manifest role is first_look and modality DE; "
           "up and down lists of one contrast are one candidate. Keep the candidates that are best on criterion 1, "
           "then among those the best on criterion 2, and so on until one is left: "
           "(1) gene-level list over miRNA, lncRNA or circRNA lists; "
           "(2) list in the study's recorded tissue class (primary_lists.json `tissue`); "
           "(3) SZ or SZ-spectrum cases over case groups that pool other diagnoses (bipolar disorder, other psychoses); "
           "(4) all cases and controls over a subset or condition (sex, age, illness duration, treatment status, "
           "time of death, clinical subtype, within-twin contrast, stimulated cells, later time point); "
           "(5) provenance published over reproduced over Gemma; "
           "(6) FDR-controlled threshold (FDR, q, Bonferroni, adjusted p) over nominal, fold-only or other; "
           "(7) prefrontal region over any other; "
           "(8) a model with added covariates over its unadjusted sibling; "
           "(9) most subjects (when every remaining candidate has a count in its manifest contrast); "
           "(10) first in extraction order (manifest.json, then manifests/*.json)."),
    "R2": ("Union within study. Apply R1 criteria (1) to (6) to the same candidates; take the union of every "
           "candidate left (regions, cell types, layers, stages, cohorts, statistical methods, and subsets when no "
           "all-subject list exists). The union list's tested universe is the union of the members' universes."),
    "A": ("Author-generous union (Paul, 2026-09-25: would the authors want it included; could they say we were too "
          "stringent?). Candidates: the study's first_look DE and DE-adjacent lists (tracker Vocabulary: lists the authors "
          "present as their SZ finding without a plain SZ vs control DE test; manifest modality DE, gene selection or other "
          "statistic). Apply R1 criteria (1) gene-level, (2) the study's "
          "tissue class and (5) the authors' own lists over reproductions and Gemma; take the union of every candidate left, "
          "whatever the case group (pooled psychosis, relatives as controls), subset or threshold. Where one contrast is "
          "reported at several thresholds the union is the least stringent. Universe: union of the members' universes."),
    "R3": ("Per list, study weighted 1/k. Apply R1 criteria (1) to (6); keep each of the k candidates left as its "
           "own list with weight 1/k, so each study contributes a total weight of 1. Study tier: a gene's recurrence "
           "is the sum of the weights of the lists that contain it, tested against the weighted Poisson-binomial "
           "null. Cohort tier: every list is a member of its cluster (a cluster counts once whatever the weights). "
           "H1: every list pair of two studies is scored, weighted 1/(k_A k_B)."),
}


def sha(p, n=64):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:n]


def data_dir(inp):
    """Directory of the data items of an inputs dir (manifest.json, manifests/, lists/, reference/, universe_sets/): the
    repo's data/ for the live analyses/ dir, the inputs dir itself for a snapshot taken in the one-folder layout."""
    return DATA if inp == GL else inp


def inputs_paths(inp):
    lineage = os.path.join(inp, "corpus", "lineage.json") if inp != GL else os.path.join(GL, "..", "corpus", "lineage.json")
    return {"pl": os.path.join(inp, "primary_lists.json"), "uni": os.path.join(inp, "universes.json"),
            "tech": os.path.join(inp, "technology.json"), "lineage": lineage,
            "hgnc": os.path.join(data_dir(inp), "reference", "hgnc_complete_set.txt"),
            # the DE prior in force (recurrence.DE_PRIOR_SOURCE): ours (live files of
            # analyses/), Crow 2019 (data/reference/DE_Prior.txt) with DE_PRIOR_SOURCE=crow2019
            "prior": {"crow2019": os.path.join(data_dir(inp), "reference", "DE_Prior.txt"),
                      "global": os.path.join(GL, "de_prior_lists", "de_prior_ranks", "1_global_de_prior.tsv")}.get(
                          os.environ.get("DE_PRIOR_SOURCE", "ours"), os.path.join(GL, "de_prior_lists", "our_prior", "1_global_our_prior.tsv"))}


def load_code(code):
    sys.path.insert(0, code)
    mods = {}
    for m in ("extract", "build_universes", "recurrence", "h1"):
        mods[m] = importlib.import_module(m)
    return mods


def manifest_specs(inp):
    m = json.load(open(os.path.join(data_dir(inp), "manifest.json")))["lists"]
    for f in sorted(glob.glob(os.path.join(data_dir(inp), "manifests", "*.json"))):
        m += json.load(open(f))["lists"]
    return {f"{s['pmid']}__{s['list_id']}": s for s in m}


def read_genes(inp, name):
    rows = [l for l in open(os.path.join(data_dir(inp), "lists", name + ".tsv"), encoding="utf-8") if not l.startswith("#")]
    return {r["symbol"] for r in csv.DictReader(rows, delimiter="\t") if r.get("symbol")}


# ---------------------------------------------------------------- per-list attributes (from manifest text)
MOLECULE = [("miRNA", r"mirna|\bmir\b|microrna"), ("lncRNA", r"lncrna"), ("circRNA", r"circrna")]
TISSUE = [("ipsc", r"ipsc|organoid|spheroid|hipsc"), ("lcl", r"lymphoblast|\blcl\b"),
          ("blood", r"blood|pbmc|lymphocyte|leukocyte|monocyte|plasma|serum|exosome"),
          ("other", r"fibroblast|olfactory|nasal"),
          ("brain", r"cortex|dlpfc|\bba ?\d|hippocamp|striat|caudate|putamen|accumbens|amygdala|cerebell|thalam|"
                    r"midbrain|choroid|habenula|subependymal|postmortem|brain|prefrontal|cingulate|temporal|neun|olig2")]
# case-group text (the contrast before its first ' vs ', or all of it when there is none) naming other diagnoses
POOLED = r"major psychosis|bipolar|\bBD\b|\bBP\b|other psycho|first-episode psychosis|psychotic disorder|^\W*psychosis\b"


def case_part(contrast):
    m = re.search(r"\bvs\.?\b", contrast)
    return contrast[:m.start()] if m else contrast
# list-id tokens that mark a subset of the subjects or a non-baseline condition
SUBSET = [("sex", r"(^|_)(male|female|men|women)(_|$)"), ("age", r"(^|_)age\d|aged"),
          ("illness duration", r"_doi\b|_doi$|duration"), ("time of death", r"(^|_)(night|day)(_|$)"),
          ("treatment status", r"drug_treated|drug_naive|treatment|(^|_)(on|off)_ap(_|$)|apnaive|apswitching|ssd_v3"),
          ("stimulated cells", r"dopamine"), ("clinical subtype", r"mrds|type2|hiscz"),
          ("within-twin contrast", r"cotwin|st_vs_ht")]
ADJUSTED = r"added as covariates|layer-adjusted|adjusted for|corrected for"
NOMINAL = r"nominal|uncorrected|unadjusted"
FDR = r"FDR|false discovery|q[- ]?value|\bq ?(<|<=|≤)|adj(usted)?[ .]?p|p[ .-]?adj|padj|Bonferroni|Benjamini|\bBH\b|Storey|corrected for multiple"
PVAL = r"\bp ?(<|<=|≤)|p value|t-test|\bP\b"
# other ways two candidates of one study differ (list id or contrast), for the characterization only
POPULATION = [("cell type, layer, nucleus or stage", r"neuron|npc|gpc|astrocyte|ngn2|(^|_)ad2|glut|gaba|migrated|sphere|"
               r"(^|_)d\d+(_|$)|endothelial|neuronal|neun|olig2|superficial|deep|infragranular|supragranular|(^|_)an_|"
               r"vpl_|carseq|layer|(^|_)dg$|ca1|ca3|microglia"),
              ("cohort or dataset", r"newcollection|cchpc|hbb_|gse\d+(?!.*gemma)|_ac_|sample|(^|_)gse"),
              ("statistical method or model", r"deseq2|edger|cash|ebayes|genemeta|fisherz|ct_adjusted|restricted|adjusted|interaction|psychencode")]


def first_match(table, text):
    for lab, rx in table:
        if re.search(rx, text, re.I):
            return lab
    return ""


def n_subjects(contrast):
    m = re.search(r"\((\d+)\s*vs\.?\s*(\d+)\)", contrast)
    if m:
        return int(m.group(1)) + int(m.group(2))
    m1 = re.search(r"\((\d+)", contrast)
    if not m1:
        return None
    rest = contrast[m1.end():]
    v = re.search(r"\bvs\.?\b", rest)
    if not v:
        return None
    m2 = re.search(r"\((\d+)", rest[v.end():])
    if not m2:
        # 'control ('reference subject role') (11)'
        m2 = re.search(r"\)\s*\((\d+)", rest[v.end():])
    return int(m1.group(1)) + int(m2.group(1)) if m2 else None


def core(t):
    """threshold text before its first parenthesis"""
    return t.split(" (")[0].split("(")[0].strip()


def threshold_kind(spec, primary_spec, study_kind, is_primary):
    """the study's recorded kind for its current primary and for lists whose threshold text (before any parenthesis) is
    the primary's; else from the threshold text before any parenthesis."""
    if is_primary:
        return study_kind
    t = spec.get("threshold", "")
    if primary_spec is not None and t and core(t) == core(primary_spec.get("threshold", "")):
        return study_kind
    t = core(t) or t
    if spec.get("provenance") == "gemma":
        return "fdr"
    if re.search(NOMINAL, t, re.I):
        return "nominal"
    if re.search(FDR, t):
        return "fdr"
    if re.search(PVAL, t, re.I):
        return "nominal"
    if re.search(r"fold", t, re.I):
        return "fold_only"
    return "other"


def unit_key(name):
    """up and down lists of one contrast are one candidate."""
    lid = name.split("__", 1)[1]
    return re.sub(r"(^|_)(up|down)(_|$)", r"\1\3", lid).strip("_")


def candidates(inp, H):
    """study -> list of candidate units (dict), and all attributes, for in-corpus studies with a primary list."""
    plfull = json.load(open(inputs_paths(inp)["pl"]))
    pl = plfull["studies"]
    specs = manifest_specs(inp)
    mindex = {n: i for i, n in enumerate(specs)}  # position in manifest.json then manifests/*.json (extraction order)
    out = {}
    for k, v in pl.items():
        if not (v.get("in_corpus", True) and v.get("primary")):
            continue
        names = list(v["primary"]) + [n for n in v.get("alternates", []) if n not in v["primary"]]
        names = [n for n in names if n not in (v.get("excluded_lists") or [])]
        prim_spec = specs.get(v["primary"][0])
        units, seen = [], {}
        for order, n in enumerate(names):
            s = specs.get(n, {})
            role, mod = s.get("role", "first_look"), s.get("modality", "DE")
            is_p = n in v["primary"]
            key = "__primary__" if is_p else unit_key(n)
            if key in seen:
                seen[key]["lists"].append(n)
                continue
            text = n.split("__", 1)[1] + " " + s.get("contrast", "")
            sub = first_match(SUBSET, n.split("__", 1)[1])
            u = {"study": k, "lists": [n], "order": order, "is_current": is_p, "role": role, "modality": mod,
                 "provenance": s.get("provenance", "published"),
                 "molecule": first_match(MOLECULE, text) or "gene",
                 "tissue": first_match(TISSUE, s.get("contrast", "") + " " + n.split("__", 1)[1]) or v["tissue"],
                 "pooled": bool(re.search(POOLED, case_part(s.get("contrast", "")))),
                 "subset": sub, "adjusted": bool(re.search(ADJUSTED, s.get("contrast", ""), re.I)),
                 "threshold_kind": threshold_kind(s, prim_spec, v["threshold_kind"], is_p),
                 "threshold": s.get("threshold", ""), "region": H.region(specs, n),
                 "population": first_match(POPULATION, re.sub(r"^gemma_gse\d+(_\d+)?_?", "", n.split("__", 1)[1], flags=re.I)),
                 "n_subjects": n_subjects(s.get("contrast", "")), "contrast": s.get("contrast", "")}
            # the study's recorded tissue applies to its current primary and to lists whose contrast names no tissue
            if is_p or not first_match(TISSUE, s.get("contrast", "")):
                u["tissue"] = v["tissue"]
            u["manifest_index"] = mindex.get(n, 10**6)
            seen[key] = u
            units.append(u)
        for u in units:
            g = set()
            for n in u["lists"]:
                g |= read_genes(inp, n)
            u["n_genes"] = len(g)
        out[k] = {"label": v["label"], "tissue": v["tissue"], "threshold_kind": v["threshold_kind"],
                  "provenance": v["provenance"], "reason": v.get("reason", ""), "confidence": v.get("confidence", ""),
                  "units": units, "eligible": [u for u in units if u["role"] == "first_look" and u["modality"] == "DE"],
                  "eligible_A": [u for u in units if u["role"] == "first_look" and u["modality"] in A_MODALITIES]}
    return out, plfull, specs


# ---------------------------------------------------------------- rules
def criteria(study):
    t = study["tissue"]
    return [
        ("1 gene-level", lambda u: u["molecule"] == "gene"),
        ("2 study tissue", lambda u: u["tissue"] == t),
        ("3 SZ-only cases", lambda u: not u["pooled"]),
        ("4 all subjects", lambda u: not u["subset"]),
        ("5 provenance", lambda u: -RANK.get(u["provenance"], 3)),
        ("6 FDR", lambda u: u["threshold_kind"] == "fdr"),
        ("7 prefrontal", lambda u: u["region"] == "prefrontal"),
        ("8 adjusted model", lambda u: u["adjusted"]),
        ("9 subjects", "n_subjects"),
        ("10 manifest order", lambda u: -u["manifest_index"]),
    ]


def criteria_genes(study):
    """R1 with most genes before manifest order (tie-break sensitivity only)."""
    c = criteria(study)
    return c[:-1] + [("10 genes", lambda u: u["n_genes"]), ("11 manifest order", lambda u: -u["manifest_index"])]


def criteria_A(study):
    """Rule A keeps only R1 criteria (1) gene-level, (2) the study's tissue class and (5) provenance (the authors' own
    lists over reproductions and Gemma); pooled diagnoses, subsets and nominal thresholds all stay."""
    c = criteria(study)
    return [c[0], c[1], c[4]]


def apply(study, upto=None, crit=None, pool="eligible"):
    """Filter the eligible units by the criteria (all, or the first `upto`); returns (units left, deciding criterion)."""
    left = list(study[pool])
    decided = ""
    for i, (name, f) in enumerate((crit or criteria)(study)):
        if upto is not None and i >= upto:
            break
        if len(left) <= 1:
            break
        if f == "n_subjects":
            if any(u["n_subjects"] is None for u in left):
                continue
            f = lambda u: u["n_subjects"]  # noqa: E731
        best = max(f(u) for u in left)
        new = [u for u in left if f(u) == best]
        if len(new) < len(left):
            decided = name
        left = new
    return left, decided


def rule_pick(study, rule):
    if rule == "A":
        if not study["eligible_A"]:
            return [], "no eligible list"
        return apply(study, crit=criteria_A, pool="eligible_A")
    if not study["eligible"]:
        return [], "no eligible list"
    if rule in ("R1", "R1g"):
        left, d = apply(study, crit=criteria if rule == "R1" else criteria_genes)
        return left[:1], d
    left, d = apply(study, upto=6)
    return left, d


# ---------------------------------------------------------------- universes for lists the current file has none for
def universe_record(name, study_key, pl_s, uni, specs, BU, extract, inp, cache):
    """Universe record for one list: the current record if universes.json has it; else a Gemma list's full table; else
    the manifest platform; else the list's own table with threshold filters removed when it covers at least 3x the
    list and 5,000 genes; else the current primary list's record (inherited)."""
    if name in uni:
        return dict(uni[name], universe_from="universes.json")
    if name in cache:
        return cache[name]
    base = uni.get(pl_s["primary"][0])
    s = specs[name]
    L = read_genes(inp, name)
    rec = {"study": study_key, "label": pl_s["label"], "provenance": s.get("provenance", "published"),
           "tissue": pl_s["tissue"], "n_tested": None, "n_tested_unit": None, "n_tested_source": None,
           "list_genes": len(L)}

    def with_set(typ, genes, fname, src, how):
        path = os.path.join(SETS, fname + ".txt")
        os.makedirs(SETS, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(f"# {how}\n# source: {src}\n# script: primary_rule.py select, run {TODAY}\n")
            fh.write("\n".join(sorted(genes)) + "\n")
        rec.update(type=typ, source=src, gene_set="rule_sets/" + fname + ".txt", n_membership=len(genes),
                   list_genes_outside=len(L - genes), n_k=len(genes), n_k_from="n_membership", evidence=how,
                   universe_from=how)
        return rec
    try:
        if s.get("provenance") == "gemma":
            g, nrows, stat = BU.fullstats_set(s, None)
            out = with_set("fullstats", g, "fullstats__" + name, s["file"], f"Gemma full results table ({nrows} rows with '{stat}')")
        elif s.get("platform") and os.path.basename(s["platform"]).split(".")[0] in BU.PLATFORM_FILES:
            gpl = os.path.basename(s["platform"]).split(".")[0]
            g = BU.platform_set(gpl)
            out = with_set("platform", g, "platform__" + gpl, gpl, f"manifest platform {gpl}")
        else:
            out = None
            if base and base["type"] == "fullstats" and s.get("file") and s.get("columns"):
                g, nrows, stat = BU.fullstats_set(s, None)
                if len(g) >= max(3 * len(L), 5000):
                    out = with_set("fullstats", g, "fullstats__" + name, s["file"],
                                   f"the list's own table with threshold filters removed ({nrows} rows with '{stat}')")
            if out is None:
                if base is None:
                    rec.update(type="exclude", source=None, gene_set=None, n_membership=None, reason="no base universe",
                               universe_from="none")
                    out = rec
                else:
                    out = dict(base, study=study_key, list_genes=len(L), universe_from=f"inherited from {pl_s['primary'][0]}")
                    if base.get("gene_set"):
                        gs = os.path.join(data_dir(inp), base["gene_set"]) if base["gene_set"].startswith("universe_sets/") else \
                            os.path.join(HERE, "sets", os.path.basename(base["gene_set"]))
                        out["list_genes_outside"] = len(L - read_set_file(gs))
    except Exception as e:  # a table that cannot be rebuilt: inherit
        out = dict(base or {}, study=study_key, list_genes=len(L), universe_from=f"inherited ({type(e).__name__}: {e})")
        if not base:
            out.update(type="exclude", gene_set=None, n_membership=None, reason="no base universe")
    cache[name] = out
    return out


def read_set_file(path):
    return {l.strip() for l in open(path, encoding="utf-8") if l.strip() and not l.startswith("#")}


def set_path(inp, gs):
    return os.path.join(data_dir(inp), gs) if gs.startswith("universe_sets/") else os.path.join(SETS, os.path.basename(gs))


def union_record(names, recs, study_key, label, inp, rule):
    """One record shared by every member of a union (recurrence.load_lists needs one universe per study)."""
    if any(not r.get("gene_set") for r in recs):
        # list_genes counts the union of the members' lists, as every other union record does (it had kept the first
        # member's count: Iwamoto 2005 67 against 76 in the two lists)
        L = set().union(*(read_genes(inp, n) for n in names))
        r0 = dict(recs[0], type="exclude", gene_set=None, n_membership=None, list_genes=len(L), list_genes_outside=None,
                  reason="a union member has no gene set", universe_from="union: member without gene set")
        return r0
    paths = [set_path(inp, r["gene_set"]) for r in recs]
    if len(set(paths)) == 1 and len({r["n_k"] for r in recs}) == 1:
        L = set().union(*(read_genes(inp, n) for n in names))
        return dict(recs[0], list_genes=len(L), list_genes_outside=len(L - read_set_file(paths[0])),
                    universe_from="union: all members share one universe (" + recs[0].get("universe_from", "") + ")")
    genes = set().union(*(read_set_file(p) for p in paths))
    fname = f"union__{rule}__{study_key}"
    os.makedirs(SETS, exist_ok=True)
    with open(os.path.join(SETS, fname + ".txt"), "w", encoding="utf-8") as fh:
        fh.write(f"# union of the universes of {', '.join(names)}\n# script: primary_rule.py select, run {TODAY}\n")
        fh.write("\n".join(sorted(genes)) + "\n")
    types = sorted({r["type"] for r in recs})
    L = set().union(*(read_genes(inp, n) for n in names))
    return {"list_genes": len(L), "list_genes_outside": len(L - genes), "study": study_key, "label": label, "type": types[0] if len(types) == 1 and types[0] != "fullstats" else "union",
            # members' sources split and de-duplicated: a member can itself be an earlier union record (promote_rule.py copies
            # rule A's into universes.json), and joining whole strings grew the source by one entry per rerun (2026-09-27)
            "source": "; ".join(sorted({s for r in recs for s in str(r.get("source")).split("; ")})), "gene_set": "rule_sets/" + fname + ".txt",
            "n_membership": len(genes), "n_k": len(genes), "n_k_from": "n_membership", "n_tested": None,
            "universe_from": "union of " + ", ".join(sorted(set(os.path.basename(p) for p in paths))),
            "evidence": "union of member universes"}


# ---------------------------------------------------------------- select
def cmd_select(inp, code):
    mods = load_code(code)
    extract, BU, H = mods["extract"], mods["build_universes"], mods["h1"]
    ip = inputs_paths(inp)
    extract.HGNC_FILE = ip["hgnc"]
    extract.HGNC = extract.load_hgnc()
    # build_universes reads platforms/ and raw/ through extract.HERE (live, read only); manifests come from `inp`
    extract.HERE = BU.HERE = GL
    cands, pl, specs = candidates(inp, H)
    uni = json.load(open(ip["uni"]))
    cache = {}
    rows_c, rows_r = [], []
    picks = {r: {} for r in ("R1", "R2", "R3", "A")}
    for k, st in sorted(cands.items(), key=lambda x: x[1]["label"]):
        for u in st["units"]:
            rows_c.append({"study": k, "label": st["label"], "lists": ";".join(u["lists"]), "current_primary": "yes" if u["is_current"] else "",
                           "role": u["role"], "modality": u["modality"], "eligible": "yes" if u in st["eligible"] else "",
                           "provenance": u["provenance"], "molecule": u["molecule"], "tissue": u["tissue"],
                           "pooled_diagnoses": "yes" if u["pooled"] else "", "subset": u["subset"],
                           "adjusted_model": "yes" if u["adjusted"] else "", "threshold_kind": u["threshold_kind"],
                           "region": u["region"], "population": u["population"],
                           "n_subjects": "" if u["n_subjects"] is None else u["n_subjects"], "n_genes": u["n_genes"],
                           "manifest_index": u["manifest_index"], "threshold": u["threshold"][:200], "contrast": u["contrast"][:200]})
        cur = [u for u in st["eligible"] if u["is_current"]]
        row = {"study": k, "label": st["label"], "tissue": st["tissue"], "n_units": len(st["units"]),
               "current_genes": cur[0]["n_genes"] if cur else "",
               "n_eligible": len(st["eligible"]), "current": ";".join(cur[0]["lists"]) if cur else "(not eligible: " + ";".join(
                   n for u in st["units"] if u["is_current"] for n in u["lists"]) + ")"}
        for r in ("R1", "R2", "R3", "A"):
            left, d = rule_pick(st, r)
            picks[r][k] = left
            names = sorted(n for u in left for n in u["lists"])
            row[r] = ";".join(names)
            row[f"{r}_decided_by"] = d
            row[f"{r}_same_as_current"] = "yes" if cur and names == sorted(cur[0]["lists"]) else "no"
            row[f"{r}_k"] = len(left)
            row[f"{r}_genes"] = len(set().union(*(read_genes(inp, n) for n in names))) if names else ""
        left, d = rule_pick(st, "R1g")
        row["R1_genes_tiebreak"] = ";".join(sorted(n for u in left for n in u["lists"]))
        row["R1_genes_tiebreak_same_as_R1"] = "yes" if row["R1_genes_tiebreak"] == row["R1"] else "no"
        rows_r.append(row)
    write_tsv(os.path.join(HERE, "candidates.tsv"), rows_c,
              [f"# primary_rule.py select, run {TODAY}; inputs {inp}",
               "# one row per candidate unit (up and down lists of one contrast merged) of every in-corpus study with a primary list; "
               "attributes parsed from the manifest list id, contrast and threshold text (regexes in primary_rule.py); eligible = role "
               "first_look and modality DE"])
    write_tsv(os.path.join(HERE, "rule_choices.tsv"), rows_r,
              [f"# primary_rule.py select, run {TODAY}",
               "# current = the proposed primary in primary_lists.json (snapshot primary_lists_snapshot.json); R1 one list; R2 union "
               "members; R3 weighted members (weight 1/R3_k each); *_decided_by = the last criterion that removed a candidate"])
    # characterization
    dims = characterize(cands)
    # per-rule primary_lists.json and universes.json
    for r in RULES:
        os.makedirs(os.path.join(HERE, r), exist_ok=True)
        newpl = json.loads(json.dumps(pl))
        newuni = json.loads(json.dumps(uni))
        changed = []
        if r != "current":
            for k, left in picks[r].items():
                st, s = cands[k], newpl["studies"][k]
                if not left:
                    continue
                names = [n for u in left for n in u["lists"]]
                if sorted(names) == sorted(s["primary"]):
                    # unchanged lists still need their shared union universe (the default primary lists are rule A's,
                    # so A finds nothing to change while build_universes.py has just written per-list records)
                    if r != "R3" and len(names) > 1:
                        recs = [universe_record(n, k, pl["studies"][k], uni["lists"], specs, BU, extract, inp, cache) for n in names]
                        ur = union_record(names, recs, k, s["label"], inp, r)
                        for n in names:
                            newuni["lists"][n] = dict(ur)
                    continue
                changed.append(k)
                allnames = list(dict.fromkeys(list(s["primary"]) + list(s.get("alternates", []))))
                s["primary"] = names
                s["alternates"] = [n for n in allnames if n not in names]
                s["provenance"] = left[0]["provenance"]
                kinds = [u["threshold_kind"] for u in left]
                s["threshold_kind"] = max(set(kinds), key=lambda x: (kinds.count(x), -kinds.index(x)))
                s["reason"] = f"{r} (primary_rule.py): " + ", ".join(names)
                s["confidence"] = "rule"
                recs = [universe_record(n, k, pl["studies"][k], uni["lists"], specs, BU, extract, inp, cache) for n in names]
                if r != "R3" and len(names) > 1:
                    ur = union_record(names, recs, k, s["label"], inp, r)
                    for n in names:
                        newuni["lists"][n] = dict(ur)
                else:
                    for n, rc in zip(names, recs):
                        newuni["lists"][n] = rc
        if r == "current":  # the default lists are rule A's (several per study): members share one union universe
            for k, s in newpl["studies"].items():
                if s.get("in_corpus", True) and len(s.get("primary", [])) > 1 and k in cands:
                    names = list(s["primary"])
                    recs = [universe_record(n, k, pl["studies"][k], uni["lists"], specs, BU, extract, inp, cache) for n in names]
                    ur = union_record(names, recs, k, s["label"], inp, "current")
                    for n in names:
                        newuni["lists"][n] = dict(ur)
        if r == "R3":  # the candidate units (up and down lists of one contrast together) behind each study's weight 1/k
            for k, left in picks[r].items():
                if left:
                    newpl["studies"][k]["rule_units"] = [u["lists"] for u in left]
        newpl["rule"] = {"rule": r, "text": RULE_TEXT.get(r, "the proposed primary lists as they are"), "studies_changed": changed,
                         "modalities": list(A_MODALITIES) if r == "A" else ["DE"],
                         "generated": TODAY, "script": "primary_rule/primary_rule.py select",
                         "source_primary_lists_sha256": sha(ip["pl"])}
        newuni["rule"] = newpl["rule"]["rule"]
        json.dump(newpl, open(os.path.join(HERE, r, "primary_lists.json"), "w"), indent=1, ensure_ascii=False)
        json.dump(newuni, open(os.path.join(HERE, r, "universes.json"), "w"), indent=1, ensure_ascii=False)
        print(r, "studies changed:", len(changed))
    # universe records for lists new to universes.json
    write_tsv(os.path.join(HERE, "new_universes.tsv"),
              [{"list": n, "type": c.get("type"), "gene_set": c.get("gene_set"), "n_membership": c.get("n_membership"),
                "n_k": c.get("n_k"), "list_genes": c.get("list_genes"), "list_genes_outside": c.get("list_genes_outside", ""),
                "universe_from": c.get("universe_from")} for n, c in sorted(cache.items())],
              [f"# primary_rule.py select, run {TODAY}: universe records built for lists that universes.json (snapshot) has no entry for"])
    # input provenance; a copy of the primary_lists.json this select read, so the rule files can be traced to it
    # after the default file is rewritten (rule A is promoted over it)
    shutil.copyfile(ip["pl"], os.path.join(HERE, "primary_lists_snapshot.json"))
    files = [ip["pl"], ip["uni"], ip["tech"], ip["lineage"], ip["hgnc"], ip["prior"], os.path.join(data_dir(inp), "manifest.json")] + \
        sorted(glob.glob(os.path.join(data_dir(inp), "manifests", "*.json"))) + \
        [os.path.join(code, f) for f in ("recurrence.py", "h1.py", "build_universes.py", "extract.py")]
    write_tsv(os.path.join(HERE, "inputs_sha256.tsv"),
              [{"file": os.path.relpath(f, inp) if f.startswith((inp, data_dir(inp))) else f, "sha256": sha(f)} for f in files],
              [f"# inputs of primary_rule.py (select and run), {TODAY}; inputs dir {inp}; code dir {code}; lists/*.tsv and "
               "universe_sets/*.txt were read from the same inputs dir"])
    return dims


def characterize(cands):
    """choice_dimensions.tsv: per study with more than one eligible candidate, which dimensions vary and where the
    current primary sits on each."""
    rows = []
    for k, st in sorted(cands.items(), key=lambda x: x[1]["label"]):
        el = st["eligible"]
        if len(el) < 2:
            continue
        cur = [u for u in el if u["is_current"]]
        if not cur:
            continue
        c = cur[0]
        el_all = el
        row = {"study": k, "label": st["label"], "tissue": st["tissue"], "n_eligible": len(el),
               "current": ";".join(c["lists"]), "confidence": st["confidence"]}
        varying = []
        # provenance across all eligible candidates; every other dimension among candidates of the current's provenance
        pv = [u["provenance"] for u in el_all]
        if len(set(pv)) > 1:
            varying.append("provenance")
            row["provenance"] = c["provenance"]
            row["provenance (current is preferred value)"] = "yes" if c["provenance"] == min(pv, key=lambda x: RANK.get(x, 3)) else "no"
        el = [u for u in el_all if u["provenance"] == c["provenance"]]
        row["n_same_provenance"] = len(el)

        def dim(name, key, pref=None):
            vals = [key(u) for u in el]
            if len(set(map(str, vals))) > 1:
                varying.append(name)
                row[name] = str(key(c))
                if pref is not None:
                    best = pref(vals)
                    row[name + " (current is preferred value)"] = "yes" if key(c) == best else "no"
        dim("molecule", lambda u: u["molecule"], lambda v: "gene" if "gene" in v else None)
        dim("tissue", lambda u: u["tissue"], lambda v: st["tissue"])
        dim("pooled diagnoses", lambda u: "pooled" if u["pooled"] else "SZ-only", lambda v: "SZ-only")
        dim("subset or condition", lambda u: u["subset"] or "all subjects", lambda v: "all subjects" if "all subjects" in v else None)
        dim("threshold kind", lambda u: u["threshold_kind"], lambda v: "fdr" if "fdr" in v else None)
        same_kind = [u for u in el if u["threshold_kind"] == c["threshold_kind"] and u["provenance"] == c["provenance"]]
        if len({u["threshold"] for u in same_kind}) > 1 and len({u["region"] for u in same_kind}) == 1 and \
                len({u["population"] for u in same_kind}) == 1 and len({u["subset"] for u in same_kind}) == 1:
            varying.append("threshold cut (same kind)")
            row["threshold cut (same kind)"] = c["threshold"][:80]
        dim("region", lambda u: u["region"], lambda v: "prefrontal" if "prefrontal" in v else None)
        # a category shared by two or more candidates means they differ on it (different cell types, cohorts or models)
        dim("cell type, layer, nucleus or stage", lambda u: u["lists"][0] if u["population"].startswith("cell") else "")
        dim("cohort or dataset", lambda u: u["lists"][0] if u["population"].startswith("cohort") else "")
        dim("statistical method or model", lambda u: u["lists"][0] if u["population"].startswith("statistical") or u["adjusted"] else "")
        ns = [u["n_subjects"] for u in el]
        if all(x is not None for x in ns) and len(set(ns)) > 1:
            varying.append("subjects")
            row["subjects"] = c["n_subjects"]
            row["subjects (current is preferred value)"] = "yes" if c["n_subjects"] == max(ns) else "no"
        if len(el) > 1:
            gs = [u["n_genes"] for u in el]
            row["genes (current is largest)"] = "yes" if c["n_genes"] == max(gs) else "no"
            row["current is first in extraction order"] = "yes" if c["manifest_index"] == min(u["manifest_index"] for u in el) else "no"
        row["dimensions varying"] = "; ".join(varying)
        cites = []
        for lab, rx in (("rule 1", r"(?i)rule 1\b"), ("rule 2", r"(?i)rule 2\b"), ("rule 4", r"(?i)rule 4\b"), ("listed first", r"listed first"),
                        ("larger list or sample", r"larger|largest|sample size"), ("title names it", r"title names"),
                        ("chosen because the paper reports most DE there", r"chosen because the paper reports most"),
                        ("the paper reports no DE in the other stratum", r"reports no DE in day"),
                        ("comparability with another study", r"comparab|as for Ketchesin"), ("latest stage", r"latest stage"),
                        ("the list the paper carries forward", r"the list the paper carries forward"),
                        ("the first count the paper states", r"the first count the paper states"),
                        ("published over Gemma", r"over the ([\w-]+ )?Gemma|beats the Gemma"),
                        ("baseline or all-subject contrast", r"(?i)baseline|combined-sex|all-SZ list|all patients|affected twins vs controls"),
                        ("gene list over non-coding lists", r"mRNA list chosen|Gene list \(|Gene list at"),
                        ("new samples over reanalysed data", r"New CCHPC|new data|first-look|First-look")):
            if re.search(rx, st["reason"]):
                cites.append(lab)
        row["reason cites"] = "; ".join(cites)
        row["reason"] = st["reason"][:300]
        rows.append(row)
    fields = ["study", "label", "tissue", "n_eligible", "n_same_provenance", "current", "confidence", "dimensions varying"]
    fields += [f for r in rows for f in r if f not in fields and f not in ("reason cites", "reason")]
    fields = list(dict.fromkeys(fields)) + ["reason cites", "reason"]
    write_tsv(os.path.join(HERE, "choice_dimensions.tsv"), rows,
              [f"# primary_rule.py select, run {TODAY}: studies with 2 or more eligible candidates (first_look, DE); provenance is "
               "compared across all of them, every other dimension among the candidates of the current primary's provenance "
               "(n_same_provenance); for each dimension on which they differ, the current primary's value and whether it is the "
               "value R1 prefers (blank where R1 has no preference); 'reason cites' = phrases found in the primary_lists.json reason"], fields)
    return rows


def write_tsv(path, rows, head, fields=None):
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    return [l for l in lines if l.startswith("#")], list(csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("select", "run", "check", "report"))
    ap.add_argument("--inputs", default=GL)
    ap.add_argument("--code", default=GL)
    ap.add_argument("--rules", default=",".join(RULES))
    a = ap.parse_args()
    if a.cmd == "select":
        cmd_select(os.path.abspath(a.inputs), os.path.abspath(a.code))
    elif a.cmd == "run":
        from primary_rule_run import cmd_run
        cmd_run(os.path.abspath(a.inputs), os.path.abspath(a.code), a.rules.split(","))
    elif a.cmd == "check":
        from primary_rule_run import cmd_check
        cmd_check(os.path.abspath(a.inputs), os.path.abspath(a.code))
    else:
        from primary_rule_report import cmd_report
        cmd_report()
