#!/usr/bin/env python3
"""Gene recurrence across studies against chance, one list per study.

Unit of analysis (tracker, "Unit of analysis for recurrence"): each study contributes its
proposed primary SZ vs control list (primary_lists.json; up and down lists of one contrast
count once). Studies that share data (same_data_as) are collapsed to one representative per
scenario, so a dataset reanalyzed five times counts once.

Left out before collapsing (both tiers): meta-lane studies (lineage lane 'meta'), and studies
whose same_data_as evidence names two or more GEO series that they share with different
partners (two partners whose shared series are disjoint): these are multi-dataset reanalyses
that behave like meta-analyses. Union-find through them chained unrelated datasets into one
group (29 studies represented by a blood study, analysis review R1, 2026-09-24). A pair that
shares two series at once (a SuperSeries and its SubSeries) is one dataset and does not count.

List role (manifest "role"; Paul, 2026-09-24: keep 'first look at the data' separate): only first_look lists
(the study's own new data, the default) can be a study's primary list. A primary list marked reanalysis (the
study's analysis of data first published elsewhere) is dropped, and a study left with no primary list is left
out. Lineage group 3 studies (no new data) are left out as before.

Representative (per scenario, among group members that pass the scenario's filter and have a
list): published over reproduced over gemma, then FDR lists, then earliest PMID. A blood study
therefore never represents a brain scenario.

Chance model (universes.json, from build_universes.py): list k has a tested gene set M_k and a
denominator n_k (the paper's stated number of genes tested when it is stated in genes and not
larger than |M_k|, else |M_k|). The list L_k is restricted to M_k, and gene g is drawn into
list k with probability p_gk = |L_k & M_k| / n_k if g is in M_k, else 0. This assumes every
tested gene is equally likely to be in a list. The number of lists containing g is
Poisson-binomial over k. Lists whose universe type is 'exclude', or 'targeted' without a known
gene set, are left out of this model.

DE-prior null (sensitivity): the same total per list (|L_k & M_k| x |M_k| / n_k), shared out in
proportion to w_g = n_g + 1, where n_g is the gene's count in the DE prior in force (DE_PRIOR_SOURCE); genes not in
the prior get the median weight; probabilities capped at 1 with the excess shared among the rest. Default since
2026-09-30 (Paul's decision): our own 2-fold build from Gemma, de_prior_lists/our_prior/1_global_our_prior.tsv (all
technologies, 8,712 contrasts), n_g = n_de_total, the number of those contrasts in which the gene is DE. From
2026-09-29 to that decision the default was the external 4-fold global prior of de_prior_lists/de_prior_ranks/
1_global_de_prior.tsv (5,785 contrasts), which stays available as DE_PRIOR_SOURCE=global. It is the analogue of the column used before, N_HitLists of the Crow et al. 2019 prior (the number
of 635 DE hit lists containing the gene; data/reference/DE_Prior.txt), which stays available as DE_PRIOR_SOURCE=crow2019
(environment variable) so the earlier results can be regenerated. Genes are keyed by NCBI gene ID through HGNC, with the
prior's own symbol as the fallback.
Use (Paul, 2026-09-29): the DE-prior null gives aggregate expectations only (expected genes in at least N lists, the
summary rows of recurrence.tsv, the gene_coverage curve). Per-gene significance is reported under the uniform and
expression-matched nulls; p_prior, q_prior_bh and q_cluster_prior_bh stay in recurrence_genes.tsv but are not evidence
for any gene. Reason: on bipolar disorder and major depression postmortem brain lists from Gemma, the global DE-prior
null calls genes the uniform null does not, nearly all of them low prior (DE_PRIOR_GENE_NOTE,
de_prior_lists/checks_2026-09-29/checks_report.md).

Sensitivity columns:
  u18_*  the same modelled lists, unrestricted, each a draw of its own size from U = 18,000 genes
         (isolates the universe model; observed counts differ from the main columns only by the
         list genes outside M_k);
  old_*  every representative list in the scenario, including lists the universe model leaves
         out (miRNA, lncRNA, cell-type-stacked, platform unknown), unrestricted, U = 18,000. Not a
         reproduction of any earlier run.

Cohort tier (conservative): each same_data_as group is assigned the named lineage lanes of ALL
its members (cmc, mssm, pitt, libd merged; stanley one cluster; other named lanes one cluster
each); groups whose members span several named clusters merge those clusters; a group with no
named lane is its own cluster. A gene counts once per cluster if any member list has it, with
P(g in cluster c) = 1 - prod over member lists k of (1 - p_gk); Poisson-binomial over clusters.

Multi-bank studies (lineage node 'banks', e.g. Mendez-Victoriano 2026, NSW TRC and Stanley; Paul, 2026-09-24): their
lane does not name or merge a cluster. In the cohort tier they contribute their per-bank lists (manifest "bank") to the
matching bank clusters; without per-bank lists they stay in the study tier and are left out of the cohort tier.

Direction: a gene in both an up and a down list of one study (or in up and down rows of one
list) has direction 'both' and is left out of the sign test.

Forking-paths switches (2026-09-30; tracker Decisions log 'forking paths', 'recurrence by tiers of analytic rigour',
'lists with nothing significant are included', 'recurrence list thresholds as tiers'). Environment variables, read at
import as DE_PRIOR_SOURCE is; each default is the behaviour before the switch existed, so the default outputs do not move.
analyses/forking_paths/ varies them in one process by setting the module globals.
  RECURRENCE_THRESHOLD_TIER  all (default) | nominal | fdr. fdr: studies whose threshold_kind is fdr, plus non-FDR studies
                             whose primary list has a computed BH table (computed_fdr_2026-09-30/, routes complete,
                             authors_q, truncated), which then count with their q < 0.05 rows in place of the paper's
                             list (such a list counts as modality DE whatever the paper's list was: Dean 2026's
                             gene selection). nominal: fdr and nominal studies as published, plus NOMINAL_TIER_OTHER.
  RECURRENCE_EMPTY_LISTS     drop (default) | keep. keep: a representative whose list holds no gene still counts (it
                             tested every gene in its universe and found none); drop: it cannot represent its group.
  PRINTED_NAME_FALLBACK      keep (default) | drop. drop: rows extract.py tagged ';tag:printed-name-fallback' (probe on the
                             list's platform, which gives it no gene; symbol from the printed name) hold no gene.
  RECURRENCE_PSEUDO_NOEXPR   keep (default) | drop. drop: HGNC pseudogenes (locus_group 'pseudogene') and genes no brain
                             source measures as expressed (EXPRESSED_FILE: LIBD DLPFC or hippocampus, PsychAD pseudobulk
                             or any Gemma brain study; lenient definition, Paul 2026-09-30) leave every list and universe.
                             Until 2026-09-30 'expressed' was a LIBD BrainSeq value (EXPR_FILE) alone.
  RECURRENCE_NONCODING       pooled (default) | drop | only. drop: HGNC locus_group 'non-coding RNA' leaves every list and
                             universe; only: they are all that is kept (non-coding RNAs analysed on their own).
                             A universe that loses genes to either filter keeps its rate: n_k is scaled by |M_k'| / |M_k|.
  RECURRENCE_ADD             comma list, default empty. meta: meta-lane studies enter as studies; reanalysis: lineage
                             group 3 studies (re-publications) enter, their reanalysis-role lists allowed. Both enter as
                             their own same_data_as groups (counted as independent, not collapsed onto their sources).

    python3 recurrence.py            -> recurrence.tsv, recurrence_genes.tsv, recurrence_units.tsv,
                                        recurrence_report.md (written from the tsv files)
    python3 recurrence.py --report   -> recurrence_report.md only, from the existing tsv files
"""
import collections
import csv
import datetime
import glob
import hashlib
import itertools
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")  # repo data/: lists, manifests, reference, raw, universe sets
U = 18000
KS = (2, 3, 5, 10)
KS_COHORT = (2, 3, 5)
LINEAGE = os.path.join(HERE, "..", "corpus", "lineage.json")
# DE prior: "ours" is the default. Each source: file, count column (the weight, n + 1),
# rank column (normalized, 1 = most often DE), ID column (NCBI gene ID) and symbol column.
DE_PRIOR_SOURCES = {
    "global": dict(file=os.path.join(HERE, "de_prior_lists", "de_prior_ranks", "1_global_de_prior.tsv"), count="n_de_total",
                   rank="de_prior_norm", id="ncbi_gene_id", symbol="gene_symbol",
                   label="global DE prior (Gemma, all technologies, 5,785 contrasts; de_prior_lists/de_prior_ranks/1_global_de_prior.tsv), n_de_total",
                   of="5,785 contrasts"),
    # our own build from Gemma (2026-09-30; de_prior_lists/our_prior/README.md): 2-fold, per-gene n_tested, genes never DE kept
    # with n_de_total 0. The default.
    "ours": dict(file=os.path.join(HERE, "de_prior_lists", "our_prior", "1_global_our_prior.tsv"), count="n_de_total",
                 rank="de_prior_norm", id="ncbi_gene_id", symbol="gene_symbol",
                 label="our global DE prior (Gemma, all technologies, 8,712 contrasts, 2-fold; de_prior_lists/our_prior/1_global_our_prior.tsv), n_de_total",
                 of="8,712 contrasts"),
    "crow2019": dict(file=os.path.join(DATA, "reference", "DE_Prior.txt"), count="N_HitLists", rank="DE_Prior_Rank",
                     id="Gene_EntrezID", symbol="Gene_Name",
                     label="Crow et al. 2019 DE prior (github.com/maggiecrow/DEprior; data/reference/DE_Prior.txt), N_HitLists",
                     of="635 DE hit lists"),
}
DE_PRIOR_SOURCE = os.environ.get("DE_PRIOR_SOURCE", "ours")
# Paul, 2026-09-29, after the checks in de_prior_lists/checks_2026-09-29/: the DE-prior null is used for aggregate expectations
# only; per-gene q under it stays in the outputs and is not cited as evidence for any gene. Numbers from checks_report.md.
DE_PRIOR_GENE_NOTE = ("the DE-prior null is used for aggregate expectations only; its per-gene p and q (p_prior, q_prior_bh, "
                      "q_cluster_prior_bh) are not evidence for any gene (Paul, 2026-09-29): on Gemma's bipolar disorder and major "
                      "depression postmortem brain lists (FDR 0.05), the global DE-prior null gives 350 genes at q < 0.1 against 1 "
                      "under the uniform null, 349 of them low prior (de_prior_lists/checks_2026-09-29/checks_report.md)")
DE_PRIOR = DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]["file"]
# forking-paths switches (module docstring); defaults = the behaviour before 2026-09-30
THRESHOLD_TIER = os.environ.get("RECURRENCE_THRESHOLD_TIER", "all")
EMPTY_LISTS = os.environ.get("RECURRENCE_EMPTY_LISTS", "drop")
PRINTED_NAME_FALLBACK = os.environ.get("PRINTED_NAME_FALLBACK", "keep")
PSEUDO_NOEXPR = os.environ.get("RECURRENCE_PSEUDO_NOEXPR", "keep")
NONCODING = os.environ.get("RECURRENCE_NONCODING", "pooled")
ADD_STUDIES = os.environ.get("RECURRENCE_ADD", "")
# n_k rule (CODE-001, code review 2026-10-11): "stated" (default) = n_k from universes.json, the paper's stated number of
# genes tested when it is below |M_k|; "membership" = n_k = |M_k|, so p_gk = |L_k & M_k| / |M_k| and each list's expected
# size equals its size. The earlier results are the default.
NK_RULE = os.environ.get("RECURRENCE_NK", "stated")
# Jaffe 2018 fork (CODE-007, code review 2026-10-11): comma-separated list names whose rows are kept only where the
# feature type (column `group`) is Gene; exons, expressed regions and junctions are dropped. Default: none.
GENELEVEL_ONLY = os.environ.get("RECURRENCE_GENELEVEL_ONLY", "")


def _out():
    """Directory the outputs are written to and read from: RECURRENCE_OUT, else this folder (the committed results)."""
    return os.environ.get("RECURRENCE_OUT") or HERE


SWITCH_VALUES = {"NK_RULE": ("stated", "membership"), "THRESHOLD_TIER": ("all", "nominal", "fdr"), "EMPTY_LISTS": ("drop", "keep"),
                 "PRINTED_NAME_FALLBACK": ("keep", "drop"), "PSEUDO_NOEXPR": ("keep", "drop"),
                 "NONCODING": ("pooled", "drop", "only")}
# computed BH q for non-FDR primary lists (c058127); 'lower_bound' (Barnes 2011) is not an exact set and is left out
COMPUTED_FDR_DIR = os.path.join(HERE, "computed_fdr_2026-09-30")
COMPUTED_FDR_ROUTES = ("complete", "authors_q", "truncated")
# the nominal tier's one study of threshold_kind 'other': Shao 2008, ANCOVA-based hits (computed_fdr_2026-09-30/README.md,
# 'Proposed tiers' 2)
NOMINAL_TIER_OTHER = {"18191109"}
HGNC_FILE = os.path.join(DATA, "reference", "hgnc_complete_set.txt")
# Expression-matched null (Paul, 2026-09-28: "a uniform random is not right - it should follow the expression level
# distributions"): LIBD BrainSeq DLPFC mean log CPM per gene (gwas_vs_de/scripts/02a_libd_expr.R)
EXPR_FILE = os.path.join(HERE, "gwas_vs_de", "results", "libd_mean_logcpm.tsv")
# 'expressed in brain' for RECURRENCE_PSEUDO_NOEXPR (lenient, Paul 2026-09-30): expressed_in_brain/expressed_in_brain.py
EXPRESSED_FILE = os.path.join(HERE, "expressed_in_brain", "expressed_in_brain.tsv")
N_EXPR_BINS = 10
# lineage lanes that name one brain bank or line collection; every other lane (brain_other,
# periph, sc) makes each group its own cohort. meta is left out.
NAMED_LANES = {"pitt", "stanley", "libd", "cmc", "mssm", "harvard", "aus", "pritzker", "ipsc"}
# libd joined 2026-09-24 (Paul: "I think we probably have to"): donor_overlap/ finds 44 of 185 HBCC sgACC donors
# (Akula 2021, lane cmc) in LIBD BrainSeq at the exact level (4.4 by chance; HBCC PMIs flagged 'Uncertain' treated as
# missing), LIBD numbers Br1030 to Br2542.
MERGE = {k: "cmc+mssm+pitt+libd" for k in ("cmc", "mssm", "pitt", "libd")}
GENE_SCENARIOS = ("published, brain", "published, brain, FDR lists", "published, brain, microarray", "published, brain, RNA-seq")
# assay technology per study (classify_technology.py; Paul 2026-09-24: "isolated for rna-seq vs microarray (excluding single cell)")
TECHNOLOGY = os.path.join(HERE, "technology.json")
GSE_RE = re.compile(r"GSE\d+")
SCENARIOS = [
    ("all provenances, all tissues", lambda s: True),
    ("published only, all tissues", lambda s: s["provenance"] == "published"),
    ("published, brain", lambda s: s["provenance"] == "published" and s["tissue"] == "brain"),
    ("published, brain, FDR lists", lambda s: s["provenance"] == "published" and s["tissue"] == "brain" and s["threshold_kind"] == "fdr"),
    ("published, brain, nominal lists", lambda s: s["provenance"] == "published" and s["tissue"] == "brain" and s["threshold_kind"] == "nominal"),
    ("published, blood", lambda s: s["provenance"] == "published" and s["tissue"] == "blood"),
    ("published, brain, microarray", lambda s: s["provenance"] == "published" and s["tissue"] == "brain" and s.get("technology") == "microarray"),
    ("published, brain, RNA-seq", lambda s: s["provenance"] == "published" and s["tissue"] == "brain" and s.get("technology") == "rnaseq"),
]
RANK = {"published": 0, "reproduced": 1, "gemma": 2}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


def list_meta():
    """list name -> {"role", "bank", "modality"} from manifest.json and manifests/*.json (role defaults to first_look,
    modality to DE)."""
    specs = json.load(open(os.path.join(DATA, "manifest.json")))["lists"]
    for f in sorted(glob.glob(os.path.join(DATA, "manifests", "*.json"))):
        specs += json.load(open(f))["lists"]
    return {f"{x['pmid']}__{x['list_id']}": {"role": x.get("role", "first_look"), "bank": x.get("bank", ""),
                                             "modality": x.get("modality", "DE")} for x in specs}


def check_switches():
    for k, ok in SWITCH_VALUES.items():
        if globals()[k] not in ok:
            raise SystemExit(f"recurrence.py: {k} = {globals()[k]!r}, not one of {ok}")
    bad = set(filter(None, ADD_STUDIES.split(","))) - {"meta", "reanalysis"}
    if bad:
        raise SystemExit(f"recurrence.py: RECURRENCE_ADD {sorted(bad)} not in meta, reanalysis")


def switches():
    """The non-default switches in force, as 'NAME=value' strings (empty at the defaults)."""
    out = [f"{k}={globals()[k]}" for k, ok in SWITCH_VALUES.items() if globals()[k] != ok[0]]
    return out + ([f"ADD_STUDIES={ADD_STUDIES}"] if ADD_STUDIES else []) + ([f"GENELEVEL_ONLY={GENELEVEL_ONLY}"] if GENELEVEL_ONLY else [])


def read_list(name):
    """symbol -> direction; a symbol with both up and down rows gets 'both'. With PRINTED_NAME_FALLBACK=drop, rows tagged
    printed-name-fallback hold no gene."""
    rows = [l for l in open(os.path.join(DATA, "lists", name + ".tsv"), encoding="utf-8") if not l.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        g = r.get("symbol")
        if g and name in GENELEVEL_ONLY.split(",") and (r.get("group") or "Gene") != "Gene":
            continue
        if g and PRINTED_NAME_FALLBACK == "drop" and "printed-name-fallback" in (r.get("map_route") or ""):
            continue
        if g:
            out[g] = merge_dir(out.get(g), r.get("direction", ""))
    return out


def computed_fdr_lists():
    """list name -> (symbol -> direction) of the rows at computed BH q < 0.05 (computed_fdr_2026-09-30/summary.tsv, routes
    COMPUTED_FDR_ROUTES); a table shared by several lists (up and down lists of one contrast) goes to each."""
    lines = [l for l in open(os.path.join(COMPUTED_FDR_DIR, "summary.tsv"), encoding="utf-8") if not l.startswith("#")]
    out = {}
    for r in csv.DictReader(lines, delimiter="\t"):
        if r["route"] not in COMPUTED_FDR_ROUTES:
            continue
        rows = [l for l in open(os.path.join(COMPUTED_FDR_DIR, r["output"]), encoding="utf-8") if not l.startswith("#")]
        d = {}
        for x in csv.DictReader(rows, delimiter="\t"):
            if x["q_lt_0.05"] == "yes" and x["symbol"]:
                d[x["symbol"]] = merge_dir(d.get(x["symbol"]), x.get("direction", ""))
        for n in (y.strip() for y in r["lists"].split(",")):
            out[n] = {**out.get(n, {}), **d}
    return out


def gene_filter_set():
    """HGNC approved symbols kept by RECURRENCE_PSEUDO_NOEXPR and RECURRENCE_NONCODING, or None when neither filters."""
    if PSEUDO_NOEXPR == "keep" and NONCODING == "pooled":
        return None
    expr = expressed_in_brain() if PSEUDO_NOEXPR == "drop" else set()
    keep = set()
    with open(HGNC_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "Approved":
                continue
            grp = r["locus_group"]
            if PSEUDO_NOEXPR == "drop" and (grp == "pseudogene" or r["symbol"] not in expr):
                continue
            if (NONCODING == "drop" and grp == "non-coding RNA") or (NONCODING == "only" and grp != "non-coding RNA"):
                continue
            keep.add(r["symbol"])
    return keep


def expressed_in_brain():
    """Symbols some brain source measures as expressed (EXPRESSED_FILE)."""
    with open(EXPRESSED_FILE, encoding="utf-8") as fh:
        return {r["symbol"] for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}


def tier_ok(k, units):
    """RECURRENCE_THRESHOLD_TIER on a study's recorded threshold_kind (primary_lists.json)."""
    kind = units["studies"][k]["threshold_kind"]
    if THRESHOLD_TIER == "nominal":
        return kind in ("fdr", "nominal") or k in NOMINAL_TIER_OTHER
    if THRESHOLD_TIER == "fdr":
        return kind == "fdr" or k in units.get("computed_fdr", {})
    return True


def usable(units, lists):
    """Which studies can represent their group: a list with a gene (or any list, RECURRENCE_EMPTY_LISTS=keep) in the
    threshold tier in force."""
    return lambda k: (bool(lists[k]) or EMPTY_LISTS == "keep") and tier_ok(k, units)


def merge_dir(a, b):
    if a is None or a == b:
        return b
    if {a, b} <= {"up", "down", "both"}:
        return "both"
    return a if a in ("up", "down", "both") else b


def read_set(path):
    return {l.strip() for l in open(os.path.join(DATA, path), encoding="utf-8") if l.strip() and not l.startswith("#")}


def poisson_binomial(ps):
    dist = [1.0]
    for q in ps:
        nd = [0.0] * (len(dist) + 1)
        for i, v in enumerate(dist):
            nd[i] += v * (1 - q)
            nd[i + 1] += v * q
        dist = nd
    return dist


def pb_matrix(P):
    """Row-wise Poisson-binomial: P is genes x units of inclusion probabilities; returns
    genes x (units + 1) probabilities of 0..units inclusions."""
    G, K = P.shape
    D = np.zeros((G, K + 1))
    D[:, 0] = 1.0
    for j in range(K):
        q = P[:, j:j + 1]
        D[:, 1:j + 2] = D[:, 1:j + 2] * (1 - q) + D[:, 0:j + 1] * q
        D[:, 0:1] = D[:, 0:1] * (1 - q)
    return D


def upper_tail(D, obs):
    """P(X >= obs) per row."""
    c = np.cumsum(D[:, ::-1], axis=1)[:, ::-1]  # c[:, k] = P(X >= k)
    return np.clip(c[np.arange(len(obs)), obs], 0.0, 1.0), c


def bh(p):
    p = np.asarray(p)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def sign_test(up, down):
    n = up + down
    if n == 0:
        return None
    k = min(up, down)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


# ---------------------------------------------------------------- units: exclusion, groups, clusters
def partner_series(s):
    return [set(GSE_RE.findall(o.get("evidence", "") if isinstance(o, dict) else "")) for o in s.get("same_data_as", [])]


def multi_series(s):
    """True when the study shares GEO series with different partners (two partners' shared series disjoint)."""
    sets = [x for x in partner_series(s) if x]
    if len(set().union(*sets) if sets else set()) < 2:
        return False
    return any(not (a & b) for a, b in itertools.combinations(sets, 2))


def load_units(pl=None, nodes=None):
    """Study selection shared by recurrence.py, h1.py, confounds/scripts/scenario_universe.py and
    effect_sizes/harmonize.py (grouping only). Returns a dict:
      all       in-corpus studies with a primary list
      studies   those left after meta-lane and multi-series studies are removed
      meta, multi, reanalysis  pmids left out
      role_dropped  pmid -> primary lists dropped because their role is not first_look
      role_only     pmids left out because every primary list was dropped for its role
      multibank     pmid -> banks (lineage 'banks', 2 or more) for studies in `studies`
      bank_lists    pmid -> {bank: [primary lists with that manifest bank]} (multi-bank studies with per-bank lists)
      groups    list of sorted member lists (same_data_as, links only among `studies`)
      group_of  pmid -> group index
      gcluster  group index -> named cluster (None when no member has a named lane)
      merged    list of cluster sets merged because one group spans them
      lanes"""
    if pl is None:
        pl = json.load(open(os.path.join(HERE, "primary_lists.json")))
    if nodes is None:
        nodes = json.load(open(LINEAGE))["nodes"]
    check_switches()
    add = set(filter(None, ADD_STUDIES.split(",")))
    lanes = {(n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")): n.get("lane") for n in nodes}  # DOI key for nodes without a PMID (Cui 2024)
    cats = {(n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")): n.get("category") for n in nodes}
    lm = list_meta()
    cf_lists = set(computed_fdr_lists()) if THRESHOLD_TIER == "fdr" else set()
    # only first_look lists represent a study (reanalysis lists are a separate category, Paul 2026-09-24)
    # only DE lists are overlapped (Paul 2026-09-24: "we should only be overlapping DE lists in the primary analysis");
    # manifest `modality` marks the others (e.g. Oldmeadow 2014's exon-splicing table)
    # a primary-list rule may widen the modalities (primary_rule rule A, Paul 2026-09-25); the default is DE only
    OK_MOD = set(pl.get("rule", {}).get("modalities", ["DE"]))
    allst, role_dropped, role_only, role_alternates, modality_dropped, modality_only = {}, {}, [], {}, {}, []
    for k, v in pl["studies"].items():
        if not (v.get("in_corpus", True) and v.get("primary")):
            continue
        alt = [n for n in v.get("alternates", []) if lm.get(n, {}).get("role", "first_look") != "first_look"]
        if alt:
            role_alternates[k] = alt
        # RECURRENCE_ADD=reanalysis: a group 3 study's reanalysis-role lists represent it
        bad = [] if ("reanalysis" in add and cats.get(k) == "3") else \
            [n for n in v["primary"] if lm.get(n, {}).get("role", "first_look") != "first_look"]
        if bad:
            role_dropped[k] = bad
            v = dict(v, primary=[n for n in v["primary"] if n not in bad])
            if not v["primary"]:
                role_only.append(k)
                continue
        # RECURRENCE_THRESHOLD_TIER=fdr: a list with a computed BH table counts as DE, since the list it is replaced by is
        # a DE test over the full table (Dean 2026's Table 2 is a gene selection; its computed list is the complete table)
        nde = [n for n in v["primary"] if lm.get(n, {}).get("modality", "DE") not in OK_MOD and n not in cf_lists]
        if nde:
            modality_dropped[k] = nde
            v = dict(v, primary=[n for n in v["primary"] if n not in nde])
            if not v["primary"]:
                modality_only.append(k)
                continue
        allst[k] = v
    tech = json.load(open(TECHNOLOGY))["studies"] if os.path.exists(TECHNOLOGY) else {}
    allst = {k: dict(v, technology=tech.get(k, {}).get("technology", "unknown")) for k, v in allst.items()}
    meta = sorted(k for k in allst if lanes.get(k) == "meta")
    multi = sorted(k for k in allst if k not in meta and multi_series(allst[k]))
    # Reanalyses (lineage group 3: no new expression data) are left out, so each dataset enters once, through
    # the study that generated it (with its reproduced or Gemma list where it published none). Single-cohort
    # reanalyses left the meta lane for display on 2026-09-24 and otherwise re-count their cohort's donors
    # (Clifton 2023 and Xia 2024 re-count CommonMind/PsychENCODE; 2026-09-24).
    reanalysis = sorted(k for k in allst if k not in meta and k not in multi and cats.get(k) == "3")
    # RECURRENCE_ADD: meta-lane studies and / or group 3 re-publications enter, each its own group (independent)
    added = set(meta if "meta" in add else []) | set(reanalysis if "reanalysis" in add else [])
    meta = [k for k in meta if k not in added]
    reanalysis = [k for k in reanalysis if k not in added]
    studies = {k: v for k, v in allst.items() if k not in meta and k not in multi and k not in reanalysis}
    nbanks = {(n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")): list(n["banks"]) for n in nodes if len(n.get("banks") or []) > 1}
    multibank = {k: nbanks[k] for k in studies if k in nbanks}
    bank_lists = {}
    for k, bs in multibank.items():
        per = {b: [n for n in studies[k]["primary"] if lm.get(n, {}).get("bank") == b] for b in bs}
        if any(per.values()):
            bank_lists[k] = {b: v for b, v in per.items() if v}
    parent = {k: k for k in studies}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for k, s in studies.items():
        for o in s.get("same_data_as", []):
            o = o if isinstance(o, str) else o.get("pmid", "")
            if o in parent and k not in added and o not in added:
                parent[find(k)] = find(o)
    gd = {}
    for k in studies:
        gd.setdefault(find(k), []).append(k)
    groups = sorted((sorted(g, key=lambda k: int(k) if k.isdigit() else 10**9) for g in gd.values()), key=lambda g: g[0])
    group_of = {k: i for i, g in enumerate(groups) for k in g}
    # named clusters from all members; merge clusters a group spans
    cpar = {}

    def cfind(x):
        cpar.setdefault(x, x)
        while cpar[x] != x:
            x = cpar[x]
        return x
    gnamed = []
    for g in groups:
        # a multi-bank study's lane does not name or merge clusters (one bridging study must not merge two banks)
        named = sorted({MERGE.get(lanes.get(k), lanes.get(k)) for k in g if lanes.get(k) in NAMED_LANES and k not in multibank})
        gnamed.append(named)
        for c in named:
            cfind(c)
        for a, b in zip(named, named[1:]):
            ra, rb = cfind(a), cfind(b)
            if ra != rb:
                cpar[rb] = ra
    comp = {}
    for c in list(cpar):
        comp.setdefault(cfind(c), set()).add(c)
    cname = {c: "+".join(sorted(members)) for members in comp.values() for c in members}
    gcluster = {i: (cname[named[0]] if named else None) for i, named in enumerate(gnamed)}
    merged = [sorted(m) for m in comp.values() if len(m) > 1]
    return {"all": allst, "studies": studies, "meta": meta, "multi": multi, "reanalysis": reanalysis, "groups": groups, "group_of": group_of,
            "gcluster": gcluster, "merged": merged, "lanes": lanes, "cname": cname, "role_dropped": role_dropped,
            "modality_dropped": modality_dropped, "modality_only": sorted(modality_only),
            "role_only": sorted(role_only), "role_alternates": role_alternates, "labels": {k: v["label"] for k, v in pl["studies"].items()},
            "multibank": multibank, "bank_lists": bank_lists, "added": sorted(added)}


def pick(members, studies, f, ok=lambda k: True):
    """Representative of one group for one scenario (None if no member passes)."""
    c = [k for k in members if f(studies[k]) and ok(k)]
    if not c:
        return None
    return min(c, key=lambda k: (RANK.get(studies[k]["provenance"], 3), studies[k]["threshold_kind"] != "fdr",
                                 int(k) if k.isdigit() else 10**9))


def cluster_of(units, k):
    """Cohort cluster of the group containing study k."""
    c = units["gcluster"][units["group_of"][k]]
    return c if c else f"{units['lanes'].get(k) or 'none'}:{k}"


def bank_cluster(units, b):
    """Cohort cluster of one bank (a multi-bank study's per-bank list goes here)."""
    c = MERGE.get(b, b)
    return units["cname"].get(c, c)


def cohort_label(units, k):
    """recurrence_units.tsv cluster column: the cluster, or how a multi-bank study enters the cohort tier."""
    if k in units["multibank"]:
        bl = units["bank_lists"].get(k)
        if bl:
            return "per-bank: " + "; ".join(f"{b} -> {bank_cluster(units, b)}" for b in bl)
        return f"left out of cohort tier (multi-bank: {'+'.join(units['multibank'][k])}; no per-bank lists)"
    return cluster_of(units, k)


def representatives(units, f, ok=lambda k: True):
    """rep -> group members, for one scenario filter."""
    out = {}
    for g in units["groups"]:
        r = pick(g, units["studies"], f, ok)
        if r:
            out[r] = g
    return out


# ---------------------------------------------------------------- lists and universes
def load_lists(units, uni):
    """lists[k] (symbol -> direction), M[k], nk[k], utype[k] for every study left in units."""
    studies = units["studies"]
    lists, M, nk, utype, nouni = {}, {}, {}, {}, {}
    cache = {}
    for k, s in studies.items():
        d = {}
        for name in s["primary"]:
            for g, v in read_list(name).items():
                d[g] = merge_dir(d.get(g), v)
        lists[k] = d
        u = [uni.get(n) for n in s["primary"]]
        if any(x is None for x in u):  # a list build_universes.py has not seen (R9)
            nouni[k] = "no universes.json entry"
            utype[k] = None
            continue
        utype[k] = u[0]["type"]
        if all(x.get("gene_set") for x in u):
            for x in u:
                if x["gene_set"] not in cache:
                    cache[x["gene_set"]] = read_set(x["gene_set"])
            if len({frozenset(cache[x["gene_set"]]) for x in u}) != 1 or len({x["n_k"] for x in u}) != 1:
                raise SystemExit(f"{k}: primary lists with different universes")
            M[k] = cache[u[0]["gene_set"]]
            nk[k] = u[0]["n_k"]
    # RECURRENCE_THRESHOLD_TIER=fdr: a non-FDR study whose primary lists have computed BH tables counts with their
    # q < 0.05 rows (members without a table drop out of its list); units["computed_fdr"] names the lists used
    units["computed_fdr"] = {}
    if THRESHOLD_TIER == "fdr":
        cf = computed_fdr_lists()
        for k, s in studies.items():
            names = [n for n in s["primary"] if n in cf]
            if s["threshold_kind"] != "fdr" and names:
                d = {}
                for n in names:
                    for g, v in cf[n].items():
                        d[g] = merge_dir(d.get(g), v)
                lists[k] = d
                units["computed_fdr"][k] = names
    # RECURRENCE_PSEUDO_NOEXPR / RECURRENCE_NONCODING: genes outside the kept set leave lists and universes; n_k keeps
    # its ratio to |M_k|; a universe left empty is not modelled
    keep = gene_filter_set()
    if keep is not None:
        for k in list(lists):
            lists[k] = {g: v for g, v in lists[k].items() if g in keep}
            if k in M:
                m0 = M[k]
                M[k] = m0 & keep
                if not M[k]:
                    del M[k], nk[k]
                    nouni[k] = "no gene of its universe left by the gene filter"
                    continue
                nk[k] = nk[k] * len(M[k]) / len(m0)
    if NK_RULE == "membership":
        for k in M:
            nk[k] = len(M[k])
    return lists, M, nk, utype, nouni


def scenario_selection(name):
    """(sel, M, lists, units) for a named scenario: the modelled representative lists. Used by
    confounds/scripts/scenario_universe.py."""
    f = dict(SCENARIOS)[name]
    units = load_units()
    uni = json.load(open(os.path.join(HERE, "universes.json")))["lists"]
    lists, M, nk, _, _ = load_lists(units, uni)
    reps = representatives(units, f, usable(units, lists))
    sel = [k for k in sorted(reps, key=lambda k: (0, int(k), "") if k.isdigit() else (1, 0, k)) if k in M]
    return sel, M, lists, units


def de_prior_table():
    """gene symbol -> (count, normalized rank) from the DE prior in force (DE_PRIOR_SOURCE; the file is DE_PRIOR). Keyed by
    NCBI gene ID through HGNC approved symbols, with the prior's own symbol as the fallback when the ID does not map; a
    symbol reached twice keeps its largest count."""
    src = dict(DE_PRIOR_SOURCES[DE_PRIOR_SOURCE])
    if DE_PRIOR != src["file"]:  # a caller pointed DE_PRIOR at another copy of the same file (primary_rule run dirs)
        src["file"] = DE_PRIOR
    ent, approved = {}, set()
    with open(HGNC_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "Approved":
                approved.add(r["symbol"])
                if r["entrez_id"]:
                    ent[r["entrez_id"]] = r["symbol"]
    out = {}
    with open(src["file"], encoding="utf-8") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            g = ent.get((r.get(src["id"]) or "").strip())
            if not g:
                sym = (r.get(src["symbol"]) or "").strip()
                g = sym if sym in approved else None
            if g:
                c, rk = int(float(r[src["count"]])), float(r[src["rank"]])
                if g not in out or c > out[g][0]:
                    out[g] = (c, rk)
    return out


def de_prior_weights():
    """gene symbol -> count of the DE prior in force (global: n_de_total; crow2019: N_HitLists); the null weight is count + 1"""
    return {g: v[0] for g, v in de_prior_table().items()}


def de_prior_ranks():
    """gene symbol -> normalized DE-prior rank (global: de_prior_norm; crow2019: DE_Prior_Rank; 1 = most often DE)"""
    return {g: v[1] for g, v in de_prior_table().items()}


def expression_levels():
    """gene symbol -> LIBD mean log CPM: DLPFC, else hippocampus. A gene in neither was removed by LIBD's low-expression
    filter (SERPINA3 and IL6 are in the hippocampus table only, at log CPM below 0), so expression_weights puts genes
    with no value in the lowest bin rather than a bin of their own (which gave them a rate of their own and inflated
    their significance, 2026-09-28)."""
    by = {"DLPFC": {}, "HIPPO": {}}
    with open(EXPR_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            if r["region"] in by and r["symbol"]:
                by[r["region"]][r["symbol"]] = max(by[r["region"]].get(r["symbol"], -1e9), float(r["mean_logcpm"]))
    return by["HIPPO"] | by["DLPFC"]


def expression_weights(genes, Lr, M, expr):
    """Per-gene weight for the expression-matched null: genes are put in N_EXPR_BINS equal-count bins of DLPFC expression
    (among the scenario's measured genes that have a value; genes without one join the lowest bin), and each bin's weight
    is its pooled rate of list membership, sum over lists of |L_k & M_k & bin| / sum of |M_k & bin|. A list then draws
    its genes in proportion to that weight, so the null follows how lists actually sample across expression levels.
    Returns (weights, bins)."""
    vals = np.array([expr.get(g, np.nan) for g in genes])
    ok = ~np.isnan(vals)
    bins = np.full(len(genes), -1)
    if ok.any():
        ranks = vals[ok].argsort().argsort()
        bins[ok] = (ranks * N_EXPR_BINS // ok.sum()).astype(int)
    bins[~ok] = 0
    gi = {g: i for i, g in enumerate(genes)}
    hit, tested = collections.Counter(), collections.Counter()
    for k in Lr:
        for g in M[k]:
            tested[bins[gi[g]]] += 1
        for g in Lr[k]:
            hit[bins[gi[g]]] += 1
    rate = {b: (hit[b] / tested[b]) if tested[b] else 0.0 for b in tested}
    w = np.array([rate.get(b, 0.0) for b in bins])
    w = np.where(w > 0, w, np.min(w[w > 0]) if (w > 0).any() else 1.0)
    return w, bins


def prior_probs(total, w):
    """p proportional to w with sum = total, capped at 1 (excess shared among the rest)."""
    p = np.zeros(len(w))
    free = np.ones(len(w), bool)
    rem = total
    while True:
        p[free] = rem * w[free] / w[free].sum()
        over = free & (p > 1)
        if not over.any():
            return p
        p[over] = 1.0
        free &= ~over
        rem = total - p[~free].sum()


# ---------------------------------------------------------------- chance model (study and cohort tier)
def study_tier(sel, M, nk, lists, hitn, expr=None):
    """The study tier for the modelled lists `sel`: each list restricted to its tested set M[k]; P[g, k] = |L_k & M_k| /
    nk[k] for g in M_k (uniform null) and the same total shared out by DE-prior weight (hitn: gene -> DE-prior count, de_prior_weights(); +1,
    missing genes the median); obs = number of lists containing g; exact Poisson-binomial upper tails.
    Returns a dict: Lr, genes, gi, wall, P, PP, meas, obs, pval, tail, pprior, tailP (tail[:, k] = P(X_g >= k))."""
    Lr = {k: set(lists[k]) & M[k] for k in sel}
    genes = sorted(set().union(*(M[k] for k in sel))) if sel else []
    gi = {g: i for i, g in enumerate(genes)}
    wall = np.array([hitn.get(g, np.nan) + 1.0 for g in genes])
    wmed = np.nanmedian(wall)
    wall = np.where(np.isnan(wall), wmed, wall)
    P = np.zeros((len(genes), len(sel)))
    PP = np.zeros((len(genes), len(sel)))
    PE = np.zeros((len(genes), len(sel)))
    wexp, ebins = expression_weights(genes, Lr, M, expr or {})
    meas = np.zeros((len(genes), len(sel)), dtype=bool)
    obs = np.zeros(len(genes), dtype=int)
    for j, k in enumerate(sel):
        p = len(Lr[k]) / nk[k]
        idx = np.array([gi[g] for g in M[k]])
        P[idx, j] = p
        PP[idx, j] = prior_probs(p * len(idx), wall[idx])
        PE[idx, j] = prior_probs(p * len(idx), wexp[idx])
        meas[idx, j] = True
        for g in Lr[k]:
            obs[gi[g]] += 1
    D = pb_matrix(P)
    pval, tail = upper_tail(D, obs)
    DP = pb_matrix(PP)
    pprior, tailP = upper_tail(DP, obs)
    pexpr, tailE = upper_tail(pb_matrix(PE), obs)
    return {"Lr": Lr, "genes": genes, "gi": gi, "wall": wall, "P": P, "PP": PP, "meas": meas, "obs": obs,
            "pval": pval, "tail": tail, "pprior": pprior, "tailP": tailP, "PE": PE, "wexp": wexp, "ebins": ebins,
            "pexpr": pexpr, "tailE": tailE}


def cohort_tier(clusters, genes, gi):
    """The cohort tier: clusters maps cluster -> member lists, each (name, p_g under the uniform null, p_g under the DE
    prior, measured mask, list genes). P(g in cluster) = 1 - prod over members (1 - p_gk); obs = clusters with a member
    list containing g; exact Poisson-binomial over clusters. Returns a dict: cl, PC, PCP, measC, obsC, pC, tailC, pCP, tailCP."""
    cl = sorted(clusters)
    PC = np.zeros((len(genes), len(cl)))
    PCP = np.zeros((len(genes), len(cl)))
    PCE = np.zeros((len(genes), len(cl)))
    measC = np.zeros((len(genes), len(cl)), dtype=bool)
    obsC = np.zeros(len(genes), dtype=int)
    for c_i, c in enumerate(cl):
        mem = clusters[c]
        PC[:, c_i] = 1 - np.prod(1 - np.column_stack([m[1] for m in mem]), axis=1)
        PCP[:, c_i] = 1 - np.prod(1 - np.column_stack([m[2] for m in mem]), axis=1)
        # members built before the expression-matched null carry 5 fields; their expression column is the uniform one
        PCE[:, c_i] = 1 - np.prod(1 - np.column_stack([m[5] if len(m) > 5 else m[1] for m in mem]), axis=1)
        measC[:, c_i] = np.column_stack([m[3] for m in mem]).any(axis=1)
        hit = np.zeros(len(genes), dtype=bool)
        for m in mem:
            hit[[gi[g] for g in m[4]]] = True
        obsC += hit
    DC = pb_matrix(PC)
    pC, tailC = upper_tail(DC, obsC)
    pCP, tailCP = upper_tail(pb_matrix(PCP), obsC)
    pCE, tailCE = upper_tail(pb_matrix(PCE), obsC)
    return {"cl": cl, "PC": PC, "PCP": PCP, "measC": measC, "obsC": obsC, "pC": pC, "tailC": tailC, "pCP": pCP, "tailCP": tailCP,
            "pCE": pCE, "tailCE": tailCE}


def cohort_members(sel, units, st_, M, nk):
    """Cohort-tier members of the study-tier selection `sel` (st_ = study_tier's result): cluster -> list of (name, p
    uniform, p DE prior, measured mask, list genes, p expression-matched), one per study, or for a multi-bank study with
    per-bank lists one per bank list; multi-bank studies without per-bank lists are returned apart (study tier only)."""
    genes, gi, wall, wexp = st_["genes"], st_["gi"], st_["wall"], st_["wexp"]
    P, PP, PE, meas, Lr = st_["P"], st_["PP"], st_["PE"], st_["meas"], st_["Lr"]
    jx = {k: j for j, k in enumerate(sel)}
    clusters, mb_out = {}, []
    for k in sel:
        if k in units["multibank"]:
            bl = units["bank_lists"].get(k)
            if not bl:
                mb_out.append(k)  # stays in the study tier only
                continue
            idx = np.array([gi[g] for g in M[k]])
            for b, names in bl.items():
                Lb = set().union(*(read_list(n) for n in names)) & M[k]
                pb = np.zeros(len(genes))
                pb[idx] = len(Lb) / nk[k]
                ppb = np.zeros(len(genes))
                ppb[idx] = prior_probs(len(Lb) / nk[k] * len(idx), wall[idx])
                peb = np.zeros(len(genes))
                peb[idx] = prior_probs(len(Lb) / nk[k] * len(idx), wexp[idx])
                clusters.setdefault(bank_cluster(units, b), []).append((f"{k}:{b}", pb, ppb, meas[:, jx[k]], Lb, peb))
            continue
        j = jx[k]
        clusters.setdefault(cluster_of(units, k), []).append((k, P[:, j], PP[:, j], meas[:, j], Lr[k], PE[:, j]))
    return clusters, mb_out


# ---------------------------------------------------------------- main
def main():
    pl = json.load(open(os.path.join(HERE, "primary_lists.json")))
    uni = json.load(open(os.path.join(HERE, "universes.json")))["lists"]
    units = load_units(pl)
    studies, lanes = units["studies"], units["lanes"]
    lists, M, nk, utype, nouni = load_lists(units, uni)
    hitn = de_prior_weights()
    expr = expression_levels()
    today = datetime.date.today().isoformat()
    n_groups = len(units["groups"])
    head = [f"# script: recurrence.py, run {today}",
            f"# input: primary_lists.json sha256 {sha(os.path.join(HERE, 'primary_lists.json'))} (rule A: DE and DE-adjacent lists, Paul 2026-09-25)",
            f"# input: universes.json sha256 {sha(os.path.join(HERE, 'universes.json'))} (build_universes.py; gene sets in universe_sets/)",
            f"# input: ../corpus/lineage.json sha256 {sha(LINEAGE)} (lanes)",
            f"# input: {os.path.relpath(DE_PRIOR, HERE)} sha256 {sha(DE_PRIOR)} ({DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['label']}; DE-prior null only)",
            f"# {len(units['all'])} in-corpus studies with a primary list; left out: {len(units['meta'])} meta-lane, "
            f"{len(units['multi'])} multi-series reanalyses ({', '.join(units['all'][k]['label'] for k in units['multi'])}), "
            f"{len(units['reanalysis'])} other reanalyses (lineage group 3, no new data; each dataset enters through the study that generated it); "
            f"{len(studies)} studies in {n_groups} same_data_as groups; representative chosen per scenario (recurrence_units.tsv)",
            f"# list modality: only DE lists represent a study; {len(units['modality_dropped'])} studies had a non-DE primary list dropped"
            + (f" ({', '.join(units['labels'][k] + ': ' + ', '.join(v) for k, v in sorted(units['modality_dropped'].items()))})" if units["modality_dropped"] else "")
            + f", {len(units['modality_only'])} of them left out for want of a DE list",
            f"# list role: only first_look lists represent a study; {len(units['role_dropped'])} studies had a primary list with role "
            f"reanalysis dropped" + (f" ({', '.join(units['labels'][k] + ': ' + ', '.join(v) for k, v in sorted(units['role_dropped'].items()))})"
                                     if units["role_dropped"] else "") + f", {len(units['role_only'])} of them left out for want of a first_look list; "
            f"{len(units['role_alternates'])} studies hold a reanalysis-role list as an alternate only"
            + (f" ({', '.join(units['labels'][k] + ': ' + ', '.join(v) for k, v in sorted(units['role_alternates'].items()))})"
               if units["role_alternates"] else ""),
            f"# multi-bank studies (lineage 'banks'): {len(units['multibank'])}"
            + (" (" + "; ".join(f"{studies[k]['label']} {'+'.join(b)}, " + ("per-bank lists" if k in units["bank_lists"] else
                                "no per-bank lists: study tier only, left out of the cohort tier") for k, b in sorted(units["multibank"].items())) + ")"
               if units["multibank"] else "") + "; their lane does not name or merge a cohort cluster",
            "# model: list k draws gene g with p_gk = |L_k & M_k| / n_k if g in M_k, else 0 (M_k tested gene set, n_k = stated genes "
            "tested if <= |M_k|, else |M_k|); every tested gene equally likely; lists restricted to M_k; count of lists containing g is "
            "Poisson-binomial; expected = sum over genes of P(X_g >= k)",
            f"# exp_expr_*: expression-matched null, same list totals shared out in proportion to the pooled list-membership rate of "
            f"the gene's expression decile (LIBD BrainSeq DLPFC mean log CPM, else hippocampus, {os.path.relpath(EXPR_FILE, HERE)}; "
            "genes in neither, removed by LIBD's low-expression filter, join the lowest decile); columns exp_expr_*, genes_q10_expr, p_expr, q_expr_bh, expr_decile, p_cluster_expr, q_cluster_expr_bh",
            f"# exp_prior_*: DE-prior null, same list totals shared out in proportion to (count + 1), count = {DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['count']} of the {DE_PRIOR_SOURCE} prior; genes not in "
            "the prior get the median weight; " + DE_PRIOR_GENE_NOTE,
            f"# u18_*: the same modelled lists, unrestricted, each a draw of its own size from U = {U} genes (universe effect alone)",
            f"# old_*: every representative list in the scenario, including lists the universe model leaves out, unrestricted, U = {U}; "
            "not a reproduction of an earlier run",
            "# cohort tier: a group takes the named lanes (pitt, stanley, libd, cmc, mssm, harvard, aus, pritzker, ipsc) of ALL its members; "
            "cmc + mssm + pitt + libd one cluster (CommonMind is built from Mount Sinai, Pittsburgh and Penn brains; HBCC donors (cmc lane) "
            "recur in LIBD BrainSeq, donor_overlap/); stanley one cluster; groups spanning several named clusters merge them"
            + (f" (merged: {'; '.join('+'.join(m) for m in units['merged'])})" if units["merged"] else " (none merged)")
            + "; a group with no named lane is its own cohort",
            "# cohort chance: P(g in cluster c) = 1 - prod over member lists k of (1 - p_gk); Poisson-binomial over clusters"]
    if switches():  # forking-paths switches away from their defaults (module docstring)
        head.append("# switches (not the defaults): " + "; ".join(switches())
                    + (f"; computed-FDR lists: {', '.join(n for v in units['computed_fdr'].values() for n in v)}" if units["computed_fdr"] else "")
                    + (f"; added as independent studies: {', '.join(units['added'])}" if units["added"] else ""))

    res, genes_out, unit_rows = [], [], []
    for name, f in SCENARIOS:
        reps = representatives(units, f, usable(units, lists))
        sel_old = sorted(reps, key=lambda k: (0, int(k), "") if k.isdigit() else (1, 0, k))
        for k in sel_old:
            unit_rows.append({"scenario": name, "representative": k, "label": studies[k]["label"],
                              "modelled": "yes" if k in M else f"no ({nouni.get(k) or utype[k]})",
                              "cluster": cohort_label(units, k), "lane": lanes.get(k) or "",
                              "group_members": ";".join(f"{studies[m]['label']} ({m}, {lanes.get(m) or 'none'})" for m in reps[k])})
        count_old = {}
        for k in sel_old:
            for g in lists[k]:
                count_old[g] = count_old.get(g, 0) + 1
        dist = poisson_binomial([min(len(lists[k]) / U, 1.0) for k in sel_old])
        sel = [k for k in sel_old if k in M]
        count18 = {}
        for k in sel:
            for g in lists[k]:
                count18[g] = count18.get(g, 0) + 1
        dist18 = poisson_binomial([min(len(lists[k]) / U, 1.0) for k in sel])
        st_ = study_tier(sel, M, nk, lists, hitn, expr)
        Lr, genes, gi = st_["Lr"], st_["genes"], st_["gi"]
        P, meas, obs = st_["P"], st_["meas"], st_["obs"]
        pval, tail, pprior, tailP = st_["pval"], st_["tail"], st_["pprior"], st_["tailP"]
        ebins, pexpr, tailE = st_["ebins"], st_["pexpr"], st_["tailE"]
        row = {"scenario": name, "tier": "study", "n_lists_old": len(sel_old), "n_lists": len(sel),
               "n_lists_dropped": len(sel_old) - len(sel),
               "median_size": int(np.median([len(Lr[k]) for k in sel])) if sel else 0,
               "genes_outside_universe": sum(len(set(lists[k]) - M[k]) for k in sel)}
        for kk in KS:
            row[f"obs_ge{kk}"] = int((obs >= kk).sum())
            row[f"exp_ge{kk}"] = round(float(tail[:, kk].sum()), 1) if kk < tail.shape[1] else 0.0
            row[f"exp_prior_ge{kk}"] = round(float(tailP[:, kk].sum()), 1) if kk < tailP.shape[1] else 0.0
            row[f"exp_expr_ge{kk}"] = round(float(tailE[:, kk].sum()), 1) if kk < tailE.shape[1] else 0.0
            row[f"u18_obs_ge{kk}"] = sum(1 for c in count18.values() if c >= kk)
            row[f"u18_exp_ge{kk}"] = round(U * sum(dist18[kk:]), 1)
            row[f"old_obs_ge{kk}"] = sum(1 for c in count_old.values() if c >= kk)
            row[f"old_exp_ge{kk}"] = round(U * sum(dist[kk:]), 1)
        q = bh(pval)
        qP = bh(pprior)
        qE = bh(pexpr)
        row["genes_q10"] = int((q < 0.1).sum())
        row["genes_q10_prior"] = int((qP < 0.1).sum())
        row["genes_q10_expr"] = int((qE < 0.1).sum())
        res.append(row)
        if name not in GENE_SCENARIOS:
            continue
        clusters, mb_out = cohort_members(sel, units, st_, M, nk)
        ct_ = cohort_tier(clusters, genes, gi)
        cl, PC, measC, obsC = ct_["cl"], ct_["PC"], ct_["measC"], ct_["obsC"]
        pC, tailC, pCP, tailCP = ct_["pC"], ct_["tailC"], ct_["pCP"], ct_["tailCP"]
        pCE, tailCE = ct_["pCE"], ct_["tailCE"]
        qC, qCP, qCE = bh(pC), bh(pCP), bh(pCE)
        crow = {"scenario": name, "tier": "cohort", "n_lists_old": "", "n_lists": sum(len(v) for v in clusters.values()),
                "n_lists_dropped": "", "median_size": "", "genes_outside_universe": "", "n_clusters": len(cl)}
        for kk in KS:
            inc = kk in KS_COHORT
            crow[f"obs_ge{kk}"] = int((obsC >= kk).sum()) if inc else ""
            crow[f"exp_ge{kk}"] = (round(float(tailC[:, kk].sum()), 1) if kk < tailC.shape[1] else 0.0) if inc else ""
            crow[f"exp_prior_ge{kk}"] = (round(float(tailCP[:, kk].sum()), 1) if kk < tailCP.shape[1] else 0.0) if inc else ""
            crow[f"exp_expr_ge{kk}"] = (round(float(tailCE[:, kk].sum()), 1) if kk < tailCE.shape[1] else 0.0) if inc else ""
            for pre in ("u18_obs", "u18_exp", "old_obs", "old_exp"):
                crow[f"{pre}_ge{kk}"] = ""
        crow["genes_q10"] = int((qC < 0.1).sum())
        crow["genes_q10_prior"] = int((qCP < 0.1).sum())
        crow["genes_q10_expr"] = int((qCE < 0.1).sum())
        crow["clusters"] = "; ".join(f"{c} ({len(clusters[c])})" for c in cl)
        crow["multibank_left_out"] = "; ".join(f"{studies[k]['label']} ({k})" for k in mb_out)
        res.append(crow)
        expd = P.sum(axis=1)
        expC = PC.sum(axis=1)
        for i in np.argsort(pval, kind="stable"):
            if obs[i] < 1:
                continue
            g = genes[i]
            ks = [k for k in sel if g in Lr[k]]
            ds = [lists[k][g] for k in ks]
            up, down, both = ds.count("up"), ds.count("down"), ds.count("both")
            st = sign_test(up, down) if q[i] < 0.1 else None
            genes_out.append({
                "scenario": name, "gene": g, "n_lists_measured": int(meas[i].sum()), "n_lists_observed": int(obs[i]),
                "expected": round(float(expd[i]), 3), "p": f"{pval[i]:.3g}", "q_bh": f"{q[i]:.3g}",
                "p_prior": f"{pprior[i]:.3g}", "q_prior_bh": f"{qP[i]:.3g}",
                "p_expr": f"{pexpr[i]:.3g}", "q_expr_bh": f"{qE[i]:.3g}", "expr_decile": int(ebins[i]) + 1,
                "de_prior_hitlists": hitn.get(g, ""),
                "up": up, "down": down, "both": both, "other_direction": len(ds) - up - down - both,
                "sign_test_p": "" if st is None else f"{st:.3g}",
                "n_clusters_measured": int(measC[i].sum()), "n_clusters_observed": int(obsC[i]),
                "expected_clusters": round(float(expC[i]), 3), "p_cluster": f"{pC[i]:.3g}", "q_cluster_bh": f"{qC[i]:.3g}",
                "q_cluster_prior_bh": f"{qCP[i]:.3g}", "p_cluster_expr": f"{pCE[i]:.3g}", "q_cluster_expr_bh": f"{qCE[i]:.3g}",
                "studies": ";".join(studies[k]["label"] for k in ks)})
    fields = list(dict.fromkeys(k for r in res for k in r))
    with open(os.path.join(_out(), "recurrence.tsv"), "w") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(res)
    with open(os.path.join(_out(), "recurrence_genes.tsv"), "w") as fh:
        fh.write("\n".join(head + [
            f"# scenarios: {'; '.join(GENE_SCENARIOS)}; genes in at least one list, sorted by p within scenario",
            "# p = exact Poisson-binomial P(X_g >= observed); q_bh = Benjamini-Hochberg over every gene measured in the scenario "
            "(including genes in no list); p_prior, q_prior_bh: the same under the DE-prior null (not evidence for any gene: see recurrence.tsv's exp_prior_* line); de_prior_hitlists: the gene's "
            f"{DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['count']} of {DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['of']} in the {DE_PRIOR_SOURCE} DE prior (blank if not in the prior); sign_test_p = two-sided exact binomial test of up vs "
            "down, for genes with q_bh < 0.1; 'both' = up and down in one study, left out of the sign test",
            "# p_cluster, q_cluster_bh, q_cluster_prior_bh: the same for the cohort tier"]) + "\n")
        w = csv.DictWriter(fh, list(genes_out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(genes_out)
    left = ([{"scenario": "(left out)", "representative": k, "label": units["all"][k]["label"], "modelled": "no (meta lane)",
              "cluster": "", "lane": lanes.get(k) or "", "group_members": ""} for k in units["meta"]]
            + [{"scenario": "(left out)", "representative": k, "label": units["all"][k]["label"],
                "modelled": "no (multi-series reanalysis: " + ", ".join(sorted(set().union(*partner_series(units["all"][k])))) + ")",
                "cluster": "", "lane": lanes.get(k) or "", "group_members": ""} for k in units["multi"]]
            + [{"scenario": "(left out)", "representative": k, "label": units["all"][k]["label"],
                "modelled": "no (reanalysis, lineage group 3)", "cluster": "", "lane": lanes.get(k) or "", "group_members": ""}
               for k in units["reanalysis"]]
            + [{"scenario": "(left out)", "representative": k, "label": units["labels"][k],
                "modelled": "no (primary list role reanalysis: " + ", ".join(units["role_dropped"][k]) + ")", "cluster": "",
                "lane": lanes.get(k) or "", "group_members": ""} for k in units["role_only"]]
            + [{"scenario": "(left out)", "representative": k, "label": units["labels"][k],
                "modelled": "no (primary list not DE: " + ", ".join(units["modality_dropped"][k]) + ")", "cluster": "",
                "lane": lanes.get(k) or "", "group_members": ""} for k in units["modality_only"]])
    with open(os.path.join(_out(), "recurrence_units.tsv"), "w") as fh:
        fh.write("\n".join(head + ["# one row per scenario representative (group_members = its whole same_data_as group, "
                                   "pmid and lane), then the studies left out before grouping"]) + "\n")
        w = csv.DictWriter(fh, list(unit_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(unit_rows + left)
    for r in res:
        print({k: v for k, v in r.items() if k != "clusters"})
    write_report()


# ---------------------------------------------------------------- report (from the tsv files only)
def read_tsv(fn):
    lines = open(os.path.join(_out(), fn), encoding="utf-8").read().splitlines()
    head = [l for l in lines if l.startswith("#")]
    rows = list(csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t"))
    return head, rows


def n_(x):
    """thousands separator for integer-like strings, else as is."""
    try:
        f = float(x)
    except (TypeError, ValueError):
        return x
    if f == int(f) and "." not in str(x):
        return f"{int(f):,}"
    return f"{f:,.1f}"


def write_report():
    head, res = read_tsv("recurrence.tsv")
    _, genes = read_tsv("recurrence_genes.tsv")
    _, unit_rows = read_tsv("recurrence_units.tsv")
    uni = json.load(open(os.path.join(HERE, "universes.json")))
    gshas = {fn: sha(os.path.join(_out(), fn)) for fn in ("recurrence.tsv", "recurrence_genes.tsv", "recurrence_units.tsv")}
    run = head[0].split("run ")[-1]
    inputs = [h[2:] for h in head if h.startswith("# input:")]
    info = [h[2:] for h in head if "in-corpus studies" in h][0]
    role_line = [h[2:] for h in head if h.startswith("# list role:")][0]
    mb_line = [h[2:] for h in head if h.startswith("# multi-bank studies")][0]
    L = ["---", "type: detail", "bib: My Library", f"modified: {run}", "---", "", "> [!note] Paul", ">", "",
         "# Recurrence against per-list gene universes", "",
         f"Generated by `recurrence.py` from `recurrence.tsv` (sha256 {gshas['recurrence.tsv']}), `recurrence_genes.tsv` "
         f"({gshas['recurrence_genes.tsv']}) and `recurrence_units.tsv` ({gshas['recurrence_units.tsv']}), run {run}. "
         "Edit the script, not this file.", "", "Inputs:", ""]
    L += [f"- {x}" for x in inputs]
    L += ["", "```", "python3 build_universes.py   # universes.json, universe_sets/*.txt",
          "python3 recurrence.py        # recurrence.tsv, recurrence_genes.tsv, recurrence_units.tsv, this report", "```", "",
          "## Model", "",
          "For list k with tested gene set M_k and list L_k (HGNC symbols), the list is restricted to M_k and gene g is drawn with "
          "p_gk = |L_k ∩ M_k| / n_k if g is in M_k, else 0. n_k is the paper's stated number of genes tested when it is stated in "
          "genes and is not larger than |M_k|; otherwise n_k = |M_k|. The number of lists containing g is Poisson-binomial over k; "
          "the expected number of genes at >= k lists is the sum over genes of P(X_g >= k). Per-gene p is the exact P(X_g >= "
          "observed); q is Benjamini-Hochberg over every gene measured in the scenario, including genes in no list.", "",
          "**Uniform-prior assumption.** The model treats every tested gene in a list's universe as equally likely to be in that "
          "list. Genes differ in how often they are differentially expressed in any comparison (the DE prior; Crow et al. 2019, "
          "[@crowPredictabilityHumanDifferential2019]). The p-values therefore measure recurrence above uniform chance, not "
          "specificity to SZ. For the aggregate expectations, the DE-prior null keeps each list's total and shares it out in proportion "
          + (f"to n_de_total + 1, where n_de_total is the number of the 8,712 Gemma contrasts (all technologies, 2-fold) in which "
             "the gene is differentially expressed (`de_prior_lists/our_prior/1_global_our_prior.tsv`, our own build; default since "
             "2026-09-30; the external 4-fold global prior of 5,785 contrasts remains available as DE_PRIOR_SOURCE=global, and the "
             "Crow et al. 2019 prior as DE_PRIOR_SOURCE=crow2019); "
             if DE_PRIOR_SOURCE == "ours" else
             f"to n_de_total + 1, where n_de_total is the number of the 5,785 Gemma contrasts (all technologies) in which the gene "
             "is differentially expressed (`de_prior_lists/de_prior_ranks/1_global_de_prior.tsv`, the external 4-fold global prior, "
             "DE_PRIOR_SOURCE=global); "
             if DE_PRIOR_SOURCE == "global" else
             "to N_HitLists + 1, where N_HitLists is the number of the paper's 635 DE hit lists that contain the gene "
             "(`reference/DE_Prior.txt`, from github.com/maggiecrow/DEprior, provenance in `reference/DE_Prior.meta.json`); ")
          + "genes not in the prior get the median weight. Its columns are `exp_prior_*`, `p_prior`, `q_prior_bh` and "
          "`q_cluster_prior_bh`. The prior was built from a compendium of many diseases and tissues, so it is one estimate of "
          "gene-level DE propensity, not a brain-specific one. **Use.** The DE-prior null gives the expected counts of genes in "
          "at least N lists (the tables below); per-gene significance is reported under the uniform and expression-matched "
          "nulls. The per-gene DE-prior columns stay in `recurrence_genes.tsv` and are not evidence for any gene (Paul, "
          "2026-09-29): on Gemma's bipolar disorder and major depression postmortem brain lists (FDR 0.05), the global DE-prior "
          "null gives 350 genes at q < 0.1 against 1 under the uniform null, 349 of them low prior "
          "(`de_prior_lists/checks_2026-09-29/checks_report.md`).", "",
          "**Expression-matched null** (Paul, 2026-09-28: \"a uniform random is not right - it should follow the expression level "
          "distributions\"). Each list keeps its total, shared out in proportion to the pooled list-membership rate of the gene's "
          "expression decile: genes are binned in ten equal bins of LIBD BrainSeq DLPFC mean log CPM (hippocampus where DLPFC has "
          "none; genes in neither, removed by LIBD's low-expression filter, join the lowest bin), and a bin's rate is the number of "
          "list memberships of its genes over the number of times they were measured, over all lists of the scenario. The null so "
          "follows how the lists actually draw genes across expression levels. The rates are estimated from the same lists, so a "
          "decile's recurrent genes contribute to their own expectation; with about 1,800 genes per decile the effect is small. "
          "Columns `exp_expr_*`, `genes_q10_expr`, `p_expr`, `q_expr_bh`, `expr_decile`, `p_cluster_expr`, `q_cluster_expr_bh`.", "",
          "**Units.** " + info + ". A study is left out as a multi-series reanalysis when its same_data_as evidence names GEO "
          "series it shares with different partners (two partners whose shared series are disjoint); a pair sharing two series at "
          "once (a SuperSeries and its SubSeries) is one dataset. No same_data_as link runs through a left-out study. The "
          "representative of a group is chosen per scenario among members that pass the scenario's filter: published over "
          "reproduced over Gemma, then FDR lists, then earliest PMID. A blood study cannot represent a brain scenario.", "",
          "**List role.** Only first_look lists (the study's own new data) can be a study's primary list; a primary list whose "
          "manifest role is reanalysis (the study's analysis of data first published elsewhere) is dropped (Paul, 2026-09-24: "
          "\"keep 'first look at the data' separate\"). This run: " + role_line.split("; ", 1)[1] + ".", "",
          "**Cohort tier.** Each same_data_as group takes the named lineage lanes of all its members, not only of its "
          "representative. cmc, mssm, pitt and libd form one cluster (CommonMind is built from Mount Sinai, Pittsburgh and Penn "
          "brains; HBCC donors, lane cmc, recur in LIBD BrainSeq, `donor_overlap/`). Stanley is one cluster: the Array Collection "
          "and Neuropathology Consortium are different subjects, but many studies do not say which, so this is conservative. A "
          "group with no named lane is its own cohort. A gene counts once per cluster if any member list has it; P(g in cluster c) "
          "= 1 - prod over member lists of (1 - p_gk); Poisson-binomial over clusters. Clusters merged because one group spans "
          "them: " + (re.search(r"\(merged: ([^)]*)\)", [h for h in head if "cohort tier:" in h][0]) or re.match("(none)", "none")).group(1) + ".",
          "", "**Multi-bank studies.** A study whose lineage node lists several banks (`banks`) does not name or merge a cohort "
          "cluster, so one bridging study cannot merge two banks. In the cohort tier it contributes its per-bank lists (manifest "
          "`bank`) to the matching bank clusters; without per-bank lists it stays in the study tier and is left out of the cohort "
          "tier. Multi-bank studies in this run: " + mb_line[len("multi-bank studies (lineage 'banks'): "):].split("; their lane")[0] + ".",
          "", "**Direction.** A gene in both an up and a down list of one study has direction 'both' and is left out of the sign "
          "test.", ""]
    # left out
    lo = [u for u in unit_rows if u["scenario"] == "(left out)"]
    ms = [u for u in lo if u["modelled"].startswith("no (multi")]
    ra = [u for u in lo if u["modelled"].startswith("no (reanalysis")]
    rr = [u for u in lo if u["modelled"].startswith("no (primary list role")]
    # studies whose primary lists are all non-DE (e.g. proteomics) are their own group, not meta-lane (code audit CA13, 2026-09-26)
    nd = [u for u in lo if u["modelled"].startswith("no (primary list not DE")]
    L += ["## Studies left out before grouping", "",
          f"{len(lo) - len(ms) - len(ra) - len(rr) - len(nd)} meta-lane studies, {len(ms)} multi-series reanalyses, {len(ra)} other reanalyses "
          "(lineage group 3: no new expression data, so each dataset enters once, through the study that generated it), "
          f"{len(nd)} studies whose primary lists are all non-DE (proteomic or other modality; list modality above) and "
          f"{len(rr)} studies whose only primary list has role reanalysis" + (f" ({', '.join(u['label'] for u in rr)})" if rr else "") + ".", "",
          "Other reanalyses: " + ", ".join(sorted(u["label"] for u in ra)) + ".", "",
          "| Study | PMID | Lane | GEO series in its same_data_as evidence |", "|---|---|---|---|"]
    L += [f"| {u['label']} | {u['representative']} | {u['lane']} | {u['modelled'][len('no (multi-series reanalysis: '):-1]} |" for u in ms]
    L += ["", "Meta-lane: " + ", ".join(sorted(u["label"] for u in lo if u not in ms and u not in ra and u not in rr and u not in nd)) + ".", "",
          "Primary lists all non-DE: " + ", ".join(sorted(u["label"] for u in nd)) + ".", ""]
    # universe types
    lists = uni["lists"]
    tc = {}
    for n, r in lists.items():
        t = tc.setdefault(r["type"], [0, set()])
        t[0] += 1
        t[1].add(r["study"])
    L += ["## Universe types", "",
          f"{len({r['study'] for r in lists.values()})} studies, {len(lists)} primary lists (`universes.json`, "
          f"primary_lists.json sha256 {uni['inputs']['primary_lists.json'][:12]}). rnaseq_default: {uni['rnaseq_default']['n_genes']:,} genes, "
          "union of " + "; ".join(f"{p['label']} ({p['n_genes']:,})" for p in uni["rnaseq_default"]["sources"]) + ".", "",
          "| Type | Lists | Studies |", "|---|---|---|"]
    L += [f"| {t} | {v[0]} | {len(v[1])} |" for t, v in sorted(tc.items(), key=lambda x: -x[1][0])]
    reasons = {}
    for n, r in lists.items():
        if not r.get("gene_set"):
            reasons.setdefault(r.get("reason") or "gene set not known", []).append(r["label"])
    L += ["", f"Left out of the universe model ({sum(len(v) for v in reasons.values())} lists):", "", "| Reason | Lists |", "|---|---|"]
    L += [f"| {k} | {len(v)}: {', '.join(sorted(v))} |" for k, v in sorted(reasons.items(), key=lambda x: -len(x[1]))]
    small = sorted(((r["label"], r["n_membership"] / r["n_k"]) for r in lists.values()
                    if r.get("gene_set") and r["n_k"] < r["n_membership"]), key=lambda x: -x[1])
    L += ["", f"{len(small)} lists use a stated n_tested smaller than |M_k| as n_k; each gene in M_k then gets |L_k| / n_k, so the "
          f"list's expected size is |L_k| x |M_k| / n_k, a factor of {small[-1][1]:.2f} ({small[-1][0]}) to {small[0][1]:.2f} "
          f"({small[0][0]}). This biases towards less recurrence." if small else "", ""]
    # study tier table
    st = [r for r in res if r["tier"] == "study"]
    L += ["## Study tier", "",
          "Each cell: observed / expected (uniform) / expected (DE prior) / expected (expression-matched). Lists: modelled lists (restricted to their universes) "
          "and, in brackets, all representative lists in the scenario.", "",
          "| Scenario | Lists | ≥ 2 | ≥ 3 | ≥ 5 | ≥ 10 | Genes q < 0.1 (uniform / expression) |", "|---|---|---|---|---|---|---|"]
    for r in st:
        L.append(f"| {r['scenario']} | {r['n_lists']} ({r['n_lists_old']}) | " + " | ".join(
            f"{n_(r[f'obs_ge{k}'])} / {n_(r[f'exp_ge{k}'])} / {n_(r[f'exp_prior_ge{k}'])} / {n_(r.get(f'exp_expr_ge{k}', ''))}" for k in KS)
            + f" | {r['genes_q10']} / {r.get('genes_q10_expr', '')} |")
    L += ["", "**Universe model alone** (`u18_*`): the same modelled lists, unrestricted, each a random draw of its own size from "
          "U = 18,000 genes; observed / expected. **All lists** (`old_*`): every representative list in the scenario, including "
          "the lists the universe model leaves out, unrestricted, U = 18,000. This is not a reproduction of any earlier run: it uses "
          "the current exclusions and collapsing.", "",
          "| Scenario | ≥ 3, U = 18,000, modelled lists | ≥ 5, same | ≥ 3, all lists | ≥ 5, all lists |", "|---|---|---|---|---|"]
    for r in st:
        L.append(f"| {r['scenario']} | {n_(r['u18_obs_ge3'])} / {n_(r['u18_exp_ge3'])} | {n_(r['u18_obs_ge5'])} / {n_(r['u18_exp_ge5'])} | "
                 f"{n_(r['old_obs_ge3'])} / {n_(r['old_exp_ge3'])} | {n_(r['old_obs_ge5'])} / {n_(r['old_exp_ge5'])} |")
    # cohort
    ch = [r for r in res if r["tier"] == "cohort"]
    L += ["", "## Cohort tier", "", "Observed / expected (uniform) / expected (DE prior) / expected (expression-matched).", "",
          "| Scenario | Clusters (lists) | ≥ 2 | ≥ 3 | ≥ 5 | Genes q < 0.1 (uniform / expression) |", "|---|---|---|---|---|---|"]
    for r in ch:
        L.append(f"| {r['scenario']} | {r['n_clusters']} ({r['n_lists']}) | " + " | ".join(
            f"{n_(r[f'obs_ge{k}'])} / {n_(r[f'exp_ge{k}'])} / {n_(r[f'exp_prior_ge{k}'])} / {n_(r.get(f'exp_expr_ge{k}', ''))}" for k in KS_COHORT)
            + f" | {r['genes_q10']} / {r.get('genes_q10_expr', '')} |")
    L.append("")
    for r in ch:
        L.append(f"Clusters, {r['scenario']}: {r['clusters']}."
                 + (f" Left out of this tier (multi-bank, no per-bank lists): {r['multibank_left_out']}." if r.get("multibank_left_out") else ""))
        L.append("")
    # genes
    for scen in GENE_SCENARIOS:
        g = [x for x in genes if x["scenario"] == scen]
        qs = [x for x in g if float(x["q_bh"]) < 0.1]
        qc = [x for x in g if float(x["q_cluster_bh"]) < 0.1]
        qe = [x for x in g if x.get("q_expr_bh") and float(x["q_expr_bh"]) < 0.1]
        qce = [x for x in g if x.get("q_cluster_expr_bh") and float(x["q_cluster_expr_bh"]) < 0.1]
        L += [f"## Genes, {scen}", "",
              f"Study tier under the uniform null, q < 0.1: {len(qs)}" + (f" ({', '.join(x['gene'] for x in qs)})" if qs else "") + ".", "",
              f"Study tier under the expression-matched null, q < 0.1: {len(qe)}" + (f" ({', '.join(x['gene'] + ' (decile ' + str(x['expr_decile']) + ')' for x in qe)})" if qe else "") + ".", "",
              f"Cohort tier under the uniform null, q < 0.1: {len(qc)}" + (f" ({', '.join(x['gene'] + ' ' + x['n_clusters_observed'] + ' clusters, q ' + x['q_cluster_bh'] for x in qc)})" if qc else "") + ".", "",
              f"Cohort tier under the expression-matched null, q < 0.1: {len(qce)}" + (f" ({', '.join(x['gene'] + ' ' + x['n_clusters_observed'] + ' clusters, q ' + x['q_cluster_expr_bh'] for x in qce)})" if qce else "") + ".", "",
              "Per-gene q under the DE-prior null is in `recurrence_genes.tsv` (q_prior_bh, q_cluster_prior_bh) and is not evidence for any gene (Model, Use).", ""]
        one = [x for x in qs if (int(x["up"]) == 0 or int(x["down"]) == 0) and int(x["both"]) == 0 and int(x["other_direction"]) == 0]
        sig = [x for x in qs if x["sign_test_p"] and float(x["sign_test_p"]) < 0.05]
        L += [f"Direction, study-tier genes at q < 0.1: {len(one)} with every list in one direction, {len(qs) - len(one)} otherwise; "
              f"sign test below 0.05 for {len(sig)}" + (f" ({', '.join(x['gene'] + ' ' + x['up'] + '/' + x['down'] for x in sig)})" if sig else "")
              + ". With 5 directional lists or fewer, a unanimous gene cannot reach a sign-test p below 0.0625.", ""]
        if scen == "published, brain":
            L += [f"Top 30 genes (study tier). Measured: lists whose universe contains the gene. DE prior: the gene's {DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['count']} of {DE_PRIOR_SOURCES[DE_PRIOR_SOURCE]['of']}.",
                  "", "| Gene | Measured | Observed | Expected | p | q | q, expression | DE prior count | Up / down (both) | Sign test p | Clusters obs / measured | q cohort |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"]
            for x in g[:30]:
                L.append(f"| {x['gene']} | {x['n_lists_measured']} | {x['n_lists_observed']} | {x['expected']} | {x['p']} | {x['q_bh']} | "
                         f"{x.get('q_expr_bh', '')} | {x['de_prior_hitlists']} | {x['up']} / {x['down']}" + (f" ({x['both']})" if x["both"] != "0" else "")
                         + (f" (+{x['other_direction']} no direction)" if x["other_direction"] != "0" else "")
                         + f" | {x['sign_test_p']} | {x['n_clusters_observed']} / {x['n_clusters_measured']} | {x['q_cluster_bh']} |")
            L.append("")
    # caveats
    outside = sorted(((r["label"], r["list_genes_outside"], r["list_genes"], r["source"]) for r in lists.values()
                      if r.get("gene_set") and r.get("list_genes_outside")), key=lambda x: -x[1])
    L += ["## Caveats that can be measured", "",
          f"- **List genes outside their universe**: {sum(x[1] for x in outside):,} of "
          f"{sum(r['list_genes'] for r in lists.values() if r.get('gene_set')):,} list genes across all lists with a gene set "
          "(including meta-lane lists), dropped by the restriction. Study-tier scenarios: " + "; ".join(f"{r['scenario']} {r['genes_outside_universe']}" for r in st)
          + ". Largest: " + "; ".join(f"{a} {b} of {c} ({s})" for a, b, c, s in outside[:6]) + ".",
          "- **rnaseq_default is brain-derived.** " + str(tc.get("rnaseq_default", [0])[0]) + " lists use it, including iPSC, blood "
          "and fibroblast RNA-seq; genes expressed only outside brain are missing from it.",
          "- **Imputed expression.** TWAS lists (Huckins 2019, Gusev 2018) use rnaseq_default because their predictor gene sets are "
          "not held.",
          "- **Single-nucleus lists do not enter.** No study in the sc lineage lane has a primary list; the cell-type-stacked lists "
          "are type exclude (table above) and enter only the all-lists columns.",
          "- **Cohort independence.** P(g in cluster) assumes the member lists of a cluster are independent draws, which "
          "overstates chance for clusters with shared donors and so is conservative.",
          "- **Per-list evidence** for each universe assignment is in `universes.json` (`evidence`, `n_tested_source`); every "
          "representative and its group are in `recurrence_units.tsv`.", ""]
    with open(os.path.join(_out(), "recurrence_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    if "--report" in sys.argv:
        write_report()
    else:
        main()
