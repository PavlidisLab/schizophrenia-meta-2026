#!/usr/bin/env python3
"""H1: pairwise agreement between studies' primary SZ vs control lists, old era vs new era.

Tracker, Step 3, H1: within-era pairwise overlap (hypergeometric against the tested universe)
and directional concordance, old-era pairs vs new-era pairs, with list size accounted for;
correlation of test statistics where full summary statistics exist (Hoffman et al.); pairs that
share samples analysed separately.

Units and chance model come from recurrence.py and universes.json:
  - one proposed primary list per study (primary_lists.json; up and down lists of one contrast
    merged); meta-lane studies and multi-series reanalyses left out, same_data_as groups
    collapsed to one representative per scenario (recurrence.load_units / representatives), so
    pairs that share data never enter;
  - list k is restricted to its tested gene set M_k; a pair is scored over M_A & M_B;
  - cohort cluster = recurrence.cluster_of (named lanes of all members of the study's group;
    cmc + mssm + pitt + libd merged); a pair is 'same' cohort when both fall in one named cluster.
    A multi-bank study (lineage 'banks') takes the clusters of its per-bank lists; with none, its pairs
    have cohort 'multi-bank' and enter only the all-cohort summaries and tests;
  - only first_look lists represent a study (recurrence.load_units drops reanalysis-role primary lists).

Per pair (A = the earlier study, by year then PMID), over U = M_A & M_B:
  overlap, expected |a||b|/|U|, one-sided hypergeometric p, Jaccard, fold = observed / expected;
  top-N overlap: N = min(|a|, |b|, 100), each list's N genes with the smallest p (or FDR when the
    list has no p), when at least 80% of the list's genes carry one; expected N^2/|U|;
  concordance: fraction of shared genes with a direction in both lists that agree; two-sided exact
    binomial p (recurrence.sign_test);
  full statistics (build_universes.fullstats_rows, the same tables the fullstats universes come
    from): rho = logFC Spearman correlation (signed effect; no table carries a t statistic) over
    genes in both tables; pi1 = 1 - qvalue(p)$pi0 (R qvalue, defaults, as in Hoffman's plot_concordance) of
    B's p-values over A's restricted list (at least 10 genes), and B's all-gene baseline over B's
    genes in M_A; and the reverse. Where qvalue's default smoother fails (short or lumpy p vectors)
    pi0 is re-estimated at lambda = 0.5; those fields are named in pi1_lambda05. '_pub' columns use a published table: the study's own, else a published table from
    another study in its same_data_as group (same data, different analysis; FDR list first, then
    earliest PMID). '_any' columns fall back to our reproduction, then a Gemma result set (the
    study's Gemma alternate or a Gemma member of its group). A table from another list is used
    only when its brain region (region(), from the manifest contrast) matches the region of the
    study's primary list; each pair also records the two tables' regions (tables_same_region_*).

Eras: old = published <= 2015, new = >= 2016 (main); sensitivity split old <= 2020, new >= 2021.
Year from ../corpus/lineage.json (26697997, not a lineage node, from its label).

    python h1.py   -> h1_pairs.tsv, h1_summary.tsv, h1_tests.tsv, h1_fold_enrichment.svg
"""
import csv
import datetime
import itertools
import json
import math
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_universes as BU  # noqa: E402
import extract  # noqa: E402
import recurrence as R  # noqa: E402

SPLITS = (2015, 2020)
N_TOP = 100
RANK_COVERAGE = 0.8
MIN_DIR = 5        # shared directional genes for a pair's concordance to enter the medians
MIN_PI1 = 10       # A-list genes with a B p-value for pi1
Q_SIG = 0.05
NPERM = 2000
SEED = 20260924
PLOT_DIR = os.path.expanduser(".claude/skills/plotting/python")

SCENARIOS = [
    ("published, brain", lambda s: s["provenance"] == "published" and s["tissue"] == "brain"),
    ("published, brain, FDR lists", lambda s: s["provenance"] == "published" and s["tissue"] == "brain"
     and s["threshold_kind"] == "fdr"),
    ("all provenances, brain", lambda s: s["tissue"] == "brain"),
    ("published, blood", lambda s: s["provenance"] == "published" and s["tissue"] == "blood"),
    ("published, brain, microarray", lambda s: s["provenance"] == "published" and s["tissue"] == "brain"
     and s.get("technology") == "microarray"),
    ("published, brain, RNA-seq", lambda s: s["provenance"] == "published" and s["tissue"] == "brain"
     and s.get("technology") == "rnaseq"),
]


def read_list_stats(name):
    """symbol -> (direction, p, fdr), keeping the row with the smallest p (then fdr) per gene."""
    rows = [l for l in open(os.path.join(R.DATA, "lists", name + ".tsv"), encoding="utf-8") if not l.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        g = r.get("symbol")
        if not g:
            continue
        p, f = extract.num(r.get("p")), extract.num(r.get("fdr"))
        rec = (r.get("direction", ""), p, f)
        cur = out.get(g)
        if cur is None or ((p if p is not None else 9), (f if f is not None else 9)) < \
                ((cur[1] if cur[1] is not None else 9), (cur[2] if cur[2] is not None else 9)):
            out[g] = rec
    return out


def full_table(name, specs):
    """symbol -> (signed effect, p) from the list's full table (row with the smallest p per gene)."""
    rows, spec, _ = BU.fullstats_rows(specs[name], BU.fullstats_override(name))
    d = spec.get("direction") or {}
    out = {}
    for r in rows:
        g = r["symbol"]
        if not g:
            continue
        e, p = extract.num(r.get("effect")), extract.num(r.get("p"))
        if e is not None and d.get("from") == "sign":
            c = d.get("center", 0)
            # no inversion here: extract.build has already turned an inverted table's effect round (direction.invert)
            e = (math.log2(e) if e > 0 else None) if c == 1 else e - c
        if p is not None and not 0 <= p <= 1:
            p = None
        key = (p if p is not None else 2.0, -(abs(e) if e is not None else 0.0))
        if g not in out or key < out[g][2]:
            out[g] = (e, p, key)
    return {g: (e, p) for g, (e, p, _) in out.items()}


REGIONS = [  # (label, regex over the manifest contrast); a contrast matching several labels is 'multiple'
    ("prefrontal", r"\bDLPFC\b|prefrontal|\bPFC\b|\bfrontal\b|\bBA ?(8|9|10|11|46)\b|BA8-9|BA9/46|BA46/10|orbitofrontal"),
    ("cingulate", r"cingulate|\bACC\b|\bsgACC\b|\bAnCg\b|\bBA ?24\b"),
    ("temporal", r"temporal|\bSTG\b|\bBA ?(21|22|42)\b"),
    ("hippocampus", r"hippocamp|dentate|\bCA1\b"),
    ("striatum", r"caudate|striat|putamen|accumbens|\bnAcc\b"),
    ("amygdala", r"amygdala"),
    ("cerebellum", r"cerebell"),
    ("thalamus", r"thalam"),
    ("midbrain", r"midbrain"),
    ("visual cortex", r"visual cortex|\bBA ?17\b"),
]


def region(specs, name):
    """Brain region of a list: from its list id when that names exactly one region (e.g. __deg_dlpfc,
    gemma_..._striatum_...), else from its manifest contrast text; one label, 'multiple: a+b', or 'unknown'."""
    lid = name.split("__", 1)[-1].replace("_", " ")
    hits = [lab for lab, rx in REGIONS if re.search(rx, lid, re.I)]
    if len(hits) == 1:
        return hits[0]
    txt = specs.get(name, {}).get("contrast", "")
    hits = [lab for lab, rx in REGIONS if re.search(rx, txt, re.I)]
    if not hits:
        return "unknown"
    return hits[0] if len(hits) == 1 else "multiple: " + "+".join(hits)


def hyper_sf(x, N, K, n):
    """P(X >= x), X ~ hypergeometric(population N, K successes, n draws)."""
    if x <= 0:
        return 1.0
    hi = min(K, n)
    if x > hi:
        return 0.0
    lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    base = lc(N, n)
    terms = [lc(K, k) + lc(N - K, n - k) - base for k in range(x, hi + 1) if n - k <= N - K]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def rankdata(v):
    o = np.argsort(v, kind="mergesort")
    r = np.empty(len(v))
    sv = v[o]
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and sv[j + 1] == sv[i]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    x, y = rankdata(np.asarray(x, float)), rankdata(np.asarray(y, float))
    return float(np.corrcoef(x, y)[0, 1])


NULL_REPS = 1000


def null_median_fold(rows, rng, reps=NULL_REPS):
    """Median per-pair fold under the null: each pair's overlap redrawn from its hypergeometric (same universe and list
    sizes). With many small expected overlaps the per-pair fold is skewed, so this median sits below 1 and the observed
    median fold is read against it, not against 1 (tests/gene_lists_core/test_glc_h1_stats.py). Returns (mean, 2.5%,
    97.5%) over reps, or (None, None, None) when no pair has an expected overlap."""
    N = np.array([r["n_universe"] for r in rows])
    a = np.array([r["size_a"] for r in rows])
    b = np.array([r["size_b"] for r in rows])
    e = a * b / np.maximum(N, 1)
    keep = e > 0
    if not keep.any():
        return None, None, None
    N, a, b, e = N[keep], a[keep], b[keep], e[keep]
    meds = np.median(rng.hypergeometric(np.tile(a, (reps, 1)), np.tile(N - a, (reps, 1)), np.tile(b, (reps, 1))) / e, axis=1)
    return float(np.mean(meds)), float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def top_genes(stats, genes, n):
    """The n genes of `genes` with the smallest p (or fdr when fewer than 80% carry a p); None if unranked."""
    genes = sorted(genes)
    for idx in (1, 2):
        have = [g for g in genes if stats[g][idx] is not None]
        if genes and len(have) >= RANK_COVERAGE * len(genes):
            key = lambda g: (stats[g][idx] if stats[g][idx] is not None else 9, g)
            return set(sorted(genes, key=key)[:n]), ("p" if idx == 1 else "fdr")
    return None, ""


def run_pi1(jobs):
    """jobs: id -> list of p. Returns id -> pi1 (None where qvalue fails)."""
    if not jobs:
        return {}
    with tempfile.TemporaryDirectory() as td:
        fin, fout = os.path.join(td, "in.tsv"), os.path.join(td, "out.tsv")
        with open(fin, "w") as fh:
            fh.write("job\tp\n")
            for j, ps in jobs.items():
                fh.writelines(f"{j}\t{p!r}\n" for p in ps)
        rs = ("suppressMessages(library(qvalue)); d <- read.delim(commandArgs(TRUE)[1], colClasses=c('character','numeric'));"
              "s <- split(d$p, d$job); f <- function(p) tryCatch(c(1 - qvalue(p)$pi0, 0), error=function(e)"
              " tryCatch(c(1 - qvalue(p, lambda=0.5)$pi0, 1), error=function(e) c(NA, NA)));"
              "r <- t(sapply(s, f)); write.table(data.frame(job=names(s), pi1=r[,1], fallback=r[,2]), commandArgs(TRUE)[2],"
              " sep='\\t', row.names=FALSE, quote=FALSE)")
        subprocess.run(["Rscript", "-e", rs, fin, fout], check=True, capture_output=True)
        out = {}
        for r in csv.DictReader(open(fout), delimiter="\t"):
            out[r["job"]] = (None, None) if r["pi1"] in ("NA", "") else (float(r["pi1"]), r["fallback"] == "1")
    return out


def med(v):
    v = [x for x in v if x is not None and not (isinstance(x, float) and math.isnan(x))]
    return (float(np.median(v)), len(v)) if v else (None, 0)


def fmt(x, d=3):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{d}g}"
    return str(x)


def era_tests(scen, rows, rng, nperm=NPERM, splits=SPLITS):
    """The era tests of one scenario (h1_tests.tsv rows): new-new minus old-old for four statistics, all pairs and
    different-cohort pairs, with era labels permuted across studies (pairs share studies, so pair-level tests are not
    valid); and the size model (least squares of log2 fold), with the new-new coefficient's p from the same permutation.
    rows: the scenario's pair rows (study_a, study_b, year_a, year_b, cohort, size_a, size_b, log2_fold, overlap,
    expected, top_overlap, top_expected, concordance, n_dir). Draws from rng in the order main() always has."""
    tests = []
    sids = sorted({r["study_a"] for r in rows} | {r["study_b"] for r in rows})
    yr = {r["study_a"]: r["year_a"] for r in rows} | {r["study_b"]: r["year_b"] for r in rows}
    ia = np.array([sids.index(r["study_a"]) for r in rows])
    ib = np.array([sids.index(r["study_b"]) for r in rows])
    coh_same = np.array([r["cohort"] == "same" for r in rows])
    smin = np.log10(np.array([max(1, min(r["size_a"], r["size_b"])) for r in rows]))
    smax = np.log10(np.array([max(1, max(r["size_a"], r["size_b"])) for r in rows]))
    lf = np.array([r["log2_fold"] for r in rows])
    ov = np.array([r["overlap"] for r in rows], float)
    ex = np.array([r["expected"] for r in rows], float)
    tov = np.array([r.get("top_overlap", np.nan) for r in rows], float)
    tex = np.array([r.get("top_expected", np.nan) for r in rows], float)
    cc = np.array([r["concordance"] if r["n_dir"] >= MIN_DIR else np.nan for r in rows], float)
    # each metric: (values that must be present, statistic over a mask)
    metrics = {
        "median log2 fold": (lf, lambda m: np.median(lf[m])),
        "pooled fold (sum overlap / sum expected)": (ex, lambda m: ov[m].sum() / ex[m].sum()),
        "pooled top-N fold": (tex, lambda m: tov[m].sum() / tex[m].sum()),
        "median concordance": (cc, lambda m: np.median(cc[m])),
    }
    for s in splits:
        new = np.array([yr[k] > s for k in sids])
        for coh in ("all", "different"):
            mask = np.ones(len(rows), bool) if coh == "all" else np.array([r["cohort"] == "different" for r in rows])
            for mname, (v, stat) in metrics.items():
                def diff(nw):
                    nn = nw[ia] & nw[ib] & mask & ~np.isnan(v)
                    oo = ~nw[ia] & ~nw[ib] & mask & ~np.isnan(v)
                    if nn.sum() == 0 or oo.sum() == 0:
                        return np.nan, int(nn.sum()), int(oo.sum()), np.nan, np.nan
                    a_, b_ = float(stat(nn)), float(stat(oo))
                    return a_ - b_, int(nn.sum()), int(oo.sum()), a_, b_
                obs, nnn, noo, vnn, voo = diff(new)
                if math.isnan(obs):
                    continue
                perm = np.array([diff(rng.permutation(new))[0] for _ in range(nperm)])
                perm = perm[~np.isnan(perm)]
                p = (1 + np.sum(np.abs(perm) >= abs(obs) - 1e-12)) / (1 + len(perm))
                tests.append({"scenario": scen, "era_split": s, "cohort": coh, "test": "new-new minus old-old",
                              "metric": mname, "new_new": vnn, "old_old": voo, "estimate": obs,
                              "n_new_new": nnn, "n_old_old": noo, "p_perm": p, "n_perm": len(perm)})
        # size model: log2 fold ~ new-new + old-new + same cohort + log10 smaller + log10 larger list size;
        # a term with no variation in the scenario (no same-cohort pairs, no new-new pairs) is left out
        nn0 = new[ia] & new[ib]
        if nn0.sum() < 3 or (~new[ia] & ~new[ib]).sum() < 3:
            continue
        names = ["intercept", "new-new", "old-new", "same cohort", "log10 smaller list", "log10 larger list"]
        use = [True, True, bool((new[ia] ^ new[ib]).any()), bool(coh_same.any() and not coh_same.all()), True, True]

        def fit(nw, ret_all=False):
            cols = [np.ones(len(rows)), (nw[ia] & nw[ib]).astype(float), (nw[ia] ^ nw[ib]).astype(float),
                    coh_same.astype(float), smin, smax]
            X = np.column_stack([c for c, u in zip(cols, use) if u])
            beta = np.linalg.lstsq(X, lf, rcond=None)[0]
            return beta if ret_all else beta[1]
        beta = fit(new, True)
        perm = np.array([fit(rng.permutation(new)) for _ in range(nperm)])
        p = (1 + np.sum(np.abs(perm) >= abs(beta[1]) - 1e-12)) / (1 + nperm)
        for nm, b in zip([n for n, u in zip(names, use) if u], beta):
            tests.append({"scenario": scen, "era_split": s, "cohort": "all", "test": "OLS log2 fold",
                          "metric": nm, "new_new": "", "old_old": "", "estimate": float(b), "n_new_new": int(nn0.sum()),
                          "n_old_old": int((~new[ia] & ~new[ib]).sum()),
                          "p_perm": p if nm == "new-new" else "", "n_perm": nperm if nm == "new-new" else ""})
    return tests


def main():
    extract.HGNC = extract.load_hgnc()
    today = datetime.date.today().isoformat()
    pl = json.load(open(os.path.join(HERE, "primary_lists.json")))
    uni = json.load(open(os.path.join(HERE, "universes.json")))["lists"]
    nodes = json.load(open(R.LINEAGE))["nodes"]
    lanes = {(n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")): n.get("lane") for n in nodes}
    years = {(n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")): n.get("year") for n in nodes}
    specs = BU.manifest_specs()
    units = R.load_units(pl, nodes)
    studies = units["studies"]
    groups = {k: units["groups"][units["group_of"][k]] for k in studies}

    def year(k):
        y = years.get(k)
        if y:
            return int(y)
        m = re.search(r"(19|20)\d\d", studies[k]["label"])
        return int(m.group(0))

    def cluster(k):
        c = R.cluster_of(units, k)
        return None if ":" in c else c

    def clusters_of(k):
        """Named clusters of study k: one, none (own cohort), or a multi-bank study's per-bank clusters;
        None for a multi-bank study without per-bank lists."""
        if k in units["multibank"]:
            bl = units["bank_lists"].get(k)
            return {R.bank_cluster(units, b) for b in bl} if bl else None
        c = cluster(k)
        return {c} if c else set()

    def relation(a, b):
        xa, xb = clusters_of(a), clusters_of(b)
        if xa is None or xb is None:
            return "multi-bank"
        return "same" if xa & xb else "different"

    # lists, universes (as recurrence.py: every primary list of a study must share one universe)
    stats, M, set_cache = {}, {}, {}
    for k in studies:
        stats[k] = {}
        for name in studies[k]["primary"]:
            for g, rec in read_list_stats(name).items():
                if g not in stats[k] or (rec[1] or 9) < (stats[k][g][1] or 9):
                    stats[k][g] = rec
        u = [uni.get(n) for n in studies[k]["primary"]]
        if all(x and x.get("gene_set") for x in u):
            for x in u:
                if x["gene_set"] not in set_cache:
                    set_cache[x["gene_set"]] = R.read_set(x["gene_set"])
            if len({frozenset(set_cache[x["gene_set"]]) for x in u}) != 1:
                raise SystemExit(f"{k}: primary lists with different universes")
            M[k] = set_cache[u[0]["gene_set"]]

    # full tables: published (own, else same_data_as group member), then any (reproduced own, Gemma).
    # A table from another list must be of the same brain region as the study's primary list (H1-3);
    # candidates from a different region are not used and are counted in region_skipped.
    def is_full(n):
        return n in uni and uni[n]["type"] == "fullstats"
    lreg = {k: region(specs, studies[k]["primary"][0]) for k in studies}

    def same_reg(k, n):
        a, b = lreg[k], region(specs, n)
        return a == b or "unknown" in (a, b)
    src_pub, src_any, region_skipped = {}, {}, []
    for k, grp in groups.items():
        own = [n for n in studies[k]["primary"] if is_full(n)]
        prov = studies[k]["provenance"]
        if own and prov == "published":
            src_pub[k] = (own[0], "own")
        else:
            # a published member's table; picked as the representative is picked: FDR list first, then earliest PMID
            allc = [(studies[m]["threshold_kind"] != "fdr", int(m), m, n) for m in grp if m != k
                    and studies[m]["provenance"] == "published" for n in studies[m]["primary"][:1] if is_full(n)]
            cands = [c for c in allc if same_reg(k, c[3])]
            region_skipped += [(k, c[3]) for c in allc if not same_reg(k, c[3])]
            if cands:
                *_, m, n = min(cands)
                src_pub[k] = (n, f"group member {studies[m]['label']}")
        if k in src_pub:
            src_any[k] = src_pub[k]
        elif own:
            src_any[k] = (own[0], f"own, {prov}")
        else:
            alts = [n for n in studies[k]["alternates"] if specs.get(n, {}).get("provenance") == "gemma"]
            gm = [studies[m]["primary"][0] for m in grp if m != k and studies[m]["provenance"] == "gemma"
                  and is_full(studies[m]["primary"][0])]
            region_skipped += [(k, n) for n in alts + gm if not same_reg(k, n)]
            alts = [n for n in alts if same_reg(k, n)]
            gm = [n for n in gm if same_reg(k, n)]
            if alts:
                src_any[k] = (alts[0], "Gemma alternate")
            elif gm:
                src_any[k] = (gm[0], "Gemma, group member")
    tables = {}

    def table(n):
        if n not in tables:
            tables[n] = full_table(n, specs)
        return tables[n]

    pairs, pi1_jobs, n_reps = [], {}, {}
    for scen, f in SCENARIOS:
        reps = R.representatives(units, f, lambda k: bool(stats[k]))
        n_reps[scen] = len(reps)
        sel = sorted((k for k in reps if k in M), key=lambda k: (year(k), (0, int(k), "") if k.isdigit() else (1, 0, k)))
        for A, B in itertools.combinations(sel, 2):
            U = M[A] & M[B]
            a = {g for g in stats[A] if g in U}
            b = {g for g in stats[B] if g in U}
            if not U:
                continue
            ov = a & b
            exp_ = len(a) * len(b) / len(U)
            ca, cb = cluster(A), cluster(B)
            row = {"scenario": scen, "study_a": A, "label_a": studies[A]["label"], "year_a": year(A),
                   "study_b": B, "label_b": studies[B]["label"], "year_b": year(B),
                   "cluster_a": ca or f"{lanes.get(A) or 'none'}:{A}", "cluster_b": cb or f"{lanes.get(B) or 'none'}:{B}",
                   "cohort": relation(A, B),
                   "region_a": lreg[A], "region_b": lreg[B],
                   "threshold_a": studies[A]["threshold_kind"], "threshold_b": studies[B]["threshold_kind"],
                   "technology_a": studies[A]["technology"], "technology_b": studies[B]["technology"],
                   "tech_pair": "-".join(sorted((studies[A]["technology"], studies[B]["technology"]))),
                   "universe_type_a": uni[studies[A]["primary"][0]]["type"], "universe_type_b": uni[studies[B]["primary"][0]]["type"],
                   "n_universe": len(U), "size_a": len(a), "size_b": len(b), "overlap": len(ov),
                   "expected": exp_, "fold": len(ov) / exp_ if exp_ else None,
                   "log2_fold": math.log2((len(ov) + 0.5) / (exp_ + 0.5)),
                   "p_hyper": hyper_sf(len(ov), len(U), len(a), len(b)),
                   "jaccard": len(ov) / len(a | b) if a | b else None}
            n = min(len(a), len(b), N_TOP)
            ta, ka = top_genes(stats[A], a, n)
            tb, kb = top_genes(stats[B], b, n)
            if ta is not None and tb is not None and n:
                e2 = n * n / len(U)
                row.update(top_n=n, top_rank_by=f"{ka}/{kb}", top_overlap=len(ta & tb), top_expected=e2,
                           top_fold=len(ta & tb) / e2)
            same = sum(1 for g in ov if stats[A][g][0] in ("up", "down") and stats[A][g][0] == stats[B][g][0])
            opp = sum(1 for g in ov if {stats[A][g][0], stats[B][g][0]} == {"up", "down"})
            row.update(n_dir=same + opp, n_same=same, concordance=same / (same + opp) if same + opp else None,
                       p_concordance=R.sign_test(same, opp))
            for tag, srcs in (("pub", src_pub), ("any", src_any)):
                if A not in srcs or B not in srcs:
                    continue
                (na, sa), (nb, sb) = srcs[A], srcs[B]
                TA, TB = table(na), table(nb)
                ra_, rb_ = region(specs, na), region(specs, nb)
                row[f"table_a_{tag}"] = f"{na} ({sa}; {ra_})"
                row[f"table_b_{tag}"] = f"{nb} ({sb}; {rb_})"
                row[f"tables_same_region_{tag}"] = "unknown" if "unknown" in (ra_, rb_) else ("yes" if ra_ == rb_ else "no")
                sh = [g for g in TA if g in TB and TA[g][0] is not None and TB[g][0] is not None]
                row[f"n_rho_{tag}"] = len(sh)
                if len(sh) >= 3:
                    row[f"rho_{tag}"] = spearman([TA[g][0] for g in sh], [TB[g][0] for g in sh])
                for d, (X, L, T) in (("ab", (A, a, TB)), ("ba", (B, b, TA))):
                    ps = [T[g][1] for g in sorted(L) if g in T and T[g][1] is not None]
                    Mx = M[A] if d == "ab" else M[B]
                    base = [T[g][1] for g in sorted(T) if g in Mx and T[g][1] is not None]
                    row[f"n_pi1_{d}_{tag}"] = len(ps)
                    if len(ps) >= MIN_PI1:
                        j = f"{scen}|{A}|{B}|{d}|{tag}"
                        pi1_jobs[j] = ps
                        row[f"_job_{d}_{tag}"] = j
                    if len(base) >= MIN_PI1:
                        jb = f"base|{(nb if d == 'ab' else na)}|{X}"
                        pi1_jobs[jb] = base
                        row[f"_base_{d}_{tag}"] = jb
            pairs.append(row)
    pi1 = run_pi1(pi1_jobs)
    for row in pairs:
        fb = []
        for d in ("ab", "ba"):
            for tag in ("pub", "any"):
                for pre, col in (("_job", "pi1"), ("_base", "pi1_base")):
                    if f"{pre}_{d}_{tag}" in row:
                        v, f_ = pi1.get(row.pop(f"{pre}_{d}_{tag}"), (None, None))
                        row[f"{col}_{d}_{tag}"] = v
                        if f_:
                            fb.append(f"{col}_{d}_{tag}")
        row["pi1_lambda05"] = ",".join(fb)
    # BH across pairs within each scenario; era labels
    for scen, _ in SCENARIOS:
        rows = [r for r in pairs if r["scenario"] == scen]
        if rows:
            for r, q in zip(rows, R.bh([r["p_hyper"] for r in rows])):
                r["q_hyper"] = float(q)
    for r in pairs:
        for s in SPLITS:
            ea, eb = ("old" if r["year_a"] <= s else "new"), ("old" if r["year_b"] <= s else "new")
            r[f"era_{s}"] = "-".join(sorted((ea, eb), key=lambda e: e != "old")) if ea != eb else f"{ea}-{ea}"

    # summaries
    summ = []
    null_rng = np.random.default_rng(SEED + 1)
    types = ("old-old", "new-new", "old-new")
    for scen, _ in SCENARIOS:
        for s in SPLITS:
            for t in types + ("all",):
                for coh in ("same", "different", "all"):
                    rows = [r for r in pairs if r["scenario"] == scen and (t == "all" or r[f"era_{s}"] == t)
                            and (coh == "all" or r["cohort"] == coh)]
                    if not rows:
                        continue
                    fold, nf = med([r["fold"] for r in rows])
                    nfold, nfold_lo, nfold_hi = null_median_fold(rows, null_rng)
                    pooled = sum(r["overlap"] for r in rows) / sum(r["expected"] for r in rows)
                    trs = [r for r in rows if r.get("top_fold") is not None]
                    pooled_top = sum(r["top_overlap"] for r in trs) / sum(r["top_expected"] for r in trs) if trs else None
                    tf, ntf = med([r.get("top_fold") for r in rows])
                    conc, nc = med([r["concordance"] for r in rows if r["n_dir"] >= MIN_DIR])
                    sd, nd = sum(r["n_same"] for r in rows), sum(r["n_dir"] for r in rows)
                    rp, nrp = med([r.get("rho_pub") for r in rows])
                    ra, nra = med([r.get("rho_any") for r in rows])
                    pa, npa = med([r.get("pi1_ab_any") for r in rows])
                    pb, npb = med([r.get("pi1_base_ab_any") for r in rows])
                    summ.append({"scenario": scen, "era_split": s, "pair_type": t, "cohort": coh, "n_pairs": len(rows),
                                 "n_studies": len({r["study_a"] for r in rows} | {r["study_b"] for r in rows}),
                                 "median_size_smaller": med([min(r["size_a"], r["size_b"]) for r in rows])[0],
                                 "median_size_larger": med([max(r["size_a"], r["size_b"]) for r in rows])[0],
                                 "median_fold": fold, "null_median_fold": nfold, "null_median_fold_lo": nfold_lo,
                                 "null_median_fold_hi": nfold_hi, "pooled_fold": pooled, "n_pairs_overlap0": sum(1 for r in rows if r["overlap"] == 0),
                                 "median_top_fold": tf, "pooled_top_fold": pooled_top, "n_top": ntf,
                                 "n_top_overlap_gt0": sum(1 for r in trs if r["top_overlap"] > 0),
                                 "median_concordance": conc, "n_concordance": nc,
                                 "pooled_concordance": sd / nd if nd else None, "pooled_n_dir": nd,
                                 "median_rho_pub": rp, "n_rho_pub": nrp, "median_rho_any": ra, "n_rho_any": nra,
                                 "median_pi1_ab_any": pa, "median_pi1_base_ab_any": pb, "n_pi1_any": npa,
                                 "n_sig_q05": sum(1 for r in rows if r["q_hyper"] < Q_SIG)})

    # tests: study-level permutation of era labels (pairs share studies, so pair-level tests are not valid)
    rng = np.random.default_rng(SEED)
    tests = []
    for scen, _ in SCENARIOS:
        rows = [r for r in pairs if r["scenario"] == scen]
        if not rows:
            continue
        tests += era_tests(scen, rows, rng)

    # write
    def sha(p):
        return R.sha(p)
    head = [f"# script: h1.py, run {today}",
            f"# input: primary_lists.json sha256 {sha(os.path.join(HERE, 'primary_lists.json'))} (rule A: DE and DE-adjacent lists, Paul 2026-09-25)",
            f"# input: universes.json sha256 {sha(os.path.join(HERE, 'universes.json'))}; lists/*.tsv; full tables via build_universes.fullstats_rows",
            f"# input: ../corpus/lineage.json sha256 {sha(R.LINEAGE)} (year, lane)",
            f"# units: recurrence.load_units ({len(units['all'])} studies with a primary list; left out {len(units['meta'])} meta-lane and "
            f"{len(units['multi'])} multi-series reanalyses, {len(units['reanalysis'])} lineage group 3 reanalyses and {len(units['role_only'])} "
            f"studies whose only primary list has role reanalysis; {len(units['groups'])} same_data_as groups; representative per scenario: "
            + ", ".join(f"{k} {v}" for k, v in n_reps.items()) + "); pairs over M_A & M_B; A is the earlier study",
            f"# list role: {len(units['role_dropped'])} studies had a primary list with role reanalysis dropped"
            + (" (" + ", ".join(units["labels"][k] for k in sorted(units["role_dropped"])) + ")" if units["role_dropped"] else "")
            + f"; multi-bank studies: {len(units['multibank'])}"
            + (" (" + ", ".join(f"{units['labels'][k]} {'+'.join(units['multibank'][k])}: " + ("per-bank lists" if k in units["bank_lists"]
                                else "no per-bank lists, its pairs have cohort 'multi-bank'") for k in sorted(units["multibank"])) + ")"
               if units["multibank"] else ""),
            "# cohort: recurrence.cluster_of (named lanes of all group members); region: brain region parsed from the list id, else "
            "the manifest contrast (region()); a full table from another list is used only when its region matches the study's list region "
            f"({len(region_skipped)} candidate tables of another region skipped: "
            + ("; ".join(f"{studies[k]['label']} not scored with {n}" for k, n in region_skipped) or "none") + ")",
            f"# era: old <= {SPLITS[0]} / new >= {SPLITS[0] + 1} (era_{SPLITS[0]}), sensitivity old <= {SPLITS[1]} (era_{SPLITS[1]})",
            "# p_hyper one-sided hypergeometric; q_hyper BH across pairs within scenario; top_*: N = min(|a|,|b|,100) genes by p (or fdr)",
            "# concordance over shared genes with a direction in both lists; p_concordance two-sided exact binomial",
            "# rho: logFC Spearman correlation (signed effect, no table carries a t statistic) over genes in both full tables; tables_same_region_*: the two tables' regions; pi1_ab = 1 - qvalue pi0 of B's p over A's list, "
            "pi1_base_ab = over B's genes in M_A (qvalue defaults; lambda = 0.5 where the default fails, named in pi1_lambda05); _pub published tables only (own or same_data_as group member), _any adds reproduced and Gemma"]
    fields = list(dict.fromkeys(k for r in pairs for k in r))
    with open(os.path.join(HERE, "h1_pairs.tsv"), "w") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", restval="")
        w.writeheader()
        for r in pairs:
            w.writerow({k: fmt(v, 4) for k, v in r.items()})
    for fn, rows, extra in (("h1_summary.tsv", summ, [f"# medians over pairs; concordance median over pairs with >= {MIN_DIR} shared directional genes; "
                                                         f"pooled = sum agreeing / sum shared directional genes; n_sig_q05 = pairs with q_hyper < {Q_SIG}"]),
                            ("h1_tests.tsv", tests, [f"# p_perm: era labels permuted across studies ({NPERM} permutations, seed {SEED}); "
                                                     "pair-level tests would treat pairs sharing a study as independent"])):
        with open(os.path.join(HERE, fn), "w") as fh:
            fh.write("\n".join(head + extra) + "\n")
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", restval="")
            w.writeheader()
            for r in rows:
                w.writerow({k: fmt(v, 4) for k, v in r.items()})
    figure([r for r in pairs if r["scenario"] == "published, brain"], today)
    for r in summ:
        if r["era_split"] == SPLITS[0] and r["scenario"] in ("published, brain", "published, blood"):
            print({k: fmt(v) for k, v in r.items()})
    print("pi1 jobs", len(pi1_jobs), "lambda 0.5 fallback", sum(1 for v in pi1.values() if v[1]),
          "failed", sum(1 for v in pi1.values() if v[0] is None))


def figure(rows, today):
    sys.path.insert(0, PLOT_DIR)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pavlab_stripchart import pavlab_stripchart
    order = ["old-old", "old-new", "new-new"]
    cnt = {t: sum(1 for r in rows if r["era_2015"] == t) for t in order}
    x = [f"{r['era_2015']}\n(n={cnt[r['era_2015']]})" for r in rows]
    col = [f"{r['cohort']} cohort" for r in rows]
    fig, ax = plt.subplots(figsize=(6, 4.8))
    xorder = [f"{t}\n(n={cnt[t]})" for t in order]
    pavlab_stripchart(x, [r["log2_fold"] for r in rows], color=col, order=xorder, ax=ax,
                      ylabel="log2 fold enrichment", xlabel="Era pair (old <= 2015)", label_size="small",
                      title="Pairwise overlap, published brain lists", origin_zero=False)
    ax.axhline(0, color="#9ca3af", lw=0.8, zorder=0)
    fig.text(0.01, 0.01, f"Fold = (obs + 0.5) / (exp + 0.5) over the pair's shared tested universe.\n"
             f"Source: h1_pairs.tsv (h1.py, {today}), scenario 'published, brain'", fontsize=8, color="#6b7280")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(HERE, "h1_fold_enrichment.svg"))


def _tsv(fn):
    lines = open(os.path.join(HERE, fn), encoding="utf-8").read().splitlines()
    return [l for l in lines if l.startswith("#")], list(csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t"))


def _f(x, d=2):
    if x in ("", None):
        return "none"
    v = float(x)
    return f"{v:.{d}f}".replace("-", "\u2212") if abs(v) >= 0.01 or v == 0 else f"{v:.3f}".replace("-", "\u2212")


def _p(x):
    v = float(x)
    return f"{v:.2g}" if v < 0.01 else f"{v:.2f}"


def tech_section(S, T, main):
    """Microarray-only and RNA-seq-only scenarios, and main-scenario pairs by technology pair (technology.json)."""
    def pooled(rows, o, e):
        rs = [r for r in rows if r.get(o) not in ("", None)]
        return sum(float(r[o]) for r in rs) / sum(float(r[e]) for r in rs) if rs else None

    def medv(rows, k):
        v = sorted(float(r[k]) for r in rows if r.get(k) not in ("", None))
        return (v[len(v) // 2] if len(v) % 2 else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2) if v else None
    L = ["## By technology: microarray vs RNA-seq", "",
         "Technology per study from `technology.json` (`classify_technology.py`: GEO platform, RNA-seq universe, or the lineage "
         "platform text, with listed overrides). Single-nucleus and imputed (TWAS) lists are outside both scenarios. Each scenario "
         "picks its own representative per same_data_as group, so the two scenarios can each hold a member of one group.", "",
         "| Scenario | Pair type | Cohort | Pairs | Median fold | Pooled fold | Pooled top-N fold | Median concordance (pairs) | "
         "Median ρ, any table (pairs) | Significant, q < 0.05 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for scen, lab in (("published, brain, microarray", "Microarray"), ("published, brain, RNA-seq", "RNA-seq")):
        first = True
        for t in ("old-old", "new-new", "old-new", "all"):
            for c in ("different", "all"):
                x = S(scen, 2015, t, c)
                if not x:
                    continue
                L.append(f"| {lab if first else ''} | {t} | {c} | {x['n_pairs']} | {_f(x['median_fold'])} | {_f(x['pooled_fold'])} | "
                         f"{_f(x['pooled_top_fold'])} | {_f(x['median_concordance'])} ({x['n_concordance']}) | "
                         + ("none" if not x["median_rho_any"] else f"{_f(x['median_rho_any'], 3)} ({x['n_rho_any']})")
                         + f" | {x['n_sig_q05']} |")
                first = False
    L += ["", "Era tests within technology (new-new minus old-old, study-level permutation p; old up to 2015):", ""]
    for scen, lab in (("published, brain, microarray", "Microarray"), ("published, brain, RNA-seq", "RNA-seq")):
        parts = []
        for m, ml in (("median log2 fold", "median log2 fold"), ("pooled top-N fold", "pooled top-N fold"),
                      ("median concordance", "concordance")):
            a = T(scen, 2015, "all", m)
            if a and a.get("p_perm") not in ("", None):
                parts.append(f"{ml} {_f(a['new_new'])} vs {_f(a['old_old'])} (p {_p(a['p_perm'])})")
        L.append(f"- {lab}: " + ("; ".join(parts) if parts else "not testable (one era has no pairs)") + ".")
    L += ["", "Main scenario (published brain) pairs by technology pair, all eras:", "",
          "| Technology pair | Cohort | Pairs | Median fold | Pooled fold | Pooled top-N fold | Median concordance | Pairs sharing no gene |",
          "|---|---|---|---|---|---|---|---|"]
    for tp in sorted({p.get("tech_pair", "") for p in main}):
        for c in ("different", "all"):
            rows = [p for p in main if p.get("tech_pair", "") == tp and (c == "all" or p["cohort"] == c)]
            if not rows:
                continue
            conc = [r for r in rows if r.get("n_dir") not in ("", None) and int(r["n_dir"]) >= MIN_DIR]
            L.append(f"| {tp or 'unrecorded'} | {c} | {len(rows)} | {_f(medv(rows, 'fold'))} | {_f(pooled(rows, 'overlap', 'expected'))} | "
                     f"{_f(pooled(rows, 'top_overlap', 'top_expected'))} | {_f(medv(conc, 'concordance'))} ({len(conc)}) | "
                     f"{sum(1 for r in rows if r['overlap'] == '0')} |")
    return L + [""]


def write_report():
    """h1_report.md from h1_summary.tsv, h1_tests.tsv and h1_pairs.tsv only."""
    head, summ = _tsv("h1_summary.tsv")
    _, tests = _tsv("h1_tests.tsv")
    _, pairs = _tsv("h1_pairs.tsv")
    run = head[0].split("run ")[-1]
    shas = {fn: R.sha(os.path.join(HERE, fn)) for fn in ("h1_pairs.tsv", "h1_summary.tsv", "h1_tests.tsv")}
    S = lambda scen, split, t, c: next((r for r in summ if r["scenario"] == scen and r["era_split"] == str(split)
                                         and r["pair_type"] == t and r["cohort"] == c), None)
    T = lambda scen, split, coh, metric, test="new-new minus old-old": next(
        (r for r in tests if r["scenario"] == scen and r["era_split"] == str(split) and r["cohort"] == coh
         and r["metric"] == metric and r["test"] == test), None)
    units = [h[2:] for h in head if h.startswith("# units:")][0]
    roleln = [h[2:] for h in head if h.startswith("# list role:")][0]
    cohort = [h[2:] for h in head if h.startswith("# cohort:")][0]
    skipped = cohort.split("skipped: ", 1)[-1].rstrip(")")
    main = [p for p in pairs if p["scenario"] == "published, brain"]
    nst = len({p["study_a"] for p in main} | {p["study_b"] for p in main})
    L = ["---", "type: detail", "bib: My Library", f"modified: {run}", "---", "", "> [!note] Paul", ">", "",
         "# H1: pairwise agreement, old era vs new era", "",
         f"Generated by `h1.py` from `h1_summary.tsv` (sha256 {shas['h1_summary.tsv']}), `h1_tests.tsv` ({shas['h1_tests.tsv']}) and "
         f"`h1_pairs.tsv` ({shas['h1_pairs.tsv']}), run {run}. Edit the script, not this file.", "",
         "Tracker [[Schizophrenia transcriptome replication — tracker]], Step 3, H1: \"Within-era pairwise overlap (hypergeometric against the tested "
         "universe) and directional concordance; old-era pairs vs new-era pairs. Overlap depends strongly on list size, so match "
         "thresholds. Where full summary statistics exist, also compute correlation of test statistics (the metric Hoffman et al. "
         "proposed). Shared-sample pairs analysed separately.\"", "", "Inputs:", ""]
    L += [f"- {h[2:].replace('input: ', '')}" for h in head if h.startswith("# input:")]
    L += ["", "```", "python h1.py   # h1_pairs.tsv, h1_summary.tsv, h1_tests.tsv, h1_fold_enrichment.svg, this report", "```", "",
          "## Method", "",
          f"**Units.** As `recurrence.py` (`recurrence.load_units`, `recurrence.representatives`): {units.split('(', 1)[1].split(')')[0]}. "
          "Pairs that share data never enter. Only lists with a gene universe are used. Cohort relation uses `recurrence.cluster_of`: "
          "a same_data_as group takes the named lanes of all its members; a pair is \"same cohort\" when both studies are in one named "
          "cluster (cmc+mssm+pitt+libd, stanley, aus, harvard, pritzker, ipsc); other studies are each their own cohort. A multi-bank "
          "study (lineage `banks`) takes the clusters of its per-bank lists; without per-bank lists its pairs have cohort \"multi-bank\" "
          "and enter only the all-cohort rows and tests. Only first_look lists (the study's own new data) represent a study; "
          "reanalysis-role lists are a separate category. This run: " + roleln[len("list role: "):] + ".", "",
          f"**Scenarios.** Main: published lists, brain ({nst} studies, {len(main)} pairs). Sensitivities: published brain FDR-controlled "
          "lists only, all provenances brain (adds reproduced and Gemma lists where a group has no published list), published blood.", "",
          "**Eras.** Old = published up to 2015, new = 2016 on, from lineage.json (`26697997`, not a lineage node, from its label). "
          "Sensitivity: old up to 2020, new 2021 on.", "",
          "**Per pair** (A = the earlier study), over U = M_A ∩ M_B, the genes both studies tested:", "",
          "- overlap, expected |a||b|/|U| with a and b the lists restricted to U, one-sided hypergeometric p, Jaccard, fold = observed / "
          "expected; BH across the pairs of a scenario.",
          "- Top-N overlap: N = min(|a|, |b|, 100); each list's N genes with the smallest p (FDR where the list has no p), when at least "
          "80% of its genes carry one; expected N²/|U|.",
          "- Concordance: fraction of shared genes with a direction in both lists that agree; two-sided exact binomial p. Medians use "
          "pairs with at least 5 such genes.",
          "- Full statistics: ρ is the logFC Spearman correlation (signed effect; none of the tables has a t statistic; Ellis 2016 has a "
          "z score) over genes in both tables. Hoffman et al.'s 0.343 is also a logFC Spearman correlation (their `plot_concordance` "
          "calls `cor.test(logFC.x, logFC.y, method = \"spearman\")`), so the comparison is like for like. π1 = 1 − `qvalue(p)$pi0` "
          "(defaults, as in Hoffman's `plot_concordance`) of B's p-values over A's list (at least 10 genes), beside B's all-gene baseline "
          "over B's genes in M_A; λ = 0.5 where the default smoother fails (`pi1_lambda05`). `_pub` columns use a published table: the "
          "study's own, else the published table of another member of its same_data_as group; `_any` adds our reproductions and Gemma "
          "result sets.",
          "- Region: each list's brain region is read from its list id, else its manifest contrast (`region()` in `h1.py`). A table from "
          "another list is used only when its region matches the region of the study's primary list. Tables are labelled by list and "
          "region in `h1_pairs.tsv` (`table_a_pub`, `table_b_pub`, ...), and `tables_same_region_*` records whether the two tables of a "
          f"pair are from one region. Candidate tables not used because of region: {skipped}.", "",
          "**Tests.** Pairs share studies, so the era comparison permutes era labels across studies (2,000 permutations) and "
          "recomputes new-new minus old-old for median log2 fold (pseudocount 0.5), pooled fold (sum of overlaps / sum of expected), "
          "pooled top-N fold and median concordance. List size: least squares of log2 fold on new-new, old-new, same cohort, log10 "
          "smaller and log10 larger list size, with the new-new coefficient's p from the same permutation.", "",
          "## Main result: published brain, old up to 2015", "",
          "Median fold is read against its null median (in brackets: each pair's overlap redrawn from its hypergeometric, "
          f"{NULL_REPS:,} draws), not against 1: many pairs expect under one shared gene, so the per-pair fold is skewed and its "
          "null median sits below 1. Pooled fold centres on 1 under the null.", "",
          "| Pair type | Cohort | Pairs | Median fold (null) | Pooled fold | Pooled top-N fold | Median concordance (pairs) | Median ρ, published "
          "tables (pairs) | Median ρ, any table (pairs) | Significant, q < 0.05 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for t in ("old-old", "new-new", "old-new", "all"):
        for c in (("same", "different", "all") if t != "all" else ("all",)):
            r = S("published, brain", 2015, t, c)
            if not r:
                continue
            rho = lambda k, n: "none" if not r[k] else f"{_f(r[k], 3)} ({r[n]})"
            L.append(f"| {t} | {c} | {r['n_pairs']} | {_f(r['median_fold'])} ({_f(r['null_median_fold'])}) | {_f(r['pooled_fold'])} | {_f(r['pooled_top_fold'])} | "
                     f"{_f(r['median_concordance'])} ({r['n_concordance']}) | {rho('median_rho_pub', 'n_rho_pub')} | "
                     f"{rho('median_rho_any', 'n_rho_any')} | {r['n_sig_q05']} |")
    r = S("published, brain", 2015, "all", "all")
    L += ["", f"{r['n_pairs_overlap0']} of the {r['n_pairs']} pairs share no gene. Median list sizes after restriction (smaller / larger "
          "list of a pair): " + ", ".join(f"{t} {_f(S('published, brain', 2015, t, 'all')['median_size_smaller'], 0)} / "
                                         f"{_f(S('published, brain', 2015, t, 'all')['median_size_larger'], 0)}" for t in ("old-old", "new-new", "old-new")) + ".", "",
          "Era tests (new-new minus old-old, study-level permutation p):", "",
          "| Statistic | All pairs: new-new vs old-old | p | Different-cohort pairs only | p |", "|---|---|---|---|---|"]
    for m, lab in (("median log2 fold", "Median log2 fold"), ("pooled fold (sum overlap / sum expected)", "Pooled fold"),
                   ("pooled top-N fold", "Pooled top-N fold"), ("median concordance", "Median concordance")):
        a, d = T("published, brain", 2015, "all", m), T("published, brain", 2015, "different", m)
        L.append(f"| {lab} | {_f(a['new_new'])} vs {_f(a['old_old'])} | {_p(a['p_perm'])} | "
                 + (f"{_f(d['new_new'])} vs {_f(d['old_old'])} | {_p(d['p_perm'])} |" if d else "none | |"))
    ols = [t for t in tests if t["scenario"] == "published, brain" and t["era_split"] == "2015" and t["test"] == "OLS log2 fold"]
    nn = next(t for t in ols if t["metric"] == "new-new")
    L += ["", f"Size model (log2 fold, {r['n_pairs']} pairs): new-new {_f(nn['estimate'])} (permutation p {_p(nn['p_perm'])}), "
          + ", ".join(f"{t['metric']} {_f(t['estimate'])}" for t in ols if t["metric"] not in ("intercept", "new-new")) + ".", "",
          "![[h1_fold_enrichment.svg]]", ""]
    L += tech_section(S, T, main)
    L += ["## Sensitivities (all cohort relations)", "",
          "| Scenario | Pair type | Pairs | Median fold | Pooled fold | Pooled top-N fold | Median concordance (pairs) | Median ρ, any table (pairs) | Significant |",
          "|---|---|---|---|---|---|---|---|---|"]
    for scen, split, lab in (("published, brain", 2020, "Published brain, old up to 2020"),
                             ("published, brain, FDR lists", 2015, "Published brain, FDR lists, old up to 2015"),
                             ("all provenances, brain", 2015, "All provenances brain, old up to 2015"),
                             ("published, blood", 2015, "Published blood, old up to 2015")):
        first = True
        for t in ("old-old", "new-new", "old-new"):
            x = S(scen, split, t, "all")
            if not x:
                continue
            L.append(f"| {lab if first else ''} | {t} | {x['n_pairs']} | {_f(x['median_fold'])} | {_f(x['pooled_fold'])} | "
                     f"{_f(x['pooled_top_fold'])} | {_f(x['median_concordance'])} ({x['n_concordance']}) | "
                     + ("none" if not x["median_rho_any"] else f"{_f(x['median_rho_any'], 3)} ({x['n_rho_any']})") + f" | {x['n_sig_q05']} |")
            first = False
    L += ["", "Permutation p for new-new minus old-old (`h1_tests.tsv`), all pairs / different-cohort pairs:", ""]
    for scen, split in (("published, brain", 2020), ("published, brain, FDR lists", 2015), ("all provenances, brain", 2015), ("published, blood", 2015)):
        parts = []
        for m, lab in (("median log2 fold", "median log2 fold"), ("pooled fold (sum overlap / sum expected)", "pooled fold"),
                       ("pooled top-N fold", "pooled top-N fold"), ("median concordance", "concordance")):
            a, d = T(scen, split, "all", m), T(scen, split, "different", m)
            if a:
                parts.append(f"{lab} {_f(a['new_new'])} vs {_f(a['old_old'])}, p {_p(a['p_perm'])}" + (f" / {_p(d['p_perm'])}" if d else ""))
        if parts:
            L.append(f"- {scen}, split {split}: " + "; ".join(parts) + ".")
    # rho / pi1
    rp = [p for p in main if p.get("rho_pub")]
    top = sorted(rp, key=lambda p: -float(p["rho_pub"]))[:3]
    pi = [p for p in main if p.get("pi1_ab_pub") and p.get("pi1_base_ab_pub")]
    above = sum(1 for p in pi if float(p["pi1_ab_pub"]) > float(p["pi1_base_ab_pub"]))
    med = lambda v: float(np.median(v)) if v else None
    nnd = [float(p["rho_pub"]) for p in rp if p["era_2015"] == "new-new" and p["cohort"] == "different"]
    nns = [float(p["rho_pub"]) for p in rp if p["era_2015"] == "new-new" and p["cohort"] == "same"]
    L += ["", "## Correlation of logFC and π1", "",
          f"Published-table ρ exists for {len(rp)} main-scenario pairs ({sum(1 for p in rp if p['era_2015'] == 'old-old')} old-old); "
          f"median {_f(med([float(p['rho_pub']) for p in rp]), 3)}. New-new different-cohort pairs: "
          + (f"median {_f(med(nnd), 3)} ({len(nnd)} pairs)" if nnd else "no pair has published tables on both sides")
          + f"; new-new same-cohort pairs: median {_f(med(nns), 3)} ({len(nns)}). Largest: "
          + "; ".join(f"{p['label_a']} ({p['table_a_pub']}) with {p['label_b']} ({p['table_b_pub']}), {_f(p['rho_pub'])}, {p['cohort']} cohort, "
                      f"tables same region: {p['tables_same_region_pub']}" for p in top)
          + ". Hoffman et al. report a logFC Spearman ρ of 0.343 between two CommonMind cohorts "
          "([[Hoffman 2021 replication statistics — pi1, lfsr, and module preservation]]).", "",
          f"π1 of the later study's p-values over the earlier study's list, published tables ({len(pi)} pairs): median "
          f"{_f(med([float(p['pi1_ab_pub']) for p in pi]))} against a median all-gene baseline of {_f(med([float(p['pi1_base_ab_pub']) for p in pi]))}; "
          f"{above} of {len(pi)} above their baseline." if pi else "No published-table π1 in the main scenario.", ""]
    low = [p for p in main if p["concordance"] and int(p["n_dir"]) >= MIN_DIR and float(p["concordance"]) < 0.5]
    cnt = {}
    for p in low:
        for lab in (p["label_a"], p["label_b"]):
            cnt[lab] = cnt.get(lab, 0) + 1
    L += ["## Direction checks", "",
          f"{len(low)} main-scenario pairs with at least {MIN_DIR} shared directional genes have concordance below 0.5. Studies in "
          "three or more of them: " + (", ".join(f"{k} ({v})" for k, v in sorted(cnt.items(), key=lambda x: -x[1]) if v >= 3) or "none")
          + ". Each list's direction handling in its manifest entry should be checked before these are read as biology.", "",
          "## Results in brief", ""]
    oo, nw = S("published, brain", 2015, "old-old", "all"), S("published, brain", 2015, "new-new", "all")
    L.append(f"Published brain lists, era split at 2015: median fold enrichment {_f(oo['median_fold'])} for old-old pairs "
             f"({oo['n_pairs']}; null median {_f(oo['null_median_fold'])}) and {_f(nw['median_fold'])} for new-new pairs "
             f"({nw['n_pairs']}; null median {_f(nw['null_median_fold'])}); pooled fold {_f(oo['pooled_fold'])} vs "
             f"{_f(nw['pooled_fold'])}; size-matched top-N pooled fold {_f(oo['pooled_top_fold'])} vs {_f(nw['pooled_top_fold'])}; median "
             f"directional concordance {_f(oo['median_concordance'])} vs {_f(nw['median_concordance'])}; {oo['n_sig_q05']} of {oo['n_pairs']} "
             f"old-old and {nw['n_sig_q05']} of {nw['n_pairs']} new-new pairs overlap significantly (q < 0.05). Study-level permutation p: "
             + ", ".join(f"{lab} {_p(T('published, brain', 2015, 'all', m)['p_perm'])}" for m, lab in (
                 ("median log2 fold", "median log2 fold"), ("pooled fold (sum overlap / sum expected)", "pooled fold"),
                 ("pooled top-N fold", "pooled top-N fold"), ("median concordance", "concordance")))
             + "; different-cohort pairs only: "
             + ", ".join(f"{lab} {_p(T('published, brain', 2015, 'different', m)['p_perm'])}" for m, lab in (
                 ("median log2 fold", "median log2 fold"), ("pooled fold (sum overlap / sum expected)", "pooled fold"),
                 ("pooled top-N fold", "pooled top-N fold"), ("median concordance", "concordance")))
             + f". Same-cohort pairs are {S('published, brain', 2015, 'new-new', 'same')['n_pairs']} of {nw['n_pairs']} new-new pairs and "
             f"{S('published, brain', 2015, 'old-old', 'same')['n_pairs']} of {oo['n_pairs']} old-old pairs. With list sizes and cohort relation "
             f"in the model the new-new term is {_f(nn['estimate'])} log2 units (permutation p {_p(nn['p_perm'])}).")
    L.append("")
    with open(os.path.join(HERE, "h1_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    if "--report" in sys.argv:
        write_report()
    else:
        main()
        write_report()
