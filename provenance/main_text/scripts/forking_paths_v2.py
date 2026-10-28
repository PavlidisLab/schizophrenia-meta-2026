#!/usr/bin/env python3
"""Redesigned forking-paths grid (2026-10-10; block 3 added 2026-10-11), 282 designed paths (312 with FP_ERROR_STUDY_CHOICE=1). Engine = forking_paths_rerun_2026-10-09b_sn_bulk_only/
forking_paths.py (recurrence.py's own code: load_units, load_lists, representatives, study_tier, cohort_members,
cohort_tier; run directories as primary_rule_run.py builds them), with these changes:

  dropped    rule R2 (rules are R1 and A), probe handling (fixed at 'printed', the 'today' setting), the one-at-a-time,
             tier, subset, add and clone designs.
  threshold  fdr | nominal | all ('any criterion'); the empty-list choice (keep | drop) varies only under fdr; elsewhere it
             is fixed at drop (identical outcomes outside fdr in the 2026-10-09b grid)
  pseudo_noexpr  drop | keep | gtex   (gtex = HGNC pseudogenes dropped, and genes without median TPM >= 1 in at least one of
                 the 13 GTEx v8 brain regions dropped; inputs/gtex_brain_expressed_tpm1.tsv)
  noncoding  drop | pooled;  single_cell  without | with
  cut        published | tight | fc026 | fc058 | cap200      (list rows kept per study, after the threshold tier and gene filters)
               tight   per list row: reported FDR column < 0.01 where the row has one, else reported p < 0.001, else the
                       row is kept; a study whose lists were replaced by computed-FDR tables (threshold fdr) uses the
                       computed BH q (q_bh) < 0.01. A gene stays if any of its rows passes.
               fc026 / fc058  |log2FC| >= 0.26 / 0.58 from effect_sizes/data/viewA.tsv.gz (largest |log2FC| over the
                       study's rows for the gene); a gene with no converted effect is kept; a study with no converted
                       effect at all keeps its list.
               cap200  the 200 genes with the smallest reported p (else smallest reported FDR when every gene lacks p,
                       computed p for computed-FDR lists); ties broken by symbol; a study without that column, or with
                       <= 200 genes, keeps its list.
  study_filter  all | n20 (total N >= 20, h2_studies.tsv n) | from2016 (year >= 2016) | cortex (every primary list's region
                in prefrontal, cingulate, temporal, visual cortex; h1.region) | no_severe_error (not among the 14
                severity-1 studies, inputs/severe_error_pmids.tsv) | covariate_adjusted (final_call 'adjusted' for >= 1 of
                13 covariates, inputs/covariate_adjusted_pmids.tsv). A study with no metadata value fails a filter that needs
                it. The filter is applied where recurrence.py picks a group's representative, as the subset fork was.

Block 1 (96): rule 2 x threshold/empty 4 (fdr keep, fdr drop, nominal, all) x pseudo_noexpr 3 x noncoding 2 x single_cell 2,
              cut published, study_filter all.
Block 2 (150; 180 with FP_ERROR_STUDY_CHOICE=1): rule 2 x threshold 3 x cut 5 x study_filter 5 (6 with the no_severe_error filter) at empty drop, pseudo_noexpr keep, noncoding pooled,
              single_cell with (the 'today' settings).
Block 3 (36, 2026-10-11): rule 2 x threshold 3 x list_level 6 (non-default levels, see LIST_LEVELS) at the 'today' settings.
Hashing: before any model run, each path's post-fork inputs are hashed (selected study ids, per-study gene sets restricted to the
universe, universes, n_k, cluster structure); only distinct hashes are run (paths_hash.tsv).

    FP_GL=rerun-full/analyses python3 forking_paths_v2.py [--prepass-only]
"""
import argparse
import collections
import csv
import datetime
import functools
import gzip
import hashlib
import itertools
import json
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.environ.get("FP_GL", "rerun-full/analyses")
REPO = os.path.dirname(GL)
sys.path[:0] = [GL, os.path.join(GL, "primary_rule")]
import primary_rule as P  # noqa: E402
import primary_rule_run as PR  # noqa: E402

MODS = P.load_code(GL)
R, H = MODS["recurrence"], MODS["h1"]
TODAY = datetime.date.today().isoformat()
RULES = ("R1", "A")
PATH_FORKS = ("rule", "threshold", "empty", "pseudo_noexpr", "noncoding", "single_cell", "cut", "study_filter", "list_level")
# list_level (2026-10-11, lists-fix flags in analyses/lists_fix_2026-10-10/subset_flags.tsv): how complete a study's primary list must be.
#   inclusive      default: every study as in primary_lists.json
#   no_partial     drop the partial DE lists: Scarr 2018 (29463818), Arion 2015 (25560755), Altar 2005 (16038679), Shao 2008 (18191109)
#   no_nonDE       drop the non-DE-style lists: Ling 2024 (38448582), Kohen 2014 (24594777)
#   no_mixed       drop the mixed-diagnosis lists: Pai 2019 (31053723), Frohlich 2024 (39227716)
#   chang_s2only   Chang 2017 (28809853) from its Table S2 any-method list only (default: Table S2 union + Table S4)
#   pai_text       Pai 2019's list replaced by the three genes named in the text (NR4A1, NR4A2, EGR1; all down, Data S5 logFC
#                  -1.43, -1.37, -1.02; FDR 0.123, 0.060, 0.123) instead of the Data S5 rows at FDR < 0.05
#   no_flagged     no_partial + no_nonDE + no_mixed together
LIST_LEVELS = ("inclusive", "no_partial", "no_nonDE", "no_mixed", "chang_s2only", "pai_text", "no_flagged")
LEVEL_DROP = {"no_partial": {"29463818", "25560755", "16038679", "18191109"}, "no_nonDE": {"38448582", "24594777"},
              "no_mixed": {"31053723", "39227716"}}
LEVEL_DROP["no_flagged"] = set().union(*LEVEL_DROP.values())
CHANG, CHANG_S2 = "28809853", "28809853__table_s2_any_method"
PAI, PAI_TEXT_GENES = "31053723", {"NR4A1": "down", "NR4A2": "down", "EGR1": "down"}
CUTS_ = ("published", "tight", "fc026", "fc058", "cap200")
# FP_ERROR_STUDY_CHOICE=1 restores the earlier design with the no_severe_error study filter (312 designed paths: 96 + 180 + 36);
# the default (unset or 0) leaves that choice out (282 designed paths: 96 + 150 + 36).
ERROR_STUDY_CHOICE = os.environ.get("FP_ERROR_STUDY_CHOICE", "0") == "1"
FILTERS = ("all", "n20", "from2016", "cortex") + (("no_severe_error",) if ERROR_STUDY_CHOICE else ()) + ("covariate_adjusted",)
CUTS = (3, 4, 5)
NULLS = [("uniform", "n/a", "ours", "tail", "pval"), ("expr", "n/a", "ours", "tailE", "pexpr"),
         ("de_prior", "ours", "ours", "tailP", "pprior"), ("de_prior", "global", "global", "tailP", "pprior")]
NULLS_C = [("uniform", "n/a", "ours", "tailC", "pC"), ("expr", "n/a", "ours", "tailCE", "pCE"),
           ("de_prior", "ours", "ours", "tailCP", "pCP"), ("de_prior", "global", "global", "tailCP", "pCP")]
PRIORS = ("ours", "global")
SPECS = P.manifest_specs(GL)
CORTEX = ("prefrontal", "cingulate", "temporal", "visual cortex")
GTEX_FILE = os.path.join(HERE, "inputs", "gtex_brain_expressed_tpm1.tsv")
DEFAULT_EXPRESSED = R.EXPRESSED_FILE
VIEWA = os.path.join(GL, "effect_sizes", "data", "viewA.tsv.gz")
H2 = os.path.join(GL, "h2", "h2_studies.tsv")


def sha(p, n=64):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:n]


def rd(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


# ---------------------------------------------------------------- memoised file readers (same values, read once)
def _memo(name):
    f = getattr(R, name)
    c = functools.lru_cache(maxsize=None)(f)
    setattr(R, name, c)


for _n in ("read_list", "read_set", "computed_fdr_lists"):
    _memo(_n)
_gfs = R.gene_filter_set


@functools.lru_cache(maxsize=None)
def _gfs_cached(pn, nc, ef):
    return _gfs()


R.gene_filter_set = lambda: _gfs_cached(R.PSEUDO_NOEXPR, R.NONCODING, R.EXPRESSED_FILE)


# ---------------------------------------------------------------- study metadata for the filters
def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


H2ROWS = {r["pmid"]: r for r in rd(H2)}
SEVERE = {r["pmid"] for r in rd(os.path.join(HERE, "inputs", "severe_error_pmids.tsv"))}
COVADJ = {r["pmid"]: int(r["n_adjusted"]) for r in rd(os.path.join(HERE, "inputs", "covariate_adjusted_pmids.tsv"))}
REGION_CACHE = {}


def regions_ok(primary):
    key = tuple(primary)
    if key not in REGION_CACHE:
        REGION_CACHE[key] = all(H.region(SPECS, n) in CORTEX for n in primary)
    return REGION_CACHE[key]


def allow_fn(filt, studies):
    """k -> bool for study_filter `filt`; a missing metadata value fails the filter."""
    if filt == "all":
        return lambda k: True
    if filt == "n20":
        return lambda k: (num(H2ROWS.get(k, {}).get("n")) or 0) >= 20
    if filt == "from2016":
        return lambda k: (num(H2ROWS.get(k, {}).get("year")) or 0) >= 2016
    if filt == "cortex":
        return lambda k: k in studies and regions_ok(studies[k]["primary"])
    if filt == "no_severe_error":
        return lambda k: k not in SEVERE
    if filt == "covariate_adjusted":
        return lambda k: COVADJ.get(k, 0) > 0
    raise ValueError(filt)


# ---------------------------------------------------------------- cut tables
class CutTables:
    """Per list name: gene -> {rows: [(p, fdr)]}; computed-FDR tables: gene -> [(p, q_bh)] ; viewA: study -> gene -> max |log2FC|."""

    def __init__(self):
        self.lists = {}
        self.cf = None
        self.fc = collections.defaultdict(dict)
        with gzip.open(VIEWA, "rt") as fh:
            next(fh)
            for l in fh:
                f = l.rstrip("\n").split("\t")
                v = num(f[4])
                if v is None:
                    continue
                d = self.fc[f[0]]
                d[f[3]] = max(d.get(f[3], 0.0), abs(v))

    def list_rows(self, name):
        if name not in self.lists:
            d = collections.defaultdict(list)
            rows = [l for l in open(os.path.join(R.DATA, "lists", name + ".tsv"), encoding="utf-8") if not l.startswith("#")]
            for r in csv.DictReader(rows, delimiter="\t"):
                if r.get("symbol"):
                    d[r["symbol"]].append((num(r.get("p")), num(r.get("fdr"))))
            self.lists[name] = d
        return self.lists[name]

    def cf_rows(self, name):
        if self.cf is None:
            self.cf = {}
            lines = [l for l in open(os.path.join(R.COMPUTED_FDR_DIR, "summary.tsv"), encoding="utf-8") if not l.startswith("#")]
            for r in csv.DictReader(lines, delimiter="\t"):
                if r["route"] not in R.COMPUTED_FDR_ROUTES:
                    continue
                d = collections.defaultdict(list)
                rows = [l for l in open(os.path.join(R.COMPUTED_FDR_DIR, r["output"]), encoding="utf-8") if not l.startswith("#")]
                for x in csv.DictReader(rows, delimiter="\t"):
                    if x["symbol"]:
                        d[x["symbol"]].append((num(x.get("p")), num(x.get("q_bh"))))
                for n in (y.strip() for y in r["lists"].split(",")):
                    for g, v in d.items():
                        self.cf.setdefault(n, collections.defaultdict(list))[g].extend(v)
        return self.cf[name]


CT = None
STATS = collections.Counter()


def apply_cut(cut, units, lists, report):
    """Cut the studies' lists in place; report[k] = 'changed' | 'unchanged'|'n/a' per study (for the log)."""
    if cut == "published":
        return
    global CT
    if CT is None:
        CT = CutTables()
    studies = units["studies"]
    for k in list(lists):
        L = lists[k]
        if not L:
            continue
        cfn = units["computed_fdr"].get(k)
        names = cfn if cfn else studies[k]["primary"]
        rows = collections.defaultdict(list)  # gene -> [(p, fdr_or_q)]
        for n in names:
            src = CT.cf_rows(n) if cfn else CT.list_rows(n)
            for g in L:
                if g in src:
                    rows[g].extend(src[g])
        keep = None
        if cut == "tight":
            if cfn:
                keep = {g for g in L if not rows.get(g) or any(q is None or q < 0.01 for _, q in rows[g])}
            else:
                def ok(p, f):
                    if f is not None:
                        return f < 0.01
                    if p is not None:
                        return p < 0.001
                    return True
                keep = {g for g in L if not rows.get(g) or any(ok(p, f) for p, f in rows[g])}
        elif cut in ("fc026", "fc058"):
            thr = 0.26 if cut == "fc026" else 0.58
            fc = CT.fc.get(k)
            if fc:
                keep = {g for g in L if g not in fc or fc[g] >= thr}
        elif cut == "cap200":
            if len(L) > 200:
                key = None
                for ix in (0, 1):
                    vals = {}
                    for g in L:
                        v = [r[ix] for r in rows.get(g, []) if r[ix] is not None]
                        if v:
                            vals[g] = min(v)
                    if len(vals) == len(L):
                        key = vals
                        break
                if key is not None:
                    order = sorted(L, key=lambda g: (key[g], g))
                    keep = set(order[:200])
                    if key[order[199]] == key[order[200]]:
                        report["_ties"] = report.get("_ties", 0) + 1
        if keep is None:
            report[k] = "n/a"
            continue
        report[k] = "changed" if len(keep) != len(L) else "unchanged"
        lists[k] = {g: v for g, v in L.items() if g in keep}


# ---------------------------------------------------------------- one path
class Ctx:
    def __init__(self):
        self.dirs = {r: PR.make_run_dir(r, GL) for r in RULES}
        PR.rebase(MODS, self.dirs["A"], GL)
        self.nodes = json.load(open(R.LINEAGE))["nodes"]
        self.hitn = {}
        keep = R.DE_PRIOR_SOURCE, R.DE_PRIOR
        for src in PRIORS:
            R.DE_PRIOR_SOURCE, R.DE_PRIOR = src, R.DE_PRIOR_SOURCES[src]["file"]
            self.hitn[src] = R.de_prior_weights()
        R.DE_PRIOR_SOURCE, R.DE_PRIOR = keep
        self.expr = R.expression_levels()

    def close(self):
        for d in self.dirs.values():
            shutil.rmtree(d, ignore_errors=True)


def set_switches(cfg):
    R.THRESHOLD_TIER = cfg["threshold"]
    R.EMPTY_LISTS = cfg["empty"]
    R.PSEUDO_NOEXPR = "keep" if cfg["pseudo_noexpr"] == "keep" else "drop"
    R.EXPRESSED_FILE = GTEX_FILE if cfg["pseudo_noexpr"] == "gtex" else DEFAULT_EXPRESSED
    R.NONCODING = cfg["noncoding"]
    R.PRINTED_NAME_FALLBACK = "keep"   # probe handling fixed at 'printed' (the 'today' setting; inert in the 2026-10-09b grid)
    R.ADD_STUDIES = ""


def prepare(cfg, ctx):
    """Everything downstream of the forks up to the model: returns dict with units, lists, M, nk, sel and the input hash."""
    d = ctx.dirs[cfg["rule"]]
    PR.rebase(MODS, d, GL)
    pl = json.load(open(os.path.join(d, "primary_lists.json")))
    uni = json.load(open(os.path.join(d, "universes.json")))["lists"]
    if cfg["single_cell"] == "without":
        pl["rule"]["modalities"] = [m for m in pl["rule"]["modalities"] if m != "other statistic"]
        for k in set(json.load(open(os.path.join(HERE, "sn_population.json")))["studies"]) | {"38448582"}:
            pl["studies"].pop(k, None)
    else:
        sn = json.load(open(os.path.join(HERE, "sn_population.json")))
        pl["studies"].update(sn["studies"])
        uni.update(sn["universes"])
    if cfg["list_level"] == "chang_s2only":
        pl["studies"][CHANG]["primary"] = [CHANG_S2]
    set_switches(cfg)
    units = R.load_units(pl, ctx.nodes)
    lists, M, nk, utype, nouni = R.load_lists(units, uni)
    if cfg["list_level"] == "pai_text" and PAI in lists:
        lists[PAI] = {g: v for g, v in PAI_TEXT_GENES.items() if g in M.get(PAI, ())}
    cutrep = {}
    apply_cut(cfg["cut"], units, lists, cutrep)
    studies = units["studies"]
    base = lambda s: s["provenance"] == "published" and s["tissue"] == "brain"   # SUBSET['brain'] of the earlier grids
    allow0 = allow_fn(cfg["study_filter"], studies)
    drop = LEVEL_DROP.get(cfg["list_level"], ())
    allow = lambda k: allow0(k) and k not in drop
    usable = R.usable(units, lists)
    reps = R.representatives(units, base, lambda k: usable(k) and allow(k))
    sel = [k for k in sorted(reps, key=lambda k: (0, int(k), "") if k.isdigit() else (1, 0, k)) if k in M]
    # hash of the post-fork inputs: selected ids, per-study list genes within the universe, universe genes and n_k, clusters
    h = hashlib.sha256()
    h.update(("|".join(sel)).encode())
    for k in sel:
        h.update(f"\n{k}\t{nk[k]:.9g}\t".encode())
        h.update(hashlib.sha256("\n".join(sorted(set(lists[k]) & M[k])).encode()).digest())
        h.update(hashlib.sha256("\n".join(sorted(M[k])).encode()).digest())
        h.update(f"\t{R.cluster_of(units, k)}\t{units['multibank'].get(k)}\t{units['bank_lists'].get(k)}".encode())
    # clusters of unselected group members cannot enter the model; cname/gcluster enter only through cluster_of above
    return {"units": units, "lists": lists, "M": M, "nk": nk, "sel": sel, "hash": h.hexdigest(), "cutrep": cutrep}


def compute(prep, ctx):
    units, lists, M, nk, sel = prep["units"], prep["lists"], prep["M"], prep["nk"], prep["sel"]
    # multi-bank studies with per-bank lists are read back from the list files by recurrence.cohort_members: not cut
    if any(k in units["multibank"] and units["bank_lists"].get(k) for k in sel):
        STATS["multibank_bank_lists_in_sel"] += 1
    res = {}
    for src in PRIORS:
        st_ = R.study_tier(sel, M, nk, lists, ctx.hitn[src], ctx.expr)
        clusters, mb_out = R.cohort_members(sel, units, st_, M, nk)
        res[src] = (st_, R.cohort_tier(clusters, st_["genes"], st_["gi"]))
    return res


def make_info(prep, res):
    units, sel = prep["units"], prep["sel"]
    st_, ct_ = res["ours"]
    st = units["studies"]
    clusters_n = sum(1 for _ in ct_["cl"])
    return {"n_studies": len(sel), "n_empty": sum(1 for k in sel if not st_["Lr"][k]),
            "n_clusters": clusters_n, "n_genes_measured": len(st_["genes"]),
            "median_list_size": int(np.median([len(st_["Lr"][k]) for k in sel])) if sel else 0,
            "studies": ";".join(sel), "empty_studies": ";".join(k for k in sel if not st_["Lr"][k])}


def summarize(cfg, info, res):
    genes = np.array(res["ours"][0]["genes"])
    out, sets = [], {}
    for tier, j, obs_key, nulls in (("study", 0, "obs", NULLS), ("cohort", 1, "obsC", NULLS_C)):
        obs = res["ours"][j][obs_key]
        gs = {}
        for n, src, run, _, pk in nulls:
            p = res[run][j][pk]
            q = R.bh(p) if len(p) else np.array([])
            gs[n if n != "de_prior" or src == "ours" else "de_prior_global"] = set(genes[q < 0.1])
        sets.update({(tier, n): g for n, g in gs.items()})
        for n, src, run, tk, _ in nulls:
            tail = res[run][j][tk]
            row = dict(cfg, unit_tier=tier, null=n, de_prior_source=src,
                       **{k: v for k, v in info.items() if k not in ("studies", "empty_studies")})
            for c in CUTS:
                o = int((obs >= c).sum())
                e = float(tail[:, c].sum()) if c < tail.shape[1] else 0.0
                row[f"obs_ge{c}"], row[f"exp_ge{c}"] = o, round(e, 1)
                row[f"ratio_ge{c}"] = round(o / e, 3) if e > 0 else ""
            row["genes_q10_uniform"], row["genes_q10_expr"] = len(gs["uniform"]), len(gs["expr"])
            row["genes_q10_prior"] = len(gs["de_prior" if src != "global" else "de_prior_global"])
            row["headline_genes"] = "" if n == "de_prior" else len(gs[n])
            out.append(row)
    return out, sets


def path_key(cfg):
    return "|".join(cfg[f] for f in PATH_FORKS)


# ---------------------------------------------------------------- designs
def designs():
    out = []
    today = dict(rule="A", threshold="all", empty="drop", pseudo_noexpr="keep", noncoding="pooled", single_cell="with",
                 cut="published", study_filter="all", list_level="inclusive")
    for rule, (thr, emp), pn, nc, sc in itertools.product(RULES, (("fdr", "keep"), ("fdr", "drop"), ("nominal", "drop"), ("all", "drop")),
                                                        ("drop", "keep", "gtex"), ("drop", "pooled"), ("without", "with")):
        out.append(("block1", dict(today, rule=rule, threshold=thr, empty=emp, pseudo_noexpr=pn, noncoding=nc, single_cell=sc)))
    for rule, thr, cut, filt in itertools.product(RULES, ("fdr", "nominal", "all"), CUTS_, FILTERS):
        out.append(("block2", dict(today, rule=rule, threshold=thr, cut=cut, study_filter=filt)))
    for rule, thr, lvl in itertools.product(RULES, ("fdr", "nominal", "all"), LIST_LEVELS[1:]):
        out.append(("block3", dict(today, rule=rule, threshold=thr, list_level=lvl)))
    return out


def header():
    files = [os.path.join(GL, "primary_rule", r, x) for r in RULES for x in ("primary_lists.json", "universes.json")]
    files += [os.path.join(GL, "recurrence.py"), os.path.join(GL, "extract.py"), __file__, R.LINEAGE, R.DE_PRIOR, R.EXPR_FILE,
              os.path.join(R.COMPUTED_FDR_DIR, "summary.tsv"), os.path.join(GL, "technology.json"),
              os.path.join(GL, "platform_type", "platform_type.json"), R.EXPRESSED_FILE, VIEWA, H2,
              os.path.join(HERE, "sn_population.json"),
              os.path.join(REPO, "data", "reference", "hgnc_complete_set.txt")]
    files += [os.path.join(HERE, "inputs", f) for f in ("severe_error_pmids.tsv", "covariate_adjusted_pmids.tsv", "gtex_brain_expressed_tpm1.tsv")]
    lists = hashlib.sha256("".join(f"{n}{sha(os.path.join(REPO, 'data', 'lists', n))}" for n in
                                   sorted(os.listdir(os.path.join(REPO, "data", "lists")))).encode()).hexdigest()
    return ([f"# script: analyses/multiverse_analysis/forking_paths_v2.py, run {TODAY}"]
            + [f"# input: {os.path.relpath(f, REPO) if f.startswith(REPO) else f} sha256 {sha(f)}" for f in files]
            + [f"# input: data/lists/*.tsv (sha256 of name + sha256 per file, sorted) sha256 {lists}"])


def write(path, rows, fields, head):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(prepass_only=False):
    ctx = Ctx()
    t0 = time.time()
    try:
        # ---- pre-pass: hash every path
        preps, order = {}, []
        for block, cfg in designs():
            key = path_key(cfg)
            t = time.time()
            prep = prepare(cfg, ctx)
            st = prep["units"]["studies"]
            pinfo = {"n_list_files": sum(len(prep["units"]["computed_fdr"].get(k) or st[k]["primary"]) for k in prep["sel"]),
                     "n_computed_fdr": sum(1 for k in prep["sel"] if k in prep["units"]["computed_fdr"]),
                     "cut_changed_studies": sum(1 for k in prep["sel"] if prep["cutrep"].get(k) == "changed"),
                     "cut_na_studies": sum(1 for k in prep["sel"] if prep["cutrep"].get(k) == "n/a"),
                     "cap200_ties_at_boundary": prep["cutrep"].get("_ties", 0)}
            order.append((block, key, cfg, prep["hash"], prep["sel"], pinfo))
            preps.setdefault(prep["hash"], (cfg, prep))
            print(f"prepass {len(order):3d} {block} {key} sel={len(prep['sel'])} {prep['hash'][:10]} {time.time()-t:.1f}s", flush=True)
            if not preps[prep["hash"]][1] is prep:
                prep.clear()   # free memory of repeated hashes (info is rebuilt below only from the first)
        distinct = len(preps)
        print(f"{len(order)} designed paths, {len({o[1] for o in order})} distinct specifications, {distinct} distinct hashes "
              f"({time.time()-t0:.0f}s)", flush=True)
        head = header()
        hash_rows = [{"path": key, "block": block, "hash": hsh, "n_studies": len(sel), "first_path_with_hash": preps[hsh][0] and path_key(preps[hsh][0])}
                     for block, key, cfg, hsh, sel, _p in order]
        write(os.path.join(HERE, "paths_hash.tsv"), hash_rows, ["path", "block", "hash", "n_studies", "first_path_with_hash"],
              head + [f"# {len(order)} designed paths ({'; '.join(f'{b} {sum(1 for o in order if o[0] == b)}' for b in ('block1', 'block2', 'block3'))}; error-study choice {'included' if ERROR_STUDY_CHOICE else 'left out'}); {len({o[1] for o in order})} distinct specifications; "
                      f"{distinct} distinct hashes of the post-fork inputs (selected study ids, per-study gene sets within the universe, "
                      "universe genes, n_k, cluster of each study); only distinct hashes were run",
                      "# path = " + "|".join(PATH_FORKS)])
        if prepass_only:
            return
        # ---- model runs, one per distinct hash
        done, runtime = {}, []
        for hsh, (cfg, prep) in preps.items():
            t = time.time()
            res = compute(prep, ctx)
            info = make_info(prep, res)
            done[hsh] = (res, info)
            runtime.append({"hash": hsh, "path": path_key(cfg), "seconds": round(time.time() - t, 2)})
            print(f"run {len(done):3d}/{distinct} {path_key(cfg)} {runtime[-1]['seconds']}s", flush=True)
        # ---- rows per designed path (the path's own cfg and info, the hash's results)
        rows, gene_rows, study_rows = [], [], []
        summ = {}
        for block, key, cfg, hsh, sel, pinfo in order:
            res, info0 = done[hsh]
            if hsh not in summ:
                summ[hsh] = summarize(cfg, info0, res)
            prow, sets = summ[hsh]
            # n_list_files / n_computed_fdr are properties of the path's own units, rebuilt for the path
            rows += [{**r, **cfg, **pinfo, "block": block, "path": key, "hash": hsh} for r in prow]
            for (tier, n), gs in sorted(sets.items()):
                if n.startswith("de_prior"):
                    continue
                for g in sorted(gs):
                    gene_rows.append({"path": key, "block": block, "unit_tier": tier, "null": n, "gene": g})
            study_rows.append({"path": key, "block": block, "hash": hsh, "studies": info0["studies"], "empty_studies": info0["empty_studies"]})
    finally:
        ctx.close()
    fields = ["block", "path", "hash"] + list(PATH_FORKS) + ["unit_tier", "null", "de_prior_source", "n_studies", "n_list_files", "n_empty",
              "n_computed_fdr", "n_clusters", "n_genes_measured", "median_list_size", "cut_changed_studies", "cut_na_studies", "cap200_ties_at_boundary"]
    fields += [f"{p}_ge{c}" for c in CUTS for p in ("obs", "exp", "ratio")]
    fields += ["headline_genes", "genes_q10_uniform", "genes_q10_expr", "genes_q10_prior"]
    write(os.path.join(HERE, "grid.tsv"), rows, fields, head + [
        "# one row per designed path x unit tier x null (x DE-prior source); path = " + "|".join(PATH_FORKS),
        "# a path's outcome columns equal those of every path with the same hash; n_list_files, n_computed_fdr and cut_* are the path's own",
        "# obs_geC: genes in >= C lists (study tier) or clusters (cohort tier); exp_geC: expected under the row's null; ratio = obs / exp; "
        "headline_genes: genes at BH q < 0.1 under the row's null (blank for de_prior)"])
    write(os.path.join(HERE, "grid_genes.tsv"), gene_rows, ["path", "block", "unit_tier", "null", "gene"], head + [
        "# genes at BH q < 0.1 per designed path, unit tier and null (uniform, expr); join to grid.tsv on block + path"])
    write(os.path.join(HERE, "grid_studies.tsv"), study_rows, ["path", "block", "hash", "studies", "empty_studies"], head + [
        "# the study-tier representatives of each path (PMIDs) and those whose list is empty"])
    write(os.path.join(HERE, "runtime.tsv"), runtime, ["hash", "path", "seconds"], head + [
        f"# total {round(time.time() - t0)} s including the pre-pass; {distinct} model runs; " + json.dumps(dict(STATS))])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepass-only", action="store_true")
    a = ap.parse_args()
    main(a.prepass_only)
