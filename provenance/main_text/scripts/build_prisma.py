#!/usr/bin/env python3
"""Project-wide search record (PRISMA 2020 / PRISMA-S): every count computed from the files of record.

Writes, in search/:
  prisma_counts.json      every flow number, nested by stream (no timestamps, so a rerun on unchanged inputs is
                          byte-identical)
  build_prisma_meta.json  sha256 of every input file (directories: sha256 of the sorted file-name list), script sha256,
                          run time
  search_log.tsv          one row per search: search_id, date, stream, source, exact_query, filters, hits, after_dedup,
                          screened, included, file_of_record, then screening_method, decided_by, count_notes. A count
                          that could not be reconstructed from a file is "NR"; count_notes gives the value stated
                          elsewhere and where
  prisma_16b.tsv          PRISMA 2020 item 16b: every study excluded at the eligibility stage (after passing a model
                          screen) whose exclusion carries a recorded reason, with the reason verbatim
  criteria.md             the screening and eligibility criteria verbatim, read from the scripts and biolit configs
and rewrites the blocks of PRISMA.md between `<!-- prisma:begin NAME -->` and `<!-- prisma:end NAME -->`.

Nothing is fetched: search/rerun_pubmed_merikangas.py saves the one network-derived input (the PubMed rerun and its
esummary metadata). Inputs outside the repo: the vault's biolit-results/ and biolit-configs/ (VAULT_RA overrides the
vault Research Assistant path, as in tests/corpus/corpus_testlib.py).

Usage:  python3 search/build_prisma.py            (write everything)
        python3 search/build_prisma.py --check    (exit 1 if any output would change; writes nothing)
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import io
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
VAULT_RA = Path(os.environ.get("VAULT_RA", str(Path.home() / "Documents/PaulsObsidianVault/Research Assistant")))
BIOLIT = VAULT_RA / "biolit-results"
CONFIGS = VAULT_RA / "biolit-configs"
CORPUS = REPO / "corpus"
GL = REPO / "analyses"
PRISMA_MD = HERE / "PRISMA.md"

# the April 2026 biolit runs (tracker "Key artefacts"; same table as tests/corpus/corpus_testlib.py)
APRIL_RUNS = {
    "pubmed_m26": "merikangas-review/run_20260404_111222",
    "pubmed_new": "merikangas-rescreen/run_20260404_133102",
    "embase_1": "merikangas-embase/run_20260404_141521",
    "embase_2": "merikangas-embase/run_20260404_143045",
    "wos": "merikangas-wos/run_20260404_151130",
}
APRIL16_RUNS = {"bousman": "merikangas-bousman/run_20260416_175227", "embase_doi": "merikangas-embase-doi/run_20260416_180012"}
BACKWARD_SCREEN = ["merikangas-backward/run_20260922_165405", "merikangas-backward/run_20260922_184508"]
BACKWARD_PASS = ["merikangas-backward/run_20260922_190415/results_part1.csv", "merikangas-backward/run_20260922_184508/results_part2.csv"]
BROADER_RUNS = [f"merikangas-backward/run_20260922_20384{i}" for i in (2, 3, 4, 5)]
APRIL_DATE = "2026/04/04"

INPUTS: dict[str, str] = {}


# ---------------------------------------------------------------- input helpers (every read is recorded with its sha256)

def _reconstruction(sid):
    """Recall and precision of the reconstructed backward filter (backward_reconstruction/reconstruction_report.tsv)."""
    import csv as _csv
    with open(Path(__file__).resolve().parent / "backward_reconstruction" / "reconstruction_report.tsv") as fh:
        r = {x["search"]: x for x in _csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}[sid]
    return (f" Filter reconstructed 2026-10-08 (search/backward_reconstruction/rebuild_backward.py): recovers {r['recovered']} of "
            f"{r['list_of_record']} papers of record (recall {r['recall']}), {r['extra']} extra (precision {r['precision']}); the list of record is authoritative.")

def _rel(p: Path) -> str:
    p = Path(p)
    for base, tag in ((REPO, ""), (VAULT_RA, "vault:Research Assistant/")):
        try:
            return tag + str(p.resolve().relative_to(base.resolve()))
        except ValueError:
            pass
    return str(p)


def use(p) -> Path:
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(p)
    if p.is_dir():
        names = "\n".join(sorted(x.name for x in p.iterdir()))
        INPUTS[_rel(p) + "/"] = hashlib.sha256(names.encode()).hexdigest()
    else:
        INPUTS[_rel(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def rjson(p):
    return json.loads(use(p).read_text(encoding="utf-8"))


def rtext(p):
    return use(p).read_text(encoding="utf-8", errors="replace")


def rtsv(p):
    return list(csv.DictReader((l for l in rtext(p).splitlines(True) if not l.startswith("#")), delimiter="\t"))


def rcsv(p):
    return list(csv.DictReader(io.StringIO(rtext(p))))


def artifact_ids(run: str) -> set[str]:
    d = use(BIOLIT / run / "artifacts")
    return {x.name.split("_")[0] for x in d.iterdir()}


def keys_of(r: dict) -> set[str]:
    """identifiers of a record: PMID and lower-case DOI (as tests/corpus/corpus_testlib.record_keys)"""
    return {x for x in ((r.get("pmid") or "").strip(), (r.get("doi") or "").strip().lower()) if x and x not in ("null", "None")}


def group_records(records: list[dict], extra_keys=None) -> list[list[int]]:
    """indices grouped into papers: records sharing any key are one paper (union-find)"""
    parent = list(range(len(records)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    seen = {}
    for i, r in enumerate(records):
        ks = keys_of(r) | (extra_keys(r) if extra_keys else set())
        for k in ks:
            if k in seen:
                parent[find(i)] = find(seen[k])
            else:
                seen[k] = i
    g = collections.defaultdict(list)
    for i in range(len(records)):
        g[find(i)].append(i)
    return list(g.values())


def norm_title(t: str) -> str:
    t = re.sub(r"<[^>]+>|\{|\}", "", t or "").lower()
    return re.sub(r"[^a-z0-9]+", " ", t).strip()


def decider(note: str) -> str:
    n = note or ""
    if re.search(r"\[decision: Paul\]|per Paul|\(Paul \d{4}-|Paul: add", n):
        return "Paul"
    if re.search(r"\[Claude", n):
        return "Claude"
    return "Paul" if "Paul" in n else "not recorded"


REASON_RULES = [  # (category, pattern) applied in order to the verbatim reason; the first match wins
    ("not a SZ vs control expression study (stock note of the cited-source pass)", r"not a sz vs control expression study"),
    ("duplicate or preprint of an included study", r"duplicate|preprint"),
    ("data reuse / reanalysis (no new data)", r"re-?use|reanaly|re-analy|same .{0,40}(data|dataset|arrays|cohort|gels|comparison)|previously reported|data first reported|first reported in|probable reanalysis|data from another study|obtained the corresponding"),
    ("meta-analysis", r"meta-analy"),
    ("non-coding RNA only", r"non-coding|circ(ular)? ?rna|lncrna|mirna|microrna|small rna"),
    ("targeted, not transcriptome-wide", r"candidate|qpcr|q-pcr|rt-pcr|targeted|not genome|glycogene|single-gene|one gene|panel|\bsrm\b|in situ|not transcriptome-wide|not a discovery|pre-selected"),
    ("cell line, animal or model system (no patient vs control samples)", r"sh-sy5y|hek293|crispr|knockdown|overexpression|differentiation|cell-line|cell line|mutation trio|macaque|ntera2|\brats?\b|mouse|mice|animal|es line"),
    ("no schizophrenia vs control result", r"no sz|no schizophrenia|not sz|not a sz|bipolar|22q11|other disorder|within-patient|no case-control|not schizophrenia|non-sz|no .{0,20}contrast|not case vs control|not a case-control|methamphetamine|alzheimer|contrast is|twins"),
    ("protein, not RNA (before proteomics were admitted)", r"protein not rna"),
    ("not relevant (not expression, method or review)", r"atlas|review|gwas|eqtl|twas|genetic|method|not relevant|not a de|imaging|off-topic|methylation|cytokine|facs"),
]


def reason_category(text: str) -> str:
    t = (text or "").lower()
    if not t.strip():
        return "no reason recorded"
    for cat, pat in REASON_RULES:
        if re.search(pat, t):
            return cat
    return "other (see reason)"


# ---------------------------------------------------------------- the lineage (final corpus)
def lineage_nodes():
    return rjson(CORPUS / "lineage.json")["nodes"]


def corpus_studies(nodes):
    """as numbers_ledger.py _studies() and corpus_testlib.studies(): not outside sources, not the review landmark"""
    return [n for n in nodes if not n.get("external") and not n.get("landmark")]


def node_index(nodes):
    """lookup by id, by id without the 'sc:' prefix, by PMID and by lower-case DOI-shaped id"""
    idx = {}
    for n in nodes:
        for k in {n["id"], n["id"].removeprefix("sc:"), n.get("pmid") or "", (n.get("doi") or "").lower()}:
            if k:
                idx.setdefault(k.lower(), n)
    return idx


def find_node(idx, key: str):
    k = (key or "").strip().lower()
    if k.startswith("doi_"):
        k = k[4:].replace("_", "/", 1)
    return idx.get(k) or idx.get(k.removeprefix("sc:"))


# ---------------------------------------------------------------- stream A: April 2026 three-database replication
def parse_embase():
    recs = []
    for part in ("part1", "part2"):
        text = rtext(HERE / "materials" / f"merikangas-embase-citation.{part}.xml")
        for rec in re.split(r"<record index=", text)[1:]:
            def field(code):
                m = re.search(rf"<F C='{code}'.*?</F>", rec, re.S)
                return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", m.group(0))).strip() if m else ""
            recs.append({"pmid": (re.search(r"\d+", field("pm")) or [""])[0], "doi": field("do").lower().rstrip("."),
                         "title": field("ti").rstrip(". "), "abstract": field("ab")})
    return recs


def parse_wos():
    text = rtext(HERE / "materials" / "merikangas-webofscience.bib")
    recs = []
    for e in re.split(r"\n(?=@)", "\n" + text)[1:]:
        def field(name):
            m = re.search(rf"\n\s*{name} = \{{(.*?)\}},?\n\s*[A-Za-z-]+ = |\n\s*{name} = \{{(.*?)\}}\}}", e, re.S)
            return re.sub(r"\s+", " ", (m.group(1) or m.group(2))).strip() if m else ""
        recs.append({"pmid": "", "doi": field("DOI").lower(), "title": field("Title"), "abstract": field("Abstract")})
    return recs


def stream_april():
    out = {}
    rr_path = sorted((HERE / "reruns").glob("pubmed_merikangas_*.json"))[-1]
    rr = rjson(rr_path)
    es = rr["esummary"]
    rerun = set(rr["pmids_all_dates"])
    indexed_by_april = sorted((p for p in rerun if (es.get(p, {}).get("mesh_date") or "9999")[:10] <= APRIL_DATE), key=int)
    art = {k: artifact_ids(v) for k, v in {**APRIL_RUNS, **APRIL16_RUNS}.items()}
    passed = {k: rcsv(BIOLIT / v / "results.csv") for k, v in {**APRIL_RUNS, **APRIL16_RUNS}.items()}
    pm_screened = art["pubmed_m26"] | art["pubmed_new"]
    out["pubmed"] = {
        "query": rr["query"], "rerun_file": _rel(rr_path), "rerun_date": rr["run_date"],
        "rerun_hits_all_dates": rr["count_all_dates"], "rerun_hits_window": rr["count_window"],
        "rerun_hits_mesh_indexed_by_2026_04_04": len(indexed_by_april),
        "rerun_after_2026_04_04": sorted(rerun - set(indexed_by_april), key=int),
        "screened_records": len(pm_screened),
        "screened_run_m26": len(art["pubmed_m26"]), "screened_run_new": len(art["pubmed_new"]),
        "passed_run_m26": len(passed["pubmed_m26"]), "passed_run_new": len(passed["pubmed_new"]),
        "screened_in_rerun": len(pm_screened & rerun),
        "screened_not_in_rerun": [{"pmid": p, "title": es.get(p, {}).get("title", "")} for p in sorted(pm_screened - rerun, key=int)],
        "rerun_not_screened": [{"pmid": p, "title": es.get(p, {}).get("title", "")} for p in sorted(set(indexed_by_april) - pm_screened, key=int)],
    }
    emb = parse_embase()
    wos = parse_wos()
    out["embase"] = {"records": len(emb), "with_pmid": sum(bool(r["pmid"]) for r in emb),
                     "with_pmid_distinct": len({r["pmid"] for r in emb if r["pmid"]}),
                     "without_pmid": sum(not r["pmid"] for r in emb),
                     "without_pmid_distinct": len(group_records([r for r in emb if not r["pmid"]], lambda r: {norm_title(r["title"])})),
                     "screened_run_1": len(art["embase_1"]), "screened_run_2": len(art["embase_2"]),
                     "screened_distinct": len(art["embase_1"] | art["embase_2"]),
                     "passed_run_1": len(passed["embase_1"]), "passed_run_2": len(passed["embase_2"]),
                     "doi_rescreen_screened": len(art["embase_doi"]), "doi_rescreen_passed": len(passed["embase_doi"])}
    kw = re.compile(r"transcriptome|RNA.seq|microarray|gene expression|differentially expressed", re.I)
    out["wos"] = {"records": len(wos), "with_doi": sum(bool(r["doi"]) for r in wos),
                  "keyword_filter_pass_recomputed": sum(bool(kw.search(r["title"] + " " + r["abstract"])) for r in wos),
                  "screened": len(art["wos"]), "passed": len(passed["wos"])}
    out["bousman"] = {"screened": len(art["bousman"]), "passed": len(passed["bousman"])}
    # global deduplication of the three exports (not done in April; computed here)
    pm_recs = [{"pmid": p, "doi": es[p]["doi"], "title": es[p]["title"], "db": "pubmed"} for p in indexed_by_april]
    allrec = pm_recs + [{**r, "db": "embase"} for r in emb] + [{**r, "db": "wos"} for r in wos]
    by_id = group_records(allrec)
    by_id_title = group_records(allrec, lambda r: {"t:" + norm_title(r["title"])} if len(norm_title(r["title"])) > 30 else set())
    overlap = collections.Counter("+".join(sorted({allrec[i]["db"] for i in g})) for g in by_id_title)
    out["dedup"] = {"records_before": len(allrec), "pubmed_records_used": len(pm_recs),
                    "unique_by_pmid_or_doi": len(by_id), "unique_by_pmid_doi_or_title": len(by_id_title),
                    "overlap_pattern": dict(sorted(overlap.items()))}
    # biolit passes across the five 2026-04-04 runs: rows and distinct papers
    rows5 = [r for k in APRIL_RUNS for r in passed[k]]
    out["biolit_2026_04_04"] = {"screened_rows": sum(len(art[k]) for k in APRIL_RUNS), "pass_rows": len(rows5),
                                "pass_distinct_papers": len(group_records(rows5))}
    # the screening pool and Paul's decisions (search/screening/)
    pool = rjson(HERE / "screening" / "papers.json")
    dec = rjson(HERE / "screening" / "decisions.json")
    groups = group_records(pool)
    gdec = [{dec[pool[i]["id"]]["decision"] for i in g} for g in groups]
    inc_groups = [g for g, d in zip(groups, gdec) if d == {"include"}]
    exc_records = [r for r in pool if dec[r["id"]]["decision"] == "exclude"]
    src_inc = collections.Counter()
    db_of = {"rescreen": "pubmed", "wos": "wos", "embase": "embase", "embase-doi": "embase", "embase-xml-no-pmid": "embase"}

    def srcs(r):
        s = r.get("sources") or []
        return s if isinstance(s, list) else ast.literal_eval(s)
    for g in inc_groups:
        for s in sorted({db_of.get(s, s) for i in g for s in srcs(pool[i])}):
            src_inc[s] += 1
    out["pool"] = {"records": len(pool), "biolit_passed_true": sum(str(r.get("biolit_passed")) == "True" for r in pool),
                   "biolit_passed_unset": sum(str(r.get("biolit_passed")) != "True" for r in pool),
                   "papers": len(groups), "duplicate_records": len(pool) - len(groups),
                   "conflicting_papers": sum(len(d) > 1 for d in gdec),
                   "include_records": sum(v["decision"] == "include" for v in dec.values()),
                   "exclude_records": len(exc_records),
                   "include_papers": len(inc_groups), "exclude_papers": sum(d == {"exclude"} for d in gdec),
                   "exclude_records_with_note": sum(bool((dec[r["id"]].get("note") or "").strip()) for r in exc_records),
                   "include_papers_by_source_database": dict(sorted(src_inc.items()))}
    rev = rjson(REPO / "review-app" / "review_decisions.json")
    key2g = {k: gi for gi, g in enumerate(groups) for i in g for k in keys_of(pool[i])}
    sp_db = collections.Counter()
    sp_unmatched = []
    for k, v in rev.items():
        if v["decision"] != "include":
            continue
        gi = key2g.get(k.lower()) if k.lower() in key2g else key2g.get(k)
        if gi is None:
            sp_unmatched.append(k)
            continue
        for s_ in sorted({db_of.get(s_, s_) for i in groups[gi] for s_ in srcs(pool[i])}):
            sp_db[s_] += 1
    out["second_pass_include_by_database"] = dict(sorted(sp_db.items()))
    out["second_pass_include_not_in_pool"] = sp_unmatched
    out["second_pass"] = {"reviewed": len(rev), "include": sum(v["decision"] == "include" for v in rev.values()),
                          "exclude": sum(v["decision"] == "exclude" for v in rev.values()),
                          "model_rows_gpt4o_mini": len(rcsv(BIOLIT / "merikangas-eligibility" / "results.csv")),
                          "model_rows_sonnet_4_6": len(rcsv(BIOLIT / "merikangas-eligibility-sonnet" / "results.csv"))}
    # Merikangas's 26: which April search files hold them
    sets = rjson(CORPUS / "lineage_sets.json")
    m26 = sets["m26"]
    cache = rjson(CORPUS / "pubmed_cache.json")
    doi = {p: ((cache.get(p) or {}).get("doi") or "").lower() for p in m26}
    wos_d, emb_d = {r["doi"] for r in wos if r["doi"]}, {r["doi"] for r in emb if r["doi"]}
    emb_p = {r["pmid"] for r in emb if r["pmid"]}
    found = {p for p in m26 if p in pm_screened or (doi[p] and (doi[p] in wos_d or doi[p] in emb_d)) or p in emb_p}
    out["m26"] = {"n": len(m26), "in_pubmed_april_runs": sum(p in pm_screened for p in m26),
                  "in_pubmed_rerun": sum(p in rerun for p in m26),
                  "in_wos_export": sum(bool(doi[p]) and doi[p] in wos_d for p in m26),
                  "in_embase_export": sum((bool(doi[p]) and doi[p] in emb_d) or p in emb_p for p in m26),
                  "in_any_april_file": len(found), "not_in_april_files": sorted(set(m26) - found),
                  "wos_only": sorted(p for p in m26 if p not in pm_screened and p not in emb_p
                                     and not (doi[p] and doi[p] in emb_d) and doi[p] and doi[p] in wos_d)}
    out["_pool_rows"] = (pool, dec, groups)
    return out


# ---------------------------------------------------------------- stream B: September 2026 review batch (backward, broader, Gemma, cited)
def stream_september_batch():
    p = rjson(REPO / "review-app" / "backward-2026-09-22" / "papers.json")
    d = rjson(REPO / "review-app" / "backward-2026-09-22" / "review_decisions.json")
    by = collections.defaultdict(lambda: collections.Counter())
    for r in p:
        b = r.get("batch") or "backward reference lists"
        v = d[r["key"]]
        by[b][v["decision"]] += 1
        by[b][f"{v['decision']}_by_{decider(v.get('note', ''))}"] += 1
    out = {"batches": {k: dict(sorted(v.items())) for k, v in sorted(by.items())}, "papers": len(p),
           "include": sum(v["decision"] == "include" for v in d.values()),
           "exclude": sum(v["decision"] == "exclude" for v in d.values())}
    b1 = rjson(BIOLIT / "merikangas-backward" / "backward_candidates.json")
    b2 = rjson(BIOLIT / "merikangas-backward" / "backward_candidates_2.json")
    scr = set().union(*(artifact_ids(r) for r in BACKWARD_SCREEN))
    pas = [row for f in BACKWARD_PASS for row in rcsv(BIOLIT / f)]
    bscr = set().union(*(artifact_ids(r) for r in BROADER_RUNS))
    bpas = [row for r in BROADER_RUNS for f in sorted(use(BIOLIT / r).glob("results_bw2_chunk_*.csv")) for row in rcsv(f)]
    b1ids = {str(x[1]) for x in b1}
    out["backward"] = {"candidates_file": len(b1), "screened": len(scr), "screened_not_in_candidates": sorted(scr - b1ids),
                       "passed": len(pas)}
    out["broader"] = {"candidates_file": len(b2), "screened": len(bscr), "passed": len(bpas)}
    g = rjson(BIOLIT / "merikangas-eligibility-gemma" / "gemma_scz.json")
    chk = rjson(BIOLIT / "merikangas-eligibility-gemma" / "gemma_scz_check.json")
    out["gemma_check"] = {"datasets": len(g), "via_uri_and_query": sum("uri" in x["via"] for x in g.values()),
                          "via_query_only": sum(x["via"] == ["query"] for x in g.values()),
                          "status_annotated": dict(sorted(collections.Counter(r[0] for r in chk if r[1]).items())),
                          "status_query_only": dict(sorted(collections.Counter(r[0] for r in chk if not r[1]).items()))}
    elig = {}
    for name in ("merikangas-eligibility-backward", "merikangas-eligibility-backward2", "merikangas-eligibility-gemma",
                 "merikangas-eligibility-cited", "merikangas-eligibility-corpus"):
        elig[name] = sum(len(rcsv(f)) for f in sorted(use(BIOLIT / name).glob("run_*/*.csv")))
    out["eligibility_model_rows"] = elig
    out["_rows"] = (p, d)
    return out


# ---------------------------------------------------------------- other streams
def header_counts(path, prefix):
    for l in rtext(path).splitlines():
        if l.startswith(prefix):
            return l
    return ""


def parse_kv_counts(line):
    return {k.strip(): int(v) for k, v in re.findall(r"([A-Za-z][A-Za-z0-9 _-]*?) (\d+)(?:;|$)", line)}


def stream_forward():
    f = CORPUS / "forward_sweep" / "candidates.tsv"
    counts = parse_kv_counts(header_counts(f, "# counts:").removeprefix("# counts:"))
    seeds = int(re.search(r"# seeds \((\d+)\)", rtext(f)).group(1))
    s2 = rjson(CORPUS / "forward_sweep" / "cache" / "screen2_summary.json")
    af = rjson(CORPUS / "added_forward.json")
    return {"seeds": seeds, "counts": counts, "sonnet": s2["counts"], "sonnet_fulltext_used": s2["fulltext_used"],
            "added": len([k for k in af if not k.startswith("_")]) + len(af.get("_removed", {})),
            "removed_at_gene_list_step": len(af.get("_removed", {}))}


def stream_single_gene():
    s = rjson(CORPUS / "single_gene_sweep" / "summary.json")
    gw = rtsv(CORPUS / "single_gene_sweep" / "genome_wide.tsv")
    ag = rjson(CORPUS / "added_genomewide_2026-09-24.json")
    ap = rjson(CORPUS / "added_proteomics_2026-09-24.json")
    rej = ag.get("_rejected", {})
    return {"keyword_queries": sum(not k.startswith(("au:", "src:")) for k in s["query_counts"]),
            "author_queries": sum(k.startswith("au:") for k in s["query_counts"]),
            "hits_by_query": s["query_counts"], "count_only": s["count_only"], "union": s["union"],
            "already_in_corpus": s["in_corpus"], "already_decided": s["already_decided"], "pool": s["pool"],
            "screened": s["screened"], "haiku": s["haiku_decisions"], "rescue_sent": s["rescue_sent"],
            "rescue_answers": s["rescue_answers"], "final": s["final_decisions"],
            "genome_wide_rows": len(gw), "genome_wide_classes": dict(sorted(collections.Counter(r["class"] for r in gw).items())),
            "batch10_transcriptome_include": len([k for k in ag if not k.startswith("_")]),
            "batch10_transcriptome_reject": sum(not v.get("duplicate_of") for v in rej.values()),
            "batch10_transcriptome_duplicate": sum(bool(v.get("duplicate_of")) for v in rej.values()),
            "batch10_pending": len(ag.get("_pending", {})),
            "batch10_proteomics_include": len([k for k in ap if not k.startswith("_")]),
            "batch10_proteomics_reject": len(ap.get("_rejected", {})),
            "queries_by_date": dict(sorted(collections.Counter(v for k, v in s["query_added"].items() if not k.startswith("src:")).items())),
            "other_sources": s.get("other_sources", {}),
            "papers_json_includes": sum(v for k, v in s["papers_json_by_batch"].items() if not k.endswith("unsure")),
            "includes_by_batch_date": s.get("includes_by_batch_date", {})}


def stream_single_cell(idx):
    sc = rjson(CORPUS / "single_cell.json")
    fate = collections.Counter()
    for grp in ("studies", "studies_by_hand", "preprints_by_hand"):
        for s in sc[grp]:
            key = s.get("pmid") or s.get("doi") or ""
            n = find_node(idx, key)
            merged = any(m.get("pmid") == key for x in idx.values() for m in x.get("merged_sources", []))
            fate[("preprint merged into its published paper" if not n and merged else "not in lineage" if not n
                  else "outside source node" if n.get("external") else n["source"])] += 1
    return {"query": sc["query"], "model": sc["model"], "candidates": sc["n_candidates"], "model_include": len(sc["studies"]),
            "by_hand": len(sc["studies_by_hand"]), "preprints_by_hand": len(sc["preprints_by_hand"]),
            "fate_in_lineage": dict(sorted(fate.items()))}


def stream_hand():
    h = rjson(CORPUS / "added_by_hand.json")
    rows = {k: v for k, v in h.items() if not k.startswith("_")}
    kind = collections.Counter("source audit 2026-09-25" if "source audit" in v.get("added_by", "")
                               else "Pritzker PubMed sweep 2026-09-26" if "Pritzker" in v.get("added_by", "")
                               else "other hand additions" for v in rows.values())
    sa = rjson(CORPUS / "source_audit_2026-09-25.json")
    return {"n": len(rows), "by_route": dict(sorted(kind.items())), "source_audit": sa["counts"],
            "source_audit_unresolved": sa["unresolved"]}


def md_section(text, start, end):
    i = text.index(start)
    j = text.index(end, i + len(start)) if end in text[i + len(start):] else len(text)
    return text[i:j]


def stream_surveys(idx):
    geo = rjson(REPO / "side-studies" / "gemma_candidates" / "geo_sz_sweep.json")
    wide = rjson(REPO / "side-studies" / "gemma_candidates" / "gemma_wide_sweep.json")
    status = rjson(REPO / "side-studies" / "gemma_candidates" / "gemma_status.json")
    harvest = rjson(REPO / "side-studies" / "gemma_candidates" / "harvest.json")
    gas = rjson(REPO / "side-studies" / "gemma_candidates" / "geo_author_search.json")
    upd = CORPUS / "forward_sweep" / "update_2026-09-27" / "candidates.tsv"
    head = rtext(upd).splitlines()[:5]
    hl = " ".join(head)
    urows = rtsv(upd)
    leads = rtext(REPO / "side-studies" / "gemma_candidates" / "new_leads_2026-09-27.md")
    link = r"\[(\d{7,9})\]\(https://pubmed"
    new_data = sorted(set(re.findall(link, md_section(leads, "## 2. Literature: new data", "## 3."))), key=int)
    reuse_sec = md_section(leads, "## 3. Literature: reuses", "## 4.")
    reuse = sorted(set(re.findall(link, reuse_sec)), key=int)
    reading = [m for m in re.findall(r"\|\s*\[(\d{7,9})\][^\n]*\n", reuse_sec) if re.search(rf"\[{m}\][^\n]*reading list", reuse_sec)]
    in_lineage = lambda ids: sum(bool(find_node(idx, p)) and not find_node(idx, p).get("external") for p in ids)
    doi_of = {r["pmid"]: (r.get("doi") or "").lower() for r in urows}
    cls = collections.Counter(r["class"] for r in urows if r["class"])
    sreuse = collections.Counter(r["sonnet_reuse"] for r in urows if r["sonnet_reuse"])
    return {
        "geo_gemma_2026_09_25": {"geo_term": geo["geo_term"], "geo_series": len(geo["geo"]), "gemma_query_datasets": len(geo["gemma"]),
                                 "harvest_targets": harvest["n_targets"], "geo_author_searches": len(gas),
                                 "gemma_status_datasets": len(status)},
        "gemma_wide_2026_09_27": {"inspected": wide["n_inspected"], "with_sz_samples": len(wide["datasets"]),
                                  "not_in_2026_09_25_survey": len(set(wide["datasets"]) - set(status))},
        "forward_update_2026_09_27": {
            "seeds": int(re.search(r"seeds (\d+)", hl).group(1)),
            "seeds_new": int(re.search(r"new since 2026-09-23: (\d+)", hl).group(1)),
            "delta_pool": int(re.search(r"delta pool (\d+)", hl).group(1)),
            "delta_prefilter_pass": int(re.search(r"delta prefilter pass (\d+)", hl).group(1)),
            "old_exclude_pool": int(re.search(r"excluded_2026-09-23 pool (\d+)", hl).group(1)),
            "class_screened_rows": len(urows), "class": dict(sorted(cls.items())),
            "sonnet_reuse": dict(sorted(sreuse.items())),
            "kept_new_data_by_hand": len(new_data), "reuse_listed": len(reuse), "reuse_from_reading_lists": len(reading),
            "kept_new_data_in_lineage": in_lineage(new_data), "reuse_in_lineage": in_lineage(reuse),
            "_lead_keys": {p: {p, doi_of.get(p, "")} - {""} for p in new_data + reuse}},
    }


# ---------------------------------------------------------------- final corpus and downstream
def final_corpus(nodes, sep):
    st = corpus_studies(nodes)
    p, d = sep["_rows"]
    batch = {r["key"]: (r.get("batch") or "backward reference lists") for r in p}
    hand = rjson(CORPUS / "added_by_hand.json")
    sub = collections.Counter()
    for n in st:
        s = n["source"]
        if s == "Added Sep 2026 (reference lists)":
            k = n["id"]
            if k in d and d[k]["decision"] == "include":
                sub[f"{s} | {batch[k]}: include"] += 1
            elif k in d:
                sub[f"{s} | {batch[k]}: excluded, re-admitted as a one-gene study (group 2)"] += 1
            else:
                sub[f"{s} | not in the review batch (group {n.get('category')})"] += 1
        elif s == "Added by hand":
            ab = (hand.get(n["id"]) or {}).get("added_by", "")
            sub[f"{s} | " + ("source audit 2026-09-25" if "source audit" in ab else "Pritzker PubMed sweep 2026-09-26"
                             if "Pritzker" in ab else "other hand additions")] += 1
        else:
            sub[f"{s} | all"] += 1
    cats = collections.Counter("sc" if n.get("source") == "Single-cell / single-nucleus" and not n.get("category") else str(n.get("category")) for n in st)
    summ = rtsv(GL / "summary.tsv")
    ex = [r for r in summ if r["status"].split(":")[0].strip() == "extracted"]
    ex_studies = sorted({r["list"].split("__")[0] for r in ex})
    idx = node_index(st)
    outside = [k for k in ex_studies if not find_node(idx, k)]
    rec = rtsv(GL / "recurrence.tsv")
    head = header_counts(GL / "recurrence.tsv", "# 2")
    m = re.search(r"# (\d+) in-corpus studies with a primary list; left out: (\d+) meta-lane, (\d+) multi-series reanalyses.*?(\d+) other reanalyses", head)
    pb = next(r for r in rec if r["scenario"] == "published, brain" and r["tier"] == "study")
    pc = next(r for r in rec if r["scenario"] == "published, brain" and r["tier"] == "cohort")
    nondel = re.search(r"(\d+) studies had a non-DE primary list dropped", rtext(GL / "recurrence.tsv"))
    return {
        "studies": len(st),
        "by_source": dict(sorted(collections.Counter(n["source"] for n in st).items())),
        "by_source_detail": dict(sorted(sub.items())),
        "by_group": dict(sorted(cats.items())),
        "outside_source_nodes": sum(bool(n.get("external")) for n in nodes),
        "gene_lists_extracted": len(ex), "studies_with_extracted_list": len(ex_studies),
        "extracted_studies_outside_corpus": outside,
        "list_provenance": dict(sorted(collections.Counter(r["provenance"] for r in ex).items())),
        "summary_rows": len(summ), "summary_studies": len({r["list"].split("__")[0] for r in summ}),
        "recurrence_primary_in_corpus": int(m.group(1)), "recurrence_left_out_meta_lane": int(m.group(2)),
        "recurrence_left_out_multi_series": int(m.group(3)), "recurrence_left_out_other_reanalyses": int(m.group(4)),
        "published_brain_candidates": int(pb["n_lists_old"]), "published_brain_non_de_dropped": int(pb["n_lists_dropped"]),
        "non_de_primary_dropped_all_scenarios": int(nondel.group(1)) if nondel else None,
        "published_brain_studies": int(pb["n_lists"]), "published_brain_in_clusters": int(pc["n_lists"]),
        "cohort_clusters": int(pc["n_clusters"]),
    }


# ---------------------------------------------------------------- item 16b
def list_16b(april, sep, nodes):
    rows = []
    pool, dec, groups = april["_pool_rows"]
    for g in groups:
        recs = [pool[i] for i in g]
        if {dec[r["id"]]["decision"] for r in recs} != {"exclude"}:
            continue
        note = "; ".join(sorted({(dec[r["id"]].get("note") or "").strip() for r in recs} - {""}))
        if not note:
            continue
        r = recs[0]
        rows.append({"id": r.get("pmid") or r.get("doi"), "label": f"{(r.get('authors') or '').split(',')[0]} {r.get('year', '')}".strip(),
                     "title": r.get("title", ""), "stream": "April 2026 databases", "stage": "title/abstract + biolit fields (screening UI, 2026-04-16)",
                     "reason": note, "decided_by": "Paul", "file_of_record": "search/screening/decisions.json"})
    rev = rjson(REPO / "review-app" / "review_decisions.json")
    rpapers = {x.get("key") or x.get("id") or x.get("pmid"): x for x in rjson(REPO / "review-app" / "papers.json")}
    for k, v in rev.items():
        if v["decision"] == "exclude":
            x = rpapers.get(k, {})
            rows.append({"id": k, "label": f"{(x.get('authors') or '').split(',')[0]} {x.get('year', '')}".strip(), "title": x.get("title", ""),
                         "stream": "April 2026 databases", "stage": f"second pass, full text (review app, {v.get('updated', '')[:10]})",
                         "reason": v.get("note") or "", "decided_by": "not recorded in the file",
                         "file_of_record": "review-app/review_decisions.json"})
    p, d = sep["_rows"]
    st_ids = {n["id"] for n in corpus_studies(nodes)}
    for r in p:
        v = d[r["key"]]
        if v["decision"] != "exclude":
            continue
        readm = r["key"] in st_ids
        rows.append({"id": r["key"], "label": f"{(r.get('authors') or '').split(',')[0]} {r.get('year', '')}".strip(), "title": r.get("title", ""),
                     "stream": "Sep 2026 " + (r.get("batch") or "backward reference lists"),
                     "stage": "abstract + model eligibility fields (review app, 2026-09-22)" + ("; later re-admitted as a one-gene study (group 2)" if readm else ""),
                     "reason": v.get("note", ""), "decided_by": decider(v.get("note", "")),
                     "file_of_record": "review-app/backward-2026-09-22/review_decisions.json"})
    for fname, stream in (("added_genomewide_2026-09-24.json", "Sep 2026 single-gene sweep, genome-wide class (batch 10)"),
                          ("added_proteomics_2026-09-24.json", "Sep 2026 single-gene sweep, proteomics (batch 10)")):
        a = rjson(CORPUS / fname)
        for k, v in a.get("_rejected", {}).items():
            rows.append({"id": k, "label": v.get("label", ""), "title": "", "stream": stream, "stage": "full text read by a worker (batch 10 verification)",
                         "reason": v.get("reason", "") + (f" (duplicate of {v['duplicate_of']})" if v.get("duplicate_of") else ""),
                         "decided_by": "Claude (worker verification)", "file_of_record": f"corpus/{fname} (_rejected)"})
        for k, v in a.get("_pending", {}).items():
            rows.append({"id": k, "label": v.get("label", ""), "title": "", "stream": stream, "stage": "pending Paul's decision",
                         "reason": v.get("worker_reason", "") + " | " + v.get("why", ""), "decided_by": "pending (Paul)",
                         "file_of_record": f"corpus/{fname} (_pending)"})
    af = rjson(CORPUS / "added_forward.json")
    for k, v in af.get("_removed", {}).items():
        rows.append({"id": k, "label": "", "title": "", "stream": "Sep 2026 forward citations", "stage": "gene-list step (after inclusion)",
                     "reason": v.get("removed", ""), "decided_by": "Paul" if "Paul:" in v.get("removed", "") else "Claude (gene-list step)",
                     "file_of_record": "corpus/added_forward.json (_removed)"})
    for r in rows:
        r["reason_category"] = reason_category(r["reason"])
    rows.sort(key=lambda r: (r["stream"], r["reason_category"], r["id"]))
    return rows


# ---------------------------------------------------------------- the search log
LOG_COLS = ["search_id", "date", "stream", "source", "exact_query", "filters", "hits", "after_dedup", "screened",
            "included", "file_of_record", "screening_method", "decided_by", "count_notes"]


def search_log(C):
    A, S, F, G, SC, H, SV = C["april"], C["september_batch"], C["forward"], C["single_gene"], C["single_cell"], C["hand"], C["surveys"]
    pm, emb, wos, pool = A["pubmed"], A["embase"], A["wos"], A["pool"]
    inc_db = A["second_pass_include_by_database"]
    b = S["batches"]
    fu = SV["forward_update_2026_09_27"]
    ovid = ("schizophrenia/ and case control study/ and human/ and article/ and (transcriptome/ or microarray analysis/ or "
            "RNA sequencing.mp. or gene expression/) [mp=title, abstract, heading word, drug trade name, original title, device "
            "manufacturer, drug manufacturer, device trade name, keyword heading word, floating subheading word, candidate term word]")
    wosq = "TS=(schizophr*) AND (TS=(Transcriptome) OR TS=(Microarray analysis) OR TS=(Transcription) OR TS=(RNA sequence)) AND TS=(case control)"
    sg_queries = " || ".join(l.split(": ", 1)[1] for l in rtext(CORPUS / "single_gene_sweep" / "candidates.tsv").splitlines()
                             if l.startswith("# query "))
    rows = [
        dict(search_id="S01", date="2026-04-04", stream="A. April 2026 database replication", source="PubMed (NCBI E-utilities)",
             exact_query=pm["query"], filters="none (Merikangas window 2000/01/01-2020/06/30 counted separately)",
             hits="NR", after_dedup=A["dedup"]["unique_by_pmid_doi_or_title"], screened=pm["screened_records"],
             included=inc_db.get("pubmed", "NR"), file_of_record=f"vault biolit-results/{APRIL_RUNS['pubmed_m26']} and /{APRIL_RUNS['pubmed_new']} (artifacts, results.csv); {pm['rerun_file']}",
             screening_method="biolit criterion screen (config merikangas-review.json; model not recorded), then Paul's manual title/abstract screen (screening UI, 2026-04-16)",
             decided_by="Paul (inclusion); model (pre-screen)",
             count_notes=(f"hits: the April hit list was not saved; the search note states 87 (73 in window). The 2026-09-27 rerun of the same query "
                          f"returns {pm['rerun_hits_all_dates']} ({pm['rerun_hits_window']} in window), {pm['rerun_hits_mesh_indexed_by_2026_04_04']} of them "
                          f"MeSH-indexed by 2026-04-04. after_dedup is the three-database total (S01-S03). included (S01-S03) = second-pass include papers whose pool "
                          f"record lists that database (a paper can list several).")),
        dict(search_id="S02", date="2026-04-04", stream="A. April 2026 database replication", source="Embase (Ovid)",
             exact_query=ovid, filters="none (no date limit)", hits=emb["records"], after_dedup=A["dedup"]["unique_by_pmid_doi_or_title"],
             screened=emb["screened_distinct"] + emb["doi_rescreen_screened"], included=inc_db.get("embase", "NR"),
             file_of_record="search/materials/merikangas-embase-citation.part1.xml, part2.xml; vault biolit-results/merikangas-embase/, merikangas-embase-doi/",
             screening_method="biolit criterion screen (two runs over the same PMIDs; DOI rescreen 2026-04-16 with GPT-4o-mini), then Paul's manual screen",
             decided_by="Paul (inclusion); model (pre-screen)",
             count_notes=f"{emb['with_pmid']} records with a PMID, {emb['without_pmid']} without; screened = distinct PMIDs screened plus DOI rescreen."),
        dict(search_id="S03", date="2026-04-04", stream="A. April 2026 database replication", source="Web of Science Core Collection",
             exact_query=wosq, filters="publication years 2000 to April 2026 (Paul's query note: 'With date limit to 2000-4/2026')",
             hits=wos["records"], after_dedup=A["dedup"]["unique_by_pmid_doi_or_title"], screened=wos["screened"], included=inc_db.get("wos", "NR"),
             file_of_record="search/materials/merikangas-webofscience.bib; vault biolit-results/merikangas-wos/",
             screening_method="keyword pre-filter on title/abstract, then biolit criterion screen, then Paul's manual screen",
             decided_by="Paul (inclusion); rule + model (pre-screen)",
             count_notes=(f"screened = artifacts of the biolit run (the search note states 168). The keyword pre-filter as written in the search note, "
                          f"recomputed on the export, passes {wos['keyword_filter_pass_recomputed']}; the April pre-filter also dropped records already in Zotero, "
                          f"which cannot be recomputed.")),
        dict(search_id="S04", date="2026-04-16", stream="A. April 2026 database replication", source="PubMed record (retrospective screen)",
             exact_query="PMID 19582768 (Bousman 2010)", filters="n/a", hits=1, after_dedup=1, screened=A["bousman"]["screened"],
             included=0, file_of_record="vault biolit-results/merikangas-bousman/", screening_method="biolit criterion screen",
             decided_by="model", count_notes="In PubMed results but not fed to the April biolit run; one of Merikangas's 26, so not a new include."),
        dict(search_id="S05", date="2026-09-22", stream="A. April 2026 database replication", source="PubMed (rerun of S01)",
             exact_query=pm["query"], filters="none", hits="NR", after_dedup="NR", screened="NR", included="NR",
             file_of_record="tracker, Step 1 item 6 (no file of the hit list)", screening_method="coverage check only", decided_by="n/a",
             count_notes="The tracker records that this rerun returned 18 of the 27 PMIDs then tabled as Merikangas's; no hit list or count was saved."),
        dict(search_id="S06", date=pm["rerun_date"], stream="A. April 2026 database replication", source="PubMed (rerun of S01, this record)",
             exact_query=pm["query"], filters="none; window count with datetype=pdat 2000/01/01-2020/06/30",
             hits=pm["rerun_hits_all_dates"], after_dedup=pm["rerun_hits_all_dates"], screened=0, included=0, file_of_record=pm["rerun_file"],
             screening_method="none (verification of S01)", decided_by="n/a",
             count_notes=f"{pm['rerun_hits_window']} in the Merikangas window; {A['m26']['in_pubmed_rerun']} of Merikangas's {A['m26']['n']} are in the hit list."),
        dict(search_id="S07", date="2026-09-22", stream="B. Reference lists (backward citation)", source="Europe PMC reference lists of included and parent studies",
             exact_query="reference lists of the included/parent studies; filter: expression-profiling MeSH or 'gene expression profil' in the title (filter code not saved; reconstructed in search/backward_reconstruction)",
             filters="schizophrenia expression papers not already retrieved or screened", hits=S["backward"]["candidates_file"],
             after_dedup=S["backward"]["candidates_file"], screened=S["backward"]["screened"], included=b["backward reference lists"].get("include", 0),
             file_of_record="vault biolit-results/merikangas-backward/backward_candidates.json, run_20260922_165405, run_20260922_184508, run_20260922_190415; review-app/backward-2026-09-22/",
             screening_method="biolit criterion screen; Sonnet 5 eligibility fields (merikangas-eligibility-lineage.json); decisions in the review app",
             decided_by=f"Paul {b['backward reference lists'].get('include_by_Paul', 0) + b['backward reference lists'].get('exclude_by_Paul', 0)}, Claude or unrecorded for the rest",
             count_notes=(f"The tracker states the reference lists of 50 studies were used; the number of studies and of references read is not in a file. "
                          f"screened includes {len(S['backward']['screened_not_in_candidates'])} papers added by hand (Mirnics 2000, Vawter 2002, Tkachev 2003); "
                          f"{S['backward']['passed']} passed biolit." + _reconstruction("S07"))),
        dict(search_id="S08", date="2026-09-22", stream="B. Reference lists (backward citation)", source="Europe PMC reference lists (broader filter)",
             exact_query="same reference lists; broader filter after the first filter dropped Fromer 2016 (filter code not saved)",
             filters="as S07", hits=S["broader"]["candidates_file"], after_dedup=S["broader"]["candidates_file"], screened=S["broader"]["screened"],
             included=b["broader reference filter"].get("include", 0),
             file_of_record="vault biolit-results/merikangas-backward/backward_candidates_2.json, run_20260922_20384[2-5]; review-app/backward-2026-09-22/",
             screening_method="biolit criterion screen; Sonnet 5 eligibility fields; decisions in the review app",
             decided_by="Claude (broader-filter pass)", count_notes=f"{S['broader']['passed']} passed biolit." + _reconstruction("S08")),
        dict(search_id="S09", date="2026-09-22", stream="B. Gemma completeness check", source="Gemma 2 REST (gemma2.msl.ubc.ca)",
             exact_query="datasets matching the text query 'schizophrenia' and datasets annotated MONDO_0005090 (exact REST calls not saved)",
             filters="human datasets", hits=S["gemma_check"]["datasets"], after_dedup=S["gemma_check"]["datasets"],
             screened=b["Gemma schizophrenia datasets"]["include"] + b["Gemma schizophrenia datasets"]["exclude"],
             included=b["Gemma schizophrenia datasets"]["include"],
             file_of_record="vault biolit-results/merikangas-eligibility-gemma/gemma_scz.json, gemma_scz_check.json; review-app/backward-2026-09-22/",
             screening_method="papers behind datasets never seen by any search, Sonnet 5 eligibility fields, then review app",
             decided_by="Claude (Gemma check)", count_notes=f"hits = datasets; screened and included = papers. {S['gemma_check']['via_uri_and_query']} datasets carry the MONDO annotation."),
        dict(search_id="S10", date="2026-09-22", stream="B. Cited sample sources", source="methods and reference lists of corpus studies (model-read, corpus/sources_llm.py)",
             exact_query="prompt in corpus/sources_llm.py (PROMPT; criteria.md)", filters="cited papers outside the corpus that are expression studies",
             hits="NR", after_dedup="NR", screened=b["cited as a data source"]["include"] + b["cited as a data source"]["exclude"],
             included=b["cited as a data source"]["include"],
             file_of_record="corpus/sources_llm.json; vault biolit-results/merikangas-eligibility-cited/; review-app/backward-2026-09-22/",
             screening_method="model reads methods and reference list; Sonnet 5 eligibility fields; review app", decided_by="Claude (cited-source pass)",
             count_notes="The tracker states 121 cited papers outside the corpus, of which 56 were expression studies; the 121 is not recomputed here."),
        dict(search_id="S11", date="2026-09-23", stream="C. Single-cell / single-nucleus", source="PubMed (NCBI E-utilities) plus DOIs cited in a vault note",
             exact_query=SC["query"], filters="NOT review[pt]", hits=SC["candidates"], after_dedup=SC["candidates"], screened=SC["candidates"],
             included=SC["fate_in_lineage"].get("Single-cell / single-nucleus", 0), file_of_record="corpus/single_cell.json",
             screening_method=f"one {SC['model']} call per abstract; plus {SC['by_hand']} study and {SC['preprints_by_hand']} preprints added by hand",
             decided_by="model; hand additions by Claude/Paul (not recorded which)",
             count_notes=f"hits = query plus vault-note DOIs (the query's own count is not saved separately). Fate of the {SC['model_include'] + SC['by_hand'] + SC['preprints_by_hand']} identified: {SC['fate_in_lineage']}."),
        dict(search_id="S12", date="2026-09-23", stream="D. Forward citation", source="Europe PMC citations (`MED/{pmid}/citations`) of corpus seed studies",
             exact_query=f"citing papers of {F['seeds']} seeds (seed list in the '# seeds' line of corpus/forward_sweep/candidates.tsv); prefilter regex in corpus/forward_sweep/README.md",
             filters="MEDLINE id, published 2002 or later, not already in lineage.json", hits=F["counts"].get("citation links"),
             after_dedup=F["counts"].get("candidate pool"), screened=F["counts"].get("prefilter pass"), included=C["final"]["by_source"].get("Added Sep 2026 (forward citations)", 0),
             file_of_record="corpus/forward_sweep/candidates.tsv, screen2_results.csv, cache/; corpus/added_forward.json",
             screening_method="regex prefilter; claude-haiku-4-5 screen; claude-sonnet-5 second pass (stricter criterion)",
             decided_by="model only (Paul, 2026-09-24: add all Sonnet includes; remove any found irrelevant at the gene-list step)",
             count_notes=f"{F['sonnet'].get('include')} Sonnet includes added, {F['removed_at_gene_list_step']} removed at the gene-list step."),
        dict(search_id="S13", date="2026-09-24", stream="E. Single-gene sweep (genome-wide and proteomic classes)", source="PubMed (NCBI E-utilities)",
             exact_query=sg_queries, filters="count-only query not screened (serum/plasma protein)", hits=sum(G["hits_by_query"].values()),
             after_dedup=G["union"], screened=G["screened"],
             included=C["final"]["by_source"].get("Added Sep 2026 (genome-wide sweep)", 0) + C["final"]["by_source"].get("Added Sep 2026 (proteomics)", 0),
             file_of_record="corpus/single_gene_sweep/candidates.tsv, summary.json, genome_wide.tsv, batch10_verification/; corpus/added_genomewide_2026-09-24.json, added_proteomics_2026-09-24.json",
             screening_method="claude-haiku-4-5 screen, claude-sonnet-5 rescue pass and genome-wide class pass, then full-text verification by worker agents",
             decided_by="Claude workers (full text); Paul for Chen 2024 (pending) and the proteomics rule",
             count_notes=(f"hits = sum over {G['keyword_queries']} keyword and {G['author_queries']} author queries (added by date: "
                          + ", ".join(f"{k} {v}" for k, v in G["queries_by_date"].items()) + "), plus PMIDs found by citation (other sources, 2026-10-01: "
                          + ", ".join(f"{k[4:]} {v}" for k, v in G["other_sources"].items()) + "; corpus/single_gene_sweep/recall_2026-10-01/). "
                          f"The {G['papers_json_includes']} single-gene includes ({G['final'].get('include')} screened, the rest corpus group 2 and review-app single-gene decisions; "
                          f"by batch date " + ", ".join(f"{k} {v}" for k, v in G["includes_by_batch_date"].items()) + ") are H4 units, not corpus studies.")),
        dict(search_id="S14", date="2026-09-25", stream="F. Cited-source audit", source="sources_class.json items, outside source nodes, pooled datasets (corpus/source_audit_2026-09-25.json)",
             exact_query="audit by hand (batch 11 fork); PubMed searches to resolve pooled 'paper' datasets (queries not saved)", filters="n/a",
             hits=H["source_audit"]["checked"], after_dedup=H["source_audit"]["checked"], screened=H["source_audit"]["checked"],
             included=H["source_audit"]["by_decision"].get("added", 0), file_of_record="corpus/source_audit_2026-09-25.json; corpus/added_by_hand.json",
             screening_method="by hand (Claude), classes i/ii/iii", decided_by="Claude", count_notes=f"decisions: {H['source_audit']['by_decision']}"),
        dict(search_id="S15", date="2026-09-25", stream="G. GEO and Gemma data survey", source="GEO DataSets (gds, E-utilities) and Gemma 2 text search",
             exact_query=SV["geo_gemma_2026_09_25"]["geo_term"] + " || Gemma /datasets query: schizophrenia, schizophrenic, schizoaffective, psychosis",
             filters="human series", hits=SV["geo_gemma_2026_09_25"]["geo_series"] + SV["geo_gemma_2026_09_25"]["gemma_query_datasets"],
             after_dedup="NR", screened="NR", included=0, file_of_record="side-studies/gemma_candidates/geo_sz_sweep.json, harvest.json, geo_author_search.json, gemma_candidates.tsv",
             screening_method="data-location survey for corpus studies (not a study search)", decided_by="Claude",
             count_notes=f"{SV['geo_gemma_2026_09_25']['geo_series']} GEO series and {SV['geo_gemma_2026_09_25']['gemma_query_datasets']} Gemma datasets; no study entered the corpus from it."),
        dict(search_id="S16", date="2026-09-26", stream="H. Pritzker collection", source="PubMed (E-utilities) and GEO DataSets",
             exact_query="Pritzker AND schizophrenia AND (expression OR transcriptome OR microarray OR RNA-seq); consortium-author queries (exact strings not saved)",
             filters="n/a", hits="NR", after_dedup="NR", screened="NR", included=H["by_route"].get("Pritzker PubMed sweep 2026-09-26", 0),
             file_of_record="corpus/pritzker_datasets_2026-09-26.md; corpus/added_by_hand.json", screening_method="by hand", decided_by="Paul and Claude",
             count_notes="corpus/pritzker_datasets_2026-09-26.md states 82 records for the author queries and 12 GEO series; not in a machine-readable file."),
        dict(search_id="S17", date="2026-09-27", stream="G. GEO and Gemma data survey", source="Gemma 2 REST, authenticated (gemma2.msl.ubc.ca)",
             exact_query="filters allCharacteristics.valueUri = MONDO_0005090; allCharacteristics.value like schizophreni / schizoaffective / psychos / psychotic; queries schizophrenia, schizophrenic, schizoaffective, psychosis; private non-GEO human imports",
             filters="human; kept when a sample characteristic matches schizo|psychos|psychotic", hits=SV["gemma_wide_2026_09_27"]["inspected"],
             after_dedup=SV["gemma_wide_2026_09_27"]["inspected"], screened=SV["gemma_wide_2026_09_27"]["inspected"], included=0,
             file_of_record="side-studies/gemma_candidates/gemma_wide_sweep.json, new_leads_2026-09-27.md", screening_method="rule (sample characteristics), then by hand",
             decided_by="Claude; admission to the corpus awaits Paul",
             count_notes=f"{SV['gemma_wide_2026_09_27']['with_sz_samples']} datasets with SZ or psychosis samples, {SV['gemma_wide_2026_09_27']['not_in_2026_09_25_survey']} not in the 2026-09-25 survey."),
        dict(search_id="S18", date="2026-09-27", stream="D. Forward citation (update)", source="Europe PMC citations of the current lineage seeds",
             exact_query=f"citing papers of {fu['seeds']} seeds ({fu['seeds_new']} new since 2026-09-23); criteria of sweep.py, screen2.py and screen_genome_wide.py",
             filters="as S12", hits=fu["delta_pool"] + fu["old_exclude_pool"], after_dedup=fu["class_screened_rows"], screened=fu["class_screened_rows"], included=0,
             file_of_record="corpus/forward_sweep/update_2026-09-27/candidates.tsv, run.log; side-studies/gemma_candidates/new_leads_2026-09-27.md",
             screening_method="claude-sonnet-5 class screen and reuse screen (abstract only); new-data class triaged by hand from titles and model fields",
             decided_by="Claude (triage); admission to the corpus awaits Paul",
             count_notes=(f"hits = {fu['delta_pool']} new citing papers + {fu['old_exclude_pool']} 2026-09-23 Haiku excludes re-classed; "
                          f"{fu['kept_new_data_by_hand']} new-data papers kept and {fu['reuse_listed']} reuse papers listed ({fu['reuse_from_reading_lists']} of them from the vault reading lists); none added to the lineage.")),
    ]
    for r in rows:
        for k in LOG_COLS:
            r.setdefault(k, "")
    return rows


# ---------------------------------------------------------------- criteria.md (verbatim)
CRITERIA_SOURCES = [
    ("corpus/forward_sweep/sweep.py", ("SYSTEM", "CRITERION")),
    ("corpus/forward_sweep/screen2.py", ("SYSTEM", "CRITERION")),
    ("corpus/single_gene_sweep/sweep_single_gene.py", ("SYSTEM", "CRITERION", "RESCUE_Q")),
    ("corpus/single_gene_sweep/screen_genome_wide.py", ("CRITERION",)),
    ("corpus/single_cell.py", ("QUERY",)),
    ("corpus/sources_llm.py", ("PROMPT",)),
]


def criteria_md():
    L = ["# Screening and eligibility criteria, verbatim", "",
         "Generated by `search/build_prisma.py` from the scripts and biolit configs named below; do not edit by hand. "
         "Each block is the string constant or config field as it stands in the current file (sha256 given); the scripts combine "
         "these with an output specification (the OUT_SPEC / OUT constants) before sending. Whether a file changed after its run "
         "is not recorded here; the scripts' own caches and run logs hold the run-time state.", ""]
    for rel, names in CRITERIA_SOURCES:
        p = REPO / rel
        tree = ast.parse(rtext(p))
        sha = INPUTS[_rel(p)][:12]
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in node.targets):
                name = next(t.id for t in node.targets if isinstance(t, ast.Name))
                try:
                    val = ast.literal_eval(node.value)
                except ValueError:
                    continue
                L += [f"## `{rel}` {name} (file sha256 {sha})", "", "```text", val.rstrip(), "```", ""]
    for cfg in ("merikangas-review.json", "merikangas-eligibility.json", "merikangas-eligibility-sonnet.json", "merikangas-eligibility-lineage.json"):
        p = CONFIGS / cfg
        d = rjson(p)
        L += [f"## biolit config `{cfg}` (vault biolit-configs; sha256 {INPUTS[_rel(p)][:12]})", "",
              f"Model: {d.get('provider', 'not set')} / {d.get('model', 'not set in the config')}.", "", "```json",
              json.dumps({k: d[k] for k in ("criterion", "fields") if k in d}, indent=1, ensure_ascii=False), "```", ""]
    return "\n".join(L)


# ---------------------------------------------------------------- PRISMA.md generated blocks
def fmt(d):
    return ", ".join(f"{k} {v}" for k, v in d.items())


def mermaid(lines):
    return "```mermaid\nflowchart TD\n" + "\n".join("  " + l for l in lines) + "\n```"


def box(i, text):
    return f'{i}["{text}"]'


def blocks(C, log, r16):
    A, S, F, G, SC, H, SV, FC = (C[k] for k in ("april", "september_batch", "forward", "single_gene", "single_cell", "hand", "surveys", "final"))
    pm, emb, wos, pool, sp, dd = A["pubmed"], A["embase"], A["wos"], A["pool"], A["second_pass"], A["dedup"]
    out = {}
    t = ["| id | date | stream | source | hits | after dedup | screened | included (in corpus) |", "|---|---|---|---|--:|--:|--:|--:|"]
    for r in log:
        t.append(f"| {r['search_id']} | {r['date']} | {r['stream']} | {r['source']} | {r['hits']} | {r['after_dedup']} | {r['screened']} | {r['included']} |")
    out["search-table"] = "\n".join(t) + "\n\nNR: not reconstructable from a file; `search_log.tsv` (`count_notes`) gives the value stated elsewhere and where. Exact query strings, filters and files of record are in `search_log.tsv`; screening criteria verbatim in `criteria.md`."
    bio = A["biolit_2026_04_04"]
    out["flow-april"] = "\n\n".join([
        mermaid([
            box("P", f"PubMed n = 87 (search note; hit list not saved)<br/>rerun {pm['rerun_date']}: {pm['rerun_hits_all_dates']}, {pm['rerun_hits_mesh_indexed_by_2026_04_04']} indexed by 2026-04-04"),
            box("E", f"Embase (Ovid) n = {emb['records']}<br/>{emb['with_pmid']} with PMID, {emb['without_pmid']} without"),
            box("W", f"Web of Science n = {wos['records']}"),
            box("D", f"Records before deduplication n = {dd['records_before']}<br/>unique by PMID/DOI {dd['unique_by_pmid_or_doi']}; by PMID/DOI/title {dd['unique_by_pmid_doi_or_title']}<br/>(computed 2026-09-27; not done in April)"),
            box("B", f"biolit criterion screen, 2026-04-04: {bio['screened_rows']} screening rows<br/>{bio['pass_rows']} pass rows = {bio['pass_distinct_papers']} distinct papers<br/>+ Embase DOI rescreen {emb['doi_rescreen_screened']} screened, {emb['doi_rescreen_passed']} passed<br/>+ Bousman 2010 {A['bousman']['passed']} passed"),
            box("Q", f"Screening pool (Paul, UI, 2026-04-16)<br/>{pool['records']} records = {pool['papers']} papers ({pool['duplicate_records']} duplicate records)"),
            box("X1", f"Excluded n = {pool['exclude_papers']} papers<br/>reason noted for {pool['exclude_records_with_note']} of {pool['exclude_records']} records"),
            box("I1", f"Included n = {pool['include_papers']} papers"),
            box("R", f"Second pass, full text (review app, 2026-04-29/30 and 2026-09-22)<br/>{sp['reviewed']} reviewed"),
            box("X2", f"Excluded n = {sp['exclude']} (reasons in prisma_16b.tsv)"),
            box("I2", f"Corpus 'Added Apr 2026' n = {FC['by_source'].get('Added Apr 2026', 0)}"),
            "P --> D", "E --> D", "W --> D", "D --> B", "B --> Q", "Q --> X1", "Q --> I1", "I1 --> R", "R --> X2", "R --> I2"]),
        (f"Merikangas's 26 (previous review): {A['m26']['in_any_april_file']} are in the April search files "
         f"(PubMed runs {A['m26']['in_pubmed_april_runs']}, Web of Science export {A['m26']['in_wos_export']}, Embase export {A['m26']['in_embase_export']}); "
         f"not in them: {', '.join(A['m26']['not_in_april_files']) or 'none'}; in the 2026-09-27 PubMed rerun: {A['m26']['in_pubmed_rerun']}. "
         f"Include papers by database (a paper can list several): first screen {fmt(pool['include_papers_by_source_database'])}; after the second pass {fmt(A['second_pass_include_by_database'])}. "
         f"Biolit pass rows per run: PubMed {pm['passed_run_m26']} + {pm['passed_run_new']}, Embase {emb['passed_run_1']} + {emb['passed_run_2']} (two runs over the same {emb['screened_distinct']} PMIDs), "
         f"Web of Science {wos['passed']} of {wos['screened']}. "
         f"The April PubMed runs screened {pm['screened_records']} PMIDs, {pm['screened_in_rerun']} of them in the rerun hit list; "
         f"{len(pm['screened_not_in_rerun'])} were not ({', '.join(x['pmid'] for x in pm['screened_not_in_rerun'])}), and "
         f"{len(pm['rerun_not_screened'])} rerun hits indexed by 2026-04-04 were not screened ({', '.join(x['pmid'] for x in pm['rerun_not_screened'])}).")])
    b = S["batches"]
    bw, br, gm, ci = (b.get(k, {}) for k in ("backward reference lists", "broader reference filter", "Gemma schizophrenia datasets", "cited as a data source"))
    readm = sum(v for k, v in FC["by_source_detail"].items() if "re-admitted" in k)
    out["flow-september"] = "\n\n".join([
        mermaid([
            box("B1", f"Reference lists, first filter<br/>{S['backward']['candidates_file']} candidates + {len(S['backward']['screened_not_in_candidates'])} by hand"),
            box("B2", f"biolit screen: {S['backward']['screened']} screened, {S['backward']['passed']} passed"),
            box("C1", f"Reference lists, broader filter<br/>{S['broader']['candidates_file']} candidates"),
            box("C2", f"biolit screen: {S['broader']['screened']} screened, {S['broader']['passed']} passed"),
            box("G1", f"Gemma datasets n = {S['gemma_check']['datasets']}<br/>({S['gemma_check']['via_uri_and_query']} MONDO-annotated)"),
            box("G2", f"papers behind unseen datasets n = {gm.get('include', 0) + gm.get('exclude', 0)}"),
            box("K1", f"Cited sample sources (sources_llm.py)<br/>expression studies n = {ci.get('include', 0) + ci.get('exclude', 0)}"),
            box("V", f"Review batch (Sonnet 5 eligibility fields + review app, 2026-09-22)<br/>{S['papers']} papers: {S['include']} include, {S['exclude']} exclude"),
            box("X", f"Excluded n = {S['exclude']}<br/>backward {bw.get('exclude', 0)}, broader {br.get('exclude', 0)}, Gemma {gm.get('exclude', 0)}, cited {ci.get('exclude', 0)}<br/>of which {readm} re-admitted later as one-gene studies"),
            box("I", f"Corpus 'Added Sep 2026 (reference lists)' n = {FC['by_source'].get('Added Sep 2026 (reference lists)', 0)}"),
            "B1 --> B2 --> V", "C1 --> C2 --> V", "G1 --> G2 --> V", "K1 --> V", "V --> X", "V --> I", "X -. re-admitted .-> I"]),
        "Composition of the corpus source 'Added Sep 2026 (reference lists)': " + "; ".join(
            f"{k.split(' | ', 1)[1]} {v}" for k, v in FC["by_source_detail"].items() if k.startswith("Added Sep 2026 (reference lists)")) + ". "
        + "Deciders, from the note tags in `review_decisions.json`: " + "; ".join(
            f"{k}: " + ", ".join(f"{kk.replace('_by_', ' by ')} {vv}" for kk, vv in v.items() if "_by_" in kk) for k, v in b.items()) + "."])
    fc = F["counts"]
    out["flow-forward"] = mermaid([
        box("F1", f"Seeds n = {F['seeds']}; citation links n = {fc.get('citation links')}<br/>non-MEDLINE links dropped {fc.get('non-MED links')}; unique citing PMIDs {fc.get('unique citing PMIDs')}"),
        box("F2", f"Candidate pool n = {fc.get('candidate pool')}<br/>(already in corpus {fc.get('already in corpus')})"),
        box("F3", f"Regex prefilter pass n = {fc.get('prefilter pass')}"),
        box("F4", f"Haiku screen: include {fc.get('include')}, unsure {fc.get('unsure')}, exclude {fc.get('exclude')}"),
        box("F5", f"Sonnet second pass: include {F['sonnet'].get('include')}, unsure {F['sonnet'].get('unsure')}, exclude {F['sonnet'].get('exclude')}"),
        box("F6", f"Added without hand review n = {F['added']}; removed at gene-list step {F['removed_at_gene_list_step']}"),
        box("F7", f"Corpus 'Added Sep 2026 (forward citations)' n = {FC['by_source'].get('Added Sep 2026 (forward citations)', 0)}"),
        "F1 --> F2 --> F3 --> F4 --> F5 --> F6 --> F7"])
    out["flow-single-gene"] = mermaid([
        box("S1", f"PubMed: {G['keyword_queries']} keyword + {G['author_queries']} author queries; citation finds {sum(G['other_sources'].values())}<br/>union n = {G['union']} (count-only query not screened: {sum(G['count_only'].values())})"),
        box("S2", f"Pool n = {G['pool']} (in corpus {G['already_in_corpus']}, already decided {G['already_decided']})"),
        box("S3", f"Haiku: {fmt(G['haiku'])}<br/>Sonnet rescue of {G['rescue_sent']}: {fmt(G['rescue_answers'])}<br/>final: {fmt(G['final'])}"),
        box("S4", f"transcriptome/proteome-wide labels n = {G['genome_wide_rows']}<br/>Sonnet class: {fmt(G['genome_wide_classes'])}"),
        box("S5", f"Full-text verification (batch 10)<br/>transcriptome: {G['batch10_transcriptome_include']} include, {G['batch10_transcriptome_reject']} reject, {G['batch10_transcriptome_duplicate']} duplicate, {G['batch10_pending']} pending<br/>proteomics: {G['batch10_proteomics_include']} include, {G['batch10_proteomics_reject']} reject"),
        box("S6", f"Corpus: genome-wide sweep {FC['by_source'].get('Added Sep 2026 (genome-wide sweep)', 0)}, proteomics {FC['by_source'].get('Added Sep 2026 (proteomics)', 0)}"),
        box("S7", f"Single-gene includes n = {G['papers_json_includes']} (2026-10-01: {G['includes_by_batch_date'].get('2026-10-01', 0)}): H4 units, not corpus studies"),
        "S1 --> S2 --> S3 --> S4 --> S5 --> S6", "S3 --> S7"])
    out["flow-single-cell"] = mermaid([
        box("C1", f"PubMed query + vault-note DOIs<br/>candidates n = {SC['candidates']}"),
        box("C2", f"{SC['model']}: include {SC['model_include']}<br/>+ {SC['by_hand']} study and {SC['preprints_by_hand']} preprints by hand"),
        box("C3", "Where they sit in the lineage: " + ", ".join(f"{k} {v}" for k, v in SC["fate_in_lineage"].items())),
        "C1 --> C2 --> C3"])
    sa = H["source_audit"]
    fu, gw = SV["forward_update_2026_09_27"], SV["gemma_wide_2026_09_27"]
    out["flow-other"] = "\n\n".join([
        mermaid([
            box("H1", f"Cited-source audit 2026-09-25: {sa['checked']} items checked<br/>classes {fmt(sa['by_class'])}"),
            box("H2", f"decisions: {fmt(sa['by_decision'])}"),
            box("H3", f"Added by hand, corpus n = {FC['by_source'].get('Added by hand', 0)}: " + ", ".join(f"{k} {v}" for k, v in H["by_route"].items())),
            "H1 --> H2 --> H3"]),
        mermaid([
            box("U1", f"Forward update 2026-09-27: {fu['seeds']} seeds<br/>new citing papers {fu['delta_pool']} (prefilter pass {fu['delta_prefilter_pass']}) + {fu['old_exclude_pool']} earlier Haiku excludes"),
            box("U2", f"Sonnet class screen n = {fu['class_screened_rows']}: {fmt(fu['class'])}"),
            box("U3", f"Sonnet reuse screen: {fmt(fu['sonnet_reuse'])}"),
            box("U4", f"Kept by hand (Claude): new data {fu['kept_new_data_by_hand']}; reuse listed {fu['reuse_listed']} ({fu['reuse_from_reading_lists']} from reading lists)<br/>in the lineage: {fu['kept_new_data_in_lineage'] + fu['reuse_in_lineage']} (awaiting Paul)"),
            box("U5", f"Gemma-wide sweep 2026-09-27: {gw['inspected']} human datasets inspected, {gw['with_sz_samples']} with SZ samples, {gw['not_in_2026_09_25_survey']} new; none added"),
            "U1 --> U2 --> U3 --> U4"])])
    out["flow-other"] += ("\n\nLeads of 2026-09-27 that already carry an eligibility-stage exclusion in `prisma_16b.tsv`: "
                          + (", ".join(fu["leads_previously_excluded"]) or "none") + ".")
    out["discrepancies"] = (
        f"13. **Discrepancies between notes and files** found here: the report's Methods gives 668 records for April, the export files give "
        f"{dd['records_before']} (the search note was corrected on 2026-09-27); the search note gives 145 Embase records with a PMID and 32 without, the "
        f"export has {emb['with_pmid']} records with a PMID ({emb['with_pmid_distinct']} distinct) and {emb['without_pmid']} without; the search note "
        f"gives 168 Web of Science records screened, the run's artifacts hold {wos['screened']}; the tracker gives 26 rejects in batch 10, "
        f"`corpus/added_genomewide_2026-09-24.json` has {G['batch10_transcriptome_reject']} rejects, {G['batch10_transcriptome_duplicate']} duplicates "
        f"and {G['batch10_pending']} pending (Chen 2024, a worker reject held by Paul).")
    out["proposed-report-text"] = (
        f"Methods, \"Corpus\", first sentence, replacing \"(668 records)\":\n\n> The Merikangas PubMed, Embase and Web of Science queries rerun on "
        f"2026-04-04 without a date limit ({dd['records_before']} records, {dd['unique_by_pmid_doi_or_title']} after deduplication by identifier and title); "
        f"biolit screening; Paul's manual screen of {pool['records']} candidates; then backward and forward citation, cited-source, single-cell, Gemma and "
        f"genome-wide sweeps, each study checked against its paper. Every search, its query, date, counts and screening method, and a PRISMA 2020 flow to "
        f"the {FC['studies']} studies, are in `search/PRISMA.md` and `search/search_log.tsv`.")
    src = FC["by_source"]
    order = ["Merikangas 2022", "Added Apr 2026", "Added Sep 2026 (reference lists)", "Added Sep 2026 (forward citations)",
             "Added Sep 2026 (genome-wide sweep)", "Added Sep 2026 (proteomics)", "Single-cell / single-nucleus", "Added by hand"]
    lines = [box(f"N{i}", f"{s} n = {src.get(s, 0)}") for i, s in enumerate(order)] + [f"N{i} --> T" for i in range(len(order))]
    pv = FC["list_provenance"]
    lines += [
        box("T", f"Total studies in the corpus (corpus/lineage.json) n = {FC['studies']}<br/>groups: " + ", ".join(f"{k} {v}" for k, v in FC["by_group"].items())),
        box("L", f"Gene lists extracted n = {FC['gene_lists_extracted']} from {FC['studies_with_extracted_list']} studies<br/>published {pv.get('published', 0)}, Gemma {pv.get('gemma', 0)}, reproduced {pv.get('reproduced', 0)}"
            + (f", review compilation {pv['review_compilation']}" if pv.get('review_compilation') else "")),
        box("U", f"In-corpus studies with a primary list n = {FC['recurrence_primary_in_corpus']}<br/>(left out: meta-lane {FC['recurrence_left_out_meta_lane']}, multi-series {FC['recurrence_left_out_multi_series']}, other reanalyses {FC['recurrence_left_out_other_reanalyses']})"),
        box("PB", f"Published brain, study tier: {FC['published_brain_candidates']} candidates, {FC['published_brain_non_de_dropped']} non-DE lists dropped<br/>= {FC['published_brain_studies']} published brain studies"),
        box("CT", f"{FC['published_brain_in_clusters']} of {FC['published_brain_studies']} in {FC['cohort_clusters']} cohort clusters"),
        "T --> L --> U --> PB --> CT"]
    out["flow-combined"] = mermaid(lines) + "\n\n" + "\n".join(
        ["| corpus source (lineage.json `source`) | studies | route |", "|---|--:|---|"]
        + [f"| {k.split(' | ')[0]} | {v} | {k.split(' | ')[1]} |" for k, v in FC["by_source_detail"].items()]
        + [f"| **all** | **{FC['studies']}** | |"]) + (
        f"\n\nExtracted lists from studies outside the corpus: {', '.join(FC['extracted_studies_outside_corpus']) or 'none'}. "
        f"`analyses/summary.tsv` has {FC['summary_rows']} list rows for {FC['summary_studies']} studies; the rows not extracted record why (no list published, targeted, blocked and others).")
    cat = collections.Counter((r["stream"].split(",")[0], r["reason_category"]) for r in r16)
    streams = sorted({s for s, _ in cat})
    cats = sorted({c for _, c in cat})
    t = ["| reason category | " + " | ".join(streams) + " | all |", "|---|" + "--:|" * (len(streams) + 1)]
    for c in cats:
        t.append(f"| {c} | " + " | ".join(str(cat.get((s, c), "")) for s in streams) + f" | {sum(cat.get((s, c), 0) for s in streams)} |")
    t.append("| **all** | " + " | ".join(str(sum(v for (s, _), v in cat.items() if s == s2)) for s2 in streams) + f" | {len(r16)} |")
    out["16b-summary"] = "\n".join(t)
    return out


def apply_blocks(text, blk):
    for name, body in blk.items():
        pat = re.compile(rf"(<!-- prisma:begin {re.escape(name)} -->\n)(?:.*?\n)?(<!-- prisma:end {re.escape(name)} -->)", re.S)
        if not pat.search(text):
            raise SystemExit(f"PRISMA.md has no block '{name}'")
        text = pat.sub(lambda m: m.group(1) + body + "\n" + m.group(2), text)
    return text


# ---------------------------------------------------------------- main
def compute():
    INPUTS.clear()
    nodes = lineage_nodes()
    idx_all = node_index(nodes)
    C = {"april": stream_april(), "september_batch": stream_september_batch(), "forward": stream_forward(),
         "single_gene": stream_single_gene(), "single_cell": stream_single_cell(idx_all), "hand": stream_hand(),
         "surveys": stream_surveys(idx_all)}
    C["final"] = final_corpus(nodes, C["september_batch"])
    r16 = list_16b(C["april"], C["september_batch"], nodes)
    C["item_16b"] = {"rows": len(r16), "by_category": dict(sorted(collections.Counter(r["reason_category"] for r in r16).items()))}
    ex_ids = {r["id"].lower(): f"{r['stream']}, {r['stage'].split(' (')[0]}" for r in r16}
    pool, dec, groups = C["april"]["_pool_rows"]
    for g in groups:
        if {dec[pool[i]["id"]]["decision"] for i in g} != {"exclude"}:
            continue
        noted = any((dec[pool[i]["id"]].get("note") or "").strip() for i in g)
        for i in g:
            for k in keys_of(pool[i]):
                ex_ids[k.lower()] = "April 2026 databases, title/abstract " + ("(reason in prisma_16b.tsv)" if noted else "(no reason recorded)")
    fu = C["surveys"]["forward_update_2026_09_27"]
    prev = {}
    for p_, ks in fu.pop("_lead_keys").items():
        for k in sorted(ks):
            if k.lower() in ex_ids:
                prev.setdefault(p_, ex_ids[k.lower()])
    fu["leads_previously_excluded"] = [f"{p_} ({v})" for p_, v in sorted(prev.items(), key=lambda x: int(x[0]))]
    log = search_log(C)
    return C, log, r16


TSV_NOTE = ("# plain tab-separated, no quoting: query strings are verbatim, double quotes included; read with quoting off "
            "(Python csv.QUOTE_NONE; pandas quoting=3, comment='#')\n")


def tsv(rows, cols):
    clean = lambda v: re.sub(r"[\t\r\n]+", " ", str(v))
    return "\t".join(cols) + "\n" + "".join("\t".join(clean(r.get(k, "")) for k in cols) + "\n" for r in rows)


def read_plain_tsv(path):
    lines = [l.rstrip("\n") for l in Path(path).read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def render(C, log, r16):
    clean = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in C.items()}
    files = {}
    files[HERE / "prisma_counts.json"] = json.dumps(clean, indent=1, ensure_ascii=False) + "\n"
    files[HERE / "search_log.tsv"] = ("# generated by search/build_prisma.py; do not edit by hand. NR = not reconstructable from a file (see count_notes).\n"
                                      + TSV_NOTE + tsv(log, LOG_COLS))
    cols = ["id", "label", "title", "stream", "stage", "reason_category", "reason", "decided_by", "file_of_record"]
    files[HERE / "prisma_16b.tsv"] = ("# PRISMA 2020 item 16b: studies excluded at the eligibility stage with a recorded reason (generated by search/build_prisma.py).\n"
                                      "# reason_category is a keyword rule over the verbatim reason (REASON_RULES in the script); the reason column is the file's own text.\n"
                                      + TSV_NOTE + tsv(r16, cols))
    files[HERE / "criteria.md"] = criteria_md() + "\n"
    files[PRISMA_MD] = apply_blocks(PRISMA_MD.read_text(encoding="utf-8"), blocks(C, log, r16))
    return files


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="exit 1 if any output would change; write nothing")
    a = ap.parse_args(argv)
    C, log, r16 = compute()
    files = render(C, log, r16)
    stale = [p for p, t in files.items() if not p.exists() or p.read_text(encoding="utf-8") != t]
    if a.check:
        for p in stale:
            print(f"stale: {_rel(p)}")
        return 1 if stale else 0
    for p, t in files.items():
        p.write_text(t, encoding="utf-8")
    meta = {"script": "search/build_prisma.py", "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "outputs": {_rel(p): hashlib.sha256(t.encode()).hexdigest() for p, t in files.items() if p != PRISMA_MD},
            "inputs": dict(sorted(INPUTS.items()))}
    (HERE / "build_prisma_meta.json").write_text(json.dumps(meta, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(files)} files; corpus {C['final']['studies']} studies, {C['final']['gene_lists_extracted']} lists, "
          f"{C['final']['published_brain_studies']} published brain studies; {len(log)} searches; {len(r16)} item-16b rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
