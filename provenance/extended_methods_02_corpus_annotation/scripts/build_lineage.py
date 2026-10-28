#!/usr/bin/env python3
"""Build the study-lineage graph for the Merikangas extension.

Reads the included studies and their relation / parent-dataset fields, and
writes lineage.json (nodes, edges, lanes) plus lineage.html (the page).

Inputs (all on disk, nothing hardcoded):
  - biolit-results/merikangas-eligibility-backward/results.csv   (73 new, 2026-09-22)
  - biolit-results/merikangas-eligibility-corpus/run_*/results_*.csv  (Merikangas 26 + April includes)
  - schizophrenia-transcriptome-replication/review-app/backward-2026-09-22/review_decisions.json
    (decisions, and Paul/Claude-reviewed relation + parent for the new 73)
  - schizophrenia-transcriptome-replication/review-app/review_decisions.json  (April decisions)
  - lineage_sets.json (which PMIDs are Merikangas's 26 / April includes)
  - milestones.json (grants and consortium announcements; drawn as the page's bottom lane, not written to lineage.json)
PubMed metadata is fetched once and cached in pubmed_cache.json.
"""
from __future__ import annotations

import csv
from collections import Counter
import glob
import hashlib
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
from vault import RA  # vault Research Assistant/ (project moved to Dev, 2026-09-25)
BACK_CSV = RA / "biolit-results/merikangas-eligibility-backward/results.csv"
CORPUS_CSVS = sorted(glob.glob(str(RA / "biolit-results/merikangas-eligibility-corpus/run_*/results_*.csv")))
BACK_DEC = RA / "schizophrenia-transcriptome-replication/review-app/backward-2026-09-22/review_decisions.json"
SETS = HERE / "lineage_sets.json"
HAND = HERE / "added_by_hand.json"          # corpus studies added by hand (reanalyses the search missed)
FORWARD = HERE / "added_forward.json"       # forward-citation sweep includes (forward_sweep/, Sonnet screen)
MILESTONES = HERE / "milestones.json"      # grants and consortium announcements, drawn as the page's bottom lane only
GENOMEWIDE = HERE / "added_genomewide_2026-09-24.json"  # genome-wide studies from single_gene_sweep/, checked against each paper
# proteome-wide SZ vs control studies from the same sweep (Paul, 2026-09-24: "proteomics is fine"; own analysis modality)
PROTEOMICS = HERE / "added_proteomics_2026-09-24.json"
# both built by build_added_batch10.py from single_gene_sweep/batch10_verification/verify_*.json
EXTRA_ROWS = [(HAND, "Added by hand"), (FORWARD, "Added Sep 2026 (forward citations)"),
              (GENOMEWIDE, "Added Sep 2026 (genome-wide sweep)"), (PROTEOMICS, "Added Sep 2026 (proteomics)")]
STANLEY = HERE / "stanley_datasets.json"    # the 12 SMRIDB array studies (Higgs 2006; Gemma stanley_*)
AUTHORS = HERE / "authors.json"             # PubMed author lists, fallback for papers not in Zotero
BIB = Path.home() / "Zotero" / "My Library.bib"
CACHE = HERE / "pubmed_cache.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SOURCES = json.loads((HERE / "sources_llm.json").read_text())["papers"] if (HERE / "sources_llm.json").exists() else {}

# Cohort / brain-bank lanes, matched against the parent-dataset text. Order is
# priority: the first match wins, so specific collections come before generic ones.
LANES = [
    ("pitt", "Pittsburgh (Lewis / Mirnics)", r"pittsburgh|mirnics|conte center|ccnmd|lewis lab"),
    ("stanley", "Stanley (SMRI) collections", r"stanley|smri"),
    ("libd", "LIBD / BrainSeq / NIMH CBDB", r"lieber|libd|brainseq|clinical brain disorders|cbdb|weinberger"),
    ("cmc", "CommonMind / PsychENCODE / HBCC", r"commonmind|\bcmc\b|psychencode|hbcc|nih neurobiobank|phs000979"),
    ("mssm", "Mount Sinai (Haroutunian / Davis)", r"mount sinai|mssm|haroutunian"),
    ("harvard", "Harvard / McLean", r"harvard|mclean|benes|hbtrc"),
    ("aus", "Australian banks (NSW TRC, Victorian, ASRB)", r"nsw|new south wales|sydney|asrb|australian|victorian|vbbn|florey"),
    # the Pritzker consortium's brains come from the UC Irvine brain bank (Hagenauer 2018 Fig. 3 "Brain Banks: UC-Irvine";
    # Zeppillo 2022 "University of California Irvine Brain Bank (UCIBB)")
    ("pritzker", "Pritzker consortium", r"pritzker|irvine brain bank|uc ?irvine brain|\bucibb\b|irvine brain donor"),
    ("ipsc", "iPSC-derived cells and organoids", r"brennand|coriell|ipsc|organoid"),  # label: 9 of 29 are Brennand-lab papers (2026-09-24)
]
OTHER_BRAIN = ("brain_other", "Other brain samples")
OTHER_PERIPH = ("periph", "Blood, LCL, fibroblast, other")
META_LANE = ("meta", "Meta-analyses and combined analyses")

REL_NORMAL = [
    ("primary", "primary"), ("data reuse", "data reuse"), ("single-gene", "single-gene follow-up"),
    ("meta-analysis", "meta-analysis / combined"), ("reanalysis", "reanalysis of public data"),
]
BRAIN_WORDS = r"cortex|cortical|dlpfc|prefrontal|hippocamp|amygdala|cerebell|striat|caudate|putamen|accumbens|thalam|gyrus|brain|entorhinal|cingulate|neuron"


def norm_rel(text: str) -> str:
    t = (text or "").strip().lower()
    for key, val in REL_NORMAL:
        if t.startswith(key):
            return val
    return "primary"


def read_rows() -> dict[str, dict]:
    rows = {}
    for path in [str(BACK_CSV), *CORPUS_CSVS]:
        with open(path) as f:
            for r in csv.DictReader(f):
                for key in (r.get("pmid"), r.get("doi")):
                    if key and key not in ("null", "None") and key not in rows:
                        rows[key] = r
    for path, _ in EXTRA_ROWS:
        for k, r in (json.loads(path.read_text()).items() if path.exists() else []):
            # a row checked against the paper (build_added_batch10.py, 'verification') replaces a model-read eligibility row:
            # Hathy 2020 and Ma 2024, excluded in the backward review and included in the 2026-09-30 iPSC round, kept the
            # backward CSV's text and lost their verification
            if not k.startswith("_") and (k not in rows or (r.get("verification") and not rows[k].get("verification"))):
                rows[k] = r
    return rows


def fetch_pubmed(ids: list[str]) -> dict[str, dict]:
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    todo = [i for i in ids if i not in cache]
    # DOIs → PMIDs
    for i in [x for x in todo if not x.isdigit()]:
        q = urllib.parse.quote(f"{i}[doi]")
        r = json.load(urllib.request.urlopen(f"{EUTILS}esearch.fcgi?db=pubmed&retmode=json&term={q}"))
        found = r["esearchresult"]["idlist"]
        cache[i] = {"pmid": found[0]} if found else {"pmid": None}
        time.sleep(0.4)
    pmids = [x for x in todo if x.isdigit()] + [cache[x]["pmid"] for x in todo if not x.isdigit() and cache[x].get("pmid")]
    for start in range(0, len(pmids), 150):
        batch = pmids[start:start + 150]
        if not batch:
            continue
        data = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(batch)}).encode()
        x = urllib.request.urlopen(f"{EUTILS}efetch.fcgi", data=data).read().decode()
        for a in x.split("<PubmedArticle>")[1:]:
            g = lambda p: (m.group(1) if (m := re.search(p, a, re.S)) else "")
            pm = g(r"<PMID[^>]*>(\d+)")
            clean = lambda s: html.unescape(re.sub("<[^>]+>", "", s)).strip()
            cache[pm] = {
                "pmid": pm,
                "title": clean(g(r"<ArticleTitle>(.*?)</ArticleTitle>")),
                "year": int(g(r"<PubDate>.*?<Year>(\d{4})") or g(r"<Year>(\d{4})") or 0),
                "first_author": g(r"<LastName>([^<]+)</LastName>"),
                "journal": clean(g(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>")),
                "doi": g(r'<ArticleId IdType="doi">([^<]+)'),
            }
        time.sleep(0.4)
    CACHE.write_text(json.dumps(cache, indent=1, ensure_ascii=False))
    out = {}
    for i in ids:
        pm = i if i.isdigit() else cache.get(i, {}).get("pmid")
        if pm and pm in cache:
            out[i] = cache[pm]
    return out


def parse_n(text: str) -> int | None:
    """Total samples from the model's tissue_and_n text: every count attached to a
    group word (cases, controls, subjects, pairs x2, 'SZ n=35', ...), summed."""
    t = (text or "").replace(",", "")
    groups = r"(?:matched\s+|drug-\w+\s+|minimally\s+medicated\s+)?(sib-pairs?|pairs?|schizo\w*|SZ|SCZ|BD|MDD|cases?|patients?|controls?|subjects?|individuals|unaffected|healthy|comparison|normal|persons|donors)"
    total = 0
    # 'SZ n=35' counts first and is blanked out, so its number is not read again as '35 controls' by the next pattern
    # ('SZ n=35, controls n=30' is 65, not 100)
    n_eq = re.compile(r"\b(SZ|SCZ|BD|MDD|schizophrenia|controls?|comparison|cases?)\s*\(?n\s*=\s*(\d{1,4})", re.I)
    for word, num in n_eq.findall(t):
        total += int(num)
    t = n_eq.sub(" ", t)
    for num, word in re.findall(r"\b(\d{1,4})\s+" + groups, t, re.I):
        total += int(num) * (2 if word.lower().endswith(("pair", "pairs")) else 1)
    if not total:
        m = re.search(r"(\d{1,4})\s*(?:/|vs\.?|versus)\s*(\d{1,4})", t)
        total = int(m.group(1)) + int(m.group(2)) if m else 0
    return total if 0 < total < 5000 else None


def lane_for(parent: str, tissue: str, rel: str, cited: str = "") -> tuple[str, str]:
    if rel == "meta-analysis / combined":
        return META_LANE
    for text, skip in ((parent or "", ()), (cited or "", ("ipsc",))):
        for key, label, pat in LANES:
            if key not in skip and re.search(pat, text, re.I):
                return key, label
    return OTHER_BRAIN if re.search(BRAIN_WORDS, tissue or "", re.I) else OTHER_PERIPH


SOURCE_PREFIXES = ("ext:", "smri:", "coll:")


# Tissue classes for the "Tissue" colouring, matched against tissue_n, then the title.
# Matched against tissue_n, then the title; the class mentioned first wins, since
# "iPSC-derived neurons" and "olfactory neuronal cells" also match brain words.
TISSUES = [
    ("ipsc", r"ipsc|induced pluripotent|hipsc|pluripotent|organoid|\bnpcs?\b|neural progenitor|stem cell"),
    ("lcl", r"lymphoblast|\blcls?\b|\blb cells"),
    ("other", r"fibroblast|skin|olfactory|nasal|turbinate|hair follicle|epitheli|saliva|muscle"),
    ("brain", BRAIN_WORDS + r"|post-?mortem|\bba\s?\d+|brodmann|white matter|gr[ae]y matter"),
    ("blood", r"blood(?![- ]brain)|pbmc|lymphocyte|leukocyte|monocyte|\bt cells?\b|serum|plasma|whole-blood"),
]


CSF_FIRST = r"\s*(csf\b|cerebrospinal)"


def tissue_class(n: dict) -> str:
    """brain / blood / lcl / ipsc / other / mixed (brain and periphery) / unknown.
    The class named first in the text wins; mixed when brain and a peripheral tissue are both named.
    A text that starts with CSF is 'other' (Harrington 1985 was 'unknown'; 2026-10-01). CSF named later, as in
    'postmortem brain (CSF used only for Western validation)', does not make a study mixed."""
    if re.match(CSF_FIRST, n.get("tissue_n") or "", re.I):
        return "other"
    for text in (n.get("tissue_n") or "", n.get("title") or ""):
        t = text.lower()
        pos = {k: m.start() for k, pat in TISSUES if (m := re.search(pat, t))}
        if not pos:
            continue
        first = min(pos, key=pos.get)
        # 'olfactory neuronal cells', 'olfactory epithelium neuronal layers': neurons of a peripheral tissue, not brain
        if first == "other" and "brain" in pos and not re.search(BRAIN_WORDS.replace("|neuron", "") + r"|post-?mortem", t):
            return "other"
        if first in ("brain", "blood", "other") and "brain" in pos and ({"blood", "other"} & set(pos)):
            return "mixed"
        return first
    if n.get("lane") == "ipsc":
        return "ipsc"
    return "unknown"


# Short tissue for the cards: the class plus the regions or cell types the eligibility text names.
REGIONS = [
    ("DLPFC", r"dlpfc|dorsolateral|\bba ?(9|46)\b|area 9\b"), ("sgACC", r"sgacc|subgenual"), ("ACC", r"\bacc\b|anterior cingulate"),
    ("OFC", r"orbitofrontal|\bofc\b"), ("frontal pole", r"frontal pole|\bba ?10\b"), ("PFC", r"prefrontal|\bpfc\b"),
    ("frontal cortex", r"(?<!pre)frontal cortex"), ("temporal cortex", r"temporal (cortex|gyrus)|\bstg\b|\bba ?2[12]\b"),
    ("entorhinal", r"entorhinal"), ("parietal", r"parietal"), ("visual cortex", r"visual cortex|occipital"),
    ("insula", r"insula"), ("hippocampus", r"hippocamp|\bca1\b|dentate"), ("amygdala", r"amygdala"),
    ("striatum", r"striatum|caudate|putamen|accumbens|\bnacc?\b"), ("thalamus", r"thalam"), ("habenula", r"habenula"),
    ("midbrain", r"midbrain|substantia nigra"), ("cerebellum", r"cerebell"), ("white matter", r"white matter"),
    ("choroid plexus", r"choroid plexus"), ("SEZ", r"subependymal"),
    ("PBMC", r"pbmc|peripheral blood mononuclear"), ("whole blood", r"whole[- ]blood"), ("lymphocytes", r"lymphocyte"),
    ("monocytes", r"monocyte"), ("leukocytes", r"leukocyte"), ("plasma/serum", r"plasma|serum|exosom|extracellular vesicle"),
    ("CSF", r"^" + CSF_FIRST), ("fibroblasts", r"fibroblast"), ("olfactory", r"olfactory|nasal|turbinate"), ("hair follicle", r"hair follicle"),  # turbinate: Tung 2023 (2026-10-01)
    ("neurons", r"neuron"), ("NPCs", r"\bnpcs?\b|neural progenitor"), ("organoids", r"organoid"),
    ("astrocytes", r"astrocyte"), ("oligodendrocytes", r"oligodendrocyte|\bopcs?\b"), ("microglia", r"microglia"),
]
TISSUE_WORD = {"brain": "Brain", "blood": "Blood", "lcl": "LCL", "ipsc": "iPSC", "other": "Peripheral", "mixed": "Brain + periphery",
               "unknown": "Not stated"}


def tissue_short(n: dict) -> str:
    """'Brain: DLPFC, hippocampus', 'Blood: PBMC', 'iPSC: neurons, NPCs', or the class alone."""
    t = (n.get("tissue_n") or "").lower()
    found = [name for name, pat in REGIONS if re.search(pat, t)]
    if not found:  # no region in the eligibility text: use the Step 2 list contrasts (e.g. Maycox 2009, BA10 and BA9)
        t = " ".join(x.get("contrast") or "" for x in (n.get("genelist") or {}).get("lists", [])).lower()
        found = [name for name, pat in REGIONS if re.search(pat, t)]
    if "sgACC" in found and "ACC" in found:
        found.remove("ACC")
    if {"DLPFC", "frontal pole", "OFC"} & set(found) and "PFC" in found:
        found.remove("PFC")
    word = TISSUE_WORD.get(n.get("tissue"), "Not stated")
    return f"{word}: {', '.join(found[:4])}{'…' if len(found) > 4 else ''}" if found else word


# Brain tissue split into neocortex / other region (Paul, 2026-09-25: "For tissue, brain stratified by neocortex and other
# brain region or something?"). Display only: n["tissue"] stays "brain", which the lane rules and analyses read.
# Entorhinal cortex counts with the medial temporal lobe (other region), as Paul's grouping puts medial temporal there.
NEOCORTEX = {"DLPFC", "sgACC", "ACC", "OFC", "frontal pole", "PFC", "frontal cortex", "temporal cortex", "parietal",
             "visual cortex", "insula", "neocortex", "motor cortex", "cingulate"}
OTHER_REGION = {"entorhinal", "hippocampus", "amygdala", "striatum", "thalamus", "habenula", "midbrain", "cerebellum",
                "white matter", "choroid plexus", "SEZ", "hypothalamus", "brainstem", "locus coeruleus"}
EXTRA_REGIONS = [("neocortex", r"neocort"), ("motor cortex", r"motor cortex"), ("cingulate", r"cingulate"),
                 ("hypothalamus", r"hypothalam"), ("brainstem", r"brain ?stem|\bpons\b|medulla"), ("locus coeruleus", r"locus co?eruleus")]
# Brodmann areas: the region each belongs to (neocortex except 27, 28, 34, 35, 36: medial temporal)
BA_REGION = {9: "DLPFC", 46: "DLPFC", 10: "frontal pole", 11: "OFC", 47: "OFC", 24: "ACC", 32: "ACC", 25: "sgACC", 33: "ACC",
             21: "temporal cortex", 22: "temporal cortex", 20: "temporal cortex", 38: "temporal cortex", 41: "temporal cortex",
             42: "temporal cortex", 17: "visual cortex", 18: "visual cortex", 19: "visual cortex", 39: "parietal", 40: "parietal",
             7: "parietal", 5: "parietal", 8: "frontal cortex", 6: "motor cortex", 4: "motor cortex", 44: "frontal cortex",
             45: "frontal cortex", 23: "cingulate", 31: "cingulate", 13: "insula", 28: "entorhinal", 34: "entorhinal",
             27: "entorhinal", 35: "entorhinal", 36: "entorhinal"}


def brain_regions(n: dict) -> tuple[str, str, list[str], str]:
    """(class, display, regions, source text) for a brain study. class: neocortex / other region / neocortex and other /
    unclassified. Regions come from the eligibility text (tissue_n), then the Step 2 list contrasts; generic 'cortex' or
    'brain' with no named area is unclassified (not guessed)."""
    sources = [("tissue_n", n.get("tissue_n") or "")] + [
        ("list contrasts", " ".join(x.get("contrast") or "" for x in (n.get("genelist") or {}).get("lists", [])))]
    for where, text in sources:
        t = text.lower()
        found = [name for name, pat in REGIONS + EXTRA_REGIONS if re.search(pat, t)]
        bas = sorted({int(b) for b in re.findall(r"\b(?:ba|brodmann(?:'s)? area|area)\s?(\d{1,2})\b", t) if 0 < int(b) < 53})
        for b in bas:
            if BA_REGION.get(b) and BA_REGION[b] not in found:
                found.append(BA_REGION[b])
        found = [f for f in found if f in NEOCORTEX | OTHER_REGION]
        if "sgACC" in found and "ACC" in found:
            found.remove("ACC")
        if {"DLPFC", "frontal pole", "OFC"} & set(found) and "PFC" in found:
            found.remove("PFC")
        if "neocortex" in found and len(found) > 1:
            found.remove("neocortex")
        if not found:
            continue
        neo, oth = bool(set(found) & NEOCORTEX), bool(set(found) & OTHER_REGION)
        cls = "neocortex and other" if neo and oth else "neocortex" if neo else "other region"
        ba_of = {}
        for b in bas:
            ba_of.setdefault(BA_REGION.get(b), []).append(f"BA{b}")
        disp = ", ".join(f + (f" ({', '.join(ba_of[f])})" if f in ba_of else "") for f in found)
        return cls, disp, found, f"{where}: {text[:300]}"
    return "unclassified", "", [], f"tissue_n: {(n.get('tissue_n') or '')[:300]}"


def enrich_collections(nodes: list[dict], edges: list[dict]) -> None:
    """Replace the bare '<consortium> (data resource)' placeholder with the curated
    description in collections.json, and link the resource to its papers that are in
    the corpus ('resource paper' edges). Papers outside the corpus are listed on the node."""
    path = HERE / "collections.json"
    if not path.exists():
        return
    cur = json.loads(path.read_text())
    by_id = {n["id"]: n for n in nodes}
    by_pmid = {n.get("pmid"): n for n in nodes if not n.get("external") and n.get("pmid")}
    for cid, c in cur.items():
        node = by_id.get(cid)
        if not isinstance(c, dict) or node is None:
            continue
        # The diamond sits at the first year a corpus study used the resource (or its data paper in the corpus, below); the
        # curated year (PsychENCODE 2015, the project start) is kept as "founded". At 2015, PsychENCODE sat before the
        # CommonMind and LIBD data it includes (first used 2018; Paul, 2026-09-27).
        node.update({"label": c["label"], "title": c["title"], "year": node["year"] or c.get("year"), "founded": c.get("year"),
                     "portal": c.get("portal", ""), "resource_source": c.get("source", "")})
        node["papers"] = []
        for p in c.get("papers", []):
            hit = by_pmid.get(p["pmid"])
            node["papers"].append({**p, "in_corpus": bool(hit)})
            if hit and hit.get("year") and hit["year"] < (node["year"] or 9999):
                node["year"] = hit["year"]  # its own data paper in the corpus (Fromer 2016 for CommonMind) counts as a first use
            if hit and not any(e["from"] == cid and e["to"] == hit["id"] for e in edges):
                edges.append({"from": cid, "to": hit["id"], "kind": "resource paper", "evidence": p["role"]})


def drop_same_subjects(nodes: list[dict], edges: list[dict]) -> list[dict]:
    """'Same subjects' links came from a model reading which references a paper cites for its
    samples; they are sparse (most papers cite the collection, not a prior paper) and some are
    wrong. Dropped as a link type (Paul, 2026-09-23): subject overlap is to be established by
    identifiers, and collection-level sharing is carried by the lanes. The citation stays in
    the study's details under cited sample sources."""
    by_id = {n["id"]: n for n in nodes}
    kept = []
    for e in edges:
        if e["kind"] != "same subjects":
            kept.append(e)
            continue
        s, t = by_id.get(e["from"]), by_id.get(e["to"])
        if s and t:
            t.setdefault("cited_sources", []).append(
                {"label": s["label"], "title": s.get("title", ""), "pmid": s.get("pmid", ""), "kind": e["kind"], "evidence": e.get("evidence", "")})
    print(f"'same subjects' links dropped: {len(edges) - len(kept)}")
    return kept


def prune_sources(nodes: list[dict], edges: list[dict]) -> tuple[list[dict], list[dict], dict]:
    """Keep an outside source node only if >= 2 corpus studies link to it, or it is itself
    an expression dataset (the `dataset` flag, or class expression_data in sources_class.json).
    A dropped source is listed on each citing study as `cited_sources`, so nothing is lost."""
    cls_path = HERE / "sources_class.json"
    cls = {k: v["class"] for k, v in json.loads(cls_path.read_text())["items"].items()} if cls_path.exists() else {}
    ovr_path = HERE / "sources_class_overrides.json"
    cls.update(json.loads(ovr_path.read_text()) if ovr_path.exists() else {})
    by_id = {n["id"]: n for n in nodes}
    citing: dict[str, list[str]] = {}
    for e in edges:
        if e["from"].startswith(SOURCE_PREFIXES) and e["to"] in by_id and e["to"] not in citing.get(e["from"], []):
            citing.setdefault(e["from"], []).append(e["to"])
    ext = [n for n in nodes if n.get("external")]
    before = {"sources": len(ext), "edges": len(edges),
              "by_citing": {k: sum(len(citing.get(n["id"], [])) == k for n in ext) for k in sorted({len(citing.get(n["id"], [])) for n in ext})}}
    drop = set()
    for n in ext:
        ids = citing.get(n["id"], [])
        n["n_citing"] = len(ids)
        n["cited_by"] = [by_id[i]["label"] for i in ids]
        n["expression_data"] = bool(n.get("dataset")) or cls.get(n["id"]) == "expression_data"
        if len(ids) < 2 and not n["expression_data"]:
            drop.add(n["id"])
    for e in edges:
        if e["from"] in drop and e["to"] in by_id:
            s = by_id[e["from"]]
            by_id[e["to"]].setdefault("cited_sources", []).append(
                {"label": s["label"], "title": s.get("title", ""), "pmid": s.get("pmid", ""), "kind": e["kind"], "evidence": e.get("evidence", "")})
    nodes = [n for n in nodes if n["id"] not in drop]
    edges = [e for e in edges if e["from"] not in drop and e["to"] not in drop]
    kept = [n for n in nodes if n.get("external")]
    after = {"sources": len(kept), "edges": len(edges),
             "kept_shared": sum(n["n_citing"] >= 2 for n in kept),
             "kept_dataset_only": sum(n["n_citing"] < 2 for n in kept)}
    print(f"cited sources: {before['sources']} before (by number of citing studies: {before['by_citing']}), "
          f"{after['sources']} after ({after['kept_shared']} cited by 2+ studies, {after['kept_dataset_only']} expression datasets cited by one); "
          f"edges {before['edges']} -> {after['edges']}; {len(drop)} moved to the citing study's details")
    return nodes, edges, {"before": before, "after": after}


DEGENES = HERE.parent / "analyses"
DATA = DEGENES.resolve().parent / "data"  # repo data/ (lists, manifests, reference, raw)
sys.path.insert(0, str(DEGENES))  # analysis_modality, classify_technology
# A manifest note records full summary statistics only in a clause that says so and does not negate it.
FULLSTATS_POS = re.compile(r"full (?:all-gene )?summary statistics|full per-(?:probeset|gene)", re.I)
FULLSTATS_NEG = re.compile(r"\bnot\b|\bpartial\b", re.I)


def genelist_status() -> dict[str, dict]:
    """Per PMID: Step 2 status, list counts, rows, unique genes, count match, full summary statistics.
    Read only on analyses (summary.tsv, manifest.json, manifests/*.json, lists/*.tsv)."""
    if not (DEGENES / "summary.tsv").exists():
        return {}
    summary = {r["list"]: r for r in csv.DictReader(open(DEGENES / "summary.tsv"), delimiter="\t")}
    lists = []
    for p in [DATA / "manifest.json", *sorted((DATA / "manifests").glob("*.json"))]:
        lists += json.loads(p.read_text()).get("lists", [])

    def bucket(text: str) -> str:
        t = text.strip().lower()
        if t.startswith("extracted"):
            return "extracted"
        if t.startswith("blocked"):
            return "blocked"
        if t.startswith("excluded"):
            return "excluded"
        if t.startswith("pending") or not t:
            return "pending"
        return "no list published"

    out: dict[str, dict] = {}
    for l in lists:
        key = f"{l['pmid']}__{l['list_id']}"
        s = summary.get(key, {})
        status_text = s.get("status") or l.get("status") or ""
        b = bucket(status_text)
        if b == "extracted" and l.get("provenance") in ("reproduced", "gemma"):
            # recomputed by us, or Gemma's own reanalysis; never counted as a published list
            b = l["provenance"]
        rec = {"list_id": l["list_id"], "contrast": l.get("contrast", ""), "status": b,
               "reason": l.get("reproduction", "") if b in ("reproduced", "gemma") else ("" if b == "extracted" else status_text),
               "rows": int(s["n_rows"]) if s.get("n_rows", "").isdigit() else None,
               "genes": int(s["n_genes"]) if s.get("n_genes", "").isdigit() else None,
               "expected_n": int(s["expected_n"]) if s.get("expected_n", "").isdigit() else None,
               "matches": {"yes": True, "NO": False}.get(s.get("matches_paper", "")), "file": key,
               # list role (manifest "role"): first_look = the study's own new data; reanalysis = its analysis of data
               # first published elsewhere (Paul, 2026-09-24); bank: one bank's part of a multi-bank study
               "role": l.get("role", "first_look"), "bank": l.get("bank", "")}
        if l.get("subset_flag"):
            # a list that is a subset or a non-DE statistic although the study is a primary SZ vs control source
            # (manifest "subset_flag", analyses/lists_fix_2026-10-10/subset_flags.tsv)
            rec["subset_flag"] = l["subset_flag"]
        st = out.setdefault(l["pmid"], {"label": l.get("label", ""), "lists": [], "full_stats": False, "full_stats_source": []})
        st["lists"].append(rec)
        for clause in re.split(r"(?<=[.;])\s+", f"{l.get('notes', '')} {s.get('status', '') or l.get('status', '')}"):
            if FULLSTATS_POS.search(clause) and not FULLSTATS_NEG.search(clause):
                st["full_stats"] = True
                src = f"{l['list_id']}: {clause.strip()[:240]}"
                if not any(clause.strip()[:80] in x for x in st["full_stats_source"]):
                    st["full_stats_source"].append(src)
    for pm, st in out.items():
        ls = st["lists"]
        ext = [x for x in ls if x["status"] == "extracted"]
        genes = set()
        for x in ext:
            f = DATA / "lists" / f"{x['file']}.tsv"
            if f.exists():
                rows = [r for r in f.read_text().splitlines() if r and not r.startswith("#")]
                hdr = rows[0].split("\t") if rows else []
                if "symbol" in hdr:
                    i = hdr.index("symbol")
                    genes |= {r.split("\t")[i] for r in rows[1:] if len(r.split("\t")) > i and r.split("\t")[i]}
        stated = [x for x in ext if x["expected_n"] is not None]
        order = ["extracted", "reproduced", "gemma", "pending", "blocked", "excluded", "no list published"]
        st.update({
            "status": next(b for b in order if any(x["status"] == b for x in ls)),
            "n_lists": len(ls), "n_extracted": len(ext),
            "n_blocked": sum(x["status"] == "blocked" for x in ls),
            "n_reproduced": sum(x["status"] == "reproduced" for x in ls),
            "n_gemma": sum(x["status"] == "gemma" for x in ls),
            "rows": sum(x["rows"] or 0 for x in ext), "unique_genes": len(genes),
            "n_stated": len(stated), "n_match": sum(bool(x["matches"]) for x in stated),
        })
        for x in ls:
            x.pop("file")
    return out


CTX_SEQ = re.compile(r"rna-?seq|sequenc|psychencode", re.I)
CTX_ARR = re.compile(r"microarray|\barrays?\b|affymetrix|genechip|beadchip", re.I)
_TECH: dict = {}


def _tech_of(c: dict) -> str:
    """technology.json's technology for a corpus study, else the keyword rule on its platform text."""
    if not _TECH:
        p = DEGENES / "technology.json"
        _TECH.update(json.loads(p.read_text())["studies"] if p.exists() else {"_": {}})
    from classify_technology import technology_of_text
    t = (_TECH.get(c.get("pmid")) or {}).get("technology")
    return t if t in ("microarray", "rnaseq") else technology_of_text(c.get("platform"))


def pick_cited(cands: list[dict], context, cite: str = "") -> tuple[dict | None, str]:
    """One corpus study for a citation matched by first author and year. Several studies can share the label (two
    'Gandal 2018' papers); linking all of them draws a study as pooling data it did not use (tests: one citation links one
    study). First rule that applies: one candidate; a candidate's PMID in the citation text; the citation context names one
    technology (sequencing / PsychENCODE, or arrays), every candidate's technology is known (technology.json) and exactly one
    is of that technology;
    otherwise none, and the citation is kept on the citing study as an ambiguous source. context: a text, or several tried
    in order, narrowest first. (A word-overlap score was tried first and linked 2 of 4 Gandal 2018 citations wrongly.)"""
    if len(cands) <= 1:
        return (cands[0] if cands else None), ""
    contexts = [context] if isinstance(context, str) else list(context)
    for c in cands:
        if c.get("pmid") and re.search(rf"\b{c['pmid']}\b", f"{cite} {' '.join(contexts)}"):
            return c, f"PMID {c['pmid']} in the citation"
    for text in contexts:
        t = f"{cite} {text}"
        seq, arr = bool(CTX_SEQ.search(t)), bool(CTX_ARR.search(t))
        if seq != arr and all(_tech_of(c) in ("microarray", "rnaseq") for c in cands):  # the others are known to differ
            want = "rnaseq" if seq else "microarray"
            hit = [c for c in cands if _tech_of(c) == want]
            if len(hit) == 1:
                w = (CTX_SEQ if seq else CTX_ARR).search(t).group(0)
                return hit[0], f"the citation context names '{w}'; the only candidate of that technology"
    return None, "ambiguous: " + ", ".join(f"{c['label']} (PMID {c.get('pmid')})" for c in cands)


def isSC_node(n: dict) -> bool:
    return bool(n.get("single_cell")) or n.get("lane") == "sc"


def study_key(n: dict) -> str:
    """The Step 2 key of a node: its PMID, or doi_<doi with / as _> for a study with none (Cui 2024, preprints)."""
    return n.get("pmid") or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")


# A stated panel size in free text: "1,128 brain-relevant genes", "custom 1,128-gene", "676 ... probe sets", ">1,800 candidate genes"
PANEL_SIZE = re.compile(r"([>~]|about |over )?(\d[\d,]*|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)[\s-]+(?:[A-Za-z][\w-]*\s+){0,3}?"
                        r"(genes?|probe sets|miRNA probes|probes|clones|miRNAs|circRNAs|lncRNAs|transcripts|targets)\b", re.I)
# Non-coding lists (coordinator, 2026-09-25): scope by molecule class, not by lineage group. A list is miRNA / circRNA / lncRNA when its
# list id names the class, or when at least half its symbols are MIR* (miRNA). lncRNA only by list id: gene-level lists carry LINC* symbols.
NONCODING_ID = [("miRNA", re.compile(r"mirna", re.I)), ("circRNA", re.compile(r"circrna", re.I)), ("lncRNA", re.compile(r"lncrna", re.I))]
MIR_SYM = re.compile(r"(?!.*HG$)(MIR\d|MIRLET)")  # not the host genes (MIR9-1HG)
# a stated count of non-coding features; the number must stand alone ("GENCODE v19", "GRCh38", "GSE17612" are not counts)
FEATURE_COUNT = re.compile(r"(?<![\w.])([>~]|about |over )?(\d[\d,]*|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)[\s-]+"
                           r"(?:[A-Za-z][\w-]*\s+){0,3}?(miRNA probes|miRNAs|circRNAs|lncRNAs|probes)\b", re.I)
# the paper says the features were pre-selected: a qPCR follow-up of screen hits, or named candidates
PRESELECTED = re.compile(r"candidate|pre-?selected|taken to qPCR", re.I)


def feature_class(list_id: str, symbols) -> str | None:
    """miRNA / circRNA / lncRNA, or None for a gene (mRNA) list."""
    for cls, rx in NONCODING_ID:
        if rx.search(list_id):
            return cls
    syms = list(symbols)
    if syms and sum(bool(MIR_SYM.match(x)) for x in syms) * 2 >= len(syms):
        return "miRNA"
    return None
GW_SUBTYPE = {"platform": "platform", "fullstats": "full table", "rnaseq_default": "RNA-seq"}


def _panel_size(*texts: str) -> tuple[str, str] | None:
    """(N with unit, the text it was read from) from the first text that states a panel size."""
    for t in texts:
        m = PANEL_SIZE.search(t or "")
        if m:
            unit = m.group(3).lower().replace("mirna", "miRNA")
            unit = "genes" if unit == "gene" else unit
            return f"{(m.group(1) or '').strip()}{m.group(2)} {unit}", t
    return None


def _feature_count(*texts: str) -> tuple[str, str] | None:
    """(N, the text) for a non-coding list: the first stated count of features, as a number only."""
    for t in texts:
        m = FEATURE_COUNT.search(t or "")
        if m:
            return f"{(m.group(1) or '').strip()}{m.group(2)}", t
    return None


def list_scope(u: dict, man: dict, summ: dict, node: dict, study_u: list[dict], tech: str | None,
               feature: str | None = None) -> tuple[str, str, str]:
    """How much of the genome a list tested: (class g genome-wide / t targeted / u unknown, label, source of the label).
    Paul, 2026-09-25: "panels are so biased, almost everything is going to sound super interesting". First rule that applies:
      1. manifest panel_genes                                   -> targeted, N = len(panel_genes)
      1n. a miRNA / circRNA / lncRNA list (feature_class)       -> targeted when the text says the features were pre-selected
                                                                   (PRESELECTED), else class n "<class>-wide: N tested", whatever the
                                                                   study's lineage group; "genome-wide" is kept for gene (mRNA) lists
      2. universes.json type targeted                           -> targeted, N = n_tested (genes), else the gene set size,
                                                                   else a size stated in its evidence
      3. universes.json type platform / fullstats / rnaseq_default -> genome-wide
      4. manifest universe or modality_note beginning 'targeted' -> targeted, N stated there
      5. a Gemma reanalysis                                     -> genome-wide (full table: Gemma tests every probe)
      6. study in group 1m (multigene panel), list not a reanalysis -> targeted, N from the study's universes.json entry or a stated size
      7. study in group 1 / 3, single-cell, or a reanalysis list of a 1m study -> genome-wide (categorize.py: transcriptome-wide)
      8. otherwise                                              -> scope unknown
    """
    uni = man.get("universe") or summ.get("universe") or ""
    if man.get("panel_genes"):
        n = len(man["panel_genes"])
        return "t", f"targeted: {n:,} genes tested", f"manifest panel_genes ({n} genes); universe: '{uni}'"
    if feature:
        texts = [uni, u.get("evidence"), *(x.get("evidence") for x in study_u), node.get("platform")]
        n = _feature_count(*texts)
        src = f"{feature} list (list id or symbols; build_lineage.py::feature_class)" + (f"; count stated in '{n[1]}'" if n else "; no count stated") + \
            (f"; manifest universe: '{uni}'" if uni and not (n and n[1] == uni) else "")
        pre = next((t for t in texts if t and PRESELECTED.search(t)), None)
        if pre:
            return "t", f"targeted: {n[0] + ' ' + feature + 's' if n else 'panel size not stated'} tested".replace("panel size not stated tested", "panel size not stated"), \
                f"pre-selected features: '{pre}'; " + src
        return "n", f"{feature}-wide: {n[0]} tested" if n else f"{feature}-wide: size not stated", src
    if u.get("type") == "targeted":
        if u.get("n_tested") and u.get("n_tested_unit") == "genes":
            return "t", f"targeted: {u['n_tested']:,} genes tested", f"universes.json n_tested ({u['n_tested_source']})"
        if u.get("n_k"):
            return "t", f"targeted: {u['n_k']:,} genes tested", f"universes.json gene set size ({u['source']}); evidence: '{u['evidence']}'"
        s = _panel_size(u.get("evidence"), uni)
        if s:
            return "t", f"targeted: {s[0]} tested", f"stated panel size: '{s[1]}' (universes.json evidence / manifest universe)"
        return "t", "targeted: panel size not stated", f"universes.json type targeted; evidence: '{u.get('evidence')}'"
    if u.get("type") in GW_SUBTYPE:
        return "g", f"genome-wide ({GW_SUBTYPE[u['type']]})", f"universes.json type {u['type']}: {u.get('evidence')}"
    note = man.get("modality_note") or ""
    if uni.lower().startswith("targeted") or note.lower().startswith("targeted"):
        s = _panel_size(uni, note)
        return "t", f"targeted: {s[0]} tested" if s else "targeted: panel size not stated", f"manifest universe / modality_note: '{uni or note}'"
    if (summ.get("provenance") or man.get("provenance")) == "gemma":
        return "g", "genome-wide (full table)", f"Gemma reanalysis; manifest universe: '{uni}'"
    role = summ.get("role") or man.get("role") or "first_look"
    cat = "sc" if node.get("single_cell") or node.get("lane") == "sc" else node.get("category")
    if cat == "1m" and role != "reanalysis":
        for x in study_u:
            if x.get("type") == "targeted" and x.get("n_tested") and x.get("n_tested_unit") == "genes":
                return "t", f"targeted: {x['n_tested']:,} genes tested", f"group 1m study; universes.json n_tested of the study's primary list ({x['n_tested_source']})"
            if x.get("type") == "targeted" and x.get("n_k"):
                return "t", f"targeted: {x['n_k']:,} genes tested", f"group 1m study; universes.json gene set size of the study's primary list ({x['source']})"
        s = _panel_size(uni, *(x.get("evidence") for x in study_u), node.get("platform"))
        return "t", f"targeted: {s[0]} tested" if s else "targeted: panel size not stated", \
            "lineage group 1m, multigene panel (corpus/categorize.py)" + (f"; stated size: '{s[1]}'" if s else "")
    if cat in ("1", "3", "sc") or cat == "1m":
        sub = "RNA-seq" if cat == "sc" else {"microarray": "platform", "rnaseq": "RNA-seq", "single_cell": "RNA-seq"}.get(tech or "")
        why = {"1": "group 1, transcriptome-wide", "3": "group 3, re-use of genome-wide data", "sc": "single-cell study",
               "1m": "reanalysis list in a group 1m study"}[cat]
        return "g", f"genome-wide ({sub})" if sub else "genome-wide", \
            f"lineage {why} (corpus/categorize.py)" + (f"; technology.json: {tech}" if tech else "") + (f"; manifest universe: '{uni}'" if uni else "")
    return "u", "scope unknown", f"no panel size, universes.json type or corpus group (universes.json: {u.get('type') or 'no entry'}; manifest universe: '{uni}')"


def platform_titles() -> dict[str, str]:
    """GPL -> platform title from the GEO annotation headers in data/platforms (Affymetrix / Illumina words trimmed)."""
    import gzip
    out = {}
    for f in sorted((DATA / "platforms").glob("GPL*")):
        opener = gzip.open if f.suffix == ".gz" else open
        try:
            with opener(f, "rt", errors="replace") as fh:
                for line in fh:
                    if line.startswith("!Annotation_platform_title"):
                        t = re.sub(r"^\[[^\]]*\]\s*", "", line.split("=", 1)[1].strip())
                        t = re.sub(r"\s*\[transcript \(gene\) version\]|\s*\(Feature Number version\)|\bAffymetrix\s+|\s+Array$|\s+expression beadchip$", "", t, flags=re.I)
                        out.setdefault(f.name.split(".")[0], t.replace("™", "").strip())
                        break
                    if not line.startswith(("^", "!")):
                        break
        except OSError:
            continue
    return out


def headline(n: dict, prim: list[str], uni: dict, titles: dict) -> dict:
    """What a study card shows under its title (Paul, 2026-09-25: "sample size and gene coverage should be headline-adjacent.
    the number of DE genes they said"): genes measured by the primary list's platform or table (universes.json), and the
    DE gene count the paper states beside the count in our list file (summary.tsv via genelist). prim: primary list keys
    (primary_lists.json), '<study>__<list_id>'."""
    out = {}
    u = next((uni[k] for k in prim if k in uni), None)
    if u:
        t = u.get("type")
        if t == "platform" and u.get("n_membership"):
            names = " + ".join(titles.get(g, g) for g in re.split(r"[+|&]", str(u.get("source") or "")) if g)
            out["coverage"] = f"{names}: {u['n_membership']:,} genes measured"
        elif t == "fullstats" and u.get("n_k"):
            out["coverage"] = f"{u['n_k']:,} genes tested (full results table)"
        elif t == "rnaseq_default":
            out["coverage"] = "RNA-seq, whole transcriptome (genes tested not published)"
        elif t == "targeted":
            k = u.get("n_tested") if u.get("n_tested_unit") == "genes" else u.get("n_k")
            out["coverage"] = f"targeted panel: {k:,} genes" if k else "targeted panel, size not stated"
        elif t == "exclude" and u.get("reason"):
            out["coverage"] = u["reason"]
        if out.get("coverage"):
            out["coverage_by"] = f"universes.json {u.get('type')}: {u.get('evidence', '')}"[:300]
    ids = {k.split("__", 1)[1] for k in prim}
    ls = [x for x in (n.get("genelist") or {}).get("lists", []) if x["list_id"] in ids]
    if ls:
        stated = [x["expected_n"] for x in ls if x.get("expected_n") is not None]
        ours = [x["genes"] for x in ls if x.get("genes") is not None]
        src = {x["status"] for x in ls}
        rows = [x["rows"] for x in ls if x.get("rows") is not None]
        out.update({"de_rows": sum(rows) if len(rows) == len(ls) else None,
                    "de_match": all(x.get("matches") for x in ls) if all(x.get("matches") is not None for x in ls) else None,
                    "de_list": ", ".join(x["list_id"] for x in ls),
                    "de_stated": sum(stated) if len(stated) == len(ls) else None,
                    "de_ours": sum(ours) if len(ours) == len(ls) else None,
                    "de_source": "gemma" if src == {"gemma"} else "reproduced" if src == {"reproduced"} else "published"})
    return out


def gene_index(nodes: list[dict]) -> tuple[dict, dict, dict]:
    """Gene symbol -> the extracted DE lists that contain it, for the page's gene search (Paul, 2026-09-25:
    searching HCN1 found nothing). Read only on analyses: lists/*.tsv, summary.tsv, manifest.json + manifests/*.json,
    primary_lists.json, analysis_modality.json. Only lists of studies whose analysis modality is DE.

    Encoding, to keep the page small:
      lists:  one row per list, [study key, node id or "", list id, contrast, threshold, role, primary, provenance, groups,
                                 scope class, scope label, scope source]
              primary is 1 / 0, or null when primary_lists.json has no entry for the study;
              groups: the list's cell types / regions (the 'group' column), or [] when the list has none;
              scope: how much of the genome the list tested (list_scope): g genome-wide (gene lists), n miRNA / circRNA /
              lncRNA-wide, t targeted, u unknown
      genes:  symbol -> tokens <list index in base 36>[.<group index in base 36>...]<direction>
              direction U up, D down, M both across the symbol's rows, N not stated
    Row fields 12-13 of lists: the effect column's name in the source table (manifest columns.effect) and the list's
    effect note (the '# effect:' header line), so the page can label an effect value in the source's own terms.

    Also returns the per-gene statistics as they appear in each list file (Paul, 2026-09-25: "when you filter on a gene,
    the statistic for the gene shows up on the card for each data set"), written by main() to lineage_gene_stats.json
    and loaded by the page on the first gene search, to keep the page small:
      symbol -> "li,gi,effect,p,fdr;li,gi,effect,p,fdr;..." one entry per row of the list file: list index and group
      index in base 36 (group "" when the list has none), values as printed in the list (rounded to 3 significant
      digits; "" when the list has none). No unit conversion.
    """
    lists_dir = DATA / "lists"
    paths = {p: DEGENES / p for p in ("summary.tsv", "primary_lists.json", "analysis_modality.json")}
    paths["manifest.json"] = DATA / "manifest.json"
    paths |= {f"manifests/{p.name}": p for p in sorted((DATA / "manifests").glob("*.json"))}
    if not lists_dir.exists() or not all(p.exists() for p in paths.values()):
        return {}, {}, {}
    summary = {r["list"]: r for r in csv.DictReader(open(paths["summary.tsv"]), delimiter="\t")}
    manifest = {}
    for k, p in paths.items():
        if k.endswith(".json") and (k == "manifest.json" or k.startswith("manifests/")):
            for l in json.loads(p.read_text()).get("lists", []):
                manifest[f"{l['pmid']}__{l['list_id']}"] = l
    primary = {k: set(v.get("primary") or []) for k, v in json.loads(paths["primary_lists.json"].read_text())["studies"].items()}
    modality = json.loads(paths["analysis_modality.json"].read_text())["studies"]
    node_of = {study_key(n): n["id"] for n in nodes if not n.get("external") and not n.get("landmark")}
    node_by = {study_key(n): n for n in nodes if not n.get("external") and not n.get("landmark")}
    # tested scope of each list (list_scope): universes.json and technology.json, read when present
    opt = {p: DEGENES / p for p in ("universes.json", "technology.json") if (DEGENES / p).exists()}
    universes = json.loads(opt["universes.json"].read_text())["lists"] if "universes.json" in opt else {}
    tech = json.loads(opt["technology.json"].read_text())["studies"] if "technology.json" in opt else {}
    study_u: dict[str, list] = {}
    for x in universes.values():
        study_u.setdefault(x["study"], []).append(x)

    rows_out, genes, files, skipped, stats = [], {}, {}, {}, {}
    per_node: dict[str, dict] = {}
    for f in sorted(lists_dir.glob("*.tsv")):
        key = f.stem
        study, list_id = key.split("__", 1)
        mo = (modality.get(study) or {}).get("modality", "not classified")
        if (node_by.get(study) or {}).get("category") == "2":
            # one-gene-at-a-time studies are left off the page (page_view), so their lists stay out of the gene search;
            # page_view stops the build otherwise (first hit: Fillman 2014, group 1m -> 2 under the 20-gene panel
            # minimum, 2026-09-30)
            skipped["category 2"] = skipped.get("category 2", 0) + 1
            continue
        if mo != "DE":
            skipped[mo] = skipped.get(mo, 0) + 1
            continue
        st = (summary.get(key) or {}).get("status", "").split(":")[0].strip()
        if st not in ("extracted", "reproduced", "gemma"):
            # e.g. de Jong 2012's Gemma reanalysis, excluded for double-logged data: GRIN1 was found through it (Paul, 2026-09-25)
            skipped[f"status {st or 'missing'}"] = skipped.get(f"status {st or 'missing'}", 0) + 1
            continue
        text = f.read_text()
        files[f"lists/{f.name}"] = hashlib.sha256(text.encode()).hexdigest()
        body = [r for r in text.splitlines() if r and not r.startswith("#")]
        if len(body) < 2:
            continue  # header only: no genes at the list's threshold (several Gemma reanalyses)
        hdr = body[0].split("\t")
        i_sym, i_dir = hdr.index("symbol"), hdr.index("direction")
        i_grp = hdr.index("group") if "group" in hdr else None
        i_stat = [hdr.index(c) if c in hdr else None for c in ("effect", "p", "fdr")]
        effect_note = next((r[len("# effect:"):].strip() for r in text.splitlines() if r.startswith("# effect:")), "")
        row_stats: list[tuple] = []
        dirs: dict[str, set] = {}
        grps: dict[str, set] = {}
        groups: list[str] = []
        for r in body[1:]:
            c = r.split("\t")
            sym = c[i_sym].strip() if len(c) > i_sym else ""
            if not sym:
                continue
            dirs.setdefault(sym, set()).add(c[i_dir] if len(c) > i_dir else "")
            g = c[i_grp].strip() if i_grp is not None and len(c) > i_grp else ""
            if g:
                if g not in groups:
                    groups.append(g)
                grps.setdefault(sym, set()).add(groups.index(g))
            row_stats.append((sym, groups.index(g) if g else -1, *(_stat(c[i]) if i is not None and len(c) > i else "" for i in i_stat)))
        if not dirs:
            continue
        m, s = manifest.get(key, {}), summary.get(key, {})
        li = len(rows_out)
        scope = list_scope(universes.get(key, {}), m, s, node_by.get(study, {}), study_u.get(study, []),
                           (tech.get(study) or {}).get("technology"), feature_class(list_id, dirs))
        # a panel is under 500 genes (Paul, 2026-09-28: '"Panels" would be <500 genes'); a larger custom array is a
        # low-coverage array, not a panel
        m_n = re.search(r"([\d,]+) (?:genes|\w+s) tested", scope[1]) if scope[0] == "t" else None
        if m_n and int(m_n.group(1).replace(",", "")) >= 500:
            scope = ("g", f"low-coverage array: {m_n.group(1)} genes tested", scope[2] + "; 500 or more genes: an array, not a panel")
        rows_out.append([study, node_of.get(study, ""), list_id, m.get("contrast", ""), m.get("threshold", ""),
                         s.get("role") or m.get("role") or "first_look",
                         (1 if key in primary[study] else 0) if study in primary else None,
                         s.get("provenance") or m.get("provenance") or "published", groups, *scope,
                         _effect_label(m), effect_note or m.get("effect_note", "")])
        for sym, gi, e, p, q in row_stats:
            if e or p or q:
                stats.setdefault(sym, []).append(",".join([_b36(li), _b36(gi) if gi >= 0 else "", e, p, q]))
        tag = _b36(li)
        for sym, ds in dirs.items():
            ds = ds - {""}
            d = "N" if not ds else "M" if len(ds) > 1 else {"up": "U", "down": "D"}.get(next(iter(ds)), "N")
            g = "".join("." + _b36(x) for x in sorted(grps.get(sym, ())))
            genes[sym] = genes.get(sym, "") + tag + g + d  # the direction letter ends each token
        pn = per_node.setdefault(study, {"lists": 0, "genes": set()})
        pn["lists"] += 1
        pn["genes"] |= set(dirs)
    for n in nodes:
        pn = per_node.get(study_key(n)) if not n.get("external") and not n.get("landmark") else None
        if pn:
            n["gene_index"] = {"lists": pn["lists"], "genes": len(pn["genes"])}
    index = {"lists": rows_out, "genes": dict(sorted(genes.items()))}
    meta = {
        "script": "corpus/build_lineage.py::gene_index",
        "scope": "extracted Step 2 lists of studies whose analysis modality is DE (analysis_modality.json); "
                 "lists with no gene rows are left out",
        "n_lists": len(rows_out), "n_symbols": len(genes), "n_studies": len({r[0] for r in rows_out}),
        "n_studies_not_in_graph": len({r[0] for r in rows_out if not r[1]}),
        "n_list_files_skipped_not_DE": skipped,
        "n_entries": sum(len(re.findall(r"[UDMN]", v)) for v in genes.values()),
        "scope_rules": "build_lineage.py::list_scope (first rule that applies; row fields 9-11: class, label, source)",
        "n_lists_by_scope": {c: sum(r[9] == c for r in rows_out) for c in ("g", "n", "t", "u")},
        "n_lists_by_label": dict(sorted(Counter(r[10] for r in rows_out).items())),
        "targeted_lists": [{"list": f"{r[0]}__{r[2]}", "label": r[10], "source": r[11]} for r in rows_out if r[9] == "t"],
        "inputs_sha256": {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in (paths | opt).items()} | files,
    }
    meta["n_stat_rows"] = sum(len(v) for v in stats.values())
    return index, meta, {k: ";".join(v) for k, v in sorted(stats.items())}


def _stat(v: str) -> str:
    """a list file's effect / p / fdr cell, as printed, rounded to 3 significant digits when numeric"""
    v = v.strip()
    if v in ("", "NA", "nan", "NaN"):
        return ""
    try:
        return f"{float(v):.3g}"
    except ValueError:
        return v[:16].replace(",", " ").replace(";", " ")


def _effect_label(m: dict) -> str:
    """the source column(s) the list's effect was read from (manifest columns.effect), e.g. 'log2FoldChange'"""
    e = (m.get("columns") or {}).get("effect") if isinstance(m.get("columns"), dict) else None
    return " | ".join(map(str, e)) if isinstance(e, list) else str(e or "")


TOUR_GENES = {"example_rare": "HCN1", "example_recurrent": "MT1X"}  # Paul, 2026-09-25: the tour's gene-search card
TOUR_SCENARIO = "published, brain"


def tour_genes() -> dict:
    """Rows of analyses/recurrence_genes.tsv (scenario 'published, brain', study tier) for the tour's example genes."""
    p = DEGENES / "recurrence_genes.tsv"
    if not p.exists():
        return {}
    rows = {r["gene"]: r for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")
            if r["scenario"] == TOUR_SCENARIO}
    keep = ("n_lists_measured", "n_lists_observed", "expected", "p", "q_bh", "p_prior", "q_prior_bh", "de_prior_hitlists",
            "up", "down", "n_clusters_measured", "n_clusters_observed", "studies")
    out = {k: {"gene": g, "scenario": TOUR_SCENARIO, "row": ({f: rows[g][f] for f in keep} if g in rows else None)} for k, g in TOUR_GENES.items()}
    out["_source"] = f"analyses/recurrence_genes.tsv sha256 {hashlib.sha256(p.read_bytes()).hexdigest()[:12]}"
    return out


# The page's "brain studies with their own published gene list, one per shared dataset": the representatives of analyses/recurrence_units.tsv, scenario 'published, brain', modelled 'yes'.
PUB_BRAIN_SCENARIO, PUB_BRAIN_N = "published, brain", 62  # 62 since 2026-10-10 (Matthews 2012 and Medina 2023 left the corpus); 64 since 2026-10-04/05: batch-13 proteomics + re-analysis studies entered the population (tracker: 59 bulk + 5 single-cell)


def published_brain_list_ids(nodes: list[dict]) -> dict:
    """Node ids of the page's published-brain-list set, with the source file's sha256. Page only (not in lineage.json).
    Stops the build when the set is not PUB_BRAIN_N studies or a representative has no node."""
    p = DEGENES / "recurrence_units.tsv"
    reps = [r["representative"] for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")
            if r["scenario"] == PUB_BRAIN_SCENARIO and r["modelled"] == "yes"]
    by_key = {study_key(n): n["id"] for n in nodes if not n.get("external")}
    missing = [k for k in reps if k not in by_key]
    if len(reps) != PUB_BRAIN_N or len(set(reps)) != len(reps) or missing:
        raise SystemExit(f"published_brain_list_ids: {len(reps)} rows ({len(set(reps))} distinct) in {p.name} scenario "
                         f"'{PUB_BRAIN_SCENARIO}' modelled 'yes', expected {PUB_BRAIN_N}; no node for {missing}")
    return {"ids": [by_key[k] for k in reps], "source": "analyses/recurrence_units.tsv",
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest()[:12]}


# One-gene papers from the single-gene sweep: the includes, as
# analyses/single_gene_study_agreement/h4_extract.py counts them (INCLUDE_BATCHES there); most are not corpus studies. From 2026-09-30 they
# include the corpus group-2 studies and the review-app single-gene decisions (added by the sweep without a screen).
SWEEP_PAPERS = HERE / "single_gene_sweep" / "papers.json"
SWEEP_INCLUDE = {"single-gene sweep: haiku include", "single-gene sweep: rescue include",
                 "single-gene: corpus group 2", "single-gene: review-app single-gene decision",
                 "single-gene: forward-cite search 2026-10-04"}


def single_gene_sweep(nodes: list[dict]) -> dict | None:
    """Count of the sweep's included papers and the node ids of those that are also corpus studies. Page only (not in
    lineage.json). Stops the build when the count differs from the rows of analyses/single_gene_study_agreement/h4_papers.tsv."""
    if not SWEEP_PAPERS.exists():
        return None
    papers = [x for x in json.loads(SWEEP_PAPERS.read_text()) if x.get("batch") in SWEEP_INCLUDE]
    h4 = DEGENES / "h4" / "h4_papers.tsv"
    if h4.exists():
        n_h4 = sum(1 for l in open(h4) if l.strip() and not l.startswith("#")) - 1  # less the header
        if n_h4 != len(papers):
            raise SystemExit(f"single_gene_sweep: {len(papers)} included papers in {SWEEP_PAPERS.name}, {n_h4} rows in {h4.name}")
    pm = {str(x.get("pmid")) for x in papers if x.get("pmid")}
    return {"papers": len(papers), "in_lineage": [n["id"] for n in nodes if not n.get("external") and str(n.get("pmid")) in pm],
            "source": "corpus/single_gene_sweep/papers.json", "sha256": hashlib.sha256(SWEEP_PAPERS.read_bytes()).hexdigest()[:12]}


_NOTE_PATH_MARKERS = re.compile(
    r"|gene-lists/|data/(?:raw|manifests)/|scratchpad/|analyses/|corpus/|side-studies/"
    r"|batch\d+_verification/|Zotero storage [A-Z0-9]{5,}"
)
_NOTE_EVIDENCE_SENT = re.compile(r"\s*Evidence:.*?\.(?=\s|\Z)")
_NOTE_TRAILING_CLAUSE = re.compile(
    r"[,;]\s*[^,;.|]*(?:|data/raw/|data/manifests/|gene-lists/|scratchpad/|analyses/|corpus/|side-studies/)[^,;.|]*"
)
_NOTE_FIELDS = ("category_note", "verification", "parent", "provisional_why", "tech_tier_by", "lane_by")
_UNUSED_FIELDS = ("field_override", "checked_by")  # internal audit trail; lineage_template.html never reads them


def _strip_path_parens(s: str) -> str:
    """Drops a parenthesised span (balanced, one or more levels) if it contains a repo path or Zotero storage id."""
    out, i, n = [], 0, len(s)
    while i < n:
        if s[i] == "(":
            depth, j = 1, i + 1
            while j < n and depth:
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                j += 1
            span = s[i:j]
            if _NOTE_PATH_MARKERS.search(span):
                i = j
                continue
            out.append(span)
            i = j
            continue
        out.append(s[i])
        i += 1
    return "".join(out)


def clean_note(s: str | None, max_segments: int = 2, max_chars: int = 320) -> str | None:
    """Strips repo/Zotero file paths out of an internal provenance note for display on the page (Paul, 2026-10-01:
    the visualizer has to be user-oriented, not a tool for us; also trims it, these notes run long)."""
    if not s:
        return s
    s = _strip_path_parens(s)
    s = _NOTE_EVIDENCE_SENT.sub("", s)
    s = _NOTE_TRAILING_CLAUSE.sub("", s)
    m = _NOTE_PATH_MARKERS.search(s)
    if m:
        s = s[: m.start()]
    s = re.sub(r"\s*\|\s*", " | ", s).strip(" |")
    segs = s.split(" | ")
    if len(segs) > max_segments:
        segs = segs[:max_segments]
    s = " | ".join(segs)
    s = re.sub(r"\s{2,}", " ", s).strip()
    s = re.sub(r"\($", "", s).rstrip()
    s = re.sub(r"\s+([.;:,])", r"\1", s)
    if len(s) > max_chars:
        cut = s.rfind(" ", 0, max_chars)
        cut = cut if cut > 200 else max_chars
        s = s[:cut].rstrip(" ,;|").rstrip() + "…"
    return s


def page_view(out: dict) -> dict:
    """lineage.json as the page gets it: without the one-gene-at-a-time studies (category 2), their links, the source
    nodes left with no link, so the page neither draws nor counts them (decided 2026-09-30). Stops if any of them has
    gene-index lists. Also with repo/Zotero file paths stripped and provenance notes shortened in the fields the page
    shows (category_note, verification, parent, provisional_why, edge evidence), and the fields the page never reads
    (field_override, checked_by) dropped -- lineage.json keeps the originals."""
    drop = {n["id"] for n in out["nodes"] if not n.get("external") and n.get("category") == "2"}
    edges = [e for e in out["edges"] if e["from"] not in drop and e["to"] not in drop]
    linked = {e["from"] for e in edges} | {e["to"] for e in edges}
    nodes = [n for n in out["nodes"] if n["id"] not in drop and not (n.get("external") and n["id"] not in linked)]
    pmids = {n.get("pmid") for n in out["nodes"] if n["id"] in drop}
    gi = out.get("gene_index")
    if gi and any(r[0] in pmids for r in gi.get("lists", [])):
        raise SystemExit("page_view: the gene index has lists from category-2 studies; drop them from 'genes' too before publishing")
    def clean_node(n):
        n = {**n, **{f: clean_note(n[f]) for f in _NOTE_FIELDS if n.get(f)}}
        if n.get("other_data"):
            n["other_data"] = [{**o, "why": clean_note(o["why"])} if o.get("why") else o for o in n["other_data"]]
        return {k: v for k, v in n.items() if k not in _UNUSED_FIELDS}

    nodes = [clean_node(n) for n in nodes]
    edges = [{**e, "evidence": clean_note(e["evidence"])} if e.get("evidence") else e for e in edges]
    view = {**out, "nodes": nodes, "edges": edges}
    if (out.get("stats") or {}).get("corpus"):
        from corpus_stats import corpus_counts
        view["stats"] = {**out["stats"], "corpus": corpus_counts(nodes)}
    return view


def _b36(i: int) -> str:
    s, digits = "", "0123456789abcdefghijklmnopqrstuvwxyz"
    while True:
        i, r = divmod(i, 36)
        s = digits[r] + s
        if not i:
            return s


SMRI = {d["short_name"]: d for d in json.loads(STANLEY.read_text())["datasets"]} if STANLEY.exists() else {}
_SMRI_META = json.loads(STANLEY.read_text()) if STANLEY.exists() else {}
OTHER_MEASURE = r"costa|cotter|fatemi|harrison|manji|o'banion|reynolds|uranova|density|immunocyto|protein levels"


# Brain-bank heads: an author from one of these lists places a study in that bank's lane.
# Vawter and Bunney are left out: they publish on Stanley, Pritzker and blood samples alike.
BANK_HEADS = {
    # UC Irvine brain bank / Pritzker (Paul, 2026-09-26: "constants" of the collection's work). Bunney WE and Vawter MP
    # are the UCI Department of Psychiatry authors on the collection's papers; Hagenauer MH submitted GSE92538 and is
    # first author of its release paper (Hagenauer 2018)
    "pritzker": {"Bunney WE", "Vawter MP", "Hagenauer MH"},
    "mssm": {"Haroutunian V", "Davis KL"},
    "pitt": {"Lewis DA", "Mirnics K"},
    "libd": {"Weinberger DR", "Kleinman JE", "Hyde TM"},
    "stanley": {"Webster MJ", "Torrey EF"},
    "aus": {"Dean B", "Weickert CS", "Shannon Weickert C"},
    "harvard": {"Benes FM"},
}


def author_lists(pmids: set[str], dois: dict[str, str]) -> dict[str, set[str]]:
    """'Last Initials' per PMID: Zotero's BibTeX export first, PubMed (authors.json) for the rest."""
    out: dict[str, set[str]] = {}
    if BIB.exists():
        by_pmid, by_doi = {}, {}
        for e in BIB.read_text().split("\n@"):
            au = re.search(r"\n\s*author = \{(.*?)\},\n", e, re.S)
            if not au:
                continue
            names = set()
            for a in re.sub(r"\s+", " ", au.group(1)).split(" and "):
                last, _, first = a.partition(",")
                ini = "".join(w[0] for w in re.split(r"[\s.\-]+", first.strip().replace("{", "").replace("}", "")) if w and w[0].isalpha()).upper()
                names.add(f"{last.strip().strip('{}')} {ini}")
            pm = re.search(r"\n\s*pmid = \{(\d+)\}", e)
            doi = re.search(r"\n\s*doi = \{([^}]*)\}", e)
            if pm:
                by_pmid[pm.group(1)] = names
            if doi:
                by_doi[doi.group(1).lower()] = names
        for p in pmids:
            hit = by_pmid.get(p) or by_doi.get((dois.get(p) or "").lower())
            if hit:
                out[p] = hit
    if AUTHORS.exists():
        pm_auth = json.loads(AUTHORS.read_text())["authors"]
        for p in pmids - set(out):
            if p in pm_auth:
                out[p] = {f"{a['last']} {a['initials']}" for a in pm_auth[p]}
    return out


def stanley_key(pm: str, ds: dict) -> tuple[str, str]:
    """Map a Stanley mention from pooled_*.json to one of the 12 SMRIDB array studies.
    Returns (kind, key): ("dataset", short_name) | ("extra", smri id) | ("collection", name) | ("other measure", name)."""
    name = (ds.get("name") or "").lower()
    num = str(ds.get("id") or "")
    if re.search(OTHER_MEASURE, name):
        return "other measure", name
    ov = _SMRI_META.get("mention_overrides", {}).get(f"{pm}|{name.split()[0] if name else ''}")
    if ov and not re.search(r"altar ?[ac]\b|collection [ac]\b|sklar ?[ab]\b", name):
        return "dataset", ov["short_name"]
    if num in _SMRI_META.get("choi_study_numbers", {}):
        return "dataset", _SMRI_META["choi_study_numbers"][num]
    for pat, short in [(r"altar ?a\b|collection a\b", "stanley_altarA"), (r"altar ?c\b|collection c\b", "stanley_altarC"),
                       (r"altar ?b\b|collection b\b", "stanley_altarB"), (r"sklar ?b\b", "stanley_sklarB"), (r"sklar", "stanley_sklarA"),
                       (r"bahn", "stanley_bahn"), (r"\bchen\b", "stanley_chen"), (r"dobrin", "stanley_dobrin"), (r"f[ei]{2}nberg", "stanley_feinberg"),
                       (r"\bkato\b", "stanley_kato"), (r"vawter", "stanley_vawter"), (r"\byoung\b", "stanley_young")]:
        if re.search(pat, name):
            return "dataset", short
    if re.search(r"array collection", name):
        return "collection", "Stanley Array Collection"
    if re.search(r"neuropathology consortium|consortium", name):
        return "collection", "Stanley Neuropathology Consortium"
    if num:
        return "extra", f"study{num}"
    return "collection", "Stanley collections"


CORPUS_EXCLUSIONS = {k: v for k, v in json.loads((HERE / "corpus_exclusions.json").read_text()).items() if not k.startswith("_")} \
    if (HERE / "corpus_exclusions.json").exists() else {}


def main() -> None:
    sets = json.loads(SETS.read_text())
    m26, april = set(sets["m26"]), set(sets["april"])
    back_dec = json.loads(BACK_DEC.read_text())
    rows = read_rows()

    cat_path = HERE / "study_categories.json"
    CATS = json.loads(cat_path.read_text())["studies"] if cat_path.exists() else {}
    # Studies added by hand after categorize.py last ran have no study_categories.json entry, and categorize.py applies
    # category_overrides.json only to studies it has an entry for; so a hand-added study's override (e.g. the batch 14
    # proteomic studies' scope 'proteome-wide', group p) never reached its node, which fell back to its relation
    # (found 2026-09-29: nine batch 14 proteomic studies counted as group 1). Take the override with categorize.assign.
    from categorize import assign as _assign
    for k, v in json.loads((HERE / "category_overrides.json").read_text()).items():
        if not k.startswith("_") and k not in CATS and isinstance(v, dict) and "scope" in v and "new_data" in v:
            CATS[k] = dict(v, pmid=k, category=_assign(v))
            CATS[k]["new_validation"] = CATS[k]["category"] == "3" and bool(v.get("new_data"))
    included = [k for k in m26 | april if k in rows] + [k for k, v in back_dec.items() if v.get("decision") == "include"]
    # Group-2 (one-gene-at-a-time) studies are part of the scheme.
    included += [k for k, v in CATS.items() if v.get("category") == "2" and k not in included and k in rows]
    extra_src = {k: label for path, label in EXTRA_ROWS if path.exists()
                 for k in json.loads(path.read_text()) if not k.startswith("_")}
    included += [k for k in sorted(extra_src) if k not in included]
    meta = fetch_pubmed(included)

    ncounts = json.loads((HERE / "n_extracted.json").read_text())["counts"] if (HERE / "n_extracted.json").exists() else {}

    fo_path = HERE / "node_field_overrides.json"  # hand corrections of model-extracted fields (2026-09-25)
    FIELD_OV = {k: v for k, v in (json.loads(fo_path.read_text()).items() if fo_path.exists() else []) if not k.startswith("_")}
    def n_for(key: str, text: str) -> tuple[int | None, int | None, int | None]:
        fo_ = FIELD_OV.get(key) or {}  # a hand / paper-check correction wins over the model extraction
        if fo_.get("cases") is not None or fo_.get("controls") is not None:
            # n only when both sides are known. Until 2026-09-29 a missing side counted as 0, so an override giving one
            # side made n equal that side (Garbett 2008 showed n 12, its 12 cases; decisions_paul_2026-09-28.json)
            ca, co = fo_.get("cases"), fo_.get("controls")
            return (ca + co if ca is not None and co is not None else None), ca, co
        c = ncounts.get(key)
        if c and (c.get("cases") or c.get("controls")):
            return (c.get("cases") or 0) + (c.get("controls") or 0), c.get("cases"), c.get("controls")
        return parse_n(text), None, None

    nodes = []
    for key in included:
        r, md = rows.get(key, {}), meta.get(key, {})
        # N is read from the model's own text: a correction that rewrites tissue_n without cases / controls is not an N
        # verdict, and parse_n summed every count in the fuller corrected text (brain + blood cohorts, validation sets,
        # a count quoted twice: Glatt 2005 46 -> 100, Huang 2020 426 -> 852), which also blocked the list-contrast N.
        n_text = r.get("tissue_and_n", "")
        fo = FIELD_OV.get(key)
        if fo:
            r = {**r, "tissue_and_n": fo.get("tissue_n", r.get("tissue_and_n", "")), "platform": fo.get("platform", r.get("platform", ""))}
        if not md.get("title"):  # no PubMed record (e.g. Cui 2024, Nat Ment Health): metadata from the hand-added row
            md = {**md, **{k: r[k] for k in ("first_author", "year", "title", "journal") if r.get(k)}}
        dec = back_dec.get(key, {})
        rel = dec.get("relation") or norm_rel(r.get("relation", ""))
        parent = dec.get("parent") or ("" if (r.get("parent_dataset", "").upper().startswith("NONE")) else r.get("parent_dataset", ""))
        cited = " ".join(f"{x.get('cite', '')} {x.get('title', '')}" for x in SOURCES.get(md.get("pmid", key), {}).get("sources", []))
        lane, lane_label = lane_for(parent, r.get("tissue_and_n", ""), rel, cited)
        source = "Merikangas 2022" if key in m26 else ("Added Apr 2026" if key in april else
                                                        extra_src.get(key, "Added Sep 2026 (reference lists)"))
        nodes.append({
            "id": key, "pmid": md.get("pmid", key if key.isdigit() else ""),
            "label": f"{html.unescape(md.get('first_author', '?'))} {md.get('year', '')}",
            "year": md.get("year") or 0, "title": md.get("title") or r.get("title", ""),
            "journal": md.get("journal", ""), "relation": rel, "lane": lane, "lane_label": lane_label,
            "parent": parent, "tissue_n": r.get("tissue_and_n", ""), "n": n_for(key, n_text)[0], "cases": n_for(key, n_text)[1], "controls": n_for(key, n_text)[2],
            "platform": r.get("platform", ""), "genes": r.get("de_genes_quote", ""), "source": source,
            "parent_by": "normalized by Claude" if dec.get("parent") else "model text",
        })
        if fo:
            for f in ("cases", "controls"):
                if f in fo:
                    nodes[-1][f] = fo[f]
            if fo.get("cases") is not None and fo.get("controls") is not None:
                nodes[-1]["n"] = fo["cases"] + fo["controls"]
            nodes[-1]["field_override"] = fo["why"]
        if r.get("verification"):
            # batch 10 rows (build_added_batch10.py): checked against the paper; provisional when the Methods were not read
            nodes[-1].update({"verification": r["verification"], "provisional": bool(r.get("provisional")),
                              "provisional_why": r.get("provisional_why", ""), "parent_by": "worker verification (bank_cohort)"})

    # Edges: parent text naming another corpus paper by first author (+ year when given).
    by_author: dict[str, list[dict]] = {}
    for n in nodes:
        by_author.setdefault(n["label"].split()[0].lower(), []).append(n)
    edges = []
    for n in nodes:
        for m in re.finditer(r"([A-Z][a-zA-Z\-]+)(?: et al\.?)?,? \(?((?:19|20)\d\d)\)?", n["parent"]):
            cands = [c for c in by_author.get(m.group(1).lower(), [])
                     if c["id"] != n["id"] and str(c["year"]) == m.group(2) and c["year"] <= n["year"]]
            cand, why = pick_cited(cands, n["parent"][max(0, m.start() - 60):m.end() + 60], m.group(0))
            if cand:
                edges.append({"from": cand["id"], "to": n["id"], "kind": "parent", **({"match": why} if why else {})})
            elif cands:
                n.setdefault("cited_sources", []).append({"label": m.group(0), "title": "", "pmid": "", "kind": "parent (ambiguous citation)", "evidence": why})
    # Edges from the papers each study cites as its sample/data source (sources_llm.py),
    # and from GEO series named in its parent text (geo_links.py). A cited source outside
    # the corpus becomes a hollow "source" node, placed in the lane of the studies citing it.
    by_pmid = {n["pmid"]: n for n in nodes if n.get("pmid")}
    external: dict[str, dict] = {}
    src_path, geo_path = HERE / "sources_llm.json", HERE / "geo_links.json"
    sources = json.loads(src_path.read_text())["papers"] if src_path.exists() else {}
    geo = json.loads(geo_path.read_text()) if geo_path.exists() else {"gse_pmids": {}, "paper_gse": {}}

    def link(src_pmid: str, cite: str, title: str, year, n: dict, how: str, evidence: str) -> None:
        if not src_pmid or src_pmid == n.get("pmid"):
            return
        if re.search(r"genome-wide association|associated genetic loci|schizophrenia loci|alleles are enriched|linkage|genome-wide scan", title or "", re.I):
            return  # GWAS / linkage papers are genotype sources, not expression data
        if src_pmid in by_pmid:
            edges.append({"from": by_pmid[src_pmid]["id"], "to": n["id"], "kind": how, "evidence": evidence})
            return
        e = external.setdefault(src_pmid, {
            "id": "ext:" + src_pmid, "pmid": src_pmid, "label": cite, "year": int(year or 0), "title": title,
            "journal": "", "relation": "source", "external": True, "cited_by": [], "lanes": [],
            "parent": "", "tissue_n": "", "n": None, "platform": "", "genes": "", "source": "Cited source (not in corpus)",
            "parent_by": ""})
        e["cited_by"].append(n["label"]); e["lanes"].append(n["lane"])
        if how == "brain collection description":
            e["collection"] = True
        edges.append({"from": e["id"], "to": n["id"], "kind": how, "evidence": evidence})

    for n in nodes:
        for s in sources.get(n.get("pmid"), {}).get("sources", []):
            m = re.search(r"((?:19|20)\d\d)", s.get("cite", ""))
            link(s.get("pmid"), s.get("cite", ""), s.get("title", ""), m.group(1) if m else 0, n, s.get("how", ""), s.get("evidence", ""))
        for g in dict.fromkeys(re.findall(r"GSE\d+", n["parent"]) + geo.get("text_gse", {}).get(n.get("pmid"), [])):
            for pm in geo["gse_pmids"].get(g, [])[:1]:
                cached = json.loads(CACHE.read_text()).get(pm, {}) if CACHE.exists() else {}
                link(pm, f"{g}", cached.get("title", f"GEO {g}"), cached.get("year", 0), n, "public dataset", f"{g} named as data source")
    # Datasets pooled by meta-analyses, transcribed from their tables (pooled_datasets.json).
    pooled_path = HERE / "pooled_datasets.json"
    pooled = {k: v for k, v in json.loads(pooled_path.read_text()).items() if not k.startswith("_")} if pooled_path.exists() else {}
    llm_path = HERE / "pooled_llm.json"
    for k, v in (json.loads(llm_path.read_text())["papers"].items() if llm_path.exists() else []):
        if k not in pooled and v.get("datasets"):
            pooled[k] = {"source": f"model reading of {v.get('where') or 'the paper'} ({v.get('text_source')})",
                         "relation": v.get("relation"), "datasets": v["datasets"]}
    CONSORTIUM_LANE = [(r"cmc|commonmind", "cmc", "CMC"), (r"psychencode", "cmc", "PsychENCODE"), (r"hbcc", "cmc", "NIMH HBCC"),
                       (r"libd|brainseq|lieber", "libd", "LIBD / BrainSeq"), (r"stanley.*array", "stanley", "Stanley Array Collection"),
                       (r"stanley.*neuropath|stanley.*consortium", "stanley", "Stanley Neuropathology Consortium"),
                       (r"stanley", "stanley", "Stanley collections"), (r"pritzker", "pritzker", "Pritzker")]
    # Only these case/control collections become nodes; any other consortium the model names
    # (GWAS consortia, GTEx/BrainSpan/Allen reference atlases, mouse phenotyping resources) is dropped.
    gemma_stanley = geo.get("gemma_stanley", {})
    stanley_nodes: dict[str, dict] = {}
    for pm, entry in pooled.items():
        n = by_pmid.get(pm)
        if not n:
            continue
        n["relation"] = entry.get("relation", n["relation"])
        if n["relation"] == "meta-analysis / combined":
            n["lane"], n["lane_label"] = META_LANE
        n["pooled_source"] = entry.get("source", "")
        for ds in entry["datasets"]:
            # an ArrayExpress copy of a GEO series (E-GEOD-21138 = GSE21138) is that GEO series: the model filed these as
            # 'paper' with the accession as the cite, so they linked to nothing (Bharadhwaj 2023: 7 of its 10 datasets)
            m_ae = re.search(r"E-GEOD-(\d+)", " ".join(str(ds.get(k) or "") for k in ("gse", "cite", "name")))
            if m_ae and ds["kind"] != "geo":
                ds = dict(ds, kind="geo", gse="GSE" + m_ae.group(1))
            ev = f"{entry.get('source', '').split(',')[0]}: {ds.get('name') or ds.get('gse') or ds.get('cite', '')} ({ds.get('n', '')})"
            if ds["kind"] == "geo":
                for src in geo["gse_pmids"].get(ds["gse"], [])[:1]:
                    cached = json.loads(CACHE.read_text()).get(src, {}) if CACHE.exists() else {}
                    link(src, ds.get("cite") or ds["gse"], cached.get("title", ds["gse"]), cached.get("year", 0), n, "pooled study", ev)
            elif ds["kind"] == "consortium":
                nm = ds.get("name", "")
                hit = next(((l, lab) for pat, l, lab in CONSORTIUM_LANE if re.search(pat, nm, re.I)), None)
                if not hit or re.search(r"gwas|genotyp|snp|exome", nm, re.I):
                    continue
                lane_id, label = hit
                cid = "coll:" + re.sub(r"\W+", "_", label.lower())
                c = stanley_nodes.setdefault(cid, {
                    "id": cid, "pmid": "", "label": label, "year": 0, "title": f"{label} (data resource)", "journal": "",
                    "relation": "source", "external": True, "collection": True, "cited_by": [], "lanes": [lane_id], "lane": lane_id,
                    "lane_label": "", "parent": "", "tissue_n": "", "n": None, "platform": "", "genes": "",
                    "source": "Cited source (not in corpus)", "parent_by": ""})
                c["cited_by"].append(n["label"])
                c["year"] = min([y for y in (c["year"], n["year"]) if y] or [n["year"]])
                edges.append({"from": cid, "to": n["id"], "kind": "pooled study", "evidence": ev})
            elif ds["kind"] == "stanley":
                kind, key = stanley_key(pm, ds)
                if kind == "other measure":
                    # Stanley neuropathology measures (cell counts, proteins): not expression data
                    n.setdefault("stanley_other", []).append(ds.get("name", ""))
                    continue
                if kind == "collection":
                    cid = "coll:" + re.sub(r"\W+", "_", key.lower())
                    c = stanley_nodes.setdefault(cid, {
                        "id": cid, "pmid": "", "label": key, "year": 0, "title": f"{key} (data resource)", "journal": "",
                        "relation": "source", "external": True, "collection": True, "cited_by": [], "lanes": ["stanley"], "lane": "stanley",
                        "lane_label": "", "parent": "", "tissue_n": "", "n": None, "platform": "", "genes": "",
                        "source": "Cited source (not in corpus)", "parent_by": ""})
                    c["cited_by"].append(n["label"])
                    c["year"] = min([y for y in (c["year"], n["year"]) if y] or [n["year"]])
                    edges.append({"from": cid, "to": n["id"], "kind": "pooled study", "evidence": ev})
                    continue
                sd = SMRI.get(key)
                paper = next((p_ for p_ in (sd or {}).get("describing_papers", []) if p_.get("in_corpus") and p_["pmid"] in by_pmid), None)
                if paper:
                    # the dataset is a corpus paper's own data: link from that paper, no separate node
                    if by_pmid[paper["pmid"]]["id"] != n["id"]:
                        edges.append({"from": by_pmid[paper["pmid"]]["id"], "to": n["id"], "kind": "pooled study", "evidence": ev})
                    continue
                sid = "smri:" + key
                if sd:
                    inv = sd["short_name"].replace("stanley_", "")
                    inv = inv[0].upper() + inv[1:-1] + inv[-1].upper() if inv[-1] in "ABC" and inv[:-1] in ("altar", "sklar") else inv.capitalize()
                    dp = sd.get("describing_papers") or []
                    s = stanley_nodes.setdefault(sid, {
                        "id": sid, "pmid": dp[0]["pmid"] if dp else "",
                        "label": f"Stanley {inv}", "year": sd.get("year") or 2006, "year_note": "" if sd.get("year") else
                        "no describing paper found; placed at the SMRIDB release (Higgs 2006, PMID 16594998)",
                        "title": f"SMRIDB study {sd['short_name']} ({sd['collection']}; {sd.get('investigator', '')}): {sd['region']}, {sd['platform']}, {sd['n_samples']} samples"
                                 + (f". Describing paper: {dp[0]['first_author']} {dp[0]['year']}, {dp[0]['title']}" if dp else ". No describing paper found."),
                        "journal": "", "relation": "source", "external": True, "collection": False, "dataset": True,
                        "cited_by": [], "lanes": ["stanley"], "lane": "stanley", "lane_label": "Stanley (SMRI) collections",
                        "parent": "", "tissue_n": f"{sd['region']}; {sd['n_samples']} samples (all diagnoses)", "n": None,
                        "platform": sd["platform"], "genes": "", "source": "Cited source (not in corpus)", "parent_by": "",
                        "gemma": [{"id": sd["gemma_id"], "short_name": sd["short_name"]}], "geo_ds": sd.get("geo", [])})
                else:
                    # SMRIDB studies beyond the 12 (e.g. Chu, Laeng, Aston in Bergon 2015)
                    s = stanley_nodes.setdefault(sid, {
                        "id": sid, "pmid": "", "label": f"Stanley study {key[5:]}" + (f" ({ds.get('name')})" if ds.get("name") else ""), "year": 2006,
                        "year_note": "not one of the 12 SMRIDB studies in Higgs 2006; year unknown, placed at 2006",
                        "title": f"SMRIDB study {ds.get('id', '')} ({ds.get('name', '')}): {ds.get('tissue', '')}, {ds.get('platform', '')}",
                        "journal": "", "relation": "source", "external": True, "collection": False, "dataset": True,
                        "cited_by": [], "lanes": ["stanley"], "lane": "stanley", "lane_label": "Stanley (SMRI) collections",
                        "parent": "", "tissue_n": ds.get("tissue", ""), "n": None, "platform": ds.get("platform", ""), "genes": "",
                        "source": "Cited source (not in corpus)", "parent_by": "", "gemma": []})
                s["cited_by"].append(n["label"])
                edges.append({"from": sid, "to": n["id"], "kind": "pooled study", "evidence": ev})
            elif ds["kind"] == "paper":
                m = re.search(r"(\d{4})", ds.get("cite", ""))
                if not m:
                    continue
                au, yr = ds["cite"].split()[0].lower(), m.group(1)
                cands = [c for c in nodes if not c.get("external") and c["label"].split()[0].lower() == au and str(c["year"]) == yr]
                cand, why = pick_cited(cands, [f"{ds.get('n', '')} {ds.get('name', '')}", entry.get("source", "")], ds["cite"])
                if cand:
                    edges.append({"from": cand["id"], "to": n["id"], "kind": "pooled study", "evidence": ev, **({"match": why} if why else {})})
                elif cands:
                    n.setdefault("cited_sources", []).append({"label": ds["cite"], "title": "", "pmid": "", "kind": "pooled study (ambiguous citation)", "evidence": f"{ev}; {why}"})
    for x in stanley_nodes.values():
        x["lane_label"] = x["lane_label"] or next((m["lane_label"] for m in nodes if m["lane"] == x["lane"]), "")
    nodes += list(stanley_nodes.values())
    # Reuse links stated in a paper but not visible from citations (pooled_datasets.json _curated_edges).
    def by_cite(c):
        au, yr = c.split()[0].lower(), re.search(r"(\d{4})", c).group(1)
        return next((x for x in nodes if not x.get("external") and x["label"].split()[0].lower() == au and str(x["year"]) == yr), None)
    def by_pm(pm):  # a PMID, where a cite string is ambiguous (two corpus nodes "Zhou 2025")
        return next((x for x in nodes if not x.get("external") and pm and x.get("pmid") == pm), None)
    for ce in json.loads(pooled_path.read_text()).get("_curated_edges", []) if pooled_path.exists() else []:
        a = by_pm(ce.get("from_pmid")) or by_cite(ce["from_cite"])
        b = by_pm(ce.get("to_pmid")) or by_cite(ce["to_cite"])
        if a and b:
            # the curated kind wins over a link read from citations for the same pair: duplicates keep the first link, so
            # a model's 'same subjects' link, dropped later as a kind, took Du 2023 -> Du 2025 with it (found 2026-09-29)
            edges[:] = [e for e in edges if not (e["from"] == a["id"] and e["to"] == b["id"])]
            edges.append({"from": a["id"], "to": b["id"], "kind": ce["kind"], "evidence": ce["evidence"]})

    for n in nodes:
        if n.get("dataset"):
            continue
        n["geo"] = geo.get("paper_gse", {}).get(n.get("pmid"), [])
        n["gemma"] = geo.get("paper_gemma", {}).get(n.get("pmid"), [])
    for e in external.values():
        e["lane"] = max(sorted(set(e["lanes"])), key=e["lanes"].count)  # sorted: deterministic ties
        e["lane_label"] = next((n["lane_label"] for n in nodes if n["lane"] == e["lane"]), "")
        if not e["year"]:
            e["year"] = min(by_label["year"] for by_label in nodes if by_label["label"] in e["cited_by"])
    REL2CAT = {"primary": "1", "single-gene follow-up": "2", "data reuse": "3", "reanalysis of public data": "3", "meta-analysis / combined": "3"}
    for n in nodes:
        if n.get("external") or n.get("lane") == "sc":
            continue
        # a study with no PMID (the round-4 bioRxiv preprint, 2026-09-30) has its facts under its DOI key; looked up by PMID
        # only, it lost its also_reuses flag and molecule
        c = CATS.get(n.get("pmid") or n["id"], {})
        n["category"] = c.get("category") or REL2CAT.get(n["relation"], "1")
        if str(n.get("pmid") or n["id"]) in CORPUS_EXCLUSIONS:  # corpus_exclusions.json: left the corpus by decision; facts unchanged
            n["category"] = "2"
            n["corpus_exclusion"] = CORPUS_EXCLUSIONS[str(n.get("pmid") or n["id"])]["date"]
        n["also_reuses"] = bool(c.get("also_reuses")) and n["category"] != "3"
        n["new_validation"] = bool(c.get("new_validation"))
        n["molecule"] = c.get("molecule", "")
        n["category_note"] = c.get("note", "")
        if n["category"] == "3" and n["lane"] not in ("meta",) and n["relation"] == "meta-analysis / combined":
            n["lane"], n["lane_label"] = META_LANE

    # Keep only cited sources classified as a real sample/cohort source or reused expression data
    # (classify_sources.py; hand overrides in sources_class_overrides.json win).
    cls_path, ovr_path = HERE / "sources_class.json", HERE / "sources_class_overrides.json"
    # Every candidate source, before filtering, is written out for classify_sources.py.
    ev_by = {}
    for e in edges:
        ev_by.setdefault(e["from"], []).append(f"{e['kind']}: {e.get('evidence', '')[:160]}")
    (HERE / "sources_all.json").write_text(json.dumps(
        # sorted: `edges` and `external` come out of set iteration, so their order (and, past three links, which three) varied
        # with the hash seed; two builds of the same inputs wrote different bytes (2026-10-08)
        [{"id": x["id"], "cite": x["label"], "title": x["title"], "links": sorted(ev_by.get(x["id"], []))[:3]}
         for x in sorted(external.values(), key=lambda x: str(x["id"]))],
        indent=1, ensure_ascii=False))
    if cls_path.exists():
        cls = {k: v["class"] for k, v in json.loads(cls_path.read_text())["items"].items()}
        cls.update(json.loads(ovr_path.read_text()) if ovr_path.exists() else {})
        # Unclassified sources are dropped until classify_sources.py has seen them.
        dropped = {k for k in external if cls.get("ext:" + k, "unclassified") not in ("sample_source", "expression_data")}
        for k in dropped:
            external.pop(k)
        edges = [e for e in edges if not (e["from"].startswith("ext:") and e["from"][4:] in dropped)]
    # Corpus papers that GEO links to the same series share data: earliest -> later ("same GEO series").
    by_gse: dict[str, list[dict]] = {}
    for pm_, gses in geo.get("paper_gse", {}).items():
        if pm_ in by_pmid:
            for g in gses:
                by_gse.setdefault(g, []).append(by_pmid[pm_])
    for g, ps in by_gse.items():
        ps = sorted(ps, key=lambda x: (x["year"], x["id"]))
        for later in ps[1:]:
            if later["year"] >= ps[0]["year"]:
                edges.append({"from": ps[0]["id"], "to": later["id"], "kind": "same dataset",
                              "evidence": f"GEO links both papers to {g}"})
    # A cited paper that is also the describing paper of an SMRIDB dataset node is that node (Ryan 2006, PMID 16894394, is
    # smri:stanley_bahn): its links move there, so one paper is one node.
    smri_by_pmid = {x["pmid"]: x for x in stanley_nodes.values() if x.get("pmid")}
    for pm_ in [k for k in external if k in smri_by_pmid]:
        tgt = smri_by_pmid[pm_]
        for e in edges:
            if e["from"] == "ext:" + pm_:
                e["from"] = tgt["id"]
        tgt["cited_by"] += [c for c in external.pop(pm_)["cited_by"] if c not in tgt["cited_by"]]
    # Links mean data re-use only (Paul, 2026-09-23: "I don't want to track mere citations"): drop a study
    # citing a brain-collection description or an unspecified "other" relation.
    edges = [e for e in edges if e["kind"] not in ("brain collection description", "other", "")]
    # Drop duplicate edges (same pair, first reason kept); only keep source nodes that anchor an edge.
    # edges whose evidence is a guess ("likely source", "typical source") are not data re-use (review 2026-09-24)
    edges = [e for e in edges if not re.search(r"\blikely\b|\btypical|presumab|\bprobabl|\bassum", e.get("evidence", ""), re.I)]
    seen, uniq = set(), []
    for e in edges:
        if (e["from"], e["to"]) not in seen:
            seen.add((e["from"], e["to"])); uniq.append(e)
    edges = uniq
    # Merikangas et al. 2022: the review this project extends. A landmark node, linked to its 26 studies.
    mer = {"id": "merikangas2022", "pmid": "35091668", "label": "Merikangas 2022", "year": 2022,
           "title": "What genes are differentially expressed in individuals with schizophrenia? A systematic review",
           "journal": "Mol Psychiatry", "relation": "systematic review", "category": "review", "landmark": True,
           "lane": "meta", "lane_label": META_LANE[1], "parent": "", "tissue_n": "26 studies, search window 2000-01 to 2020-06",
           "n": None, "platform": "", "genes": "6,771 genes reported DE in any study; 160 in 3+ studies; 22 in 4+ (Supp Table 1)",
           "source": "Merikangas 2022", "parent_by": "", "geo": [], "gemma": [], "window": [2000.0, 2020.5]}
    nodes.append(mer)
    for x in nodes:
        if x.get("source") == "Merikangas 2022" and x is not mer:
            edges.append({"from": x["id"], "to": mer["id"], "kind": "reviewed", "evidence": "one of the 26 studies in Merikangas 2022"})
    anchored = {e["from"] for e in edges}
    nodes += [x for x in external.values() if x["id"] in anchored]
    edges = drop_same_subjects(nodes, edges)
    nodes, edges, prune_stats = prune_sources(nodes, edges)
    enrich_collections(nodes, edges)

    # Step 2 gene-list status per study (analyses, read only).
    gl = genelist_status()
    for n in nodes:
        if n.get("external") or n.get("landmark"):
            continue
        n["genelist"] = gl.get(n.get("pmid") or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")) or {"status": "not attempted", "n_lists": 0}  # DOI key: no PMID (Cui 2024)

    for n in nodes:
        for nm in n.pop("stanley_other", []):
            n.setdefault("cited_sources", []).append(
                {"label": nm, "title": "Stanley Neuropathology Consortium measure, not expression data", "kind": "other data on the same samples", "pmid": ""})

    # Cohort from authorship: a brain bank's director on the author list is strong evidence the
    # samples are that bank's (Paul, 2026-09-23). Applied only to studies whose parent text named
    # no cohort, and only when exactly one bank's heads are authors.
    authors = author_lists({n["pmid"] for n in nodes if n.get("pmid") and not n.get("external")},
                           {k: v.get("doi") for k, v in json.loads(CACHE.read_text()).items()} if CACHE.exists() else {})
    lane_label = {k: lab for k, lab, _ in LANES}
    for n in nodes:
        if n.get("external") or n.get("landmark") or n["lane"] not in (OTHER_BRAIN[0], OTHER_PERIPH[0], "meta"):
            continue
        if n["lane"] == OTHER_PERIPH[0] and "pritzker" in BANK_HEADS and any(
                h in authors.get(n.get("pmid"), set()) for h in BANK_HEADS["pritzker"]):
            # Vawter's blood and cell-line family studies (Vawter 2004, Martin 2009) are not UC Irvine brain bank samples
            n["bank_authors"] = {}
            continue
        au = authors.get(n.get("pmid"), set())
        banks = {b: sorted(h for h in heads if h in au) for b, heads in BANK_HEADS.items()}
        banks = {b: h for b, h in banks.items() if h}
        n["bank_authors"] = banks
        if len(banks) == 1 and n["lane"] != "meta":
            b, h = next(iter(banks.items()))
            n["lane"], n["lane_label"], n["lane_by"] = b, lane_label[b], "author: " + ", ".join(h)

    # Meta-analysis lane only for pooling across cohorts (Paul, 2026-09-23: "re-uses by a consortium
    # doesn't count"). A cohort is a named lane, or one input study/dataset from an unnamed lane.
    by_id = {n["id"]: n for n in nodes}
    named = {k for k, _, _ in LANES}
    for n in nodes:
        if n["lane"] != "meta" or n.get("landmark") or n.get("external"):
            continue
        units = {k for k, _, pat in LANES if re.search(pat, n["parent"], re.I)}
        for e in edges:
            if e["to"] == n["id"] and e["from"] in by_id and e["kind"] != "reviewed":
                src = by_id[e["from"]]
                if src.get("external") and src["lane"] == "meta":
                    continue  # a database/resource paper (e.g. SNCID) is not a cohort
                units.add(src["lane"] if src["lane"] in named else src["id"])
        n["pools"] = sorted(u for u in units if u in named)
        n["n_cohorts"] = len(units)
        # Move out only when exactly one cohort is identified and the parent text does not describe several
        # datasets; zero identified inputs means unknown, not single-cohort (review 2026-09-24: 41 moved wrongly).
        several = len(set(re.findall(r"GSE\d+", n["parent"]))) >= 2 or re.search(
            r"\b([2-9]|\d\d+|two|three|four|five|six|seven|eight|nine|ten|multiple|several)\s+([\w-]+\s+){0,3}"
            r"(datasets|data sets|cohorts|studies|brain banks|collections)\b", n["parent"] + " " + n["tissue_n"], re.I)
        if len(units) == 1 and not several:
            home = next((k for k in units if k in named), None)
            if not home and len(n.get("bank_authors", {})) == 1:
                home = next(iter(n["bank_authors"]))
                n["lane_by"] = "author: " + ", ".join(n["bank_authors"][home])
            if home:
                n["lane"], n["lane_label"] = home, lane_label[home]
            else:
                n["lane"], n["lane_label"] = lane_for(n["parent"], n["tissue_n"], "primary")
            n["meta_note"] = f"combined analysis within one cohort ({n['lane_label']}); not drawn as a cross-cohort meta-analysis"
    # Sources follow the lane of the studies citing them, which may have moved.
    for e in external.values():
        cit = [by_id[x["to"]]["lane"] for x in edges if x["from"] == e["id"] and x["to"] in by_id]
        if cit:
            e["lane"] = max(sorted(set(cit)), key=cit.count)  # sorted: ties break the same way every run
            e["lane_label"] = lane_label.get(e["lane"]) or dict([OTHER_BRAIN, OTHER_PERIPH, META_LANE]).get(e["lane"], "")

    # Single-cell / single-nucleus studies (single_cell.py). Analysed as cell-type-resolved lists
    # (Paul: bring them into the analysed corpus). Each sits in the lane of its brain bank or collection
    # (sc_banks.json; Paul: "the single-cell studies have to go in the lane for their sources, or near them, and link up
    # to the consortium nodes"), flagged single_cell so the page can mark and filter them.
    sc_banks = json.loads((HERE / "sc_banks.json").read_text()) if (HERE / "sc_banks.json").exists() else {"studies": {}, "edges": []}
    lane_lab_all = {k: lab for k, lab, _ in LANES} | {OTHER_BRAIN[0]: OTHER_BRAIN[1], OTHER_PERIPH[0]: OTHER_PERIPH[1], META_LANE[0]: META_LANE[1]}
    sc_path = HERE / "single_cell.json"
    if sc_path.exists():
        sc = json.loads(sc_path.read_text())
        have = {n.get("pmid") for n in nodes if not n.get("external")}
        seen_titles = set()
        for st in sc["studies"] + sc.get("studies_by_hand", []) + sc.get("preprints_by_hand", []):
            key = st.get("pmid") or st.get("doi")
            t = st["title"].lower()[:60]
            if key in have or t in seen_titles:
                continue
            # a single-cell study already drawn as a cited source (e.g. Batiuk 2022, pooled by meta-analyses)
            # becomes the corpus node; its edges move to it
            ext_id = "ext:" + key
            if any(x["id"] == ext_id for x in nodes):
                for e in edges:
                    if e["from"] == ext_id:
                        e["from"] = "sc:" + key
                    if e["to"] == ext_id:
                        e["to"] = "sc:" + key
                nodes = [x for x in nodes if x["id"] != ext_id]
            seen_titles.add(t)
            n_tot = (st.get("cases") or 0) + (st.get("controls") or 0)
            bk = sc_banks["studies"].get(key, {})
            sc_lane = bk.get("lane", "sc")
            nodes.append({"id": "sc:" + key, "pmid": st.get("pmid", ""), "label": st.get("label") or f"{html.unescape(st['first_author'])} {st['year']}",
                          "year": int(st["year"]), "title": st["title"], "journal": "", "relation": "single-cell / single-nucleus",
                          "single_cell": True, "bank_names": bk.get("bank_names", []),
                          **({"banks": bk["banks"]} if len(bk.get("banks") or []) > 1 else {}),
                          **({"lane_by": "sc_banks.json: " + bk["evidence"]} if bk.get("evidence") else {}),
                          "lane": sc_lane, "lane_label": lane_lab_all.get(sc_lane, "Single-cell / single-nucleus"), "parent": "",
                          "tissue_n": st.get("tissue", ""), "n": n_tot or None, "cases": st.get("cases"), "controls": st.get("controls"),
                          "platform": "single-cell / single-nucleus RNA-seq", "genes": "", "source": "Single-cell / single-nucleus",
                          "parent_by": "", "geo": geo.get("paper_gse", {}).get(st.get("pmid"), []),
                          "gemma": geo.get("paper_gemma", {}).get(st.get("pmid"), [])})

    for n in nodes:
        if n.get("single_cell"):
            # preprints have no PMID; Step 2 keys them as doi_<doi with / as _> (manifests/batch9_sc.json)
            key = n.get("pmid") or "doi_" + n["id"][3:].replace("/", "_")
            n["genelist"] = gl.get(key) or {"status": "not attempted", "n_lists": 0}

    # Hand-curated data re-use and subject links (curated_edges.json, e.g. Hoffman 2021; sc_banks.json edges for the
    # single-cell studies). An end is a corpus PMID / DOI key, or a consortium node id (coll:...), created here if no
    # meta-analysis pooled it. Kinds follow the tracker Vocabulary (subjects, samples, regions).
    by_key = {}
    for x in nodes:
        if x.get("external") or x.get("landmark"):
            continue
        by_key.setdefault(x.get("pmid") or x["id"].split(":", 1)[-1], x)
        by_key.setdefault(x["id"], x)
    coll_lane = {"coll:cmc": "cmc", "coll:psychencode": "cmc", "coll:nimh_hbcc": "cmc", "coll:libd_brainseq": "libd",
                 "coll:stanley_array_collection": "stanley", "coll:stanley_neuropathology_consortium": "stanley",
                 "coll:stanley_collections": "stanley", "coll:pritzker": "pritzker"}
    ce_path = HERE / "curated_edges.json"
    ce_all = json.loads(ce_path.read_text()) if ce_path.exists() else {}
    curated = ce_all.get("edges", []) + sc_banks.get("edges", [])
    n_curated, missing = 0, []
    for ce in curated:
        a_key, b_key = ce.get("from_id") or ce.get("from_pmid") or ce.get("from"), ce.get("to_pmid") or ce.get("to_id") or ce.get("to")
        b = by_key.get(b_key)
        if b is None and str(b_key).startswith("coll:"):  # one collection's data included in another (PsychENCODE)
            b = next((x for x in nodes if x["id"] == b_key), None)
        if a_key.startswith(("ext:", "smri:")):  # a drawn source node (a GEO series' paper, an SMRIDB dataset)
            a = next((x for x in nodes if x["id"] == a_key), None)
        elif a_key.startswith("coll:"):
            a = next((x for x in nodes if x["id"] == a_key), None)
            if a is None and a_key in coll_lane:
                a = {"id": a_key, "pmid": "", "label": a_key[5:], "year": 0, "title": f"{a_key[5:]} (data resource)", "journal": "",
                     "relation": "source", "external": True, "collection": True, "cited_by": [], "lanes": [coll_lane[a_key]],
                     "lane": coll_lane[a_key], "lane_label": lane_lab_all.get(coll_lane[a_key], ""), "parent": "", "tissue_n": "",
                     "n": None, "platform": "", "genes": "", "source": "Cited source (not in corpus)", "parent_by": "", "n_citing": 0}
                nodes.append(a)
        else:
            a = by_key.get(a_key)
        if not a or not b:
            missing.append(f"{a_key} -> {b_key}")
            continue
        if any(e["from"] == a["id"] and e["to"] == b["id"] for e in edges):
            edges[:] = [e for e in edges if not (e["from"] == a["id"] and e["to"] == b["id"])]  # the curated kind wins
        edges.append({"from": a["id"], "to": b["id"], "kind": ce["kind"], "evidence": ce["evidence"], "curated": True})
        n_curated += 1
        # a curated edge settles an ambiguous citation that listed this source as a candidate (Hess 2016: Katsel 2005,
        # de Jong 2012; 2026-09-30)
        if a.get("pmid") and b.get("cited_sources"):
            b["cited_sources"] = [c for c in b["cited_sources"] if not ("ambiguous citation" in c.get("kind", "")
                                                                         and f"(PMID {a['pmid']})" in c.get("evidence", ""))]
        if a.get("external") and not b.get("external"):
            if b["label"] not in a.setdefault("cited_by", []):
                a["cited_by"].append(b["label"])
            a["n_citing"] = len(a["cited_by"])
            if not a.get("year") or (b.get("year") and b["year"] < a["year"]):
                a["year"] = b.get("year") or a.get("year")
    print(f"curated edges: {n_curated} added" + (f"; {len(missing)} with an end not in the graph: {missing}" if missing else ""))
    # Data a study re-used that no drawn node holds (curated_edges.json "other_data", keyed by study PMID / node id): listed
    # under the study's "Other data used", like the folded sources below (2026-09-29)
    od_missing = []
    for key, entries in ce_all.get("other_data", {}).items():
        if key.startswith("_"):
            continue
        n = by_key.get(key)
        if n is None:
            od_missing.append(key)
            continue
        for o in entries:
            n.setdefault("other_data", []).append({k: o[k] for k in ("label", "pmid", "gse", "why", "kind") if k in o})
    if od_missing:
        print(f"curated other_data: studies not in the graph: {od_missing}")
    lane_order = [l[0] for l in LANES] + [OTHER_BRAIN[0], OTHER_PERIPH[0], META_LANE[0], "sc"]
    lane_labels = {l[0]: l[1] for l in LANES} | {OTHER_BRAIN[0]: OTHER_BRAIN[1], OTHER_PERIPH[0]: OTHER_PERIPH[1], META_LANE[0]: META_LANE[1],
                                                 "sc": "Single-cell / single-nucleus"}
    # Step 2 records a tissue per study from its lists (analyses/primary_lists.json); it fills in where the
    # eligibility text named none. Without it, a study with no tissue text fell into the blood/other lane
    # (Maycox 2009, a brain study, 2026-09-24).
    cites = json.loads((HERE / "citations.json").read_text()) if (HERE / "citations.json").exists() else {}  # citation_counts.py
    th = json.loads((HERE / "themes.json").read_text()) if (HERE / "themes.json").exists() else {}  # themes.py
    for n in nodes:
        t = th.get("studies", {}).get(n["id"])
        if t and not n.get("external"):
            n["themes"] = [x["theme"] for x in t["themes"]]
            n["theme_evidence"] = {x["theme"]: x["quote"] for x in t["themes"]}
            if t.get("title_only"):
                n["themes_title_only"] = True
    nfill = json.loads((HERE / "n_from_lists.json").read_text())["studies"] if (HERE / "n_from_lists.json").exists() else {}
    for n in nodes:  # N stated in a Step 2 list contrast, where nothing else gave one (fill_n.py)
        f = nfill.get(n.get("pmid") or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_"))
        if f and not n.get("external") and not n.get("n"):
            n["n"], n["cases"], n["controls"] = f["n"], f["cases"], f["controls"]
            n["n_by"] = f"Step 2 list contrast ({f['list']}): \"{f['quote']}\""
    for n in nodes:  # Europe PMC times cited, studies and outside sources alike
        k = str(n.get("pmid") or "") if n.get("external") else study_key(n)  # DOI key for a study with no PMID (citation_counts.py)
        if k in cites.get("counts", {}):
            n["times_cited"] = cites["counts"][k]  # not "cited_by": that is the list of citing corpus studies
    # outside sources: N from a stated sample count in their description (the SMRIDB datasets: "44 samples")
    for n in nodes:
        if n.get("external") and not n.get("n"):
            m = re.search(r"(\d[\d,]*) samples", f"{n.get('title') or ''} {n.get('tissue_n') or ''}")
            if m:
                n["n"], n["n_by"] = int(m.group(1).replace(",", "")), "sample count in the dataset description"
    mo_path = DEGENES / "analysis_modality.json"
    modality = json.loads(mo_path.read_text())["studies"] if mo_path.exists() else {}
    rk_path = DEGENES / "result_kinds.json"
    result_kinds = json.loads(rk_path.read_text())["studies"] if rk_path.exists() else {}
    import sys
    sys.path.insert(0, str(DEGENES))
    from analysis_modality import KIND_TO_MODALITY
    pl_path = DEGENES / "primary_lists.json"
    pl_tissue = {k: v.get("tissue") for k, v in json.loads(pl_path.read_text())["studies"].items()} if pl_path.exists() else {}
    for n in nodes:
        if not n.get("external"):
            n["tissue"] = tissue_class(n)
            if n["tissue"] == "unknown" and pl_tissue.get(n.get("pmid")) in ("brain", "blood", "lcl", "ipsc", "other"):
                n["tissue"], n["tissue_by"] = pl_tissue[n["pmid"]], "gene lists"
            n["tissue_short"] = tissue_short(n)
            # analysis modality (analyses/analysis_modality.py): a study whose result is not gene-level DE shows
            # "not DE" rather than "no list published" (Paul 2026-09-24, on Kiss 2026: cell composition)
            mo = modality.get(study_key(n))  # DOI key for a study with no PMID (Cui 2024, preprints; 2026-10-01)
            if mo:
                n["modality"], n["de_step"] = mo["modality"], mo["de_step"]
                # DE only through Gemma's reanalysis while the paper's own analysis (result_kinds) is something else:
                # the detail panel shows both (2026-09-25, de Jong 2012: WGCNA modules; Gemma GSE38485 lists)
                own = KIND_TO_MODALITY.get((result_kinds.get(study_key(n)) or {}).get("kind"))
                sts = {x.get("status") for x in (n.get("genelist") or {}).get("lists", [])}
                if mo["modality"] == "DE" and own and own != "DE" and "gemma" in sts and not sts & {"extracted", "reproduced"}:
                    n["paper_modality"] = own
                g = n.get("genelist")
                if g and g.get("status") == "no list published" and mo["modality"] != "DE":
                    g["status"] = "not DE"
    # Card headline (headline()), technology tier (analyses/classify_technology.py: technology.json where the study has a
    # primary list, else the same keyword rule on the platform text) and citations per year since publication.
    from classify_technology import array_generation, tier, technology_of_text, TIERS
    uni_lists = json.loads((DEGENES / "universes.json").read_text())["lists"] if (DEGENES / "universes.json").exists() else {}
    tech_path = DEGENES / "technology.json"
    techs = json.loads(tech_path.read_text())["studies"] if tech_path.exists() else {}
    pl_all = json.loads(pl_path.read_text())["studies"] if pl_path.exists() else {}
    titles = platform_titles()
    fetched_year = int(str(cites.get("fetched") or datetime.now().year)[:4])
    # Field corrections for every study node, including the single-cell nodes built from single_cell.json (the corpus rows
    # got theirs when built, above; applying again is a no-op). Keyed by node id, or PMID for "sc:" ids.
    def fov(n):
        return FIELD_OV.get(n["id"]) or FIELD_OV.get(n.get("pmid") or "") or FIELD_OV.get(n["id"].removeprefix("sc:"))
    for n in nodes:
        f_ = None if n.get("external") or n.get("landmark") else fov(n)
        if not f_:
            continue
        if f_.get("platform"):
            n["platform"] = f_["platform"]
        if f_.get("tissue_n") and f_["tissue_n"] != n.get("tissue_n"):
            n["tissue_n"] = f_["tissue_n"]
            if not n.get("tissue_by"):  # the tissue class follows the corrected text (set above from the old one)
                n["tissue"] = tissue_class(n)
            n["tissue_short"] = tissue_short(n)
        if f_.get("cases") is not None or f_.get("controls") is not None:
            n["cases"], n["controls"] = f_.get("cases"), f_.get("controls")
            # n only when both sides are known (see n_for)
            n["n"] = n["cases"] + n["controls"] if n["cases"] is not None and n["controls"] is not None else None
            n.pop("n_by", None)  # an N filled from a list contrast (n_from_lists.json) above is replaced by the override
        # cases_sz: the schizophrenia count of a pooled psychosis case group, whose cases field counts every case in the
        # contrast (N audit convention 2, decisions_paul_2026-09-28.json); distinct_subjects.py uses it for SZ people
        if f_.get("cases_sz") is not None:
            n["cases_sz"] = f_["cases_sz"]
        # n_note: a documented conflict between two counts in the paper; the page marks N with an asterisk
        if f_.get("n_note"):
            n["n_note"] = f_["n_note"]
        n["checked_by"] = f_.get("why", "")
    # The GEO datasets each re-use study pooled, as GEO describes them (pooled_geo.py; Paul, 2026-09-27: "go by what the
    # data are" and, on meta-analyses that pooled a wrong dataset, "we should be consistent"): shown on the study, and the
    # ones with no schizophrenia in the GEO record or its samples join "Other data used".
    pg_path = HERE / "pooled_geo.json"
    PG = json.loads(pg_path.read_text())["studies"] if pg_path.exists() else {}
    for n in nodes:
        pg = PG.get(n.get("pmid") or n["id"]) if not n.get("external") else None
        if not pg:
            continue
        n["geo_datasets"] = [{"gse": x["gse"], "title": x["title"], "tier": x["tier"], "sz_data": x["sz_data"],
                              "platforms": "; ".join(f"{p_['gpl']} {p_['title'] or ''}".strip() for p_ in x["platforms"])} for x in pg["datasets"]]
        have = {o["label"] for o in n.get("other_data", [])} | {o.get("gse") for o in n.get("other_data", [])}
        for x in pg["datasets"]:
            if x["gse"] in pg["not_sz_data"] and x["gse"] not in have:
                n.setdefault("other_data", []).append({"label": x["gse"], "gse": x["gse"], "pmid": "", "kind": "pooled (GEO)",
                    "why": (x["title"] or "")[:90] + ("" if x["human"] else f" ({x['taxon']})") + "; no schizophrenia in its GEO record or samples"})
    for n in nodes:  # one entry per dataset (a source linked twice, e.g. pooled and public dataset, was listed twice)
        if n.get("other_data"):
            seen_o, kept_o = set(), []
            for o in n["other_data"]:
                if o["label"] not in seen_o:
                    seen_o.add(o["label"])
                    kept_o.append(o)
            n["other_data"] = kept_o
    tier_conflicts = []
    for n in nodes:
        if n.get("times_cited") is not None and n.get("year"):
            n["cites_per_year"] = round(n["times_cited"] / max(1, fetched_year - int(n["year"])), 1)
        if n.get("external") or n.get("landmark"):
            continue
        key = study_key(n)
        prim = [p for p in (pl_all.get(key) or {}).get("primary", [])]
        n.update(headline(n, prim, uni_lists, titles))
        te = techs.get(key)
        if isSC_node(n):
            t, gen, by = "single_cell", None, "single-cell / single-nucleus study"
        elif te:
            t, gen, by = te["technology"], te.get("array_generation"), f"technology.json: {te['evidence']}" + (f"; {te.get('array_generation_evidence')}" if te.get("array_generation_evidence") else "")
        else:
            # The eligibility extraction read an abstract or excerpt and often names no platform ("Not stated in provided
            # excerpt" for Katsel 2008, which used HG-U133A/B); the Step 2 gene-list record, read from the paper, names it
            # in its contrast. Use that when the platform text names no technology or no array (Paul, 2026-09-27).
            ptxt = n.get("platform") or ""
            ltxt = " | ".join(l.get("contrast", "") for l in (n.get("genelist") or {}).get("lists", []))
            # first looks only: a re-use study's list text names its inputs, and those studies' tiers are inherited below
            if ltxt and str(n.get("category")) in ("1", "1m") and (technology_of_text(ptxt) == "unknown" or (technology_of_text(ptxt) == "microarray" and array_generation(ptxt)[0] == "not stated")) \
                    and technology_of_text(ltxt) != "unknown":
                n["platform_from_list"], ptxt = ltxt, ltxt
                if not re.search(r"\w", n.get("platform") or "") or re.search(r"not (stated|specified|given)|none stated|unspecified|excerpt|specific pla", n.get("platform") or "", re.I):
                    parts = [q.strip() for q in re.split(r"[,;(|]", ltxt) if technology_of_text(q) != "unknown"]
                    if parts:
                        n["platform"] = parts[0] + " (from our gene-list record)"
            t = "proteomics" if n.get("modality") == "proteomics" or n.get("category") == "p" else \
                "single_cell" if re.search(r"\b(sn|sc)RNA-?seq|single[- ](cell|nucle)", ptxt, re.I) else technology_of_text(ptxt)
            gen, gev = array_generation(ptxt) if t == "microarray" else (None, "")
            if n.get("platform_from_list"):
                gev = (gev + "; " if gev else "") + "platform from the gene-list record (Step 2 contrast), the extracted platform text named none"
            by = f"platform text ({'proteomics study' if t == 'proteomics' else 'no primary list in technology.json'})" + (f"; {gev}" if gev else "")
        if t == "unknown" and re.search(r"\bpcr\b|qpcr|in situ|northern|nanostring|differential display|immunohisto|western", n.get("platform") or "", re.I):
            t, by = "other", "platform text: a targeted method (qPCR, in situ, NanoString, ...)"
        n["tech_tier"], n["tech_tier_by"] = tier(t, gen), by
        pg = PG.get(n.get("pmid") or n["id"])
        if not te and str(n.get("category")) == "3" and pg and pg.get("tier_by_data"):
            n["tech_tier"] = pg["tier_by_data"]
            n["tech_tier_by"] = "GEO platforms of the pooled datasets (pooled_geo.json): " + ", ".join(pg["dataset_tiers"])
            continue
        tov = (fov(n) or {}).get("tech_tier")
        if tov and tov != n["tech_tier"]:
            if te and n["tech_tier"] not in ("array_unstated", "unknown"):
                # technology.json (used by the analyses) named a specific tier from the list's own probe IDs or universe;
                # a disagreement with the paper check is listed for review, not silently diverged (Kuzman 2009's list
                # has U133 Plus 2.0-only probe sets though the paper says U133A: build_universes.py)
                tier_conflicts.append({"id": n["id"], "label": n["label"], "technology_json": n["tech_tier"], "paper_check": tov,
                                       "why": fov(n).get("why", "")})
            else:
                n["tech_tier"], n["tech_tier_by"] = tov, "paper check: " + fov(n).get("why", "")[:200]
    (HERE / "tech_tier_conflicts.json").write_text(json.dumps(
        {"_about": "Studies whose technology tier from analyses/technology.json (used by the analyses) differs from the tier the "
                   "paper check (node_field_overrides.json) gives. The page follows technology.json, which is set by the data: the "
                   "platform of the study's own list or of the datasets it pooled (Paul, 2026-09-27: \"go by what the data are\").",
         "conflicts": tier_conflicts}, indent=1, ensure_ascii=False))
    if tier_conflicts:
        print(f"technology tier conflicts (technology.json vs paper check): {len(tier_conflicts)}; see tech_tier_conflicts.json")
    # A re-use study (meta-analysis, reanalysis) with no technology of its own takes that of the corpus studies whose data
    # it re-used (links pooled study / public dataset / same dataset / same samples); several technologies -> mixed.
    by_id_t = {x["id"]: x for x in nodes}
    for n in nodes:
        if n.get("external") or n.get("landmark") or n.get("tech_tier") != "unknown":
            continue
        ins = {by_id_t[e["from"]].get("tech_tier") for e in edges if e["to"] == n["id"] and e["from"] in by_id_t
               and e["kind"] in ("pooled study", "public dataset", "same dataset", "same samples (reanalysis)")} - {None, "unknown"}
        if ins:
            n["tech_tier"] = next(iter(ins)) if len(ins) == 1 else "mixed"
            n["tech_tier_by"] = "technology of the corpus studies whose data it re-used: " + ", ".join(sorted(ins))
    tier_labels = dict(TIERS)
    # iPSC studies belong in the ipsc lane whatever bank their donors came from; non-iPSC studies do not
    # (review 2026-09-24: 16 iPSC studies outside it, 2 non-iPSC inside it).
    lab = {k: l for k, l, _ in LANES} | {OTHER_BRAIN[0]: OTHER_BRAIN[1], OTHER_PERIPH[0]: OTHER_PERIPH[1]}
    for n in nodes:
        if n.get("external") or n.get("landmark") or n["lane"] in ("meta", "sc") or n.get("single_cell"):
            continue  # single-cell studies: lane from sc_banks.json
        if n.get("tissue") == "ipsc" and n["lane"] != "ipsc":
            n["lane"], n["lane_label"], n["lane_by"] = "ipsc", lab["ipsc"], "tissue: iPSC-derived"
        elif n.get("tissue") in ("blood", "lcl", "other") and n["lane"] not in (OTHER_PERIPH[0], "ipsc"):
            # a peripheral sample is not from a brain bank, whatever cohort the text names (Paul, 2026-09-27, on Gardiner
            # 2012: "why is it in the brain bank lane not blood?"; the Australian Schizophrenia Research Bank collects
            # blood from living participants)
            n["lane"], n["lane_label"], n["lane_by"] = OTHER_PERIPH[0], OTHER_PERIPH[1], f"tissue: {n['tissue']} (not a brain-bank sample)"
        elif n["lane"] == OTHER_PERIPH[0] and n.get("tissue") == "brain":
            n["lane"], n["lane_label"], n["lane_by"] = OTHER_BRAIN[0], OTHER_BRAIN[1], "tissue: brain"
        elif n["lane"] == "ipsc" and n.get("tissue") not in ("ipsc", "mixed"):
            k = next((k for k, _, pat in LANES if k != "ipsc" and re.search(pat, n["parent"], re.I)), None)
            k = k or (OTHER_BRAIN[0] if n.get("tissue") == "brain" else OTHER_PERIPH[0])
            n["lane"], n["lane_label"], n["lane_by"] = k, lab[k], "tissue: not iPSC"
    ov_path = HERE / "lane_overrides.json"
    for pm, ov in (json.loads(ov_path.read_text()).items() if ov_path.exists() else []):
        for n in nodes:
            if not pm.startswith("_") and n.get("pmid") == pm and not n.get("external"):
                n["lane"], n["lane_label"], n["lane_by"] = ov["lane"], lab.get(ov["lane"], META_LANE[1] if ov["lane"] == "meta" else ov["lane"]), "override: " + ov["why"]
                if ov.get("relation"):
                    n["relation"] = ov["relation"]
                if ov.get("also_banks"):
                    n["banks"] = [ov["lane"]] + ov["also_banks"]  # study spans several banks
                if ov["lane"] == "meta":
                    n.pop("meta_note", None)
    # Meta-analysis lane in the middle of the populated lanes (Paul, 2026-09-23).
    populated = {n["lane"] for n in nodes if not n.get("external")}
    # Bottom of the figure (Paul, 2026-09-24): other brain, then iPSC below it, single-cell above blood/other.
    tail = [k for k in (OTHER_BRAIN[0], "ipsc", "sc", OTHER_PERIPH[0]) if k in populated]
    base = [k for k in lane_order if k not in ["meta"] + tail and k in populated] + tail
    pos = len(base) // 2
    lane_order = base[:pos] + ["meta"] + base[pos:]
    lanes = [{"id": k, "label": lane_labels[k], "count": sum(n["lane"] == k and not n.get("external") for n in nodes)} for k in lane_order]
    lanes = [l for l in lanes if l["count"]]

    # Formal bank names with acronyms for lanes and consortium nodes (bank_names.json, each name sourced; Paul, 2026-09-25:
    # "put HBTRC, we need those kinds of things there"). Lane keys are unchanged, so the analyses are unaffected.
    bn_path = HERE / "bank_names.json"
    bn_raw = json.loads(bn_path.read_text()) if bn_path.exists() else {}
    bn = {**bn_raw.get("lanes", {}), **bn_raw.get("collections", {})}
    for l in lanes:
        if bn.get(l["id"], {}).get("label"):
            l["label"], l["short"], l["detail"] = bn[l["id"]]["label"], bn[l["id"]].get("short", ""), bn[l["id"]].get("detail", "")
    for n in nodes:
        if n.get("lane") in bn and bn[n["lane"]].get("label"):
            n["lane_label"] = bn[n["lane"]]["label"]
        if n["id"] in bn and bn[n["id"]].get("label"):
            n["label"], n["bank_detail"] = bn[n["id"]]["label"], bn[n["id"]].get("detail", "")

    # Brain region class for display (brain_regions); the mapping per study is written out so it can be checked.
    reg_out = {}
    for n in nodes:
        if n.get("external") or n.get("landmark") or n.get("tissue") != "brain":
            continue
        cls, disp, found, src = brain_regions(n)
        n["brain_class"], n["brain_regions"] = cls, disp
        reg_out[n["id"]] = {"label": n["label"], "class": cls, "regions": found, "display": disp, "source": src}
    from collections import Counter as _C
    (HERE / "tissue_regions.json").write_text(json.dumps({
        "_about": "Brain studies split by region for the graph's tissue colouring (build_lineage.py::brain_regions, 2026-09-25). "
                  "Classes: neocortex; other region (hippocampus / medial temporal incl. entorhinal, striatum, thalamus, cerebellum, "
                  "midbrain, white matter, other); neocortex and other; unclassified (no named region in the eligibility text or list "
                  "contrasts). 'source' is the text the regions were read from. Display only; analyses read analyses/primary_lists.json.",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "counts": dict(_C(v["class"] for v in reg_out.values())), "studies": reg_out}, indent=1, ensure_ascii=False))
    print("brain region classes:", dict(_C(v["class"] for v in reg_out.values())))

    # gene symbol search index over the extracted DE lists (Step 2); nodes get "gene_index" counts
    gindex, gindex_meta, gstats = gene_index(nodes)
    (HERE / "lineage_gene_stats.json").write_text(json.dumps(gstats, ensure_ascii=False, separators=(",", ":")))
    # Technology shown as "targeted panel" when every gene list of the study tests fewer than 500 pre-selected genes
    # (display only; analyses/technology.json, used by the analyses, keeps the assay's technology)
    scopes_by_node = {}
    for r in gindex.get("lists", []):
        scopes_by_node.setdefault(r[1], set()).add(r[9] or "u")
    for n in nodes:
        if scopes_by_node.get(n["id"]) == {"t"} and str(n.get("category")) in ("1", "1m"):  # primary studies only: a re-use study's data are not a panel
            n["tech_tier"], n["tech_tier_by"] = "panel", "every gene list tests fewer than 500 pre-selected genes (targeted panel)"
    print(f"gene index: {gindex_meta.get('n_lists', 0)} DE lists, {gindex_meta.get('n_symbols', 0)} symbols, "
          f"{gindex_meta.get('n_entries', 0)} symbol-list entries; skipped non-DE list files {gindex_meta.get('n_list_files_skipped_not_DE')}")

    # source_folds.json (hand-curated 2026-09-27): an alias source is the data of a node already in the graph, so its
    # links move there; a folded source is not schizophrenia case-control expression data. Folded nodes and their links
    # stay in lineage.json (h3.py counts a meta-analysis's inputs from them); each gets "folded" and each study that used
    # one lists it under "other_data", and the page draws neither the node nor its links.
    sf_path = HERE / "source_folds.json"
    sf = json.loads(sf_path.read_text()) if sf_path.exists() else {"alias": {}, "fold": {}}
    ids_now = {x["id"]: x for x in nodes}
    for src, a in sf["alias"].items():
        tgt = ids_now.get(a["to"])
        if src not in ids_now or not tgt:
            print(f"source_folds alias {src} -> {a['to']}: an end is not in the graph")
            continue
        for e in edges:
            for end in ("from", "to"):
                if e[end] == src:
                    e[end] = tgt["id"]
                    e["evidence"] = (e.get("evidence", "") + f" [cited source {ids_now[src]['label']} is this node: {a['why']}]").strip()
        tgt.setdefault("merged_sources", []).append({"id": src, "pmid": ids_now[src].get("pmid", ""), "label": ids_now[src]["label"]})
        nodes.remove(ids_now.pop(src))
    seen_pair, kept = set(), []
    for e in edges:
        if (e["from"], e["to"]) not in seen_pair and e["from"] != e["to"]:
            seen_pair.add((e["from"], e["to"]))
            kept.append(e)
    edges[:] = kept
    for src, f in sf["fold"].items():
        s = ids_now.get(src)
        if not s:
            print(f"source_folds fold {src}: not in the graph")
            continue
        s["folded"] = f["why"]
        for e in edges:
            if e["from"] == src and e["to"] in ids_now:
                ids_now[e["to"]].setdefault("other_data", []).append(
                    {"label": s["label"], "pmid": s.get("pmid", ""), "why": f["why"], "kind": e["kind"]})
    linked = {e["from"] for e in edges} | {e["to"] for e in edges}
    loose = sorted(x["label"] for x in nodes if x["id"].startswith("ext:") and not x.get("folded") and x["id"] in linked)
    if loose:
        print(f"cited-paper source nodes not in source_folds.json (drawn as sources): {loose}")

    # A source node whose only links were dropped above ('same subjects' links, drop_same_subjects) anchors nothing; the
    # citation is already on the citing study under cited_sources.
    anchored = {e["from"] for e in edges} | {e["to"] for e in edges}
    orphans = [x["id"] for x in nodes if x.get("external") and x["id"] not in anchored]
    nodes = [x for x in nodes if x["id"] not in set(orphans)]
    if orphans:
        print(f"source nodes with no link left, dropped: {orphans}")
    # every author, for the page's search (Paul, 2026-09-26: searching "mccarroll" did not find Ling 2024)
    all_au = author_lists({n["pmid"] for n in nodes if n.get("pmid") and not n.get("external")},
                          {k: v.get("doi") for k, v in json.loads(CACHE.read_text()).items()} if CACHE.exists() else {})
    for n in nodes:
        if all_au.get(n.get("pmid")):
            # BibTeX accents ({\'a}, {\"o}) to plain letters, for search
            n["authors"] = "; ".join(sorted(re.sub(r"[{}]", "", re.sub(r"\{\\.\s*\{?(\w)\}?\}", r"\1", a)) for a in all_au[n["pmid"]]))
    # Edges in a fixed order, so a rebuild with unchanged inputs writes the same file.
    edges = sorted(edges, key=lambda e: (e["from"], e["to"], e["kind"], e.get("evidence", "")))

    inputs = [BACK_CSV, BACK_DEC, SETS, *map(Path, CORPUS_CSVS)] + [HERE / f for f in ("n_extracted.json", "sources_llm.json", "geo_links.json", "sources_class.json") if (HERE / f).exists()] \
        + [p for p in [DEGENES / "summary.tsv", DATA / "manifest.json", *sorted((DATA / "manifests").glob("*.json"))] if p.exists()]
    # No build timestamp in lineage.json: its content (and sha256) changes only when the graph does, so the analyses that
    # record its hash do not look stale after a rebuild. The time goes to lineage_meta.json and into the page.
    built = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out = {
        "script": "schizophrenia-transcriptome-replication/corpus/build_lineage.py",
        "inputs": {(str(p.relative_to(HERE.parent)) if p.is_relative_to(HERE.parent) else str(p.relative_to(RA)) if p.is_relative_to(RA) else str(p)): hashlib.sha256(p.read_bytes()).hexdigest()[:12] for p in inputs},
        "citations": {k: cites.get(k) for k in ("fetched", "source")} if cites else None,
        "themes": {k: th.get(k) for k in ("generated", "model", "prompt_sha256", "vocabulary")} if th else None,
        # corpus description for the page's "Corpus" panel (analyses/corpus_stats.py)
        "stats": json.loads((DEGENES / "corpus_stats.json").read_text()) if (DEGENES / "corpus_stats.json").exists() else None,
        "tech_tiers": [[k, tier_labels[k]] for k, _ in TIERS],
        "array_early_max": __import__("classify_technology").ARRAY_EARLY_MAX,
        "gene_index_meta": gindex_meta,
        "tour_genes": tour_genes(),
        "prune": prune_stats, "lanes": lanes, "nodes": sorted(nodes, key=lambda n: (n["year"], n["label"])), "edges": edges,
        "gene_index": gindex,
    }
    body = json.dumps(out, indent=1, ensure_ascii=False)
    (HERE / "lineage.json").write_text(body)
    (HERE / "lineage_meta.json").write_text(json.dumps({
        "generated": built, "script": out["script"], "lineage.json sha256": hashlib.sha256(body.encode()).hexdigest(),
        "note": "build time of lineage.json, kept out of it so its hash follows its content"}, indent=1) + "\n")
    template = (HERE / "lineage_template.html").read_text()
    # Bottom lane of the page: major grants and consortium announcements (milestones.json; Paul, 2026-09-27). Page only:
    # lineage.json stays as it is, since downstream analyses record its sha256.
    ms = json.loads(MILESTONES.read_text()) if MILESTONES.exists() else None
    if ms:
        ms = {**ms, "sha256": hashlib.sha256(MILESTONES.read_bytes()).hexdigest()[:12],
              "milestones": sorted(ms.get("milestones", []), key=lambda m: str(m.get("date", "")))}
    # GEO release dates of studies whose data were public before their paper (geo_release.py; Paul, 2026-09-28): page only
    gr = HERE / "geo_release.json"
    geo_rel = json.loads(gr.read_text())["studies"] if gr.exists() else {}
    # the page's published-brain-list filter (published_brain_list_ids): page only, like geo_release
    pub_brain = published_brain_list_ids(out["nodes"])
    single_gene_sweep(out["nodes"])  # keeps its papers.json / h4_papers.tsv agreement check; nothing goes to the page
    (HERE / "lineage.html").write_text(template.replace("/*__DATA__*/null", json.dumps({"generated": built, **page_view(out), "milestones": ms, "geo_release": geo_rel,
                                                                                        "published_brain_list_ids": pub_brain}, ensure_ascii=False)))
    print(f"{sum(not n.get('external') and not n.get('single_cell') and not n.get('landmark') for n in nodes)} studies + {sum(bool(n.get('single_cell')) for n in nodes)} single-cell + {sum(bool(n.get('external')) for n in nodes)} source nodes, {len(edges)} edges, "
          f"{len(lanes)} lanes; {sum(1 for n in nodes if n['n'])} with N")


if __name__ == "__main__":
    main()
