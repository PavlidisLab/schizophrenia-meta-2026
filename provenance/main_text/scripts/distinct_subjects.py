#!/usr/bin/env python3
"""How many distinct subjects (donors) the studies represent: the sum of each study's own n against a range for the
number of distinct people once shared donors are collapsed. Counts SUBJECTS (tracker, Vocabulary, 2026-09-25): a
study's n is taken as its subjects; subject overlap between studies is collapsed; samples and regions never enter.

Paul, 2026-09-27: "the executive summary of the data (and in the report) has to say how many different subjects are
represented in the [data]."

Study sets
  published brain   the 47 units of recurrence_units.tsv, scenario 'published, brain', modelled 'yes' (the report's
                    "47 published brain studies"); counts are the representative study's n / cases / controls in
                    h2/h2_studies.tsv, the basis of the report's naive total (report/numbers.tsv).
  first-look        the 204 first-look DE studies of h2/h2_studies.tsv (every tissue); studies that share one data set
                    (h2 column group, from primary_lists.json same_data_as) count once, at the group's largest count.
  A study with no recorded count for a measure is left out of that measure (listed in the meta file).

Measures: all (n), SZ (cases, or cases_sz where the study has one), CTRL (controls). In h2_studies.tsv n = cases +
  controls wherever all three are recorded. Where a list's case group pools psychosis diagnoses (N audit convention 2,
  analyses/list_extraction_audit/decisions_paul_2026-09-28.json), cases (and so n) counts every case in the
  contrast, and the schizophrenia count is cases_sz (corpus/node_field_overrides.json, carried onto the node by
  corpus/build_lineage.py); SZ reads cases_sz from corpus/lineage.json for those studies. The meta file lists them
  (sz_counts_from_cases_sz).

Cohort clusters
  published brain: the recurrence_units.tsv cluster column (recurrence.py cohort tier). first-look: the same rule from
  the lineage lane (recurrence.NAMED_LANES, recurrence.MERGE), split by tissue (brain vs each other tissue); a study
  with no named lane is its own cluster. Multi-bank studies (recurrence load_units 'multibank') have no cluster.

Lower end ("maximal overlap within a bank"). Within a cluster, every smaller study is assumed to be a subset of the
  largest, except where banks are known to hold different people: a cluster is a tree of sub-collections (BANK_TREE),
  and lower(node) = max(largest study placed at the node itself, sum of lower(child) over its children). Children
  are different banks or collections: Stanley Array Collection vs Neuropathology Consortium; within
  cmc+mssm+pitt+libd, CommonMind's banks (Mount Sinai, Penn, Pittsburgh) vs LIBD+HBCC, and Pittsburgh vs Mount Sinai.
  Studies spanning children stay at the parent unless their paper states the split (STATED_SPLITS). Clusters count as
  subject-separate (the cohort tier's rule); clusters joined by a measured overlap are merged and counted at the
  larger one. A multi-bank study adds nothing (it may lie inside the clusters it draws from).
Upper end ("collapse only measured overlap"). Sum of the counts, minus (a) the size of a fixed Stanley collection
  where the studies placed in it add up to more than it holds (overlap_graph.stanley_sizes; a study whose count
  exceeds its labelled collection is treated as collection not stated), then the stanley cluster as a whole is capped
  at the two collections together (Array Collection 35 SZ / 35 CTRL plus Neuropathology Consortium 15 / 15: 100, 50
  SZ, 50 CTRL; Paul, 2026-09-28: "Stanley collections were just those brains. 100 total or whatever."), whether or not
  a study states its collection, plus the donors a study states come from another Stanley collection (STANLEY_OUTSIDE,
  checked against the paper's text at run time), and (b) measured subject overlap between
  counted studies (subject_overlap_edges.tsv, identifier and fingerprint tiers), collapsed with Hunter's bound
  |union| <= sum |A_i| - sum over a spanning forest of |A_i & A_j| (maximum-weight forest). A fingerprint overlap
  enters as its excess over chance (k minus the chance mean of its source row), an identifier overlap as k; a pair
  measured in two banks adds the two; a pair measured in several region pairs takes the largest. Per diagnosis, only
  overlaps whose source row splits k by diagnosis (donor_overlap/overlap_summary.tsv) are subtracted. The LIBD
  BrainSeq donor table ('box:libd') stands for the union of the LIBD-bank studies counted; exact-crosswalk overlap
  between two LIBD-bank studies (donor_overlap/libd_crosswalk.tsv) is collapsed inside that union, by its own forest.

Exact crosswalks: every donor_overlap/*_crosswalk.tsv (columns biobank_id, dataset_or_study, study_pmid_if_known,
  diagnosis; at present katsel_crosswalk.tsv, Mount Sinai brain-bank IDs) is read as sets of donor IDs per study. For
  two counted studies a and b with ID sets S_a, S_b (a study's subjects lie within its set; slack = |S| - count when
  the set is larger than the study's cohort): shared IDs less both slacks are subtracted as an identifier overlap at the
  upper end; where b's count equals |S_b|, n_a + |S_b - S_a| is a lower end for the two together. A crosswalk that
  reaches CommonMind or PsychAD IDs would tighten the Mount Sinai part of cmc+mssm+pitt+libd; the Katsel IDs do not
  (katsel_crosswalk_report.md), so their effect is limited to Katsel 2005 and Wu 2021.

    opt/anaconda3/bin/python3 analyses/shared_subjects_between_datasets/distinct_subjects.py
      -> subject_relations/distinct_subjects.tsv, subject_relations/distinct_subjects_meta.json
"""
import csv
import datetime
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.path.dirname(HERE)
REPO = os.path.dirname(GL)
sys.path.insert(0, GL)
sys.path.insert(0, HERE)
import overlap_graph as OG  # noqa: E402  load(): units, Banks; stanley_sizes()
import recurrence as R  # noqa: E402  NAMED_LANES, MERGE

F = {
    "recurrence_units": os.path.join(GL, "recurrence_units.tsv"),
    "h2_studies": os.path.join(GL, "h2", "h2_studies.tsv"),
    "edges": os.path.join(HERE, "subject_overlap_edges.tsv"),
    "subject_relations": os.path.join(HERE, "subject_relations.tsv"),
    "overlap_summary": os.path.join(GL, "donor_overlap", "overlap_summary.tsv"),
    "provenance": os.path.join(REPO, "corpus", "provenance.json"),
    "lineage": os.path.join(REPO, "corpus", "lineage.json"),
    "primary_lists": os.path.join(GL, "primary_lists.json"),
}
CROSSWALKS = os.path.join(GL, "donor_overlap")  # *_crosswalk.tsv: exact donor-ID crosswalks (see docstring)
OUT_TSV = os.path.join(HERE, "distinct_subjects.tsv")
OUT_META = os.path.join(HERE, "distinct_subjects_meta.json")
MEASURES = (("all", "n"), ("SZ", "cases"), ("CTRL", "controls"))

# Sub-collections known to hold different people, per cluster: bank (subject_relations.Banks) -> path in the tree.
# CommonMind's DLPFC cohort is Mount Sinai, Penn and Pittsburgh brains; LIBD and NIMH HBCC share donors (Akula 2021,
# donor_overlap/), so they are one node; Pittsburgh and Mount Sinai are different banks (subject_relations.py rule 5).
BANK_TREE = {
    "cmc+mssm+pitt+libd": {"pitt": ("CommonMind banks", "Pittsburgh"), "mssm": ("CommonMind banks", "Mount Sinai"),
                           "libd": ("LIBD+HBCC",), "hbcc": ("LIBD+HBCC",)},
    "stanley": {},  # filled from Banks.collection: ('Array Collection',) / ('Neuropathology Consortium',)
}
STANLEY_NODE = {"array": "Array Collection", "npc": "Neuropathology Consortium"}
STANLEY_SIZE_KEY = {"Array Collection": "stanley:array", "Neuropathology Consortium": "stanley:npc"}

# Splits stated in the paper, checked against the paper's text at run time.
STATED_SPLITS = {
    "27668389": {  # Fromer 2016
        "text": "data/raw/27668389/fulltext/pmc_fulltext.txt",
        "parts": [
            # the 258 / 279 CommonMind DLPFC subjects by bank: Supplementary Table 1, "Brain Bank" rows (SCZ N, Control N)
            {"path": ("CommonMind banks", "Mount Sinai"), "label": "CommonMind, Mount Sinai", "SZ": 144, "CTRL": 158,
             "text": "data/raw/27668389/fulltext/supp_moesm35_layout.txt", "quote": "MSSM 144 0.56 158 0.57"},
            {"path": ("CommonMind banks", "Penn"), "label": "CommonMind, Penn", "SZ": 57, "CTRL": 37,
             "text": "data/raw/27668389/fulltext/supp_moesm35_layout.txt", "quote": "Penn 57 0.22 37 0.13"},
            {"path": ("CommonMind banks", "Pittsburgh"), "label": "CommonMind, Pittsburgh", "SZ": 57, "CTRL": 84,
             "text": "data/raw/27668389/fulltext/supp_moesm35_layout.txt", "quote": "Pitt 57 0.22 84 0.30"},
            {"path": ("LIBD+HBCC",), "label": "NIMH HBCC", "SZ": 131, "CTRL": 176,
             "quote": "from 131 SCZ cases and 176 controls"},
        ],
    },
    "27622934": {  # Darby 2016: two Stanley collections, SCZ and controls of each (the BPD and MDD groups are not counted)
        "text": "data/raw/27622934/fulltext/pmc_fulltext.txt",
        "parts": [
            # Table 1, hippocampus (cohort 1, Array Collection): SCZ 35, BPD 33, controls 32 (N = 100 after two outliers)
            {"path": ("Array Collection",), "label": "Array Collection, hippocampus", "SZ": 35, "CTRL": 32,
             "quote": "Number 35 33 32 13 14 15 15"},
            # Table 1, orbitofrontal cortex (replication cohort, Neuropathology Consortium): SCZ 13, BPD 14, MDD 15, controls 15
            {"path": ("Neuropathology Consortium",), "label": "Neuropathology Consortium, OFC", "SZ": 13, "CTRL": 15,
             "quote": "Number 35 33 32 13 14 15 15"},
        ],
    },
    "21183010": {  # Moreau 2011: Stanley Array Collection, BA9: 35 SCZ, 33 controls, 33 BPD (BPD not counted)
        "text": "data/raw/21183010/fulltext/pmc_fulltext.txt",
        "parts": [
            {"path": ("Array Collection",), "label": "Array Collection, BA9", "SZ": 35, "CTRL": 33,
             "quote": "Statistical analysis was performed using 35, 33, and 33 samples from the schizophrenia, control, and bipolar groups, respectively"},
        ],
    },
    "24475125": {  # Smalheiser 2014: Neuropathology Consortium, 15 per group; one schizophrenia and one bipolar sample removed
        "text": "data/raw/24475125/fulltext/pmc_fulltext.txt",
        "parts": [
            {"path": ("Neuropathology Consortium",), "label": "Neuropathology Consortium, BA10", "SZ": 14, "CTRL": 15,
             "quote": "Two samples were removed from consideration as unsuitable (one bipolar and one schizophrenia sample"},
        ],
    },
    "15018814": {  # Kuromitsu 2001: Neuropathology Consortium, 15 subjects in each of four groups (the 4 + 4 array pools are subsets)
        "text": "data/raw/15018814/fulltext/pdf_layout.txt",  # pdftotext -layout of the publisher PDF (data/raw/15018814/handfetch/)
        "parts": [
            {"path": ("Neuropathology Consortium",), "label": "Neuropathology Consortium, frontal cortex", "SZ": 15, "CTRL": 15,
             "quote": "15 samples of each frontal or temporal cortex of four matched groups of schizophrenia, bipolar disorder, major depression and normal controls"},
        ],
    },
}
# Stanley donors a study states come from outside the Array Collection and the Neuropathology Consortium: added to the
# stanley cluster's cap (upper end). Found by a grep of the lineage sample descriptions for other Stanley sources.
STANLEY_OUTSIDE = {
    "27898074": {  # Kim 2016, choroid plexus: the replication set
        "text": "data/raw/27898074/fulltext/pmc_fulltext.txt", "label": "New Stanley Collection", "SZ": 19, "CTRL": 19,
        "quote": "we sequenced the CP transcriptome of 19 subjects with schizophrenia and 19 unaffected controls from the New Stanley Collection",
    },
}
# Studies of one bank whose analysis also uses another bank's data: kept at the cluster root (not placed in a child).
SPANS_BANKS = {
    "30050107": {"text": "data/raw/30050107/fulltext/pmc_fulltext.txt",  # Jaffe 2018
                 "quote": "Meta-analysis between discovery (LIBD) and CMC"},
}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_tsv(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def num(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def rel(p):
    return os.path.relpath(p, GL)


def check_quotes():
    out = []
    for k, spec in list(STATED_SPLITS.items()) + list(SPANS_BANKS.items()) + list(STANLEY_OUTSIDE.items()):
        # a part may name its own source file (a supplement table); whitespace runs compare as one space, so a table row
        # extracted with pdftotext -layout matches its quote
        for f, q in [(p.get("text", spec["text"]), p["quote"]) for p in spec.get("parts", [])] or [(spec["text"], spec["quote"])]:
            text = " ".join(open(os.path.join(REPO, f), encoding="utf-8", errors="replace").read().split())
            if " ".join(q.split()) not in text:
                sys.exit(f"distinct_subjects: quote for {k} not found in {f}: {q!r}")
            out.append({"study": k, "file": f, "quote": q})
    return out


# ---------------------------------------------------------------------------------------------- study sets
def sz_counts():
    """{pmid: cases_sz} for the lineage nodes that carry cases_sz (pooled psychosis case groups; see the docstring)."""
    with open(F["lineage"], encoding="utf-8") as fh:
        nodes = json.load(fh)["nodes"]
    return {str(n.get("pmid") or n["id"]): n["cases_sz"] for n in nodes
            if n.get("cases_sz") is not None and not n.get("external")}


def count(s, col, sz):
    """A study's count for one measure: the SZ measure takes cases_sz where the study has one."""
    if col == "cases" and s.get("pmid") in sz:
        return sz[s["pmid"]]
    return num(s.get(col))


def published_brain(h2, sz):
    units = [r for r in read_tsv(F["recurrence_units"]) if r["scenario"] == "published, brain" and r["modelled"] == "yes"]
    items = []
    for u in units:
        k = u["representative"]
        s = h2.get(k, {"pmid": k})
        multi = u["cluster"].startswith("left out of cohort tier")
        items.append({"node": k, "pmid": k, "pmids": [k], "label": u["label"], "cluster": None if multi else u["cluster"],
                      "multibank": multi, "n": count(s, "n", sz), "cases": count(s, "cases", sz),
                      "controls": count(s, "controls", sz), "naive_rows": [k]})
    return items


def first_look(h2, D, sz):
    groups = {}
    for k, s in h2.items():
        groups.setdefault(s["group"] or k, []).append(s)
    multibank = D["units"]["multibank"]
    items = []
    for g, rows in sorted(groups.items()):
        rep = max(rows, key=lambda s: num(s["n"]) or -1)
        lane = rep["lane"]
        tissue = rep["tissue"] or "unknown"
        multi = any(s["pmid"] in multibank for s in rows)
        if lane in R.NAMED_LANES:
            cl = R.MERGE.get(lane, lane)
            cl = cl if tissue == "brain" else f"{cl} [{tissue}]"
        else:
            cl = f"own:{g}"
        it = {"node": rows[0]["pmid"] if len(rows) == 1 else g, "pmid": rep["pmid"], "pmids": [s["pmid"] for s in rows], "label": " + ".join(s["label"] for s in rows),
              "cluster": None if multi else cl, "multibank": multi, "naive_rows": [s["pmid"] for s in rows]}
        for _, col in MEASURES:
            vals = [count(s, col, sz) for s in rows if count(s, col, sz) is not None]
            it[col] = max(vals) if vals else None
            it["naive_" + col] = sum(vals) if vals else None
            it["rows_" + col] = len(vals)
        items.append(it)
    return items


# ---------------------------------------------------------------------------------------------- tree placement
def place(it, D, sizes):
    """Placements [(path, label, counts)] of one study in its cluster's tree; a study stated to split gives several."""
    k = it["pmid"]
    cl = it["cluster"]
    if k in STATED_SPLITS and cl in BANK_TREE:
        return [(p["path"], p["label"], {"cases": p["SZ"], "controls": p["CTRL"], "n": p["SZ"] + p["CTRL"]})
                for p in STATED_SPLITS[k]["parts"]]
    counts = {col: it[col] for _, col in MEASURES}
    if cl not in BANK_TREE or k in SPANS_BANKS or len(it["pmids"]) > 1 and cl == "cmc+mssm+pitt+libd":
        return [((), "", counts)]
    bank = D["banks"].bank(k)
    if cl == "stanley":
        node = STANLEY_NODE.get(D["banks"].collection(k))
        if node:
            szc, ctc, _ = sizes[STANLEY_SIZE_KEY[node]]
            if (counts["cases"] or 0) <= szc and (counts["controls"] or 0) <= ctc:
                return [((node,), "", counts)]
        return [((), "", counts)]
    return [(BANK_TREE[cl].get(bank, ()), "", counts)]


def lower_of(placements, col):
    """lower(node) = max(largest count at the node, sum over children of lower(child)); placements: [(path, count)]
    (an exact pair lower end enters as one more count at the node holding both studies)."""
    def rec(prefix):
        here = [c for p, c in placements if p == prefix]
        kids = sorted({p[len(prefix)] for p, _ in placements if len(p) > len(prefix) and p[:len(prefix)] == prefix})
        return max(max(here, default=0), sum(rec(prefix + (ch,)) for ch in kids))
    return rec(())


# ---------------------------------------------------------------------------------------------- measured overlap
def source_weights(src, SR, OS):
    """(bank key, k, excess over chance, {SZ, CTRL} excess or None) for one 'sources' entry of subject_overlap_edges.tsv."""
    m = re.match(r"subject_relations\.tsv row \d+ \((.+) / (.+)\)$", src.strip())
    if m:
        r = SR.get((m.group(1), m.group(2))) or SR.get((m.group(2), m.group(1)))
        if not r or num(r["k_shared"]) is None:
            return None
        k = num(r["k_shared"])
        lvl = r["fp_level"] or "exact"
        mu = float(r.get(f"fp_chance_mean_{lvl}") or 0) if r["evidence"] == "fingerprint" else 0.0
        return ("datasets", k, max(0.0, k - mu), None)
    m = re.match(r"overlap_summary\.tsv row \d+ \((\S+) / (\S+)\)$", src.strip())
    if m:
        r = OS.get((m.group(1), m.group(2))) or OS.get((m.group(2), m.group(1)))
        if not r:
            return None
        nid = num(r["same_id_pairs"]) or 0
        k = num(r["exact_shared"]) or 0
        banks = tuple(sorted(x.split(":")[1] for x in (r["unit_a"], r["unit_b"])))
        if nid > 0:
            return (banks, nid, float(nid), None)
        mu = float(r["exact_chance_mean"] or 0)
        exc = max(0.0, k - mu)
        f = exc / k if k else 0.0
        return (banks, k, exc, {"SZ": f * (num(r["exact_SZ"]) or 0), "CTRL": f * (num(r["exact_CTRL"]) or 0)})
    return None  # study_relations.tsv rows repeat a subject_relations.tsv / overlap_summary.tsv row


def load_edges():
    SR = {(r["dataset_a"], r["dataset_b"]): r for r in read_tsv(F["subject_relations"])}
    OS = {(r["unit_a"], r["unit_b"]): r for r in read_tsv(F["overlap_summary"])}
    rows = [r for r in read_tsv(F["edges"]) if r["tier"] in ("identifier", "fingerprint")]
    pair = {}
    for r in rows:
        a, b = sorted((r["a"], r["b"]))
        ws = [w for w in (source_weights(s, SR, OS) for s in re.split(r" \| |; ", r["sources"])) if w]
        byb = pair.setdefault((a, b), {})
        for bank, k, exc, dx in ws:  # several region pairs of one bank pair: the largest
            if bank not in byb or exc > byb[bank]["all"]:
                byb[bank] = {"all": exc, "k": k, "SZ": None if dx is None else dx["SZ"], "CTRL": None if dx is None else dx["CTRL"],
                             "tier": r["tier"], "text": r["text"]}
    out = []
    for (a, b), byb in pair.items():  # different bank pairs of one study pair: added
        out.append({"a": a, "b": b, "all": sum(v["all"] for v in byb.values()),
                    "SZ": sum(v["SZ"] or 0 for v in byb.values()), "CTRL": sum(v["CTRL"] or 0 for v in byb.values()),
                    "k": sum(v["k"] for v in byb.values()), "tiers": sorted({v["tier"] for v in byb.values()}),
                    "text": " + ".join(v["text"] for v in byb.values())})
    return out


def crosswalk_files():
    return sorted(os.path.join(CROSSWALKS, f) for f in os.listdir(CROSSWALKS) if f.endswith("_crosswalk.tsv"))


def crosswalk_sets():
    """{pmid: {dx: set of donor IDs}} from the exact crosswalks; 'all' is the union. IDs are prefixed with the file so
    two banks' numbers never meet."""
    out = {}
    for f in crosswalk_files():
        tag = os.path.basename(f)
        for r in read_tsv(f):
            for k in filter(None, re.split(r"[;,\s]+", r["study_pmid_if_known"])):
                d = out.setdefault(k, {"all": set(), "SZ": set(), "CTRL": set()})
                i = f"{tag}:{r['biobank_id']}"
                d["all"].add(i)
                if r["diagnosis"] in d:
                    d[r["diagnosis"]].add(i)
    return out


def crosswalk_pairs(items, cw):
    """Per counted pair with crosswalk sets: guaranteed shared subjects (upper end) and exact pair lower ends."""
    MC = {"all": "n", "SZ": "cases", "CTRL": "controls"}
    have = [it for it in items if len(it["pmids"]) == 1 and it["pmid"] in cw]
    edges, pair_lower = [], []
    for i, a in enumerate(have):
        for b in have[i + 1:]:
            e = {"a": a["pmid"], "b": b["pmid"], "tiers": ["identifier (crosswalk)"], "k": 0, "text": "exact crosswalk"}
            for m, col in MC.items():
                Sa, Sb = cw[a["pmid"]][m], cw[b["pmid"]][m]
                na, nb = a[col], b[col]
                if na is None or nb is None or not Sa or not Sb:
                    e[m] = 0.0
                    continue
                slack = max(0, len(Sa) - na) + max(0, len(Sb) - nb)
                e[m] = float(max(0, len(Sa & Sb) - slack))
                if m == "all":
                    e["k"] = len(Sa & Sb)
                for x, y, Sx, Sy, nx, ny in ((a, b, Sa, Sb, na, nb), (b, a, Sb, Sa, nb, na)):
                    if ny == len(Sy) and nx <= len(Sx):  # y's subjects are exactly Sy; x's lie within Sx
                        pair_lower.append({"a": x["node"], "b": y["node"], "measure": m, "value": nx + len(Sy - Sx)})
            edges.append(e)
    return edges, pair_lower


def forest(nodes, edges, key):
    """Maximum-weight spanning forest (Kruskal); returns the edges used."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    used = []
    for e in sorted(edges, key=lambda e: -e[key]):
        if e[key] <= 0:
            continue
        ra, rb = find(e["na"]), find(e["nb"])
        if ra != rb:
            parent[ra] = rb
            used.append(e)
    return used


def map_edges(items, edges, D):
    """Edges between counted super-nodes: a study maps to its item; box:libd to the union of LIBD-bank items."""
    of = {}
    for it in items:
        for p in it["pmids"]:
            of[p] = it["node"]
    libd = {it["node"] for it in items if it["cluster"] and len(it["pmids"]) == 1 and D["banks"].bank(it["pmid"]) == "libd"}
    box = "box:libd" if libd else None
    out = []
    for e in edges:
        na, nb = (box if x == "box:libd" or of.get(x) in libd else of.get(x) for x in (e["a"], e["b"]))
        if not na or not nb or na == nb:
            continue
        out.append(dict(e, na=na, nb=nb))
    return out, libd


# ---------------------------------------------------------------------------------------------- compute
def compute(set_name, items, D, sizes, edges, cw):
    rows, left_out = [], {}
    cw_edges, pair_lower = crosswalk_pairs(items, cw)
    emap, libd = map_edges(items, edges + [e for e in cw_edges if e["all"] > 0], D)
    # exact-crosswalk overlap between two LIBD-bank studies: both map to box:libd, so map_edges drops it; it is collapsed
    # inside the box instead (donor_overlap/libd_crosswalk.tsv)
    of_node = {p: it["node"] for it in items for p in it["pmids"]}
    in_box = [dict(e, na=of_node[e["a"]], nb=of_node[e["b"]], in_box=True) for e in cw_edges
              if e["all"] > 0 and of_node.get(e["a"]) in libd and of_node.get(e["b"]) in libd]
    clusters = sorted({it["cluster"] for it in items if it["cluster"]})
    cl_of = {it["node"]: it["cluster"] for it in items}
    for it in items:
        if it["node"] in libd:
            cl_of["box:libd"] = it["cluster"]
    for m, col in MEASURES:
        counted = [it for it in items if it[col] is not None]
        left_out[m] = sorted(it["label"] for it in items if it[col] is None)
        ids = {it["node"] for it in counted}
        sn = lambda n: "box:libd" if n in libd else n  # noqa: E731  super-node
        naive_col = "naive_" + col
        per_cl = {}
        for c in clusters + [None]:
            mem = [it for it in counted if it["cluster"] == c]
            if not mem:
                continue
            pl = [(path, cnt[col]) for it in mem for path, _, cnt in place(it, D, sizes)]
            pos1 = {it["node"]: place(it, D, sizes)[0][0] for it in mem if len(place(it, D, sizes)) == 1}
            pl += [(pos1[p["a"]], p["value"]) for p in pair_lower if p["measure"] == m and p["a"] in pos1
                   and pos1.get(p["b"]) == pos1[p["a"]]]
            lower = lower_of(pl, col) if c else 0
            naive = sum(it.get(naive_col, it[col]) for it in mem)
            dedup = sum(it[col] for it in mem)  # stated same-data groups at their largest count
            cap_saved = 0
            if c == "stanley":
                for node, key in STANLEY_SIZE_KEY.items():
                    at = [cnt for path, cnt in pl if path == (node,)]
                    szc, ctc, _ = sizes[key]
                    cap = {"cases": szc, "controls": ctc, "n": szc + ctc}[col]
                    cap_saved += max(0, sum(at) - cap)
            capped = {node for node in STANLEY_SIZE_KEY} if c == "stanley" else set()
            nodes_c = {sn(it["node"]) for it in mem}
            ee = [e for e in emap if e["na"] in nodes_c and e["nb"] in nodes_c and e["na"] in {sn(i) for i in ids}
                  and e["nb"] in {sn(i) for i in ids}]
            if capped:  # an overlap inside one capped collection is already inside its cap
                pos = {it["node"]: path for it in mem for path, _, _ in place(it, D, sizes)}
                ee = [e for e in ee if not (pos.get(e["na"]) and pos.get(e["na"]) == pos.get(e["nb"]) and pos[e["na"]][0] in capped)]
            used = forest(nodes_c, ee, m)
            box_ids = {it["node"] for it in mem if it["node"] in libd}
            used += forest(box_ids, [e for e in in_box if e["na"] in box_ids and e["nb"] in box_ids], m)
            ov = sum(e[m] for e in used)
            upper = dedup - cap_saved - ov
            if c == "stanley":  # the cluster's union cannot exceed the two collections, plus donors stated to be from outside
                whole = sum({"cases": sizes[k][0], "controls": sizes[k][1], "n": sizes[k][0] + sizes[k][1]}[col]
                            for k in STANLEY_SIZE_KEY.values())
                whole += sum({"cases": o["SZ"], "controls": o["CTRL"], "n": o["SZ"] + o["CTRL"]}[col]
                             for it in mem for k, o in STANLEY_OUTSIDE.items() if k in it["pmids"])
                if upper > whole:
                    cap_saved += upper - whole
                    upper = whole
            per_cl[c] = {"naive": naive, "lower": lower, "upper": upper, "members": mem, "used": used, "eligible": ee,
                         "cap_saved": cap_saved, "same_data": naive - dedup}
            rows.append({"set": set_name, "level": "cluster", "cluster": c or "multi-bank (no cluster)", "measure": m,
                         "n_studies": sum(it.get("rows_" + col, 1) for it in mem), "naive": naive, "lower": lower,
                         "upper": round(upper), "same_data_collapsed": naive - dedup, "capacity_collapsed": cap_saved,
                         "overlap_collapsed": round(ov, 1),
                         "overlaps": "; ".join(f"{e['na']}~{e['nb']} {e[m]:.1f}" for e in used),
                         "members": "; ".join(f"{it['label']} ({it[col]})" for it in mem)})
        # totals: clusters joined by a measured overlap are merged for the lower end (counted at the larger)
        par = {c: c for c in per_cl if c}

        def find(x):
            while par[x] != x:
                x = par[x]
            return x
        cross = [e for e in emap if cl_of.get(e["na"]) and cl_of.get(e["nb"]) and cl_of[e["na"]] != cl_of[e["nb"]]
                 and e["na"] in {sn(i) for i in ids} and e["nb"] in {sn(i) for i in ids}]
        for e in cross:
            if cl_of[e["na"]] in par and cl_of[e["nb"]] in par:
                par[find(cl_of[e["na"]])] = find(cl_of[e["nb"]])
        comps = {}
        for c in par:
            comps.setdefault(find(c), []).append(c)
        lower_total = sum(max(per_cl[c]["lower"] for c in cs) for cs in comps.values())
        # upper: one forest over every measured edge between counted studies (within clusters, across clusters and to
        # multi-bank studies); the within-cluster forests are part of it
        live = {sn(i) for i in ids}
        within = [e for c in per_cl.values() for e in c["used"] if not e.get("in_box")]  # box-internal: already in their cluster
        inside = {(e["na"], e["nb"]) for c in per_cl.values() for e in c.get("eligible", [])}
        other = [e for e in emap if e["na"] in live and e["nb"] in live and cl_of.get(e["na"]) != cl_of.get(e["nb"])
                 or e["na"] in live and e["nb"] in live and cl_of.get(e["na"]) is None]
        used_all = forest(live, within + [e for e in other if (e["na"], e["nb"]) not in inside], m)
        ov_cross = sum(e[m] for e in used_all) - sum(e[m] for e in within)
        naive_total = sum(v["naive"] for v in per_cl.values())
        upper_total = sum(v["upper"] for v in per_cl.values()) - ov_cross
        rows.append({"set": set_name, "level": "total", "cluster": "all", "measure": m,
                     "n_studies": sum(it.get("rows_" + col, 1) for it in counted), "naive": naive_total,
                     "lower": lower_total, "upper": round(upper_total),
                     "same_data_collapsed": sum(v["same_data"] for v in per_cl.values()),
                     "capacity_collapsed": sum(v["cap_saved"] for v in per_cl.values()),
                     "overlap_collapsed": round(sum(sum(e[m] for e in v["used"]) for v in per_cl.values()) + ov_cross, 1),
                     "overlaps": f"cross-cluster {ov_cross:.1f}; merged for the lower end: "
                                 + ("; ".join("+".join(cs) for cs in comps.values() if len(cs) > 1) or "none"),
                     "members": f"{len(per_cl) - (1 if None in per_cl else 0)} clusters"
                                + (f" + {len(per_cl[None]['members'])} multi-bank" if None in per_cl else "")})
    return rows, left_out, cw_edges, pair_lower


COLS = ["set", "level", "cluster", "measure", "n_studies", "naive", "lower", "upper", "same_data_collapsed",
        "capacity_collapsed", "overlap_collapsed", "overlaps", "members"]


def main():
    quotes = check_quotes()
    D = OG.load()
    sizes = OG.stanley_sizes(D)
    h2 = {r["pmid"]: r for r in read_tsv(F["h2_studies"])}
    sz = sz_counts()
    edges = load_edges()
    cw = crosswalk_sets()
    rows, left, cwr = [], {}, {}
    for name, items in (("published brain", published_brain(h2, sz)), ("first-look", first_look(h2, D, sz))):
        r, lo, ce, pl = compute(name, items, D, sizes, edges, cw)
        rows += r
        left[name] = lo
        cwr[name] = {"pairs": [{k: e[k] for k in ("a", "b", "k", "all", "SZ", "CTRL")} for e in ce], "pair_lower": pl}
    inputs = {rel(p): sha(p) for p in F.values()}
    for f in crosswalk_files():
        inputs[rel(f)] = sha(f)
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    with open(OUT_TSV, "w", encoding="utf-8", newline="") as fh:
        fh.write(f"# built {stamp} by subject_relations/distinct_subjects.py; method in the script docstring\n")
        fh.write("# naive = sum of each study's own count; lower = maximal overlap within each bank or collection; "
                 "upper = only measured overlap and fixed collection sizes collapsed. Subjects, not samples.\n")
        for p, h in inputs.items():
            fh.write(f"# input: {p} sha256 {h}\n")
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    meta = {"script": "subject_relations/distinct_subjects.py", "run": stamp, "inputs": inputs,
            "outputs": {rel(OUT_TSV): sha(OUT_TSV)},
            "stanley_sizes": {k: {"SZ": v[0], "CTRL": v[1], "source": v[2]} for k, v in sizes.items()},
            "stated_splits_and_spans_checked": quotes,
            "left_out_no_count": left,
            "sz_counts_from_cases_sz": {k: {"label": h2[k]["label"], "cases": num(h2[k]["cases"]), "cases_sz": v}
                                        for k, v in sorted(sz.items()) if k in h2},
            "exact_crosswalks": {"files": [rel(f) for f in crosswalk_files()], "studies": sorted(cw),
                                 "counted_pairs": cwr}}
    json.dump(meta, open(OUT_META, "w"), indent=1)
    for r in rows:
        if r["level"] == "total" or r["set"] == "published brain":
            print(f"{r['set'][:15]:15} {r['cluster'][:40]:40} {r['measure']:4} n_st={r['n_studies']:<4} naive={r['naive']:<6} "
                  f"lower={r['lower']:<6} upper={r['upper']:<6} cap={r['capacity_collapsed']} ov={r['overlap_collapsed']}")


if __name__ == "__main__":
    main()
