#!/usr/bin/env python3
"""Subject and sample overlap among the first-look SZ vs control brain studies, drawn as one node-edge figure.

Paul, 2026-09-25: "for the sample/subject overlap, we need a graph (nodes-edges) or something. something to explain the
situation across all these 'primary' studies". Vocabulary: tracker, "Vocabulary: subjects, samples, regions".

Nodes
  - the first-look primary studies of recurrence.load_units() (studies left after meta-lane, multi-series and reanalysis
    studies are removed) whose primary list is from brain;
  - single-nucleus first looks of the lineage (relation 'single-cell / single-nucleus', brain, not the meta lane), which
    have no primary list here but carry most of the donor-table evidence (square markers);
  - any other study that a measured relation joins to one of the above (open markers, e.g. a reanalysis whose dataset
    shares subjects with a first look).
  Blood, LCL, iPSC and other first-look units are summarized in the report table, not drawn.

Bank (box) of a node: subject_relations.Banks (the project's one bank classifier: lineage lane, Akula 2021 as HBCC,
Stanley split by sub-collection); for nodes Banks leaves without a bank, the donor table's single bank
(donor_overlap/tables/INDEX.tsv) when there is exactly one; else the lineage lane when it names a bank; else the box
'no single named bank'. Node outline: solid = bank from the lane as stated in the paper; dashed = lane inferred from an
author (lineage lane_by 'author: ...'); dotted = bank read from the donor table.

Edges, one colour per evidence tier, never pooled:
  identifier   the same subject identifiers in both datasets (subject_relations.tsv, study_relations.tsv evidence
               'identifier'; donor_overlap/overlap_summary.tsv same_id_pairs > 0)
  fingerprint  demographic fingerprint above chance (subject_relations.tsv / study_relations.tsv evidence 'fingerprint';
               overlap_summary.tsv rows passing subject_relations.fp_overlap without its permutation p, which that
               file does not carry: k >= 3, k > chance 95th percentile, k - chance mean >= 3, k >= 2 x chance mean)
  stated       stated in a paper or curated by hand: study_relations.tsv evidence 'stated', primary_lists.json
               same_data_as between drawn studies, lineage.json edges marked curated with a subject-vocabulary kind
  citation     lineage.json 'same dataset' edges between two drawn studies that are not marked
               curated (grey, dashed; not verified as a subject relation)
An edge whose evidence is a bank's donor table rather than one study (the LIBD BrainSeq table) points at the bank box.
A red '+Bank' tag after a label names another bank the study also draws donors from (lineage 'banks', donor-table
banks, curated 'same collection' edges from a consortium node).

Pairs of drawn studies with no measured relation are counted, not drawn: same bank, untested (overlap unknown); bank not
stated for one or both, untested (overlap unknown); different named banks, untested (presumed separate by bank only).

    opt/anaconda3/bin/python3 analyses/shared_subjects_between_datasets/overlap_graph.py
      -> subject_overlap_graph.svg / .png, subject_overlap_graph_meta.json, subject_overlap_edges.tsv,
         subject_overlap_pairs.tsv, subject_overlap_report.md (all in this directory)
"""
import collections
import csv
import datetime
import hashlib
import itertools
import json
import math
import os
import re
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.path.dirname(HERE)
REPO = os.path.dirname(GL)
sys.path.insert(0, GL)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser(".claude/skills/architecture-figures/python"))
import recurrence as R  # noqa: E402
import build_universes as BU  # noqa: E402
import subject_relations as S  # noqa: E402  Banks, NAMED, BANK_SAME, fingerprint thresholds
import relations as SR  # noqa: E402  SEPARATE

TODAY = datetime.date.today().isoformat()
F = {
    "primary_lists": os.path.join(GL, "primary_lists.json"),
    "lineage": R.LINEAGE,
    "study_relations": os.path.join(HERE, "study_relations.tsv"),
    "subject_relations": os.path.join(HERE, "subject_relations.tsv"),
    "units": os.path.join(HERE, "subject_relations_units.tsv"),
    "donor_summary": os.path.join(GL, "donor_overlap", "overlap_summary.tsv"),
    "donor_index": os.path.join(GL, "donor_overlap", "tables", "INDEX.tsv"),
    # donor overlap of a study whose donor table cannot be committed (Nersisyan 2026, not to be redistributed, Paul
    # 2026-10-01): counts only, same rules (subject_relations.fingerprint / fp_overlap); see SIDE_DONOR_STUDY
    "nersisyan_donors": os.path.join(REPO, "side-studies", "nersisyan_2026", "donor_overlap_summary.tsv"),
    "provenance": os.path.join(REPO, "corpus", "provenance.json"),
    "bank_names": os.path.join(REPO, "corpus", "bank_names.json"),
}
OUT = {k: os.path.join(HERE, v) for k, v in {
    "svg": "subject_overlap_graph.svg", "png": "subject_overlap_graph.png", "meta": "subject_overlap_graph_meta.json",
    "edges": "subject_overlap_edges.tsv", "pairs": "subject_overlap_pairs.tsv", "report": "subject_overlap_report.md"}.items()}

# donor_overlap table -> study key (PMID) it documents. gemma_brainscope is Gemma's copy of Emani 2024's sample file;
# libd_benjamin2024 is the LIBD BrainSeq phenotype table (Benjamin 2024), drawn as the LIBD bank box; gse224683 is the
# GEO series of Schulmann 2023 (a reanalysis).
TABLE_STUDY = {"psychad2025": "42778699", "ruzicka2024": "38781388", "ling2024": "38448582", "batiuk2022": "36223459",
               "emani2024": "38781369", "gemma_brainscope": "38781369", "pietersen2014_GSE37981": "24702465",
               "pietersen2014_GSE46509": "24628518", "glatt2005": "16223876", "pai2019": "31053723",
               "mendizabal2019": "31288836", "hbcc_sgacc": "33558674", "libd_benjamin2024": "box:libd",
               "gse224683": "36932057"}
# side-study donor summaries (F key -> the study they describe); rows: evidence 'identifier' (other_pmid, k_shared) or
# 'fingerprint (...)' (other_pmid or the donor table in brackets, overlap_rule = subject_relations.fp_overlap with p)
SIDE_DONOR_STUDY = {"nersisyan_donors": "41639050"}
# subject_relations.tsv datasets with no study in subject_relations_units.tsv: PsychAD pseudobulk donors (the PsychAD
# single-nucleus data, SZ vs control DE in Lee 2025)
DATASET_STUDY = {"PsychAD_HBCC": "42778699", "PsychAD_MSSM": "42778699"}
# donor-table bank names -> bank box
TABLE_BANK = {"McLean": "harvard", "MSSM": "mssm", "HBCC": "hbcc", "LIBD": "libd", "Pittsburgh": "pitt"}
# lineage lanes that name a bank box (subject_relations.NAMED) and box titles
BOX_TITLE = {"stanley:array": "Stanley Array Collection", "stanley:npc": "Stanley Neuropathology Consortium",
             "stanley:unknown": "Stanley, sub-collection not stated", "aus": "NSW TRC / Victoria (Australia)",
             "pritzker": "Pritzker Consortium", "pitt": "University of Pittsburgh", "cmc": "CommonMind / PsychENCODE",
             "hbcc": "NIMH HBCC", "mssm": "Mount Sinai (MSSM)", "libd": "LIBD / BrainSeq",
             "harvard": "McLean (HBTRC)", "none": "No single named bank (not stated, or several banks)"}
# box layout: rows of (box, columns of nodes)
BANK_SHORT = {"stanley": "Stanley", "aus": "NSW/Victoria", "pritzker": "Pritzker", "pitt": "Pittsburgh", "cmc": "CommonMind",
              "hbcc": "HBCC", "mssm": "Mount Sinai", "libd": "LIBD", "harvard": "McLean"}
LAYOUT = [[("pritzker", 1), ("stanley:array", 2), ("stanley:npc", 2), ("stanley:unknown", 2), ("aus", 2)],
          [("pitt", 2), ("cmc", 1), ("mssm", 1), ("hbcc", 1), ("libd", 2)],
          [("harvard", 2), ("none", 4)]]
SUBJECT_KINDS = {"same subjects", "same subjects, same region, different samples", "same subjects, different cell type",
                 "same subjects, different tissue", "partial subject overlap", "same samples (reanalysis)"}
TIERS = ["identifier", "fingerprint", "stated", "citation"]
TIER_LABEL = {"identifier": "identifier match", "fingerprint": "demographic fingerprint above chance",
              "stated": "stated in a paper or curated", "citation": "lineage 'same dataset' link, not curated"}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def key_of(n):
    return n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fix_label(label):
    """lineage labels title-case every word; roman-numeral name suffixes (Enwright III) are upper case"""
    return re.sub(r"\b(Ii|Iii|Iv)\b", lambda m: m.group(1).upper(), label)


def stanley_sizes(D):
    """{box: (n_SZ, n_control, source)} for the two fixed Stanley collections, read from repo files; a box is left out
    when no file states its size. Array Collection: the Puvogel 2022 methods sentence in corpus/provenance.json.
    Neuropathology Consortium: the curated lineage edge from coll:stanley_neuropathology_consortium ('15 per group')."""
    out = {}
    prov = json.load(open(F["provenance"]))
    for k, recs in prov.items():
        for r in recs if isinstance(recs, list) else []:
            m = re.search(r"(\d+) schizophrenia and (\d+) controls of the SMRI Array Collection", r.get("sentence", ""))
            if m and "stanley:array" not in out:
                lab = D["node"].get(k, {}).get("label", k)
                out["stanley:array"] = (int(m.group(1)), int(m.group(2)),
                                        f"{lab}, methods: \"{m.group(0)}\" (corpus/provenance.json, {k})")
    for e in D["lin"]["edges"]:
        if e["from"] == "coll:stanley_neuropathology_consortium" and e.get("curated"):
            m = re.search(r"(\d+) per group", e.get("evidence", ""))
            if m:
                n = int(m.group(1))
                out["stanley:npc"] = (n, n, f"{e['evidence']} (corpus/lineage.json, curated edge to {D['id2key'][e['to']]})")
    return out


def fp_pass(k, mu, q):
    """subject_relations.fp_overlap without the permutation p (overlap_summary.tsv does not carry it)."""
    return k >= S.MIN_K and k > q and k - mu >= S.MIN_EXCESS and k >= S.MIN_RATIO * mu


# ------------------------------------------------------------------------------------------------ data
def load():
    pl = json.load(open(F["primary_lists"]))
    lin = json.load(open(F["lineage"]))
    nodes = lin["nodes"]
    node = {key_of(n): n for n in nodes}
    id2key = {n["id"]: key_of(n) for n in nodes}
    units = R.load_units(pl, nodes)
    ctx = types.SimpleNamespace(units=units, node=node, pl=pl, specs=BU.manifest_specs())
    return dict(pl=pl, lin=lin, node=node, id2key=id2key, units=units, banks=S.Banks(ctx),
                sr=read_tsv(F["study_relations"]), dr=read_tsv(F["subject_relations"]), su=read_tsv(F["units"]),
                ds=read_tsv(F["donor_summary"]), di=read_tsv(F["donor_index"]),
                side={k: read_tsv(F[k]) for k in SIDE_DONOR_STUDY if os.path.exists(F[k])})


def build(D):
    node, units, B = D["node"], D["units"], D["banks"]
    studies = units["studies"]
    V = {}  # key -> {label, year, kind, bank, coll, bank_src, extra_banks}

    def add(k, kind):
        if k in V or k not in node:
            return
        n = node[k]
        V[k] = {"key": k, "label": fix_label(n["label"]), "year": int(n.get("year") or 0), "kind": kind,
                "lane": n.get("lane"), "lane_by": n.get("lane_by", ""), "relation": n.get("relation", ""),
                "extra": set(), "subject_data": False}

    brain_units = [k for k, v in studies.items() if v.get("tissue") == "brain"]
    for k in brain_units:
        add(k, "unit")
    sn = [key_of(n) for n in D["lin"]["nodes"] if n.get("relation") == "single-cell / single-nucleus"
          and n.get("tissue") == "brain" and n.get("lane") != "meta"]
    for k in sn:
        add(k, "single-nucleus")

    # ---- donor tables: banks per study
    table_banks = collections.defaultdict(set)
    for r in D["di"]:
        st = TABLE_STUDY.get(r["table"])
        if st and not st.startswith("box:"):
            table_banks[st].add(r["bank"])

    # ---- edges
    E = []

    def edge(a, b, tier, relation, text, source, k=""):
        if a == b:
            return
        k = str(int(float(k))) if str(k).replace(".", "", 1).isdigit() else ""
        E.append({"a": a, "b": b, "tier": tier, "relation": relation, "text": text, "source": source, "k": k})

    unit_study = {r["unit"]: r["study"] for r in D["su"] if r["kind"] == "dataset" and r["study"]}
    unit_study.update(DATASET_STUDY)
    for r in D["su"]:
        if r["kind"] == "dataset" and r["study"]:
            if r["study"] in V:
                V[r["study"]]["subject_data"] = True
    sep = SR.SEPARATE
    unknown = "same collection, subject overlap unknown"
    tested_sep = []  # (a, b, evidence, source)
    # study_relations.tsv
    for i, r in enumerate(D["sr"]):
        a, b = r["unit_a"].split("__")[0], r["unit_b"].split("__")[0]
        src = f"study_relations.tsv row {i + 1} ({r['label_a']} / {r['label_b']})"
        if r["relation"] in (sep, unknown):
            if r["evidence"] in ("identifier", "fingerprint") and a != b:
                tested_sep.append((a, b, r["evidence"], src, r["evidence_detail"]))
            continue
        edge(a, b, r["evidence"] if r["evidence"] in TIERS else "stated", r["relation"], r["relation_text"], src, r["k_shared"])
    # subject_relations.tsv (dataset pairs)
    for i, r in enumerate(D["dr"]):
        a, b = unit_study.get(r["dataset_a"]), unit_study.get(r["dataset_b"])
        if not a or not b or a == b:
            continue
        src = f"subject_relations.tsv row {i + 1} ({r['dataset_a']} / {r['dataset_b']})"
        if r["relation"] in (sep, unknown):
            if r["evidence"] in ("identifier", "fingerprint"):
                tested_sep.append((a, b, r["evidence"], src, r["evidence_detail"]))
            continue
        tier = r["evidence"] if r["evidence"] in TIERS else ("identifier" if r["evidence"] == "same samples" else "stated")
        edge(a, b, tier, r["relation"], f"{r['relation_text']} [{r['dataset_a']} / {r['dataset_b']}]", src, r["k_shared"])
    # donor_overlap/overlap_summary.tsv
    for i, r in enumerate(D["ds"]):
        ta, ba = r["unit_a"].split(":")
        tb, bb = r["unit_b"].split(":")
        a, b = TABLE_STUDY.get(ta), TABLE_STUDY.get(tb)
        if not a or not b or a == b:
            continue
        k, mu, q, nid = num(r["exact_shared"]), num(r["exact_chance_mean"]), num(r["exact_chance_q95"]), num(r["same_id_pairs"])
        src = f"overlap_summary.tsv row {i + 1} ({r['unit_a']} / {r['unit_b']})"
        detail = f"{bb if a.startswith('box:') else ba} vs {bb}: {int(k)} matched, chance mean {mu:g}, 95th percentile {q:g}"
        if nid > 0:
            edge(a, b, "identifier", "subject overlap", f"{int(nid)} shared identifiers ({ba} / {bb}); fingerprint {int(k)} matched", src, nid)
        elif fp_pass(k, mu, q):
            edge(a, b, "fingerprint", "partial subject overlap", f"{ba} / {bb}: {int(k)} matched, chance mean {mu:g}", src, k)
        elif not (a.startswith("box:") or b.startswith("box:")):
            tested_sep.append((a, b, "fingerprint (donor table)", src, detail))
    # side-study donor summaries (no committed donor table)
    for fk, st in SIDE_DONOR_STUDY.items():
        for i, r in enumerate(D["side"].get(fk, [])):
            ev = r["evidence"]
            src = f"{os.path.relpath(F[fk], REPO)} row {i + 1} ({r['other']})"
            if ev == "identifier" and r["other_pmid"]:
                if num(r["k_shared"]) > 0:
                    edge(st, r["other_pmid"], "identifier", "partial subject overlap",
                         f"{int(num(r['k_shared']))} shared CommonMind individual IDs ({r['nersisyan_unit']})", src, r["k_shared"])
            elif ev.startswith("fingerprint") and not r["relation"].startswith(("null", "bank unknown")):
                tb = r["other"].rsplit("[", 1)[-1].split(":")[0]
                b = r["other_pmid"] or TABLE_STUDY.get(tb)
                if not b:
                    continue
                if r["overlap_rule"] == "yes":
                    edge(st, b, "fingerprint", "partial subject overlap",
                         f"{r['nersisyan_unit']} / {r['other'].rsplit('[', 1)[-1].rstrip(']')}: {r['k_shared']} matched, chance mean {r['chance_mean']}",
                         src, r["k_shared"])
                elif not b.startswith("box:"):
                    tested_sep.append((st, b, "fingerprint (donor table)", src,
                                       f"{r['k_shared']} matched, chance mean {r['chance_mean']}, 95th percentile {r['chance_q95']}"))
    # primary_lists same_data_as
    for k, v in D["pl"]["studies"].items():
        for o in v.get("same_data_as", []):
            ok = o if isinstance(o, str) else o.get("pmid", "")
            ev = "" if isinstance(o, str) else o.get("evidence", "")
            edge(k, ok, "stated", "same data", ev, f"primary_lists.json {k} same_data_as {ok}")
    # lineage edges
    for e in D["lin"]["edges"]:
        a, b = D["id2key"].get(e["from"]), D["id2key"].get(e["to"])
        if e["from"].startswith("coll:"):
            continue
        if e.get("curated") and e["kind"] in SUBJECT_KINDS | {"same dataset", "same collection"}:
            if e["kind"] == "same collection":
                continue  # a collection relation, not a subject relation between the two studies
            edge(a, b, "stated", e["kind"], e.get("evidence", ""), f"lineage.json edge {e['from']} -> {e['to']} ({e['kind']}, curated)")
        elif e["kind"] == "same dataset":
            edge(a, b, "citation", "same dataset", e.get("evidence", ""), f"lineage.json edge {e['from']} -> {e['to']} (same dataset)")

    # keep edges that touch a drawn node; add the other end as an 'other' node when a measured or stated relation needs it
    # measured relations (identifier, fingerprint) may bring in another study as an open node; stated and citation
    # relations are drawn only between first looks (a reanalysis or meta-analysis of the same data is not a first look)
    first = lambda y: y.startswith("box:") or (y in V and V[y]["kind"] != "other")
    for x in E:
        if x["tier"] in ("identifier", "fingerprint") and (first(x["a"]) or first(x["b"])):
            for y in (x["a"], x["b"]):
                if not y.startswith("box:") and y not in V:
                    add(y, "other")
    keep = []
    for x in E:
        a, b = x["a"], x["b"]
        if x["tier"] in ("identifier", "fingerprint"):
            ok = (a in V or a.startswith("box:")) and (b in V or b.startswith("box:"))
        else:
            ok = first(a) and first(b) and not (a.startswith("box:") or b.startswith("box:"))
        if ok:
            keep.append(x)
    # labels unique among drawn studies: a, b by key when two share an author-year label
    seen = collections.defaultdict(list)
    for k in sorted(V, key=lambda k: (int(k) if k.isdigit() else 10 ** 9)):
        seen[V[k]["label"]].append(k)
    for labl, ks in seen.items():
        if len(ks) > 1:
            for i, k in enumerate(ks):
                V[k]["label"] = labl + "abcdefgh"[i]
    # donor tables mark subject-level data
    for st in table_banks:
        if st in V:
            V[st]["subject_data"] = True

    # ---- banks
    for k, v in V.items():
        b = B.bank(k)
        src = "lane"
        if b is None and v["kind"] == "other":
            db = {r["bank"] for r in D["su"] if r["kind"] == "dataset" and r["study"] == k} - {"not stated", ""}
            if len(db) == 1 and next(iter(db)) in S.NAMED:
                b, src = next(iter(db)), "dataset table"
        if b is None:
            allt = table_banks.get(k, set())
            if len(allt) == 1 and next(iter(allt)) in TABLE_BANK:
                b, src = TABLE_BANK[next(iter(allt))], "donor table"
        if b is not None and src == "lane" and v["lane_by"].startswith("author"):
            src = "author"
        if b == "stanley":
            v["box"] = f"stanley:{B.collection(k)}"
        else:
            v["box"] = b if b in BOX_TITLE else "none"
        v["bank"], v["bank_src"] = b, src
        extra = {x for x in (node[k].get("banks") or []) if x in S.NAMED}
        extra |= {TABLE_BANK[x] for x in table_banks.get(k, ()) if x in TABLE_BANK}
        v["extra"] = {x for x in extra if x != b}
        v["table_banks"] = sorted(table_banks.get(k, ()))
    for e in D["lin"]["edges"]:
        if e["from"].startswith("coll:") and e.get("curated") and e["kind"] == "same collection":
            k = D["id2key"].get(e["to"])
            coll = {"coll:cmc": "cmc", "coll:psychencode": "cmc", "coll:nimh_hbcc": "hbcc", "coll:libd_brainseq": "libd",
                    "coll:stanley_neuropathology_consortium": "stanley"}.get(e["from"])
            if k in V and coll and coll != V[k]["bank"]:
                V[k]["extra"].add(coll)
    # collapse: one study joined by the same evidence text to 2+ nodes of one box -> an edge to the box
    grp = collections.defaultdict(list)
    for x in keep:
        if x["b"] in V and x["a"] in V:
            grp[(x["a"], x["tier"], x["text"], V[x["b"]]["box"])].append(x)
            grp[(x["b"], x["tier"], x["text"], V[x["a"]]["box"])].append(dict(x, a=x["b"], b=x["a"]))
    collapsed, gone = [], set()
    for (a, tier, text, box), xs in grp.items():
        if len({x["b"] for x in xs}) >= 2 and V[a]["box"] != box:
            collapsed.append({"a": a, "b": "box:" + box, "tier": tier, "relation": xs[0]["relation"],
                              "text": text + f" (study_relations lists it for {len(xs)} studies of this bank)",
                              "source": "; ".join(sorted({x['source'] for x in xs})), "k": xs[0]["k"]})
            gone |= {frozenset((x["a"], x["b"])) for x in xs}
    keep = [x for x in keep if not (frozenset((x["a"], x["b"])) in gone and x["a"] in V and x["b"] in V)] + collapsed
    # one edge per unordered pair and tier; texts merged
    M = {}
    for x in keep:
        key = (frozenset((x["a"], x["b"])), x["tier"])
        if key in M:
            if x["text"] and x["text"] not in M[key]["texts"]:
                M[key]["texts"].append(x["text"])
            M[key]["sources"].append(x["source"])
            if x["k"] and x["k"] not in M[key]["ks"]:
                M[key]["ks"].append(x["k"])
        else:
            M[key] = dict(x, texts=[x["text"]] if x["text"] else [], sources=[x["source"]], ks=[x["k"]] if x["k"] else [])
    edges = sorted(M.values(), key=lambda x: (TIERS.index(x["tier"]), x["a"], x["b"]))
    # a citation edge between a pair that already has a stronger tier adds nothing
    strong = {frozenset((x["a"], x["b"])) for x in edges if x["tier"] != "citation"}
    edges = [x for x in edges if x["tier"] != "citation" or frozenset((x["a"], x["b"])) not in strong]
    return V, edges, tested_sep


# ------------------------------------------------------------------------------------------------ pair accounting
PAIR_CLASSES = ["overlap: identifier", "overlap: fingerprint", "same data or subjects: stated", "'same dataset': lineage link, not curated",
                "tested, no overlap found", "same bank, untested", "bank not stated, untested",
                "different named banks, untested (presumed separate)"]


def pair_class(a, b, V, E, tested):
    tiers = E.get(frozenset((a, b)), set())
    for t, c in (("identifier", 0), ("fingerprint", 1), ("stated", 2), ("citation", 3)):
        if t in tiers:
            return PAIR_CLASSES[c]
    if frozenset((a, b)) in tested:
        return PAIR_CLASSES[4]
    va, vb = V[a], V[b]
    ba, bb = {va["bank"]} | va["extra"] - {None}, {vb["bank"]} | vb["extra"] - {None}
    ba.discard(None)
    bb.discard(None)
    if not ba or not bb:
        return PAIR_CLASSES[6]
    if "stanley" in ba & bb:
        ca, cb = va["box"].split(":")[-1], vb["box"].split(":")[-1]
        if len(ba & bb) == 1 and "unknown" not in (ca, cb) and ca != cb:
            return PAIR_CLASSES[7]
        return PAIR_CLASSES[5]
    if ba & bb or any(frozenset((x, y)) in S.BANK_SAME for x in ba for y in bb):
        return PAIR_CLASSES[5]
    return PAIR_CLASSES[7]


def pairs(V, edges, tested_sep):
    E = collections.defaultdict(set)
    for x in edges:
        E[frozenset((x["a"], x["b"]))].add(x["tier"])
    tested = {frozenset((a, b)) for a, b, *_ in tested_sep}
    keys = sorted(k for k in V if V[k]["kind"] != "other")
    return [(a, b, pair_class(a, b, V, E, tested)) for a, b in itertools.combinations(keys, 2)]


# ------------------------------------------------------------------------------------------------ figure
def draw(V, edges, pngtoo=True, sizes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Circle, RegularPolygon, FancyArrowPatch
    from matplotlib.lines import Line2D
    from pavlab_arch.style import apply_rcparams
    from pavlab_arch import palette as P
    from matplotlib.colors import to_rgba
    PLATE = to_rgba(P.SOFT_BG, 0.78)
    from pavlab_arch.layout import svg_safe
    apply_rcparams()
    COL = {"identifier": P.ACCENT, "fingerprint": P.ACCENT_3, "stated": P.ACCENT_5, "citation": P.SUBTLE}
    FS, COLW, PITCH, HEAD, PAD, GAP = 7.0, 1.62, 0.235, 0.42, 0.14, 0.34
    boxes = collections.defaultdict(list)
    for k, v in V.items():
        boxes[v["box"]].append(k)
    for b in boxes:
        boxes[b].sort(key=lambda k: (V[k]["year"], V[k]["label"]))
    # geometry: rows top to bottom
    geo, y = {}, 0.0
    rows_h = []
    for row in LAYOUT:
        x, hmax = 0.0, 0.0
        for b, nc in row:
            n = len(boxes.get(b, []))
            nr = max(1, math.ceil(n / nc))
            w, h = nc * COLW + 2 * PAD, HEAD + nr * PITCH + PAD
            geo[b] = [x, None, w, h, nc]
            x += w + GAP
            hmax = max(hmax, h)
        rows_h.append((row, hmax, x - GAP))
    W = max(r[2] for r in rows_h)
    LEG_H = 1.7
    H = sum(r[1] for r in rows_h) + GAP * (len(rows_h) - 1) + LEG_H + 0.75
    ytop = H - 0.55
    for row, hmax, wrow in rows_h:
        # stretch each row to the full width: spread the slack over its boxes
        slack = (W - wrow) / len(row)
        x = 0.0
        for b, nc in row:
            g = geo[b]
            g[0], g[2] = x, g[2] + slack
            g[1], g[3] = ytop - hmax, hmax
            x += g[2] + GAP
        ytop -= hmax + GAP
    pos = {}
    for b, ks in boxes.items():
        x0, y0, w, h, nc = geo[b]
        colw = (w - 2 * PAD) / nc
        nr = math.ceil(len(ks) / nc)
        for i, k in enumerate(ks):
            c, r = divmod(i, nr)
            pos[k] = (x0 + PAD + c * colw + 0.08, y0 + h - HEAD - (r + 0.5) * PITCH)
    fig = plt.figure(figsize=(W + 0.3, H))
    ax = fig.add_axes([0.15 / (W + 0.3), 0, W / (W + 0.3), 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.text(0.012, 1 - 0.18 / H, "Subject and sample overlap among first-look SZ vs control brain studies",
             fontsize=13, color=P.TEXT, va="top")
    fig.text(0.012, 1 - 0.40 / H, "Boxes are brain banks; a line joins two studies only where a relation was measured or stated. "
             "Studies in one box without a line: overlap not known.", fontsize=8.5, color=P.SUBTLE, va="top")
    head_arts = []
    for b, (x0, y0, w, h, nc) in geo.items():
        ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=0.08", fc=P.SOFT_BG,
                                    ec=P.DET, lw=1.0, zorder=0))
        n = sum(V[k]["kind"] != "other" for k in boxes.get(b, []))
        head_arts.append(ax.text(x0 + PAD, y0 + h - 0.1, BOX_TITLE[b], fontsize=8.2, fontweight="bold", color=P.TEXT,
                                 va="top", zorder=5, bbox=dict(fc=PLATE, ec="none", pad=0.4)))
        ax.text(x0 + w - PAD, y0 + h - 0.1, f"{n}", fontsize=8, color=P.SUBTLE, va="top", ha="right", zorder=5)
        inc = sorted(BANK_SHORT[x] for pair in S.BANK_SAME if b in pair for x in pair if x != b)
        sub = ("includes donors of " + ", ".join(inc)) if inc else ""
        if sizes and b in sizes:
            sub = f"fixed collection, {sizes[b][0]} SZ and {sizes[b][1]} controls"
        if sub:
            head_arts.append(ax.text(x0 + PAD, y0 + h - 0.27, sub, fontsize=6.8, color=P.SUBTLE, va="top", style="italic",
                                     zorder=5, bbox=dict(fc=PLATE, ec="none", pad=0.4)))
    def box_anchor(b, toward):
        x0, y0, w, h, _ = geo[b]
        cx, cy = x0 + w / 2, y0 + h / 2
        tx, ty = toward
        dx, dy = tx - cx, ty - cy
        if dx == 0 and dy == 0:
            return cx, cy
        s = min((w / 2) / abs(dx) if dx else 1e9, (h / 2) / abs(dy) if dy else 1e9)
        return cx + dx * s, cy + dy * s
    # nodes and labels first: labels are opaque, so lines pass under them; their extents steer lines and badges
    RM, HALO = 0.058, 0.105
    lab_ext = {}
    fig.canvas.draw()
    inv = ax.transData.inverted()
    ext = lambda art: inv.transform(art.get_window_extent(fig.canvas.get_renderer()).get_points())
    for k, (x, y) in pos.items():
        v = V[k]
        fc = P.TEXT if v["kind"] == "unit" else (P.GRID if v["kind"] == "other" else P.SUBTLE)
        ec, lw = (P.SUBTLE if v["kind"] == "other" else P.TEXT), 0.9
        if v["bank_src"] == "author":  # hollow: bank inferred from an author name
            fc, ec, lw = "white", P.TEXT, 1.6
        if v["kind"] == "single-nucleus":
            m = RegularPolygon((x, y), 4, radius=0.075, orientation=math.pi / 4, fc=fc, ec=ec, lw=lw, zorder=4)
        else:
            m = Circle((x, y), RM, fc=fc, ec=ec, lw=lw, zorder=4)
        ax.add_patch(m)
        if v["bank_src"] in ("donor table", "dataset table"):  # ring: bank read from the study's own subject table
            ax.add_patch(Circle((x, y), HALO, fc="none", ec=P.SUBTLE, lw=0.9, zorder=4))
        lab = v["label"]
        if v["kind"] == "other":
            lab += {"meta-analysis / combined": " (meta)", "reanalysis of public data": " (reanalysis)",
                    "data reuse": " (data reuse)"}.get(v["relation"], "")
        t = ax.text(x + 0.13, y, lab, fontsize=FS, color=P.TEXT if v["kind"] != "other" else P.SUBTLE, va="center",
                    zorder=5, bbox=dict(fc=PLATE, ec="none", pad=0.4))
        arts = [t]
        if v["extra"]:
            arts.append(ax.annotate(" +" + ", ".join(BANK_SHORT[b] for b in sorted(v["extra"])), xy=(1, 0.5), xycoords=t,
                                    va="center", fontsize=6.3, color=P.ACCENT_4, zorder=5,
                                    bbox=dict(fc=PLATE, ec="none", pad=0.4)))
        fig.canvas.draw()
        pts = [ext(a) for a in arts]
        lab_ext[k] = (min(p[0][0] for p in pts), min(p[0][1] for p in pts), max(p[1][0] for p in pts), max(p[1][1] for p in pts))
    obstacles = [tuple(ext(a).ravel()) for a in head_arts] + list(lab_ext.values()) + [(x - HALO, y - HALO, x + HALO, y + HALO) for x, y in pos.values()]

    def anchor(k, toward):
        """a study's end of a line: its marker, or the right end of its label when the other end lies to the right"""
        x, y = pos[k]
        return (lab_ext[k][2] + 0.03, y) if toward[0] > lab_ext[k][2] + 0.2 else (x, y)

    def bez(pa, pb, rad, t):
        c = ((pa[0] + pb[0]) / 2 + rad * (pb[1] - pa[1]), (pa[1] + pb[1]) / 2 - rad * (pb[0] - pa[0]))
        return tuple((1 - t) ** 2 * pa[i] + 2 * (1 - t) * t * c[i] + t ** 2 * pb[i] for i in (0, 1))

    def free(bx):
        return not any(bx[0] < o[2] and o[0] < bx[2] and bx[1] < o[3] and o[1] < bx[3] for o in obstacles)

    # edges: one line per pair and measured tier; a stated line is left out where the pair has a measured one
    measured = {frozenset((x["a"], x["b"])) for x in edges if x["tier"] in ("identifier", "fingerprint")}
    # edges to draw: a stated/citation line is left out where the pair has a measured one
    draw = [x for x in edges if not (x["tier"] in ("stated", "citation") and frozenset((x["a"], x["b"])) in measured)]
    # bundle inter-box stated/citation edges per box pair (Paul 2026-10-02: the stated mesh is unreadable): each box
    # pair gets its own corridor offset and its edges fan within it; measured tiers keep their own dodged arcs
    bundle = collections.defaultdict(list)
    for x in draw:
        ba = V[x["a"]]["box"] if x["a"] in V else x["a"]
        bb = V[x["b"]]["box"] if x["b"] in V else x["b"]
        if x["tier"] in ("stated", "citation") and ba != bb:
            bundle[tuple(sorted((ba, bb)))].append(x)
    binfo = {}
    for gi_, key in enumerate(sorted(bundle)):
        n_ = len(bundle[key])
        for i_, x in enumerate(bundle[key]):
            binfo[id(x)] = (0.55 + 0.5 * (gi_ % 3), (-1) ** (gi_ % 2), i_, n_)
    nth = collections.Counter()
    for x in draw:
        a, b = x["a"], x["b"]
        pk = frozenset((a, b))
        nth[pk] += 1
        ca = pos[a] if a in pos else box_anchor(a[4:], pos[b])
        cb = pos[b] if b in pos else box_anchor(b[4:], pos[a])
        pa = anchor(a, cb) if a in pos else ca
        pb = anchor(b, pa) if b in pos else cb
        if a in pos:
            pa = anchor(a, pb)
        same_box = a in V and b in V and V[a]["box"] == V[b]["box"]
        ba_ = V[a]["box"] if a in V else a
        bb_ = V[b]["box"] if b in V else b
        bundled = x["tier"] in ("stated", "citation") and not same_box and ba_ != bb_ and id(x) in binfo
        lw = {"identifier": 2.0, "fingerprint": 1.6, "stated": 1.3, "citation": 1.0}[x["tier"]]
        al = {"identifier": 0.95, "fingerprint": 0.9, "stated": 0.5, "citation": 0.6}[x["tier"]]
        if bundled:
            # ribbon through the box pair's corridor: pass point at perpendicular offset d from the chord midpoint
            D, sgn, i_, n_ = binfo[id(x)]
            chord = (pb[0] - pa[0], pb[1] - pa[1])
            L = math.hypot(*chord) or 1.0
            mult = 1 if (ba_, bb_) == min((ba_, bb_), (bb_, ba_)) else -1
            d = sgn * mult * (D + 0.13 * (i_ - (n_ - 1) / 2))
            rad = 2 * d / L
            ax.add_patch(FancyArrowPatch(pa, pb, arrowstyle="-", connectionstyle=f"arc3,rad={rad}",
                                         color=COL[x["tier"]], lw=lw, zorder=2, alpha=al, shrinkA=0, shrinkB=0))
            continue
        rad = (0.55 if same_box else 0.18) * (1 if nth[pk] == 1 else -1)
        if same_box and abs(pb[0] - pa[0]) < 0.05:
            # same column: bow to the left, into the gap before the markers, away from the labels
            rad = -abs(rad) * (1 if pb[1] > pa[1] else -1)
        else:
            # bow to whichever side crosses fewer labels and box headings (the two ends' own labels aside)
            own = [lab_ext[k] for k in (a, b) if k in lab_ext]
            hits = lambda r: sum(any(o[0] < q[0] < o[2] and o[1] < q[1] < o[3] for o in obstacles if o not in own)
                                 for q in (bez(pa, pb, r, t / 40) for t in range(2, 39)))
            if hits(-rad) < hits(rad):
                rad = -rad
        ls = (0, (4, 2.5)) if x["tier"] == "citation" else "-"
        ax.add_patch(FancyArrowPatch(pa, pb, arrowstyle="-", connectionstyle=f"arc3,rad={rad}", color=COL[x["tier"]],
                                     lw=lw, linestyle=ls, zorder=2, alpha=al, shrinkA=0, shrinkB=0))
        if x["tier"] in ("identifier", "fingerprint") and x["ks"]:
            ks = sorted({int(k) for k in x["ks"]})
            txt = str(ks[0]) if len(ks) == 1 else f"{ks[0]}-{ks[-1]}"
            bw, bh = 0.045 * len(txt) + 0.09, 0.14
            # badge on the curve (data units are inches, aspect 1:1), at the first point clear of labels and badges
            cands = [(tt, 0.0) for tt in (0.5, 0.42, 0.58, 0.34, 0.66, 0.26, 0.74, 0.18, 0.82)]
            cands += [(tt, off) for off in (-0.1, 0.1, -0.18, 0.18) for tt in (0.5, 0.4, 0.6, 0.3, 0.7)]
            for tt, off in cands + [(0.5, 0.0)]:
                mx, my = bez(pa, pb, rad, tt)
                mx -= off  # off the curve sideways, toward the free side of a crowded column
                bx = (mx - bw / 2, my - bh / 2, mx + bw / 2, my + bh / 2)
                if free(bx):
                    break
            obstacles.append(bx)
            ax.text(mx, my, txt, fontsize=6.3, color=COL[x["tier"]], ha="center", va="center", zorder=6,
                    fontweight="bold", bbox=dict(fc="white", ec=COL[x["tier"]], lw=0.5, boxstyle="round,pad=0.15"))
    # legend
    lx, ly = 0.1, LEG_H - 0.05
    items = [("line", COL[t], t) for t in TIERS]
    for i, (_, c, t) in enumerate(items):
        yy = ly - 0.2 - i * 0.24
        ax.add_line(Line2D([lx, lx + 0.45], [yy, yy], color=c, lw=2.0 if t != "citation" else 1.0,
                           ls="-" if t != "citation" else (0, (4, 2.5))))
        ax.text(lx + 0.55, yy, f"{t}: {TIER_LABEL[t]}", fontsize=7.5, color=P.TEXT, va="center")
    lx2 = lx + 4.6
    marks = [("unit", "circle", P.TEXT, P.TEXT, 0.9, "first-look study with a primary list (recurrence unit)"),
             ("sn", "square", P.SUBTLE, P.TEXT, 0.9, "single-nucleus first look (no primary list here)"),
             ("other", "circle", P.GRID, P.SUBTLE, 0.9, "other study joined by a measured relation (grey label)"),
             ("author", "circle", "white", P.TEXT, 1.6, "hollow: bank inferred from an author (lineage lane_by)"),
             ("table", "circle", P.TEXT, P.TEXT, 0.9, "ringed: bank read from the study's own donor or dataset table")]
    for i, (key, shape, fc, ec, lw, t) in enumerate(marks):
        yy = ly - 0.2 - i * 0.24
        if shape == "square":
            m = RegularPolygon((lx2, yy), 4, radius=0.075, orientation=math.pi / 4, fc=fc, ec=ec, lw=lw)
        else:
            m = Circle((lx2, yy), RM, fc=fc, ec=ec, lw=lw)
        ax.add_patch(m)
        if key == "table":
            ax.add_patch(Circle((lx2, yy), HALO, fc="none", ec=P.SUBTLE, lw=0.9))
        ax.text(lx2 + 0.15, yy, t, fontsize=7.5, color=P.TEXT, va="center")
    lx3 = lx2 + 4.2
    ax.text(lx3, ly - 0.2, "+Bank (red, after a label): the study also draws donors from that bank", fontsize=7.5,
            color=P.TEXT, va="center")
    ax.text(lx3, ly - 0.44, "Line ending on a box: the evidence is that bank's donor table, not one study", fontsize=7.5,
            color=P.TEXT, va="center")
    ax.text(lx3, ly - 0.68, "No line between two studies: overlap not measured (unknown within a bank)", fontsize=7.5,
            color=P.TEXT, va="center")
    fig.text(0.012, 0.012, f"Source: analyses/shared_subjects_between_datasets/overlap_graph.py ({TODAY}) from study_relations.tsv, "
             "subject_relations.tsv, donor_overlap/overlap_summary.tsv, primary_lists.json, corpus/lineage.json.",
             fontsize=8, color=P.SUBTLE)
    svg_safe(ax)
    fig.savefig(OUT["svg"], format="svg")
    if pngtoo:
        fig.savefig(OUT["png"], dpi=110)
    plt.close(fig)
    return pos


# ------------------------------------------------------------------------------------------------ report
def write_outputs(D, V, edges, tested_sep, PR):
    lab = lambda k: V[k]["label"] if k in V else BOX_TITLE.get(k[4:], k)
    with open(OUT["edges"], "w") as fh:
        fh.write(f"# built {TODAY} by subject_relations/overlap_graph.py; one row per edge (unordered pair x evidence tier); the figure leaves out a stated or lineage line where the pair also has a measured one\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["a", "b", "label_a", "label_b", "tier", "relation", "k_shared", "text", "sources"])
        for x in edges:
            w.writerow([x["a"], x["b"], lab(x["a"]), lab(x["b"]), x["tier"], x["relation"], "/".join(x["ks"]), " | ".join(x["texts"]),
                        " | ".join(x["sources"])])
    with open(OUT["pairs"], "w") as fh:
        fh.write(f"# built {TODAY} by subject_relations/overlap_graph.py; every pair of drawn first-look studies (open "
                 "markers left out) with its class\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["a", "b", "label_a", "label_b", "box_a", "box_b", "banks_a", "banks_b", "class"])
        bk = lambda k: ",".join(sorted(({V[k]["bank"]} - {None}) | V[k]["extra"]))
        for a, b, c in PR:
            w.writerow([a, b, lab(a), lab(b), V[a]["box"], V[b]["box"], bk(a), bk(b), c])


def md_table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def report(D, V, edges, tested_sep, PR, inputs, sizes=None):
    lab = lambda k: V[k]["label"] if k in V else BOX_TITLE.get(k[4:], k)
    boxes = [b for row in LAYOUT for b, _ in row]
    drawn = [k for k in V if V[k]["kind"] != "other"]
    nk = collections.Counter(V[k]["kind"] for k in V)
    pc = collections.Counter(c for *_, c in PR)
    within = collections.defaultdict(collections.Counter)
    for a, b, c in PR:
        if V[a]["box"] == V[b]["box"]:
            within[V[a]["box"]][c] += 1
    tier_n = collections.Counter(x["tier"] for x in edges)
    measured = pc[PAIR_CLASSES[0]] + pc[PAIR_CLASSES[1]]
    unknown = pc[PAIR_CLASSES[5]] + pc[PAIR_CLASSES[6]]
    cross = [x for x in edges if x["tier"] in ("identifier", "fingerprint") and
             ((x["a"] in V and x["b"] in V and V[x["a"]]["box"] != V[x["b"]]["box"]) or x["a"].startswith("box:") or x["b"].startswith("box:"))]
    boxname = lambda k: BOX_TITLE[k[4:]] if k.startswith("box:") else BOX_TITLE[V[k]["box"]]
    short = lambda t: t.split(" (")[0]
    cross_boxes = collections.Counter(tuple(sorted((short(boxname(x["a"])), short(boxname(x["b"]))))) for x in cross)
    other_units = {k: v for k, v in D["units"]["studies"].items() if v.get("tissue") != "brain"}
    sr_ns = [r for r in D["sr"] if r["relation"] == SR.SEPARATE and r.get("bank_relation") == "bank not stated"]
    dr_ns = [r for r in D["dr"] if r["relation"] == SR.SEPARATE and r.get("bank_relation") == "bank not stated"]
    dr_ns_coll = [r for r in dr_ns if r["evidence"] == "collection"]
    sr_ns_coll = [r for r in sr_ns if r["evidence"] == "collection"]
    sr_db = [r for r in D["sr"] if r["relation"] == SR.SEPARATE and r["evidence"] == "collection"
             and r.get("bank_relation", "").startswith("different banks")]
    a4s = read_tsv(os.path.join(GL, "fairness_lists", "a4_summary.tsv"))
    a4_n = lambda part: next((r["n_pairs"] for r in a4s if r["part"] == part and r["region_relation"] == "all"
                              and r["donor_relation"] == SR.SEPARATE + " (bank not stated)"), "?")
    untested = {frozenset((r["dataset_a"], r["dataset_b"])) for r in dr_ns_coll}
    sl_files, sl_hits = [], 0
    for f in ("a_cross_by_pair.tsv", "b_pairs.tsv", "dataset_pairs.tsv"):
        pth = os.path.join(GL, "fairness_sample_level", "results", f)
        if os.path.exists(pth):
            sl_files.append(f)
            sl_hits += sum(frozenset((r.get("dataset_a"), r.get("dataset_b"))) in untested for r in read_tsv(pth))
    L = ["---", "type: detail", "bib: My Library", f"modified: {TODAY}", "---", "", "> [!note] Paul", ">", "",
         "# Subject and sample overlap among the first-look brain studies", ""]
    L += ["## Executive summary", "",
          f"- **Question.** Which of the first-look SZ vs control brain studies share subjects (donors) or samples, and for how many pairs is that known?",
          f"- **Data.** {nk['unit']} first-look brain studies with a primary list (the recurrence units), {nk['single-nucleus']} single-nucleus "
          f"first looks, and {nk['other']} other studies that a measured relation joins to them; relations from the project's subject "
          "relation tables, the donor-table fingerprints, and relations stated in papers.",
          f"- **Known.** Of {len(PR)} pairs of drawn first-look studies, {measured} have a measured subject overlap (identifier "
          f"{pc[PAIR_CLASSES[0]]}, fingerprint {pc[PAIR_CLASSES[1]]}), {pc[PAIR_CLASSES[2]]} share data or subjects by a stated or "
          f"curated relation, and {pc[PAIR_CLASSES[4]]} were tested with no overlap found.",
          f"- **Not known.** {pc[PAIR_CLASSES[5]]} pairs are in the same bank and untested, and {pc[PAIR_CLASSES[6]]} involve a study "
          f"whose bank is not stated; together {unknown} pairs ({unknown / len(PR):.0%}) whose overlap is unknown. "
          f"{pc[PAIR_CLASSES[7]]} pairs are in different named banks and untested; they are presumed separate on the bank alone. "
          "HBCC and LIBD are different banks and share subjects by fingerprint, so a bank difference alone does not "
          "establish separation.",
          f"- **Cross-bank.** {len(cross)} measured-overlap lines cross a box boundary (listed below), by pair of boxes: "
          + "; ".join(f"{a} with {b} {n}" for (a, b), n in sorted(cross_boxes.items(), key=lambda x: (-x[1], x[0]))) + ". "
          "A study in one box can share subjects with another box because it draws donors from several banks.",
          f"- **Label in the relation tables.** `study_relations.tsv` labels {len(sr_ns)} study pairs 'subject-separate' "
          f"when a bank is not stated, {len(sr_ns_coll)} of them with no identifier or fingerprint test; "
          f"`subject_relations.tsv` does the same for {len(dr_ns)} dataset pairs, {len(dr_ns_coll)} untested. Here the "
          "untested ones are counted as unknown; see [[#The 'bank not stated' label]].",
          ""]
    L += ["![[subject_overlap_graph.svg]]", "",
          "## How to read the figure", "",
          "- **Boxes** are brain banks. Stanley is split by sub-collection where the paper names it (Array Collection vs "
          "Neuropathology Consortium). The last box holds studies with no stated bank and studies drawing on several banks.",
          "- **Markers.** Filled circle: a first-look study with a primary list (a recurrence unit). Grey square: a "
          "single-nucleus first look, which has no primary list here. Light grey circle, grey label: another study (a "
          "reanalysis, meta-analysis or data reuse) that a measured relation joins to a drawn study. A hollow marker was "
          "placed in its bank from an author name (lineage lane_by), not a stated bank; a ringed marker, from the bank "
          "named in the study's own donor or dataset table.",
          "- **Lines** join two studies only where a relation was measured or stated. Colour gives the evidence: blue, the "
          "same subject identifiers; amber, a demographic fingerprint (diagnosis, sex, age, PMI) above chance; violet, "
          "stated in a paper or curated by hand (same data, same subjects); grey dashed, a lineage 'same dataset' link "
          "that is not marked curated (its evidence sentence is in the table below). A line ending on a box border: the evidence is that bank's donor table "
          "(the LIBD BrainSeq phenotype table), not one study.",
          "- **Red +Bank after a label**: the study also draws donors from that bank (lineage 'banks', its donor table, "
          "or a curated collection link). The CommonMind, Mount Sinai and Pittsburgh boxes note that CommonMind includes "
          "donors of the other two.",
          "- **Numbers on lines**: subjects shared (k), or the range of k when several dataset pairs of the two studies "
          "were compared.",
          "- **No line** between two studies in one box means their overlap is unknown, not that they are separate.",
          "- Lines join studies, not subjects: one line can stand for 4 shared subjects or 300. The counts are in the table "
          "of cross-bank overlaps below and in the edge file.", ""]
    rows = []
    for b in boxes:
        ks = [k for k in drawn if V[k]["box"] == b]
        if not ks:
            continue
        sd = sum(V[k]["subject_data"] for k in ks)
        e_touch = sum(1 for x in edges if x["tier"] in ("identifier", "fingerprint") and
                      ((x["a"] in V and V[x["a"]]["box"] == b) or (x["b"] in V and V[x["b"]]["box"] == b) or
                       x["a"] == "box:" + b or x["b"] == "box:" + b))
        w = within[b]
        tot = sum(w.values())
        rows.append([BOX_TITLE[b], len(ks), sd, e_touch, tot, w[PAIR_CLASSES[0]] + w[PAIR_CLASSES[1]],
                     w[PAIR_CLASSES[2]], w[PAIR_CLASSES[4]], w[PAIR_CLASSES[5]] + w[PAIR_CLASSES[6]]])
    if sizes:
        bn = json.load(open(F["bank_names"]))["collections"]
        bn_detail = "; ".join(f"{bn[c]['short']}, \"{bn[c]['detail']}\"" for c in
                              ("coll:stanley_array_collection", "coll:stanley_neuropathology_consortium") if c in bn)
        outside = [f"{V[k]['label']} ({re.search(r'New Stanley[^;,.)]*', D['node'][k].get('tissue_n', '')).group(0)})"
                   for k in drawn if V[k]["box"] in sizes and re.search(r"New Stanley", D["node"][k].get("tissue_n", ""))]
        L += ["## The fixed Stanley collections", "",
              "The Stanley Array Collection and the Neuropathology Consortium are fixed collections of donors "
              f"(`corpus/bank_names.json`: {bn_detail}). The sizes in the box subtitles come from these repo records:", ""]
        L += [f"- **{BOX_TITLE[b]}**: {n1} SZ and {n2} controls. Source: {src}." for b, (n1, n2, src) in sorted(sizes.items())]
        L += ["", "Every SZ vs control study placed in one of these boxes therefore samples the same donor pool. Which "
              "donors each study used is not recorded here, so the number shared by two studies of one box is not known; "
              "those pairs are counted as 'same bank, untested'." + (" Studies whose lineage sample description names "
              "donors from outside the collection: " + "; ".join(outside) + "." if outside else ""), ""]
    L += ["## Per bank", "",
          "Studies: drawn first-look studies (recurrence units and single-nucleus first looks) in the box. Subject-level "
          "data: studies with a per-subject table here (a dataset in the subject relation tables or a donor table). Measured "
          "overlaps: identifier or fingerprint lines touching the box, including lines to other banks. Pairs: pairs of "
          "studies within the box, by what is known.", ""]
    L += md_table(["Bank", "Studies", "With subject-level data", "Measured overlaps (lines)", "Pairs within bank",
                   "overlap measured", "stated same data or subjects", "tested, none found", "unknown"], rows)
    L += ["", "All pairs of drawn studies, by class:", ""]
    L += md_table(["Class", "Pairs"], [[c, pc[c]] for c in PAIR_CLASSES])
    L += ["", "## Cross-bank overlaps", "",
          "Measured overlaps (identifier or fingerprint) between studies placed in different boxes, or between a study and a "
          "bank's donor table.", ""]
    L += md_table(["Study", "Study or bank", "Evidence", "What was matched"],
                  [[lab(x["a"]), lab(x["b"]), x["tier"], "; ".join(x["texts"])] for x in cross])
    stated = [x for x in edges if x["tier"] in ("stated", "citation")]
    L += ["", "## Stated and citation-only relations", "",
          "Violet lines are relations stated in a paper or curated by hand; grey dashed lines are lineage 'same dataset' "
          "links not marked curated. Their evidence sentences differ in what they show (one is a row-for-row check of "
          "two lists, another a sentence about reprocessing data for comparison), so they are drawn apart and counted as "
          "their own class, not as overlap.", ""]
    L += md_table(["Study", "Study", "Tier", "Relation", "Evidence"],
                  [[lab(x["a"]), lab(x["b"]), x["tier"], x["relation"], "; ".join(t[:220] for t in x["texts"])] for x in stated])
    by_lane = collections.Counter((v.get("tissue"), D["node"][k].get("lane")) for k, v in other_units.items())
    groups_nb = [g for g in D["units"]["groups"] if len(g) > 1 and all(k in other_units for k in g)]
    L += ["", "## Studies not drawn: blood, LCL, iPSC and other tissue", "",
          f"{len(other_units)} first-look units are not from brain. Brain banks do not apply to them; the only subject "
          "relations recorded for them are shared data (primary_lists same_data_as).", ""]
    L += md_table(["Tissue", "Lineage lane", "Studies"], [[t, l, n] for (t, l), n in sorted(by_lane.items())])
    L += ["", "Pairs recorded as the same data: " + "; ".join(" = ".join(D["units"]["labels"].get(k, k) for k in g) for g in groups_nb) + "."]
    L += ["", "## The 'bank not stated' label", "",
          "`subject_relations.py` (`bank_relation`) returns 'subject-separate' when either study's bank is not stated, "
          "and the same label for different named banks. Neither is a measurement. Counts in the current tables:", "",
          f"- `study_relations.tsv`: {len(sr_ns)} study pairs 'subject-separate' with bank_relation 'bank not stated', "
          f"and {len(sr_db)} more with 'different banks' and no identifier or fingerprint test.",
          f"- `subject_relations.tsv`: {len(dr_ns)} dataset pairs 'subject-separate' with bank_relation 'bank not stated'; "
          f"{len(dr_ns_coll)} of them have no identifier or fingerprint test (evidence 'collection').",
          "",
          "Who reads the label: `fairness_lists/a4_donors.py` reports 'subject-separate (bank not stated)' as its own class "
          f"({a4_n('A published')} study pairs in part A, {a4_n('B Gemma')} result-set pairs in part B of `a4_summary.tsv`) and "
          "uses 'subject-separate (different banks)' as the reference for its bootstrap differences, so the bank-not-stated "
          "pairs are kept apart from that reference but still carry the word 'separate' in the fairness report's tables. "
          f"The sample-level pair files checked (`fairness_sample_level/results/` {', '.join(sl_files)}) contain "
          f"{sl_hits} of the {len(dr_ns_coll)} untested bank-not-stated dataset pairs; the tested ones (a fingerprint without "
          "excess) can enter their 'subject-separate' comparison sets, where a coarse fingerprint is a weak test. In the "
          "figure and the pair counts above, untested pairs with a bank not stated are unknown. HBCC and LIBD are different "
          "banks and share subjects by fingerprint; `bank_relation` handles that pair of banks as a special case.", ""]
    ns_box = [k for k in drawn if V[k]["box"] == "none"]
    nodata = collections.Counter(V[k]["box"] for k in drawn if not V[k]["subject_data"])
    L += ["## What would resolve the unknowns", "",
          f"- **Per-subject tables** (diagnosis, sex, age, PMI; identifiers where given) for the {sum(nodata.values())} drawn "
          "studies with none here: " + ", ".join(f"{BOX_TITLE[b]} {n}" for b, n in sorted(nodata.items(), key=lambda x: -x[1])) + ". "
          "With such a table the existing fingerprint method (`donor_overlap/match.py`) runs unchanged.",
          f"- **The bank** for the {len(ns_box)} studies in the last box, from their methods sections.",
          "- **Stanley sub-collection** for the studies in the 'sub-collection not stated' box, and SMRI subject "
          "identifiers for any Stanley study.",
          "- **CommonMind metadata** (Synapse syn3354385, controlled access; `donor_overlap_report.md`) for CommonMind against "
          "HBCC, Pittsburgh and Mount Sinai studies.", ""]
    L += ["## Methods", "",
          "Built by `analyses/shared_subjects_between_datasets/overlap_graph.py` (edit the script, not this file). Nodes: "
          "`recurrence.load_units()` studies with a brain primary list; lineage nodes with relation 'single-cell / "
          "single-nucleus', tissue brain, lane not 'meta'; and any other study a measured or stated relation joins to one "
          "of these. Bank: `subject_relations.Banks` (lane, Akula 2021 as HBCC, Stanley sub-collection from the list "
          "contrast or lineage parent text); else the donor table's bank when the study's table (`donor_overlap/tables/"
          "INDEX.tsv`) has exactly one; else the lineage lane when it names a bank and the node has no lineage 'banks' list. "
          "Edges: `study_relations.tsv` and `subject_relations.tsv` rows whose relation is not 'subject-separate' or 'same "
          "collection, subject overlap unknown' (dataset rows mapped to studies through `subject_relations_units.tsv`; "
          "PsychAD datasets to Lee 2025); `donor_overlap/overlap_summary.tsv` rows with shared identifiers (identifier) or "
          "passing the fingerprint rule of `subject_relations.fp_overlap` without its permutation p, which that file does "
          "not carry (fingerprint); `primary_lists.json` same_data_as and curated `corpus/lineage.json` edges with a subject "
          "kind (stated); uncurated lineage 'same dataset' edges between two drawn studies (citation). One study joined by "
          "the same text to two or more studies of one bank is drawn as one line to the bank box. Pair classes: the "
          "strongest tier of any edge; else 'tested, no overlap found' when an identifier or fingerprint comparison of the "
          "two exists without overlap (including donor tables below the fingerprint rule, which for coarse ages is weak); "
          "else by bank (same bank or CommonMind with a source bank: unknown; a bank not stated: unknown; different named "
          "banks, or different Stanley sub-collections: presumed separate).", "",
          "Files: `subject_overlap_graph.svg` (and `.png`), `subject_overlap_edges.tsv` (every drawn edge with its source "
          "rows), `subject_overlap_pairs.tsv` (every pair with its class), `subject_overlap_graph_meta.json` (input sha256).",
          "", "Inputs:", ""]
    L += [f"- `{os.path.relpath(p, GL)}` sha256 {h[:12]}" for p, h in inputs.items()]
    open(OUT["report"], "w").write("\n".join(L) + "\n")


def main(png=True):
    D = load()
    V, edges, tested_sep = build(D)
    PR = pairs(V, edges, tested_sep)
    sizes = stanley_sizes(D)
    draw(V, edges, pngtoo=png, sizes=sizes)
    inputs = {p: sha(p) for p in F.values()}
    write_outputs(D, V, edges, tested_sep, PR)
    report(D, V, edges, tested_sep, PR, inputs, sizes)
    meta = {"artifact": os.path.basename(OUT["svg"]), "sha256": sha(OUT["svg"]), "script": "analyses/shared_subjects_between_datasets/overlap_graph.py",
            "generated": TODAY, "inputs": {os.path.relpath(p, REPO): h for p, h in inputs.items()},
            "outputs": {os.path.basename(p): sha(p) for k, p in OUT.items() if k not in ("meta",) and os.path.exists(p)},
            "counts": {"nodes": collections.Counter(v["kind"] for v in V.values()),
                       "edges": collections.Counter(x["tier"] for x in edges),
                       "pairs": collections.Counter(c for *_, c in PR)}}
    json.dump(meta, open(OUT["meta"], "w"), indent=1, sort_keys=True, default=dict)
    return V, edges, PR


if __name__ == "__main__":
    main()
