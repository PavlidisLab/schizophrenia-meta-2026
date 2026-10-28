"""P10 META15 inputs: adaptation of analyses/mab_numbers_2026-10-12/task_b_meta_pooled.py to the core-80 base (80 studies, 72 contributing lists).
Changes: the 63 'published, brain' studies are replaced by the 80 base studies (base80_studies.tsv); same-data groups and lanes from the B80
recurrence_units.tsv (72 non-empty lists); primary_lists.json / universes.json are the B80 versions; meta_link_evidence.tsv is the 2026-10-12 file
(copied unchanged to in_meta_link_evidence.tsv: its confirmed links are re-checked against the 80, see link_check table). Output only in the private folder."""
import collections, csv, hashlib, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = HERE + "/out"; os.makedirs(OUT, exist_ok=True)
REPO = "sztr-rebuild-infra/repo"; AN = REPO + "/analyses"
WS = "base78_2026-10-16/B80/ws"; B80 = "base78_2026-10-16/B80/out"
MAIL = "agent-mail/sz-transcriptome"
def read_tsv(path): return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()
H3 = os.path.join(AN, "h3")
UNITS, INV, DS = (os.path.join(H3, f) for f in ("h3_units.tsv", "h3_inventory.tsv", "h3_datasets.tsv"))
RU = B80 + "/recurrence_units.tsv"; BASE80 = MAIL + "/recurrence_on_80_studies/base80_studies.tsv"; BASE81 = MAIL + "/base81_2026-10-14/base81.tsv"
LIN = os.path.join(REPO, "corpus", "lineage.json")
PL = WS + "/analyses/primary_lists.json"
UNI = WS + "/analyses/universes.json"
EXC = os.path.join(REPO, "corpus", "corpus_exclusions.json")
GEO = os.path.join(REPO, "corpus", "geo_links.json")
STAN = os.path.join(REPO, "corpus", "stanley_datasets.json")
EVID = os.path.join(HERE, "in_meta_link_evidence.tsv")
SRC = [UNITS, INV, DS, RU, BASE80, BASE81, LIN, PL, UNI, EXC, GEO, STAN, EVID]

units = [r for r in read_tsv(UNITS) if r["tissue_set"] == "brain"]
inv = {(r["h3_unit"]): r for r in read_tsv(INV) if r["h3_tissue_set"] == "brain" and r["in_pairs"] == "yes"
       and r["h3_class"].startswith("cross-cohort")}
dsrows = {r["dataset"]: r for r in read_tsv(DS)}
ru = read_tsv(RU)
b80 = read_tsv(BASE80); assert len(b80) == 80
empty80 = {r["pmid"] for r in b80 if r["empty_list_dropped"] == "True"}; assert len(empty80) == 8
rep72 = {r["representative"]: r for r in ru if r["scenario"] == "published, brain" and r["modelled"] == "yes"}
assert len(rep72) == 72 and len(units) == 15, (len(rep72), len(units))
reps = {r["pmid"]: dict(label=r["label"], lane=(rep72.get(r["pmid"]) or {}).get("lane", ""), cluster=(rep72.get(r["pmid"]) or {}).get("cluster", "")) for r in b80}
assert set(rep72) <= set(reps), set(rep72) - set(reps)
member_of = {}
for k, r in rep72.items():
    for m in re.findall(r"\((\d+),", r["group_members"]):
        member_of.setdefault(m, k)
left = {}
notmod = {r["representative"]: r for r in ru if r["scenario"] == "published, brain" and r["modelled"] != "yes"}
lin = json.load(open(LIN))
nodes = {n["id"]: n for n in lin["nodes"]}
pl = json.load(open(PL))["studies"]
uni = json.load(open(UNI))["lists"]
exc = json.load(open(EXC))
geo = json.load(open(GEO))
stan = {d["short_name"]: d for d in json.load(open(STAN))["datasets"]}
lab = lambda n: nodes[n]["label"] if n in nodes else (nodes["sc:" + n]["label"] if "sc:" + n in nodes else n)

# SMRIDB key -> primary representative (describing paper in a primary group; Altar 2005 for A/B/C as prior_meta_use.py does)
smri_primary = {}
for k, d in stan.items():
    for p in d.get("describing_papers", []):
        if p.get("in_corpus") and p.get("pmid") in reps:
            smri_primary[k] = p["pmid"]
for k in ("stanley_altarA", "stanley_altarB", "stanley_altarC"):
    smri_primary[k] = "16038679"  # Altar 2005, the Array Collection's own paper (prior_meta_use.py)

lane_cmc = sorted(r["label"] for r in reps.values() if r["lane"] == "cmc")


# dataset keys whose h3.py label names a corpus study or consortium, but whose cohort_evidence does not cite a lineage node
NAME_MAP = {"ruzicka2024": ("38781388", "dataset label 'Ruzicka 2024 snRNA-seq' (h3.py DATASETS) names the study"),
            "ling2024": ("38448582", "dataset label 'Ling 2024 snRNA-seq' (h3.py DATASETS) names the study"),
            "frohlich2024": ("39227716", "dataset label 'Frohlich 2024 snRNA-seq' (h3.py DATASETS) names the study"),
            "cmc_hbcc": ("coll:cmc", "dataset label 'CMC-HBCC' (h3.py DATASETS)"),
            "libd_jaffe": ("coll:libd_brainseq", "dataset label 'LIBD DLPFC RNA-seq' (h3.py DATASETS)")}


evid = {r["dataset"]: r for r in read_tsv(EVID)}
assert all(r["status"] in ("confirmed", "unresolved") for r in evid.values())
for r in evid.values():  # every evidence file named must exist
    for f in r["evidence_files"].split("; "):
        assert os.path.exists(os.path.join(REPO, f)), (r["dataset"], f)
    assert r["dataset"] in dsrows, r["dataset"]
    if r["status"] == "confirmed":
        assert r["link_to_node"] in nodes or "sc:" + r["link_to_node"] in nodes, r["link_to_node"]


def node_of(ds):
    ev = dsrows[ds]["cohort_evidence"]
    if ds in evid and evid[ds]["status"] == "confirmed":
        return evid[ds]["link_to_node"], "meta_link_evidence.tsv: " + evid[ds]["basis"]
    if ds in NAME_MAP:
        return NAME_MAP[ds][0], NAME_MAP[ds][1]
    m = re.search(r"lineage node (\S+) \(", ev)
    if m:
        return m.group(1), "h3_datasets.tsv cohort_evidence"
    if ds.startswith("GSE"):
        c = [p for p in geo["gse_pmids"].get(ds, []) if p in nodes]
        if len(c) == 1:
            return c[0], "corpus/geo_links.json gse_pmids"
    return None, ev


def classify(ds):
    """-> (is_primary, primary_unit, category, resolved, reason_recorded, evidence, settle)"""
    nid, how = node_of(ds)
    if nid is None:
        extra_note = ""
        if ds == "cohen2012":
            extra_note = "; lineage has a node labelled 'Cohen 2012' (PMID 23062752, lane harvard) but the match to Hess 2016's citation is not verified"
        return (False, "", "not recorded", "no lineage node",
                f"dataset has no corpus node ({how}){extra_note}; whether its source paper reports its own SZ list is not recorded",
                "h3.py DATASETS entry", "the dataset's source paper / BioProject record (accession or the citing meta-analysis's methods table) and a check of its supplementary tables")
    if nid.startswith("sc:"):
        nid = nid[3:]
    if nid.startswith("smri:"):
        key = nid[5:]
        if key in smri_primary:
            return (True, lab(smri_primary[key]), "", f"SMRIDB {key}", "", "", "")
        d = stan.get(key, {})
        dp = d.get("describing_papers", [])
        if dp:
            p = dp[0]
            return (False, "", "other", f"Stanley Online Genomics Database study {key}",
                    f"Stanley database study, not a corpus study; recorded describing paper {p['first_author']} {p['year']} (PMID {p['pmid']}) is not in the corpus (in_corpus {p['in_corpus']}) and is a bipolar disorder study",
                    "corpus/stanley_datasets.json describing_papers", "")
        return (False, "", "not recorded", f"Stanley Online Genomics Database study {key}",
                "Stanley database study with no describing paper recorded; lineage node is 'Cited source (not in corpus)'; whether a published SZ list exists is not recorded",
                "corpus/stanley_datasets.json (describing_papers empty); lineage.json node",
                "the Stanley Online Genomics Database study page and PubMed for the study's investigator")
    if nid.startswith("coll:"):
        return (False, "", "other", lab(nid) + " (consortium)",
                f"consortium data set, not a single study; lineage node 'Cited source (not in corpus)'; the primary studies on this data are the lane-cmc studies ({len(lane_cmc)} of 80: " + "; ".join(lane_cmc) + ")",
                "lineage.json node; recurrence_units.tsv lane", "")
    if nid in reps:
        linked = ds in evid and evid[ds]["status"] == "confirmed"
        return (True, lab(nid), "", lab(nid) + (f" (linked by meta_link_evidence.tsv: {evid[ds]['basis']})" if linked else ""), "",
                "meta_link_evidence.tsv" if linked else "", "")
    if nid in member_of:
        return (False, lab(member_of[nid]), "reanalysis only", f"{lab(nid)} (same_data_as group of {lab(member_of[nid])})",
                f"same data as primary study {lab(member_of[nid])} (same_data_as group)", "recurrence_units.tsv group_members", "")
    L = nodes.get(nid) or nodes.get("sc:" + nid, {})
    label = L.get("label", nid)
    pe = pl.get(nid)
    gl = (L.get("genelist") or {})
    gl_none = bool(gl.get("lists")) and all(x.get("status") == "no list published" for x in gl["lists"])
    ex = exc.get(nid)
    if ex:
        return (False, "", "failed eligibility", label, f"corpus_exclusions.json: {ex['reason']} ({ex['date']})", "corpus/corpus_exclusions.json", "")
    if nid in left:
        t = left[nid]
        cat = "reanalysis only" if "reanalysis" in t else "other"
        return (False, "", cat, label, f"recurrence_units.tsv (left out): {t}; lineage category {L.get('category')}, relation '{L.get('relation')}'",
                "analyses/recurrence_units.tsv", "")
    if nid in notmod:
        ul = [v for k, v in uni.items() if v.get("study") == nid]
        why = ul[0]["reason"] + "; " + ul[0]["evidence"] if ul else "universe entry not found"
        return (False, "", "other", label, f"'published, brain' study left out of the modelled units (universe type 'exclude'): {why}",
                "recurrence_units.tsv; universes.json", "")
    if pe and pe["provenance"] != "published":
        return (False, "", "no published SZ list of its own", label,
                f"primary_lists.json provenance '{pe['provenance']}': {pe['reason']}", "analyses/primary_lists.json", "")
    if gl_none:
        return (False, "", "no published SZ list of its own", label, "lineage genelist status 'no list published': " + gl["lists"][0]["reason"],
                "corpus/lineage.json genelist", "")
    ext = [x for x in gl.get("lists", []) if x.get("status") in ("extracted", "gemma", "reproduced")]
    if ext and not pe:
        return (False, "", "no primary_lists.json entry", label,
                f"corpus study (lineage category {L.get('category')}, relation '{L.get('relation')}') with {len(ext)} extracted SZ list(s) ({', '.join(x['list_id'] for x in ext)}) and no entry in primary_lists.json, so it was never a candidate for the 80; no exclusion recorded",
                "corpus/lineage.json genelist; analyses/primary_lists.json (absent); corpus/corpus_exclusions.json (absent)", "")
    return (False, "", "not recorded", label, "no decision-table entry found", "", "recurrence_units.tsv / primary_lists.json / corpus_exclusions.json")



# ================= core-80 tail: coverage classes, counts, subjects, overlap =================
b81 = {r["pmid"]: r for r in read_tsv(BASE81)}
def num(x):
    try: return int(float(x))
    except Exception: return None
def base_pmid(ds):
    nid, _ = node_of(ds)
    if nid is None: return None, None
    if nid.startswith("sc:"): nid = nid[3:]
    if nid.startswith("smri:"):
        k = smri_primary.get(nid[5:]); return (k, "direct") if k in reps else (None, None)
    if nid.startswith("coll:"): return None, None
    if nid in reps: return nid, "direct"
    if nid in member_of and member_of[nid] in reps: return member_of[nid], "same_data_as"
    return None, None
def outside_corpus(ds):
    nid, _ = node_of(ds)
    if nid is None or nid.startswith("smri:") or nid.startswith("coll:"): return True
    n = nodes.get(nid) or nodes.get("sc:" + nid)
    return n is None or bool(n.get("external")) or n.get("source") == "Cited source (not in corpus)"
def subjects(ds, bp):
    """(cases, controls, source) known, else (None, None, 'unknown'); no estimates"""
    nid, _ = node_of(ds)
    if bp and bp in b81 and num(b81[bp]["n_SZ"]) is not None:
        return num(b81[bp]["n_SZ"]), num(b81[bp]["n_control"]), "base81.tsv (authors' recorded counts, before QC) for study " + b81[bp]["study"]
    if nid:
        n = nodes.get(nid) or nodes.get("sc:" + nid)
        if n and num(n.get("cases")) is not None:
            return num(n["cases"]), num(n.get("controls")), "lineage.json node cases/controls"
    return None, None, "unknown"
CLASS = {"i": "(i) among the 80 base studies", "ii": "(ii) in the corpus, not among the 80", "iii": "(iii) outside the corpus"}
rows = []
for u in units:
    r = inv[u["unit"]]
    for ds in r["input_datasets"].split(";"):
        isp, punit, cat, res, why, ev, settle = classify(ds)
        if ds in evid and evid[ds]["status"] == "unresolved":
            why += f"; link to {lab(evid[ds]['link_to_node'])} NOT applied (unresolved): {evid[ds]['evidence']}"
            ev += ("; " if ev else "") + "in_meta_link_evidence.tsv"
        bp, how = base_pmid(ds)
        if bp: cl = "i"
        elif outside_corpus(ds): cl = "iii"
        else: cl = "ii"
        if bp:
            cat2, why2 = "", ""
        else:
            cat2, why2 = cat, why
        c, k, src = subjects(ds, bp)
        n80 = reps[bp]["label"] if bp else ""
        rows.append({"meta_analysis": u["label"], "meta_unit": u["unit"], "dataset": ds, "dataset_label": dsrows[ds]["label"],
                     "coverage_class": CLASS[cl], "base_study_pmid": bp or "", "base_study": n80, "link_type": how or "",
                     "base_study_has_list_among_72": ("no (empty derived list, dropped)" if bp in empty80 else "yes") if bp else "",
                     "cohort_cluster": reps[bp]["cluster"] if bp else "",
                     "resolved_to": res, "reason_category_if_not_in_80": cat2, "reason_recorded": why2,
                     "cases": "" if c is None else c, "controls": "" if k is None else k, "subjects_source": src,
                     "evidence_source": ev})
F = list(rows[0])
TODAY = "2026-10-16"
def write(name, rows, claim, extra=()):
    head = [f"# figure/claim: {claim}", "# script: task_b_core80.py (assembled from task_b_core80_head.py + task_b_core80_tail.py; adapted from figure_data/mab_numbers_2026-10-12/task_b_meta_pooled.py)",
            "# sources (sha256): " + "; ".join(f"{os.path.relpath(p, '$HOME')} {sha(p)}" for p in SRC), f"# date: {TODAY}",
            "# parameters: no random draws; base = 80 studies / 72 contributing lists (8 empty derived lists dropped); 15 brain meta-analyses = brain units of h3_units.tsv"] + [f"# {e}" for e in extra]
    with open(OUT + "/" + name, "w", newline="") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
note = ["Table 5 (data needed) register is not on this host: class (ii) is split by the reason recorded in the repository (reason_category_if_not_in_80), not by Table 5 membership. Ketchesin 2023 and Pai 2019 (D33 rows) are among the 80 (derived lists; Pai's is empty and dropped) and are reached by no meta-analysis dataset",
        "subject counts: only where recorded (base81.tsv n_SZ / n_control for studies among the 80, else lineage node cases/controls); otherwise blank = unknown; none estimated"]
write("meta15_input_datasets.tsv", rows, "META15 part 4 / B: each dataset pooled by each of the 15 brain meta-analyses, with its coverage class against the 80 base studies", note)

# per meta-analysis
per = []
ucl = {}  # unit -> set clusters of its input cohorts
for u in units:
    rr = [x for x in rows if x["meta_unit"] == u["unit"]]
    coh = [c for c in inv[u["unit"]]["input_cohorts"].split(";") if c]
    ucl[u["unit"]] = set(coh)
    studies = {x["base_study_pmid"] for x in rr if x["base_study_pmid"]}
    n_dir = len({x["base_study_pmid"] for x in rr if x["link_type"] == "direct"})
    cluster80 = {reps[p]["cluster"] for p in reps if reps[p]["cluster"]}
    cc = collections.Counter(x["coverage_class"] for x in rr)
    def sm(cls, key):
        s = [x for x in rr if x["coverage_class"] == CLASS[cls]]
        kn = [x for x in s if x[key] != ""]
        return sum(int(x[key]) for x in kn), len(s) - len(kn)
    ci, ui = sm("i", "cases"); cn, _ = sm("i", "controls")
    c2, u2 = sm("ii", "cases"); k2, _ = sm("ii", "controls")
    c3, u3 = sm("iii", "cases"); k3, _ = sm("iii", "controls")
    reas = collections.Counter(x["reason_category_if_not_in_80"] for x in rr if x["coverage_class"] == CLASS["ii"])
    N = len(rr)
    per.append({"meta_analysis": u["label"], "meta_unit": u["unit"], "n_input_datasets": N,
                "i_datasets_among_80": cc[CLASS["i"]], "i_share": f"{cc[CLASS['i']] / N:.2f}", "i_distinct_base_studies_pooled": len(studies),
                "i_of_which_direct": n_dir, "i_datasets_reached_via_same_data_group": sum(1 for x in rr if x["link_type"] == "same_data_as"),
                "i_studies_contributing_one_of_72_lists": len(studies - empty80), "i_studies_with_empty_list_dropped": len(studies & empty80),
                "ii_datasets_in_corpus_not_80": cc[CLASS["ii"]], "ii_share": f"{cc[CLASS['ii']] / N:.2f}", "ii_by_reason": "; ".join(f"{k}: {v}" for k, v in sorted(reas.items())),
                "iii_datasets_outside_corpus": cc[CLASS["iii"]], "iii_share": f"{cc[CLASS['iii']] / N:.2f}",
                "i_cases_known": ci, "i_controls_known": cn, "i_datasets_subjects_unknown": ui,
                "ii_cases_known": c2, "ii_controls_known": k2, "ii_datasets_subjects_unknown": u2,
                "iii_cases_known": c3, "iii_controls_known": k3, "iii_datasets_subjects_unknown": u3,
                "n_input_cohort_clusters": len(coh), "cohort_clusters_shared_with_80": len(set(coh) & cluster80),
                "cohort_clusters": ";".join(coh)})
for p_ in per:
    others = set().union(*[ucl[x] for x in ucl if x != p_["meta_unit"]])
    p_["cohort_clusters_shared_with_another_meta"] = len(ucl[p_["meta_unit"]] & others)
write("meta15_input_coverage.tsv", per, "META15 part 4 / B: per meta-analysis, how many of its input datasets are among the 80 base studies, in the corpus but not among the 80, or outside the corpus; studies pooled; subject counts where recorded; cohort clusters", note +
      ["the 80 base studies = base80_studies.tsv; 72 contribute a list (8 derived lists are empty and dropped); cohort clusters are the input_cohorts of h3_inventory.tsv and the 'cluster' column of the B80 recurrence_units.tsv (empty-list studies have no cluster there)"])

# pooled over metas / distinct datasets
dist = {}
for x in rows: dist.setdefault(x["dataset"], x)
dc = collections.Counter(x["coverage_class"] for x in dist.values()); pc = collections.Counter(x["coverage_class"] for x in rows)
bs = {x["base_study_pmid"] for x in rows if x["base_study_pmid"]}
reasd = collections.Counter(x["reason_category_if_not_in_80"] for x in dist.values() if x["coverage_class"] == CLASS["ii"])
summ = [{"quantity": "meta-analyses", "value": len(units)}, {"quantity": "dataset slots (meta-analysis x dataset pairs)", "value": len(rows)},
        {"quantity": "distinct pooled datasets", "value": len(dist)}]
for k_, v in CLASS.items():
    summ.append({"quantity": f"slots {v}", "value": pc[v]}); summ.append({"quantity": f"distinct datasets {v}", "value": dc[v]})
summ += [{"quantity": "distinct base studies reached", "value": len(bs)}, {"quantity": "of which with an empty list (no contribution among the 72)", "value": len(bs & empty80)},
         {"quantity": "of which contributing one of the 72 lists", "value": len(bs - empty80)}]
for k_, v in sorted(reasd.items()): summ.append({"quantity": f"distinct datasets in (ii), reason '{k_}'", "value": v})
write("meta15_input_counts.tsv", summ, "META15 part 4 / B: input datasets of the 15 meta-analyses against the 80 base studies (counts)", note)

# link check of the confirmed evidence rows against the 80
chk = []
for d_, e in evid.items():
    nid = e["link_to_node"]; in80 = nid in reps or (nid in member_of and member_of[nid] in reps)
    chk.append({"dataset": d_, "status_2026_10_12": e["status"], "link_to_node": nid, "node_label": lab(nid), "node_among_80_base": "yes" if in80 else "no",
                "node_in_80_has_list_among_72": ("yes" if nid not in empty80 else "no (empty)") if in80 else "", "basis": e["basis"]})
write("meta15_link_check_vs_80.tsv", chk, "META15 B: the 2026-10-12 meta_link_evidence.tsv linkings re-checked against the 80 base studies (status unchanged; only membership in the 80 is recomputed)", note)

# pairwise overlap counts between metas (shared base studies)
pw = []
st_of = {u["unit"]: {x["base_study_pmid"] for x in rows if x["meta_unit"] == u["unit"] and x["base_study_pmid"]} for u in units}
ds_of = {u["unit"]: {x["dataset"] for x in rows if x["meta_unit"] == u["unit"]} for u in units}
for a in units:
    for b in units:
        if a["unit"] < b["unit"]:
            pw.append({"meta_a": a["label"], "meta_b": b["label"], "shared_datasets": len(ds_of[a["unit"]] & ds_of[b["unit"]]),
                       "shared_base_studies": len(st_of[a["unit"]] & st_of[b["unit"]]), "shared_cohort_clusters": len(ucl[a["unit"]] & ucl[b["unit"]])})
write("meta15_pairwise_shared_inputs.tsv", pw, "META15 part 4: counts of datasets, base studies and cohort clusters shared by each pair of the 15 meta-analyses", note)
for s_ in summ: print(s_["quantity"], "=", s_["value"])
for p_ in per: print(p_["meta_analysis"], p_["n_input_datasets"], p_["i_datasets_among_80"], p_["i_distinct_base_studies_pooled"], p_["i_studies_contributing_one_of_72_lists"], p_["ii_datasets_in_corpus_not_80"], p_["iii_datasets_outside_corpus"], p_["cohort_clusters_shared_with_80"], "/", p_["n_input_cohort_clusters"])
