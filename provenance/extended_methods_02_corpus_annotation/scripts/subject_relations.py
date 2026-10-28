#!/usr/bin/env python3
"""Canonical subject relations between the sample-level datasets of the project (tracker, Vocabulary, 2026-09-25).

    opt/anaconda3/bin/python3 subject_relations/subject_relations.py      (from analyses/)

A *subject* is a person (donor); a *sample* is one tissue piece or RNA extract from one subject and one region. Every
analysis that needs the relation of two datasets (or two studies' lists) reads it from the outputs here; no other
script classifies subject relations.

Datasets (units):
  registry   the 28 datasets of controls_age_sex/data/datasets.tsv (used by controls_age_sex and fairness_sample_level);
             samples = the GSMs those analyses use (controls_age_sex/results/sex_check_samples.tsv), fields from the
             registry covariate files
  LIBD       LIBD_DLPFC, LIBD_HIPPO: adult BrainSeq Phase 2 samples of libd_brainseq/analysis/results/samples.tsv
  PsychAD    PsychAD_HBCC, PsychAD_MSSM: psychad_cellxgene/donors.tsv, donors kept by the PsychAD pseudobulk rule
             (bank H or M; SZ, or control: dx_gemma CTRL / CTRL_SOMATIC or paper_dx CTRL); age, sex, PMI from Supp. Table 1
  Gemma sets result sets of fairness_lists a4 part B that are not registry datasets: two composites of registry datasets
             (GSE12679_all, GSE107638_all) and GSE12654, GSE12649, GSE62191, GSE78936 (3 regions), GSE93987 (L3, L5),
             GSE87610 (L3, L5) from their GEO series matrices (SZ and control samples; no age/sex/PMI except GSE62191 age)
Only SZ and control samples are counted.

Subjects: samples are collapsed to subjects by a subject identifier where the source carries one (ID_RULES below; each
rule was checked: every identifier shared across regions of one series carries identical diagnosis, sex, age, PMI and pH).
Without an identifier each sample counts as one subject.

Per dataset pair (subject_relations.tsv), evidence in this order:
  1. same samples     the smaller dataset's GSMs are all in the other (composite vs its component)
  2. identifier       both datasets carry identifiers from one identifier space (same series, or GSE17612/GSE21935 whose
                      sample titles share the CCHPC subject codes): k = identifiers in common
  3. fingerprint      donor_overlap/match.py (imported): same dx class and sex, age and PMI intervals overlap at each
                      source's stated precision (intervals from donor_overlap/build_tables.py age_interval/pmi_interval,
                      read from that file's source); k = maximum one-to-one assignment; chance = B's ages and PMIs permuted
                      within dx x sex, 200 times (mean, 95th percentile, p). Level 'exact' when both datasets record PMI,
                      else 'no_pmi' (dx + sex + age, weak). Counts as overlap only if k >= 3, k > chance 95th
                      percentile, permutation p <= 0.01, k - chance mean >= 3 and k >= 2 x chance mean (fp_overlap). The
                      last three guard against the ~560 pairs tested: the one pair with identifiers showing no shared subject
                      and a fingerprint (PsychAD HBCC vs MSSM, 0 shared donor_id) has 14 fingerprint matches against a
                      chance mean of 6.9, p 0.01.
  4. stated           the submitter states the samples come from the same subjects (STATED)
  5. collection       presumption from the brain bank (bank per study: lineage lane, overrides as fairness_lists a4 had
                      them; Stanley split by sub-collection from the list contrast or lineage parent text)
Relation (column relation; Vocabulary labels):
  same samples | same subjects, same region | same subjects, different region | partial subject overlap |
  same collection, subject overlap unknown | subject overlap unknown, bank not stated | subject-separate
  'same subjects' when k >= 90% of the smaller dataset's subjects (SAME_FRAC), else 'partial subject overlap'; k, n_a and
  n_b are always given (relation_text states 'k of n'). Without identifier, fingerprint or stated overlap: same bank ->
  'same collection, subject overlap unknown'; different banks -> 'subject-separate'; a bank not stated ->
  'subject overlap unknown, bank not stated', or 'subject-separate' when a fingerprint test found no excess.
The fingerprint columns are filled for every pair where both datasets have diagnosis, sex and age, also where an
identifier decided the relation (a check of the fingerprint against identifiers).

Outputs (each with '#' provenance lines):
  subject_relations.tsv        one row per unordered dataset pair
  subject_relations_units.tsv  one row per dataset (n samples, n subjects, identifier space, bank, region, platform) and
                               per published-list study used by fairness_lists (bank, Stanley sub-collection, dataset)
  study_relations.tsv          one row per pair of studies (fairness_lists a4 part A, published lists) and per extra
                               list pair: the dataset-pair relation where both map to datasets here (STUDY_DATASET,
                               LIST_DATASET), else documented facts (STUDY_FACTS), else the collection presumption
"""
import ast
import collections
import csv
import datetime
import gzip
import hashlib
import itertools
import json
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GL = os.path.dirname(HERE)
TODAY = datetime.date.today().isoformat()
sys.path.insert(0, os.path.join(GL, "donor_overlap"))
import match as M  # noqa: E402  donor_overlap fingerprint method (level, matched_sets, max_matching, permuted, q95)

SEED = 20260925
SAME_FRAC = 0.9
MIN_K, MAX_P, MIN_EXCESS, MIN_RATIO = 3, 0.01, 3, 2
CAS = os.path.join(GL, "controls_age_sex")
REG = os.path.join(CAS, "data", "datasets.tsv")
USED = os.path.join(CAS, "results", "sex_check_samples.tsv")
LIBD = os.path.join(GL, "libd_brainseq", "analysis", "results", "samples.tsv")
PSYCHAD = os.path.join(GL, "psychad_cellxgene", "donors.tsv")
GEO = {"GSE12654": "confounds/data/geo/GSE12654_series_matrix.txt.gz",
       "GSE12649": "confounds/data/geo/GSE12649_series_matrix.txt.gz",
       "GSE62191": "confounds/data/geo/GSE62191_series_matrix.txt.gz",
       "GSE93987": "confounds/data/geo/GSE93987_series_matrix.txt.gz",
       "GSE87610": "confounds/data/geo/GSE87610_series_matrix.txt.gz",
       "GSE78936": "subject_relations/raw/geo/GSE78936_series_matrix.txt.gz",
       "GSE53987": "confounds/data/geo/GSE53987_series_matrix.txt.gz",
       "GSE138082": "subject_relations/raw/geo/GSE138082_series_matrix.txt.gz"}
INPUTS = []


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def rel(p):
    return os.path.relpath(p, GL)


def use(p):
    INPUTS.append((rel(p), sha(p)))
    return p


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


# interval rules of donor_overlap/build_tables.py, taken from its source (the module writes tables when imported)
_bt = os.path.join(GL, "donor_overlap", "build_tables.py")
_ns = {"re": re}
for node in ast.parse(open(_bt).read()).body:
    if isinstance(node, ast.FunctionDef) and node.name in ("blank", "ndec", "age_interval", "pmi_interval"):
        exec(compile(ast.Module([node], []), _bt, "exec"), _ns)
age_interval, pmi_interval, blank = _ns["age_interval"], _ns["pmi_interval"], _ns["blank"]

# ---------------------------------------------------------------- dataset facts
# region (anatomical; different cell populations of one region are the same region) and compartment
REGION = {
    "GSE21138": "DLPFC", "GSE35977": "parietal cortex", "GSE35974": "cerebellum", "GSE53987_BA46": "DLPFC",
    "GSE53987_hippocampus": "hippocampus", "GSE53987_striatum": "associative striatum", "GSE17612": "anterior PFC (BA10)",
    "GSE21935": "superior temporal cortex (BA22)", "GSE92538_GPL570": "DLPFC", "GSE92538_GPL96": "DLPFC",
    "GSE80655_DLPFC": "DLPFC", "GSE80655_AnCg": "anterior cingulate", "GSE80655_nAcc": "nucleus accumbens",
    "GSE208338": "DLPFC", "GSE202537_NAc": "nucleus accumbens", "GSE202537_caudate": "caudate",
    "GSE202537_putamen": "putamen", "GSE224683": "DLPFC", "GSE26927_temporal": "temporal lobe",
    "GSE138082_CA1": "hippocampus CA1", "GSE138082_CA3": "hippocampus CA3", "GSE138082_DG": "dentate gyrus",
    "GSE37981": "superior temporal gyrus", "GSE12679_neuron": "DLPFC", "GSE12679_endothelial": "DLPFC",
    "GSE112523": "frontal cortex (not specified)", "GSE107638_NeuN": "DLPFC", "GSE107638_OLIG2": "DLPFC",
    "LIBD_DLPFC": "DLPFC", "LIBD_HIPPO": "hippocampus", "PsychAD_HBCC": "DLPFC", "PsychAD_MSSM": "DLPFC",
    "GSE12679_all": "DLPFC", "GSE107638_all": "DLPFC", "GSE12654": "anterior PFC (BA10)", "GSE12649": "DLPFC",
    "GSE62191": "frontal cortex (not specified)", "GSE78936_BA9": "DLPFC", "GSE78936_BA11": "orbitofrontal cortex (BA11)",
    "GSE78936_BA24": "anterior cingulate", "GSE93987_L3": "DLPFC", "GSE93987_L5": "DLPFC", "GSE87610_L3": "DLPFC",
    "GSE87610_L5": "DLPFC"}
# study (PMID) whose data the dataset is, for the bank; bank overrides where the lineage lane is not a bank
STUDY = {"GSE21138": "18778695", "GSE35977": "23147385", "GSE35974": "23147385", "GSE53987": "31123247",
         "GSE17612": "19255580", "GSE21935": "21538462", "GSE92538": "30016334", "GSE80655": "28754123",
         "GSE208338": "40253403", "GSE202537": "36302706", "GSE224683": "36932057", "GSE26927": "25119539",
         "GSE138082": "32152472", "GSE37981": "24702465", "GSE12679": "19088852", "GSE112523": "31053723",
         "GSE107638": "31288836", "BrainSeq Phase 2": "31174959", "PsychAD": "", "GSE12654": "14743183",
         "GSE12649": "15563509", "GSE62191": "25981335", "GSE78936": "27661005", "GSE93987": "25560755",
         "GSE87610": "28476208"}
# PsychAD by source bank; GSE92538 (lineage lane brain_other) is the Pritzker consortium DLPFC series (registry cluster
# 'pritzker', controls_age_sex/scripts/build_covariates.py)
BANK_DATASET = {"PsychAD_HBCC": "hbcc", "PsychAD_MSSM": "mssm", "GSE92538_GPL570": "pritzker", "GSE92538_GPL96": "pritzker"}
# submitter statements of same subjects
STATED = {frozenset(("GSE78936_BA9", "GSE78936_BA24")):
          "GEO GSE78936 overall design: '19 BA24 and 19 BA9 samples from the same subjects'"}
# identifier rules per series: (id space, rule text); checked in check_ids()
ID_RULES = {
    "GSE92538": ("GSE92538 'subject id'", "GEO/Gemma characteristic 'subject id' (registry covariate 'subject')"),
    "GSE107638": ("GSE107638 'subject id'", "GEO/Gemma characteristic 'subject id' (registry covariate 'subject')"),
    "GSE202537": ("GSE202537 'patients id'", "GEO/Gemma characteristic 'patients id' (registry covariate 'subject')"),
    "GSE53987": ("GSE53987 title index", "GEO sample title '<dx>_<region>_<n>': <dx>_<n>"),
    "GSE80655": ("GSE80655 title prefix", "GEO sample title prefix 'X<n>' (e.g. X3079_nAcc_B_SL7392)"),
    "GSE17612": ("CCHPC code", "GEO sample title prefix (C014, S014), shared by GSE17612 and GSE21935"),
    "GSE21935": ("CCHPC code", "GEO sample title prefix (C014, S014), shared by GSE17612 and GSE21935"),
    "GSE35977": ("GSE35977 sample code", "Gemma sample name prefix (CM123); replicate arrays share it"),
    "GSE35974": ("GSE35974 sample code", "Gemma sample name prefix (AC49, AC101_rep1); replicate arrays share it"),
    "GSE21138": ("GSE21138 sample name", "Gemma sample name (Scz-28)"),
    "GSE208338": ("GSE208338 'sample'", "GEO characteristic 'sample' (9-105)"),
    "GSE224683": ("GSE224683 sample name", "Gemma sample name (DLPFC_Sample17)"),
    "GSE26927": ("GSE26927 sample name", "Gemma sample name code (BNE_MS_Control_CMS-RZ84)"),
    "GSE37981": ("GSE37981 title", "GEO sample title (Brain_Control_subject2)"),
    "GSE112523": ("GSE112523 'project.sampleid'", "GEO characteristic 'project.sampleid'"),
    "GSE93987": ("GSE93987 title ID", "GEO sample title '<id>-<CR|SZ>-<layer>-<replicate>'"),
    "GSE87610": ("GSE87610 title ID", "GEO sample title '<id>-<dx>-<layer>-<replicate>' ('subject <id>' in the description)"),
    "GSE12654": ("GSE12654 title number", "GEO sample title '<dx>, BA10-<n>'"),
    "GSE12649": ("GSE12649 title number", "GEO sample title '<dx>, BA46-<n>'"),
    "GSE62191": ("GSE62191 title", "GEO sample title (SZ28)"),
    "BrainSeq Phase 2": ("LIBD BrNum", "LIBD Br number"),
    "PsychAD": ("PsychAD donor_id", "PsychAD donor_id (Donor_n)"),
}
# no subject identifier: GSE138082 (GEO 'original_id' differs per subfield sample), GSE12679 (titles '<cell>_rep<n>' do not
# pair across cell types: rep numbers carry different diagnosis, sex, age or PMI), GSE78936 (title numbers per region)
NO_ID = {"GSE138082": "no subject identifier: GEO 'original_id' is per sample (different in each subfield)",
         "GSE12679": "no subject identifier: titles '<cell type>_rep<n>' do not pair across cell types (9 of 11 shared rep numbers carry different dx, sex, age or PMI)",
         "GSE78936": "no subject identifier: title numbers are per region"}


PAIR_NOTE = {"GSE138082_DG": "GEO GSE138082 (series matrix, read 2026-09-25): each of the 13 DG schizophrenia samples carries the "
                            "gender, age, PMI and RIN of one DG control sample, so DG schizophrenia subjects cannot be matched by fingerprint"}


class Dataset:
    def __init__(self, name, series, source, compartment, platform, bank_study):
        self.name, self.series, self.source, self.compartment, self.platform = name, series, source, compartment, platform
        self.study = bank_study
        self.region = REGION[name]
        self.samples = {}      # gsm -> (subject key, record)
        self.id_space = ID_RULES[series][0] if series in ID_RULES else ""
        self.note = ""

    def add(self, gsm, key, dx, sex, age, pmi, ph=""):
        dxc = {"schizophrenia": "SZ", "SZ": "SZ", "control": "CTRL", "CTRL": "CTRL"}.get(dx)
        if dxc is None:
            return
        sexc = {"M": "M", "F": "F", "male": "M", "female": "F", "Male": "M", "Female": "F"}.get(str(sex).strip(), "")
        alo, ahi, _ = age_interval(age) if not blank(age) else (None, None, None)
        try:
            plo, phi = pmi_interval(pmi) if not blank(pmi) else (None, None)
        except ValueError:
            plo, phi = None, None
        self.samples[gsm] = (key or gsm, {"dx": dxc, "sex": sexc, "age_lo": alo, "age_hi": ahi, "pmi_lo": plo, "pmi_hi": phi,
                                          "age_raw": str(age), "pmi_raw": str(pmi), "ph": str(ph), "gsm": gsm})

    def subjects(self):
        """subject key -> record (first sample's fields); n samples per subject."""
        out = {}
        for gsm, (k, r) in self.samples.items():
            out.setdefault(k, r)
        return out

    def has_fp(self):
        S = self.subjects().values()
        return bool(S) and all(r["sex"] and r["age_lo"] is not None for r in S)

    def has_pmi(self):
        return any(r["pmi_lo"] is not None for r in self.subjects().values())


def gemma_json(path):
    return {s["accession"]["accession"]: s for s in json.load(open(path))["data"]}


def gsm_chars(s):
    return {c.get("category"): c.get("value") for c in (s.get("sample", {}).get("characteristics") or [])}


def registry_key(gse, s, cov):
    """Subject key for a registry sample (ID_RULES), or '' when the series has none."""
    name = s["name"]
    if gse in ("GSE92538", "GSE107638", "GSE202537"):
        return cov.get("subject", "")
    if gse == "GSE53987":
        m = re.fullmatch(r"(\w+?)_(\w+?)_(\d+)", name)
        return f"{m.group(1)}_{m.group(3)}"
    if gse in ("GSE80655", "GSE17612", "GSE21935", "GSE35977", "GSE35974"):
        return name.split("_")[0]
    if gse in ("GSE21138", "GSE224683"):
        return name.split("_")[-1]
    if gse == "GSE208338":
        return gsm_chars(s).get("sample", "")
    if gse == "GSE26927":
        return name.split("_")[-1]
    if gse == "GSE37981":
        return name
    if gse == "GSE112523":
        return gsm_chars(s).get("project.sampleid", "")
    return ""


def platform_of(sj, gsms):
    p = set()
    for g in gsms:
        s = sj[g]
        p.add((s.get("originalPlatform") or s["arrayDesign"])["shortName"])
    return "+".join(sorted(p))


def load_registry():
    reg = read_tsv(use(REG))
    used = collections.defaultdict(list)
    for r in read_tsv(use(USED)):
        used[r["dataset"]].append(r["gsm"])
    cp = lambda p: p if p.startswith("/") else os.path.normpath(os.path.join(CAS, p))
    out = {}
    for r in reg:
        sj = gemma_json(use(cp(r["samples_json"])))
        cov = {x["gsm"]: x for x in read_tsv(use(cp(r["cov"])))}
        d = Dataset(r["dataset"], r["gse"], "registry", r["compartment"], platform_of(sj, used[r["dataset"]]), STUDY[r["gse"]])
        d.cluster = r["cluster"]
        for g in used[r["dataset"]]:
            c = cov[g]
            key = registry_key(r["gse"], sj[g], c) if r["gse"] not in NO_ID else ""
            d.add(g, key, c["dx"], c["sex"], c["age"], c.get("PMI", ""), c.get("pH", ""))
        d.note = NO_ID.get(r["gse"], "")
        out[d.name] = d
    return out


def load_libd():
    rows = read_tsv(use(LIBD))
    out = {}
    for reg, name in (("DLPFC", "LIBD_DLPFC"), ("HIPPO", "LIBD_HIPPO")):
        d = Dataset(name, "BrainSeq Phase 2", "libd_brainseq", "bulk", "RNA-seq, RiboZero (BrainSeq Phase 2)", STUDY["BrainSeq Phase 2"])
        d.cluster = "cmc+mssm+pitt+libd"
        for r in rows:
            if r["region"] == reg:
                d.add(r["RNum"], r["BrNum"], {"Schizo": "SZ", "Control": "CTRL"}[r["dx"]], r["sex"], r["age"], r["PMI"])
        d.note = "adult (age > 17) Phase 2 samples; PMI from the Phase 1 CSV only"
        out[name] = d
    return out


def load_psychad():
    if not os.path.exists(PSYCHAD):
        return {}
    rows = read_tsv(use(PSYCHAD))
    out = {}
    for bank, name in (("H", "PsychAD_HBCC"), ("M", "PsychAD_MSSM")):
        d = Dataset(name, "PsychAD", "psychad_cellxgene", "snRNA-seq pseudobulk", "10x snRNA-seq (PsychAD)", "")
        d.cluster = "cmc+mssm+pitt+libd"
        seen = set()
        for r in rows:
            if r["source_bank"] != bank or r["donor_id"] in seen:
                continue
            if r["dx_gemma"] == "SZ":
                dx = "SZ"
            elif r["dx_gemma"] in ("CTRL", "CTRL_SOMATIC") or r["paper_dx"] == "CTRL":
                dx = "CTRL"
            else:
                continue
            seen.add(r["donor_id"])
            d.add(r["donor_id"], r["donor_id"], dx, r["paper_sex"], r["paper_age"], r["paper_pmi"])
        d.note = "donors of the PsychAD pseudobulk rule (scripts/build_pseudobulk.py); age, sex, PMI from Supp. Table 1 (paper_*)"
        out[name] = d
    return out


def geo_samples(gse):
    p = use(os.path.join(GL, GEO[gse]))
    rows = [l.rstrip("\n").split("\t") for l in gzip.open(p, "rt") if l.startswith("!Sample_")]
    n = len(rows[0]) - 1
    S = [collections.defaultdict(list) for _ in range(n)]
    for r in rows:
        for i, v in enumerate(r[1:]):
            S[i][r[0]].append(v.strip('"'))
    return S


def load_geo_sets():
    out = {}
    # GSE12654, GSE12649: dx in characteristics; title '<dx>, BA10-<n>'
    for gse, reg in (("GSE12654", "BA10"), ("GSE12649", "BA46")):
        d = Dataset(gse, gse, "fairness_lists a4 (GEO)", "bulk", "", STUDY[gse])
        for s in geo_samples(gse):
            d.platform = s["!Sample_platform_id"][0]
            dx = s["!Sample_characteristics_ch1"][1]
            key = s["!Sample_title"][0].split("-")[-1]
            d.add(s["!Sample_geo_accession"][0], key, dx, "", "", "")
        d.note = "GEO series matrix: diagnosis only (no age, sex or PMI)"
        out[gse] = d
    d = Dataset("GSE62191", "GSE62191", "fairness_lists a4 (GEO)", "bulk", "", STUDY["GSE62191"])
    for s in geo_samples("GSE62191"):
        d.platform = s["!Sample_platform_id"][0]
        ch = dict(x.split(": ", 1) for x in s["!Sample_characteristics_ch1"] if ": " in x)
        dx = {"schizophrenia": "SZ", "healthy control": "CTRL", "control": "CTRL"}.get(ch.get("disease state", "").lower(), "")
        d.add(s["!Sample_geo_accession"][0], s["!Sample_title"][0], dx, "", ch.get("age", "").replace(" yr", ""), "")
    d.note = "GEO series matrix: diagnosis and age only (no sex or PMI)"
    out["GSE62191"] = d
    for gse in ("GSE93987", "GSE87610"):
        for layer in ("L3", "L5"):
            d = Dataset(f"{gse}_{layer}", gse, "fairness_lists a4 (GEO)", f"LCM layer {layer[1]} pyramidal neurons", "", STUDY[gse])
            for s in geo_samples(gse):
                t = s["!Sample_title"][0].split("-")
                if t[2] != layer:
                    continue
                d.platform = s["!Sample_platform_id"][0]
                dx = {"CR": "CTRL", "UC": "CTRL", "SZ": "SZ"}.get(t[1], "")
                d.add(s["!Sample_geo_accession"][0], t[0], dx, "", "", "")
            d.note = "GEO series matrix: identifier and diagnosis only; replicate arrays collapsed by identifier"
            out[d.name] = d
    for reg, lab in (("BA9", "BA9"), ("BA11", "BA11"), ("BA24", "BA24")):
        d = Dataset(f"GSE78936_{lab}", "GSE78936", "fairness_lists a4 (GEO)", "bulk", "", STUDY["GSE78936"])
        for s in geo_samples("GSE78936"):
            t = s["!Sample_title"][0].split("_")
            if t[1] != reg:
                continue
            d.platform = s["!Sample_platform_id"][0]
            d.add(s["!Sample_geo_accession"][0], "", {"SZ": "SZ", "Control": "CTRL"}.get(t[0], ""), "", "", "")
        d.note = NO_ID["GSE78936"] + "; GEO series matrix: diagnosis only"
        out[d.name] = d
    return out


def composite(name, parts, compartment):
    """Union of registry datasets of one series. With an identifier: union by identifier. Without one: samples of the
    parts matched one-to-one at the exact fingerprint level count as one subject (the parts' own pair evidence)."""
    a, b = parts
    d = Dataset(name, a.series, "composite of " + " + ".join(p.name for p in parts), compartment,
                "+".join(sorted({a.platform, b.platform})), a.study)
    d.cluster = a.cluster
    if a.id_space:
        for p in parts:
            for g, (k, r) in p.samples.items():
                d.samples[g] = (k, r)
        d.note = "union of the parts by identifier"
    else:
        SA, SB = list(a.subjects().items()), list(b.subjects().items())
        m = M.max_matching(M.matched_sets([r for _, r in SA], [r for _, r in SB])["exact"])
        link = {SB[j][0]: SA[i][0] for i, j in m.items()}
        for g, (k, r) in a.samples.items():
            d.samples[g] = (k, r)
        for g, (k, r) in b.samples.items():
            d.samples[g] = (link.get(k, k), r)
        d.note = f"union of the parts; {len(link)} samples of {b.name} linked to {a.name} subjects by exact one-to-one fingerprint match"
    return d


def check_ids(ds):
    """Every identifier within a dataset or shared across datasets of one identifier space carries one dx, sex, age, PMI."""
    out = []
    by_space = collections.defaultdict(lambda: collections.defaultdict(set))
    for d in ds.values():
        if not d.id_space or d.source.startswith("composite"):
            continue
        for g, (k, r) in d.samples.items():
            by_space[d.id_space][k].add((r["dx"], r["sex"], r["age_raw"] if r["age_lo"] is not None else "", r["pmi_raw"]))
    for sp, ids in sorted(by_space.items()):
        bad = {k: v for k, v in ids.items() if len({x[:2] for x in v}) > 1 or len({x[2] for x in v if x[2]}) > 1 and sp != "LIBD BrNum"}
        out.append((sp, len(ids), len(bad), sorted(bad)[:5]))
    return out


# ---------------------------------------------------------------- banks (moved from fairness_lists/a4_donors.py)
ARRAY = re.compile(r"Array Collection|Microarray Collection|SMRI Array|Array cohort|Collection A\b", re.I)
NPC = re.compile(r"Neuropathology Consortium|\bSNC\b", re.I)
BANK_OVERRIDE = {"33558674": "hbcc", "31288836": None}
NAMED = {"pitt", "libd", "stanley", "aus", "harvard", "mssm", "pritzker", "cmc", "hbcc"}
BANK_SAME = {frozenset(("cmc", "pitt")), frozenset(("cmc", "mssm"))}


class Banks:
    """Bank of a study = its lineage lane (pitt, libd, stanley, aus, harvard, mssm, pritzker, cmc), with Akula 2021 as
    'hbcc' (NIMH HBCC donors, donor_overlap/) and Mendizabal 2019 as unknown (its donor table lists McLean and other
    banks). Stanley studies split by collection: 'Array Collection' / 'Microarray Collection' vs 'Neuropathology
    Consortium', read from the primary list's manifest contrast, else the lineage 'parent' text."""

    def __init__(self, ctx):
        self.ctx = ctx

    def bank(self, k):
        if not k:
            return None
        if k in BANK_OVERRIDE:
            return BANK_OVERRIDE[k]
        if k in self.ctx.units["multibank"]:
            return None
        lane = self.ctx.node.get(k, {}).get("lane")
        return lane if lane in NAMED else None

    def collection(self, k):
        if self.bank(k) != "stanley":
            return ""
        prim = self.ctx.pl["studies"].get(k, {}).get("primary", [])
        for txt in ([self.ctx.specs.get(n, {}).get("contrast", "") for n in prim][:1] + [str(self.ctx.node.get(k, {}).get("parent", ""))]):
            a, n = bool(ARRAY.search(txt)), bool(NPC.search(txt))
            if a != n:
                return "array" if a else "npc"
        return "unknown"


UNKNOWN_NO_BANK = "subject overlap unknown, bank not stated"


def bank_relation(ba, ca, bb, cb):
    """(bank_relation, relation when nothing else is known)."""
    if ba is None or bb is None:
        # untested and no bank to go on: unknown, not separate (the first version called these subject-separate)
        return "bank not stated", UNKNOWN_NO_BANK
    if ba == bb == "stanley":
        if "unknown" in (ca, cb):
            return "Stanley, sub-collection not stated", "same collection, subject overlap unknown"
        if ca == cb:
            return f"Stanley, same sub-collection ({ca})", "same collection, subject overlap unknown"
        return "Stanley, different sub-collections (Array vs Neuropathology Consortium)", "subject-separate"
    if ba == bb:
        return f"same bank ({ba})", "same collection, subject overlap unknown"
    if frozenset((ba, bb)) in BANK_SAME:
        return f"CommonMind with a source bank ({ba}, {bb})", "same collection, subject overlap unknown"
    if frozenset((ba, bb)) == frozenset(("hbcc", "libd")):
        return "HBCC with LIBD (donor_overlap: fingerprint overlap)", "same collection, subject overlap unknown"
    return f"different banks ({ba}, {bb})", "subject-separate"


# ---------------------------------------------------------------- pair evidence
CACHE = os.path.join(HERE, "fingerprint_cache.json")
FP_METHOD = "age and PMI shuffled independently within dx x sex (match.permuted, the original null)"


def fp_key(A, B, seed):
    # the null's method is part of the key: a cached result from an earlier null must not be reused
    return hashlib.sha256(json.dumps([A, B, seed, FP_METHOD], sort_keys=True, default=str).encode()).hexdigest()


def fingerprint(args):
    """donor_overlap/match.py on subject records: k at 'exact' and 'no_pmi', chance mean / q95 / p for each. Results are
    cached by a hash of the two record lists and the seed (fingerprint_cache.json), so reruns only recompute changed pairs."""
    A, B, seed = args
    m = M.matched_sets(A, B)
    obs = {lv: len(M.max_matching(m[lv])) for lv in ("exact", "no_pmi")}
    M.random.seed(seed)
    null = {"exact": [], "no_pmi": []}
    for _ in range(M.N_PERM):
        mp = M.matched_sets(A, M.permuted(B))
        for lv in null:
            null[lv].append(len(M.max_matching(mp[lv])))
    out = {}
    for lv in ("exact", "no_pmi"):
        out[lv] = (obs[lv], statistics.mean(null[lv]), M.q95(null[lv]), (1 + sum(x >= obs[lv] for x in null[lv])) / (M.N_PERM + 1))
    return out


def fp_overlap(k, mu, q, p):
    """Fingerprint overlap rule: k >= 3, k above the chance 95th percentile, permutation p <= 0.01 (at most 1 of 200
    permutations reaches k), at least 3 above the chance mean and at least twice the chance mean."""
    return k >= MIN_K and k > q and p <= MAX_P and k - mu >= MIN_EXCESS and k >= MIN_RATIO * mu


def same_region(a, b):
    if a == b:
        return "y"
    if "not specified" in a + b and ("cortex" in a or "PFC" in a) and ("cortex" in b or "PFC" in b):
        return "unknown"
    return "n"


def classify(a, b, fp, banks, seed_rel):
    Sa, Sb = a.subjects(), b.subjects()
    na, nb = len(Sa), len(Sb)
    row = {"k_shared": "", "evidence": "", "evidence_detail": ""}
    row["note"] = "; ".join(PAIR_NOTE[n] for n in (a.name, b.name) if n in PAIR_NOTE and (a.series != b.series or a.name != b.name))
    ga, gb = set(a.samples), set(b.samples)
    small, big = (ga, gb) if len(ga) <= len(gb) else (gb, ga)
    sreg = same_region(a.region, b.region)
    same_subj = lambda k: k >= SAME_FRAC * min(na, nb)
    subj_label = lambda: "same subjects, same region" if sreg == "y" else "same subjects, different region"
    if small and small <= big:
        row.update(k_shared=min(na, nb), evidence="same samples",
                   evidence_detail=f"all {len(small)} samples of the smaller dataset are in the other")
        return row, "same samples"
    if a.id_space and a.id_space == b.id_space:
        k = len(set(Sa) & set(Sb))
        row.update(k_shared=k, evidence="identifier", evidence_detail=f"{a.id_space}: {k} identifiers in common")
        if k == 0:
            br, _ = seed_rel
            row["evidence_detail"] += f"; {br}"
            return row, "subject-separate"
        return row, subj_label() if same_subj(k) else "partial subject overlap"
    if fp:
        lv = fp["level"]
        k, mu, q, p = fp[lv]
        row.update(k_shared=k, evidence="fingerprint",
                   evidence_detail=f"{lv}: {k} matched one-to-one, chance mean {mu:.1f}, 95th percentile {q}, p {p:.3f}")
        if fp_overlap(k, mu, q, p):
            return row, subj_label() if same_subj(k) else "partial subject overlap"
    st = STATED.get(frozenset((a.name, b.name)))
    if st:
        row.update(k_shared=min(na, nb), evidence="stated", evidence_detail=st)
        return row, subj_label()
    br, rel_ = seed_rel
    if not row["evidence"]:
        row.update(k_shared="", evidence="collection", evidence_detail=f"no identifier or fingerprint; {br}")
    else:
        row["evidence_detail"] += f" (no excess); {br}"
        if rel_ == UNKNOWN_NO_BANK:
            # a fingerprint test with no excess is the only evidence here, and it says separate
            rel_ = "subject-separate"
    return row, rel_


def relation_group(relation, same_series):
    """Reporting group: the relation, with same-subject relations split by whether the two datasets are one series
    (same lab, protocol, usually batch) or two (different studies)."""
    if relation.startswith("same subjects") or relation == "partial subject overlap":
        return f"{relation}, {'same series' if same_series == 'y' else 'across series'}"
    return relation


def fmt(x):
    if isinstance(x, float):
        return f"{x:.4g}"
    return "" if x is None else str(x)


def write(path, rows, head, fields=None):
    fields = fields or list(rows[0])
    with open(path, "w", newline="") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(v) for k, v in r.items()})


# ---------------------------------------------------------------- study level (published lists; fairness_lists a4)
# study -> dataset whose subjects its primary list (primary_lists.json) was computed on
STUDY_DATASET = {"18778695": "GSE21138", "21538462": "GSE21935", "24702465": "GSE37981", "31123247": "GSE53987_BA46",
                 "36302706": "GSE202537_caudate", "32152472": "GSE138082_DG", "31174959": "LIBD_DLPFC",
                 "31288836": "GSE107638_NeuN", "19255580": "GSE17612"}
# published lists of the a4 'not in H1' pairs -> dataset
LIST_DATASET = {"31123247__sz_pfc": "GSE53987_BA46", "31123247__sz_hip": "GSE53987_hippocampus",
                "31123247__sz_str": "GSE53987_striatum", "36302706__caudate": "GSE202537_caudate",
                "36302706__nac": "GSE202537_NAc", "36302706__putamen": "GSE202537_putamen",
                "32152472__dg": "GSE138082_DG", "32152472__ca3": "GSE138082_CA3", "32152472__ca1": "GSE138082_CA1",
                "31174959__dlpfc": "LIBD_DLPFC", "31174959__hippo": "LIBD_HIPPO",
                "31288836__neun_szdeg": "GSE107638_NeuN", "31288836__olig2_szdeg": "GSE107638_OLIG2",
                "19255580__cchpc_ba10": "GSE17612", "21538462__scz_ba22": "GSE21935"}
LIST_PAIRS = [("31123247__sz_pfc", "31123247__sz_hip"), ("31123247__sz_pfc", "31123247__sz_str"),
              ("31123247__sz_hip", "31123247__sz_str"), ("36302706__caudate", "36302706__nac"),
              ("36302706__caudate", "36302706__putamen"), ("36302706__nac", "36302706__putamen"),
              ("32152472__dg", "32152472__ca3"), ("32152472__dg", "32152472__ca1"), ("32152472__ca3", "32152472__ca1"),
              ("31174959__dlpfc", "31174959__hippo"), ("31288836__neun_szdeg", "31288836__olig2_szdeg"),
              ("19255580__cchpc_ba10", "21538462__scz_ba22"), ("16223876__si_table4_dlpfc", "19255580__hbb_ba9")]
# documented facts for units without datasets here
FACTS = {
    frozenset(("16223876__si_table4_dlpfc", "19255580__hbb_ba9")):
        ("same samples", "stated", "", "Maycox 2009 reanalyzed Glatt 2005's arrays (corpus/pooled_datasets.json; tracker Step 1, item 10)"),
}
HBCC_LIBD_FACT = ("partial subject overlap", "fingerprint", "44",
                  "donor_overlap_report.md: 44 of 185 HBCC sgACC subjects (Akula 2021 table) match a LIBD BrainSeq subject "
                  "(Benjamin 2024 table) at the exact fingerprint level, 4.4 expected by chance")


def main():
    ds = {}
    ds.update(load_registry())
    ds.update(load_libd())
    ds.update(load_psychad())
    ds.update(load_geo_sets())
    ds["GSE12679_all"] = composite("GSE12679_all", (ds["GSE12679_neuron"], ds["GSE12679_endothelial"]), "LCM neurons + endothelial cells")
    ds["GSE107638_all"] = composite("GSE107638_all", (ds["GSE107638_NeuN"], ds["GSE107638_OLIG2"]), "sorted NeuN+ and OLIG2+ nuclei")
    for d in ds.values():
        if d.name.startswith("GSE") and not d.platform:
            d.platform = "?"
    for sp, n, nbad, ex in check_ids(ds):
        print(f"identifier check {sp}: {n} identifiers, {nbad} with conflicting dx/sex/age" + (f" {ex}" if nbad else ""))

    sys.path.insert(0, os.path.join(GL, "fairness_lists"))
    import fl_common as C  # noqa: E402
    ctx = C.Ctx()
    use(ctx.pl_path)
    B = Banks(ctx)

    def dbank(d):
        return BANK_DATASET.get(d.name) or B.bank(d.study)

    def dcoll(d):
        return B.collection(d.study) if dbank(d) == "stanley" else ""

    names = list(ds)
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    task = {}
    for x, y in itertools.combinations(names, 2):
        a, b = ds[x], ds[y]
        if a.has_fp() and b.has_fp():
            args = (list(a.subjects().values()), list(b.subjects().values()),
                    SEED + int(hashlib.sha256(f"{x}|{y}".encode()).hexdigest()[:6], 16))
            k = fp_key(*args)
            task[(x, y)] = k
            if k not in cache:
                cache[k] = args
    todo = [k for k, v in cache.items() if isinstance(v, tuple)]
    print(f"fingerprints: {len(task)} pairs, {len(todo)} to compute", flush=True)
    import multiprocessing as mp
    with mp.Pool(int(os.environ.get("NCORES", "6"))) as pool:
        for n, (k, res) in enumerate(zip(todo, pool.imap(fingerprint, [cache[k] for k in todo], chunksize=4))):
            cache[k] = res
            if n % 50 == 0:
                print("fingerprint", n, flush=True)
    json.dump({k: v for k, v in cache.items() if k in set(task.values())}, open(CACHE, "w"))
    rows = []
    for i, (x, y) in enumerate(itertools.combinations(names, 2)):
        a, b = ds[x], ds[y]
        fp = None
        if a.has_fp() and b.has_fp():
            fp = {lv: tuple(v) for lv, v in cache[task[(x, y)]].items()}
            fp["level"] = "exact" if a.has_pmi() and b.has_pmi() else "no_pmi"
        br = bank_relation(dbank(a), dcoll(a), dbank(b), dcoll(b))
        r, relation = classify(a, b, fp, B, br)
        na, nb = len(a.subjects()), len(b.subjects())
        k = r["k_shared"]
        if relation == "partial subject overlap":
            text = f"partial subject overlap ({k} of {na} and {nb})"
        elif relation.startswith("same subjects"):
            text = f"{relation} ({k} of {na} and {nb})"
        else:
            text = relation
        row = {"dataset_a": x, "dataset_b": y, "series_a": a.series, "series_b": b.series,
               "region_a": a.region, "region_b": b.region, "compartment_a": a.compartment, "compartment_b": b.compartment,
               "platform_a": a.platform, "platform_b": b.platform,
               "n_samples_a": len(a.samples), "n_samples_b": len(b.samples), "n_subjects_a": na, "n_subjects_b": nb,
               "k_shared": k, "relation": relation, "relation_group": relation_group(relation, "y" if a.series == b.series else "n"),
               "relation_text": text, "evidence": r["evidence"],
               "evidence_detail": r["evidence_detail"], "note": r["note"],
               "same_series": "y" if a.series == b.series else "n",
               "same_platform": "y" if a.platform == b.platform else "n",
               "same_region": same_region(a.region, b.region),
               "bank_a": dbank(a) or "not stated", "bank_b": dbank(b) or "not stated", "bank_relation": br[0],
               "fp_level": fp["level"] if fp else "",
               "fp_k_exact": fp["exact"][0] if fp else "", "fp_chance_mean_exact": fp["exact"][1] if fp else "",
               "fp_chance_q95_exact": fp["exact"][2] if fp else "", "fp_p_exact": fp["exact"][3] if fp else "",
               "fp_k_no_pmi": fp["no_pmi"][0] if fp else "", "fp_chance_mean_no_pmi": fp["no_pmi"][1] if fp else "",
               "fp_chance_q95_no_pmi": fp["no_pmi"][2] if fp else ""}
        rows.append(row)
        if i % 50 == 0:
            print(i, x, y, relation, k, flush=True)


    units = []
    for n, d in ds.items():
        S = d.subjects()
        units.append({"unit": n, "kind": "dataset", "source": d.source, "series": d.series, "study": d.study,
                      "region": d.region, "compartment": d.compartment, "platform": d.platform,
                      "n_samples": len(d.samples), "n_subjects": len(S),
                      "n_SZ": sum(r["dx"] == "SZ" for r in S.values()), "n_CTRL": sum(r["dx"] == "CTRL" for r in S.values()),
                      "id_space": d.id_space, "id_rule": ID_RULES.get(d.series, ("", ""))[1] if d.id_space else "",
                      "fingerprint_fields": ("dx,sex,age" + (",pmi" if d.has_pmi() else "")) if d.has_fp() else "",
                      "bank": dbank(d) or "not stated", "stanley_collection": dcoll(d), "dataset": "", "note": d.note})

    # ---- study level
    base = ctx.reps(C.MAIN)
    for k in base:
        units.append({"unit": k, "kind": "study", "source": "primary_lists.json (published, brain)", "series": "",
                      "study": k, "region": "", "compartment": "", "platform": "", "n_samples": "", "n_subjects": "",
                      "n_SZ": "", "n_CTRL": "", "id_space": "", "id_rule": "", "fingerprint_fields": "",
                      "bank": B.bank(k) or "not stated", "stanley_collection": B.collection(k),
                      "dataset": STUDY_DATASET.get(k, ""), "note": ctx.label(k)})
    drow = {}
    for r in rows:
        drow[(r["dataset_a"], r["dataset_b"])] = r
        drow[(r["dataset_b"], r["dataset_a"])] = r

    def from_dataset(da, db):
        r = drow[(da, db)]
        return {"relation": r["relation"], "relation_group": r["relation_group"], "relation_text": r["relation_text"], "evidence": r["evidence"],
                "evidence_detail": f"dataset pair {da} / {db}: {r['evidence_detail']}", "k_shared": r["k_shared"],
                "n_subjects_a": r["n_subjects_a"] if r["dataset_a"] == da else r["n_subjects_b"],
                "n_subjects_b": r["n_subjects_b"] if r["dataset_a"] == da else r["n_subjects_a"],
                "same_series": r["same_series"], "bank_relation": r["bank_relation"]}

    srows = []
    for ka, kb in itertools.combinations(base, 2):
        br = bank_relation(B.bank(ka), B.collection(ka), B.bank(kb), B.collection(kb))
        da, db = STUDY_DATASET.get(ka), STUDY_DATASET.get(kb)
        if da and db:
            s = from_dataset(da, db)
        elif {B.bank(ka), B.bank(kb)} == {"hbcc", "libd"}:
            s = {"relation": HBCC_LIBD_FACT[0], "relation_group": relation_group(HBCC_LIBD_FACT[0], "n"), "relation_text": "partial subject overlap (44 of 185 HBCC sgACC subjects)",
                 "evidence": HBCC_LIBD_FACT[1], "evidence_detail": HBCC_LIBD_FACT[3], "k_shared": HBCC_LIBD_FACT[2],
                 "n_subjects_a": "", "n_subjects_b": "", "same_series": "n", "bank_relation": br[0]}
        else:
            s = {"relation": br[1], "relation_group": br[1], "relation_text": br[1], "evidence": "collection",
                 "evidence_detail": f"no subject-level data for both studies here; {br[0]}", "k_shared": "",
                 "n_subjects_a": "", "n_subjects_b": "", "same_series": "n", "bank_relation": br[0]}
        srows.append(dict({"unit_a": ka, "unit_b": kb, "kind": "study", "label_a": ctx.label(ka), "label_b": ctx.label(kb),
                           "dataset_a": da or "", "dataset_b": db or "", "bank_a": B.bank(ka) or "not stated",
                           "bank_b": B.bank(kb) or "not stated"}, **s))
    for la, lb in LIST_PAIRS:
        f = FACTS.get(frozenset((la, lb)))
        da, db = LIST_DATASET.get(la), LIST_DATASET.get(lb)
        if da and db:
            s = from_dataset(da, db)
        else:
            s = {"relation": f[0], "relation_group": f[0], "relation_text": f[0], "evidence": f[1], "evidence_detail": f[3], "k_shared": f[2],
                 "n_subjects_a": "", "n_subjects_b": "", "same_series": "", "bank_relation": ""}
        srows.append(dict({"unit_a": la, "unit_b": lb, "kind": "list", "label_a": ctx.label(la.split("__")[0]),
                           "label_b": ctx.label(lb.split("__")[0]), "dataset_a": da or "", "dataset_b": db or "",
                           "bank_a": B.bank(la.split("__")[0]) or "not stated", "bank_b": B.bank(lb.split("__")[0]) or "not stated"}, **s))

    head = [f"# built {TODAY} by subject_relations/subject_relations.py (tracker, Vocabulary, 2026-09-25); rules in the script docstring and README.md",
            f"# fingerprint: donor_overlap/match.py (N_PERM {M.N_PERM}), intervals of donor_overlap/build_tables.py; overlap counted when k >= {MIN_K}, k > chance 95th percentile, permutation p <= {MAX_P}, k - chance mean >= {MIN_EXCESS} and k >= {MIN_RATIO} x chance mean; 'same subjects' when k >= {SAME_FRAC:.0%} of the smaller dataset's subjects"]
    head += [f"# input: {p} sha256 {h}" for p, h in sorted(set(INPUTS))]
    write(os.path.join(HERE, "subject_relations.tsv"), rows, head)
    write(os.path.join(HERE, "subject_relations_units.tsv"), units, head)
    write(os.path.join(HERE, "study_relations.tsv"), srows, head + [
        "# study pairs: fairness_lists a4 part A representatives (published, brain); list pairs: a4 'not in H1' pairs. Relation from subject_relations.tsv where both map to a dataset (STUDY_DATASET, LIST_DATASET), else documented facts, else the collection presumption"],
        fields=["unit_a", "unit_b", "kind", "label_a", "label_b", "dataset_a", "dataset_b", "bank_a", "bank_b", "relation",
                "relation_group", "relation_text", "evidence", "evidence_detail", "k_shared", "n_subjects_a", "n_subjects_b", "same_series",
                "bank_relation"])
    c = collections.Counter(r["relation_group"] for r in rows)
    for rl, n in sorted(c.items()):
        print(f"{n:5d}  {rl}")
    print(collections.Counter(r["relation"] for r in srows))


if __name__ == "__main__":
    main()
