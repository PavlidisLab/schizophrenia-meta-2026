"""Assay technology of each study's primary list: microarray, rnaseq, single_cell, imputed, other, proteomics (batch 10).

Paul, 2026-09-24: "we need to do an analysis that is isolated for rna-seq vs microarray (excluding single
cell for now)". recurrence.py and h1.py read technology.json to build per-technology scenarios.

Order of evidence, first that applies:
  1. OVERRIDES below (each with its evidence).
  2. lineage node single_cell (lane "sc" before 2026-09-25) -> single_cell.
  3. universes.json list type: "platform" (a GEO GPL annotation) -> microarray; "rnaseq_default" -> rnaseq.
  4. keywords in the lineage node's platform text (sequencing words vs array words); both or neither -> unknown.

Microarray studies also get array_generation (early / whole / noncoding / not stated; rule at ARRAY_EARLY_MAX), which the
lineage page reads for its technology tiers (TIERS; build_lineage.py imports array_generation, tier, technology_of_text).

    python classify_technology.py   # writes technology.json
"""
import datetime
import hashlib
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
LINEAGE = os.path.join(HERE, "..", "corpus", "lineage.json")

SEQ = re.compile(r"rna-?seq|sequenc|hiseq|nextseq|novaseq|solid|genome analyzer|solexa|ion proton|ampliseq", re.I)
PROT = re.compile(r"mass spectrometry|lc-?(ms|srm)|proteomic", re.I)
ARR = re.compile(r"affymetrix|array|beadchip|genechip|agilent|codelink|u133|hg-u|illumina ht-12|humanref", re.I)

OVERRIDES = {
    # TWAS: imputed expression, not an assay
    "30911161": ("imputed", "Huckins 2019: TWAS on CommonMind-trained predictors"),
    "29632383": ("imputed", "Gusev 2018: TWAS"),
    # platform text names both technologies; the primary list is from the one given here
    "19088852": ("microarray", "Harris 2008: reproduced list from U133 Plus 2.0 / CodeLink arrays"),
    "36379214": ("microarray", "Namkung 2023: primary list from GSE215985, U133 Plus 2.0"),
    "31073119": ("microarray", "Bowen 2019: primary list on GPL10558 (HT-12 v4)"),
    "30016334": ("microarray", "Hagenauer 2018: Gemma GSE92538, U133A / U133 Plus 2.0"),
    "23147385": ("microarray", "Chen 2013: Gemma GSE35974, Human Gene 1.0 ST"),
    "22832356": ("microarray", "Tang 2011: Gemma GSE21138, U133 Plus 2.0"),
    "22397921": ("microarray", "Tang 2012: Gemma GSE21138, U133 Plus 2.0"),
    "30576541": ("microarray", "Piras 2019: Human Gene 1.1 ST / Exon 1.0 ST arrays"),
    "36932057": ("rnaseq", "Schulmann 2023: primary list from CMC_HBCC RNA-seq"),
    "41210297": ("rnaseq", "Mas 2025: Ion AmpliSeq Transcriptome (targeted whole-transcriptome sequencing)"),
    # no platform text on the lineage node
    "25119539": ("microarray", "Durrenberger 2015: Gemma GSE26927 on GPL6255 (Illumina humanRef-8 v2), technologyType ONECOLOR"),
    "25487697": ("microarray", "de Baumont 2015: Gemma GSE62191 on GPL4133 (Agilent 4x44K), technologyType DUALMODE"),
    "31399567": ("rnaseq", "Seney 2019: CommonMind / Pittsburgh RNA-seq"),
    "36652379": ("rnaseq", "Clifton 2023: RNA-seq (CommonMind/PsychENCODE)"),
    "38781321": ("rnaseq", "Xia 2024: RNA-seq (CommonMind/PsychENCODE)"),
    "30545856": ("rnaseq", "Gandal 2018 (Science): PsychENCODE RNA-seq"),
    "29439242": ("microarray", "Gandal 2018 (Science, cross-disorder): primary list scz_meta_fdr05 is the microarray meta-analysis (lineage text: 'Published microarray studies (refs 3-11); RNAseq replication')"),
    "34957416": ("rnaseq", "Jin 2021: DESeq2 bulk RNA-seq"),
    "41714609": ("rnaseq", "Smith 2026: sgACC RNA-seq"),
    # platform text matched both keyword sets but names only sequencing
    "33585819": ("rnaseq", "Gatta 2021: Illumina TruSeq stranded RNA-seq"),
    "36319771": ("rnaseq", "Benjamin 2022: Illumina HiSeq 3000 total RNA-seq"),
    "31053723": ("rnaseq", "Pai 2019: Gemma GSE112523, NextSeq 500 RNA-seq"),
    # batch 10 (2026-09-24)
    "17464717": ("other", "Dean 2007: restriction fragment differential display (displayPROFILE), neither array nor RNA-seq"),
    "42693176": ("single_cell", "Wu 2026: public frontal cortex snRNA-seq datasets, pseudobulk (paper check 2026-09-27)"),
    "37031222": ("other", "Medina 2023: the primary lists are TaqMan qPCR array cards, 115 targets (own BA10 data; paper check 2026-09-27)"),
    "36192459": ("single_cell", "Puvogel 2022: 10x single-nucleus RNA-seq (paper check 2026-09-27: 'Single-nucleus cDNA libraries were constructed ...')"),
    "41757081": ("other", "Kwon 2026: Visium spatial transcriptomics (Visium-SPG), pseudobulk DE; neither bulk RNA-seq nor array"),
    # 2026-10-01: the eight studies left 'unknown' by the keyword rule, settled from the paper or the GEO record
    # (quotes checked with verify_quote.py: corpus_unknowns_2026-10-01/quotes_checked.json)
    "16103141": ("microarray", "Frank 2005: 'we performed a comprehensive microarray-based analysis of HERV transcriptional activity in human brains' (custom retrovirus pol chip)"),
    "19721432": ("microarray", "Beveridge 2010: 'Microarrays were prepared using anti-sense Locked Nucleic Acid (LNA) oligonucleotides (Exiqon ...) ... containing 322 human miRNA sequences'"),
    "21183010": ("other", "Moreau 2011: 'Real-time qPCR was performed using the TaqMan® ... single tube human miRNA assay panel (Sanger miRBase v9.2)'"),
    "24475125": ("other", "Smalheiser 2014: the primary list (Table 1) is from the TaqMan Low Density Array: 'The TLDA platform version 2.0 measures miRNA expression over two different 384-well plates; we employed the A plate'"),
    "26183697": ("other", "Wei 2015: the primary list (Table 2) is qRT-PCR of eight candidate miRNAs: 'The qRT-PCR validation showed that of these eight miRNAs identified, only miR-130b and miR-193a-3p showed significant differences (Table 2)'"),
    "26710095": ("microarray", "Zheutlin 2016: 'Technical replicates were run on all samples using the Illumina Human WG6 v3.0 chip'"),
    "27117414": ("rnaseq", "Topol 2016: the primary lists are Gemma reanalyses of GSE63734, GSE63738 and GSE80170, each 'Expression profiling by high throughput sequencing' on GPL16791 (GEO, read 2026-10-01)"),
    "27219343": ("rnaseq", "Ellis 2016: 'RNA-seq, re-used SMRI Array Collection BA24 data from Zhao et al. 2015' (lineage platform text; 'Array Collection' is the Stanley cohort's name, not a microarray)"),
    "doi_10.64898_2026.08.14.744920": ("rnaseq", "Casas 2026 (bioRxiv): 'a standard RNA sequencing service, PolyA, 30M reads per sample, paired-end' (Supplementary methods, media-1.pdf p. 3)"),
}
# What the assay is, for every study whose technology is 'other' (the Corpus panel names it instead of a bare "other"; 2026-10-01)
ASSAY = {
    "17464717": "differential display",
    "37031222": "qPCR array (TaqMan cards)",
    "41757081": "spatial transcriptomics (Visium)",
    "21183010": "qPCR panel (TaqMan miRNA assays)",
    "24475125": "qPCR array (TaqMan Low Density Array)",
    "26183697": "qPCR of candidates",
}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


# ---------------------------------------------------------------- array generation (Paul, 2026-09-25: "I don't see a way to do it
# by technology (old microarray, better microarrays, rnaseq)"). A microarray is an EARLY array when it measures fewer than
# ARRAY_EARLY_MAX genes (partial-genome designs: spotted / cDNA / membrane arrays, HuGeneFL, U95, U133A alone, CodeLink UniSet
# 20K), and a WHOLE-GENOME array otherwise (U133A+B, U133 Plus 2.0, X3P, U219, Gene / Exon ST, Clariom, Illumina HumanRef-8 /
# Human-6 / HT-12, Agilent Whole Human Genome, CodeLink Whole Genome). The cut sits in the gap between the platform gene
# counts in universe_sets/ (U133A 13,569 and CodeLink 20K 14,476 genes; next are HumanRef-8 v2 17,143 and HT-12 v3 18,978).
# First rule that applies:
#   1. the gene count of the platform gene set of the study's primary list (universes.json type platform, n_membership);
#   2. a platform named in the text (ARRAY_PLATFORMS, each placed by its design);
#   3. a stated gene count in the platform text ("12,000 genes", "1127 brain-relevant genes"; probe or probe-set counts are
#      not gene counts and are not read);
#   4. custom, spotted, cDNA or membrane arrays with no count -> early (these were partial-genome designs);
#   5. otherwise the generation is not stated.
# miRNA and lncRNA arrays are their own class ("noncoding"): they do not measure genes.
ARRAY_EARLY_MAX = 16000
ARRAY_PLATFORMS = [  # (pattern, generation); first match wins, so the specific names come before the family names
    (r"u133\s*a\s*(and|\+|/|&)\s*(hu|hg-?)?u?133\s*b|hu133a and hu133b|u133 set a and b", "whole"),
    (r"u133\s*(\+|p\b|plus)|u133plus|plus[\s_]*2|u133_x3p|x3p|hg-?u219|\bu219\b|u133\+ pm", "whole"),
    (r"gene[\s-]*(1\.0|1\.1|2\.0)[\s-]*st|hugene-|exon[\s-]*1\.0|exon arrays?|clariom|primeview", "whole"),
    (r"ht-?12|humanht|humanref-?8|humanref8|human-?6\b|human 6v2|ht-?8 v3|whole human genome|4x44k|codelink[^;]*whole genome|whole genome bioarray", "whole"),
    (r"hugenefl|hu6800|u95|hg-?u133\s*a\b|hg-?u133a|hu-?133\s*a\b|u133\s*a\b|u133a\b|hgu133a|codelink[^;]*(uniset|20k)", "early"),
    (r"compugen[^;]*19\s*k|\b19\s*k\b", "whole"),  # Compugen 19K oligonucleotide set: over 16,000 genes (Bowden 2008, McCurdy 2006)
    (r"atlas|clontech|incyte|unigem|nytran|nylon|membrane|neuroarray|macroarray", "early"),
]
# Mixed generations (Paul, 2026-09-27: "we just need to be consistent. I would go by what the data are. u133a and u95 with
# whole genome is indeed mixed"): a study whose data span early and whole-genome arrays is "mixed". The platforms come from
# the GPL ids of its universe; U133A with U133B (GPL96 + GPL97) is the A+B set, one whole-genome design.
GPL_GEN = {"GPL8300": "early", "GPL91": "early", "GPL96": "early", "GPL80": "early", "GPL97": "whole-with-GPL96",
           "GPL570": "whole", "GPL13158": "whole", "GPL5175": "whole", "GPL5188": "whole", "GPL6244": "whole", "GPL6883": "whole",
           "GPL6947": "whole", "GPL10558": "whole", "GPL6255": "whole", "GPL6104": "whole", "GPL2507": "whole", "GPL6102": "whole"}


def gpl_generation(source):
    """'mixed', 'early', 'whole' or None from the GPL ids named in a universe source."""
    g = set(re.findall(r"GPL\d+", str(source or "")))
    if not g or not g <= set(GPL_GEN):
        return None
    if {"GPL96", "GPL97"} <= g:  # the U133 A+B set counts as one whole-genome design
        g = (g - {"GPL96", "GPL97"}) | {"GPL570"}
    gens = {GPL_GEN[x] for x in g} - {"whole-with-GPL96"}
    return "mixed" if len(gens) > 1 else (gens.pop() if gens else None)


# Generation by what the data are, where the list's universe cannot tell (paper check 2026-09-27,
# corpus/corpus_node_paper_checks/checked_*.tsv; Paul: "go by what the data are").
GEN_OVERRIDES = {
    "30016334": ("mixed", "Hagenauer 2018: the Pritzker data (GSE92538) are U133A and U133 Plus 2.0 arrays; both Gemma subsets are primary lists"),
    "27450777": ("mixed", "Hess 2016: the mega-analysis pools U133A-only datasets with U133 Plus 2.0, Gene/Exon ST and HT-12 datasets"),
    "26285829": ("mixed", "Bergon 2015: Table 1 datasets are U95Av2, U133A, U133 Plus 2.0, HumanRef-8 and HT-12"),
    "29439242": ("mixed", "Gandal 2018 (Science): the microarray meta-analysis pools Human Gene 1.0 ST (GSE35978), U133A (GSE12649) and U133 Plus 2.0 (GSE17612, GSE21138) sets"),
    "17326821": ("noncoding", "Perkins 2007: the study's own data are a custom miRNA array"),
    "16103141": ("noncoding", "Frank 2005: a custom chip of HERV pol probes (retroviral elements, not genes; 2026-10-01)"),
    "26733343": ("whole", "Xu 2016: the list is mRNA from the coding-transcript probes (8,594 mRNAs detected)"),
}
NONCODING_ARRAY = re.compile(r"mirna|microrna|lncrna|circrna", re.I)
MRNA_ARRAY = re.compile(r"mrna|u133|gene 1\.|exon|ht-12|humanref", re.I)
GENE_COUNT = re.compile(r"(\d{1,2},?\d{3}|\d{3,4})\s*(?:[\w-]+\s+){0,3}?(genes|gene transcripts|nuclear genes|unigene clusters)\b", re.I)
EARLY_WORDS = re.compile(r"custom|spotted|cdna|in-?house|candidate gene array", re.I)


def array_generation(platform_text, primary_universes=()):
    """(generation, evidence) for a microarray study: 'early', 'whole', 'noncoding' or 'not stated'. primary_universes: the
    universes.json entries of the study's primary lists. Rules in the block comment above."""
    gg = {gpl_generation(u.get("source")) for u in primary_universes if u.get("type") in ("platform", "union")} - {None}
    if "mixed" in gg or {"early", "whole"} <= gg:
        return "mixed", "platforms of different generations in the list's universe: " + "; ".join(sorted({str(u.get("source")) for u in primary_universes if u.get("source")}))[:200]
    for u in primary_universes:
        if u.get("type") == "platform" and u.get("n_membership"):
            g = u["n_membership"]
            return ("early" if g < ARRAY_EARLY_MAX else "whole"), f"platform gene set {u.get('source')}: {g:,} genes (universes.json)"
    t = platform_text or ""
    if NONCODING_ARRAY.search(t) and not MRNA_ARRAY.search(t):
        return "noncoding", "miRNA / lncRNA array (platform text)"
    for pat, gen in ARRAY_PLATFORMS:
        m = re.search(pat, t, re.I)
        if m:
            return gen, f"platform named in the text: '{m.group(0)}'"
    m = GENE_COUNT.search(t)
    if m:
        g = int(m.group(1).replace(",", ""))
        return ("early" if g < ARRAY_EARLY_MAX else "whole"), f"gene count stated in the platform text: '{m.group(0)}'"
    m = EARLY_WORDS.search(t)
    if m:
        return "early", f"custom / spotted / cDNA array, no gene count: '{m.group(0)}'"
    return "not stated", "platform text names no array: " + (t[:80] or "(none)")


# Technology tiers for the lineage page (Colour by: Technology). Keys are stable; labels are what the page shows.
TIERS = [
    ("array_early", "Low-coverage microarrays (500 to 16,000 genes: two-channel cDNA; one-channel HuGeneFL, U95, U133A alone)"),
    ("array_whole", "High-coverage microarrays (16,000 genes or more: U133A+B, U133 Plus 2.0, Gene/Exon ST, HT-12)"),
    ("array_unstated", "Microarray, generation not stated"),
    ("array_noncoding", "miRNA / lncRNA arrays"),
    ("panel", "Targeted panels (under 500 pre-selected genes)"),
    ("rnaseq", "Bulk RNA-seq"),
    ("single_cell", "Single-cell / single-nucleus"),
    ("proteomics", "Proteomics"),
    ("other", "Other (qPCR, in situ, differential display, TWAS, spatial)"),
    ("mixed", "Mixed: data from several technologies or array coverages (e.g. U133A or U95 with high-coverage arrays)"),
    ("unknown", "Not stated"),
]


def tier(technology, generation=None):
    """The page's technology tier from technology.json's technology and array_generation."""
    if technology == "microarray":
        return {"early": "array_early", "whole": "array_whole", "noncoding": "array_noncoding", "mixed": "mixed"}.get(generation, "array_unstated")
    if technology == "mixed":
        return "mixed"
    if technology in ("rnaseq", "single_cell", "proteomics"):
        return technology
    if technology in ("imputed", "other"):
        return "other"
    return "unknown"


def technology_of_text(text):
    """Rule 4 of classify() on a platform text alone (lineage nodes with no primary list): 'rnaseq', 'microarray' or 'unknown'."""
    seq, arr = bool(SEQ.search(text or "")), bool(ARR.search(text or ""))
    return ("rnaseq" if seq else "microarray") if seq != arr else "unknown"


def data_tier(k, node):
    """(technology, generation, evidence) from the GEO platforms of a re-use study's pooled datasets (corpus/pooled_geo.py;
    Paul, 2026-09-27: "go by what the data are"), or None. Explicit OVERRIDES / GEN_OVERRIDES win."""
    dt = DATA_TIER.get(k) if str(node.get("category")) == "3" and k not in OVERRIDES and k not in GEN_OVERRIDES else None
    if not (dt and dt.get("tier_by_data")):
        return None
    if dt["tier_by_data"] == "mixed":
        t, g = ("microarray", "mixed") if all(x.startswith("array_") for x in dt["dataset_tiers"]) else ("mixed", None)
    else:
        t, g = DATA_TO_TECH[dt["tier_by_data"]]
    return t, g, "GEO platforms of the pooled datasets (corpus/pooled_geo.json): " + ", ".join(dt["dataset_tiers"])


def classify(k, s, node, uni):
    if k in OVERRIDES:
        return OVERRIDES[k][0], "override: " + OVERRIDES[k][1]
    d = data_tier(k, node)
    if d:
        return d[0], d[2]
    if node.get("single_cell") or node.get("lane") == "sc":  # single_cell flag: sc studies sit in bank lanes since 2026-09-25
        return "single_cell", "lineage lane sc"
    ut = [uni.get(n, {}).get("type") for n in s["primary"]]
    if ut and all(uni.get(n, {}).get("reason") == "proteomics list" for n in s["primary"]):
        return "proteomics", "universe: proteomics list (batch 10; manifest modality 'proteomics')"
    if ut and all(t == "platform" for t in ut):
        return "microarray", "universe type platform (" + ", ".join(str(uni[n].get("source")) for n in s["primary"]) + ")"
    if ut and all(t == "rnaseq_default" for t in ut):
        return "rnaseq", "universe type rnaseq_default"
    text = node.get("platform") or ""
    if node.get("modality") == "proteomics" and PROT.search(text):  # 2026-10-01: four proteomic studies were left 'unknown'
        return "proteomics", "lineage modality proteomics, platform text: " + text[:80]
    seq, arr = bool(SEQ.search(text)), bool(ARR.search(text))
    if seq != arr:
        return ("rnaseq" if seq else "microarray"), "lineage platform text: " + text[:80]
    return "unknown", "lineage platform text: " + (text[:80] or "(none)")


POOLED_GEO = os.path.join(HERE, "..", "corpus", "pooled_geo.json")
DATA_TIER = json.load(open(POOLED_GEO))["studies"] if os.path.exists(POOLED_GEO) else {}
DATA_TO_TECH = {"array_early": ("microarray", "early"), "array_whole": ("microarray", "whole"), "array_noncoding": ("microarray", "noncoding"),
                "array_unstated": ("microarray", "not stated"), "rnaseq": ("rnaseq", None), "single_cell": ("single_cell", None)}


def main():
    pl = json.load(open(os.path.join(HERE, "primary_lists.json")))
    uni = json.load(open(os.path.join(HERE, "universes.json")))["lists"]
    nodes = {n["pmid"]: n for n in json.load(open(LINEAGE))["nodes"]}
    out = {}
    for k, s in pl["studies"].items():
        if not (s.get("in_corpus", True) and s.get("primary")):
            continue
        t, ev = classify(k, s, nodes.get(k, {}), uni)
        d = data_tier(k, nodes.get(k, {}))
        if d:
            out[k] = {"label": s["label"], "technology": d[0], "evidence": d[2]}
            if d[0] == "microarray":
                out[k].update({"array_generation": d[1], "array_generation_evidence": d[2]})
            continue
        out[k] = {"label": s["label"], "technology": t, "evidence": ev}
        if t == "other":
            out[k]["assay"] = ASSAY.get(k) or "not named (add it to ASSAY)"
        if t == "microarray":
            # an override's evidence names the platform where the lineage text does not (Durrenberger 2015: GPL6255)
            text = "; ".join(x for x in (nodes.get(k, {}).get("platform"), OVERRIDES.get(k, (0, ""))[1]) if x)
            g, gev = GEN_OVERRIDES[k] if k in GEN_OVERRIDES else array_generation(text, [uni.get(n, {}) for n in s["primary"]])
            if k in GEN_OVERRIDES:
                gev = "override: " + gev
            out[k].update({"array_generation": g, "array_generation_evidence": gev})
    doc = {"generated": datetime.date.today().isoformat(), "script": "classify_technology.py",
           "inputs": {"primary_lists.json": sha(os.path.join(HERE, "primary_lists.json")),
                      "universes.json": sha(os.path.join(HERE, "universes.json")), "lineage.json": sha(LINEAGE)},
           "studies": out}
    json.dump(doc, open(os.path.join(HERE, "technology.json"), "w"), indent=1, ensure_ascii=False)
    from collections import Counter
    print(Counter(v["technology"] for v in out.values()))
    print("array generation:", Counter(v.get("array_generation") for v in out.values() if v["technology"] == "microarray"))
    for k, v in out.items():
        if v["technology"] == "unknown" or k not in ASSAY and v["technology"] == "other":
            print(v["technology"], k, v["label"], v["evidence"])


if __name__ == "__main__":
    main()
