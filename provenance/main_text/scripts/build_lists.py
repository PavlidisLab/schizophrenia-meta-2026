"""Group lists of human Gemma experiments with a DEA, for the revised DE priors (Paul, 2026-09-26; built by Neera).

Reads cache/datasets.json and cache/analyses.json (fetch_gemma_human.py) and the project corpus
(../../corpus/lineage.json). Writes:

    experiments.tsv                     one row per human experiment with a DEA: classes, flags, groups
    group_<name>.tsv                    gemma_id, short_name  (all technologies)
    group_<name>_no_microarray.tsv      the same without microarray datasets
    lists_meta.json                     counts, rules, inputs with sha256

Groups (Neera's numbering): 1 global, 2 brain_postmortem, 3 blood, 4 global minus (2 + 3); plus
2b brain_postmortem_no_sz_banks and brain_single_cell (single-cell kept out of 1 to 4).
Decisions (Paul, 2026-09-26): hold out whole GEO series that contain SZ; group 2 is postmortem brain only
(surgical, fetal and organoid brain go to 4); non-SZ brain datasets from the SZ cohort's banks are flagged and
listed without them as 2b. Troubled datasets are left out of every group (flagged in experiments.tsv).
"""
import csv
import hashlib
import json
import re
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DS = json.loads((HERE / "cache/datasets.json").read_text())
AN = json.loads((HERE / "cache/analyses.json").read_text())
LIN = ROOT / "corpus/lineage.json"

SZ_RE = re.compile(r"schizophreni|schizoaffective|MONDO_0005090|DOID_5419|MONDO_0005091|HP_0100753", re.I)
BRAIN_RE = re.compile(
    r"\bbrain\b|cortex|cortical|gyrus|hippocamp|amygdala|cerebell|striatum|caudate|putamen|accumbens|thalam|"
    r"substantia nigra|midbrain|\bpons\b|medulla oblongata|brain ?stem|white matter|gr[ae]y matter|frontal|temporal lobe|"
    r"parietal|occipital|cingulate|entorhinal|locus coeruleus|insula|dentate|prefrontal|olfactory bulb|brodmann|"
    r"corpus callosum|\bBA ?\d+\b|globus pallidus|nucleus basalis|claustrum|septal|choroid plexus|ventricular zone|neocortex|cerebr", re.I)
NOT_BRAIN_RE = re.compile(r"kidney|renal|adrenal|thymus|ovar|lens|hair|bone|skin|cortex of|lymph node|tooth|dental", re.I)
BLOOD_RE = re.compile(
    r"\bblood\b|peripheral blood|mononuclear|\bPBMC|leukocyte|lymphocyte|monocyte|neutrophil|granulocyte|"
    r"\bT[- ]cell|\bB[- ]cell|\bNK[- ]cell|natural killer|erythro|reticulocyte|platelet|buffy coat|\bCD[0-9]+\+", re.I)
CELL_LINE_RE = re.compile(r"cell line|lymphoblastoid|\bLCL\b|immortali[sz]ed|hek ?293|hela|jurkat|k562|thp-1|u937", re.I)
NONTISSUE_RE = re.compile(r"organoid|iPSC|induced pluripotent|stem cell[- ]derived|in vitro|cultured|culture of|xenograft|"
                          r"differentiated from|neurospheres?|spheroid|primary (human )?(astrocyte|neuron|microglia|oligodendrocyte|cell)s?|"
                          r"infection of", re.I)
NONBLOOD_CELL_RE = re.compile(r"neur|progenitor|stem cell|fibroblast|epitheli|endotheli|hepatocyte|keratinocyte|astrocyte|glia", re.I)
SURGICAL_RE = re.compile(r"surg|resect|epilep|tumou?r|glioma|glioblastoma|astrocytoma|medulloblastoma|ependymoma|"
                         r"meningioma|oligodendroglioma|biops|neurosurg|cancer|carcinoma|lobectomy", re.I)
FETAL_RE = re.compile(r"fetal|foetal|fetus|prenatal|gestation|embryonic|trimester|post-conception|\bpcw\b", re.I)
POSTMORTEM_RE = re.compile(r"post-?mortem|post mortem|autops|brain bank|tissue bank|decedent|deceased|\bPMI\b|"
                           r"agonal|donor brains?|neurobiobank|brain tissue resource", re.I)
# Banks of the project's SZ cohort (corpus lanes); text cues in the dataset name and description only.
SZ_BANKS = {
    "stanley": r"stanley",
    "libd": r"lieber|\bLIBD\b|brainseq",
    "hbcc": r"\bHBCC\b|human brain collection core",
    "mssm": r"mount sinai|\bMSSM\b|bronx va|jj peters",
    "pitt": r"pittsburgh",
    "penn": r"university of pennsylvania|\bpenn brain",
    "cmc": r"commonmind|psychencode",
    "hbtrc": r"harvard brain|\bHBTRC\b|mclean",
    "nsw": r"new south wales|\bNSW\b|victorian brain bank|australian schizophrenia research bank|\bASRB\b",
    "pritzker": r"pritzker|uc ?irvine|university of california,? irvine",
}
SC_RE = re.compile(r"single[- ]?(cell|nucle[iu]s?) (RNA|transcriptom|sequenc)|\bscRNA|\bsnRNA|single[- ]?nucle(us|i)\b|drop-?seq|sci-rna|"
                   r"10x genomics|chromium single", re.I)
OTHER_BANK_CUES = {"neurobiobank": r"neurobiobank|\bNBB\b", "maryland": r"maryland brain collection|university of maryland brain"}


def base_gse(short):
    m = re.match(r"(GSE\d+)", short or "")
    return m.group(1) if m else (short or "")


def fv_terms(analyses):
    """(category, value, uri) of every factor value used by the dataset's analyses."""
    out = []
    for a in analyses:
        fvu = a.get("factorValuesUsed") or {}
        for fvs in (fvu.values() if isinstance(fvu, dict) else [fvu]):
            for fv in fvs or []:
                for c in fv.get("characteristics") or []:
                    out.append((c.get("category") or "", c.get("value") or "", c.get("valueUri") or ""))
                if fv.get("factorValue"):
                    out.append(("", fv["factorValue"], ""))
    return out


def corpus_gses():
    L = json.loads(LIN.read_text())
    g = set()
    for n in L["nodes"] + L.get("hidden_nodes", []):
        if n.get("external"):
            continue
        g.update(base_gse(x) for x in (n.get("geo") or []))
        g.update(base_gse(x.get("short_name")) for x in (n.get("gemma") or []))
    g.update(p.name for p in (ROOT / "data/gemma").iterdir() if p.name.startswith("GSE"))
    g.discard("")
    return g


def main():
    corpus = corpus_gses()
    rows = []
    for d in DS:
        an = AN.get(str(d["id"]))
        if not an or isinstance(an, dict):
            continue
        chars = [(c.get("category") or "", c.get("value") or "", c.get("valueUri") or "") for c in d.get("characteristics") or []]
        fvs = fv_terms(an)
        text = f"{d.get('name') or ''} {d.get('description') or ''}"
        parts = [v for c, v, u in chars + fvs if c.lower() in ("organism part", "cell type", "biological sex", "")]
        org = sorted({v for c, v, u in chars + fvs if c.lower() == "organism part"})
        celltypes = sorted({v for c, v, u in chars + fvs if c.lower() == "cell type"})
        cellline = any(c.lower() == "cell line" for c, v, u in chars + fvs) or bool(CELL_LINE_RE.search(" ".join(org + celltypes)))
        tissue_txt = " ".join(org + celltypes)
        brain_parts = [o for o in org if BRAIN_RE.search(o) and not NOT_BRAIN_RE.search(o)]
        # brain only when every organism part is brain (multi-tissue atlases go to "other"); with no organism part,
        # the name must look like brain and the text must say postmortem
        non_brain = [o for o in org if o not in brain_parts and o.lower() not in ("control", "reference")]
        brain = (bool(brain_parts) and not non_brain) or \
            (not org and bool(BRAIN_RE.search(d.get("name") or "")) and bool(POSTMORTEM_RE.search(text)))
        non_blood = [o for o in org if not BLOOD_RE.search(o) and o.lower() not in ("control", "reference")] + \
            [c for c in celltypes if NONBLOOD_CELL_RE.search(c)]
        blood = bool(BLOOD_RE.search(tissue_txt)) and not brain and not non_blood  # blood-only (no multi-tissue sets)
        nontissue = bool(NONTISSUE_RE.search(text)) or cellline
        sz_hits = sorted({v for c, v, u in chars + fvs if SZ_RE.search(v) or SZ_RE.search(u)})
        if brain and not nontissue:
            if SURGICAL_RE.search(text) and not POSTMORTEM_RE.search(text):
                tissue, pm = "brain_other", "surgical or tumour"
            elif FETAL_RE.search(text):
                tissue, pm = "brain_other", "fetal or developmental"
            else:
                tissue, pm = "brain_postmortem", ("stated" if POSTMORTEM_RE.search(text) else "assumed (no surgical, tumour or fetal cue)")
        elif blood and not nontissue:
            tissue, pm = "blood", ""
        else:
            tissue, pm = ("in_vitro_or_cell_line" if nontissue else "other"), ""
        banks = sorted(k for k, r in SZ_BANKS.items() if re.search(r, text, re.I)) if tissue.startswith("brain") else []
        other_banks = sorted(k for k, r in OTHER_BANK_CUES.items() if re.search(r, text, re.I)) if tissue.startswith("brain") else []
        platforms = ";".join(p.get("shortName", "") for p in d.get("platforms") or [])
        tech = d.get("technologyType") or ""
        # RNA-seq datasets carry technologyType GENELIST in Gemma, so the assay comes from libraryStrategies
        strategies = sorted({v.get("value") or "" for v in d.get("libraryStrategies") or []} - {""})
        microarray = any(v.startswith("MICROARRAY") for v in strategies) if strategies else tech in ("ONECOLOR", "TWOCOLOR", "DUALMODE")
        # isSingleCell is false for every human dataset in the listing (2026-09-26), so text cues are used as well
        sc_basis = "isSingleCell" if d.get("isSingleCell") else ("text" if SC_RE.search(text) else "")
        rows.append({
            "gemma_id": d["id"], "short_name": d.get("shortName") or "", "geo_series": base_gse(d.get("shortName")),
            "technology_type": tech, "library_strategies": ";".join(strategies), "is_microarray": microarray,
            "is_single_cell": bool(sc_basis), "single_cell_basis": sc_basis,
            "platforms": platforms, "n_samples": d.get("numberOfBioAssays") or d.get("bioAssayCount") or "",
            "is_public": d.get("isPublic"), "troubled": bool(d.get("troubled")),
            "n_analyses": len(an), "organism_parts": "; ".join(org), "cell_types": "; ".join(celltypes[:10]),
            "tissue_class": tissue, "postmortem_basis": pm,
            "sz_terms": "; ".join(sz_hits), "in_sz_corpus": base_gse(d.get("shortName")) in corpus,
            "sz_bank_cues": ";".join(banks), "other_bank_cues": ";".join(other_banks),
            "name": (d.get("name") or "")[:200],
        })
    sz_series = {r["geo_series"] for r in rows if r["sz_terms"] or r["in_sz_corpus"]}
    for r in rows:
        why = []
        if r["in_sz_corpus"]:
            why.append("in SZ corpus")
        if r["sz_terms"]:
            why.append("SZ factor or annotation")
        if not why and r["geo_series"] in sz_series:
            why.append("same GEO series as an SZ dataset")
        r["held_out"] = bool(why)
        r["hold_out_reason"] = "; ".join(why)
        ok = not r["held_out"] and not r["troubled"]
        bulk = ok and not r["is_single_cell"]
        g = []
        if bulk:
            g.append("1_global")
            if r["tissue_class"] == "brain_postmortem":
                g.append("2_brain_postmortem")
                if not r["sz_bank_cues"]:
                    g.append("2b_brain_postmortem_no_sz_banks")
            elif r["tissue_class"] == "blood":
                g.append("3_blood")
            else:
                g.append("4_global_minus_brain_blood")
        if ok and r["is_single_cell"] and r["tissue_class"].startswith("brain"):
            g.append("brain_single_cell")
        r["groups"] = ";".join(g)

    cols = list(rows[0].keys())
    with open(HERE / "experiments.tsv", "w", newline="") as f:
        f.write(f"# Human Gemma experiments with >= 1 differential expression analysis; built {date.today()} by build_lists.py\n")
        f.write("# Source: gemma2 REST (fetch_gemma_human.py, authenticated read); corpus: corpus/lineage.json\n")
        w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    groups = ["1_global", "2_brain_postmortem", "2b_brain_postmortem_no_sz_banks", "3_blood", "4_global_minus_brain_blood", "brain_single_cell"]
    counts = {}
    for gname in groups:
        for suffix, keep in (("", lambda r: True), ("_no_microarray", lambda r: not r["is_microarray"])):
            sel = [r for r in rows if gname in r["groups"].split(";") and keep(r)]
            with open(HERE / f"group_{gname}{suffix}.tsv", "w") as f:
                f.write("gemma_id\tshort_name\n")
                f.writelines(f"{r['gemma_id']}\t{r['short_name']}\n" for r in sorted(sel, key=lambda r: r["gemma_id"]))
            counts[gname + suffix] = len(sel)
    with open(HERE / "held_out.tsv", "w") as f:
        f.write("gemma_id\tshort_name\thold_out_reason\n")
        f.writelines(f"{r['gemma_id']}\t{r['short_name']}\t{r['hold_out_reason']}\n" for r in sorted(rows, key=lambda r: r["gemma_id"]) if r["held_out"])
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    meta = {
        "built": str(date.today()), "script": "analyses/de_prior_and_recurrence_lists/build_lists.py",
        "inputs": {p: sha(HERE / p) for p in ("cache/datasets.json", "cache/analyses.json")} | {"corpus/lineage.json": sha(LIN)},
        "human_datasets_listed": len(DS), "with_dea": len(rows),
        "held_out": sum(r["held_out"] for r in rows), "troubled": sum(r["troubled"] for r in rows),
        "single_cell": sum(r["is_single_cell"] for r in rows),
        "tissue_class": {t: sum(r["tissue_class"] == t for r in rows) for t in sorted({r["tissue_class"] for r in rows})},
        "postmortem_assumed": sum(r["postmortem_basis"].startswith("assumed") for r in rows if "2_brain_postmortem" in r["groups"]),
        "group_counts": counts,
    }
    (HERE / "lists_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
