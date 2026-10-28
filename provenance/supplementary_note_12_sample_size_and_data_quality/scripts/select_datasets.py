#!/usr/bin/env python3
"""Choose external human brain datasets (non-psychiatric studies) from Gemma for the age comparison. READ ONLY (GET).

    opt/anaconda3/bin/python3 scripts/select_datasets.py         (run from controls_age_external/)

Host: ../gemma_host.py (keychain GEMMA_BASE_URL, currently https://gemma2.msl.ubc.ca; refuses localhost). GET requests
only, authenticated when credentials exist, through controls_age_sex/scripts/fetch_gemma.py get() (imported, not copied).

Selection, in this order (every screened dataset is a row of data/candidates.tsv with the step that removed it):
  S1  universe: Gemma human datasets annotated to brain (filter allCharacteristics.valueUri = UBERON_0000955; Gemma
      includes datasets annotated to its subclasses), plus Gemma human datasets with a sample characteristic 'age'
      whose name or description matches BRAIN_TEXT (the brain annotation misses some, e.g. GSE53890). Troubled datasets
      are dropped.
  S2  title screen (dataset name only; the description mentions other disorders in passing): psychiatric or substance
      studies (PSYCH), in-vitro, tumour, fetal or developmental, single-cell and non-tissue studies (NOT_TISSUE).
  S3  sample annotations (GET /datasets/{id}/samples; cached in data/cache/samples/): per sample age, sex, diagnosis,
      region, subject, bank (parse_samples). A dataset any of whose samples carries a psychiatric diagnosis is removed.
  S4  per dataset x region: non-disease samples (dx_class control) with age >= 18 (AGE_MIN) and sex recorded (the age model
      is ~ age + sex), not surgical (SURGICAL: resected
      or biopsy tissue, or an epilepsy / surgery annotation). A region is eligible when
      it has >= 20 such subjects (MIN_N) and their ages span >= 40 years (MIN_SPAN).
  S5  cortex first: eligible regions matching CORTEX; one per series (PICK: frontal before other cortex, then most
      subjects). These are the primary external datasets. Eligible non-cortical regions are listed, not used.
Subject separation from the SZ-study datasets is checked afterwards (scripts/subject_check.py).

Writes data/candidates.tsv (one row per dataset in S1), data/candidate_regions.tsv (one row per dataset x region reached
S4), data/selected.tsv (the S5 picks), data/cache/listing_*.json with their _meta.json.
"""
import csv, datetime, hashlib, json, os, re, sys, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
GL = os.path.join(ROOT, "..")
sys.path.insert(0, GL)
sys.path.insert(0, os.path.join(GL, "controls_age_sex", "scripts"))
import gemma_host  # noqa: E402
import fetch_gemma as FG  # noqa: E402  (get: GET with retries on 502)

CACHE = os.path.join(ROOT, "data", "cache")
TODAY = datetime.date.today().isoformat()
AGE_MIN, MIN_N, MIN_SPAN = 18, 20, 40
BRAIN_URI = "http://purl.obolibrary.org/obo/UBERON_0000955"
BRAIN_TEXT = re.compile(r"brain|cortex|cortical|hippocamp|gyrus|cerebr|frontal|temporal lobe|entorhinal|striatum|caudate|putamen|"
                        r"substantia nigra|cerebell|amygdala|thalam|Brodmann|white matter|gr[ae]y matter|post-?mortem", re.I)
PSYCH = re.compile(r"schizo|bipolar|depress|autis|\bASD\b|alcohol|cocaine|heroin|opioid|suicid|psychiatr|PTSD|post-?traumatic|"
                   r"obsessive|addict|substance|psychos|psychotic|antidepress|antipsychotic|smok|drug abuse|anxiety|mood", re.I)
NOT_TISSUE = re.compile(r"organoid|iPSC|iPS\b|pluripotent|stem cell|progenitor|cell line|cultur|in vitro|xenograft|hESC|"
                        r"glioma|glioblastoma|meningioma|medulloblastoma|astrocytoma|neurocytoma|tumou?r|neoplas|metasta|cancer|"
                        r"carcinoma|fetal|foetal|fetus|embryo|prenatal|gestation|trimester|developing|development|"
                        r"single[- ]cell|single[- ]nucle|snRNA|scRNA|nuclei|FANS|laser[- ]captur|microdissect|blood|leukocyte|"
                        r"lymphocyte|fibroblast|monocyte|adrenal|muscle|retina|olfactory|islet|pancrea|placenta|chimpanzee|"
                        r"macaque|mouse|microglia|astrocytes\b|neurons from|sorted|exosom|artery|endothelial|cell[- ]type[- ]specific|pyramidal neurons", re.I)
# region of datasets whose samples carry no region annotation and no single BioSource naming one, from the GEO series text
REGION_OVERRIDE = {"GSE53890": ("frontal cortex", "GEO GSE53890 summary: 'adult human brain samples from frontal cortical regions'; overall design: "
                                "'12 young (<40yr), 9 middle aged (40-70yr), 16 normal aged (70-94yr), and 4 extremely aged (95-106yr)'")}
CORTEX = re.compile(r"cortex|cortical|gyrus|\bBA ?\d|Brodmann|frontal|temporal|parietal|occipital|cingulate|insula|FCTX|TCTX|OCTX|lobe", re.I)
NOT_CORTEX = re.compile(r"hippocamp|entorhinal|white matter|cerebell|dentate|subiculum", re.I)
FRONTAL = re.compile(r"frontal|prefrontal|\bBA ?(8|9|10|11|44|45|46|47)\b|FCTX|DLPFC|orbito", re.I)
AGE_CAT = re.compile(r"^(donor |expired |)age( at death| of death|_at_death|_death|atdeath)?( \(.*\)| \[.*\])?$", re.I)
SEX_CAT = re.compile(r"^(sex|gender|biological sex|donor sex)$", re.I)
DX_CAT = re.compile(r"^(disease|disease state|disease status|diagnosis|dx|condition|clinical diagnosis|group|phenotype|"
                    r"neuropathological diagnosis|pathological diagnosis|disease group|subject group|cohort|status|clinical status)$", re.I)
REGION_CAT = re.compile(r"^(organism part|tissue|brain region|region|brain area|tissue region|anatomical region|cortical region)$", re.I)
SUBJ_CAT = re.compile(r"^(individual|individual id|subject|subject id|donor|donor id|donor_id|patient|patient id|brain id|brain code|case id|individual identifier)$", re.I)
SURGICAL = re.compile(r"resect|surgical|surgery|biops|neurosurg", re.I)  # sample-level: tissue not taken post mortem
BANK_CAT = re.compile(r"^(tissuebank|tissue bank|brainbank|brain bank|bank|source bank)$", re.I)
PMI_CAT = re.compile(r"^(pmi|pmi \(.*\)|post-?mortem interval.*|postmortem interval.*|post mortem delay.*|post-mortem delay.*)$", re.I)
PH_CAT = re.compile(r"^(ph|brain ph|tissue ph)$", re.I)
RIN_CAT = re.compile(r"^(rin|rna integrity number)$", re.I)
CONTROL = re.compile(r"^\s*(control|controls|normal|healthy|healthy control|unaffected|non-?demented|non-?demented control|nd|ndc|ctrl|ctl|"
                     r"neurologically normal|reference subject role|no disease|none|nci|no cognitive impairment|normal control|"
                     r"non-?disease|non-?neurological control|non-?ad|con|c|neurologically healthy|neurotypical|"
                     r"age-matched control|normal aging|normal ageing|not applicable|n/?a)\s*$", re.I)


def gemma_get(path, **params):
    url = gemma_host.rest_base() + path + ("?" + urllib.parse.urlencode(params) if params else "")
    body = FG.get(url, "application/json")
    return url, body, json.loads(body)


def listing(name, filt):
    """All datasets matching a Gemma filter; cached as data/cache/listing_<name>.json with a _meta.json."""
    fp = os.path.join(CACHE, f"listing_{name}.json")
    if os.path.exists(fp):
        return json.load(open(fp))
    out, off, meta = [], 0, {"filter": filt, "pages": []}
    while True:
        url, body, d = gemma_get("/datasets", filter=filt, limit=100, offset=off, sort="+id")
        out += d["data"]
        meta["pages"].append({"url": url, "sha256_body": hashlib.sha256(body).hexdigest(),
                              "retrieved": datetime.datetime.now(datetime.timezone.utc).isoformat()})
        off += 100
        if off >= d["totalElements"]:
            break
    json.dump(out, open(fp, "w"))
    json.dump(meta, open(fp[:-5] + "_meta.json", "w"), indent=1)
    return out


def samples(gid):
    fp = os.path.join(CACHE, "samples", f"{gid}.json")
    if os.path.exists(fp):
        return json.load(open(fp)), json.load(open(fp[:-5] + "_meta.json"))
    url, body, d = gemma_get(f"/datasets/{gid}/samples")
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    open(fp, "wb").write(body)
    meta = {"url": url, "sha256_body": hashlib.sha256(body).hexdigest(), "bytes": len(body),
            "retrieved": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    json.dump(meta, open(fp[:-5] + "_meta.json", "w"), indent=1)
    return d, meta


def age_years(v):
    """Age in years from a Gemma characteristic value; None when absent or a range. Months, weeks, days and
    post-conception or gestational units are converted (they only ever put a sample below AGE_MIN)."""
    if v is None:
        return None
    s = str(v).strip().lower()
    if re.search(r"\d\s*(-|to)\s*\d", s) or not re.search(r"\d", s):
        return None
    m = re.search(r"-?\d+(\.\d+)?", s)
    x = float(m.group(0))
    if re.search(r"pcw|gw|gestation|post-?concep|fetal|embryo", s):
        return -1.0
    if re.search(r"month|\bmo\b|mos\b", s):
        return x / 12
    if re.search(r"week|wk", s):
        return x / 52
    if re.search(r"day", s):
        return x / 365
    return x


def name_age(s):
    """GSE1572-style sample names ('70 year old male')."""
    m = re.search(r"(\d+(\.\d+)?)\s*(-\s*)?(yo\b|y\.o\.|years?\b|yrs?\b)", s or "", re.I)
    return float(m.group(1)) if m else None


def parse_samples(sj, short_name=""):
    """One record per sample: gsm, name, age, sex, dx_values, dx_class (control / disease / none), psych, region,
    subject, bank, pmi, ph, rin, outlier. dx_class 'none' = no diagnosis-like annotation on the sample."""
    out = []
    for s in sj:
        ch = {}
        for c in s["sample"].get("characteristics") or []:
            ch.setdefault((c.get("category") or "").strip(), c.get("value"))
        fv = {}
        for f in s["sample"].get("fvbasicVOs") or []:
            cat = (f.get("experimentalFactorCategory") or {}).get("category") or ""
            fv.setdefault(cat.strip(), f.get("summary") or f.get("value"))
        pick = lambda rx, d: next((v for k, v in d.items() if rx.match(k)), None)
        age = next((age_years(v) for k, v in ch.items() if AGE_CAT.match(k) and age_years(v) is not None), None)
        if age is None:
            age = age_years(fv.get("age"))
        if age is None:
            age = name_age(s["name"]) or name_age(s["sample"].get("description")) or name_age(ch.get("BioSource"))
        sex = (pick(SEX_CAT, ch) or fv.get("biological sex") or fv.get("sex") or "").strip().lower()
        sex = "M" if sex in ("m", "male") else "F" if sex in ("f", "female") else ""
        if not sex:  # GSE1572-style sample names ('70 year old male')
            m = re.search(r"\b(female|male)\b", s["name"], re.I)
            sex = {"male": "M", "female": "F"}[m.group(1).lower()] if m else ""
        dxv = [str(v) for k, v in ch.items() if DX_CAT.match(k) and v is not None] + [str(fv[k]) for k in ("disease", "disease staging") if fv.get(k)]
        if not dxv:
            dxc = "none"
        elif all(CONTROL.match(v) for v in dxv):
            dxc = "control"
        else:
            dxc = "disease"
        allv = " ".join(str(v) for v in list(ch.values()) + list(fv.values()) if v)
        region = fv.get("organism part") or fv.get("cell type") or pick(REGION_CAT, ch) or ""
        if not region and short_name.split(".")[0] in REGION_OVERRIDE:
            region = REGION_OVERRIDE[short_name.split(".")[0]][0]
        if not region:  # fall back to the BioSource text (two-colour arrays also carry the reference RNA's BioSource)
            bs = [c.get("value") for c in s["sample"].get("characteristics") or [] if c.get("category") == "BioSource"
                  and not re.search(r"reference", c.get("value") or "", re.I)]
            region = bs[0] if bs else ""
        surg = bool(SURGICAL.search(allv)) or any(re.search(r"epilep|surgery", k, re.I) for k in ch)
        out.append({"gsm": s["accession"]["accession"] if s.get("accession") else s["name"], "name": s["name"], "surgical": surg,
                    "bioassay_id": s["id"], "age": age, "sex": sex, "dx_values": "|".join(dxv), "dx_class": dxc,
                    "psych": bool(PSYCH.search(" ".join(dxv))), "region": str(region).strip(),
                    "subject": str(pick(SUBJ_CAT, ch) or ""), "bank": str(pick(BANK_CAT, ch) or ""),
                    "pmi": str(pick(PMI_CAT, ch) or ""), "ph": str(pick(PH_CAT, ch) or ""), "rin": str(pick(RIN_CAT, ch) or ""),
                    "outlier": bool(s.get("outlier")), "all_values": allv})
    return out


def region_summary(recs):
    """per region: adult non-disease samples (dx_class control, or none when no sample of the dataset has a diagnosis)."""
    any_dx = any(r["dx_class"] != "none" for r in recs)
    by = {}
    for r in recs:
        by.setdefault(r["region"], []).append(r)
    rows = []
    for reg, rs in by.items():
        ok = [r for r in rs if not r["outlier"] and not r["surgical"] and r["sex"] in ("M", "F") and r["age"] is not None and r["age"] >= AGE_MIN
              and (r["dx_class"] == "control" or (r["dx_class"] == "none" and not any_dx))]
        subj = {r["subject"] or r["gsm"] for r in ok}
        ages = [r["age"] for r in ok]
        rows.append({"region": reg, "n_samples": len(rs), "n_age": sum(r["age"] is not None for r in rs),
                     "n_adult_control_samples": len(ok), "n_adult_control_subjects": len(subj),
                     "age_min": min(ages) if ages else "", "age_max": max(ages) if ages else "",
                     "n_male": sum(r["sex"] == "M" for r in ok), "n_female": sum(r["sex"] == "F" for r in ok),
                     "dx_values_used": "|".join(sorted({r["dx_values"] or "(none)" for r in ok}))[:300],
                     "banks": "|".join(sorted({r["bank"] for r in ok if r["bank"]}))[:200]})
    return rows


def write(path, rows, head):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})


def main():
    os.makedirs(CACHE, exist_ok=True)
    A = listing("brain_uri", f"taxon.commonName = human and allCharacteristics.valueUri = {BRAIN_URI}")
    Bl = listing("age_category", "taxon.commonName = human and bioAssays.sampleUsed.characteristics.category = age")
    uni = {}
    for x in A:
        uni[x["id"]] = (x, "brain annotation")
    for x in Bl:
        if x["id"] not in uni and BRAIN_TEXT.search(x["name"] + " " + (x.get("description") or "")):
            uni[x["id"]] = (x, "age characteristic + brain text")
    cand, regions = [], []
    for gid in sorted(uni):
        x, how = uni[gid]
        row = {"gemma_id": gid, "short_name": x["shortName"], "n_bioassays": x.get("numberOfBioAssays"), "name": x["name"][:160],
               "universe": how, "step_removed": "", "reason": ""}
        if x.get("troubled"):
            row.update(step_removed="S1", reason="troubled in Gemma")
        elif PSYCH.search(x["name"]):
            row.update(step_removed="S2", reason=f"title: psychiatric or substance study ({PSYCH.search(x['name']).group(0)})")
        elif NOT_TISSUE.search(x["name"]):
            row.update(step_removed="S2", reason=f"title: not adult bulk brain tissue ({NOT_TISSUE.search(x['name']).group(0)})")
        cand.append(row)
        if row["step_removed"]:
            continue
        sj, meta = samples(gid)
        recs = parse_samples(sj["data"], x["shortName"])
        row["samples_sha256"] = meta["sha256_body"][:12]
        psych = [r for r in recs if r["psych"]]
        if psych:
            row.update(step_removed="S3", reason=f"{len(psych)} samples with a psychiatric diagnosis ({psych[0]['dx_values'][:60]})")
            continue
        if not any(r["age"] is not None for r in recs):
            row.update(step_removed="S3", reason="no age parsed from the sample annotations")
            continue
        rs = region_summary(recs)
        best = None
        for r in rs:
            r.update({"gemma_id": gid, "short_name": x["shortName"]})
            span = (r["age_max"] - r["age_min"]) if r["n_adult_control_subjects"] else 0
            r["eligible"] = r["n_adult_control_subjects"] >= MIN_N and span >= MIN_SPAN
            r["cortex"] = bool(CORTEX.search(r["region"])) and not NOT_CORTEX.search(r["region"])
            r["frontal"] = bool(FRONTAL.search(r["region"]))
            regions.append(r)
            if r["eligible"] and r["cortex"]:
                if best is None or (r["frontal"], r["n_adult_control_subjects"]) > (best["frontal"], best["n_adult_control_subjects"]):
                    best = r
        if best is None:
            el = [r for r in rs if r["eligible"]]
            row.update(step_removed="S4" if not el else "S5",
                       reason=(f"no region with >= {MIN_N} non-disease subjects aged >= {AGE_MIN} with sex recorded, ages spanning >= {MIN_SPAN} years "
                               f"(best: {max((r['n_adult_control_subjects'] for r in rs), default=0)} subjects)") if not el else
                       f"eligible regions not cortex: {', '.join(r['region'] for r in el)}")
            continue
        best["picked"] = True
        row.update(reason=f"picked region '{best['region']}': {best['n_adult_control_subjects']} subjects aged {best['age_min']:.0f} to {best['age_max']:.0f}")
    head = [f"# scripts/select_datasets.py run {TODAY}; Gemma host {gemma_host.host()[0]} ({gemma_host.host()[1]}); GET only",
            f"# S1 universe: human datasets annotated to brain ({BRAIN_URI}, Gemma term inference) plus human datasets with sample characteristic "
            "'age' whose name or description names a brain part; S2 title screen (psychiatric / substance; not adult bulk brain tissue); "
            "S3 sample annotations (a psychiatric diagnosis on any sample; no age); "
            f"S4 a region with >= {MIN_N} non-disease, non-surgical subjects aged >= {AGE_MIN} with sex recorded, ages spanning >= {MIN_SPAN} years; S5 cortex, one region per series",
            "# step_removed blank = selected; the regexes are in the script"]
    write(os.path.join(ROOT, "data", "candidates.tsv"), cand, head)
    write(os.path.join(ROOT, "data", "candidate_regions.tsv"), regions,
          [f"# scripts/select_datasets.py run {TODAY}; one row per dataset x region of datasets that reached S4; non-disease = dx_class control "
           "(every diagnosis-like annotation reads as control) or no diagnosis annotation in the whole dataset; outliers flagged by Gemma left out"])
    sel = [r for r in regions if r.get("picked")]
    write(os.path.join(ROOT, "data", "selected.tsv"), sel, [f"# scripts/select_datasets.py run {TODAY}; S5 picks (one cortical region per series)"])
    from collections import Counter
    print(Counter(r["step_removed"] or "selected" for r in cand))
    for r in sel:
        print(r["gemma_id"], r["short_name"], r["region"], r["n_adult_control_subjects"], r["age_min"], r["age_max"], r["dx_values_used"][:80], r["banks"][:60])


if __name__ == "__main__":
    main()
