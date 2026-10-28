#!/usr/bin/env python3
"""Per-sample covariates and the dataset registry for the age/sex positive-control analysis.

    python3 scripts/build_covariates.py

Reads (never edits):
  ../confounds/data/covariates/*.tsv       GEO-derived covariates of the 14 confounds dataset-regions
                                           (../confounds/scripts/build_covariates.R)
  data/gemma/<id>_samples.json             Gemma sample annotations of the extra datasets (scripts/fetch_gemma.py)
Writes:
  data/covariates/<gemma_id>.tsv           gsm, dx, age, sex, pH, PMI, RIN, region, subject, bank for the extra datasets
  data/datasets.tsv                        one row per dataset-region: expression file, covariate file, region filter,
                                           compartment, collection cluster(s), n SZ / control with age and sex
  data/usable_gemma_ids.txt                extra Gemma ids with >= 4 SZ and >= 4 controls carrying age and sex
  data/excluded_datasets.tsv               Gemma brain SZ datasets in ../../data/gemma/_stats.json not used, and why

Collection clusters: the lineage lane of the study that generated the series (../corpus/lineage.json, field lane),
with recurrence.py's MERGE (cmc, mssm, pitt, libd one cluster). Exceptions, from project files: GSE17612 and
GSE21935 share donors (../confounds/confounds_report.md, "Donor overlap"), so both get cluster 'GSE17612+GSE21935';
GSE92538 (lineage lanes meta / brain_other only) gets 'pritzker', as the same report lists GSE80655 with GSE92538;
GSE112523 and GSE107638 list their banks per sample (Gemma 'tissuebank' / 'brainbank'); a bank named Harvard/HBTC(R)
adds cluster 'harvard' (McLean), the lane of GSE37981.
"""
import csv, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
GL = os.path.join(ROOT, "..")
DATA = os.path.normpath(os.path.join(GL, "..", "data"))  # repo data/ (lists, manifests, reference, raw)
sys.path.insert(0, GL)
import recurrence  # noqa: E402  (MERGE, sha)

CONF = os.path.join(GL, "confounds")

# name, gemma id, covariate source, region filter (value of covariate 'region'), compartment, GSE
CONFOUNDS = [
    ("GSE21138", 1880, "GSE21138", None), ("GSE35977", 5949, "GSE35977", None), ("GSE35974", 5939, "GSE35974", None),
    ("GSE53987_BA46", 8359, "GSE53987", "Pre-frontal cortex (BA46)"), ("GSE53987_hippocampus", 8359, "GSE53987", "hippocampus"),
    ("GSE53987_striatum", 8359, "GSE53987", "Associative striatum"), ("GSE17612", 1381, "GSE17612", None),
    ("GSE21935", 4150, "GSE21935", None), ("GSE92538_GPL570", 18742, "GSE92538_GPL10526", None),
    ("GSE92538_GPL96", 18741, "GSE92538_GPL17027", None), ("GSE80655_DLPFC", 13014, "GSE80655", "DLPFC"),
    ("GSE80655_AnCg", 13014, "GSE80655", "AnCg"), ("GSE80655_nAcc", 13014, "GSE80655", "nAcc"), ("GSE208338", 26880, "GSE208338", None),
]
# extra: name, gemma id, region filter, compartment
EXTRA = [
    ("GSE202537_NAc", 25749, "nucleus accumbens", "bulk"), ("GSE202537_caudate", 25749, "caudate nucleus", "bulk"),
    ("GSE202537_putamen", 25749, "putamen", "bulk"), ("GSE224683", 28105, None, "bulk"),
    ("GSE26927_temporal", 2684, "left temporal lobe", "bulk"),
    ("GSE138082_CA1", 23579, "CA1 field of hippocampus", "bulk"), ("GSE138082_CA3", 23579, "CA3 field of hippocampus", "bulk"),
    ("GSE138082_DG", 23579, "dentate gyrus of hippocampal formation", "bulk"),
    ("GSE37981", 6488, None, "LCM pyramidal neurons"),
    ("GSE12679_neuron", 946, "neuron", "LCM neurons"), ("GSE12679_endothelial", 946, "cerebral cortex endothelial cell", "LCM endothelial cells"),
    ("GSE112523", 16773, None, "sorted NeuN+ nuclei"),
    ("GSE107638_NeuN", 15479, "RBFOX3 [human] RNA binding fox-1 homolog 3  RBFOX3-positive cell", "sorted NeuN+ nuclei"),
    ("GSE107638_OLIG2", 15479, "OLIG2 [human] oligodendrocyte transcription factor 2  OLIG2-positive cell", "sorted OLIG2+ nuclei"),
]
NOT_USED = {  # gemma id: reason, from the sample annotations fetched by fetch_gemma.py --samples-only
    1205: "GSE12649 (Stanley Array Collection): no age or sex in the Gemma or GEO sample annotations",
    1204: "GSE12654 (Stanley Consortium Collection): no age or sex in the Gemma or GEO sample annotations",
    8995: "GSE62191: sex recorded for 1 of 88 samples",
    12983: "GSE78936: no age or sex in the Gemma sample annotations",
    31972: "GSE235055: every sample annotated schizophrenia (no controls) and no age or sex",
    85317: "GSE244782: schizophrenia only (clozapine vs not), no controls",
    10878: "GSE87610 (LCM L3/L5 pyramidal): no age or sex in the Gemma or GEO sample annotations",
    10879: "GSE93987 (LCM L3/L5 pyramidal): no age or sex in the Gemma or GEO sample annotations",
    51178: "CMC (brainSCOPE): single-nucleus data, per-cell-type result sets; not fetched",
    51180: "SZBDMulti-Seq (brainSCOPE): single-nucleus data, per-cell-type result sets; not fetched",
    51181: "MultiomeBrain (brainSCOPE): single-nucleus data, per-cell-type result sets; not fetched",
}


def lanes_by_gse():
    lin = json.load(open(os.path.join(GL, "..", "corpus", "lineage.json")))
    out = {}
    for n in lin["nodes"]:
        for g in n.get("geo") or []:
            if n.get("lane") and n.get("lane") != "meta":
                out.setdefault(g, set()).add(n["lane"])
    return out


def cluster_of(gse, lanes, banks=()):
    ls = {recurrence.MERGE.get(l, l) for l in lanes.get(gse, set())}
    if gse in ("GSE17612", "GSE21935"):
        ls = {"GSE17612+GSE21935"}
    if gse == "GSE92538":  # its lineage nodes are lane meta / brain_other; confounds_report.md: GSE80655 shares donors with GSE92538
        ls = {"pritzker"}
    for b in banks:
        if re.search(r"harvard|hbtc", b, re.I):
            ls.add("harvard")
    if not ls or ls == {"brain_other"}:
        ls = {gse}
    ls.discard("brain_other")
    return "|".join(sorted(ls))


def num(x):
    if x is None:
        return ""
    m = re.search(r"-?\d+(\.\d+)?", str(x))
    return m.group(0) if m and not re.search(r"\bto\b", str(x)) else ""


def sexmap(x):
    x = (x or "").strip().lower()
    return "M" if x in ("m", "male") else "F" if x in ("f", "female") else ""


def dxmap(x):
    x = (x or "").strip().lower()
    if re.search(r"schizo", x) and "affective" not in x:
        return "schizophrenia"
    if re.search(r"control|unaffected|reference subject role|healthy", x):
        return "control"
    return "other" if x else ""


DX_CATS = ("disease state", "psychiatric diagnosis", "diagnosis", "diagnostic group", "dist.dx")
AGE_CATS = ("age", "age (years)", "age at death (in years)", "age at death")
SEX_CATS = ("sex", "gender")
PMI_CATS = ("pmi", "pmi (hours)", "post-mortem delay (in hours)", "post mortem interval (pmi)", "pmi (hrs)")
PH_CATS = ("ph", "tissue ph")
SUBJ_CATS = ("patients id", "subject id")
BANK_CATS = ("tissuebank", "brainbank")


def gemma_covariates(gid):
    sj = json.load(open(os.path.join(ROOT, "data", "gemma", f"{gid}_samples.json")))["data"]
    rows = []
    for s in sj:
        ch = {}
        for c in s["sample"]["characteristics"]:
            ch.setdefault(c["category"].strip().lower(), c["value"])
        fv = {}
        for f in s["sample"]["fvbasicVOs"]:
            fv.setdefault(f["experimentalFactorCategory"]["category"], f["summary"])
        dxv = next((ch[c] for c in DX_CATS if c in ch), None)
        dx = dxmap(dxv) if dxv is not None else dxmap(fv.get("disease"))
        first = lambda cats: next((ch[c] for c in cats if c in ch), None)
        rows.append({"gsm": s["accession"]["accession"] if s.get("accession") else s["name"], "bioassay_id": s["id"],
                     "dx": dx, "age": num(first(AGE_CATS)), "sex": sexmap(first(SEX_CATS)), "pH": num(first(PH_CATS)),
                     "PMI": num(first(PMI_CATS)), "RIN": num(ch.get("rin")),
                     "region": fv.get("organism part") or fv.get("cell type") or "", "subject": first(SUBJ_CATS) or "",
                     "bank": first(BANK_CATS) or "", "gemma_outlier": "yes" if s.get("outlier") else ""})
    return rows


def read_conf_cov(name):
    with open(os.path.join(CONF, "data", "covariates", f"{name}.tsv")) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def count(rows, region):
    r = [x for x in rows if (region is None or x["region"] == region) and x["dx"] in ("control", "schizophrenia") and x["age"] != "" and x["sex"] in ("M", "F")
         and x.get("gemma_outlier", "") != "yes"]
    n = lambda dx, sx=None: sum(1 for x in r if x["dx"] == dx and (sx is None or x["sex"] == sx))
    return n("schizophrenia"), n("control"), sum(1 for x in r if x["sex"] == "M"), sum(1 for x in r if x["sex"] == "F")


def main():
    lanes = lanes_by_gse()
    reg = []
    for name, gid, cov, region in CONFOUNDS:
        rows = read_conf_cov(cov)
        gse = name.split("_")[0]
        nsz, nct, nm, nf = count(rows, region)
        reg.append({"dataset": name, "gse": gse, "gemma_id": gid, "expr": f"../confounds/data/gemma/{gid}_processed.tsv.gz",
                    "samples_json": f"../confounds/data/gemma/{gid}_samples.json", "cov": f"../confounds/data/covariates/{cov}.tsv",
                    "region": region or "", "compartment": "bulk", "cluster": cluster_of(gse, lanes), "source": "confounds",
                    "n_sz": nsz, "n_ctl": nct, "n_M": nm, "n_F": nf})
    usable = set()
    for name, gid, region, comp in EXTRA:
        rows = gemma_covariates(gid)
        fp = os.path.join(ROOT, "data", "covariates", f"{gid}.tsv")
        with open(fp, "w") as fh:
            fh.write(f"# Gemma dataset {gid}; source data/gemma/{gid}_samples.json (sha256 {recurrence.sha(os.path.join(ROOT, 'data', 'gemma', f'{gid}_samples.json'))}); "
                     "dx from the first of the sample characteristics 'disease state', 'psychiatric diagnosis', 'diagnosis', 'diagnostic group', 'dist.dx', "
                     "else the Gemma 'disease' factor value; region = Gemma 'organism part' or 'cell type' factor value\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        gse = name.split("_")[0]
        banks = sorted({r["bank"] for r in rows if r["bank"]})
        nsz, nct, nm, nf = count(rows, region)
        reg.append({"dataset": name, "gse": gse, "gemma_id": gid, "expr": f"data/gemma/{gid}_processed.tsv.gz",
                    "samples_json": f"data/gemma/{gid}_samples.json", "cov": f"data/covariates/{gid}.tsv", "region": region or "",
                    "compartment": comp, "cluster": cluster_of(gse, lanes, banks), "source": "gemma samples",
                    "n_sz": nsz, "n_ctl": nct, "n_M": nm, "n_F": nf})
        if nsz >= 4 and nct >= 4:
            usable.add(gid)
    with open(os.path.join(ROOT, "data", "datasets.tsv"), "w") as fh:
        fh.write(f"# built by scripts/build_covariates.py; lanes from ../corpus/lineage.json sha256 {recurrence.sha(os.path.join(GL, '..', 'corpus', 'lineage.json'))}; "
                 "n_* = SZ and control samples with age and sex (not Gemma outliers), before Gemma's missing-value filter and subject averaging\n")
        w = csv.DictWriter(fh, fieldnames=list(reg[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(reg)
    with open(os.path.join(ROOT, "data", "usable_gemma_ids.txt"), "w") as fh:
        fh.write("# extra Gemma ids with >= 4 SZ and >= 4 controls carrying age and sex (scripts/build_covariates.py)\n")
        fh.writelines(f"{g}\n" for g in sorted(usable))
    stats = json.load(open(os.path.join(DATA, "gemma", "_stats.json")))["stats"]
    used = {r["gemma_id"] for r in reg}
    seen = {}
    for v in stats.values():
        seen.setdefault(v["gemma_id"], v["short_name"])
    with open(os.path.join(ROOT, "data", "excluded_datasets.tsv"), "w") as fh:
        fh.write("gemma_id\tshort_name\treason\n")
        for gid, sn in sorted(seen.items()):
            if gid in used:
                continue
            fh.write(f"{gid}\t{sn}\t{NOT_USED.get(gid, 'not postmortem brain tissue by its Gemma title (iPSC-derived cells or organoids, blood, lymphocytes, fibroblasts or olfactory tissue); not fetched')}\n")
    for r in reg:
        print(r["dataset"], r["cluster"], r["compartment"], r["n_sz"], r["n_ctl"], r["n_M"], r["n_F"])


if __name__ == "__main__":
    main()
