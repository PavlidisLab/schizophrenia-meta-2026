#!/usr/bin/env python3
"""One row per person with a brain pH for the 14 dataset-regions (9 studies) of the pH panel; public GEO sample characteristics only (RA, 2026-10-19).
Selection is the one of analyses/supplement_number_tasks/taskE_ph_overlap.py (samples kept by the confounds run: present in the Gemma matrix with >= 50% non-missing values, one array per subject, dx control or schizophrenia, pH present)."""
import importlib.util, os, sys, re, json, hashlib, datetime
import pandas as pd, numpy as np
H = os.path.expanduser("~")
SRC = f"{H}/Dev/srep-mauney-integrated/analyses/supplement_number_tasks/taskE_ph_overlap.py"
ROOT = f"{H}/Dev/srep-mauney-integrated/analyses/ph_and_confound_analyses"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
txt = open(SRC).read()
ns = {"__file__": SRC}
exec(compile(txt.split("rows=[]; srcs=[]")[0], SRC, "exec"), ns)   # imports, DATASETS, gemma_cols
DATASETS, gemma_cols = ns["DATASETS"], ns["gemma_cols"]
rows, src = [], []
for name, gid, cov, reg in DATASETS:
    p = f"{ROOT}/data/covariates/{cov}.tsv"; gp = f"{ROOT}/data/gemma/{gid}_processed.tsv.gz"
    d = pd.read_csv(p, sep="\t", comment="#", dtype=str); d.columns = [c.strip() for c in d.columns]
    hdr = [l.strip() for l in open(p) if l.startswith("#")]
    if reg: d = d[d.region == reg]
    ok, inmat = gemma_cols(gid)
    d = d[d.gsm.isin(inmat) & d.dx.isin(["control", "schizophrenia"])]; d = d[d.gsm.isin(ok)]
    if d.subject.notna().any() and d.subject.duplicated().any(): d = d[~d.subject.duplicated() | d.subject.isna()]
    d = d.assign(dataset_region=name, study=name.split("_")[0] if not name.startswith("GSE92538") else "GSE92538")
    d["pH"] = pd.to_numeric(d.pH, errors="coerce"); d = d.dropna(subset=["pH"])
    d["group"] = np.where(d.dx == "schizophrenia", "schizophrenia", "control")
    for c in ("age", "PMI"): d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.rename(columns={"gsm": "geo_sample", "PMI": "pmi_hours"})
    rows.append(d[["dataset_region", "study", "geo_sample", "group", "pH", "age", "sex", "pmi_hours"]])
    srcline = next((l for l in hdr if l.lower().startswith("# source")), ""); shaline = next((l for l in hdr if "sha256" in l.lower()), "")
    src.append(f"# {name}: covariate file analyses/ph_and_confound_analyses/data/covariates/{cov}.tsv sha256 {sha(p)[:12]}; {srcline[2:].strip()}; {shaline[2:].strip()}")
S = pd.concat(rows, ignore_index=True)
def inrange(g):
    co = g[g.group == "control"].pH; ca = g[g.group == "schizophrenia"].pH
    return int(((ca >= co.min()) & (ca <= co.max())).sum()), len(ca)
tot_in = tot_n = 0
for n, g in S.groupby("dataset_region"):
    a, b = inrange(g); tot_in += a; tot_n += b
hd = ["# figure: pH distribution, cases against controls, 14 dataset-regions (9 studies); strip plots by dataset-region",
 "# script: ph_per_person_14_dataset_regions_2026-10-19/build_ph_per_person.py (reuses the selection of analyses/supplement_number_tasks/taskE_ph_overlap.py); run %s" % datetime.date.today(),
 "# rows: one per person with a pH; geo_sample = the public GEO sample accession of the array (no subject identifiers); pH as recorded in the GEO sample characteristics; age in years, pmi_hours = postmortem interval in hours, as recorded (blank where not recorded); sex as recorded. Only values from public GEO series matrices; no controlled-access data. GSE208338 shares donors with GSE21138, so persons are not independent across those two dataset-regions.",
 f"# pooled check: {tot_in} of {tot_n} schizophrenia cases ({tot_in/tot_n:.1%}) lie within the [min, max] pH of the controls of the same dataset-region (inclusive); regions {S.dataset_region.nunique()}, studies {S.study.nunique()}, persons {len(S)} ({(S.group=='schizophrenia').sum()} schizophrenia, {(S.group=='control').sum()} control)"] + src
with open("ph_per_person_14_dataset_regions.tsv", "w") as f:
    f.write("\n".join(hd) + "\n"); S.round({"pH": 3, "age": 1, "pmi_hours": 1}).to_csv(f, sep="\t", index=False)
print(tot_in, tot_n, S.dataset_region.nunique(), S.study.nunique(), len(S)); print(S.groupby(["dataset_region","group"]).size().unstack().T.to_string()); print(S[["age","sex","pmi_hours"]].notna().mean().round(2).to_dict())
