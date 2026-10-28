#!/usr/bin/env python3
"""Per-gene PsychAD HBCC control cell-class specificity for the 383 recurrent genes of the core 80-study base (RA, 2026-10-19). Aggregate class means only."""
import csv, gzip, hashlib, os, datetime
import pandas as pd
H = os.path.expanduser("~"); REPO = f"{H}/Dev/schizophrenia-transcriptome-replication"
SPEC = f"{REPO}/analyses/interneuron_ph_2026-10-01/results/specificity_HBCC.tsv.gz"
REC = f"{H}/out/recurrence_genes.tsv"
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
rows = [r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"]
S383 = [r["gene"] for r in rows if float(r["q_bh"]) < 0.1]; assert len(S383) == 383
sp = pd.read_csv(SPEC, sep="\t", comment="#", index_col=0)
fine = ["Astro", "EN", "IN_other", "LAMP5", "Micro", "OPC", "Oligo", "PVALB", "SST", "VIP", "Vascular"]
cols = ["spec_coarse_IN", "spec_coarse_EN"] + [f"spec_{c}" for c in fine] + [f"cpm_{c}" for c in fine] + ["cpm_coarse_IN"]
d = sp.reindex(S383)[cols]
d["n_classes_measured"] = (sp.reindex(S383)[[f"cpm_{c}" for c in fine]] > 0).sum(axis=1).where(d["spec_coarse_IN"].notna())
d["spec_IN_fine_sum"] = d[["spec_IN_other", "spec_LAMP5", "spec_PVALB", "spec_SST", "spec_VIP"]].sum(axis=1, skipna=False)
d.index.name = "gene"
hd = ["# figure: Fig. 2c interneuron strip (proposed); companion of fig02_b_gene_agonal_axis.tsv",
 "# script: interneuron_specificity_per_gene_383_2026-10-19/build_interneuron_specificity_383.py, run %s" % datetime.date.today(),
 f"# input: analyses/interneuron_ph_2026-10-01/results/specificity_HBCC.tsv.gz sha256 {sha(SPEC)} (made by analyses/interneuron_ph_2026-10-01/scripts/specificity.py, run 2026-10-01 from PsychAD HBCC_Cohort control donor-class pseudobulks, 95 control donors, donor-class pseudobulks with >= 20 nuclei; aggregate class means only)",
 f"# input: recurrence_on_80_studies/out/recurrence_genes.tsv sha256 {sha(REC)} (scenario 'published, brain', q_bh < 0.1 = 383 genes)",
 "# definition: spec_<class> = mean CPM in the class / sum of mean CPM over the 11 fine classes (EWCE specificity). Fine classes: Astro, EN (all excitatory), IN_other (ADARB2), LAMP5, Micro, OPC, Oligo, PVALB, SST, VIP, Vascular (Endo+PC+SMC+VLMC); PVM and Adaptive left out. Shares sum to 1 over the 11.",
 "# spec_coarse_IN / spec_coarse_EN = the same share with all interneuron subclasses (7) pooled into one class IN (and EN pooled), over the 7 coarse classes Astro, EN, IN, Micro, OPC, Oligo, Vascular, from the file's own columns; cpm_coarse_IN = its mean CPM. spec_IN_fine_sum = sum of the five fine interneuron shares (over 11 classes), for comparison.",
 "# n_classes_measured = fine classes (of 11) with mean CPM > 0; blank where the gene is not in the file (removed by the producer's expression filter).",
 "# No private data."]
with open("interneuron_specificity_383.tsv", "w") as f:
    f.write("\n".join(hd) + "\n"); d.round(5).to_csv(f, sep="\t")
print(d.spec_coarse_IN.notna().sum(), "of 383 in file;", "median spec_coarse_IN", d.spec_coarse_IN.median(), "; IN_fine_sum vs coarse corr", d.spec_IN_fine_sum.corr(d.spec_coarse_IN))
