#!/usr/bin/env python3
"""Per-gene GTEx slow-minus-fast-death axis for the 383 recurrent genes of the core 80-study base (RA, 2026-10-19).
Axis = mean over Brain-Frontal Cortex (BA9), Brain-Anterior cingulate cortex (BA24), Brain-Cortex of the per-region t (slow death, Hardy 4, minus fast natural death, Hardy 2;
OLS log2(TPM+1) ~ group + ischemic time + age + sex), exactly as ad_383.py builds 'GTEx agonal t'. Positive t = higher in slow death.
"""
import csv, gzip, hashlib, os, math, sys
import numpy as np
from scipy import stats
H = os.path.expanduser("~")
REPO = f"{H}/Dev/schizophrenia-transcriptome-replication"
HARDY = f"{REPO}/analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
REC = f"{H}/out/recurrence_genes.tsv"
CNT = f"{REPO}/figure_data/gtex_death_speed_all_tissues/donor_counts_by_tissue.tsv"
REG = ["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]
SHORT = {REG[0]: "BA9", REG[1]: "BA24", REG[2]: "Cortex"}
def sha(p):
    h = hashlib.sha256(); h.update(open(p, "rb").read()); return h.hexdigest()
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
# donors per region in the slow (4) vs fast (2) contrast (donor_counts_by_tissue.tsv; one sample per donor per region); df = n - 5 (intercept, group, ischemic time, age, sex)
n4 = {}; n2 = {}
for r in csv.DictReader((l for l in open(CNT) if not l.startswith("#")), delimiter="\t"):
    if r["tissue"] in REG: n4[r["tissue"]] = int(r["n_hardy_4"]); n2[r["tissue"]] = int(r["n_hardy_2"])
rows = list(csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t"))
rec = {r["gene"]: r for r in rows if r["scenario"] == "published, brain"}
S383 = [g for g, r in rec.items() if float(r["q_bh"]) < 0.1]; assert len(S383) == 383
T = {}
with gzip.open(HARDY, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")[1:]
    ix = [hdr.index(r) for r in REG]
    for l in f:
        s = l.rstrip("\n").split("\t"); v = []
        for i in ix:
            try: v.append(float(s[1 + i]))
            except (ValueError, IndexError): v.append(float("nan"))
        T[s[0]] = v
# axis-wide q: intersection-union p = largest regional two-sided p (valid whatever the correlation between regions, which share donors); needs the same sign in all three
def pvals(v):
    return [2 * stats.t.sf(abs(t), n4[r] + n2[r] - 5) if not math.isnan(t) else float("nan") for t, r in zip(v, REG)]
allg = [g for g, v in T.items() if not any(math.isnan(x) for x in v)]
pa = {}
for g in allg:
    v = T[g]; p = pvals(v)
    pa[g] = max(p) if (all(x > 0 for x in v) or all(x < 0 for x in v)) else 1.0
qa = dict(zip(allg, bh([pa[g] for g in allg])))
out = []
for g in S383:
    v = T.get(g, [float("nan")] * 3); sg = int(rec[g]["up"]) - int(rec[g]["down"]); sz = "up" if sg > 0 else "down" if sg < 0 else "none"
    ok = not any(math.isnan(x) for x in v)
    ax = float(np.mean(v)) if ok else float("nan")
    cons = ("TRUE" if (np.sign(ax) == np.sign(sg)) else "FALSE") if ok and sg != 0 else ""
    out.append([g, sz, rec[g]["up"], rec[g]["down"]] + [("%.4f" % x if not math.isnan(x) else "") for x in v] + [("%.4f" % ax) if ok else "", ("%.3g" % qa[g]) if ok else "", cons, "3" if ok else str(3 - sum(math.isnan(x) for x in v))])
recs = f"{REC}"
hdrs = [
 "# figure: Fig. 2c agonal strip (proposed); per-gene companion of fig07_a2_a4_gene_scatter.tsv (ph_confound_consistent, q_ph_controls)",
 "# script: agonal_axis_per_gene_383_2026-10-19/build_agonal_axis_per_gene_383.py, run 2026-10-19",
 f"# input: analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz sha256 {sha(HARDY)} (GTEx v8 brain, per-region t of slow death (Hardy 4) minus fast natural death (Hardy 2), OLS log2(TPM+1) ~ group + ischemic time + age + sex, genes with median TPM >= 1 in the region)",
 f"# input: recurrence_on_80_studies/out/recurrence_genes.tsv sha256 {sha(REC)} (scenario 'published, brain', q_bh < 0.1 = 383 genes; core 80-study base)",
 f"# input: figure_data/gtex_death_speed_all_tissues/donor_counts_by_tissue.tsv sha256 {sha(CNT)} (donors per Hardy group, used for degrees of freedom)",
 "# sign convention: t > 0 = higher expression in slow death than in fast natural death. agonal_t_axis = mean of the t of BA9, BA24 and Cortex (as in ad_383.py 'GTEx agonal t'); blank where the gene is not measured (median TPM < 1) in any of the three.",
 "# sz_direction = up if n_up > n_down lists, down if n_down > n_up, none if equal (same sign rule as ad_383.py / sz_sign).",
 "# agonal_consistent = TRUE when sign(agonal_t_axis) == sign of the SZ direction, i.e. the schizophrenia change goes the way slow death goes (the analogue of ph_confound_consistent, which is sign(Z_pH_controls) == -sz_sign, i.e. SZ goes the way low pH goes). Blank where no SZ direction or axis.",
 "# q_agonal = BH q across all genes measured in the three regions, from an intersection-union p = the largest of the three two-sided regional p (t distribution, df = slow + fast donors - 5), set to 1 when the three t differ in sign. It is conservative and does not assume independence of the regions, which share donors.",
 f"# donors (slow/fast) behind the t: " + "; ".join(f"{SHORT[r]} {n4[r]}/{n2[r]}" for r in REG) + ". Donor counts are all samples with a Hardy score in the region; samples lacking an ischemic time (dropped by the producer) are not subtracted, so df can be up to a few too large.",
 "# n_regions = regions (of 3) in which the gene was measured. No private data (GTEx open-access summary t only)."]
with open("agonal_axis_per_gene_383.tsv", "w") as f:
    f.write("\n".join(hdrs) + "\n")
    f.write("\t".join(["gene", "sz_direction", "n_up_lists", "n_down_lists", "t_BA9", "t_BA24", "t_Cortex", "agonal_t_axis", "q_agonal", "agonal_consistent", "n_regions"]) + "\n")
    for o in out: f.write("\t".join(o) + "\n")
m = [o for o in out if o[7] != ""]; c = [o for o in out if o[9] == "TRUE"]; d = [o for o in out if o[9] == "FALSE"]
print("383 genes; with axis", len(m), "; consistent", len(c), "; inconsistent", len(d), "; q<0.05", sum(float(o[8]) < 0.05 for o in m), "; q<0.1", sum(float(o[8]) < 0.1 for o in m))
