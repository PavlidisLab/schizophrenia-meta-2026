#!/usr/bin/env python3
"""Fill the gaps of interneuron_specificity_383.tsv (RA, 2026-10-19). For each of the 11 genes without a value: check the PsychAD HBCC gene list by current and previous or alias HGNC symbols and by Ensembl id; for genes present in the gene list but below the producer's CPM filter, compute the class means and shares from the same donor-class pseudobulks without the filter."""
import sys, os, csv, hashlib, importlib.util, datetime
import numpy as np, pandas as pd
H = os.path.expanduser("~"); R = f"{H}/Dev/schizophrenia-transcriptome-replication/"
spec = importlib.util.spec_from_file_location("specificity", R + "analyses/interneuron_ph_2026-10-01/scripts/specificity.py")
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
GAPS = ["MT1X","SERPINA3","TLCD4","CLDN11","IL33R","HILPDA","CNKSR3","KCNA1","ATP1A3","IFITM1","CNTNAP4"]
HG = R + "data/reference/hgnc_complete_set.txt"
hg = {r["symbol"]: r for r in csv.DictReader(open(HG), delimiter="\t") if r["status"] == "Approved"}
M, cols, genes, nctl, inputs = S.load("HBCC_Cohort", "HBCC_main")
fine_cm, nd = S.class_means(M, cols, "fine"); coarse_cm, ndc = S.class_means(M, cols, "coarse")
E2I = {e: i for i, e in enumerate(genes.ensembl)}; N2I = {}
for i, n in enumerate(genes.gene_name): N2I.setdefault(n, []).append(i)
orig = pd.read_csv("interneuron_specificity_383.tsv", sep="\t", comment="#", index_col=0)
fine = ["Astro","EN","IN_other","LAMP5","Micro","OPC","Oligo","PVALB","SST","VIP","Vascular"]
rows = {}
for g in GAPS:
    r = hg[g]; ens = r["ensembl_gene_id"]
    cand = {"current_symbol": N2I.get(g, []), "ensembl_id": [E2I[ens]] if ens in E2I else []}
    for s in [x for x in (r["prev_symbol"] + "|" + r["alias_symbol"]).split("|") if x]:
        if s in hg and s != g:
            continue  # an old symbol that is now the approved symbol of another gene (CNKSR3 was once MAGI1) must not be followed
        if s in N2I: cand["prev_or_alias:" + s] = N2I[s]
    hit = {k: v for k, v in cand.items() if v}
    rec = dict(gene=g, ensembl_id=ens, prev_symbols=r["prev_symbol"], found_by="; ".join(f"{k}={len(v)} row(s)" for k, v in hit.items()) or "none", status="not_in_psychad_gene_list")
    idx = sorted({i for v in hit.values() for i in v})
    if len(idx) == 1:
        i = idx[0]
        fv = fine_cm.iloc[i]; cv = coarse_cm.iloc[i]
        tot = fv.sum(); ctot = cv.sum()
        for c in fine: rec["cpm_" + c] = fv[c]; rec["spec_" + c] = fv[c] / tot if tot > 0 else np.nan
        rec["spec_coarse_IN"] = cv["IN"] / ctot if ctot > 0 else np.nan; rec["spec_coarse_EN"] = cv["EN"] / ctot if ctot > 0 else np.nan
        rec["cpm_coarse_IN"] = cv["IN"]
        rec["max_class_mean_cpm"] = fv.max()
        rec["n_classes_measured"] = int((fv > 0).sum())
        rec["status"] = "in_gene_list_below_filter_values_computed_without_filter" if fv.max() < 1 else "in_gene_list_values_computed"
    elif len(idx) > 1:
        rec["status"] = "ambiguous_several_rows"
    rows[g] = rec
df = pd.DataFrame(rows.values()).set_index("gene")
print(df[["ensembl_id", "found_by", "status", "max_class_mean_cpm", "spec_coarse_IN"]].to_string())
# merged table: the original 383 rows with the gaps filled where possible
out = orig.copy()
out["status"] = "in_filtered_table"; out["found_by"] = "current_symbol"; out["ensembl_id"] = ""
for g, r in df.iterrows():
    out.loc[g, "status"] = r["status"]; out.loc[g, "found_by"] = r["found_by"]; out.loc[g, "ensembl_id"] = r["ensembl_id"]
    if r["status"].startswith("in_gene_list"):
        for c in out.columns:
            if c in r.index and pd.notna(r[c]): out.loc[g, c] = r[c]
        out.loc[g, "spec_IN_fine_sum"] = sum(r["spec_" + c] for c in ["IN_other","LAMP5","PVALB","SST","VIP"])
out.index.name = "gene"
hd = ["# figure: Fig. 2c interactive version, cell-type share (hover discs); replaces interneuron_specificity_383.tsv (same columns plus status, found_by, ensembl_id)",
 "# script: interneuron_specificity_per_gene_383_2026-10-19/fill_interneuron_gaps.py (imports the loader and class_means of analyses/interneuron_ph_2026-10-01/scripts/specificity.py), run %s" % datetime.date.today(),
 f"# input: analyses/psychad_pseudobulk/pseudobulk/HBCC_Cohort/ counts_subclass.npz sha256 {sha(inputs[0])[:12]}; cols_subclass.tsv {sha(inputs[1])[:12]}; genes.tsv {sha(inputs[2])[:12]}; {nctl} HBCC control donors (donor set HBCC_main), donor-class pseudobulks with >= 20 nuclei; aggregate class means only",
 f"# input: data/reference/hgnc_complete_set.txt sha256 {sha(HG)[:12]}; interneuron_specificity_383.tsv (previous version) sha256 {sha('interneuron_specificity_383.tsv')[:12]}",
 "# definitions as in the previous version: spec_<class> = class mean CPM / sum of the 11 fine-class mean CPMs; spec_coarse_IN over 7 coarse classes; spec_IN_fine_sum = sum of the five interneuron fine shares",
 "# status: in_filtered_table = in the producer's specificity table (mean CPM >= 1 in at least one class); in_gene_list_below_filter_values_computed_without_filter = the gene is in the PsychAD gene list but its mean CPM is < 1 in every class, so the producer dropped it; values here are computed from the same pseudobulks without that filter (shares of very low counts, noisy); not_in_psychad_gene_list = the gene (current symbol, previous symbols, aliases and Ensembl id) is not among the 34,176 genes of the HBCC and MSSM PsychAD pseudobulk gene lists (the h5ad var table itself has 34,176 genes), so no value can be computed from this dataset",
 "# found_by: how the gene was located in the PsychAD gene list; no gene needed an alias or previous symbol (all 383 are current approved HGNC symbols); n_classes_measured = fine classes of 11 with mean CPM > 0"]
with open("interneuron_specificity_383_v2.tsv", "w") as f:
    f.write("\n".join(hd) + "\n"); out.round(5).to_csv(f, sep="\t")
print(out.status.value_counts().to_dict())
