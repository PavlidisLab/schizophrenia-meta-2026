#!/usr/bin/env python3
"""Supplementary Table 1 on the 516-study corpus (Mauney 2015 added). Adapted from code001_nk_correction_2026-10-11/supp_table1/export_supp_table1_final.py:
recurrence columns come from the in-place analyses/recurrence_genes.tsv (RECURRENCE_NK=membership, 64 studies); in_analysis_set is read from
confounds/recurrent_sets.tsv (set study_brain, the 209 recurrent genes of the freeze, not rerun); rows are the genes with q_uniform < 0.1 now plus any of the 209 that are not.
Original description: final-state version: same rows, values and sort as export_supp_table1.py (after_membership run),
with in_analysis_set (yes/no) in place of in_176, no new_vs_176 column, and header lines that describe the final definitions only.
in_analysis_set = gene has q_bh < 0.1 in the 209-gene recurrence file analyses/recurrence_genes.tsv ('published, brain'), the set used in all analyses.
Run: python3 export_supp_table1_final.py   (writes figure_data/supplementary_tables/supp_table1_*.tsv and supp_table1_meta.json)"""
import hashlib, json, os, gzip, io
from pathlib import Path
import numpy as np, pandas as pd
A = Path(__file__).resolve().parents[1]; R = A.parent; OUT = Path(os.environ["OUT"]); OUT.mkdir(parents=True, exist_ok=True); AXV2 = Path(os.environ["AXIS_V2"])
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
G = A / "recurrence_genes.tsv"; SET = A / "confounds/recurrent_sets.tsv"
P = A / "confounds/results/combined_controls_pH.tsv"; S = R / "figure_data/fig07/fig07_a2_a4_gene_scatter.tsv"
rd = lambda f: (lambda d: d[d.scenario == "published, brain"])(pd.read_csv(f, sep="\t", comment="#"))
g = rd(G).set_index("gene"); aset = set(pd.read_csv(SET, sep="\t", comment="#").query("set == 'study_brain'").gene); assert len(aset) == 383, len(aset)
p = pd.read_csv(P, sep="\t", comment="#"); p = p[p.factor == "pH"].drop_duplicates("gene").set_index("gene")
s = pd.read_csv(S, sep="\t", comment="#").set_index("gene")
z2 = pd.read_csv(AXV2, sep="\t", comment="#").drop_duplicates("gene").set_index("gene").z_stouffer
def row(gn):
    r = g.loc[gn] if gn in g.index else None
    f = lambda v, c=lambda x: x: None if r is None else c(r[v])
    return dict(gene=gn, lists_naming=f("n_lists_observed", int), lists_measured=f("n_lists_measured", int), n_up=f("up", int), n_down=f("down", int),
        q_uniform=f("q_bh"), q_expression_matched=f("q_expr_bh"), q_de_prior=f("q_prior_bh"),
        n_clusters_observed=f("n_clusters_observed", int), q_cluster_bh=f("q_cluster_bh"),
        in_analysis_set="yes" if gn in aset else "no",
        ph_z_controls=z2.get(gn, np.nan), sz_z_raw=s.z_sz_raw.get(gn, np.nan), sz_z_ph_adjusted=s.z_sz_ph_adjusted.get(gn, np.nan))
main = pd.DataFrame([row(x) for x in sorted(set(g[g.q_bh < 0.1].index) | aset, key=lambda k: (-g.loc[k].n_lists_observed, k))])
out_set = main[main.in_analysis_set == "no"]; out_txt = ", ".join(f"{r.gene} ({r.q_uniform:.4f})" for r in out_set.itertuples()); lost = sorted(aset - set(g[g.q_bh < 0.1].index)); lost_txt = ", ".join(f"{x} ({g.loc[x].q_bh:.4f})" if x in g.index else x for x in lost)
assert set(g[g.q_bh < 0.1].index) == aset, "q<0.1 rows differ from the analysis set"
prov = ["# figure: supplementary gene table",
 "# script: notes_tables_v2_2026-10-15/L02/g1/export_supp_table1_383.py (edit of analyses/mauney_integrated_2026-10-12/export_supp_table1_integrated.py), run 2026-10-16",
 f"# input: analyses/recurrence_genes.tsv sha256 {sha(G)[:12]} (remote environment repo, final run; scenario 'published, brain', core 80-study base, 72 lists contributing)",
 f"# input: analyses/confounds/recurrent_sets.tsv sha256 {sha(SET)[:12]} (set study_brain: 383 genes at q_bh < 0.1)",
 f"# input: inputs_ph_axis_v2/axis_v2_controls.tsv.gz sha256 {sha(AXV2)[:12]} (pH axis, controls, column z_stouffer over 10 datasets, genes in at least 8; column ph_z_controls)",
 f"# input: figure_data/fig07/fig07_a2_a4_gene_scatter.tsv sha256 {sha(S)[:12]} (L02 rebuild (export_fig02_fig07.py) from recurrent_genes_pH.tsv; SZ z raw and after pH adjustment, five datasets, 5+ list set only)",
 f"# rows: the {len(main)} genes with q_uniform < 0.1 (= the 383 of recurrent_sets.tsv set study_brain); q columns are BH within the scenario",
 "# q_cluster_bh: BH q of the cohort-tier Poisson-binomial, uniform null, same BH family as q_uniform; n_clusters_observed = clusters with a member list naming the gene"]
nsz = int(main.sz_z_raw.notna().sum()); nph2 = int(main.ph_z_controls.notna().sum())
hdr = [f"# Supplementary Table 1. Recurrent genes: the {len(main)} genes at q < 0.1 against random draws in the lists of the core 80-study base (72 lists contributing), counting each study once (in_analysis_set = yes for all {len(main)}). The schizophrenia z columns are blank where not computed ({nsz} of {len(main)} genes have a value).",
 "# Columns: gene; lists_naming, lists that report the gene; lists_measured, lists from studies that measured it; n_up, n_down, lists reporting it up and down; q_uniform, q_expression_matched, q_de_prior, Benjamini-Hochberg q values against random, expression-matched and DE-prior-weighted draws of genes (the recurrent set is defined by q_uniform; the DE-prior q values are for reference); n_clusters_observed, cohort clusters with a list naming the gene; q_cluster_bh, q value counting each cohort cluster once; in_analysis_set, yes for the recurrent genes used in all analyses; ph_z_controls, association with pH in controls on the pH axis (combined z over the control collections; above 0 means higher expression at higher pH); sz_z_raw and sz_z_ph_adjusted, schizophrenia effect (combined z) before and after adjustment for each study's measured pH (genes in 5 or more lists only).",
 "# code: <repository link to be added on release>"]
hdr_c = ["# Supplementary Table 1 (continued). 38 further genes often discussed in the literature, recurrent or not, with the same columns as the recurrent-gene table.", hdr[1], hdr[2]]
cands = ["C4A", "C4B", "DRD2", "GRIN2A", "GRIN2B", "SETD1A", "CACNA1C", "TCF4", "NRGN", "GAD1", "GAD2", "SST", "PVALB", "VIP", "BDNF", "VGF", "NPTX2", "SLC6A1", "RELN", "DISC1", "COMT", "AKT1", "HTR2A", "GRM3", "MIR137", "ZNF804A", "CHRNA7", "SHANK3", "NRXN1", "CUL1", "GABRA1", "SYT1", "SNAP25", "CHI3L1", "SERPINA3", "MT1X", "GADD45B", "HSPB1"]
cand = pd.DataFrame([row(x) for x in cands])
for name, df in (("recurrent_genes", main), ("named_candidates", cand)):
    with open(OUT / f"supp_table1_{name}.tsv", "w") as f:
        f.write("\n".join(hdr if name == "recurrent_genes" else hdr_c) + "\n"); df.round(5).to_csv(f, sep="\t", index=False)
open(OUT / "supp_table1_provenance.txt", "w").write("\n".join(prov) + "\n")
json.dump({"script": "L02/g1/export_supp_table1_383.py", "run": "2026-10-15",
 "inputs": {"analyses/recurrence_genes.tsv": sha(G), "analyses/confounds/recurrent_sets.tsv": sha(SET),
  "analyses/confounds/results/combined_controls_pH.tsv": sha(P), "inputs_ph_axis_v2/axis_v2_controls.tsv.gz": sha(AXV2), "figure_data/fig07/fig07_a2_a4_gene_scatter.tsv": sha(S)},
 "rows": {"recurrent_genes": len(main), "analysis_set_genes": int((main.in_analysis_set == "yes").sum()), "named_candidates": len(cand)}},
 open(OUT / "supp_table1_meta.json", "w"), indent=1)
print(len(main), (main.in_analysis_set == "yes").sum(), len(cand), (cand.in_analysis_set == "yes").sum())
