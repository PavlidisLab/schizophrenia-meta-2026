#!/usr/bin/env python3
"""Provenance (2026-10-19; no randomness; base: core-80, axis v2 where stated).
(1) Mann-Whitney P for Note 9: single-nucleus cohort-pair medians (median over cell classes within each Endresz cohort pair) vs bulk dataset-pair rho, both restricted to 11-30 donors per group.
    Core-80 input: ext_fig14_agreement_data/fig09/fig09_a_agreement_pairs.tsv. Validation: the same code on the pre-core-80 pairs file (figure_data/fig09, 153 bulk pairs in the bin) reproduces the stored 0.76.
(2) Frohlich (GSE254569) change in Pearson r after adding donor pH, 7 classes, model age_sex, gene set study_brain (label study_brain_209 in the file; n_genes up to 296 = the 383 set), min and max over classes, and classes whose bootstrap interval of the change excludes zero.
    Core-80 inputs: rerun_on_80_studies/sn_ph_and_hidden_factor_reruns/item10_sn_donor_ph/axis_v1 and axis_v2 (bootstrap seed 20261012, 1000 donor resamples, run by L09). Validation: L09 README: axis_v1 run on the earlier inputs reproduced 6 of 6 stored tables; the axis_v1 range here (0.11 to 0.33) equals the text."""
import pandas as pd,os,hashlib
from scipy.stats import mannwhitneyu
H=os.path.expanduser("")
def sh(p): return hashlib.sha256(open(p,'rb').read()).hexdigest()[:12]
print("part\tinput\tinput_sha256_12\tn_sn_pairs_11_30\tn_bulk_pairs_11_30\tmedian_sn\tmedian_bulk\tP_two_sided")
for lab,p in [("core80 (new)",H+"ext_fig14_agreement_data/fig09/fig09_a_agreement_pairs.tsv"),("pre-core80 (validation: stored 0.76)",H+"figure_data/fig09/fig09_a_agreement_pairs.tsv")]:
    d=pd.read_csv(p,sep="\t",comment="#"); dt=d.columns[0]
    b=d[d[dt].str.startswith("bulk")]; s=d[~d[dt].str.startswith("bulk")]
    med=s.groupby(["unit_a","unit_b"]).agg(rho=("rho_t","median"),n=("n_per_group","first")).reset_index()
    bb=b[b.n_per_group.between(11,30)].rho_t; m1=med[med.n.between(11,30)].rho
    print(f"MW\t{lab}\t{sh(p)}\t{len(m1)}\t{len(bb)}\t{m1.median():.3f}\t{bb.median():.3f}\t{mannwhitneyu(m1,bb).pvalue:.3f}")
print()
print("part\taxis\tinput_sha256_12\tgene_set\tmin_n_genes_over_classes\tmin_r_change\tmax_r_change\tclasses_interval_excludes_zero\tn_classes")
for ax in ("axis_v1","axis_v2"):
    p=H+f"rerun_on_80_studies/sn_ph_and_hidden_factor_reruns/item10_sn_donor_ph/{ax}/A_sz_effect_vs_ZpH_with_without_pH.tsv"
    d=pd.read_csv(p,sep="\t",comment="#"); g=d[(d.cohort=="GSE254569")&(d.model=="age_sex")]
    for gs,x in g.groupby("gene_set"):
        print(f"Frohlich\t{ax}\t{sh(p)}\t{gs}\t{x.n_genes.min()}\t{x.r_change.min():.3f}\t{x.r_change.max():.3f}\t{(x.r_change_lo>0).sum()}\t{len(x)}")
