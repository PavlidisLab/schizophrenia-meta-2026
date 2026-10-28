#!/usr/bin/env python3
"""Writes SUMMARY_for_writeup.md from tables/*.tsv (no typed numbers). Usage: python3 07_summary.py <out dir>"""
import sys, re, datetime
sys.dont_write_bytecode = True
import pandas as pd
OUT = sys.argv[1]
rd = lambda f: pd.read_csv(f"{OUT}/tables/{f}", sep="\t", comment="#")
AL, AS, CO, OV, DA, DEF = rd("axis_correlation_all_genes.tsv"), rd("axis_correlation_recurrent_sets.tsv"), rd("direction_concordance.tsv"), rd("overlap_processing_sets.tsv"), rd("dachet_axis_correlation.tsv"), rd("contrast_definitions.tsv")
BO, BA, BG = rd("bag_vs_recurrent.tsv"), rd("bag_pattern_axis_correlation.tsv"), rd("bag_go_themes_overlap.tsv")
SD = pd.read_csv(f"{OUT}/tables/surgical_sample_design.tsv", sep="\t")
PH, AG, S3, S1 = "pH axis v2", "agonal axis (GTEx slow-minus-fast)", "383 recurrent (q_bh<0.1)", "107 in >=10 lists"
nd = SD.donor.nunique()
def c(r, v="spearman", lo="spearman_boot95_low", hi="spearman_boot95_high"):
    return f"{r[v]:+.2f} [{r[lo]:+.2f}, {r[hi]:+.2f}]"
def sset(eff, ax, gs):
    r = AS[(AS.gene_set == gs) & (AS.effect == eff) & (AS.axis == ax) & (AS.statistic == "spearman")].iloc[0]
    return c(r, "value", "boot95_low", "boot95_high")
def conc(eff, gs):
    r = CO[(CO.gene_set == gs) & (CO.effect == eff)].iloc[0]
    return f"{r.frac_concordant:.2f} [{r.wilson95_low:.2f}, {r.wilson95_high:.2f}] vs null {r.null_mean:.2f}"
def ovl(nm, gs):
    r = OV[(OV.gene_set == gs) & (OV.processing_set == nm)].iloc[0]
    return f"{r.overlap}, {r.fold_vs_matched_null}x [{r.fold_null_q025}, {r.fold_null_q975}]"
CON = [("S_6h4C_vs_0h", "6 h 4 C vs frozen", "6h4C", "S_6h4C_vs_0h"), ("S_24h4C_vs_0h", "24 h 4 C", "24h4C", None), ("S_36h4C_vs_0h", "36 h 4 C", "36h4C", None),
       ("S_6h20C_vs_0h", "6 h 20 C", "6h20C", None), ("S_24h20C_vs_0h", "24 h 20 C", "24h20C", None), ("S_36h20C_vs_0h", "36 h 20 C", "36h20C", None),
       ("S_20C_minus_4C_6h", "20 minus 4 C, 6 h", None, None), ("S_20C_minus_4C_24h", "20 minus 4 C, 24 h", None, None), ("S_20C_minus_4C_36h", "20 minus 4 C, 36 h", None, None),
       ("A_SPMI_vs_CTRL", "Autopsy short PMI vs CTRL", None, None), ("A_LPMI_vs_CTRL", "Autopsy long PMI vs CTRL", None, None)]
L = []; w = L.append
w("# Yaqubi 2026 processing conditions against the pH axis v2, the agonal axis and the recurrent genes: summary (NEW, 2026-10-18)\n")
w(f"Generated {datetime.date.today().isoformat()} by scripts/07_summary.py from tables/*.tsv; details, file list and sources in README.md. Effects are condition minus control (positive = higher after processing/death). Intervals: bootstrap over genes (rho), Wilson (shares), matched-null 95% range (folds). pH axis v2: positive = higher at higher pH. Agonal axis: GTEx slow-minus-fast-death t (positive = higher in slow death).\n")
w("## Design and definitions\n")
w(f"{nd} surgical donors (cortex from tumour-resection access, frozen within ~31 min = CTRL); {len(SD)} libraries (donor-paired, two duplicate 0 h libraries); six portions per donor at 4 or 20 C for 6, 24 or 36 h. Autopsy: SPMI (n 9, mean 6 h, Netherlands Brain Bank, superior temporal gyrus) and LPMI (n 18, mean 36 h, UK banks, mid temporal gyrus) against the same CTRL (n 17). Re-fit from the released raw counts: OLS on log2(median-of-ratios count + 1), donor fixed effects for the holds, group + age + sex + RIN for autopsy; genes with mean normalised count below 10 dropped. Matched null: expression decile x DE-prior quintile, 10,000 sets (5,000 for correlations). 383 = study_brain, 107 = study_brain_ge10.\n")
w("## Per contrast\n")
w("| contrast | n | pH: all genes | pH: 383 | pH: 107 | agonal: all genes | agonal: 383 | agonal: 107 | SZ sign = effect sign (383) | DEG overlap with 383, up / down |\n|---|---|---|---|---|---|---|---|---|---|")
for eff, lab, cond, _ in CON:
    d = DEF[DEF.contrast == eff].iloc[0]
    n = f"{d.n_samples} libraries" if eff.startswith("S_") else f"{d.n_samples} ({d.definition.split('groups n = ')[1]})"
    a = AL[AL.effect == "t:" + eff]
    pa, ga = a[a.axis == PH].iloc[0], a[a.axis == AG].iloc[0]
    dg = ""
    if eff.startswith("S_") and "minus" not in eff:
        dg = f"{ovl(eff + ' DEG up', S3)} / {ovl(eff + ' DEG down', S3)}" if ((OV.processing_set == eff + ' DEG up').any()) else "no DEGs"
    elif eff.startswith("A_"):
        k = eff[2:6] + " vs CTRL"
        dg = f"{ovl(k + ' DEG up', S3)} / {ovl(k + ' DEG down', S3)}"
    w(f"| {lab} | {n} | {c(pa)} | {sset('t:' + eff, PH, S3)} | {sset('t:' + eff, PH, S1)} | {c(ga)} | {sset('t:' + eff, AG, S3)} | {sset('t:' + eff, AG, S1)} | {conc('t:' + eff, S3)} | {dg} |")
for h in sorted(DA.hours.unique()):
    p = DA[(DA.hours == h) & (DA.scope == "all genes") & (DA.axis == PH)].iloc[0]; g = DA[(DA.hours == h) & (DA.scope == "all genes") & (DA.axis == AG)].iloc[0]
    def dsc(ax, gs):
        r = DA[(DA.hours == h) & (DA.scope == gs) & (DA.axis == ax)].iloc[0]; return c(r, "spearman", "boot95_low", "boot95_high")
    w(f"| Dachet {h} h vs 0 h (GSE153669) | 1 patient, 7 time points | {c(p, 'spearman', 'boot95_low', 'boot95_high')} | {dsc(PH, S3)} | {dsc(PH, S1)} | {c(g, 'spearman', 'boot95_low', 'boot95_high')} | {dsc(AG, S3)} | {dsc(AG, S1)} | not computed | not computed |")
w("\nDEG overlap = genes of the supplement's DEG list (|log2FC| above 1, adjusted p below 0.05) inside the 383: count, fold over the matched null, null 95% range as a fold. 6 h 4 C has no DEGs. Released-log2FC versions of the correlations (`relLFC:`) are in the tables.\n")
w("## Brain Artifact Genes (BAGs; Supplementary Data 4; common to SPMI and LPMI vs CTRL)\n")
nb = BO[BO.bag_set == "BAG all"].iloc[0]
LG = re.search(r"BAG released rows (\d+) (\d+) mapped unique (\d+) up (\d+) down (\d+)", open(f"{OUT}/logs/06_bag.log").read()).groups()
w(f"Released: {LG[0]} up and {LG[1]} down; after HGNC mapping of the released symbols and removal of duplicates the BAG set is {LG[2]} genes ({LG[3]} up, {LG[4]} down; the rest did not map or were repeated). Of these {int(BO[(BO.gene_set==S3)&(BO.bag_set=='BAG up')].iloc[0].bag_genes_in_universe)} up and {int(BO[(BO.gene_set==S3)&(BO.bag_set=='BAG down')].iloc[0].bag_genes_in_universe)} down are in the published-brain gene universe of the recurrence analysis.\n")
w("| BAG subset | set | overlap | fold over matched null [null range] | BAG direction = SZ sign | null share [range] |\n|---|---|---|---|---|---|")
for gs in (S3, S1):
    for b in ("BAG up", "BAG down", "BAG all"):
        r = BO[(BO.gene_set == gs) & (BO.bag_set == b)].iloc[0]
        w(f"| {b} | {gs} | {r.overlap} of {r.bag_genes_in_universe} | {r.fold}x [{r.fold_null_q025}, {r.fold_null_q975}] | {r.concordant} of {r.n_with_SZ_sign} = {r.frac_concordant:.2f} [{r.wilson95_low:.2f}, {r.wilson95_high:.2f}] | {r.null_frac_mean:.2f} [{r.null_frac_q025:.2f}, {r.null_frac_q975:.2f}] |")
w("\nBAG log2FC (mean of the SPMI and LPMI values) against the axes, Spearman [bootstrap 95% over genes]:\n")
w("| genes | n (pH / agonal) | pH axis v2 | agonal axis |\n|---|---|---|---|")
for sc in ("all BAG genes", "BAG genes in the 383", "BAG genes in the 107"):
    p = BA[(BA.scope == sc) & (BA.axis == PH) & (BA.bag_lfc == "lfc_mean")].iloc[0]; g = BA[(BA.scope == sc) & (BA.axis == AG) & (BA.bag_lfc == "lfc_mean")].iloc[0]
    w(f"| {sc} | {p.n_genes} / {g.n_genes} | {c(p)} | {c(g)} |")
def gx(subset, other):
    r = BG[(BG.bag_subset.str.startswith(subset)) & (BG.compared_with.str.startswith(other))].iloc[0]
    return f"{r.overlap} of {r.n_subset_in_axis_universe} (expected {r.expected_if_independent})"
ox = [s for s in BG.bag_subset.unique() if "oxidative" in s and s.startswith("BAG up")][0]
lp = [s for s in BG.bag_subset.unique() if "lipopolysaccharide" in s and s.startswith("BAG up")][0]
w(f"\nGO themes the paper names (annotation from analyses/fairness_lists/go_bp_sets.tsv.gz; BAG-up genes annotated to the oxidative phosphorylation / ATP synthesis terms, and to response to lipopolysaccharide). Overlap, with the expectation if independent over the {int(BG.iloc[0].background_genes)} genes of the pH axis (hypergeometric, not DE-prior matched):\n")
w("| BAG-up subset | pH-axis top 200 (higher at higher pH) | pH-axis bottom 200 (higher at lower pH) | oxphos term group, GO members | oxphos term group, list hits | 383 |\n|---|---|---|---|---|---|")
for s in ("BAG up (all)", ox, lp):
    w(f"| {s} | {gx(s, 'pH-axis top')} | {gx(s, 'pH-axis bottom')} | {gx(s, 'oxphos term group, GO')} | {gx(s, 'oxphos term group, genes')} | {gx(s, '383')} |")
w("\nThe oxphos term group is the one rebuilt in rerun_core80_2026-10-16/M09 (14 GO BP terms; `term_group_genes.tsv`: 216 GO members, of which the 'terms' count above lie in the pH-axis universe, and 130 list hits). A 166-gene list of that name was not found on disk; these are the current files. Mitochondrially encoded genes the paper names (COX2, COX3, ATP6, ATP8) are mostly absent from the axis universe, so the oxphos BAG subset is small.\n")
w("## Limits\n")
w("- Own t statistics are an OLS refit from released counts (not DESeq2, no shrinkage); agreement with the released log2FC (Spearman, own vs released): " + "; ".join(f"{r.contrast} {r.validation_spearman:.2f}" for r in DEF.dropna(subset=['validation_spearman']).itertuples()) + ".")
w(f"- {nd} surgical donors, one lab, one batch; tumour-resection access tissue; holds in covered PBS with temperatures set by the investigators; unbalanced (some donors lack conditions).")
w("- Surgical tissue is not agonal-state tissue: no terminal hypoxia, acidosis or illness precedes the hold, so the holds can speak to time and temperature only.")
w("- Autopsy contrasts differ from CTRL in cohort (brain bank, region, age, disease status) as well as processing; BAGs share this confound (the paper says so).")
w("- Dachet et al. is one patient with seven time points; no donor variance.")
w("- Genes are not independent; intervals over genes are descriptive; BAG-in-383 and BAG-in-107 subsets are small (tens of genes); overlaps and GO subsets are not DE-prior matched except where stated. The paper's 237-gene artifact-responsive set and its 15.8% and 11.3% overlaps were not reproduced exactly (263 genes recomputed).")
open(f"{OUT}/SUMMARY_for_writeup.md", "w").write("\n".join(L) + "\n")
