---
modified: 2026-10-18
---
"# Figure captions

One entry per figure. Kept in sync with the figure on every edit.

## Figures

One entry per figure, in the numbering of `Figures_ph.md` and `supplement.md`; output files are named by content, not
number. Each figure is one matplotlib figure at print size (180 mm, 7 pt; `figures/scripts/figkit.py`). The full legends
are in `Figures_ph.md` and `supplement.md`.

### Fig. 1 — The corpus (`fig_corpus`)

What: project phases with publications screened, corpus studies and primary brain studies (a; counts read from figures/data/funnel/funnel_table.tsv); lineage of all corpus studies by year and source with reuse links (b); (cohort clusters, Extended Fig. 2f: 22 clusters, 17 of one study each; a study that pools several banks counts in the cluster of each bank). How to read: b, each point is a study, circled if in the 80; Extended Fig. 2f, grey bars are studies no meta-analysis or reanalysis has used. Shows: 28 of the 80 were never reused; reuse flows into meta-analyses from a few banks.
Source: `figures/scripts/build_fig_corpus.py`; data `figures/data/fig01/` (export_fig01_fig11.py, 2026-10-08).

### Extended Fig. 8 — Single-nucleus data (`fig_single_nucleus`)

What: one row of three panels. Agreement between cohorts with no subjects in common against donors per group, bulk dataset pairs vs single-nucleus cohort pairs with bin medians (a); PsychAD HBCC vs Mount Sinai per cell class at 20 donors per group, SZ and age, with within-bank split halves (b); ashr s.d. of true SZ effects by data type, studies and cell types (c). How to read: a, open circles at or below the blue bars mean single-nucleus pairs agree no better than bulk at the same N; b, the median row compares the four series; c, higher points mean larger true effects. Shows: single-nucleus medians 0.07 vs bulk 0.11 (0.16 vs 0.14 at 11 to 30 donors); across PsychAD banks SZ 0.04 vs age 0.17; true-effect s.d. 0.021 in single-nucleus vs 0.035 to 0.048 in bulk. Source: `figures/scripts/build_fig_single_nucleus.py`; data `figures/data/fig09/` (export_fig08_fig09_fig10.py, 2026-10-08).

### Extended Fig. 2 — The corpus by year, technology and sample size (`supp_corpus`)

What: corpus studies per year by tissue (a); technology by era (b); cases plus controls per study by technology (c); reported versus distinct subjects in the 80 (d); covariate handling in the 80 by era (e). How to read: c, one point is one study; d, bars are the sum of reported N and the lower and upper bounds on distinct people. Shows: the corpus grew after 2016 and moved to RNA-seq; 6,073 reported subjects are 1,858 to 4,401 people. Drawn by `figures/scripts/build_fig_corpus.py`.

### Supplementary Note 11 Figure 1 — Sample size and heterogeneity (`supp_sample_size`)

What: agreement vs subjects per group, split halves vs different cohorts, diagnosis and age (a); percentile of the 149 published different-cohort pairs among split halves of their size (b); Gemma reanalysis agreement by subject relation (c); agreement vs case-mix difference and comparable vs random subsets (d); donor identification vs diagnosis across dataset pairs (e); reliability of the GTEx slow-death axis (donors split within one dataset) and the control-pH axis (donors split within one collection, and between collections) (f). How to read: a, solid lines rise and dashed do not; b, a histogram piled left of 0.5 means published pairs sit below split halves; c, top rows share donors; d right, a blue dot inside the grey bar means no better than random; e, blue above grey means the person, not the diagnosis, is what carries over; f, further right is more repeatable, hollow points are an overstated value and a single collection. Shows: within-cohort agreement grows with n while cross-cohort agreement does not; published pairs sit at the 31st percentile of split halves; agreement tracks shared donors; case mix does not explain disagreement; donors are identifiable across datasets (0.76 to 0.99) while diagnosis is near chance; with all 285 GTEx donors the slow-death axis has reliability 0.94, but subsampled to the pH collections' size (50 and 45 donors) it is 0.64 and 0.59, close to the pH axis within one collection (0.73 and 0.55); between collections the pH axis is 0.42 to 0.55. Source: `figures/scripts/build_supp_sample_size.py`; data `figures/data/fig03/` (export_fig03_fig05.py, 2026-10-08); panel f `figures/data/figures_data_quality/panel_c.tsv` `split_half_ph_axis.tsv` and `split_half_gtex_at_ph_size.tsv`.

### Supplementary Note 11 Figure 2 — Positive controls (`supp_positive_controls`)

What: ρ(t) and top-200 overlap for diagnosis and age over the 192 pairs of bulk datasets with no subjects in common, full data and matched group sizes (a); per-dataset age agreement within the SZ-study data, against 9 external frontal-cortex datasets, and among the external datasets (b). How to read: a, one point per dataset pair, bars are medians, open points at the floor are pairs with zero top-200 overlap; b, one point per dataset, bars and numbers are medians over pairs. Shows: age reproduces better than diagnosis in the same data (0.14 vs 0.11), better still against external brain data (0.21), and sex not at all (0.00). Source: `figures/scripts/build_supp_positive_controls.py`; data `figures/data/fig05/`.

### Extended Fig. 3 — The recurrent genes (`supp_recurrent_genes`)

What: the recurrent genes by direction (a), majority-direction share (b), DE-prior AUC (c). How to read: b, blue far right of grey means hits agree on direction; c, bars near 1 mean the DE prior picks the set out almost perfectly. Shows: 85% vs 63% for the recurrent genes; AUC 0.88 and 0.92. Source: `figures/scripts/build_fig_recurrent_genes.py`; data `figures/data/fig06/`.

### Supplementary Note 10 Figure 1 — Antipsychotic treatment (`supp_medication`)

What: direction score of the 383 recurrent genes in 15 chronic drug-treated rodent experiments against matched random genes (a); leave-one-experiment-out pooled P (b); genes changed in the same or opposite direction as in SZ in rodent, human neural cell, and dose-within-patients data (c); each gene's mean drug z by SZ direction (d). How to read: a, right of 0 = drug changes the genes in the same direction as in SZ, open circle = matched null; c, blue share = same direction as SZ. Shows: a slight pooled trend in the SZ direction in chronic rodents (0.10 vs 0.04, P 0.046) that depends on a few experiments, absent or reversed elsewhere. Source: `figures/scripts/build_supp_medication.py`; data `figures/data/supp_fig_medication/`.

### Extended Fig. 7 — Other disorders (`supp_other_disorders`)

What: the 383 among DE genes of bipolar, autism, Alzheimer and depression postmortem brain tables against matched genes (a); share in the SZ direction among DE and all measured genes (b); disorder effect against pH and agonal axes, z vs matched genes (c); Hagihara 2018 pH case minus control by dataset, SZ and BD (d). Rows grouped: series sharing no donors with SZ data (main), tables that share donors (comparison), SZ and depression reference. Shows: autism DE genes include more of the 209 than matched genes; the 7 bipolar lists without shared donors have no DE genes but carry the 209 in the SZ direction (169 of 208); both correlate with the pH and agonal axes; depression only weakly; pH is as low in BD as in SZ (difference -0.006). Small dots are individual studies and larger markers the pooled table (per-study files panel_a/b/c_per_study.tsv); panel c is on pH axis v2 for all rows. Source: `figures/scripts/build_supp_other_disorders.py`; data `figures/data/supp_fig_other_disorders/`.

### Supplementary Note 5 Figure 1 — Pathways, coexpression modules and SNAP across studies (`supp_pathways_coexpression`)

What: agreement of two lists in top pathways and in genes over three nulls (a); control-pH association of shared term groups and of the recurrent genes (b); best match of Ling's SNAP factor in other datasets (c); schizophrenia effect on the SNAP score per donor set, and SNAP cell types minus other cell types (d); literature modules on the pH and agonal axes (e); the 196 modules by their schizophrenia difference in other cohorts, without and with the pH and agonal axis scores in the model (f). How to read: a, bars are observed over null, whiskers 95% intervals; c and d, one point is one dataset or donor set; e, one point is one paper; f, the stacked bar counts modules in the originally reported direction, opposite, or unclear. Shows: pathways agree more than genes and the excess sits in oxidative phosphorylation and myelination terms; SNAP re-derives in four of five datasets; the modules' schizophrenia differences are mostly unclear in other cohorts, and with the pH and agonal scores added 14 of 196 are in the reported direction and 12 opposite. Source: `supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_pathways_coexpression.py`; data `supplementary_notes/pathways_coexpression_modules_and_snap/data/`.

### Fig. 4 — The recurrent genes, the circumstances of death, and evidence for an acute effect (`fig_dying`)

What: GTEx slow-death vs pH by region (a); agreement after removing pH/agonal axes (b); per-study
agonal score (c); cell-class marker pH association (e).
How to read: a, negative ρ means genes lower at low pH rise with slow death. Shows: the pH signature
matches slow-death changes; cell type (neurons vs glia) predicts pH sensitivity. Source: `figures/scripts/build_fig_dying.py`; data `figures/data/fig08/` and
`figures/data/fig_interneuron/fig5g_marker_sets_vif.tsv`.

Panel d: other disorders' postmortem brain, correlation of each study's effects on the 383 with the pH and agonal axes minus matched random genes; one dot = one study, bar = median. Same data as Extended Fig. 7c (`figures/data/supp_fig_other_disorders_corrected/panel_c_per_study.tsv`); script `figures/scripts/build_fig_dying.py::panel_survey`.

Panel f: pH association in NSW single-nucleus controls at each pooling level against the bulk pH axis with the composition-only simulation; proportion of nuclei in each cell class against pH (crumblr); data `supplementary_notes/composition_and_axes/data/levels_summary_main_rows.tsv`, `sim_fraction_of_ceiling.tsv` and `composition_ph_summary_table.tsv`; script `figures/scripts/build_fig_dying.py::panel_levels, panel_fractions`.

### Fig. 2 — Agreement between studies (`fig_agreement`)

What: recurrence curve against three nulls (a); the AUROC of the 383 recurrent genes in the lists of 15 published brain meta-analyses and in all 15 ranked by recurrence, against the DE-prior and uniform expectations (b); gene-by-study calls
for the 80 most recurrent genes, studies in staircase (oncoprint memo-sort) order with strips for cohort cluster and N
(right, genes tested in grey and genes reported DE in black, log; dot above a gene symbol, top 10% of the DE prior; strips under the bars: association with pH, with a slow death and relative specificity to interneurons) (c).
How to read: a, the black line above a dashed null means more recurrence than that null predicts; b, a dot to the right of its grey bar means the 383 rank higher among that
list's tested genes than the DE prior predicts (top row: all 15 lists, genes ranked by the number of lists containing them). Shows: recurrence is above the DE-prior expectation from 3 lists (6,796
against 6,436) and, at 10 or more lists, is 107 genes against 44.3 expected; the pooled AUROC of the 383 is 0.87 against 0.73 expected from the DE prior, and 13 of 15 lists are above their DE-prior range alone (0.50 to 0.74). Data `figures/data/meta15/`.
Source: `figures/scripts/build_fig_agreement.py`; data `figures/data/fig02/` (export_fig02_fig07.py, 2026-10-08); panel a
from `analyses/code001_nk_correction_2026-10-11/gene_coverage/recurrence_curve.tsv`.

### Supplementary Note 4 Figure 1 — Meta-analyses agree in proportion to shared inputs (`supp_study_pairs`)

What: overlap of each pair of 15 brain meta-analyses over the DE-prior expectation, against the share of input datasets
the two have in common (Jaccard index). How to read: one point is one pair; blue, same technology; amber, different.
Shows: pairs sharing more datasets overlap more (Spearman 0.48 over 104 pairs). Source:
`figures/scripts/build_supp_study_pairs.py`; data `figures/data/fig02/fig02_e_meta_analysis_pairs.tsv`.

### Supplementary Note 5 Figure 3 — SNAP's astrocyte and neuron components against pH and agonal state (`supp_snap_ph`)

What: correlation of SNAP-a, SNAP-n and SNAP total with the pH-axis score, measured pH in controls and the estimated agonal score (a); schizophrenia effect on each score with age and sex in the model and with the agonal score or the pH-axis score added (b). How to read: one dot is one estimate with a 95% interval; a negative correlation with a pH measure means higher at low pH. Shows: SNAP-a is higher at low pH and SNAP-n lower, SNAP-n is lower in more agonal donors, the total shows little relation to pH, and the schizophrenia effect sits on SNAP-n and falls with the agonal score. Source: `supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_snap_ph.py`; data `supplementary_notes/pathways_coexpression_modules_and_snap/data/snap_ph_parts/`.

### Fig. 3 — The recurrent genes track pH (`fig_ph`)

What: pH difference per dataset-region (a); gene-level control-pH association against the schizophrenia effect (b); pH difference by cohort (c). How to read: b, a negative slope means genes lower at low pH are lower in schizophrenia. Extended Fig. 4a and b (agreement as covariates are added; the pH score of each study's list by pH handling; `supp_ph_covariates`) are drawn by the same script.
Source: `figures/scripts/build_fig_ph.py`; data `figures/data/fig07/` (export_fig02_fig07.py, 2026-10-08).

### Extended Fig. 5 — Slow death in every GTEx tissue (`supp_gtex_tissues`)

What: per GTEx tissue (13 brain regions, 32 others), the Spearman rho of the slow (Hardy 4) vs fast natural (Hardy 2) death effect with the brain slow-death axis (a), the low-pH pattern (b, control pH association sign-reversed) and the schizophrenia direction of the 176 recurrent genes (c); the same for Hardy 3 vs 2 (length of dying) against Hardy 4 vs 3 (long illness) (d); share of agreement between the 192 bulk SZ dataset pairs removed by each brain axis (e).
How to read: a-d, one point is one tissue; blue brain, grey other tissues; open circles in a are the three regions the brain axis is built from (circular).
Shows: slow death resembles the brain slow-death effect in 31 of 32 other tissues and is associated with the 176 genes differing in the schizophrenia direction in all 45, but resembles low brain pH only in brain; the 176 genes correlate more with the length of dying than with long illness, and in age- and sex-adjusted models that axis removes the most agreement (48% vs 39% and 20%; with all covariates 26% vs 37% and 8.5%).

### Extended Fig. 1 — Literature flow, PRISMA 2020 (`supp_prisma`)

What: records identified per route, deduplicated, screened on title/abstract (two model passes), corpus studies included on title and abstract, full texts sought, retrieved and assessed, exclusions by data use and scope, reports not added (corpus frozen), the 536 taken into the corpus, 21 targeted-panel exclusions, the 515, and the funnel to the 62.
How to read: main column top to bottom; boxes to the right are removals; blue boxes lead into the corpus.
Shows: 14,788 records screened, 1,139 full texts assessed, 665 studies in the corpus, 80 in the main analysis; final counts on the frozen corpus, arithmetic asserted in the script.
Source: `figures/scripts/build_supp_prisma.py`; data `figures/data/prisma_flow/prisma_flow.tsv`.

### Supplementary Note 3 Figure 1 — Recurrence over 248 distinct ways of building the lists (`supp_forking_paths_v2`)

What: the choices crossed (a); genes in 10 or more lists over expected under four chance models, one point per distinct analysis, published lists (b); the DE-prior ratio for published lists by studies included (c); signed pH score of each analysis's recurrent genes (d).
How to read: b, c log/linear ratio, dashed line 1, bars medians; d, positive means the recurrent genes are associated with low pH.
Shows: for published lists the median ratio at 10 or more lists is 3.23 against the DE prior and 4.6 to 13.7 by studies included; the pH signature is positive in every analysis.
Source: `figures/scripts/build_supp_forking_paths_v2.py`; data `figures/data/forking_paths_v2/`.

### Extended Fig. 6 — The recurrent genes across Gemma experiments (`supp_gemma_scan`)

What: for one random contrast per public human Gemma experiment, whether the 383 recurrent genes rank above DE-prior-matched expectation (a, by tissue; c, schizophrenia and oxygen-loss experiments), and how closely each contrast correlates with the low-pH and slow-death patterns (b). How to read: a contrast counts when the 383 score above the 97.5th percentile of that contrast's own 200 DE-prior-matched random sets; bars are shares with Wilson 95% intervals; dotted line, the random sets' own rate (3%). Shows: the 383 stand out mainly in postmortem brain (45%), and its contrasts correlate with low pH and slow death.
Source: `figures/scripts/build_supp_gemma_scan.py`; data a, c `figures/data/gemma_209_scan_calibrated/`; b `figures/data/gemma_176_scan_figure/panel_b.tsv`.

### Supplementary Note 5 Figure 2 — SNAP in bulk tissue (`supp_snap_bulk`)

What: per bulk dataset, coexpression of SNAP's top 1,000 genes over expression-matched genes (a); schizophrenia effect on the bulk SNAP score with age and sex, and with the estimated agonal score added (b).
How to read: points with 95% intervals; diamonds pooled over 22 brain-bank datasets (cohort-cluster resampling).
Shows: SNAP genes covary above matched genes in every bulk dataset; the bulk schizophrenia effect shrinks toward zero when the estimated agonal score is in the model.
Source: `supplementary_notes/pathways_coexpression_modules_and_snap/build_supp_snap_bulk.py`; data `figures/data/snap_what_is_it/`.


### Supplementary Note 7 Figure 1 — Held and autopsy tissue against the pH and agonal axes (`fig_yaqubi_axes`)

What: Spearman correlation of the gene-level effect of each handling condition (surgical tissue held at 4 or 20 C for 6, 24 or 36 h; autopsy tissue) and of the Dachet ex vivo time course with the pH and agonal axes, over all genes (a) and the 383 recurrent genes (b), and the share of the 383 whose schizophrenia sign equals the effect sign (c).
How to read: dots, observed value with bootstrap (a, b) or Wilson (c) 95% interval over genes; translucent bars, 95% range in matched gene sets; open squares, Dachet; dashed lines in c, the same share for the pH and agonal axes.
Shows: see Supplementary Note 7.
Source: `supplementary_notes/death_signatures/build_supp_death_signatures.py`; data `supplementary_notes/death_signatures/data/yaqubi_*.tsv`.

### Supplementary Note 7 Figure 2 — Living against postmortem cortex and the axes (`fig_livpm_vs_axes_383`)

What: A, Spearman correlation across genes of each living-minus-postmortem signature with the pH and agonal axes; B, share of the 383 recurrent genes with a direction for which the signature points the same way as schizophrenia (grey box, 2.5 to 97.5% of 10,000 matched gene sets); C, mean signature t of the 100 most cell-class-specific genes of each class minus matched genes.
How to read: dots with 95% bootstrap intervals over genes; a dot above its grey box means more agreement with schizophrenia than matched genes.
Shows: see Supplementary Note 7.
Source: `supplementary_notes/death_signatures/build_supp_death_signatures.py`; data `supplementary_notes/death_signatures/data/liharska_*.tsv`.

### Supplementary Note 8 Figure 1 — Cell composition and the pH and agonal axes (`fig_composition_key`)

What: Spearman correlation across genes between the pH association in 32 NSW control donors and the bulk pH axis (a) or agonal axis (b) at four pseudobulk levels, with the composition-only simulation as a band; the cell-class composition effects of pH and schizophrenia (c).
How to read: a, b, points are correlations with 95% donor-bootstrap intervals; the orange band is the simulated composition-only correlation; c, bars are 95% intervals of the difference in centred log-ratio composition.
Shows: the pH association's correlation with the bulk pH axis falls from 0.23 to about 0.05 when composition is held constant and the simulation gives 0.15 to 0.18; the agonal axis keeps more (−0.27 to −0.16); pH and schizophrenia shift composition in the same direction.
Source: `supplementary_notes/composition_and_axes/build_supp_composition.py`; data `supplementary_notes/composition_and_axes/data/`.
