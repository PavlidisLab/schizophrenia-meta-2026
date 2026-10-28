# figures/data: index

Generated 2026-10-14 from `READS.tsv` (which figure script reads which file) and `MANIFEST.tsv` (sha256, size and producing script of each file). The folder names are those of the delivery folder the tables arrive in (`sync_figure_data.py` keeps each file's relative path), so `figNN` is the number in the export script's own numbering, **not** the number of a figure in the manuscript. Use the table below to find a file's figure.

| Folder | Files | Read by (figure script) |
|---|--:|---|
| `fig01/` | 7 | `fig_agreement`, `fig_corpus` |
| `fig03/` | 7 | `supp_sample_size` |
| `fig05/` | 3 | `supp_positive_controls` |
| `fig08/` | 25 | `fig_dying` |
| `fig09/` | 10 | `fig_single_nucleus` |
| `fig_interneuron/` | 1 | `fig_dying` |
| `figures_data_quality/` | 3 | `supp_sample_size` |
| `forking_paths_v2/` | 4 | `supp_forking_paths_v2` |
| `funnel/` | 1 | `fig_corpus` |
| `gemma_176_scan_figure/` | 1 | `supp_gemma_scan` |
| `gemma_209_scan_calibrated/` | 2 | `supp_gemma_scan` |
| `gtex_death_speed_all_tissues/` | 3 | `supp_gtex_tissues` |
| `ph_by_cohort/` | 1 | none (not read by any figure) |
| `prisma_flow/` | 1 | `fig_corpus`, `supp_prisma` |
| `supp_fig_medication/` | 6 | `supp_medication` |
| `supp_fig_other_disorders/` | 1 | `supp_other_disorders` |
| `supp_fig_other_disorders_corrected/` | 3 | `supp_other_disorders` |

## Figure scripts

| Script (`build_<name>.py`) | First line of its docstring |
|---|---|
| `fig_agreement` | Figure 2 -- agreement between published brain studies is near chance. |
| `fig_corpus` | Figure 1 -- the corpus. |
| `fig_dying` | Figure 5 -- the pH signature resembles the dying process; which genes are pH-sensitive. |
| `fig_ph` | Figure 7 -- the recurrent genes track brain pH in controls. |
| `fig_recurrent_genes` | Supplementary Figure 6 (the former Fig. 4) -- the recurrent genes: direction, sources, co-occurrence, DE prior. |
| `fig_single_nucleus` | Fig. 3 -- single-nucleus data do not agree better than bulk and carry the same pH signature. |
| `supp_forking_paths_v2` | Supplementary Fig. 2 (v2): the recurrence result over the distinct ways of building the lists. |
| `supp_gemma_scan` | Extended Fig. 6 -- where the 209 recurrent genes stand out across Gemma experiments. |
| `supp_gtex_tissues` | Supplementary Fig. 9 -- a slow death in every GTEx tissue: like brain's slow-death effect, like the schizophrenia |
| `supp_medication` | Supplementary Fig. 7 -- antipsychotic treatment and the 176 recurrent genes. |
| `supp_other_disorders` | Extended Fig. 7 -- the 209 recurrent schizophrenia genes in other disorders' postmortem brain. |
| `supp_positive_controls` | Supplementary Fig. 3: positive control (age) beside diagnosis, and age against external brain data. |
| `supp_prisma` | Supplementary Fig. (PRISMA) -- literature flow, PRISMA 2020 style, from the search to the corpus and on to the primary brain studies. |
| `supp_proteomics` | Supplementary Figure 4 -- the brain pH signature in postmortem brain proteomics. |
| `supp_sample_size` | Supplementary Fig. 1: agreement against sample size, and differences between patients. |
| `supp_study_pairs` | Supplementary Note 4 Figure 1 -- overlap of pairs of brain meta-analyses against the share of input datasets they have in common. |
| `supp_snap_bulk` | Supplementary Note 5 Figure 2 -- SNAP in bulk tissue. |
| `supp_pathways_coexpression` | Supplementary Note 6 figure -- pathways and coexpression modules: no stronger agreement than genes; what agrees tracks pH. |

## Files by figure script

- `fig_agreement`: fig01/fig01_b_corpus_studies.tsv, fig02/fig02_a_recurrence_curve.tsv, fig02/fig02_b_gene_by_study_long.tsv, fig02/fig02_b_gene_de_prior.tsv, fig02/fig02_b_study_annotation.tsv, fig02/fig02_c_gemma_reanalysis_pairs.tsv, fig02/fig02_c_gemma_reanalysis_pooled.tsv, fig02/fig02_c_pairwise_fold_by_era.tsv, fig02/fig02_c_pooled_by_era.tsv, fig02/fig02_d_hoffman82_by_group.tsv, fig02/fig02_d_hoffman82_by_study.tsv, fig02/fig02_e_meta_analysis_pairs.tsv, fig02/fig2a_unit_band_centred.tsv, fig02/fig2a_unit_bootstrap.tsv, fig02/list_nk_vs_M.tsv, fig07/fig07_a2_a4_gene_scatter.tsv
- `fig_corpus`: fig01/fig01_b_corpus_studies.tsv, fig01/fig01_b_distinct_subjects.tsv, fig01/fig01_b_reporting_completeness.tsv, fig01/fig01_c_lineage_edges.tsv, fig01/fig01_c_lineage_nodes.tsv, fig01/fig01_d_neglected_datasets.tsv, fig01/fig01_f_covariate_calls_64.tsv, fig02/fig02_b_study_annotation.tsv, funnel/funnel_table.tsv, prisma_flow/prisma_flow.tsv
- `fig_dying`: fig08/fig08_a_alignment_by_time.tsv, fig08/fig08_a_dachet_gene_scatter.tsv, fig08/fig08_a_recurrent_direction.tsv, fig08/fig08_b_alignment_medians_by_target_set.tsv, fig08/fig08_b_axis_correlations.tsv, fig08/fig08_b_dataset_axis_alignment.tsv, fig08/fig08_b_pair_agreement_gtex_axis.tsv, fig08/fig08_b_pair_agreement_single_axis_same_genes.tsv, fig08/fig08_b_pair_medians_by_target_set.tsv, fig08/fig08_b_pairs_remote.tsv, fig08/fig08_c_class_direction.tsv, fig08/fig08_c_gene_ph_specificity.tsv, fig08/fig08_c_reference_set_tests.tsv, fig08/fig08_d_interneuron_score_by_dataset.tsv, fig08/fig08_d_interneuron_score_by_model.tsv, fig08/fig08_d_single_gene_papers_by_ph_handling.tsv, fig08/fig08_e_gtex_by_region.tsv, fig08/fig08_f_expression_summary.tsv, fig08/fig08_f_gse92538_subjects.tsv, fig08/fig08_f_ph_by_manner_of_death.tsv, fig08/fig08_f_ph_effect_by_dataset.tsv, fig08/fig08_f_ph_effect_pooled.tsv, fig08/fig08_f_ph_models.tsv, fig08/fig08_f_ph_natural_death_only.tsv, fig08/fig08_g_hitlist_agonal_score_by_study.tsv, fig_interneuron/fig5g_marker_sets_vif.tsv
- `fig_ph`: fig07/fig07_a1_ph_difference.tsv, fig07/fig07_a2_a4_gene_scatter.tsv, fig07/fig07_a3_cohort_sets_corrected_null.tsv, fig07/fig07_a3_cohort_sets_genes.tsv, fig07/fig07_a3_observed_vs_expected.tsv, fig07/fig07_a_covariate_series_pH_RIN_steps.tsv, fig07/fig07_a_covariate_series_pH_RIN_steps_summary.tsv, fig07/fig07_a_covariate_series_with_pH_step.tsv, fig07/fig07_b_covariate_series_genes.tsv, fig07/fig07_b_covariate_series_pairs.tsv, fig07/fig07_b_levels_by_dataset.tsv, fig07/fig07_c_ph_score_by_adjustment.tsv, fig07/fig07_c_ph_score_by_handling_v2.tsv, fig07/fig07_c_ph_score_by_handling_v2_summary.tsv
- `fig_single_nucleus`: fig09/fig09_a_agreement_by_n_bin.tsv, fig09/fig09_a_agreement_pairs.tsv, fig09/fig09_a_single_nucleus_by_class.tsv, fig09/fig09_b_medians_over_8_classes.tsv, fig09/fig09_b_psychad_20_per_group.tsv, fig09/fig09_c_ph_vs_sn_effect_by_subclass.tsv, fig09/fig09_d_effect_sd_by_study.tsv, fig09/fig09_d_effect_sd_single_nucleus_cell_types.tsv, fig09/fig09_d_endresz_cohorts.tsv, fig09/fig09_d_summary_by_type.tsv
- `supp_forking_paths_v2`: fig02/fig02_b_study_annotation.tsv, forking_paths_v2/forking_v2_a_design.tsv, forking_paths_v2/forking_v2_b_spec_curve.tsv, forking_paths_v2/forking_v2_c_ph_paths.tsv, forking_paths_v2/forking_v2_e_paths.tsv
- `supp_gemma_scan`: gemma_176_scan_figure/panel_b.tsv, gemma_209_scan_calibrated/panel_a.tsv, gemma_209_scan_calibrated/panel_c.tsv
- `supp_gtex_tissues`: fig06/fig06_b_direction_summary.tsv, gtex_death_speed_all_tissues/pair_agreement_hardy3_axis_removal.tsv, gtex_death_speed_all_tissues/rho_by_tissue_and_contrast.tsv, gtex_death_speed_all_tissues/rho_hardy3_contrasts_by_tissue.tsv
- `supp_medication`: supp_fig_medication/a_unit_tests.tsv, supp_fig_medication/b_leave_one_out.tsv, supp_fig_medication/c_gene_effects.tsv, supp_fig_medication/e_units.tsv, supp_fig_medication/f_fromer2016_rhesus.tsv, supp_fig_medication/h_dose_within_sz.tsv
- `supp_other_disorders`: fig06/fig06_b_direction_summary.tsv, supp_fig_other_disorders/panel_d_hagihara_ph_by_diagnosis.tsv, supp_fig_other_disorders_corrected/panel_a_enrichment_vs_null.tsv, supp_fig_other_disorders_corrected/panel_b_concordance_vs_null.tsv, supp_fig_other_disorders_corrected/panel_c_axis_effect.tsv
- `supp_positive_controls`: fig05/fig05_a_pairs_dx_age_sex.tsv, fig05/fig05_b_age_external_per_dataset.tsv, fig05/fig05_b_age_external_summary.tsv
- `supp_prisma`: prisma_flow/prisma_flow.tsv
- `supp_recurrent_genes`: fig02/fig02_b_study_annotation.tsv, fig06/fig06_a_gene_direction.tsv, fig06/fig06_a_gene_direction_211.tsv, fig06/fig06_b_direction_summary.tsv, fig06/fig06_b_direction_summary_both_versions.tsv, fig06/fig06_c_strata.tsv, fig06/fig06_d_cooccurrence_cohort_tier.tsv, fig06/fig06_d_cooccurrence_cohort_tier_200k.tsv, fig06/fig06_d_cooccurrence_study_tier.tsv, fig06/fig06_d_cooccurrence_study_tier_200k.tsv, fig06/fig06_e_de_prior_auc_by_set.tsv, fig06/fig06_e_de_prior_per_list_auroc.tsv, fig06/fig06_e_de_prior_top_fraction.tsv, fig06/recurrent_set_membership_209_211_115_150.tsv
- `supp_sample_size`: fig03/fig03_a_c_d_agreement_summary.tsv, fig03/fig03_a_c_d_agreement_units.tsv, fig03/fig03_b_published_vs_split_halves.tsv, fig03/fig03_e_gemma_rho_by_relation.tsv, fig03/fig03_f_casemix_pairs.tsv, fig03/fig03_f_comparable_vs_random.tsv, fig03/fig03_g_identification_by_relation.tsv, figures_data_quality/panel_c.tsv, figures_data_quality/split_half_gtex_at_ph_size.tsv, figures_data_quality/split_half_ph_axis.tsv
- `supp_study_pairs`: fig02/fig02_e_meta_analysis_pairs.tsv
- `supp_snap_bulk`: snap_what_is_it/snap_bulk_coherence.tsv, snap_what_is_it/snap_estimated_agonal.tsv
- `supp_pathways_coexpression`: suppnote5_coexpr_null/panel_b_recurrent_row_combined.tsv, suppnote5_coexpr_null/suppnote5_a_coexpr_null_current.tsv, suppnote6_pathways_coexpression/panel_k_snap_best_match_chu_vs_sz.tsv, suppnote6_pathways_coexpression/panel_l_snap_specificity.tsv, suppnote6_pathways_coexpression/panel_o_deans_convergent_genes_recurrence.tsv, suppnote6_pathways_coexpression/suppnote6_d_pH_test_term_groups.tsv, suppnote6_pathways_coexpression/suppnote6_e_coexpression_claim_classes.tsv, suppnote6_pathways_coexpression/suppnote6_f_retested_claims.tsv, suppnote6_pathways_coexpression/suppnote6_g_snap_other_snrnaseq.tsv, suppnote6_pathways_coexpression/suppnote6_h_modules_signed_axes_by_paper.tsv, suppnote6_pathways_coexpression/suppnote6_i_gandal_copa_axes.tsv

## Files read by no figure script

- `ph_by_cohort/ph_by_cohort.tsv`
