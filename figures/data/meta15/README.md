# META15: recurrent genes in the 15 brain meta-analysis lists (core-80 base), 2026-10-16

Job P10 for RA_TO_SZW_2026_10_16_META15_RECURRENT_GENES.md plus the coordinator addendum (13 vs 15; input coverage). Base is 80 studies / 72 contributing lists. Sets: 383 (study tier, q<0.1) and 107 (in 10 or more lists).

Parameters: NPERM 5000, seed 20261016 (joint draws: all 15 lists drawn in every replicate). Nulls: (a) DE prior (h3.py `expected_deprior` / `recurrence.prior_probs`, prior file `analyses/de_prior_lists/our_prior/1_global_our_prior.tsv`, weight n_de_total+1, median for missing genes, inclusion probability capped at 1, applied within each unit's universe; here a list is drawn as independent Bernoulli(p_g), mean size = the list's size); (b) uniform: the same number of genes without replacement from the unit's universe. Universe = the unit's own tested universe (universes.json gene set) for 12 units; Li 2022, Nakatsuka 2025 and Endresz 2026 have none (partner's universe used in h3), so they use the union of tested sets of the 72 non-empty base-80 lists (34,365 genes). Each list is restricted to its universe (Hess 2016 92 to 90, Priyanka 2025 605 to 444; others unchanged); Chen 2023's own universe is the 1,275-gene targeted panel.

## Inputs (sha256, full list in `sources_sha256.tsv` and in each table header)
- recurrent_sets.tsv (B80, 383 and 107): B80/out/recurrent_sets.tsv  sha256 3c42f67cdd7acb751592af24bae4c2c1a917bc5a880b8fe74d2d14cd8ef85403
- base80_studies.tsv: base80_recurrence_2026-10-16/base80_studies.tsv  sha256 09dc9e52c590a731f906881604d85132b058a9bddd445899ea766678f8a6187a
- set_gene_differences.tsv (read for the comparison; sets taken from recurrent_sets.tsv, which it was built from): sha256 bee99ba028aec06b779826144b6f6293ebdbf9e321fa6da8aa709d717024eff4
- h3_units.tsv, h3_inventory.tsv, h3_datasets.tsv (repo analyses/h3), 15 meta-analysis list files data/lists/<id>.tsv, DE prior file, universes.json and primary_lists.json (B80 ws), lineage.json, meta_link_evidence.tsv (2026-10-12, copied unchanged as `in_meta_link_evidence.tsv`) and meta_link_audit.md: sha256 in `sources_sha256.tsv` / table headers.

## The 15 versus 13 question
All 15 units are cross-cohort, pool 2 or more cohorts (2 to 7), publish a list and are in_pairs = yes, so no criterion on those drops any. No single rule is declared; two rules from the unit definitions give exactly 13 (`meta15_units_rule.tsv`):
- Variant A, at least 4 input cohorts: drops Gandal 2018 (RNA-seq) (2 cohorts) and Nakatsuka 2025 (3 cohorts).
- Variant B, bulk tissue only (drops the single-nucleus units): Nakatsuka 2025 and Endresz 2026 bioRxiv (also the only preprint, and the only list with role "reanalysis").
Other single criteria drop one unit (list of fewer than 10 genes: Bharadhwaj 2023; preprint: Endresz 2026). Every number is given for all 15, variant A and variant B;

## Per meta-analysis hits (observed vs expected, 95% interval of the simulation)
### 383 set
| meta-analysis | list genes in universe | universe | obs | uniform exp (95%) | DE-prior exp (95%) |
|---|---|---|---|---|---|
| Pérez-Santiago 2012 | 131 | own 13634 | 49 | 3.190 (0.0-7.0) | 5.519 (2.0-11.0) |
| Mistry 2013 | 98 | own 13634 | 17 | 2.386 (0.0-6.0) | 4.129 (1.0-9.0) |
| Bergon 2015 | 265 | own 21383 | 51 | 4.722 (1.0-9.0) | 8.918 (4.0-15.0) |
| Hess 2016 | 90 | own 21383 | 17 | 1.604 (0.0-4.0) | 3.029 (0.0-7.0) |
| Qin 2016 | 607 | own 13634 | 100 | 14.781 (8.0-22.0) | 25.574 (16.0-35.0) |
| Gandal 2018 (RNA-seq) | 4538 | own 21576 | 223 | 80.345 (65.0-96.0) | 144.899 (128.0-163.0) |
| Gandal 2018 (arrays) | 2028 | own 12211 | 157 | 53.478 (41.0-66.0) | 88.039 (73.0-104.0) |
| Li 2022 | 36 | fallback 34365 | 18 | 0.401 (0.0-2.0) | 0.970 (0.0-3.0) |
| Sadeghi 2022 | 75 | own 17309 | 26 | 1.642 (0.0-4.0) | 2.745 (0.0-6.0) |
| Bharadhwaj 2023 | 7 | own 21569 | 0 | 0.124 (0.0-1.0) | 0.234 (0.0-1.0) |
| Chen 2023 | 220 | own 1275 | 34 | 9.318 (4.0-15.0) | 14.729 (9.0-21.0) |
| Priyanka 2025 | 444 | own 25137 | 50 | 6.765 (2.0-12.0) | 13.288 (7.0-21.0) |
| Rebouças 2025 | 531 | own 13634 | 59 | 12.930 (7.0-20.0) | 22.372 (14.0-31.0) |
| Nakatsuka 2025 | 98 | fallback 34365 | 3 | 1.092 (0.0-3.0) | 2.640 (0.0-6.0) |
| Endresz 2026 (bioRxiv) | 4032 | fallback 34365 | 185 | 44.937 (33.0-57.0) | 108.612 (92.0-126.0) |

### 107 set
| meta-analysis | list genes in universe | universe | obs | uniform exp (95%) | DE-prior exp (95%) |
|---|---|---|---|---|---|
| Pérez-Santiago 2012 | 131 | own 13634 | 24 | 0.942 (0.0-3.0) | 1.800 (0.0-5.0) |
| Mistry 2013 | 98 | own 13634 | 8 | 0.704 (0.0-3.0) | 1.347 (0.0-4.0) |
| Bergon 2015 | 265 | own 21383 | 26 | 1.326 (0.0-4.0) | 2.775 (0.0-6.0) |
| Hess 2016 | 90 | own 21383 | 6 | 0.450 (0.0-2.0) | 0.943 (0.0-3.0) |
| Qin 2016 | 607 | own 13634 | 37 | 4.363 (1.0-9.0) | 8.343 (3.0-14.0) |
| Gandal 2018 (RNA-seq) | 4538 | own 21576 | 70 | 22.505 (15.0-31.0) | 45.104 (36.0-55.0) |
| Gandal 2018 (arrays) | 2028 | own 12211 | 50 | 15.611 (9.0-23.0) | 28.164 (20.0-37.0) |
| Li 2022 | 36 | fallback 34365 | 9 | 0.112 (0.0-1.0) | 0.301 (0.0-2.0) |
| Sadeghi 2022 | 75 | own 17309 | 16 | 0.455 (0.0-2.0) | 0.843 (0.0-3.0) |
| Bharadhwaj 2023 | 7 | own 21569 | 0 | 0.035 (0.0-1.0) | 0.073 (0.0-1.0) |
| Chen 2023 | 220 | own 1275 | 15 | 3.451 (1.0-7.0) | 5.427 (2.0-9.0) |
| Priyanka 2025 | 444 | own 25137 | 24 | 1.890 (0.0-5.0) | 4.123 (1.0-8.0) |
| Rebouças 2025 | 531 | own 13634 | 25 | 3.817 (1.0-8.0) | 7.298 (3.0-13.0) |
| Nakatsuka 2025 | 98 | fallback 34365 | 2 | 0.305 (0.0-2.0) | 0.819 (0.0-3.0) |
| Endresz 2026 (bioRxiv) | 4032 | fallback 34365 | 56 | 12.554 (6.0-19.0) | 33.697 (25.0-43.0) |

(`meta15_per_meta_analysis_hits.tsv` has set genes in universe, folds and P). Every list except Bharadhwaj 2023 (7 genes, 0 hits) exceeds both expectations; Nakatsuka 2025 exceeds the uniform null with P 0.10 for the 383 and 0.04 for the 107 and the DE-prior null by only a little (see table).

## Pooled hits (`meta15_pooled_hits.tsv`; sum of per-list hits)
| set | units | observed | uniform expected (95%) | DE-prior expected (95%) |
|---|---|---|---|---|
| 383 | all 15 | 989 | 237.7 (211-265) | 445.7 (412-480) |
| 383 | A 13 | 763 | 156.3 (134-179) | 298.2 (270-329) |
| 383 | B 13 | 801 | 191.7 (167-217) | 334.4 (306-364) |
| 107 | all 15 | 368 | 68.5 (54-83) | 141.1 (123-161) |
| 107 | A 13 | 296 | 45.7 (34-58) | 95.1 (79-112) |
| 107 | B 13 | 310 | 55.7 (43-69) | 106.5 (90-123) |

## Genes in 5 or more of the lists (`meta15_genes_by_k.tsv`, `meta15_genes_k_ge3.tsv`)
| units | genes in 5+ | in the 383 | in the 107 | expected 5+ (uniform / DE prior, 95% upper) | expected in 383 (uniform / DE prior) | expected in 107 (uniform / DE prior) |
|---|---|---|---|---|---|---|
| all 15 | 147 | 65 | 35 | 0.75 (3) / 12.6 (20) | 0.03 / 1.57 (4) | 0.01 / 0.61 (2) |
| A 13 | 80 | 42 | 25 | 0.10 / 1.85 | 0.00 / 0.25 | 0.00 / 0.10 |
| B 13 | 93 | 48 | 29 | 0.15 / 2.35 | 0.00 / 0.33 | 0.00 / 0.13 |

## Inputs overlap and coverage (`meta15_input_coverage.tsv`, `meta15_input_datasets.tsv`, `meta15_input_counts.tsv`, `meta15_pairwise_shared_inputs.tsv`, `meta15_link_check_vs_80.tsv`)
Adapted from task_b_meta_pooled.py (2026-10-12 classification, meta_link_evidence.tsv unchanged: 7 confirmed links and 1 unresolved, Bergon's Stanley study 18 to Aston 2004). All eight evidence nodes are among the 80 (Chu 2009 is a base study), so Chu 2009 counts as covered.
- Of 91 dataset slots (42 distinct datasets): 69 slots (27 distinct) are among the 80 base studies, 7 slots (3 distinct) are in the corpus but not among the 80, 15 slots (12 distinct) are outside the corpus. The datasets reach 25 distinct base studies, 21 contributing one of the 72 lists; 4 have an empty derived list (Iwamoto 2004, Iwamoto 2005, Maschietto 2015, Abdolmaleky 2019).
- Class (ii) reasons: no published SZ list of its own (2 datasets, Katsel 2005 and GSE92538 Gemma-derived primary), no decision-table entry (1, Glatt 2005). The Supplementary Table 5 (data needed) register is not on this host, so class (ii) is not split by Table 5 membership; Ketchesin 2023 and Pai 2019 (D33 rows) are among the 80 and no meta-analysis pools a dataset linked to them.
- Subject counts (cases/controls) are given only where recorded (the 80-study base table for studies among the 80, lineage node otherwise); unknown otherwise, none estimated; the number of datasets with unknown subjects is in the table. The 80-study base table counts are the per-study counts of the 81-base file (authors' counts before QC).
- Cohort clusters: per meta-analysis the count of its input cohort clusters and how many are shared with the clusters of the 80 studies and with another meta-analysis.

## Files
META15 tables `meta15_*.tsv`; scripts in `scripts/` (meta15_recurrent.py; task_b_core80.py); `in_*` copies of the 2026-10-12 evidence/audit; `sources_sha256.tsv`; SHA256SUMS. Run: nice, 1 thread, ~10 min. No files outside the private folder P10 and this landing folder were written.
