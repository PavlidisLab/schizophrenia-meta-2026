# Cohort clusters for the seven multi-bank studies

Run in a remote environment in a copy of the core-80 workspace (the stored inputs were left untouched).

## Rule
- A study that pools donors from several banks is a member of the cohort cluster of each of its banks (`multibank_cluster_sets.tsv`). The bank sets come from the lineage node `banks`, checked against the papers.
- A bank outside the named lanes (Charing Cross, Dallas, Miami, Sepulveda, Maudsley) gives that study its own one-study cluster.
- A spanning study never merges the clusters it touches.
- Cohort tier: the study's one list counts once in each of its clusters; a gene counts once per cluster. The Poisson-binomial across clusters treats clusters as independent although a spanning study sits in several, so the tier is slightly anticonservative.
- Sensitivity (`recurrence/sensitivity_own_cluster/`): each spanning study is one separate cluster of its own instead (counted once, independent of every bank).
- Pai 2019 has no gene list, so it enters no recurrence run; it is in the cluster counts over the 80 studies.
- None of the four unresolved papers reports results per bank (Mendizabal, Pai, Mendez-Victoriano, Saetre), so no per-bank lists exist and the one pooled list is used.

## Values
Cohort tier, scenario "published, brain" (`recurrence.tsv`):

| Quantity | Sets rule | Sensitivity (own cluster) |
|---|---|---|
| Studies in a cluster, of 80 | 80 | 80 |
| Clusters over the 80 studies | 22 | 25 |
| Lists in the cohort tier / clusters in the run (72 studies with a list) | 78 / 21 | 72 / 24 |
| Genes recurring with each cluster counted once, q < 0.1, uniform null | 140 | 117 |
| same, expression-matched null | 37 | 29 |
| same, DE-prior null | 83 | 48 |
| Genes in 5 or more clusters, observed | 1,414 | 1,150 |
| expected, uniform / DE prior / expression-matched | 275.3 / 1,026.0 / 839.1 | not tabulated here |
| Genes in 3 or more clusters, observed | 6,236 | |
| Study tier (72 lists), genes q < 0.1 | 383 | 383 |

Leave one study out, cohort tier, uniform null (`loso_sets_summary.tsv`, all 72 lists): 140 -> 88 without Voisey 2017, 46 without Xiao 2014, 79 without Maycox 2009, 100 without Zhu 2024, 107 without Saetre 2007; every other omission leaves 117 to 141.
Leave one cluster out (`loco_sets_summary.tsv`, 21 clusters of the run; the cluster is dropped from the cohort tier only): cmc+mssm+pitt+libd -> 42, harvard -> 44, stanley -> 46, Xiao's cluster -> 46, aus -> 60.

Distinct subjects (S80 run folder, spanning studies read as multi-bank): 1,858 to 4,401 (SZ 741 to 1,797; CTRL 864 to 2,151). By the rule of that table a spanning study adds nothing to the lower end because it may lie inside the clusters it draws from, and it counts at its reported N at the upper end. Only a paper-stated split could change it (Mendizabal Table S2 and Pai Supp. Data 1 give per-bank counts that do not match the counts in the corpus table: 29 SZ in Table S2 against 23 recorded), so none was applied.

Pairs of datasets with no subjects in common (`subject_relations/subject_relations.tsv`; bank sets added for the datasets of Maycox GSE17612, Pai GSE112523 and Mendizabal GSE107638, `code/patch3.py`):
- 192 bulk pairs of the 22 dataset-regions.
- 19 single-nucleus cohort pairs (they rest on donor data of the Endresz cohorts, not on banks).
- 26 bulk datasets of Extended Fig. 11a/14a: 229 pairs. Median rho(t) 0.1456; at 11 to 30 donors per group, 173 bulk pairs, median 0.146 against 0.163 for the 8 single-nucleus pairs, Mann-Whitney P 0.890 (`fig09_panel_a/`).

Fig 3b pairs of published lists (`h1/`): different-cohort pairs 2,080; same-cohort 476.
Published meta-analysis pairs (`h3/h3_pairs.tsv`, 246 pairs): 44 pairs share no cluster; of the 105 pairs of `fig02_e_meta_analysis_pairs`, 102 share a cluster and 3 do not.
