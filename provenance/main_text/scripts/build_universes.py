#!/usr/bin/env python3
"""Tested gene universe for every primary list (primary_lists.json, in_corpus studies).

    python3 build_universes.py   -> universes.json, universe_sets/*.txt

Each primary list gets a type:
  platform        array: the HGNC genes the GPL's probes map to (GEO .annot, or the GEO
                  platform table where no .annot exists), in platforms/
  fullstats       a table covering every tested gene (published, Gemma's, or our
                  reproduction): the gene set is its rows, with the list's threshold
                  filters removed and its selection filters (region, feature type) kept
  rnaseq_default  RNA-seq with no full table: the union of the RNA-seq full-statistics
                  tables named in RNASEQ_DEFAULT
  targeted        a custom panel; a gene set only when the panel's genes are known
  exclude         no usable gene universe (non-gene molecules, cell-type-stacked lists,
                  platform unknown)
n_tested is the number of genes the paper says it tested after its own filters, with the
wording it comes from. Only counts stated in genes are used as n_k in recurrence.py; probe,
transcript and mixed counts are recorded but not used.

Every symbol goes through extract.normalize_gene (the HGNC chokepoint).
"""
import copy
import csv
import datetime
import glob
import gzip
import hashlib
import json
import os

import extract

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")  # repo data/: lists, manifests, reference, raw, universe sets
SETS = os.path.join(DATA, "universe_sets")
PL = "platforms/"

# RNA-seq full-statistics tables whose union is the default RNA-seq universe (brain,
# gene level, every tested gene): list ids whose fullstats gene sets are pooled.
RNASEQ_DEFAULT = [
    "30545856__dge_scz",          # Gandal 2018 PsychENCODE, Table S1 DGE
    "38781321__scz_combined",     # Xia 2024, Data file S1
    "36652379__st2_fdr05",        # Clifton 2023, Supplementary Table 2
    "36319771__caudate_genes",    # Benjamin 2022, Supplementary Data 8 (gene rows)
    "32203495__table_s10_sz_fdr10",  # Jaffe 2020, Supplementary Table 10
    "34957416__deseq2_bulk",      # Jin 2021, CARseq_pipelines DESeq2 bulk table
    "36090817__scz",              # Sadeghi 2022, Supplementary Data 1 (Scz rows)
]

G, F, R, T, X = "platform", "fullstats", "rnaseq_default", "targeted", "exclude"

# list id: (type, source, n_tested, unit, where n_tested is stated, evidence)
# source: platform -> list of GPL files ("A|B" = intersection of A and B, "A+B" = union);
#         fullstats -> None (the list's own table) or a dict overriding the table spec;
#         targeted -> platform file or None (gene set not known); exclude -> reason
SPEC = {
 "11296301__table1": (G, "GPL80", None, None, None, "manifest platform GPL80 (HuGeneFL); 'over 6,000 human genes' is approximate, not a tested count"),
 "19358997__suppl_table": (G, "GPL570", None, None, None, "manifest platform GPL570; the list's probe IDs include 1560228_at, a U133 Plus 2.0-only probe set (lineage platform text says U133A)"),
 "19088852__endothelial_reproduced": (F, None, None, None, None, "our reproduction, reproduce/19088852/full_results_endothelial.tsv"),
 "27846841__up": (F, {"file": "raw/27846841/publisher/12918_2016_366_MOESM6_ESM.xls", "reader": "xls", "sheet": "Sheet1", "header_row": 2,
                      "columns": {"gene_reported": "symbol", "id_reported": "EnsemblID", "p": "pvalue"}, "filters": []},
                  14549, "transcripts", "manifest universe: '14,549 expressed transcripts (Additional file 6)'",
                  "Additional file 6 (MOESM6, 'Normalized counts of RNA-seq data and statistics for differential expression') has every expressed gene with pvalue and padj"),
 "27846841__down": (F, "=27846841__up", 14549, "transcripts", "manifest universe: '14,549 expressed transcripts (Additional file 6)'", "as 27846841__up"),
 "25658856__table2": (G, "GPL11532", None, None, None, "lineage geo GSE62333; GEO: GSE62333 is on GPL11532 (Human Gene 1.1 ST)"),
 "30077719__all_scz": (R, None, None, None, None, "RNA-seq (lineage: Illumina HiSeq2000); the supplement holds the DE rows only"),
 "31174959__dlpfc": (R, None, None, None, None, "RNA-seq; Table S11 holds genes at FDR < 10% only (632 DLPFC rows), not every tested gene"),
 "29740122__genes_fdr05": (R, None, 21312, "genes", "manifest universe: '21,312 Ensembl genes tested (plus 1,387 eRNAs)'", "RNA-seq; the supplement holds the DE rows only"),
 "16038679__tables2-4": (G, "GPL876", 9384, "genes", "manifest universe: '9,384 genes above the intensity cutoff in both cohorts'", "manifest platform GPL876 (Agilent Human 1 cDNA, layout B)"),
 "18985160__table3": (R, None, None, None, None, "mRNA sequencing (Illumina Genome Analyzer; GEO GSE12297 is GPL7110); '~33,200 RefSeq transcripts detected' is a detection count, not a gene count"),
 "19582768__san_diego_reproduced": (F, None, None, None, None, "our reproduction, reproduce/19582768/full_results.tsv (GPL5175 transcript clusters)"),
 "21727898__table2_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "22503327__suppl_table2": (X, "platform unknown", None, None, None, "main text not retrieved; the supplement does not name the array"),
 "22441714__s1_up": (G, "GPL6104", None, None, None, "manifest universe and paper: Illumina HumanRef-8 v2; GEO GPL6104 'Illumina humanRef-8 v2.0 expression beadchip'"),
 "22441714__s2_down": (G, "GPL6104", None, None, None, "as 22441714__s1_up"),
 "24287731__suppl1_s1": (G, "GPL6104", 15375, "genes and gene variants", "manifest universe: '15,375 after the intensity filter'",
                         "paper Methods: 'Illumina whole-genome HumanRef8 v2 BeadChip'; supplement Sequence Codes are ILMN_ IDs; GEO GPL6104"),
 "22558445__table_s1": (R, None, None, None, None, "RNA-seq (Illumina Genome Analyzer II); the table holds DE rows only"),
 "24169640__suppl_table4": (R, None, None, None, None, "RNA-seq (Illumina GA); the table holds DE rows only"),
 "23904455__table_s1": (G, "GPL10558", 9208, "genes", "paper Results: '14 588 transcripts (representing 9208 genes) were then selected for further analyses'",
                        "paper: 'Illumina HT12v4 microarrays'; GEO GPL10558"),
 "26183697__table2_qpcr": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "25985888__suppl_table3_array": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "25429046__table4_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "27650396__suppl_table1_lncrna": (X, "lncRNA list", None, None, None, "lncRNA IDs, no HGNC genes"),
 "28418402__table_s1": (R, None, 21146, "genes", "manifest universe: '21,146 genes with non-zero FPKM in >= 80% of subjects'", "RNA-seq (LCL); the table holds DE rows only"),
 "29391398__genes": (R, None, None, None, None, "RNA-seq (amygdala); the table holds p < 0.05 rows only"),
 "30911161__table1": (R, None, None, None, None, "imputed expression (S-PrediXcan) from CMC DLPFC and 12 GTEx predictors; the predictor gene sets are not held, so the RNA-seq default stands in; the stated 10,929 covers the CMC predictor only"),
 "31113746__scz_vs_hc": (F, None, None, None, None, "mmc2.xls sheet 'Scz vs HC' has 4,062 probes with adj.P.Val from 2.5e-07 to 0.999 (every probe tested, not only the DE rows); Illumina HT-12 v4"),
 "33558674__stable4_fdr05": (R, None, None, None, None, "RNA-seq; Supplementary Table 4 holds nominal p < 0.05 genes only (1,373)"),
 "33585819__p05_reproduced": (F, None, None, None, None, "our reproduction, reproduce/33585819/full_results.tsv"),
 "33677214__migrated_cins": (R, None, None, None, None, "RNA-seq (Kallisto); no full table"),
 "32152472__dg": (R, None, None, None, None, "RNA-seq; no full table"),
 "33288841__neurons": (R, None, None, None, None, "RNA-seq; no full table"),
 "36319771__caudate_genes": (F, None, None, None, None, "Supplementary Data 8 (MOESM10), Type == Gene rows: every tested gene"),
 "34936134__st_vs_ctrl": (R, None, 34783, "genes", "manifest universe: '34,783 genes (text: 8461/34,783 in the sex comparison)'", "RNA-seq; Supplementary Data 1 holds DE rows only"),
 "32542700__exo_circrna": (X, "circRNA list", None, None, None, "circRNAs"),
 "35701597__gn_suppl_table7": (R, None, None, None, None, "RNA-seq (Kallisto); no full table"),
 "36192459__mast_bbb_celltypes": (X, "cell-type-stacked list", None, None, None, "seven cell types stacked; per cell type only 80 to 565 pre-selected genes were tested (S. Table 6)"),
 "37117634__microglia_suppl_table2": (R, None, None, None, None, "RNA-seq (CeGaT); no full table"),
 "36302706__caudate": (F, None, 15041, "transcripts", "manifest universe: '15,041 transcripts tested (rows in Dataset S4 Caudate All)'", "Dataset S4 'Caudate All' has every tested transcript"),
 "36379214__olfactory_microarray_reproduced": (F, None, None, None, None, "our reproduction, reproduce/36379214/full_results.tsv (GPL570 probe sets)"),
 "37267675__mrna_table_s2": (R, None, None, None, None, "RNA-seq (StringTie); no full table"),
 "37661009__d150": (R, None, None, None, None, "RNA-seq; Table S5 holds DEGs only"),
 "39307399__ngn2": (R, None, None, None, None, "RNA-seq; Table S2 holds |log2FC| > 2, p < 0.05 rows only"),
 "39028640__gpc_reproduced": (F, None, None, None, None, "our reproduction, reproduce/39028640/full_results_gpc.tsv"),
 "41210297__text_4genes": (R, None, None, None, None, "AmpliSeq Transcriptome (Ion Torrent) targeted whole-transcriptome panel; the panel's gene list is not held, so the RNA-seq default stands in"),
 "41567988__bulk_table_s9": (R, None, None, None, None, "bulk RNA-seq; no full table"),
 "38704507__dg_neurons_reproduced": (F, None, None, None, None, "our reproduction, reproduce/38704507/full_results.tsv (HOMER gene rows)"),
 "15563509__fig2_sz_mito": (T, None, None, None, None, "only 676 mitochondria-related probe sets were tested (of 966 selected by GO keyword); the probe set list is in Supplementary Tables S1-S3, not retrieved"),
 "15892139__table1_top80": (G, "GPL96", None, None, None, "manifest platform GPL96 (HG-U133A)"),
 "16236523__tables2-3": (X, "platform unknown", 4311, "genes", "manifest universe: '4,311 expressed in >= 4 of 7 pairs'", "two-colour oligonucleotide arrays 'representing 19,000 gene transcripts'; no GEO record for the print"),
 "18778695__s1a_short_doi": (G, "GPL570", 14438, "probe sets", "manifest universe: '14,438 probe sets after presence, redundancy and variance filtering'", "manifest platform GPL570; lineage geo GSE21138 (GPL570)"),
 "20034564__sm2_infragranular": (G, "GPL570", None, None, None, "manifest platform GPL570; lineage platform text HG-U133 Plus 2.0"),
 "21091092__table_s1_down": (G, "GPL6104", None, None, None, "same data as Schmitt 2012 (manifest notes: every gene name and p value matches its tables), which used Illumina HumanRef-8 v2 (GPL6104)"),
 "21091092__table_s2_up": (G, "GPL6104", None, None, None, "as 21091092__table_s1_down"),
 "21490598__si_table3": (G, "GPL6244", None, None, None, "lineage geo GSE25673; GEO: GPL6244 (Human Gene 1.0 ST)"),
 "21738743__fig1_8mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "22393424__table_s2": (G, "GPL570", None, None, None, "manifest platform GPL570; lineage geo GSE27383 (GPL570)"),
 "23282246__table_s1": (R, None, 4715, "transcripts", "manifest universe: '4,715 RefSeq transcripts with FPKM > 5.0 in all 6 subjects'", "RNA-seq; Table S1 holds DE rows only"),
 "23164819__stable1_padj05": (F, None, 22278, "genes", "manifest universe: '22,278 human genes analyzed (text)'", "Supplementary Table 1 has every gene"),
 "24594777__stable1_s_vs_c": (R, None, 6314, "transcripts/genes", "manifest universe: '6,314 transcripts/genes expressed above background (of 22,075)'", "RNA-seq (SOLiD); Supplementary Table 1 has 141 genes over all comparisons only"),
 "25113377__table_s1_scz": (R, None, None, None, None, "RNA-seq; Table S1 holds q <= 0.1 rows only"),
 "25560755__tables2ab_top75": (G, "GPL13158", 19667, "probe sets", "manifest universe: '19,667 probe sets after expression and variance filtering of 54,715'", "lineage geo GSE93987; GEO: GPL13158 (HT HG-U133+ PM)"),
 "25981335__brain_reproduced": (F, None, None, None, None, "our reproduction, reproduce/25981335/full_results.tsv (GPL4133 probes)"),
 "27117414__sit1_mirna_npc": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "27668389__sdf3_genes": (R, None, 16423, "genes", "manifest universe: '16,423 genes (Ensembl models) expressed at levels sufficient for analysis'", "RNA-seq; Supplementary Data 3 holds DE genes only"),
 "27732689__etable3_diagnosis": (R, None, 13903, "genes", "manifest universe: '13,903 genes after normalization'", "RNA-seq; eTable 3 holds DE genes only"),
 "29263384__sd8_npc_fdr30": (F, None, 18910, "genes", "manifest universe: '18,910 genes (rows in Supplementary Data 8)'", "Supplementary Data 8 has every tested gene"),
 "30115913__s2_baseline_fdr05": (F, None, 21043, "genes", "manifest universe: '21,043 genes with RPKM > 0 in >= 50% of samples (rows in Table S2)'", "Table S2 has every tested gene"),
 "30664768__suppfig6f": (R, None, None, None, None, "RNA-seq; no full table"),
 "30781971__table2_circrna": (X, "circRNA list", None, None, None, "circRNAs"),
 "31123247__sz_pfc": (G, "GPL570", None, None, None, "lineage geo GSE53987 (GPL570); the supplement's Ingenuity export (3,173 genes) is not every tested gene"),
 "31959813__table_s2": (R, None, None, None, None, "RNA-seq (fibroblasts); no full table"),
 "32203495__table_s10_sz_fdr10": (F, None, 21460, "genes", "manifest universe: '21,460 genes (rows in Supplementary Table 10)'", "Supplementary Table 10 has every gene"),
 "15334603__table2": (G, "GPL91", 6452, "probe sets", "manifest universe: \"'approximately 12,000' transcripts; 6,452 met the robust-expression criteria\"", "manifest platform GPL91 (HG-U95A)"),
 "16406496__mmc5": (X, "platform unknown", None, None, None, "Compugen human 19K oligo set spotted at the Ramaciotti Centre; no GEO record for that print"),
 "18081029__table1": (G, "GPL96", None, None, None, "manifest platform GPL96 (HU133A)"),
 "14743183__scz_reproduced": (F, None, None, None, None, "our reproduction, reproduce/14743183/full_results.tsv, probe sets passing the presence and sex filters (sex_filter_pass)"),
 "18191109__suppl_table3_shared": (G, "GPL2891", None, None, None, "supplementary methods: Codelink UniSet Human 20K I Bioarray; GEO GPL2891 'CodeLink UniSet Human 20K I Bioarray'"),
 "15920292__qpcr_confirmed": (T, None, None, None, None, "3,200-clone cDNA array hybridized with pooled RNA; the clone list is not published"),
 "18445270__tables2-3": (X, "platform unknown", 8737, "genes", "manifest universe: '8,737 genes detected in >= 4 of 7 pairs'", "Compugen 19K two-colour oligonucleotide arrays (lineage); no GEO record for the print"),
 "19255580__cchpc_ba10": (G, "GPL570", 29464, "probe sets", "manifest universe: '29,464 after the intensity and intronic/antisense filters'", "manifest platform GPL570; lineage geo GSE17612 (GPL570)"),
 "18771015__table2_scz": (G, "GPL96+GPL97", None, None, None, "paper: HU133A and HU133B; union of GPL96 and GPL97"),
 "20695885__table_s1": (G, "GPL96+GPL97", 11483, "probe sets", "manifest universe: '11,051 to 11,904 probe sets per region passed the initial filter (Table S2; 11,483 in BA22)'", "probe IDs from the HG-U133 A and B arrays; union of GPL96 and GPL97"),
 "21111402__fig1b_mirna_ba46": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "21538462__scz_ba22": (G, "GPL570", 24715, "probe sets", "manifest universe: '24,715 after intensity/intron/antisense filtering'", "manifest platform GPL570; lineage geo GSE21935 (GPL570)"),
 "22315408__table1_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "22869038__suppl_table2": (R, None, None, None, None, "RNA-seq (SOLiD); Supplementary Table 2 holds DE transcripts only"),
 "23218666__table2": (G, "GPL6947", None, None, None, "manifest universe and lineage: Illumina HT-12 V3; GEO GPL6947"),
 "24507884__table1_splicing": (G, "GPL5175", None, None, None, "lineage: Affymetrix Human Exon 1.0 ST; GEO GPL5175 (transcript clusters)"),
 "24702465__table_s3": (G, "GPL1352", None, None, None, "lineage geo GSE37981; GEO: GPL1352 (Human X3P)"),
 "16604350__lcl_named_4q": (G, "GPL96", 11282, "UniGene clusters", "manifest universe: 'HG-U133A with a custom chip definition file: 11,282 UniGene clusters'", "HG-U133A (GPL96); the custom CDF is not available, so the stock platform stands in"),
 "20083392__text_sz_genes_up": (T, None, 19121, "probes", "manifest universe: '19,121 filtered probes'", "platform probe list not retrieved (batch 11)"),
 "20083392__fig3b_ann_predictors": (T, None, 19121, "probes", "manifest universe: '19,121 filtered probes'", "platform probe list not retrieved (batch 11)"),
 "22161352__text_shared_sz_bd": (T, None, None, None, None, "reanalysis list (batch 11); universe not stated"),
 "35339428__table2_hub_genes": (T, None, None, None, None, "reanalysis list (batch 11); universe not stated"),
 "42623758__text_top_genes": (T, None, None, None, None, "reanalysis list (batch 11); universe not stated"),
 "38448582__snap_a_top100": (F, None, None, None, None, "every gene with a SNAP-a loading (Ling 2024 Supplementary Table, gene loadings)"),
 "38448582__snap_n_top100": (F, None, None, None, None, "every gene with a SNAP-n loading (Ling 2024 Supplementary Table, gene loadings)"),
 "24628518__supp_table3": (G, "GPL1352", None, None, None, "GSE46509; Human X3P (GPL1352), as the sister pyramidal-neuron paper"),
 "15563509__supp3_ph64_sz": (T, None, None, None, None, "only 676 mitochondria-related probe sets were tested (as 15563509__fig2_sz_mito)"),
 "24686136__si_table3": (G, "GPL6244", None, None, None, "lineage geo GSE40102; GEO: GPL6244 (Human Gene 1.0 ST)"),
 "26173148__de_mirnas": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "26733343__table_s1": (X, "platform annotation has no genes", 8518, "genes", "manifest universe: '8,518 genes after intensity filtering'",
                        "lineage geo GSE54913 is on GPL15314 (Arraystar Human LncRNA V2.0), but the GEO platform table carries no gene symbol for its mRNA probes (4 of 60,756 rows have an Ensembl ID); GPL13825, the same design, has none either"),
 "27727244__suppl_table1": (G, "GPL570", None, None, None, "manifest platform GPL570; lineage geo GSE73129 (GPL570)"),
 "28094170__ncc_rnaseq_reproduced": (F, {"filters": [{"col": "expressed_all8", "equals": "yes"}]}, None, None, None,
                                     "the authors' Cuffdiff output (reproduce/28094170/full_results.tsv), genes expressed in all eight samples (expressed_all8)"),
 "28736215__table_s2": (R, None, None, None, None, "RNA-seq; no full table"),
 "28809853__table_s2_all3": (R, None, None, None, None, "RNA-seq; no full table"),
 "28809853__table_s2_any_method": (R, None, None, None, None, "RNA-seq; no full table (as table_s2_all3)"),
 "28809853__table_s4_761_nominal": (R, None, None, None, None, "RNA-seq; no full table (as table_s2_all3)"),
 "30050107__table_s9_replicated": (R, None, 24122, "genes", "manifest universe: '24,122 genes, 420,022 exons, 61,950 transcripts ... tested'", "RNA-seq; full statistics are on eqtl.brainseq.org, not retrieved"),
 "30050107__discovery_fdr10_qsva": (F, {"file": "raw/30050107/brainseq/expressed_de_Gene.tsv.gz", "reader": "tsv", "header_row": 1,
                                      "columns": {"gene_reported": "Symbol", "id_reported": "EnsemblGeneID", "effect": "log2FC_qsva", "p": "pval_qsva", "fdr": "fdr_qsva"},
                                      "filters": []},
                                  24122, "genes", "manifest universe: '24,122 genes ... tested' (Methods: FDR over the expressed features of this sample set); same count as the table's Gene rows",
                                  "gene-level rows of the authors' public 'SCZD vs Control: All Expressed Features' table (eqtl.brainseq.org/phase1, expressed_de_features.rda; derived expressed_de_Gene.tsv.gz): 24,122 expressed genes, qSVA statistics"),
 "30870714__table_s7_all_genes": (R, None, None, None, None, "RNA-seq; the sheet holds FDR-significant rows only"),
 "31288836__neun_szdeg": (R, None, None, None, None, "RNA-seq (sorted nuclei); no full table"),
 "31540999__glut_neurons_reproduced": (F, None, None, None, None, "our reproduction, reproduce/31540999/full_results_glutamatergic.tsv (Clariom S transcript clusters)"),
 "31019265__suppl_table1_cins": (R, None, None, None, None, "RNA-seq; no full table"),
 "32186681__organoids": (R, None, None, None, None, "RNA-seq (Cuffdiff); the sheet holds significant genes only"),
 "32566292__an_sz_nc": (G, "GPL570", None, None, None, "manifest platform GPL570; lineage platform text HG-U133 Plus 2.0"),
 "28641886__pbmc": (G, "GPL10558", 11818, "genes", "paper Methods: 'Probes with a detection p-value >0.01 were removed (22,599), resulting in 11,818 genes for further analysis'",
                    "paper Methods: 'Human HT-12 v4 Expression BeadChip Array (Illumina)'; GEO GPL10558"),
 "16223876__si_table4_dlpfc": (G, "GPL96", None, None, None, "manifest notes: the probe IDs in Table 4 are HG-U133A IDs; GPL96"),
 "23442539__table_s2_before_treatment": (G, "GPL6947", 10207, "genes", "paper Results: '48803 transcripts (comprising 25202 genes) from PBMC revealed 10207 genes as being expressed'",
                                         "paper Methods: 'Illumina HT-12_V3 beadchips'; GEO GPL6947"),
 "18184693__sam_2mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "11576761__drug_treated": (T, None, None, None, None, "custom NIA neuroarray of 1,128 genes; its gene list is not held"),
 "12363385__table2_10pct": (T, None, None, None, None, "custom NIA neuroarray of 1,127 genes; its gene list is not held"),
 "14741323__table3": (T, None, None, None, None, "custom 1,128-gene Neuroarray; its gene list is not held"),
 "16414245__table2_brain": (T, None, None, None, None, "custom 6,000-gene cDNA array (Ramaciotti Centre); its gene list is not held"),
 "17568569__table2": (T, None, None, None, None, "custom NimbleGen array of >1,800 candidate genes; its gene list is not held"),
 "17326821__table1_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 # batch 12 (2026-09-28): six studies added by hand (manifests/batch12_added_2026_09_28.json)
 "38260577__table_s3_diffgenes": (R, None, None, None, None, "iPSC-derived neuron RNA-seq (DESeq2); Table S3 holds DE rows only; the tested gene count is not stated"),
 "38260577__table_s3_diffmicrornas": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "38260577__table_s3_dapa_ipsc": (T, None, None, None, None, "dAPA over genes with more than one polyadenylation site in iNs (about 4,400 by the paper's 14.29%); that gene set is not published"),
 "38260577__table_s3_dapa_postmortem": (T, None, None, None, None, "genes with two 3'UTR regions in the iN polyadenylation-site library; that gene set is not published"),
 "31642882__s2_homogenate_adjp05": (T, "own_table", 402, "proteins", "manifest universe: '402 protein targets passing QC in homogenate (478 assayed by LC-SRM; preprint Table S18)'",
                                    "LC-SRM targeted panel (478 assayed, under the 500 line); Supplemental Table 2 lists the 402 targets passing QC"),
 "31642882__s3_synaptosome_adjp05": (T, "own_table", 155, "proteins", "manifest universe: '155 synaptosome-enriched protein targets analysed (348 quantified in synaptosomes of 478 assayed by LC-SRM)'",
                                     "LC-SRM targeted panel (478 assayed); Supplemental Table 3 lists the 155 synaptosome targets analysed"),
 "25433904__table2_p01": (T, None, 155, "proteins", "manifest universe: '155 selected synaptic protein targets (223 peptides from 155 proteins passed QC of 289 peptides from 180 proteins assayed)'",
                          "LC-SRM targeted panel; the 155 proteins are not listed"),
 "42748139__sd01_isoforms_lfsr05_all_chains": (F, None, 7271, "isoform groups", "manifest universe: '7,271 isoforms or isoform groups (from 9,169 transcripts matched to 79,040 peptides) in Dataset S1'",
                                               "Dataset S1 has every isoform or isoform group"),
 "28485405__supp_table3b_dlpfc_fdr05": (F, None, 561, "genes", "manifest universe: '561 immune genes in 20 pathways (Supplementary Tables 1 and 3B)'",
                                        "Supplementary Table 3 lists all 561 genes; above the 500-gene panel line (Paul, 2026-09-28), so not a panel"),
 # batch 13 (2026-09-28): the two 'not attempted' corpus studies (manifests/batch13_notattempted_2026_09_28.json)
 "32513544__supp2_down_p01": (R, None, 13849, "genes", "manifest universe: '13,849 expressed genes (CPM > 5 in at least two samples; Methods)'",
                              "DLPFC RNA-seq (edgeR); supplementary data 2 holds the p < 0.01 rows only"),
 "32513544__supp2_up_p01": (R, None, 13849, "genes", "manifest universe: '13,849 expressed genes (CPM > 5 in at least two samples; Methods)'",
                            "as 32513544__supp2_down_p01"),
 "41611717__sd2_layers123_fdr05": (F, None, 5243, "protein groups", "manifest universe: '5,243 protein groups quantified (Supplementary Data 2, every row)'",
                                   "Supplementary Data 2 has every quantified protein group with its authors' HGNC symbol and both Layer sets' statistics"),
 "41611717__sd2_layers456_fdr05": (F, None, 5243, "protein groups", "manifest universe: '5,243 protein groups quantified (Supplementary Data 2, every row)'",
                                   "as 41611717__sd2_layers123_fdr05"),
 # batch 14 (2026-09-29): proteomic studies the proteomics side study found missing (manifests/batch14_proteomics_2026_09_29.json)
 "16491132__english2011_s1": (X, "proteomics list", None, None, None, "2-DE; spots matched not in the material read; list is English 2011's compilation"),
 "21136665__english2011_s2": (X, "proteomics list", 423, "spots", "abstract: 'From the 423 matched spots between the two groups'",
                              "2-DE spots, not proteins; list is English 2011's compilation"),
 "21136626__english2011_s2": (X, "proteomics list", None, None, None, "2-DE; spot count not in the material read; list is English 2011's compilation"),
 "19003868__english2011_s1": (X, "proteomics list", None, None, None, "2-D DIGE; spot count not in the material read; list is English 2011's compilation"),
 "19562803__english2011_s2": (X, "proteomics list", None, None, None, "2-D DIGE; spot count not in the material read; list is English 2011's compilation"),
 "18268500__suppt1_dige_scz": (X, "proteomics list", None, None, None, "2D-DIGE of pooled samples; spot count not in the material read"),
 "18268500__suppt3_gelc_scz": (X, "proteomics list", 56, "proteins", "Supplementary Table 3 'Spectral counts for all proteins identified by GeLC/MS-MS': 56 rows",
                               "the 56 identified proteins are in the transcribed image only as far as the 17 listed rows"),
 "19110265__table2": (X, "proteomics list", 750, "spots", "Results: 'On average, 750 spots were detected in each 2-DE'", "2-DE spots, not proteins"),
 "25048004__st3_p05": (F, None, 727, "proteins", "Results: '... left 727 proteins' (Supplementary Table 3, every row)",
                       "Supplementary Table 3 has every tested protein with its ANCOVA p and FDR"),
 "25048004__st3_fdr05": (F, None, 727, "proteins", "Results: '... left 727 proteins' (Supplementary Table 3, every row)", "as 25048004__st3_p05"),
 "25728835__supp_ph6_11_scz": (X, "proteomics list", None, None, None, "pH 6-11 2D-DIGE; spot count not in the material read"),
 "19721432__suppl_table2_dlpfc": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "20675101__suppl_table7_sz": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "22675524__atlas_array_pooled": (T, "GPL138", None, None, None, "BD Clontech Atlas Human Neurobiology Array (593 genes); GEO GPL138 'Atlas Human Neurobiology Array (Cat. #7736-1)', genes by Entrez ID"),
 "18992145__additional_file1": (G, "GPL96+GPL570+GPL8300", 19502, "genes", "manifest universe: '19,502 genes across HG-U133A, HG-U133 Plus 2.0 and HG-U95Av2'", "union of the three named platforms"),
 "22954356__st2": (G, "GPL96", None, None, None, "HG-U133A and HG-U133 Plus 2.0 studies combined; only probe sets on both (the HG-U133A set) can be combined, so GPL96"),
 "22212594__meta_up": (G, "GPL96", None, None, None, "manifest platform GPL96; lineage: 'Affymetrix HG-U133A (primary analysis platform)'"),
 "22212594__meta_down": (G, "GPL96", None, None, None, "as 22212594__meta_up"),
 "24569695__suppT2_hiSCZ": (G, "GPL570", None, None, None, "manifest platform GPL570; manifest universe '54,675 probe sets (three HG-U133A studies and one HG-U133 Plus 2.0 study)'"),
 "26285829__brain_up": (G, "GPL570", 8655, "genes", "manifest universe: '8,655 genes common to the 10 brain datasets'",
                        "the genes common to the 10 datasets are a subset of each; four are GEO series on GPL570 (GSE17612, GSE21138, GSE21935, GSE53987), so GPL570 is a superset"),
 "26285829__brain_down": (G, "GPL570", 8655, "genes", "manifest universe: '8,655 genes common to the 10 brain datasets'", "as 26285829__brain_up"),
 "27219343__scz_z2.2": (F, None, None, None, None, "Supplementary Table 3: 48,197 Ensembl rows, of which 8,856 carry a SCZ Z-score (the rows kept)"),
 "27450777__brain_s3a": (G, "GPL570", 20767, "genes", "manifest universe: '20,767 genes'",
                         "PFC mega-analysis of Affymetrix HG-U133 family studies (lineage platform text; Lanz GSE53987, Maycox GSE17612 and Narayan GSE21138 are GPL570); GPL570 as the superset"),
 "30545856__dge_scz": (F, None, None, None, None, "Table S1 DGE has every gene"),
 "29632383__twas_s3": (R, None, None, None, None, "imputed expression (TWAS; CMC DLPFC, NTR and YFS blood, METSIM adipose weights); the weight gene sets are not held, so the RNA-seq default stands in"),
 "31073119__tableB1_all": (G, "GPL10558", None, None, None, "lineage platform text: 'Illumina HumanHT-12 v4 expression array (NIMH dataset)'; GEO GPL10558"),
 "30576541__genemeta_s4": (F, None, 2641, "genes", "manifest universe: '2,641 genes common to the three datasets'", "Supplementary Table S4 has every gene"),
 "34969953__table_s1": (G, "GPL570", None, None, None, "manifest platform GPL570; contrast names GEO GSE73129 (GPL570)"),
 "34010744__cotwin_neurons": (F, None, None, None, None, "the authors' DESeq2 table has every Ensembl gene (rows with a p value kept)"),
 "36932057__scz_dge_human": (R, None, None, None, None, "RNA-seq; Table S3 holds P < 0.05 genes only"),
 "22373040__af2_scz": (G, "GPL8300", None, None, None, "manifest universe: 'HG-U95Av2 probes (GPL8300)'; contrast names GEO GSE12654 (GPL8300)"),
 "30976116__s1_age40to60_replicated": (R, None, None, None, None, "RNA-seq (CMC and psychENCODE); no full table"),
 "33562405__table_s1": (G, "GPL5175|GPL570", 11112, "genes", "manifest universe: '11,112 genes shared among the datasets'", "meta-analysis of GSE18312 (GPL5175) and GSE27383 (GPL570): intersection of the two platforms"),
 "33805976__st2_p05": (F, None, 27100, "genes and lncRNAs", "manifest universe: '27,100 genes and lncRNAs in Supplementary Table 2'", "Supplementary Table 2 has every gene"),
 "34952924__s3_ba46": (G, "GPL570", None, None, None, "manifest platform GPL570; contrast names GEO GSE53987 (GPL570)"),
 "36652379__st2_fdr05": (F, None, 16278, "genes", "manifest universe: '16,278 protein-coding genes expressed at > 5 CPM in >= 10 samples'", "Supplementary Table 2 has every gene"),
 "38377413__deseq2_scz": (R, None, 14753, "genes", "manifest universe: '14,753 genes after low-count filtering'", "RNA-seq (GSE80655); Supplementary Table 1 holds DE genes only"),
 "38730231__male_dlpfc": (R, None, None, None, None, "RNA-seq (BrainSeq); Data S10 holds FDR < 0.05 features only"),
 "38843579__scz_caudate": (R, None, None, None, None, "RNA-seq; the sheets hold DE transcripts only"),
 "39972407__intersection": (G, "GPL570|GPL6947|GPL6883", None, None, None, "genes DE in all three of GSE27383 (GPL570), GSE38484 (GPL6947) and GSE38481 (GPL6883); intersection of the three platforms"),
 "23554570__gse12649_ba46": (G, "GPL96", None, None, None, "manifest platform GPL96; GEO GSE12649 (GPL96)"),
 "31399567__night_corrp05": (F, None, 13910, "genes", "manifest universe: '13,910 genes'", "Supplementary Data 9 has every gene"),
 "31920755__table1_all3regions": (G, "GPL570", None, None, None, "manifest universe: HG-U133 Plus 2.0 (GSE53987); GEO GSE53987 (GPL570)"),
 "35855328__fig1a": (R, None, None, None, None, "RNA-seq (GSE121376); no table published"),
 "35983039__pfc": (G, "GPL570", None, None, None, "contrast names GEO GSE53987 (GPL570)"),
 "36162734__st3_genes": (R, None, None, None, None, "RNA-seq (ACC); no full table"),
 "36405929__demrna": (G, "GPL570", None, None, None, "contrast names GEO GSE17612 (GPL570)"),
 "36462601__meta_sz": (G, "GPL96+GPL570+GPL13158", None, None, None, "vote counting over GSE12649 (GPL96), GSE21138, GSE21935, GSE53987 (GPL570) and GSE93987 (GPL13158); union of the three platforms"),
 "38145985__tableS1_genes": (R, None, 19383, "genes", "manifest universe: '19,383 expressed genes in the joint analysis (Results)'", "RNA-seq (CMC); Supplementary Table 1 holds FDR < 5% rows only"),
 "40727846__tableS3_scz": (R, None, None, None, None, "RNA-seq (pooled public brain RNA-seq); Table S3 holds padj < 0.05 rows only"),
 "41422210__suppl5_degs": (G, "GPL96|GPL570", 10829, "genes", "manifest universe: '10,829 genes common to the 4 discovery datasets'",
                           "GSE17612 and GSE53987 (GPL570), GSE12649 (GPL96), GSE92538 (GPL17027 and GPL10526, HG-U133A and Plus 2 with Entrez CDFs): intersection of GPL96 and GPL570"),
 "41714609__tableS4_nominal": (F, None, 18677, "genes", "manifest universe: '18,677 genes in Table S4'", "Supplementary Table 4 has every gene"),
 "42463512__deg_dlpfc": (F, None, 22542, "genes", "manifest universe: '22,542 genes in the region's network'", "Supplementary Data 16, DLPFC rows: every gene"),
 "42693176__tableS16_celltypes": (X, "cell-type-stacked list", None, None, None, "20 snRNA-seq cell types stacked; no per-gene tested universe published"),
 "38781321__scz_combined": (F, None, None, None, None, "Data file S1 has every gene"),
 "25522158__sczcgs": (R, None, None, None, None, "RNA-seq reanalysis of Wu 2012; Additional file 3 keys transcripts by RefSeq accession only"),
 "26818902__table_s2": (G, "GPL96", 22277, "probe sets", "manifest universe: '22,277 HG-U133A probe sets'", "manifest platform GPL96"),
 "29439242__scz_meta_fdr05": (F, None, 12384, "genes", "manifest universe: '12,384 genes with SCZ statistics (26,647 rows in the table)'", "Data Table S1: rows with a SCZ P value"),
 "34957416__deseq2_bulk": (F, None, 20788, "genes", "manifest universe: '20,788 genes'", "the authors' DESeq2 bulk table has every gene"),
 "36090817__scz": (F, None, 18480, "genes", "manifest universe: '18,480 genes with Scz statistics'", "Supplementary Data 1, Scz rows: every gene"),
 "36424395__scz_dirg_fdr05": (T, "own_table", 1275, "genes", "manifest universe: '1,275 curated immune-related genes (IRGs)'", "restricted to 1,275 curated immune-related genes; Table S3 has all of them"),
 "37031222__ba10_meta_fdr10": (F, None, 8398, "genes", "manifest universe: '8,398 genes (Entrez) in the meta-analysis'", "Supplementary Table S8 has every gene"),
 "37031222__qpcr_sz_fdr10": (T, "own_table", 111, "genes", "manifest universe: '111 target genes (Supplementary Table S5, MOESM6)'",
                             "custom Taqman qPCR cards (115 targets); Supplementary Table S5 has all 111 target genes with SCHIZ statistics"),
 "39727940__pfc_cash001": (G, "GPL570", None, None, None, "manifest platform GPL570; GSE53987 (GPL570)"),
 "40796563__sumrank_scz_broad": (X, "cell-type-stacked list", 17420, "genes", "manifest universe: '17,420 genes common to the 3 SCZ datasets (Methods), per cell type'", "seven snRNA-seq cell types stacked"),
 "41443001__scz_vs_control_brain": (G, "GPL570", None, None, None, "contrast names GEO GSE53987 (GPL570)"),
 "41748837__scz_edger_adjp05": (F, None, 63145, "Ensembl genes", "manifest universe/notes: '63,145 Ensembl genes'", "MOESM1 sheet 'DEG result of Schizophrenia' has every gene"),
 "19245676__supp_t3_early_scz": (G, "GPL570", None, None, None, "manifest platform GPL570; reanalysis of GSE21138 (GPL570)"),
 # batch 10 (2026-09-24): genome-wide sweep studies
 "27801890__supp_table3_sz": (G, "GPL5175", None, None, None, "manifest: Affymetrix Human Exon 1.0 ST, core transcripts at gene level; GPL5175 transcript clusters"),
 "29463818__supp_table1_ba10_interactome": (G, "GPL5175", None, None, None, "manifest: Affymetrix Human Exon 1.0 ST (the list is the IPA-interactome subset of the BA 10 DE genes)"),
 "34440415__table_s2_padj05": (F, None, None, None, None, "Table S2 has every gene passing the count filter (23,640 rows; 2,292 with padj NA)"),
 "27898074__cp_ac_supp_table3": (R, None, None, None, None, "poly-A RNA-seq (HiSeq 2000); Supplementary Table 3 holds DE rows only"),
 "30786269__total_rna_mrna_fdr05": (R, None, None, None, None, "total RNA-seq; MOESM5 holds DE rows only; number tested not stated"),
 "34059796__sez_supp_table_s9_q05": (R, None, 13387, "genes", "manifest universe: '13,387 genes (of 60,433 mappable) with >= 5 CPM in at least 10 cases'", "RNA-seq; Table S9 holds DE rows only"),
 "39526318__amygdala_scz_p05": (F, None, 15940, "genes", "manifest universe: '15,940 genes after filtering'", "preprint supplement media-1.xlsx has full DESeq2 statistics"),
 "29112193__suppl_table3_pv_cells": (G, "GPL13667", None, None, None, "manifest universe: Affymetrix HG-U219 (GPL13667); the paper's expression-filtered count is not stated"),
 "34272489__superficial_p05": (G, "GPL6244", None, None, None, "Affymetrix Human Gene 1.0 ST (GSE145554 is on GPL6244)"),
 "39237723__scz_celltypes_fdr20": (X, "cell-type-stacked list", None, None, None, "5 LCM neuron types stacked; per-type universes of about 24,000 to 26,000 genes"),
 "41174894__hb_bulk_fdr10": (F, None, 22756, "genes", "manifest universe: '22,756 genes tested (rows of ETable 7)'", "ETable 7 has every tested gene"),
 "41757081__layer_adjusted_fdr10": (F, None, 13977, "genes", "manifest universe: '13,977 genes detected by Visium (rows of Supplementary Table 4)'", "Supplementary Table 4 has every tested gene"),
 "31431609__table1_mirna_fdr05": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "19772658__suppl_table3_dx_by_probeset": (X, "exon-level probeset interaction list", None, None, None, "manifest modality 'other statistic'; not a gene-level DE list"),
 "37147304__table_s3_scz_hc_fdr05": (F, None, None, None, None, "Table S3 sheet SCZ_HC has every probe (47,323 rows)"),
 "42033308__table_s2_ssd_v1": (R, None, None, None, None, "whole-blood RNA-seq; Table S2 holds DE rows only"),
 "37861850__suppmat_table3_scz_fdr10": (F, None, 10710, "genes", "manifest universe: 'text: 10,710 genes after filtering < 10 counts' (11,627 Ensembl genes in the table)", "Supplementary Table 3 has every tested gene"),
 "28822801__suppl_table1_lcl_p01": (X, "platform unknown", None, None, None, "main text paywalled; the retrieved supplements do not name the array"),
 "34630179__suppl_table_s3_pbmc": (R, None, 8911, "transcripts", "manifest universe: '8,911 transcripts with FPKM > 1 in more than 80% of samples in all four data sets'", "PBMC RNA-seq; Table S3 holds DE rows only"),
 "37149105__table_s5_deseq2": (R, None, None, None, None, "monocyte RNA-seq; Supplementary Table 5 holds DE rows only"),
 "38988837__table_s7_degs": (R, None, None, None, None, "leukocyte RNA-seq; Table S7 holds DE rows only"),
 "41125054__table_s1_fep_vs_hc": (X, "miRNA list", None, None, None, "serum EV miRNAs"),
 "40267851__table2_q05": (R, None, 13861, "genes", "manifest universe: '13,861 genes with non-zero expression in at least 90% of samples'", "whole-blood RNA-seq; Table S1 holds genes at q < 0.10 only"),
 "41788806__table_s13_ssd_on_ap": (R, None, None, None, None, "whole-blood RNA-seq; Table S13 holds DE rows only"),
 "37605121__mrna_fc2_p05": (R, None, None, None, None, "plasma exosome total RNA-seq; Table S1 holds DE rows only"),
 "doi_10.1038_s44220-024-00306-1__table_s2_fdr05": (R, None, None, None, None, "Cui 2024 blood RNA-seq; Table S2 holds DE rows only (added 2026-09-24)"),
 "40770299__mirna_fc2_p05": (X, "miRNA list", None, None, None, "plasma exosome miRNAs"),
 "30786269__circrna_enriched_p05": (X, "circRNA list", None, None, None, "circRNAs by host gene (circRNA-enriched RNA-seq)"),
 "17464717__table4_rfdd": (X, "platform gene set unknown", None, None, None, "restriction fragment differential display of about 12,500 genes; no gene-level platform table"),
 "18226501__table2_down_8p_13q_22q": (X, "locus-filtered subset", None, None, None, "only the down-regulated genes on 8p, 13q and 22q are printed (6 of 64)"),
 "28822801__table2_lcl_fdr05": (G, "GPL13667", None, None, None, "Methods: Affymetrix HG-U219 (GPL13667)"),
 "40015618__table_s2_degs": (R, None, None, None, None, "leukocyte RNA-seq (NovaSeq); mmc8 holds DE rows only"),
 "35660966__rra_fig2": (X, "meta-analysis rank list", None, None, None, "RRA over 10 datasets on different platforms; no single tested universe (group 3 reanalysis)"),
 "36125368__m2_male_vs_control": (R, None, None, None, None, "RNA-seq reanalysis of GSE107638 (NeuN+ nuclei); SDC6 holds DE rows only"),
 "31157371__signature_table_s2": (X, "classifier features", None, None, None, "RFE-SVM features (manifest modality 'gene selection')"),
 "26697997__additional_file3": (X, "platform unknown", None, None, None, "Illumina BeadChip; version not stated in the retrieved material (main text paywalled)"),
 # batch 15 (2026-09-30): the 46 studies of the 2026-09-30 genome-wide round (manifests/batch15_genomewide_2026_09_30.json;
 # primary lists from tools/b15_primary_lists.py; proteomic lists take 'exclude' in main())
 "39397771__table_s13_clusters": (X, "cell-type-stacked list", None, None, None, "snRNA-seq neuronal clusters stacked, genes pre-selected per cluster (Puvogel 2022 criteria), as 36192459"),
 "39397771__table_s14_subclusters": (X, "cell-type-stacked list", None, None, None, "as 39397771__table_s13_clusters, per subcluster"),
 "37498306__table_s1": (R, None, None, None, None, "BA10 RNA-seq (Cuffdiff); Table S1 holds DE rows only; the tested gene count is not stated"),
 "27622934__table2_hippocampus_sz": (R, None, 21861, "genes", "manifest universe: '21,861 genes with detected expression in the hippocampus (Results)'",
                                     "hippocampus RNA-seq (DESeq2); Table 2 holds the significant genes only"),
 "30053206__supp_table5_te_signature": (X, "transposable-element loci", None, None, None, "TE loci (modality gene selection), not genes"),
 "33184259__figshare_noiseq_prob08": (F, None, None, None, None, "the authors' figshare NOISeq table has every tested gene (18,226 rows) with its probability"),
 "28123834__table2_both_regions": (G, "GPL5175", None, None, None,
                                   "the supplementary expression file's 21,980 identifiers are all Human Exon 1.0 ST transcript clusters (GPL5175; "
                                   "17,868 of category 'main', matching the Methods' 'over 17,868 NCBI Reference Sequence (RefSeq) transcripts'), "
                                   "none is an HG-U133 Plus 2.0 probe set, though the Methods name that array; no GEO series (PubMed-GEO link "
                                   "and GEO author search, 2026-09-30); as 27801890"),
 "28123834__table2_dlpfc_both_hemispheres": (G, "GPL5175", None, None, None, "as 28123834__table2_both_regions"),
 "28123834__table2_mofc_both_hemispheres": (G, "GPL5175", None, None, None, "as 28123834__table2_both_regions"),
 "28123834__table3_dlpfc_left": (G, "GPL5175", None, None, None, "as 28123834__table2_both_regions"),
 "28123834__table4_mofc_left": (G, "GPL5175", None, None, None, "as 28123834__table2_both_regions"),
 "28123834__table4_mofc_right": (G, "GPL5175", None, None, None, "as 28123834__table2_both_regions"),
 "41639050__s3a_sncrna_hbcc_fdr05": (X, "sncRNA list", None, None, None, "small RNAs (miRNAs, isomiRs, tRNA and other fragments), not genes"),
 "41639050__s3d_sncrna_mssm_fdr05": (X, "sncRNA list", None, None, None, "as 41639050__s3a_sncrna_hbcc_fdr05"),
 "38239571__table_s3_dels": (X, "lncRNA list", 768, "lincRNAs", "manifest universe: '861 lincRNAs detected; 768 entered the ANOVA'", "lincRNAs, as 27650396"),
 "27916288__sfvscf_red": (X, "sncRNA list", None, None, None, "small RNAs (miRNAs, piRNAs, snoRNAs, snRNAs), female-only selection"),
 "32006661__table_s1_degs": (G, "GPL6947", None, None, None, "manifest: Illumina HumanHT-12 v3 (GPL6947); the tested count is not stated"),
 "41942462__sd3a_dopa_fdr05": (F, None, None, None, None, "Supplementary Data 3 sheet A has every tested gene (30,828 rows) with its statistics"),
 "41942462__remacorr_sd5a": (X, "transcript-usage list", None, None, None, "remaCorr transcript-usage statistic over GENCODE transcripts (modality other statistic), not a gene-level DE test (decided 2026-09-30: exclude)"),
 "24475125__table1_mirna_p05": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "24475125__table3_mirna_fc40": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "27940260__suppl_table2_fdr05": (R, None, None, None, None, "frontal cortex RNA-seq (TruSeq RNA Access, edgeR); Supplementary Table 2 holds DE rows only"),
 "16103141__fig3_herv_scz_vs_control": (X, "HERV element", None, None, None, "one HERV pol element on a custom 50-probe chip; no HGNC gene"),
 "21183010__table2_mirna_dx_pp95": (X, "miRNA list", None, None, None, "mature miRNAs and snoRNAs (targeted TaqMan panel)"),
 "18683247__fig1b_top12_short_doi": (T, None, 783, "genes", "manifest universe: '783 genes passing the Present-call filter (Methods, journal p. 280) of 1,900 transcripts on the GLYCOv2 array'",
                                     "custom GLYCOv2 glycogene array; its gene list is not held"),
 "18683247__fig3a_qpcr_short_doi": (T, None, 7, "genes", "manifest universe: '7 genes tested by qPCR (the six listed and ASAH1, p = 0.09)'", "qPCR panel; ASAH1 is the seventh gene"),
 "30770930__qpcr_validation": (X, "miRNA list", None, None, None, "serum exosome miRNAs"),
 "42034251__mmc3_65_dems_names": (X, "miRNA list", None, None, None, "blood miRNAs"),
 "42034251__mmc2_core13": (X, "miRNA list", None, None, None, "blood miRNAs (gene selection)"),
 "41691208__table_s4_sdc_degs": (R, None, None, None, None, "blood RNA-seq (DESeq2); Table S4 holds DE rows only"),
 "33782720__table_s3_sz_circrna": (X, "circRNA list", None, None, None, "circRNAs"),
 "33782720__table_s3_sz_clr": (X, "circRNA list", None, None, None, "circular-to-linear ratios of circRNAs (other statistic)"),
 "32015466__table_s3_de_lncrna": (X, "lncRNA list", None, None, None, "lncRNAs"),
 "40263528__st1_discovery_psychotic": (X, "miRNA list", 495, "miRNAs", "manifest universe: '495 expressed miRNAs'", "whole-blood miRNAs"),
 "40263528__st1_validation_psychotic": (X, "miRNA list", 495, "miRNAs", "manifest universe: '495 expressed miRNAs'", "whole-blood miRNAs"),
 "41988836__table_s10_mrna": (R, None, None, None, None, "leukocyte total RNA-seq (edgeR); Table S10 holds DE rows only"),
 "33000265__table2_mirna": (X, "miRNA list", None, None, None, "PBMC miRNAs"),
 "32038049__table_s1": (R, None, None, None, None, "leukocyte total RNA-seq (edgeR); Table S1 holds DE rows only"),
 # batch 16 (2026-09-30): the iPSC-lane round (round 3) and round-4 includes (manifests/batch16_genomewide_2026_09_30.json;
 # primary lists from tools/b16_primary_lists.py; proteomic lists take 'exclude' in main())
 "24776767__table_s5_fc_ba9": (R, None, None, None, None, "BA9 poly(A) RNA-seq; Table S5 holds the > 1.5-fold genes only"),
 "24776767__table_s5_fc_ba24": (R, None, None, None, None, "BA24 poly(A) RNA-seq; Table S5 holds the > 1.5-fold genes only"),
 "31988434__table_s3_circrna": (X, "circRNA list", None, None, None, "circRNAs by host gene (Arraystar circRNA array)"),
 "35301253__appendix1_table_s1": (R, None, None, None, None, "whole-blood RNA-seq (BGIseq-500); Appendix 1 Table S1 holds DE rows only"),
 "35898501__supp_table2_monocyte_degs": (T, None, 79, "genes", "manifest universe: '79 monocyte subset signature genes (Results)'",
                                         "DE genes among 79 monocyte genes chosen from the literature (modality gene selection); the 79 are not listed"),
 "32109548__table_s4_il1b_interaction": (R, None, 17880, "genes", "manifest universe: '17,880 genes examined after pre-filtering (< 1 CPM in more than 3 samples removed) and DESeq2 independent filtering (Methods 2.8)'",
                                         "iPSC-astrocyte RNA-seq; Table S4 holds the nominal interaction genes only (modality other statistic)"),
 "33571628__table_s3_baseline_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "33571628__table_s7_il1b_interaction_mirna": (X, "miRNA list", None, None, None, "mature miRNAs"),
 "34357384__suppl_table6_subtype_astrocytes": (R, None, None, None, None, "iPSC-astrocyte RNA-seq (DESeq2); Suppl. Table 6 holds DE rows only; the tested count is not stated"),
 "34716291__suppl_table6_day0": (R, None, 28564, "genes", "manifest universe: '28,564 genes retained after pre-filtering (<= 3 reads in more than 50% of samples removed) and DESeq2 independent filtering (Supplementary Methods, 'Differential expression analysis')'",
                                 "iPSC astroglial differentiation RNA-seq; Suppl. Table 6 holds DE rows only"),
 "34716291__suppl_table6_day7": (R, None, 28564, "genes", "as 34716291__suppl_table6_day0", "as 34716291__suppl_table6_day0"),
 "34716291__suppl_table6_day14": (R, None, 28564, "genes", "as 34716291__suppl_table6_day0", "as 34716291__suppl_table6_day0"),
 "34716291__suppl_table6_day30": (R, None, 28564, "genes", "as 34716291__suppl_table6_day0", "as 34716291__suppl_table6_day0"),
 "34716291__suppl_table6_day40": (R, None, 28564, "genes", "as 34716291__suppl_table6_day0", "as 34716291__suppl_table6_day0"),
 "39617069__lpc_interaction_mmc4": (R, None, None, None, None, "iPSC-astrocyte RNA-seq; mmc4 holds the interaction genes at FDR < 0.1 only (modality other statistic)"),
 "39372772__ext_data_table2d_sz_vs_ct": (R, None, None, None, None, "iPSC GABAergic neuron poly-A RNA-seq (DESeq2); Extended Data Table 2d holds DE rows only"),
 "41388144__table_s4_aln_neurons": (R, None, None, None, None, "iPSC dopaminergic neuron RNA-seq (DESeq2); Table S4 holds DE rows only"),
 "33246498__supp_table1_up": (R, None, None, None, None, "iPSC NPC RNA-seq; Supplementary Table 1 holds the genes up against both parents only"),
 "33246498__supp_table2_down": (R, None, None, None, None, "iPSC NPC RNA-seq; Supplementary Table 2 holds the genes down against both parents only"),
 "33246498__supp_table3_khsrp_targets": (X, "gene selection", None, None, None, "KHSRP targets among the DE genes (modality gene selection), as 31157371"),
 "40311511__table_s3_mrna": (R, None, None, None, None, "iPSC whole-transcriptome sequencing; Table S3 holds the mRNAs called by all three methods only"),
 "37507766__additional_file3_degs": (R, None, None, None, None, "iPSC cortical neuron total RNA-seq (Cuffdiff); Additional file 3 holds DE rows only"),
 "37507766__table1_gwas_exome": (T, None, 136, "genes", "manifest universe: '106 prioritized protein-coding GWAS genes and 30 exome genes at FDR < 0.05 (Results), tested within the RNA-seq'",
                                 "DEGs among SZ GWAS and exome genes (modality gene selection); those gene sets are not held"),
 "34685556__table_s1_fdr01_base100": (F, None, None, None, None, "Table S1 sheet DESeq2_out has every gene of the DESeq2 output (30,142 Ensembl genes)"),
 "34685556__table_s1_fc2_up_down": (F, {"sheet": "DESeq2_out", "sheets": None, "header_row": 2, "filters": []}, None, None, None,
                                    "the up and down sheets hold the 2-fold genes only; their universe is the DESeq2_out sheet of the same workbook, as 34685556__table_s1_fdr01_base100"),
 "38684795__table2": (R, None, None, None, None, "iPSC cortical differentiation RNA-seq (DESeq2 LRT); Table 2 holds DE rows only"),
 "27801899__table_s4_22q11_genes": (X, "locus-filtered subset", None, None, None, "genes in the 22q11.2 deleted region only (modality gene selection), as 18226501"),
 "27801899__fig_s4c_top15": (G, "GPL570", None, None, None, "manifest universe: Affymetrix Human Genome U133 Plus 2.0 (GPL570); top 15 up and down of 386 DEGs, as 15892139"),
 "40947439__fig_s7b_shared": (X, "platform unknown", None, None, None, "main text not open access; test, threshold and tested genes not known"),
 "40947439__fig_s6b_labelled": (X, "platform unknown", None, None, None, "as 40947439__fig_s7b_shared"),
 "35705634__fig1_rtqpcr": (T, None, 11, "genes", "manifest universe: 'targeted: 11 genes in Figure 1a-b'", "RT-qPCR panel; its 11 genes are named in the manifest but not held as a gene set"),
 "doi_10.64898_2026.08.14.744920__text_ipsc_degs": (R, None, None, None, None, "iPSC poly(A) RNA-seq (DESeq2); the two DEGs are named in the text only"),
 "doi_10.64898_2026.08.14.744920__supp_table4_syngo_nsc": (X, "gene selection", None, None, None, "SynGO-annotated NSC DEGs (modality gene selection)"),
 "doi_10.64898_2026.08.14.744920__supp_table5_nsc_commonmind": (X, "gene selection", None, None, None, "NSC DEGs shared with CommonMind DLPFC DEGs (modality gene selection)"),
 # batch 17 (2026-09-30): the round-5 lead screen includes (manifests/batch17_genomewide_2026_09_30.json; primary lists from
 # tools/b17_primary_lists.py; proteomic lists take 'exclude' in main())
 "29507296__table2b_scz_hc": (T, None, 141, "genes", "manifest universe: 'targeted: 141 KEGG Wnt pathway genes on the HumanHT-12 v4 (Methods)'",
                              "SZ vs HC tested for 141 KEGG Wnt pathway genes only (arrays are genome-wide); the 141 are not listed"),
 "29467462__fig2_nsc_qpcr": (T, None, 8, "genes", "manifest universe: 'targeted: 8 genes (Supplementary Table 3 primers)'",
                             "qPCR panel; its 8 genes are named in the manifest but not held as a gene set"),
 # batch 19 (2026-09-30): lists re-extracted from the files fetched by hand and the round-7 includes
 # (manifests/batch19_genomewide_2026_09_30.json; primary lists from tools/b19_primary_lists.py; proteomic lists take 'exclude'
 # in main()).
 "42461238__suppl_table4_de47": (F, {"file": "raw/42461238/figshare/iwbp_a_2696844_sm2740.xlsx", "reader": "xlsx", "sheet": "Supplementary Table 2", "header_row": 4,
                                    "columns": {"gene_reported": "Gene Symbol", "id_reported": "Transcript Cluster id", "effect": "Fold Sz / Control", "p": "p"},
                                    "filters": []},
                                17304, "RNAs", "manifest universe: 'A signal of quantifiable intensity was obtained for 17,304 RNAs in BA 9' (Results p. 722)",
                                "Supplementary Table 2 of the figshare supplement (probe intensities, fold Sz / Control and p for every RNA with a quantifiable signal; Affymetrix Human Exon 1.0 ST)"),
 "42461238__table2_connectome_changed": (G, "GPL5175", 17304, "RNAs", "manifest universe: 'A signal of quantifiable intensity was obtained for 17,304 RNAs in BA 9' (Results p. 722)",
                                         "Affymetrix Human Exon 1.0 ST; the list is the IPA-connectome subset of the DE RNAs (gene selection), as 29463818"),
 "34004381__suppl_s1_acc": (G, "GPL6104", None, None, None, "Illumina HumanRef-8 v2 (GEO GPL6104), as the same group's STC arrays 22441714 / 24287731"),
 "15018814__tables1-2": (G, "GPL80", None, None, None, "manifest platform GPL80 (Hu6800 / HuFL); 'about 6800 human genes' is approximate, as 11296301"),
 "25665552__table2_mirna_array": (X, "miRNA list", None, None, None, "PBMC miRNAs (Affymetrix miRNA 3.0)"),
 "30770930__table_s5_pooled": (X, "miRNA list", 315, "miRNAs", "manifest universe: the pooled table has 315 rows", "serum exosome miRNAs"),
 "26710095__s3_sz_vs_control_sweden": (X, "gene selection", 76, "probes", "manifest universe: diagnosis tested on the 76 CVLT-associated probes only",
                                       "diagnosis tested only on probes selected by association with CVLT (modality gene selection), as 31157371"),
 "32015466__mrna_table_s4": (R, None, 11920, "mRNAs", "manifest universe: '11,920 unique mRNAs' (Results p. 4517)",
                             "white-blood-cell RNA-seq (edgeR); Supplementary Table 4 holds every tested mRNA but the list needs P < 0.05 in two comparisons, as blood RNA-seq 41691208 is rnaseq_default"),
 "20699480__table_s1_ons": (G, "GPL6104", 10515, "probes", "manifest universe: '10,515 probes detected on the microarrays of the SZ and control ONS cells' (Results)",
                            "Illumina HumanRef-8 v2 (GEO GPL6104)"),
 "38493364__fig1_men_partial": (R, None, None, None, None, "LCL RNA-seq (DESeq2); the list is the genes Fig. 1 names of the 48 men-only DEGs"),
 "38493364__women_hla_dqa2": (R, None, None, None, None, "LCL RNA-seq (DESeq2); one DEG in women"),
 "38493364__all_sex_adjusted": (R, None, None, None, None, "LCL RNA-seq (DESeq2); one DEG in the whole sample"),
 "21451437__table3_selected": (G, "GPL5175", None, None, None, "Affymetrix Human Exon 1.0 ST gene-level summaries; the 10 validated genes (gene selection)"),
 # batch 20 (2026-09-30): the round-8 include Kaurani 2024 (manifests/batch20_genomewide_2026_09_30.json; primary list from
 # tools/b20_primary_lists.py)
 "38528182__ev5_pfc_mirna": (X, "miRNA list", None, None, None, "prefrontal cortex miRNAs (small RNA-seq, DESeq2)"),
 "38528182__ev4_blood_mirna": (X, "miRNA list", None, None, None, "PsyCourse whole-blood miRNAs (small RNA-seq, DESeq2)"),
}
PLATFORM_FILES = {
    "GPL80": "GPL80.geo.annot.gz", "GPL570": "GPL570.geo.annot.gz", "GPL96": "GPL96.geo.annot.gz", "GPL97": "GPL97.geo.annot.gz",
    "GPL91": "GPL91.geo.annot.gz", "GPL8300": "GPL8300.geo.annot.gz", "GPL876": "GPL876.geo.annot.gz", "GPL6244": "GPL6244.geo.annot.gz",
    "GPL11532": "GPL11532.geo.annot.gz", "GPL1352": "GPL1352.geo.annot.gz", "GPL13158": "GPL13158.geo.annot.gz",
    "GPL6947": "GPL6947.geo.annot.gz", "GPL10558": "GPL10558.geo.annot.gz", "GPL6104": "GPL6104.geo.annot.gz",
    "GPL6883": "GPL6883.geo.annot.gz", "GPL2891": "GPL2891.geo.annot.gz", "GPL5175": "GPL5175.geo.txt", "GPL138": "GPL138.geo.txt",
    # no GEO .annot for GPL13667 (HG-U219): ID and Gene Symbol cut from the GEO platform table (header lines give the source)
    "GPL13667": "GPL13667.fromtable.annot.gz",
}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def norm(sym):
    # platform annotations write official symbols, so a previous symbol there is the previous-symbol gene even where
    # it is also another gene's alias (extract.symbol_to_hgnc official=True; 2026-09-30)
    return extract.symbol_to_hgnc(sym, official=True)[0] or extract.normalize_gene({"gene_reported": sym}, {})[0]


ENTREZ = None


def platform_genes(gpl):
    path = PL + PLATFORM_FILES[gpl]
    genes = set()
    if path.endswith(".annot.gz"):
        for syms in extract.platform(path).values():
            genes.update(filter(None, (norm(s) for s in syms)))
    elif gpl == "GPL5175":
        # SOFT table: gene_assignment is 'accession // SYMBOL // name // ... /// ...'
        with open(os.path.join(DATA, path), encoding="utf-8", errors="replace") as fh:
            rows = csv.DictReader((l for l in fh if not l.startswith(("#", "!", "^"))), delimiter="\t")
            for r in rows:
                for part in (r.get("gene_assignment") or "").split("///"):
                    f = [x.strip() for x in part.split("//")]
                    if len(f) > 1 and f[1] and f[1] != "---":
                        s = norm(f[1])
                        if s:
                            genes.add(s)
    elif gpl == "GPL138":
        # SOFT table: GENE is the Entrez Gene ID
        global ENTREZ
        if ENTREZ is None:
            ENTREZ = {r["entrez_id"]: r["symbol"] for r in csv.DictReader(open(extract.HGNC_FILE, encoding="utf-8"), delimiter="\t")
                      if r["status"] == "Approved" and r["entrez_id"]}
        with open(os.path.join(DATA, path), encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith(("#", "!", "^"))), delimiter="\t"):
                s = ENTREZ.get((r.get("GENE") or "").strip())
                if s and norm(s):
                    genes.add(norm(s))
    return genes


def platform_set(expr):
    if "|" in expr:
        parts = [platform_set(p) for p in expr.split("|")]
        return set.intersection(*parts)
    if "+" in expr:
        return set.union(*(platform_set(p) for p in expr.split("+")))
    return platform_genes(expr)


def manifest_specs():
    m = json.load(open(os.path.join(DATA, "manifest.json")))["lists"]
    for f in sorted(glob.glob(os.path.join(DATA, "manifests", "*.json"))):
        m += json.load(open(f))["lists"]
    return {f"{s['pmid']}__{s['list_id']}": s for s in m}


def fullstats_override(name):
    """The table spec override for a fullstats list: SPEC's dict, or another list's ('=name'), else None."""
    src = SPEC.get(name, (None, None))[1]
    if isinstance(src, str) and src.startswith("="):
        return SPEC[src[1:]][1]
    return src if isinstance(src, dict) else None


def fullstats_rows(spec, override):
    """Rows of the list's own table with its threshold filters removed; selection filters
    (equals / nonempty) kept; rows without a numeric p (or fdr, or effect) dropped.
    Also returns the merged spec (for direction handling) and the column that decided a row was kept."""
    s = copy.deepcopy(spec)
    s.update(copy.deepcopy(override or {}))
    if not override or "filters" not in override:
        s["filters"] = [f for f in s.get("filters", []) if "op" not in f]
    cols = s["columns"]
    stat = cols.get("p") or cols.get("fdr") or cols.get("effect")
    if stat:
        s["filters"].append({"col": stat, "op": ">=", "value": -1e300})
    return extract.build(s), s, stat


def fullstats_set(spec, override):
    rows, _, stat = fullstats_rows(spec, override)
    return {r["symbol"] for r in rows if r["symbol"]}, len(rows), stat


def list_genes(name):
    rows = [l for l in open(os.path.join(DATA, "lists", name + ".tsv"), encoding="utf-8") if not l.startswith("#")]
    return {r["symbol"] for r in csv.DictReader(rows, delimiter="\t") if r.get("symbol")}


def write_set(key, genes, header):
    os.makedirs(SETS, exist_ok=True)
    path = os.path.join(SETS, key + ".txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(f"# {h}\n" for h in header))
        fh.write("\n".join(sorted(genes)) + "\n")
    return os.path.relpath(path, DATA)


def main():
    extract.HGNC = extract.load_hgnc()
    today = datetime.date.today().isoformat()
    pl = json.load(open(os.path.join(HERE, "primary_lists.json")))["studies"]
    specs = manifest_specs()
    prim = [(k, v, n) for k, v in pl.items() if v.get("in_corpus", True) and v.get("primary") for n in v["primary"]]
    # BU_ONLY=<list>,<list>: rebuild only these lists' records and keep every other record of the current universes.json
    # (a scoped rerun after a few lists changed; the full run is the default). The RNA-seq default set is read from the
    # file the current universes.json names instead of being rebuilt from its seven source tables.
    only = {n for n in os.environ.get("BU_ONLY", "").split(",") if n}
    if only:
        unknown = only - {n for _, _, n in prim}
        if unknown:
            raise SystemExit(f"BU_ONLY names lists that are not primary: {sorted(unknown)}")
        prim_all = prim
        prim = [x for x in prim if x[2] in only]
    for _, v, n in prim:
        # Gemma's own analyses: the tested single-gene probes of our full results table
        if n not in SPEC and specs[n].get("provenance") == "gemma":
            SPEC[n] = (F, None, None, None, None, f"Gemma full results table {specs[n]['file']}: single-gene probes with a contrast p")
        # proteomic lists (batch 10) are their own modality and never enter the DE overlaps
        if n not in SPEC and specs[n].get("modality") == "proteomics":
            SPEC[n] = (X, "proteomics list", None, None, None, "manifest modality 'proteomics' (protein abundance, not RNA)")
    # lists with no curated SPEC entry (rule A, Paul 2026-09-25, made former alternates primary): keep the record the
    # current universes.json already has for them (primary_rule.py select built it: universe_record / union_record)
    upath = os.path.join(HERE, "universes.json")
    prior = json.load(open(upath))["lists"] if os.path.exists(upath) else {}
    carried = {n for _, _, n in prim if n not in SPEC and n in prior}
    missing = [n for _, _, n in prim if n not in SPEC and n not in prior]
    if missing:
        raise SystemExit(f"no universe assigned for: {missing}")
    hgnc = f"hgnc_complete_set.txt sha256 {sha(extract.HGNC_FILE)}"
    cache, out = {}, {}

    def get_full(name, override):
        key = ("full", name, json.dumps(override, sort_keys=True))
        if key not in cache:
            genes, nrows, stat = fullstats_set(specs[name], override)
            src = override.get("file") if override and override.get("file") else specs[name]["file"]
            path = write_set("fullstats__" + name, genes,
                             [f"fullstats gene set for {name}", f"source: {src} sha256 {sha(os.path.join(DATA, src))}",
                              f"rows: {nrows} (threshold filters removed; rows with a numeric '{stat}' kept)", f"hgnc: {hgnc}",
                              f"script: build_universes.py, run {today}"])
            cache[key] = (genes, path, src)
        return cache[key]

    # the RNA-seq default: union of the named full-statistics tables
    rs = set()
    rs_parts = []
    if only:
        prior_doc = json.load(open(upath))
        rs_path = prior_doc["rnaseq_default"]["gene_set"]
        rs_parts = prior_doc["rnaseq_default"]["sources"]
        rs = {l.strip() for l in open(os.path.join(DATA, rs_path), encoding="utf-8") if l.strip() and not l.startswith("#")}
    else:
        for name in RNASEQ_DEFAULT:
            g, _, src = get_full(name, None)
            rs |= g
            rs_parts.append({"list": name, "label": pl[name.split("__")[0]]["label"], "file": src, "n_genes": len(g)})
        rs_path = write_set("rnaseq_default", rs, ["rnaseq_default: union of the gene sets of these full-statistics tables:"]
                            + [f"  {p['label']}: {p['file']} ({p['n_genes']} genes)" for p in rs_parts]
                            + [f"hgnc: {hgnc}", f"script: build_universes.py, run {today}"])

    for study, v, name in prim:
        if name in carried:
            out[name] = dict(prior[name], carried_from="universes.json before this run (no curated SPEC entry)")
            continue
        typ, src, n_tested, unit, n_src, evidence = SPEC[name]
        rec = {"study": study, "label": v["label"], "provenance": v["provenance"], "tissue": v["tissue"], "type": typ,
               "n_tested": n_tested, "n_tested_unit": unit, "n_tested_source": n_src, "evidence": evidence}
        genes, set_path = None, None
        if typ == G or (typ == T and src and src != "own_table"):
            key = ("platform", src)
            if key not in cache:
                g = platform_set(src)
                files = [PL + PLATFORM_FILES[p] for p in src.replace("|", "+").split("+")]
                p = write_set("platform__" + src.replace("|", "_AND_").replace("+", "_OR_"), g,
                              [f"platform gene set {src} ('|' intersection, '+' union): HGNC genes the probes map to"]
                              + [f"source: {f} sha256 {sha(os.path.join(DATA, f))}" for f in files]
                              + [f"hgnc: {hgnc}", f"script: build_universes.py, run {today}"])
                cache[key] = (g, p, None)
            genes, set_path, _ = cache[key]
            rec["source"] = src
        elif typ == F or (typ == T and src == "own_table"):
            genes, set_path, file_ = get_full(name, fullstats_override(name))
            rec["source"] = file_
        elif typ == R:
            genes, set_path = rs, rs_path
            rec["source"] = "rnaseq_default"
        else:
            rec["source"] = None
            rec["reason"] = src if typ == X else "gene set not known"
        L = list_genes(name)
        rec["list_genes"] = len(L)
        if genes is not None:
            rec["gene_set"] = set_path
            rec["n_membership"] = len(genes)
            rec["list_genes_outside"] = len(L - genes)
            use = n_tested is not None and unit == "genes" and n_tested <= len(genes)
            rec["n_k"] = n_tested if use else len(genes)
            rec["n_k_from"] = "n_tested" if use else "n_membership"
        else:
            rec["gene_set"] = None
            rec["n_membership"] = None
        out[name] = rec
        print(name, typ, rec.get("n_membership"), rec.get("n_k"), rec.get("list_genes"), rec.get("list_genes_outside"), flush=True)

    if only:
        # every other list keeps its current record; lists no longer primary (replaced by the lists rebuilt here) are dropped
        studies = {k for k, _, n in prim_all if n in only}
        now_primary = {n for _, _, n in prim_all}
        merged = {n: r for n, r in prior_doc["lists"].items()
                  if n not in only and not (r.get("study") in studies and n not in now_primary)}
        merged.update(out)
        out = merged
    counts = {}
    for r in out.values():
        counts[r["type"]] = counts.get(r["type"], 0) + 1
    doc = {
        "generated": today, "script": "build_universes.py",
        "note": ("Tested gene universe per primary list (primary_lists.json, in_corpus). n_k is the denominator recurrence.py uses: "
                 "n_tested when stated in genes and not larger than n_membership, else n_membership. Counts stated in probe sets, "
                 "transcripts or mixed units are recorded but not used. 'targeted' lists without a gene set and 'exclude' lists "
                 "are left out of the universe model."),
        "inputs": {"primary_lists.json": sha(os.path.join(HERE, "primary_lists.json")), "hgnc": sha(extract.HGNC_FILE)},
        "rnaseq_default": {"gene_set": rs_path, "n_genes": len(rs), "sources": rs_parts},
        "type_counts_lists": counts,
        "lists": out,
    }
    with open(os.path.join(HERE, "universes.json"), "w") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
    print(counts, "rnaseq_default", len(rs))


if __name__ == "__main__":
    main()
