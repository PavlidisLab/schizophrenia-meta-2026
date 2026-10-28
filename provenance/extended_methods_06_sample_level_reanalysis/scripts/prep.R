# Build the analysis matrices and sample tables from the downloaded BrainSeq Phase 2 objects.
#
# Run from libd_brainseq/analysis/:  Rscript scripts/prep.R
# Inputs (read only): ../raw/rse_gene_unfiltered.Rdata (rse_gene), ../raw/methprop_pd.Rdata (pd; mean_mitoRate),
#   ../raw/phenotype_annotated_szControlEqtl_DLPFC_withDBGAP.csv (Phase 1 DLPFC phenotypes; donor-level PMI, manner of
#   death, toxicology and smoking joined by BrNum); ../hippo_kit_exclusion.tsv (43 hippocampus samples left out).
# Outputs:
#   data/libd_<region>.rds        list(E = logCPM genes x samples, cv = sample table, genes = gene table)  [gitignored]
#   results/phase2_donors.tsv     every Phase 2 sample (all ages): RNum, BrNum, Region, Dx, Age, Sex
#   results/samples.tsv           one row per adult sample used: RNum, BrNum, region, dx, covariates, Phase 1 fields
#   results/prep_genes.tsv        per region: genes at each filtering step
#   results/libd_gene_map.tsv.gz  every RSE Ensembl gene: the HGNC gene it is mapped to (or why none), and by which route
#   results/race_snpPC1.tsv       snpPC1 by recorded race (per region)
SCRIPT <- "scripts/prep.R"

# Ensembl gene -> HGNC gene: map_libd_genes() in scripts/libd_common.R (rule of 2026-09-30; moved there 2026-09-30 so that
# ../../split_halves/scripts/prep_libd.R maps its LIBD datasets with the same function).
source("scripts/libd_common.R")
suppressMessages({library(SummarizedExperiment); library(edgeR)})

load("../raw/rse_gene_unfiltered.Rdata")
load("../raw/methprop_pd.Rdata")
stopifnot(identical(colnames(rse_gene), rownames(pd)))
cd <- as.data.table(as.data.frame(pd)[, c("RNum", "BrNum", "Region", "mean_RIN", "Age", "Sex", "Race", "Dx", paste0("snpPC", 1:10),
                                         "mean_mitoRate", "rRNA_rate", "totalAssignedGene")])
p1 <- fread("../raw/phenotype_annotated_szControlEqtl_DLPFC_withDBGAP.csv", na.strings = c("", "NA", "NULL"))
p1 <- p1[, .(BrNum, p1_Dx = Dx, p1_Age = Age, PMI, Manner_Of_Death, Antipsychotics, Antidepressants, Cotinine, Nicotine, SmokingEither)]
stopifnot(!any(duplicated(p1$BrNum)))
cv <- merge(cd, p1, by = "BrNum", all.x = TRUE, sort = FALSE)
cv <- cv[match(cd$RNum, cv$RNum)]
cat("Phase 1 join: samples", nrow(cv), "with Phase 1 row", sum(!is.na(cv$p1_Dx)), "dx disagree", sum(!is.na(cv$p1_Dx) & cv$p1_Dx != cv$Dx),
    "max |age diff|", max(abs(cv$Age - cv$p1_Age), na.rm = TRUE), "\n")
cv[, `:=`(gsm = RNum, region = Region, dx = Dx, sz = as.integer(Dx == "Schizo"), age = Age, male = as.integer(Sex == "M"),
          race = Race, mito = mean_mitoRate, RIN = mean_RIN)]
stopifnot(is.numeric(cv$RIN), is.numeric(cv$mito), !anyNA(cv[, .(age, RIN, mito, snpPC1)]))
# donor-level Phase 1 fields as model covariates (NA stays NA; complete cases are taken where they are used)
cv[, manner := factor(Manner_Of_Death)]
cv[, antipsych := ifelse(Antipsychotics %in% c("Positive", "Negative"), as.integer(Antipsychotics == "Positive"), NA_integer_)]
cv[, antidep := ifelse(Antidepressants %in% c("Positive", "Negative"), as.integer(Antidepressants == "Positive"), NA_integer_)]
cv[, smoking := ifelse(SmokingEither %in% c("Yes", "No"), as.integer(SmokingEither == "Yes"), NA_integer_)]

# 43 adult hippocampus samples LIBD excluded for a different sequencing kit (Ribo-Zero Gold instead of HMR; Collado-Torres
# et al. 2019 p. 207 and p. e2), identified as the adult HIPPO samples absent from LIBD's public hippocampus qSVA design
# (../hippo_kit_exclusion.tsv, built by scripts/derive_kit_exclusion.R; provenance ../hippo_kit_exclusion_meta.json).
# Dropped from every analysis here (adversarial review round 1, F4); all 43 are controls.
kit <- fread("../hippo_kit_exclusion.tsv", skip = "RNum\t")
stopifnot(nrow(kit) == 43, all(kit$RNum %in% cv[Region == "HIPPO", RNum]))
cat("hippocampus kit exclusion:", nrow(kit), "samples (", sha256("../hippo_kit_exclusion.tsv"), ")\n")
ad <- cv$Age > 17 & !(cv$RNum %in% kit$RNum)
write_tsv(cv[, .(RNum, BrNum, Region, Dx, Age, Sex)], "results/phase2_donors.tsv", "every BrainSeq Phase 2 sample (all ages): RNum, BrNum, region, diagnosis, age, sex (../raw/methprop_pd.Rdata)")
samples <- cv[ad]
write_tsv(samples[, .(RNum, BrNum, region, dx, age, sex = Sex, race, RIN, mito, snpPC1, snpPC2, PMI, Manner_Of_Death, Antipsychotics,
                      Antidepressants, SmokingEither, in_phase1 = !is.na(p1_Dx))], "results/samples.tsv",
          c("adult (Age > 17) BrainSeq Phase 2 samples, less the 43 hippocampus samples LIBD excluded for a different sequencing kit (../hippo_kit_exclusion.tsv); Phase 1 fields joined by BrNum from ../raw/phenotype_annotated_szControlEqtl_DLPFC_withDBGAP.csv",
            sprintf("rse sha256 %s; pd sha256 %s; phase1 csv sha256 %s", sha256("../raw/rse_gene_unfiltered.Rdata"),
                    sha256("../raw/methprop_pd.Rdata"), sha256("../raw/phenotype_annotated_szControlEqtl_DLPFC_withDBGAP.csv"))))

rd <- as.data.table(as.data.frame(rowData(rse_gene)))
chr <- as.character(seqnames(rowRanges(rse_gene)))
HGNC <- "../../../data/reference/hgnc_complete_set.txt"
hgf <- fread(HGNC, select = c("hgnc_id", "symbol", "entrez_id", "ensembl_gene_id"), quote = "")
gm <- map_libd_genes(rd, hgf)
gm[, `:=`(gene_type = rd$gene_type, chr = chr)]
write_tsv(gm[, .(gencodeID, ensemblID, chr, gene_type, rse_symbol, rse_entrez, hgnc_id, symbol, entrez, source, note)],
          "results/libd_gene_map.tsv.gz",
          c(sprintf("every gene of ../raw/rse_gene_unfiltered.Rdata (sha256 %s): HGNC gene by map_libd_genes() (rule in scripts/libd_common.R); HGNC %s sha256 %s",
                    sha256("../raw/rse_gene_unfiltered.Rdata"), HGNC, sha256(HGNC)),
            "source: ensembl = HGNC ensembl_gene_id column; entrez / rse_symbol = fallback for Ensembl IDs HGNC does not list; note: why a gene has no HGNC gene"))
cat("gene map:", paste(names(table(gm$source, useNA = "ifany")), table(gm$source, useNA = "ifany"), collapse = ", "), "\n")
# counts summed per HGNC gene (only _PAR_Y copies share one); Ensembl genes with no HGNC gene stay as their own rows
# through filterByExpr and TMM (as before 2026-09-30) and are left out of E
key <- ifelse(is.na(gm$hgnc_id), gm$gencodeID, gm$hgnc_id)
counts <- rowsum(assay(rse_gene, "counts"), key, reorder = FALSE)
grow <- gm[, .(ensembl = paste(gencodeID, collapse = ";"), n_ensembl = .N, hgnc_id = hgnc_id[1], symbol = symbol[1], entrez = entrez[1],
               chr = chr[1], gene_type = gene_type[1], symbol_source = source[1]), by = .(key = key)]
stopifnot(identical(grow$key, rownames(counts)))
gstep <- list(); race_rows <- list()
for (r in REGIONS) {
  s <- which(cv$Region == r & ad)
  d <- cv[s]
  y <- DGEList(counts = counts[, s], group = d$dx)
  keep <- filterByExpr(y, group = d$dx)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lc <- cpm(y, log = TRUE, prior.count = 3)
  gi <- grow[match(rownames(lc), key)][, mean_logcpm := rowMeans(lc)]
  n_sym <- sum(!is.na(gi$hgnc_id))
  gi <- gi[!is.na(hgnc_id)][order(symbol)]
  stopifnot(!anyDuplicated(gi$symbol), !anyDuplicated(gi$hgnc_id))
  E <- lc[gi$key, , drop = FALSE]; rownames(E) <- gi$symbol
  gi[, key := NULL]
  colnames(E) <- d$gsm
  gi[, autosomal := chr %in% paste0("chr", 1:22)]
  saveRDS(list(E = E, cv = d, genes = gi, norm_factors = setNames(y$samples$norm.factors, d$gsm)), sprintf("data/libd_%s.rds", r))
  gstep[[r]] <- data.table(region = r, n_samples = length(s), n_sz = sum(d$sz), n_ctl = sum(1 - d$sz), genes_total = nrow(rse_gene),
                           genes_filterByExpr = sum(keep), genes_with_hgnc_symbol = n_sym, genes_symbol_from_ensembl = sum(gi$symbol_source == "ensembl"),
                           genes_symbol_from_entrez = sum(gi$symbol_source == "entrez"),
                           genes_symbol_from_rse_symbol = sum(gi$symbol_source == "rse_symbol"), genes_unique_symbol = nrow(E),
                           genes_autosomal = sum(gi$autosomal), genes_chrX = sum(gi$chr == "chrX"), genes_chrY = sum(gi$chr == "chrY"))
  race_rows[[r]] <- d[, .(n = .N, snpPC1_median = median(snpPC1), snpPC1_min = min(snpPC1), snpPC1_max = max(snpPC1)), by = .(region, race)]
  # separation of AA from CAUC by snpPC1 (AUC)
  aa <- d[race == "AA", snpPC1]; ca <- d[race == "CAUC", snpPC1]
  race_rows[[r]][, auc_AA_vs_CAUC := mean(outer(aa, ca, ">")) + 0.5 * mean(outer(aa, ca, "=="))]
  cat(r, "samples", length(s), "genes", nrow(E), "\n")
}
write_tsv(rbindlist(gstep), "results/prep_genes.tsv",
          c("Ensembl gene -> HGNC gene through HGNC's ensembl_gene_id column (../../../data/reference/hgnc_complete_set.txt); for Ensembl IDs HGNC does not list, the RSE NCBI gene ID, else the RSE Symbol, only to an HGNC gene the RSE has no Ensembl ID for (rule and exclusions: map_libd_genes() in scripts/libd_common.R, per gene: results/libd_gene_map.tsv.gz); counts summed per HGNC gene (only _PAR_Y copies share one); no two genes folded",
            "then edgeR filterByExpr(group = dx) on the region's adult samples (genes_filterByExpr counts HGNC genes plus unmapped Ensembl genes); TMM; logCPM prior.count 3; unmapped genes left out after normalization; genes_with_hgnc_symbol = genes_unique_symbol = rows of E; genes_symbol_from_* by mapping route"))
write_tsv(rbindlist(race_rows), "results/race_snpPC1.tsv", "snpPC1 by recorded race; auc_AA_vs_CAUC = P(snpPC1 of an AA sample > that of a CAUC sample)")

# donors with both regions (adults)
both <- intersect(cv[ad & Region == "DLPFC", BrNum], cv[ad & Region == "HIPPO", BrNum])
bd <- unique(cv[BrNum %in% both, .(BrNum, dx)])
cat("adult donors with both regions:", length(both), "SZ", sum(bd$dx == "Schizo"), "control", sum(bd$dx == "Control"), "\n")
print(samples[, .(.N, female = sum(male == 0), with_PMI = sum(!is.na(PMI)), l4_complete = sum(!is.na(PMI) & !is.na(manner) & !is.na(antipsych) & !is.na(antidep) & !is.na(smoking))), by = .(region, dx)])
cat("done\n")
