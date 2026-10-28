# Analysis matrices for the LIBD datasets of split_halves/, with the pipeline of ../libd_brainseq/analysis/scripts/prep.R.
#
# Run from split_halves/:  Rscript scripts/prep_libd.R      (after scripts/fetch_libd.py)
# Inputs (read only): raw/<LIBD gene RSE>.rda (scripts/fetch_libd.py, raw/download_meta.json);
#   ../libd_brainseq/analysis/data/libd_DLPFC.rds, libd_HIPPO.rds (BrainSeq Phase 2, built by that folder's prep.R; read,
#   not rebuilt); ../../data/reference/hgnc_complete_set.txt.
# Pipeline (copied from prep.R): adults (Age > 17); SZ and controls only (caudate bipolar samples left out); edgeR
#   Ensembl genes -> HGNC genes by map_libd_genes() (../libd_brainseq/analysis/scripts/libd_common.R, rule of 2026-09-30:
#   HGNC's Ensembl column, fallbacks only to genes the RSE lacks, no folding of distinct genes; counts summed per HGNC
#   gene, only _PAR_Y copies share one; genes with no HGNC gene kept through filterByExpr and TMM, then left out);
#   filterByExpr(group = diagnosis) on the dataset's adult samples; TMM; logCPM (prior.count 3). Before 2026-09-30:
#   the NCBI gene ID, else the RSE Symbol, one gene per symbol (highest mean logCPM), which folded distinct genes.
#   Differential expression later: limma-trend (libd_common.R fit_coef). Per-gene map: results/libd_gene_map_<dataset>.tsv.gz.
# Covariates per dataset (a covariate with one value in the samples at hand is dropped by fit_coef):
#   Phase 1 DLPFC  age + sex + snpPC1 + RIN + mitoRate       (the Phase 2 base model, ../libd_brainseq libd_common.R BASE_COVS)
#   caudate        age + sex + race (AA vs CAUC, 0/1) + RIN + mitoRate   (the RSE has no genotype PCs; race has two
#                  recorded values here, so the indicator takes one degree of freedom as snpPC1 does)
#   habenula       age + sex + snpPC1 + RIN + mitoRate       (all 68 samples male and CAUC: sex drops out)
#   Phase 2 DLPFC / HIPPO: the Phase 2 base model, from the Phase 2 matrices as built there.
# Outputs: data/<dataset>.rds = list(E, cv, genes, covs, trend = TRUE, family = "libd", name); results/prep_libd.tsv
SCRIPT <- "scripts/prep_libd.R"
SH <- normalizePath(".")
local({ owd <- setwd("../libd_brainseq/analysis"); on.exit(setwd(owd)); source("scripts/libd_common.R", local = FALSE) })
setwd(SH)
suppressMessages({library(SummarizedExperiment); library(edgeR)})
dir.create("data", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)

num1 <- function(x) if (is.list(x)) vapply(x, function(v) mean(as.numeric(v)), 0) else as.numeric(x)
HGNC <- "../../data/reference/hgnc_complete_set.txt"
hgf <- fread(HGNC, select = c("hgnc_id", "symbol", "entrez_id", "ensembl_gene_id"), quote = "")
# rowData in map_libd_genes' shape: gencodeID (the RSE row name), unversioned ensemblID, EntrezID, Symbol (Phase 1's
# hg19 RSE names rows by unversioned Ensembl ID and has no gencodeID / ensemblID columns)
rse_genes <- function(rse) {
  rd <- as.data.frame(rowData(rse))
  sym_col <- intersect(c("Symbol", "symbol"), names(rd))[1]
  data.table(gencodeID = rownames(rse), ensemblID = sub("\\..*$", "", rownames(rse)), EntrezID = as.character(rd$EntrezID),
             Symbol = as.character(rd[[sym_col]]))
}

build <- function(name, file, covs, dxcol = "Dx", agecol = "Age", mito = "mitoRate", race_indicator = FALSE) {
  e <- new.env(); nm <- load(file.path("raw", file), envir = e); rse <- e[[nm[1]]]
  cd <- as.data.frame(colData(rse))
  cv <- data.table(gsm = as.character(cd$RNum), BrNum = as.character(cd$BrNum), dx = as.character(cd[[dxcol]]), age = as.numeric(cd[[agecol]]),
                   male = as.integer(cd$Sex == "M"), race = as.character(cd$Race), RIN = num1(cd$RIN), mito = num1(cd[[mito]]),
                   snpPC1 = if ("snpPC1" %in% names(cd)) as.numeric(cd$snpPC1) else NA_real_)
  cv[, race_aa := as.integer(race == "AA")]
  keep_s <- which(cv$dx %in% c("Control", "Schizo") & cv$age > 17)
  cv <- cv[keep_s]; cv[, sz := as.integer(dx == "Schizo")]
  stopifnot(!anyDuplicated(cv$gsm), !anyDuplicated(cv$BrNum), !anyNA(cv[, c("age", "male", "RIN", "mito"), with = FALSE]))
  gm <- map_libd_genes(rse_genes(rse), hgf)
  gm[, chr := as.character(seqnames(rowRanges(rse)))]
  write_tsv(gm[, .(gencodeID, ensemblID, chr, rse_symbol, rse_entrez, hgnc_id, symbol, entrez, source, note)],
            sprintf("results/libd_gene_map_%s.tsv.gz", name),
            c(sprintf("every gene of raw/%s (sha256 %s): HGNC gene by map_libd_genes() (../libd_brainseq/analysis/scripts/libd_common.R); HGNC %s sha256 %s",
                      file, sha256(file.path("raw", file)), HGNC, sha256(HGNC))))
  # counts summed per HGNC gene (only _PAR_Y copies share one); Ensembl genes with no HGNC gene stay as their own rows
  # through filterByExpr and TMM and are left out of E (as prep.R)
  key <- ifelse(is.na(gm$hgnc_id), gm$gencodeID, gm$hgnc_id)
  counts <- rowsum(assay(rse, "counts")[, keep_s], key, reorder = FALSE)
  grow <- gm[, .(ensembl = paste(gencodeID, collapse = ";"), n_ensembl = .N, hgnc_id = hgnc_id[1], symbol = symbol[1], entrez = entrez[1],
                 rse_symbol = rse_symbol[1], chr = chr[1], symbol_source = source[1]), by = .(key = key)]
  stopifnot(identical(grow$key, rownames(counts)))
  y <- DGEList(counts = counts, group = cv$dx)
  keep <- filterByExpr(y, group = cv$dx)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lc <- cpm(y, log = TRUE, prior.count = 3)
  gi <- grow[match(rownames(lc), key)][, mean_logcpm := rowMeans(lc)]
  gi <- gi[!is.na(hgnc_id)][order(symbol)]
  stopifnot(!anyDuplicated(gi$symbol), !anyDuplicated(gi$hgnc_id))
  E <- lc[gi$key, , drop = FALSE]; rownames(E) <- gi$symbol; colnames(E) <- cv$gsm
  gi[, key := NULL]
  gi[, autosomal := chr %in% c(paste0("chr", 1:22), as.character(1:22))]
  x <- list(E = E, cv = cv, genes = gi, covs = covs, trend = TRUE, family = "libd", name = name)
  saveRDS(x, sprintf("data/%s.rds", name))
  cat(name, "samples", ncol(E), "SZ", sum(cv$sz), "control", sum(1 - cv$sz), "genes", nrow(E), "\n")
  data.table(dataset = name, source = file, n_sz = sum(cv$sz), n_ctl = sum(1 - cv$sz), n_male = sum(cv$male), n_female = sum(1 - cv$male),
             genes_total = nrow(rse), genes_filterByExpr = sum(keep), genes_unique_symbol = nrow(E), genes_autosomal = sum(gi$autosomal),
             covariates = paste(covs, collapse = " + "), trend = TRUE)
}

rows <- list(
  build("LIBD_Phase1_DLPFC", "rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda", c("age", "male", "snpPC1", "RIN", "mito")),
  build("LIBD_Caudate", "caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda", c("age", "male", "race_aa", "RIN", "mito")),
  build("LIBD_Habenula", "rse_gene_Habenula_Pilot.rda", c("age", "male", "snpPC1", "RIN", "mito"), dxcol = "PrimaryDx", agecol = "AgeDeath"))

# BrainSeq Phase 2, as built by ../libd_brainseq/analysis/scripts/prep.R
for (r in c("DLPFC", "HIPPO")) {
  f <- sprintf("../libd_brainseq/analysis/data/libd_%s.rds", r)
  x <- readRDS(f)
  x$cv[, race_aa := as.integer(race == "AA")]
  x$covs <- BASE_COVS; x$trend <- TRUE; x$family <- "libd"; x$name <- paste0("LIBD_Phase2_", r)
  saveRDS(x, sprintf("data/%s.rds", x$name))
  rows[[length(rows) + 1]] <- data.table(dataset = x$name, source = f, n_sz = sum(x$cv$sz), n_ctl = sum(1 - x$cv$sz), n_male = sum(x$cv$male),
                                          n_female = sum(1 - x$cv$male), genes_total = NA, genes_filterByExpr = NA, genes_unique_symbol = nrow(x$E),
                                          genes_autosomal = sum(x$genes$autosomal), covariates = paste(BASE_COVS, collapse = " + "), trend = TRUE)
}
fwrite(rbindlist(rows), "results/prep_libd.tsv", sep = "\t")
cat("done\n")
