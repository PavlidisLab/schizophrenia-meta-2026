# Shared loaders, donor groups and the DE pipeline for the PsychAD count pseudobulk analyses (sourced by de.R,
# split_half.R, libd_compare.R, st4_compare.R). Run the scripts from analyses/psychad_pseudobulk/.
#
# Reuse: ../fairness_sample_level/scripts/common.R is sourced (from its own directory, which it needs for its relative
# paths) for agree() (Spearman of t and of log2 FC, top-N overlap fold, top-100 sign agreement, BH-list fold, pi1),
# pi0_safe() and TOPN. Its Gemma loaders are not used.
#
# Pseudobulks: pseudobulk/<bank file>/counts_<level>.npz (scipy CSC, genes x pseudobulks, int32 UMI sums; built on remote
# by scripts/build_pseudobulk.py) with cols_<level>.tsv and genes.tsv.
#
# Pipeline (every fit): pseudobulks with >= MIN_NUC nuclei; edgeR filterByExpr(group = diagnosis); TMM; voom; limma
# lmFit + eBayes. Model: ~ sz + male + age + PMI + log(nuclei). Ensembl ids are reported as the h5ad var gene_name; when
# several Ensembl ids share a name the one with the highest voom average expression is kept.
suppressMessages({library(edgeR); library(limma); library(data.table); library(parallel); library(Matrix)})

PC <- normalizePath(".")
FSL <- normalizePath("../fairness_sample_level")
local({ owd <- setwd(FSL); on.exit(setwd(owd)); source("scripts/common.R", local = FALSE) })
setwd(PC)

NCORES <- as.integer(Sys.getenv("NCORES", "4"))
SEED <- 20260925
MIN_NUC <- 20      # nuclei per donor x cell-type pseudobulk
MIN_DONORS <- 20   # donors per group (SZ, control) for a cell type to be analysed
COVS <- c("male", "age", "pmi", "lognuc")
PB <- "pseudobulk"
BANK_FILE <- c(HBCC = "HBCC_Cohort", MSSM = "MSSM_Cohort")

sha256 <- function(f) sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE))

write_tsv <- function(x, path, header = character(0)) {
  gz <- grepl("\\.gz$", path)
  hdr <- paste("#", c(sprintf("built %s by psychad_cellxgene/%s", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT), header))
  if (gz) {
    tmp <- sub("\\.gz$", "", path)
    writeLines(hdr, tmp); fwrite(x, tmp, sep = "\t", append = TRUE, col.names = TRUE)
    system(paste("gzip -f", shQuote(tmp)))
  } else {
    writeLines(hdr, path); fwrite(x, path, sep = "\t", append = TRUE, col.names = TRUE)
  }
}

# ---- .npz (scipy sparse) reader --------------------------------------------------------------------------------
read_npy <- function(npz, name) {
  con <- pipe(paste("unzip -p", shQuote(npz), paste0(name, ".npy")), "rb")
  on.exit(close(con))
  readBin(con, "raw", 6)
  ver <- readBin(con, "integer", 2, size = 1)
  hl <- if (ver[1] == 1) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little") else readBin(con, "integer", 1, size = 4, endian = "little")
  hdr <- rawToChar(readBin(con, "raw", hl))
  descr <- sub(".*'descr': *'([^']+)'.*", "\\1", hdr)
  shp <- trimws(strsplit(sub(".*'shape': *\\(([^)]*)\\).*", "\\1", hdr), ",")[[1]])
  n <- prod(as.numeric(shp[nzchar(shp)]))
  if (descr %in% c("<i4", "<u4")) return(readBin(con, "integer", n, size = 4, endian = "little"))
  if (descr == "<i8") { v <- readBin(con, "integer", 2 * n, size = 4, endian = "little"); return(v[c(TRUE, FALSE)]) }
  if (descr == "<f8") return(readBin(con, "double", n, size = 8, endian = "little"))
  stop("unhandled npy dtype ", descr)
}
read_npz_csc <- function(npz) {
  shape <- read_npy(npz, "shape")
  i <- read_npy(npz, "indices"); p <- read_npy(npz, "indptr"); x <- read_npy(npz, "data")
  new("dgCMatrix", i = i, p = p, x = as.double(x), Dim = as.integer(shape))
}

# ---- donors -----------------------------------------------------------------------------------------------------
# donors.tsv (Gemma sample characteristics joined to Supp. Table 1; build_donors.py). Age and PMI from Supp. Table 1
# (paper_age, '89+' read as 89; paper_pmi), sex from paper_sex. dementia: the CELLxGENE disease field names dementia or
# Alzheimer disease.
load_donors <- function() {
  d <- fread("donors.tsv", colClasses = "character")
  d[, age := as.numeric(sub("\\+", "", paper_age))]
  d[, pmi := as.numeric(paper_pmi)]
  d[, male := as.integer(paper_sex == "M")]
  d[, dementia := grepl("dementia|alzheimer", disease, ignore.case = TRUE)]
  d[, bank := c(H = "HBCC", M = "MSSM", R = "RADC")[source_bank]]
  d
}

# Donor sets. Each returns data.table(donor_id, sz) for one bank.
#   HBCC_main        SZ (dx_gemma SZ) vs controls whose CELLxGENE disease is 'normal', 18 or older
#   HBCC_ctrlpaper   SZ vs the paper's controls (Supp. Table 1 CTRL: normal plus somatic-only), 18 or older
#   MSSM_main        SZ vs 'normal' controls
#   MSSM_nodementia  SZ without a dementia / Alzheimer term vs 'normal' controls
#   MSSM_ctrlpaper   SZ vs the paper's controls (adds the dementia-only donors, drops 12 'normal' donors the paper calls OTHER)
SETS <- c("HBCC_main", "HBCC_ctrlpaper", "MSSM_main", "MSSM_nodementia", "MSSM_ctrlpaper")
donor_set <- function(D, set) {
  b <- sub("_.*", "", set)
  x <- D[D$bank == b]
  sz <- x[dx_gemma == "SZ"]
  if (set == "MSSM_nodementia") sz <- sz[dementia == FALSE]
  ctl <- if (grepl("ctrlpaper", set)) x[paper_dx == "CTRL"] else x[dx_gemma == "CTRL"]
  if (b == "HBCC") ctl <- ctl[age >= 18]
  rbind(sz[, .(donor_id, sz = 1L)], ctl[, .(donor_id, sz = 0L)])
}

# ---- pseudobulk loading ------------------------------------------------------------------------------------------
load_level <- function(bank, level) {
  dir <- file.path(PB, BANK_FILE[[bank]])
  M <- read_npz_csc(file.path(dir, sprintf("counts_%s.npz", level)))
  cols <- fread(file.path(dir, sprintf("cols_%s.tsv", level)))
  genes <- fread(file.path(dir, "genes.tsv"))
  dimnames(M) <- list(genes$ensembl, cols$pb_id)
  list(M = M, cols = cols, genes = genes)
}

# samples (pseudobulks) of one cell type for one donor set, with covariates; NULL when a group has < MIN_DONORS
cell_samples <- function(L, D, set, ct, min_nuc = MIN_NUC) {
  s <- donor_set(D, set)
  cc <- L$cols[cell_type == ct & n_nuclei >= min_nuc]
  d <- merge(cc, s, by = "donor_id")
  d <- merge(d, D[, .(donor_id, male, age, pmi, dementia)], by = "donor_id")
  d <- d[!is.na(age) & !is.na(pmi) & !is.na(male)]
  d[, lognuc := log(n_nuclei)]
  setorder(d, -sz, donor_id)
  d
}
eligible <- function(d) !is.null(d) && sum(d$sz == 1) >= MIN_DONORS && sum(d$sz == 0) >= MIN_DONORS

# ---- DE ----------------------------------------------------------------------------------------------------------
fit_pb <- function(M, genes, d, covs = COVS, full = FALSE) {
  covs <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
  Y <- as.matrix(M[, d$pb_id, drop = FALSE])
  Y <- Y[rowSums(Y) > 0, , drop = FALSE]   # genes with no counts in these samples (would fail filterByExpr anyway)
  y <- DGEList(Y)
  keep <- filterByExpr(y, group = d$sz)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  des <- model.matrix(as.formula(paste("~ sz +", paste(covs, collapse = " + "))), data = d)
  if (qr(des)$rank < ncol(des)) return(NULL)
  v <- voom(y, des)
  f <- eBayes(lmFit(v, des))
  j <- which(colnames(des) == "sz")
  gname <- genes$gene_name[match(rownames(y), genes$ensembl)]
  o <- order(gname, -f$Amean)
  o <- o[!duplicated(gname[o])]
  r <- list(gene = gname[o], t = f$t[o, j], logfc = f$coefficients[o, j], p = f$p.value[o, j],
            n_sz = sum(d$sz), n_ctl = sum(1 - d$sz), n_genes = length(o), df = f$df.total[1], covs = covs)
  if (full) {
    r$ensembl <- rownames(y)[o]; r$se <- (f$stdev.unscaled[, j] * sqrt(f$s2.post))[o]; r$ave <- f$Amean[o]
    r$median_nuclei <- median(d$n_nuclei); r$median_lib <- median(d$lib_size)
  }
  r
}
res_table <- function(r) data.table(gene = r$gene, ensembl = r$ensembl, logfc = signif(r$logfc, 5), se = signif(r$se, 5),
                                    t = signif(r$t, 5), p = signif(r$p, 5), q = signif(p.adjust(r$p, "BH"), 5), ave_logcpm = signif(r$ave, 5))
