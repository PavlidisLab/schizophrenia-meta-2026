# Shared settings and helpers for libd_brainseq/analysis/ (sourced by every R script here).
# Run the scripts from libd_brainseq/analysis/.
#
# Reuse: ../../fairness_sample_level/scripts/common.R is sourced (from its own directory, which it needs for its relative
# paths) for agree() (agreement metrics: Spearman of t and log2 FC, top-N fold, top-100 sign agreement, BH-list fold, pi1),
# pi0_safe(), TOPN and ent2sym (NCBI gene ID -> HGNC symbol, ../../../data/reference/hgnc_complete_set.txt). Its Gemma loaders are
# not used here.
#
# Pipeline (fixed for every analysis in this folder; prep.R): adults (Age > 17), per region; edgeR filterByExpr with
# group = Dx on the region's adult samples; TMM (calcNormFactors) on those samples; logCPM = edgeR::cpm(log = TRUE,
# prior.count = 3); Ensembl genes mapped to HGNC genes by map_libd_genes() below (rule of 2026-09-30: HGNC's Ensembl
# column first, no folding of distinct genes; before that, through their NCBI gene ID, one gene per symbol). Differential
# expression: limma-trend (lmFit, eBayes(trend = TRUE)) on logCPM.
#
# Base model (BASE_COVS): diagnosis + age + sex + snpPC1 + RIN + mean_mitoRate. snpPC1 stands in for race (see prep.R
# output data/race_snpPC1.tsv for how snpPC1 separates the recorded race categories).
suppressMessages({library(limma); library(data.table); library(parallel)})

AN <- normalizePath(".")
FSL <- normalizePath("../../fairness_sample_level")
GL <- normalizePath("../..")
local({ owd <- setwd(FSL); on.exit(setwd(owd)); source("scripts/common.R", local = FALSE) })
setwd(AN)

NCORES <- as.integer(Sys.getenv("NCORES", "4"))
SEED <- 20260925
BASE_COVS <- c("age", "male", "snpPC1", "RIN", "mito")
# limma-trend in fit_coef (default TRUE, every analysis here); ../../split_halves sets FALSE for Gemma datasets, whose
# fairness_sample_level fit_sz uses plain eBayes
FIT_TREND <- TRUE
REGIONS <- c("DLPFC", "HIPPO")

# Ensembl gene -> HGNC gene (rule of 2026-09-30). Until then genes were named through the RSE's NCBI gene ID, which
# disagrees with HGNC's own Ensembl column for 167 RSE genes (CALM3's Ensembl gene carries 805, CALM2's ID) and, with
# one row kept per symbol, folded distinct genes into one. Rule:
#   1. the unversioned Ensembl gene ID through HGNC's ensembl_gene_id column. An Ensembl ID that HGNC lists for two or
#      more genes names no single gene and is left out.
#   2. only for an Ensembl ID that HGNC does not list at all: the RSE NCBI gene ID, else the RSE Symbol when it is an
#      HGNC approved symbol, and only to an HGNC gene whose own Ensembl ID (HGNC column) is blank or absent from the RSE.
#      An HGNC gene whose Ensembl ID is in the RSE is that RSE gene; another RSE gene carrying its NCBI ID or symbol is
#      a different gene and is left out. An HGNC gene reached this way by two or more RSE genes is left out (ambiguous).
#   3. several RSE rows -> one HGNC gene happens only for one Ensembl ID listed twice (the GENCODE _PAR_Y copy of an X/Y
#      pseudoautosomal gene); their counts are summed. Any other many-to-one stops the script.
# rd: RSE rowData (gencodeID, ensemblID, EntrezID, Symbol); hg: HGNC complete set (hgnc_id, symbol, entrez_id,
# ensembl_gene_id). Returns one row per RSE gene, in rd order: hgnc_id, symbol (HGNC approved), entrez (HGNC), source
# (ensembl / entrez / rse_symbol / NA), note (why unmapped).
map_libd_genes <- function(rd, hg) {
  hg <- as.data.table(hg)
  hg[is.na(ensembl_gene_id), ensembl_gene_id := ""]
  he <- hg[ensembl_gene_id != ""]
  amb_e <- he[, .N, by = ensembl_gene_id][N > 1, ensembl_gene_id]
  he <- he[!ensembl_gene_id %in% amb_e]
  m <- data.table(gencodeID = rd$gencodeID, ensemblID = rd$ensemblID, rse_entrez = rd$EntrezID, rse_symbol = rd$Symbol,
                  hgnc_id = NA_character_, source = NA_character_, note = NA_character_)
  i <- match(m$ensemblID, he$ensembl_gene_id)
  m[!is.na(i), `:=`(hgnc_id = he$hgnc_id[i[!is.na(i)]], source = "ensembl")]
  m[ensemblID %in% amb_e, note := "ensembl_id_on_several_hgnc_genes"]
  # fallback targets: HGNC genes with no Ensembl ID in HGNC, or one the RSE does not have
  open_g <- hg[ensembl_gene_id == "" | !ensembl_gene_id %in% m$ensemblID]
  oz <- open_g[!is.na(entrez_id)]
  fb <- is.na(m$hgnc_id) & !(m$ensemblID %in% hg$ensembl_gene_id)
  iz <- match(m$rse_entrez, oz$entrez_id); iz[is.na(m$rse_entrez)] <- NA
  k <- fb & !is.na(iz)
  m[k, `:=`(hgnc_id = oz$hgnc_id[iz[k]], source = "entrez")]
  is_ <- match(m$rse_symbol, open_g$symbol); is_[is.na(m$rse_symbol) | m$rse_symbol == ""] <- NA
  k <- fb & is.na(m$hgnc_id) & !is.na(is_)
  m[k, `:=`(hgnc_id = open_g$hgnc_id[is_[k]], source = "rse_symbol")]
  # an RSE gene that names, by NCBI ID or symbol only, an HGNC gene the RSE has under another Ensembl ID
  taken <- fb & is.na(m$hgnc_id) & ((!is.na(m$rse_entrez) & m$rse_entrez %in% hg[ensembl_gene_id %in% m$ensemblID, entrez_id]) |
                                     (m$rse_symbol %in% hg[ensembl_gene_id %in% m$ensemblID, symbol]))
  m[taken, note := "names_hgnc_gene_held_by_other_ensembl_id"]
  # several RSE genes reaching one HGNC gene by the fallback: ambiguous
  dz <- m[source %in% c("entrez", "rse_symbol"), .N, by = hgnc_id][N > 1, hgnc_id]
  m[hgnc_id %in% dz & source %in% c("entrez", "rse_symbol"), `:=`(hgnc_id = NA_character_, source = NA_character_,
                                                                  note = "several_ensembl_ids_one_hgnc_gene_by_fallback")]
  m[is.na(hgnc_id) & is.na(note), note := ifelse(ensemblID %in% hg$ensembl_gene_id, "not_mapped", "no_hgnc_gene")]
  # many-to-one left after the rules: only copies of one Ensembl ID
  mo <- m[!is.na(hgnc_id), .(n_ens = uniqueN(ensemblID)), by = hgnc_id]
  stopifnot(all(mo$n_ens == 1))
  m[, symbol := hg$symbol[match(hgnc_id, hg$hgnc_id)]]
  m[, entrez := hg$entrez_id[match(hgnc_id, hg$hgnc_id)]]
  m[]
}

sha256 <- function(f) sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE))

load_region <- function(r) readRDS(sprintf("data/libd_%s.rds", r))

# limma-trend fit of one coefficient; covariates with a single value in the samples at hand are dropped (as common.R
# fit_sz does). E: genes x samples logCPM with colnames = RNum; d: data.table with column gsm (RNum) and the covariates.
fit_coef <- function(E, d, var = "sz", covs = BASE_COVS, extra_design = NULL, genes = NULL, trend = FIT_TREND) {
  if (var != "sz") covs <- c("sz", covs)   # age and sex contrasts keep diagnosis in the model
  covs <- setdiff(covs, var)
  covs <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
  des <- model.matrix(as.formula(paste("~", paste(c(var, covs), collapse = " + "))), data = d)
  if (!is.null(extra_design)) des <- cbind(des, extra_design)
  if (qr(des)$rank < ncol(des) || nrow(des) - ncol(des) < 1) return(NULL)
  Y <- E[if (is.null(genes)) seq_len(nrow(E)) else genes, d$gsm, drop = FALSE]
  f <- eBayes(lmFit(Y, des), trend = trend)
  j <- which(colnames(des) == var)
  list(gene = rownames(Y), t = f$t[, j], logfc = f$coefficients[, j], p = f$p.value[, j],
       se = f$stdev.unscaled[, j] * sqrt(f$s2.post), sigma = f$sigma, df = f$df.total[1], covs = covs,
       n = nrow(des), n_sz = sum(d$sz), n_ctl = sum(1 - d$sz))
}

# numbered TSV writer with '#' header lines
write_tsv <- function(x, path, header = character(0)) {
  # a .gz path is written plain then gzipped whole: otherwise the plain header lines sit in front of fwrite's gzip stream and the file is not valid gzip
  out <- sub("\\.gz$", "", path)
  fh <- file(out, "w")
  writeLines(paste("#", c(sprintf("built %s by libd_brainseq/analysis/%s", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT), header)), fh)
  close(fh)
  fwrite(x, out, sep = "\t", append = TRUE, col.names = TRUE)
  if (out != path) system(paste("gzip -f", shQuote(out)))
}

# restrict a fit_coef result to a gene subset
sub_res <- function(x, genes) {
  if (is.null(x)) return(NULL)
  k <- x$gene %in% genes
  x$gene <- x$gene[k]; x$t <- x$t[k]; x$logfc <- x$logfc[k]; x$p <- x$p[k]; x$se <- x$se[k]; x$sigma <- x$sigma[k]
  x
}

# Sampling units for a contrast. var 'sz': groups SZ / control; 'male': groups male / female; 'age': one pool (2n drawn).
# Returns a list of index vectors (row numbers of cv) per group.
groups_of <- function(cv, var) {
  if (var == "sz") return(list(which(cv$sz == 1), which(cv$sz == 0)))
  if (var == "male") return(list(which(cv$male == 1), which(cv$male == 0)))
  list(seq_len(nrow(cv)))
}
# n per group per half; the 'age' pool draws 2n per half so every contrast uses the same number of samples
draw_halves <- function(G, n, seed, k_per_group = if (length(G) == 1) 2L * n else n) {
  set.seed(seed)
  h1 <- integer(0); h2 <- integer(0)
  for (g in G) {
    s <- g[sample.int(length(g))]
    h1 <- c(h1, s[seq_len(k_per_group)]); h2 <- c(h2, s[k_per_group + seq_len(k_per_group)])
  }
  list(h1, h2)
}
max_n <- function(G) if (length(G) == 1) floor(length(G[[1]]) / 4) else floor(min(lengths(G)) / 2)

# donors (BrNum) with an adult sample in both regions, and their diagnosis / sex
both_region_donors <- function(D) {
  a <- D$DLPFC$cv; b <- D$HIPPO$cv
  br <- intersect(a$BrNum, b$BrNum)
  a[match(br, BrNum), .(BrNum, sz, male)]
}

# Split-half engine shared by a_split_half.R (SZ) and d_controls.R (age, sex).
# For one repeat r and one level n (numeric, or "max" = each unit's own largest n), rows for:
#   split_halves                   two halves of one region, no shared subjects (n per group each)
#   same_subjects_diff_region      DLPFC vs hippocampus on the same n per group subjects (subjects with both regions)
#   separate_subjects_diff_region  DLPFC on one half and hippocampus on the other, subject-separate halves of the subjects with both regions
# gene_sets: named list of gene vectors (NULL = all genes of the fit) on which agree() is computed.
splithalf_rep <- function(D, var, r, n, gene_sets = list(all = NULL), level_id = 0L, kinds = c("split_halves", "same_subjects_diff_region", "separate_subjects_diff_region")) {
  out <- list()
  add <- function(kind, region, nn, a, b) {
    if (is.null(a) || is.null(b)) return()
    for (gs in names(gene_sets)) {
      g <- gene_sets[[gs]]
      aa <- if (is.null(g)) a else sub_res(a, g); bb <- if (is.null(g)) b else sub_res(b, g)
      out[[length(out) + 1]] <<- cbind(data.table(contrast = var, gene_set = gs, kind = kind, region = region, level = as.character(n),
                                                  n_per_group = nn, rep = r), suppressWarnings(agree(aa, bb)))
    }
  }
  seed0 <- SEED + 1000000L * level_id + 1000L * r
  if ("split_halves" %in% kinds) for (k in seq_along(REGIONS)) {
    x <- D[[REGIONS[k]]]; G <- groups_of(x$cv, var)
    nn <- if (n == "max") max_n(G) else as.integer(n)
    if (nn > max_n(G)) next
    h <- draw_halves(G, nn, seed0 + k)
    add("split_halves", REGIONS[k], nn, fit_coef(x$E, x$cv[h[[1]]], var), fit_coef(x$E, x$cv[h[[2]]], var))
  }
  # split halves only (D need not hold DLPFC and HIPPO; ../../split_halves passes one dataset)
  if (!any(c("same_subjects_diff_region", "separate_subjects_diff_region") %in% kinds)) return(rbindlist(out))
  B <- both_region_donors(D)
  GB <- groups_of(B, var)
  dl <- D$DLPFC; hp <- D$HIPPO
  on_region <- function(x, br) x$cv[match(br, x$cv$BrNum)]
  if ("same_subjects_diff_region" %in% kinds) {
    nmax <- if (length(GB) == 1) floor(length(GB[[1]]) / 2) else min(lengths(GB))
    nn <- if (n == "max") nmax else as.integer(n)
    if (nn <= nmax) {
      h <- draw_halves(GB, nn, seed0 + 11)
      br <- B$BrNum[h[[1]]]
      add("same_subjects_diff_region", "DLPFC_vs_HIPPO", nn, fit_coef(dl$E, on_region(dl, br), var), fit_coef(hp$E, on_region(hp, br), var))
    }
  }
  if ("separate_subjects_diff_region" %in% kinds) {
    nn <- if (n == "max") max_n(GB) else as.integer(n)
    if (nn <= max_n(GB)) {
      h <- draw_halves(GB, nn, seed0 + 12)
      add("separate_subjects_diff_region", "DLPFC_vs_HIPPO", nn, fit_coef(dl$E, on_region(dl, B$BrNum[h[[1]]]), var),
          fit_coef(hp$E, on_region(hp, B$BrNum[h[[2]]]), var))
    }
  }
  rbindlist(out)
}
