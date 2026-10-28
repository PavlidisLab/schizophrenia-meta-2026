# Q3/Q4: hemoglobin and blood-contamination markers in the 28 sample-level datasets.
# Run from anywhere:  Rscript analyses/hemoglobin_gene_analyses/02_sample_level_blood.R
# Loads the datasets with fairness_sample_level/scripts/common.R (read only), bottom-tenth expression filter OFF so low-expressed
# erythroid genes stay measurable. Outputs: tsv files next to this script.
args <- commandArgs(FALSE); me <- sub("--file=", "", args[grep("--file=", args)])
HERE <- normalizePath(dirname(me)); REPO <- normalizePath(file.path(HERE, "..", ".."))
setwd(file.path(REPO, "analyses/sample_level_reanalysis_of_datasets"))
source("scripts/common.R")
set.seed(20261011)
HB   <- c("HBB","HBA1","HBA2","HBD","HBG1","HBG2","HBM","HBQ1","HBE1","HBZ")
ERY  <- c("HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA","GYPB","AHSP","CA1","BPGM","SLC25A37","EPB42")
NONHB <- setdiff(ERY, HB)    # erythroid genes that are not hemoglobin chains
LEUK <- c("S100A8","S100A9","S100A12","FCGR3B","CSF3R","CXCR2","FPR1","MNDA","CD177","MMP9")
sets <- fread(file.path(REPO, "analyses/ph_and_confound_analyses/recurrent_sets.tsv"), skip = "set\t")[set == "study_brain"]
R209 <- sets$gene; SIGN <- setNames(ifelse(sets$sz_direction == "up", 1, -1), sets$gene)
cat("209 set:", length(R209), "\n")

zrow <- function(E, g) { g <- intersect(g, rownames(E)); if (!length(g)) return(NULL); t(scale(t(E[g, , drop = FALSE]))) } # genes x samples z
# score = mean of per-gene z over the available genes (>= 3)
score_mean <- function(E, g, min_g = 3) { Z <- zrow(E, g); if (is.null(Z) || nrow(Z) < min_g) return(NULL); colMeans(Z, na.rm = TRUE) }
# PC1 of the genes, oriented to agree with the mean z
score_pc1 <- function(E, g, min_g = 3) {
  Z <- zrow(E, g); if (is.null(Z) || nrow(Z) < min_g) return(NULL)
  Z <- Z[apply(Z, 1, function(v) all(is.finite(v))), , drop = FALSE]; if (nrow(Z) < min_g) return(NULL)
  p <- prcomp(t(Z), center = TRUE, scale. = FALSE); s <- p$x[, 1]
  if (cor(s, colMeans(Z)) < 0) s <- -s
  list(s = s / sd(s), varexp = p$sdev[1]^2 / sum(p$sdev^2), n = nrow(Z))
}
lmfit <- function(y, X) {   # coefficient table for the first non-intercept column
  ok <- complete.cases(cbind(y, X)); y <- y[ok]; X <- X[ok, , drop = FALSE]
  if (length(y) < 8) return(NULL)
  f <- tryCatch(lm(y ~ X), error = function(e) NULL); if (is.null(f)) return(NULL)
  s <- summary(f)$coefficients; if (nrow(s) < 2) return(NULL)
  c(coef = s[2, 1], se = s[2, 2], t = s[2, 3], p = s[2, 4], n = length(y))
}
cohens <- function(sc, d) {   # SZ minus control in SD of control score, with age + sex
  covs <- c("age", "male"); covs <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, TRUE)]
  X <- as.matrix(d[, c("sz", covs), with = FALSE]); r <- lmfit(sc / sd(sc[d$sz == 0]), X); r
}

REGR <- REG$dataset
out_assoc <- list(); out_cc <- list(); out_adj <- list(); out_209 <- list(); out_meas <- list(); out_cor <- list()
for (nm in REGR) {
  D <- as.list(REG[dataset == nm]); L <- load_dataset(D, bottom_filter = FALSE); E <- L$E; cv <- L$cv
  cat(nm, ncol(E), "samples", nrow(E), "genes\n")
  cvd <- as.data.frame(cv)
  gs <- c(HB, NONHB, LEUK)
  out_meas[[nm]] <- data.table(dataset = nm, gene = gs, measured = gs %in% rownames(E))
  # scores
  s_ery  <- score_mean(E, ERY);   s_hb <- score_mean(E, HB[1:3], 2)
  p_non  <- score_pc1(E, NONHB);  p_ery <- score_pc1(E, ERY); s_leu <- score_mean(E, LEUK)
  scores <- list(ery_mean_z = s_ery, hbABB_mean_z = s_hb, nonHb_ery_PC1 = if (!is.null(p_non)) p_non$s, ery_PC1 = if (!is.null(p_ery)) p_ery$s, leuk_mean_z = s_leu)
  scores <- scores[!vapply(scores, is.null, TRUE)]
  if (!is.null(p_non)) out_cor[[nm]] <- data.table(dataset = nm, compartment = L$compartment, n_nonHb_genes = p_non$n, nonHb_PC1_var_explained = p_non$varexp,
      cor_HBB = if ("HBB" %in% rownames(E)) cor(E["HBB", ], p_non$s) else NA, cor_HBA1 = if ("HBA1" %in% rownames(E)) cor(E["HBA1", ], p_non$s) else NA,
      cor_HBA2 = if ("HBA2" %in% rownames(E)) cor(E["HBA2", ], p_non$s) else NA, cor_leuk_score = if (!is.null(s_leu)) cor(s_leu, p_non$s) else NA)
  # --- association of each HB / erythroid gene and each score with covariates (controls only; all samples with dx)
  traits <- list(); for (v in c("pH","PMI","RIN","agonal","age")) traits[[v]] <- cv[[v]]; traits[["male"]] <- cv$male
  items <- c(lapply(intersect(c(HB, NONHB), rownames(E)), function(g) E[g, ]), scores); names(items) <- c(intersect(c(HB, NONHB), rownames(E)), names(scores))
  for (it in names(items)) {
    y <- as.numeric(items[[it]]); names(y) <- NULL
    for (v in names(traits)) {
      x <- traits[[v]]; if (all(is.na(x)) || length(unique(na.omit(x))) < 4 && v != "male") next
      if (v == "agonal" && sum(x != min(x, na.rm = TRUE), na.rm = TRUE) < 5) next
      for (scope in c("controls", "all")) {
        keep <- if (scope == "controls") cv$sz == 0 else rep(TRUE, nrow(cv))
        cov <- c(if (v != "age") "age", if (v != "male" && length(unique(cv$male[keep])) > 1) "male", if (scope == "all") "sz")
        X <- cbind(x = x, as.matrix(cv[, ..cov])); r <- lmfit(y[keep], X[keep, , drop = FALSE])
        if (!is.null(r)) out_assoc[[length(out_assoc) + 1]] <- data.table(dataset = nm, compartment = L$compartment, item = it, trait = v, scope = scope, t(r))
      }
    }
  }
  # --- case vs control differences in scores (SD of controls, age + sex), and with quality covariates
  for (sn in names(scores)) {
    sc <- scores[[sn]]
    r0 <- cohens(sc, cv)
    qc <- intersect(c("pH", "PMI", "RIN"), names(cv)[colSums(!is.na(cv[, c("pH","PMI","RIN")])) > 0.8 * nrow(cv)])
    cv2 <- copy(cv); cv2 <- cv2[complete.cases(cv2[, ..qc])]; sc2 <- sc[match(cv2$gsm, cv$gsm)]
    r1 <- if (length(qc) && nrow(cv2) > 15 && length(unique(cv2$sz)) == 2) { X <- as.matrix(cv2[, c("sz", "age", "male", qc), with = FALSE]); lmfit(sc2 / sd(sc2[cv2$sz == 0]), X) } else NULL
    if (!is.null(r0)) out_cc[[length(out_cc) + 1]] <- data.table(dataset = nm, compartment = L$compartment, cluster = L$cluster, score = sn, adj = "age+sex", n_sz = sum(cv$sz), n_ctl = sum(1 - cv$sz), t(r0))
    if (!is.null(r1)) out_cc[[length(out_cc) + 1]] <- data.table(dataset = nm, compartment = L$compartment, cluster = L$cluster, score = sn, adj = paste0("age+sex+", paste(qc, collapse = "+")), n_sz = sum(cv2$sz), n_ctl = sum(1 - cv2$sz), t(r1))
  }
  # --- SZ effect of HB genes and of the 209 recurrent genes, base vs adjusted for blood score
  adjs <- list(base = NULL)
  if (!is.null(p_non)) adjs$nonHb_PC1 <- matrix(p_non$s, ncol = 1, dimnames = list(cv$gsm, "nonHbPC1"))
  if (!is.null(p_ery)) adjs$ery_PC1 <- matrix(p_ery$s, ncol = 1, dimnames = list(cv$gsm, "eryPC1"))
  qc <- intersect(c("pH", "PMI"), names(cv)[colSums(!is.na(cv[, c("pH","PMI")])) > 0.8 * nrow(cv)])
  fits <- list()
  for (an in names(adjs)) {
    ex <- adjs[[an]]
    fits[[an]] <- fit_sz(E, cv, extra_design = ex)
  }
  if (length(qc)) {   # pH+PMI adjusted on the complete-case samples (reported for comparison)
    cc <- cv[complete.cases(cv[, ..qc])]; if (nrow(cc) > 15 && length(unique(cc$sz)) == 2) {
      cvq <- copy(cc); for (q in qc) cvq[[paste0("q_", q)]] <- cvq[[q]]
      fits$quality_base <- fit_sz(E, cc); 
      ex <- as.matrix(cc[, ..qc]); fits$quality_pHPMI <- fit_sz(E, cc, extra_design = ex)
      if (!is.null(p_non)) { ex2 <- cbind(ex, nonHbPC1 = p_non$s[match(cc$gsm, cv$gsm)]); fits$quality_pHPMI_nonHbPC1 <- fit_sz(E, cc, extra_design = ex2) }
    }
  }
  for (an in names(fits)) {
    f <- fits[[an]]; if (is.null(f)) next
    names(f$t) <- names(f$logfc) <- names(f$p) <- f$gene
    gg <- intersect(HB, f$gene)
    if (length(gg)) out_adj[[length(out_adj) + 1]] <- data.table(dataset = nm, compartment = L$compartment, model = an, gene = gg, logfc = f$logfc[gg], t = f$t[gg], p = f$p[gg], n_sz = f$n_sz, n_ctl = f$n_ctl)
    g9 <- intersect(R209, f$gene); sg <- SIGN[g9]
    out_209[[length(out_209) + 1]] <- data.table(dataset = nm, compartment = L$compartment, model = an, n_genes = length(g9), n_sz = f$n_sz, n_ctl = f$n_ctl,
      sum_signed_logfc = sum(sg * f$logfc[g9]), mean_signed_t = mean(sg * f$t[g9]), n_same_dir_p05 = sum(sg * f$t[g9] > 0 & f$p[g9] < 0.05), n_opp_dir_p05 = sum(sg * f$t[g9] < 0 & f$p[g9] < 0.05))
  }
  # null: adjusting for the PC1 of 9 random genes (expression-matched is not attempted), 200 draws; signed 209 logFC sum ratio
  if (!is.null(fits$base) && !is.null(p_non)) {
    g9 <- intersect(R209, fits$base$gene); sg <- SIGN[g9]
    b <- fits$base; names(b$logfc) <- b$gene; base_sum <- sum(sg * b$logfc[g9])
    pool <- setdiff(rownames(E), c(R209, ERY, HB)); rr <- numeric(200)
    for (i in 1:200) {
      pcr <- score_pc1(E, sample(pool, length(NONHB)), 3)
      f <- fit_sz(E[g9, , drop = FALSE], cv, extra_design = matrix(pcr$s, ncol = 1))
      rr[i] <- sum(sg * f$logfc) / base_sum
    }
    fb <- fits$nonHb_PC1; names(fb$logfc) <- fb$gene
    out_209[[length(out_209) + 1]] <- data.table(dataset = nm, compartment = L$compartment, model = "null_randomPC1_ratio_quantiles", n_genes = length(g9), n_sz = b$n_sz, n_ctl = b$n_ctl,
      sum_signed_logfc = base_sum, mean_signed_t = quantile(rr, 0.5), n_same_dir_p05 = quantile(rr, 0.025), n_opp_dir_p05 = quantile(rr, 0.975))
    out_209[[length(out_209)]]$obs_ratio_nonHbPC1 <- sum(sg * fb$logfc[g9]) / base_sum
  }
}
w <- function(x, fn, hdr_extra = "") {
  inputs <- c("analyses/age_and_sex_control_comparisons/data/datasets.tsv", "analyses/ph_and_confound_analyses/recurrent_sets.tsv")
  sh <- function(p) sub(" .*", "", system(paste("shasum -a 256", shQuote(file.path(REPO, p))), intern = TRUE))
  con <- file(file.path(HERE, fn), "w")
  writeLines(c("# script: analyses/hemoglobin_gene_analyses/02_sample_level_blood.R", "# run: 2026-10-11",
    paste0("# input: ", inputs, " sha256 ", substr(vapply(inputs, sh, ""), 1, 12)), "# datasets: the 28 of controls_age_sex/data/datasets.tsv (Gemma processed expression, log2 scale as delivered; bottom-tenth filter off)", hdr_extra), con)
  write.table(rbindlist(x, fill = TRUE), con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
}
w(out_meas, "samplelevel_gene_measured.tsv")
w(out_cor, "samplelevel_score_validity.tsv", "# Pearson correlation of the non-hemoglobin erythroid PC1 (ALAS2 SLC4A1 GYPA GYPB AHSP CA1 BPGM SLC25A37 EPB42, those measured) with HBB, HBA1, HBA2 and the neutrophil score")
w(out_assoc, "samplelevel_assoc_hb_and_scores.tsv", "# lm(y ~ trait + age + sex [+ dx in scope 'all']); coef per unit of trait; t = Wald t")
w(out_cc, "samplelevel_blood_score_case_control.tsv", "# SZ minus control in SD of the control score (age + sex; second row per score adds the pH/PMI/RIN available)")
w(out_adj, "samplelevel_hb_sz_effect_by_adjustment.tsv", "# limma SZ vs control log2FC of each hemoglobin gene; models: base (age+sex), nonHb_PC1/ery_PC1 added, quality_* on complete cases")
w(out_209, "samplelevel_209_signed_effect_by_adjustment.tsv", "# 209 recurrent genes (study_brain): signed by SZ direction; null row = ratio of adjusted/base signed sum over 200 random-gene PC1 adjustments (median in mean_signed_t, 2.5% in n_same_dir_p05, 97.5% in n_opp_dir_p05)")
