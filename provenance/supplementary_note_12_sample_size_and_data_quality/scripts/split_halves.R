# Split-half benchmark and positive controls for every dataset in data/, and cross-dataset pairs at matched n.
# The engine is ../libd_brainseq/analysis/scripts/libd_common.R (splithalf_rep, draw_halves, groups_of, max_n, fit_coef,
# sub_res) with agree() of ../fairness_sample_level/scripts/common.R; the designs are those of ../libd_brainseq/analysis/
# scripts/a_split_half.R (SZ) and d_controls.R (SZ, age, sex).
#
# Run from split_halves/:
#   Rscript scripts/split_halves.R within <dataset> [<dataset> ...]
#   Rscript scripts/split_halves.R cross  <pair_id> [<pair_id> ...]      (pairs in PAIRS below)
# Env: A_REPS (default 100), D_REPS (default 50), NCORES (default 4), PARTS (default "sz,controls").
#
# within, part 'sz' (a_split_half.R design): SZ and control subjects of one dataset shuffled and cut into two disjoint
#   halves of n per group; levels 5, 8, 10, 12, 15, 20, 25, 30, 50, 75 per group where the dataset has 2n per group, and
#   'max' (floor(min(n_SZ, n_control) / 2)); A_REPS repeats. Output results/within_sz/<dataset>.tsv.gz.
# within, part 'controls' (d_controls.R design): contrasts sz (n SZ + n controls per half), male (n males + n females per
#   half; model ~ male + sz + other covariates), age (2n subjects drawn at random per half, age continuous; ~ age + sz +
#   other covariates); levels 10, 20, 30, 50 where available; D_REPS repeats; agree() on autosomal genes ('auto') and all
#   genes ('all'). Output results/within_controls/<dataset>.tsv.gz.
# Model per dataset (data/<dataset>.rds, fields covs and trend): LIBD, limma-trend with the dataset's covariates
#   (scripts/prep_libd.R); Gemma, ~ contrast + age + sex with plain eBayes (fairness_sample_level fit_sz).
# cross: two datasets at the same n per group.
#   kind 'same_subjects' / 'separate_subjects' (LIBD pairs; subjects matched by BrNum; libd_common.R splithalf_rep kinds
#     same_subjects_diff_region / separate_subjects_diff_region): n + n subjects with data in both datasets fitted in both
#     (same subjects), or two subject-separate halves of those subjects, one per dataset. Both sides use the pair's covs.
#   kind 'separate_cohorts' (no shared subjects by design: different brain banks): n + n drawn independently from each
#     dataset, each fitted with its own model (as fairness_sample_level a_split_half.R 'cross').
#   Parts 'sz' (levels as within, A_REPS) and 'controls' (as within, D_REPS). Output results/cross_{sz,controls}/<pair_id>.tsv.gz.
SH <- normalizePath(".")
local({ owd <- setwd("../libd_brainseq/analysis"); on.exit(setwd(owd)); source("scripts/libd_common.R", local = FALSE) })
setwd(SH)
args <- commandArgs(TRUE)
MODE <- args[1]; WHICH <- args[-1]
A_REPS <- as.integer(Sys.getenv("A_REPS", "100"))
D_REPS <- as.integer(Sys.getenv("D_REPS", "50"))
PARTS <- strsplit(Sys.getenv("PARTS", "sz,controls"), ",")[[1]]
LEVELS_A <- c(5, 8, 10, 12, 15, 20, 25, 30, 50, 75)
LEVELS_D <- c(10, 20, 30, 50)
CONTRASTS <- c("sz", "age", "male")
LIBD_CAUDATE_COVS <- c("age", "male", "race_aa", "RIN", "mito")
PAIRS <- list(
  P1_P2_DLPFC = list(a = "LIBD_Phase1_DLPFC", b = "LIBD_Phase2_DLPFC", kinds = c("same_subjects", "separate_subjects"), covs = BASE_COVS),
  P2_DLPFC_Caudate = list(a = "LIBD_Phase2_DLPFC", b = "LIBD_Caudate", kinds = c("same_subjects", "separate_subjects"), covs = LIBD_CAUDATE_COVS),
  P2_DLPFC_GSE208338 = list(a = "LIBD_Phase2_DLPFC", b = "GSE208338", kinds = "separate_cohorts"),
  P2_DLPFC_GSE35977 = list(a = "LIBD_Phase2_DLPFC", b = "GSE35977", kinds = "separate_cohorts"),
  P2_DLPFC_GSE21138 = list(a = "LIBD_Phase2_DLPFC", b = "GSE21138", kinds = "separate_cohorts"),
  P2_DLPFC_GSE80655_DLPFC = list(a = "LIBD_Phase2_DLPFC", b = "GSE80655_DLPFC", kinds = "separate_cohorts"),
  P2_DLPFC_GSE17612 = list(a = "LIBD_Phase2_DLPFC", b = "GSE17612", kinds = "separate_cohorts"),
  P1_DLPFC_GSE208338 = list(a = "LIBD_Phase1_DLPFC", b = "GSE208338", kinds = "separate_cohorts"),
  Caudate_GSE202537_caudate = list(a = "LIBD_Caudate", b = "GSE202537_caudate", kinds = "separate_cohorts"),
  GSE208338_GSE35977 = list(a = "GSE208338", b = "GSE35977", kinds = "separate_cohorts"))

ld <- function(nm) readRDS(file.path("data", paste0(nm, ".rds")))
autos <- function(x) x$genes$symbol[x$genes$autosomal]
use_model <- function(x) { BASE_COVS <<- x$covs; FIT_TREND <<- x$trend }
run_levels <- function(levels, reps, f, tag) {
  res <- list()
  for (i in seq_along(levels)) {
    t0 <- Sys.time()
    rr <- mclapply(seq_len(reps), function(r) f(r, levels[i], i), mc.cores = NCORES, mc.preschedule = FALSE)
    bad <- vapply(rr, function(x) inherits(x, "try-error"), logical(1))
    if (any(bad)) { print(rr[bad][[1]]); stop("repeat failed") }
    res[[i]] <- rbindlist(rr, fill = TRUE)
    cat(tag, levels[i], nrow(res[[i]]), "rows", format(Sys.time() - t0), "\n")
  }
  rbindlist(res, fill = TRUE)
}
out <- function(x, dir, nm) { dir.create(file.path("results", dir), showWarnings = FALSE, recursive = TRUE); fwrite(x, file.path("results", dir, paste0(nm, ".tsv.gz")), sep = "\t") }

feasible <- function(G, levels) { m <- max_n(G); c(as.character(levels[levels <= m]), "max") }

if (MODE == "within") for (nm in WHICH) {
  x <- ld(nm); use_model(x)
  REGIONS <- nm
  D <- setNames(list(x), nm)
  cat(nm, "SZ", sum(x$cv$sz), "control", sum(1 - x$cv$sz), "covs", paste(x$covs, collapse = " + "), "trend", x$trend, "\n")
  if ("sz" %in% PARTS) {
    lv <- feasible(groups_of(x$cv, "sz"), LEVELS_A)
    A <- run_levels(lv, A_REPS, function(r, n, i) splithalf_rep(D, "sz", r, n, level_id = match(n, c(as.character(LEVELS_A), "max")), kinds = "split_halves"), paste(nm, "sz"))
    A[, `:=`(dataset_a = nm, dataset_b = nm)]
    out(A, "within_sz", nm)
  }
  if ("controls" %in% PARTS) {
    GS <- list(auto = autos(x), all = NULL)
    res <- list()
    for (v in CONTRASTS) {
      if (v == "male" && length(unique(x$cv$male)) < 2) { cat(nm, "male: one sex only, skipped\n"); next }
      lv <- intersect(as.character(LEVELS_D), feasible(groups_of(x$cv, v), LEVELS_D))
      if (!length(lv)) next
      res[[v]] <- run_levels(lv, D_REPS, function(r, n, i) splithalf_rep(D, v, r, n, gene_sets = GS, kinds = "split_halves",
                                                                        level_id = 100L + 10L * match(v, CONTRASTS) + match(n, as.character(LEVELS_D))), paste(nm, v))
    }
    C <- rbindlist(res, fill = TRUE); C[, `:=`(dataset_a = nm, dataset_b = nm)]
    out(C, "within_controls", nm)
  }
}

# independent draws from two datasets, each fitted with its own model
cross_rep <- function(A, B, var, r, n, gene_sets, level_id) {
  seed0 <- SEED + 1000000L * level_id + 1000L * r
  ga <- groups_of(A$cv, var); gb <- groups_of(B$cv, var)
  nn <- as.integer(n)
  ha <- draw_halves(ga, nn, seed0 + 21)[[1]]; hb <- draw_halves(gb, nn, seed0 + 22)[[1]]
  fa <- fit_coef(A$E, A$cv[ha], var, covs = A$covs, trend = A$trend)
  fb <- fit_coef(B$E, B$cv[hb], var, covs = B$covs, trend = B$trend)
  if (is.null(fa) || is.null(fb)) return(NULL)
  rbindlist(lapply(names(gene_sets), function(gs) {
    g <- gene_sets[[gs]]
    cbind(data.table(contrast = var, gene_set = gs, kind = "separate_cohorts", region = paste(A$name, B$name, sep = "_vs_"), level = as.character(n),
                     n_per_group = nn, rep = r), suppressWarnings(agree(if (is.null(g)) fa else sub_res(fa, g), if (is.null(g)) fb else sub_res(fb, g))))
  }))
}

if (MODE == "cross") for (pid in WHICH) {
  p <- PAIRS[[pid]]; A <- ld(p$a); B <- ld(p$b)
  cat(pid, p$a, "x", p$b, paste(p$kinds, collapse = ","), "\n")
  GS <- list(auto = intersect(autos(A), autos(B)), all = NULL)
  if (identical(p$kinds, "separate_cohorts")) {
    maxn <- function(v) min(max_n(groups_of(A$cv, v)), max_n(groups_of(B$cv, v)))
    one <- function(v, levels, reps, gs, lid0) {
      lv <- levels[levels <= maxn(v)]
      if (!length(lv)) return(NULL)
      run_levels(lv, reps, function(r, n, i) cross_rep(A, B, v, r, n, gs, lid0 + i), paste(pid, v))
    }
    if ("sz" %in% PARTS) { X <- one("sz", LEVELS_A, A_REPS, list(all = NULL), 300L); X[, `:=`(dataset_a = p$a, dataset_b = p$b, pair = pid)]; out(X, "cross_sz", pid) }
    if ("controls" %in% PARTS) {
      X <- rbindlist(lapply(CONTRASTS, function(v) if (v == "male" && (length(unique(A$cv$male)) < 2 || length(unique(B$cv$male)) < 2)) NULL else
        one(v, LEVELS_D, D_REPS, GS, 400L + 10L * match(v, CONTRASTS))), fill = TRUE)
      X[, `:=`(dataset_a = p$a, dataset_b = p$b, pair = pid)]; out(X, "cross_controls", pid)
    }
  } else {
    # LIBD pair: libd_common.R splithalf_rep between-'region' kinds, with the two datasets in the DLPFC / HIPPO slots and the
    # pair's covariates on both sides
    BASE_COVS <- p$covs; FIT_TREND <- TRUE
    D <- list(DLPFC = A, HIPPO = B)
    B0 <- both_region_donors(D)
    kmap <- c(same_subjects_diff_region = "same_subjects", separate_subjects_diff_region = "separate_subjects")
    kinds <- names(kmap)[kmap %in% p$kinds]
    cat(pid, "subjects in both:", nrow(B0), "SZ", sum(B0$sz), "\n")
    lv_of <- function(v, levels) { G <- groups_of(B0, v); nmax <- if (length(G) == 1) floor(length(G[[1]]) / 2) else min(lengths(G)); as.character(levels[levels <= nmax]) }
    fix <- function(X) { X[, kind := kmap[kind]]; X[, `:=`(dataset_a = p$a, dataset_b = p$b, pair = pid, region = paste(p$a, p$b, sep = "_vs_"))]; X }
    if ("sz" %in% PARTS) {
      X <- run_levels(c(lv_of("sz", LEVELS_A), "max"), A_REPS, function(r, n, i) splithalf_rep(D, "sz", r, n, kinds = kinds, level_id = 500L + i), paste(pid, "sz"))
      out(fix(X), "cross_sz", pid)
    }
    if ("controls" %in% PARTS) {
      X <- rbindlist(lapply(CONTRASTS, function(v) run_levels(lv_of(v, LEVELS_D), D_REPS, function(r, n, i)
        splithalf_rep(D, v, r, n, gene_sets = GS, kinds = kinds, level_id = 600L + 10L * match(v, CONTRASTS) + i), paste(pid, v))), fill = TRUE)
      out(fix(X), "cross_controls", pid)
    }
  }
}
cat("done\n")
