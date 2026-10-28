# B1 item "Cause and manner of death": c16 manner-of-death x score on the recurrent set. Usage: B_OUT=<dir> B_SETF=<recurrent_sets.tsv> Rscript b1_09_c16_manner_ph_score.R
# C16 / ADV-002 donor-level test: manner of death x pH x per-donor 209 recurrent-gene score (set study_brain) in the datasets that record manner per donor.
# Run: Rscript analyses/clinical_reruns_209_2026-10-11/c16_manner_ph_score.R
# Units and loaders as analyses/hbb_manner_of_death_2026-10-11/02_manner_of_death_sz.R (GSE92538 GPL570/GPL96, GSE208338, GSE202537 NAc/caudate/putamen, LIBD DLPFC/HIPPO).
# Per-donor score = mean over the 209 genes present of (SZ direction sign, +1 up / -1 down, 'mixed' dropped) x within-unit gene z-score.
# Tests: (a) score ~ exposure + covariates, within SZ only / controls only / all (dx in model); (b) pH ~ exposure (+age, male) where pH is recorded;
# (c) SZ coefficient on the score with and without the exposure on the same samples. Exposures: suicide, nonnatural. Effects in units of the control SD of the score.
suppressMessages({library(data.table)})
args <- commandArgs(FALSE); me <- sub("--file=", "", args[grep("--file=", args)])
REPO <- normalizePath("schizophrenia-transcriptome-replication"); HERE <- Sys.getenv("B_OUT"); dir.create(HERE, showWarnings = FALSE, recursive = TRUE)
FSL <- file.path(REPO, "analyses/fairness_sample_level"); setwd(FSL); source("scripts/common.R")
EXTRA <- file.path(REPO, "analyses/confounds/ph_adjustment_2026-10-01/results/extra_covariates.tsv")
SETF <- Sys.getenv("B_SETF", file.path(REPO, "analyses/confounds/recurrent_sets.tsv"))
ex <- fread(EXTRA, skip = "dataset_source\t"); ex <- ex[!duplicated(gsm)]
S <- fread(cmd = paste("grep -v '^#'", shQuote(SETF))); S <- S[set == "study_brain" & sz_direction %in% c("up", "down")]
sg <- setNames(ifelse(S$sz_direction == "up", 1, -1), S$gene)
LIBD <- function(r) file.path(REPO, sprintf("analyses/libd_brainseq/analysis/data/libd_%s.rds", r))
units <- list()
for (nm in c("GSE92538_GPL570", "GSE92538_GPL96", "GSE208338", "GSE202537_NAc", "GSE202537_caudate", "GSE202537_putamen")) {
  D <- as.list(REG[dataset == nm]); L <- load_dataset(D, bottom_filter = FALSE)
  cv <- merge(L$cv, ex[, .(gsm, suicide, death_violent, death_nonviolent, manner_nonnatural)], by = "gsm", all.x = TRUE, sort = FALSE)
  cv <- cv[match(colnames(L$E), cv$gsm)]
  cv[, nonnatural := fifelse(!is.na(manner_nonnatural), as.integer(manner_nonnatural),
                      fifelse(!is.na(death_violent) & !is.na(death_nonviolent), as.integer(death_violent == 1 | death_nonviolent == 1), NA_integer_))]
  units[[nm]] <- list(E = L$E, cv = cv, covs = c("age", "male"), cluster = L$cluster)
}
for (r in c("DLPFC", "HIPPO")) {
  x <- readRDS(LIBD(r)); E <- x$E; rownames(E) <- x$genes$symbol; cv <- copy(x$cv)
  cv[, suicide := fifelse(is.na(manner), NA_integer_, as.integer(manner == "Suicide"))]
  cv[, nonnatural := fifelse(is.na(manner), NA_integer_, as.integer(manner != "Natural"))]
  colnames(E) <- cv$gsm
  E <- E[!duplicated(rownames(E)), , drop = FALSE]
  units[[paste0("LIBD_", r)]] <- list(E = E, cv = cv, covs = c("age", "male", "snpPC1", "RIN", "mito"), cluster = "libd")
}
fit1 <- function(y, d, term, covs) {
  X <- as.data.frame(d[, c(term, covs), with = FALSE]); X$y <- y; X <- X[complete.cases(X), ]
  if (nrow(X) < 8 || length(unique(X[[term]])) < 2 || sum(X[[term]] == 1) < 3) return(NULL)
  covs <- covs[vapply(covs, function(v) length(unique(X[[v]])) > 1, TRUE)]
  f <- lm(as.formula(paste("y ~", paste(c(term, covs), collapse = " + "))), data = X)
  s <- summary(f)$coefficients; if (!term %in% rownames(s)) return(NULL)
  tq <- qt(.975, f$df.residual); list(b = s[term, 1], se = s[term, 2], lo = s[term, 1] - tq * s[term, 2], hi = s[term, 1] + tq * s[term, 2], p = s[term, 4], n = nrow(X), nexp = sum(X[[term]] == 1))
}
W <- list(); P <- list(); C <- list(); SC <- list()
for (n in names(units)) { u <- units[[n]]; cv <- u$cv; E <- u$E
  g <- intersect(names(sg), rownames(E)); Z <- t(scale(t(E[g, , drop = FALSE]))); Z <- Z[complete.cases(Z), , drop = FALSE]; g <- rownames(Z)
  cv[, score := colMeans(Z * sg[g])]; sdc <- sd(cv[sz == 0, score])
  SC[[n]] <- data.table(unit = n, genes_scored = length(g), n = nrow(cv), n_sz = sum(cv$sz), sd_control_score = sdc, mean_sz_minus_control = mean(cv[sz == 1, score]) - mean(cv[sz == 0, score]), ph_recorded = sum(!is.na(cv$pH)), pmi_recorded = sum(!is.na(cv$PMI)))
  for (exn in c("suicide", "nonnatural")) { if (all(is.na(cv[[exn]]))) next
    for (grp in c("SZ only", "controls only", "all, dx in model")) {
      sel <- switch(grp, "SZ only" = cv$sz == 1, "controls only" = cv$sz == 0, rep(TRUE, nrow(cv))) & !is.na(cv[[exn]]); d <- cv[sel]
      for (sc in c("usual", "usual+PMI+pH")) {
        covs <- c(u$covs, if (grp == "all, dx in model") "sz")
        if (sc != "usual") { c2 <- covs; for (v in c("PMI", "pH")) if (sum(!is.na(d[[v]])) >= 8) c2 <- c(c2, v); if (identical(c2, covs)) next; covs <- c2 }
        r <- fit1(d$score, d, exn, covs); if (is.null(r)) next
        W[[length(W) + 1]] <- data.table(unit = n, exposure = exn, group = grp, covariates = sc, n = r$n, n_exposed = r$nexp, est_score = r$b, lo = r$lo, hi = r$hi, p = r$p, est_in_control_SD = r$b / sdc, lo_SD = r$lo / sdc, hi_SD = r$hi / sdc)
      }
      if (sum(!is.na(d$pH)) >= 8) { dd <- d[!is.na(pH)]; r <- fit1(dd$pH, dd, exn, c("age", "male", if (grp == "all, dx in model") "sz"))
        if (!is.null(r)) P[[length(P) + 1]] <- data.table(unit = n, exposure = exn, group = grp, n = r$n, n_exposed = r$nexp, est_pH = r$b, lo = r$lo, hi = r$hi, p = r$p) }
    }
    d <- cv[!is.na(get(exn))]
    for (sc in c("usual", "usual+PMI+pH")) { covs <- u$covs
      if (sc != "usual") { c2 <- covs; for (v in c("PMI", "pH")) if (sum(!is.na(d[[v]])) >= 8) c2 <- c(c2, v); if (identical(c2, covs)) next; covs <- c2 }
      X <- as.data.frame(d[, unique(c("sz", exn, covs, "score")), with = FALSE]); X <- X[complete.cases(X), ]
      if (length(unique(X$sz)) < 2 || sum(X[[exn]] == 1) < 3) next
      cok <- covs[vapply(covs, function(v) length(unique(X[[v]])) > 1, TRUE)]
      f0 <- lm(as.formula(paste("score ~ sz +", paste(cok, collapse = " + "))), X); f1 <- lm(as.formula(paste("score ~ sz +", exn, "+", paste(cok, collapse = " + "))), X)
      s0 <- summary(f0)$coefficients["sz", ]; s1 <- summary(f1)$coefficients["sz", ]
      C[[length(C) + 1]] <- data.table(unit = n, exposure = exn, covariates = sc, n = nrow(X), n_sz = sum(X$sz), n_exposed = sum(X[[exn]]), sz_base = s0[1], sz_base_se = s0[2], sz_adj = s1[1], sz_adj_se = s1[2], exposure_coef_in_adj = coef(f1)[exn], pct_sz_removed = 100 * (1 - s1[1] / s0[1])) }
  }
}
W <- rbindlist(W); P <- rbindlist(P); C <- rbindlist(C); SC <- rbindlist(SC)
dl <- function(b, se) { k <- length(b); w <- 1 / se^2; fe <- sum(w * b) / sum(w); Q <- sum(w * (b - fe)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + tau2); re <- sum(wr * b) / sum(wr); sr <- sqrt(1 / sum(wr))
  data.table(k = k, re = re, lo = re - 1.96 * sr, hi = re + 1.96 * sr, p = 2 * pnorm(-abs(re / sr)), I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0) }
INDEP <- c("GSE92538_GPL570", "GSE208338", "GSE202537_caudate", "LIBD_DLPFC")
W[, se := (hi - lo) / (2 * 1.96)]
M <- rbindlist(lapply(c("suicide", "nonnatural"), function(e) rbindlist(lapply(c("SZ only", "controls only", "all, dx in model"), function(g) {
  w <- W[unit %in% INDEP & exposure == e & group == g & covariates == "usual"]; if (nrow(w) < 2) return(NULL)
  r <- dl(w$est_in_control_SD, w$se / w$sd); cbind(exposure = e, group = g, units = paste(w$unit, collapse = ","), r) }))), fill = TRUE)
# se in SD units
W[, se_SD := (hi_SD - lo_SD) / (2 * 1.96)]
M <- rbindlist(lapply(c("suicide", "nonnatural"), function(e) rbindlist(lapply(c("SZ only", "controls only", "all, dx in model"), function(g) {
  w <- W[unit %in% INDEP & exposure == e & group == g & covariates == "usual"]; if (nrow(w) < 2) return(NULL)
  cbind(exposure = e, group = g, units = paste(w$unit, collapse = ","), dl(w$est_in_control_SD, w$se_SD)) }))), fill = TRUE)
PG <- list()
for (n in names(units)) { u <- units[[n]]; cv <- u$cv; if (sum(!is.na(cv$pH)) < 8) next
  for (exn in c("suicide", "nonnatural")) { if (all(is.na(cv[[exn]]))) next
    X <- as.data.frame(cv[, c("pH", "sz", exn, "age", "male"), with = FALSE]); X <- X[complete.cases(X), ]; if (sum(X[[exn]] == 1) < 3) next
    f0 <- lm(pH ~ sz + age + male, X); f1 <- lm(as.formula(paste("pH ~ sz +", exn, "+ age + male")), X)
    c0 <- confint(f0)["sz", ]; c1 <- confint(f1)["sz", ]
    PG[[length(PG) + 1]] <- data.table(unit = n, exposure = exn, n = nrow(X), n_sz = sum(X$sz), share_exposed_sz = mean(X[[exn]][X$sz == 1]), share_exposed_ctl = mean(X[[exn]][X$sz == 0]),
      pH_gap_sz_minus_ctl = coef(f0)["sz"], gap_lo = c0[1], gap_hi = c0[2], pH_gap_adj_for_exposure = coef(f1)["sz"], adj_lo = c1[1], adj_hi = c1[2], exposure_coef_pH = coef(f1)[exn]) } }
PG <- rbindlist(PG)
sha256f <- function(f) substr(sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE)), 1, 12)
ins <- c(EXTRA, SETF, LIBD("DLPFC"), LIBD("HIPPO"), file.path(FSL, "scripts/common.R"), file.path(REPO, "analyses/controls_age_sex/data/datasets.tsv"))
hdr <- function(note) c("# script: factcheck_resolution_2026-10-19/B/scripts/b1_09_c16_manner_ph_score.R (copy of analyses/clinical_reruns_209_2026-10-11/c16_manner_ph_score.R; only REPO, HERE and SETF changed)", "# run: 2026-10-19; no random draws in this script (seed n/a); base: set per B_SETF",
  sprintf("# input: %s sha256 %s", sub(paste0(REPO, "/"), "", ins), vapply(ins, sha256f, "")), "# score: mean over the study_brain genes of B_SETF (mixed dropped) of SZ-direction sign x within-unit gene z; expression log2 (LIBD TMM logCPM)", paste0("# ", note))
put <- function(x, f, note) { writeLines(hdr(note), file.path(HERE, f)); fwrite(x[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)], file.path(HERE, f), sep = "\t", append = TRUE, col.names = TRUE) }
put(SC, "c16_unit_scores.tsv", "per-unit score summary: genes scored, n, SD of control scores, SZ minus control mean score, donors with pH / PMI recorded")
put(W, "c16_manner_score_within_group.tsv", "score difference exposed minus not exposed (lm; est_score in score units, *_SD in units of the control SD of the score); covariates usual or usual+PMI+pH where >= 8 recorded; n_exposed = donors with the exposure")
put(P, "c16_manner_ph.tsv", "pH difference exposed minus not exposed (lm with age, male [, dx]) in units with recorded pH")
put(C, "c16_sz_coefficient_with_without_manner.tsv", "SZ-minus-control score coefficient on the same samples without (sz_base) and with (sz_adj) the exposure")
put(PG, "c16_ph_gap_with_without_manner.tsv", "SZ-minus-control pH gap (lm with age, male) on the same donors without and with the exposure; share_exposed = fraction of the diagnosis group with the exposure; GSE202537 regions share donors")
put(M, "c16_manner_meta.tsv", "DerSimonian-Laird pooled exposure effect (control-SD units) over independent units GSE92538_GPL570, GSE208338, GSE202537_caudate, LIBD_DLPFC (one per donor cluster); these are not the same donors as one another by cluster assignment")
options(width = 250); print(SC); print(W[covariates == "usual"][order(exposure, group, unit)]); print(P); print(PG); print(C[covariates == "usual"]); print(M)
