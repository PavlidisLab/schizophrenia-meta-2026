# Do recurrent "DE in schizophrenia" genes track brain pH / PMI / RIN / agonal state in public data?
#
# Run from confounds/:  Rscript scripts/within_dataset_confounds.R
# Inputs (all read only):
#   data/gemma/<id>_processed.tsv.gz  Gemma processed expression (anonymous GET, scripts/fetch_gemma_processed.py)
#   data/covariates/<dataset>.tsv     per-sample covariates from GEO (scripts/build_covariates.R)
#   recurrent_sets.tsv                recurrent genes and their SZ direction (scripts/recurrent_sets.R)
#   ../../data/reference/hgnc_complete_set.txt  NCBI gene id -> current HGNC symbol
# Outputs:
#   results/within_controls_<dataset>.tsv   per gene: factor t, p (controls only, ~ factor + age + sex)
#   results/sz_adjust_<dataset>.tsv         per gene: SZ log2FC / t / p without and with pH (and with all quality covariates)
#   results/within_dataset_summary.tsv      one row per dataset x factor x recurrent set
#   results/group_differences.tsv           SZ vs control difference in pH, PMI, RIN, agonal, age
#
# Methods, in short. Per dataset (and per region where one GEO series holds several), probes annotated to one
# gene are mapped to HGNC via NCBI gene id; per gene the probe with the highest mean is kept; genes in the
# bottom tenth of mean expression are dropped. Arrays from the same subject (GSE92538) are averaged.
# Controls only: limma, expression ~ factor + age + sex, one model per factor (pH, PMI, RIN, agonal).
# Enrichment: count of recurrent genes at nominal p < 0.05 vs the genome-wide rate (Fisher's exact test,
# recurrent vs all other genes) and vs 10,000 random gene sets matched on mean-expression decile.
# Direction: for pH, a recurrent gene is "confound-consistent" when its pH coefficient has the opposite sign to
# its SZ direction (up in SZ and higher at low pH, or down in SZ and lower at low pH), which is what a lower pH in
# cases would produce. Signed score = mean over recurrent genes of (-sz_sign x t_pH); matched null as above.
# Cases vs controls: limma, expression ~ dx + age + sex, then + pH, then + all available quality covariates.
# Signed attenuation: signed_ratio_* = sum over the set's genes of (SZ sign x log2FC after adjustment) / the same sum
# before adjustment (a pooled, sign-aware ratio; median_ratio_pH, the median of per-gene ratios, is kept for comparison).
# A factor with 4 or fewer distinct values (agonal) is tested only when at least 5 samples differ from its most
# common value (GSE92538 GPL570 controls: 2 of 56 non-zero, not testable).
suppressMessages({library(limma); library(data.table)})
set.seed(20260924)
# CONFOUNDS_SETS / CONFOUNDS_OUT: another gene-set file (recurrent_sets.tsv layout) and output directory, for side runs
# such as careful_core_2026-10-01/; unset, the main run's recurrent_sets.tsv and results/ are used.
SETS_FILE <- Sys.getenv("CONFOUNDS_SETS", "recurrent_sets.tsv")
OUT <- Sys.getenv("CONFOUNDS_OUT", "results")
dir.create(OUT, showWarnings = FALSE)
# CONFOUNDS_FACTORS: the control-only factors, comma-separated (unset: pH,PMI,RIN,agonal, the main run). Factor age is fitted
# as expr ~ age + sex, without age again as a covariate (age_genes_2026-10-01/). The signed score keeps its pH convention,
# mean(-SZ sign x t), so for age a negative score means SZ-up genes rise with age.
FACTORS <- strsplit(Sys.getenv("CONFOUNDS_FACTORS", "pH,PMI,RIN,agonal"), ",")[[1]]
ctl_design <- function(d, fac) {
  covs <- c(if (fac != "age") "age", if (length(unique(d$sex)) > 1) "sex")
  model.matrix(as.formula(paste(c("~ d[[fac]]", covs), collapse = " + ")), data = d)
}
NPERM <- 10000
# Sets other than these four draw their random numbers from a stream of their own (seeded from set, dataset and factor)
# and hand the global stream back as they found it, so adding a set (study_brain_5plus, 2026-09-29) leaves the matched
# and permutation p values of the earlier sets unchanged.
LEGACY_SETS <- c("cohort_brain", "study_brain", "cohort_brainFDR", "study_brainFDR")
own_stream <- function(S, ...) {
  if (S %in% LEGACY_SETS) return(NULL)
  saved <- get(".Random.seed", envir = globalenv())
  key <- utf8ToInt(paste(S, ..., sep = "|"))
  set.seed(sum(key * seq_along(key)) %% 2147483647)
  saved
}
restore_stream <- function(saved) if (!is.null(saved)) assign(".Random.seed", saved, envir = globalenv())

hg <- fread("../../data/reference/hgnc_complete_set.txt", select = c("symbol", "entrez_id"), quote = "")
hg <- hg[!is.na(entrez_id)]
ent2sym <- setNames(hg$symbol, as.character(hg$entrez_id))

rs <- fread(SETS_FILE, skip = "set\t")
rs[, sz_sign := ifelse(sz_direction == "up", 1, ifelse(sz_direction == "down", -1, 0))]
SETS <- unique(rs$set)

DATASETS <- list(
  list(name = "GSE21138", gemma = 1880, cov = "GSE21138"),
  list(name = "GSE35977", gemma = 5949, cov = "GSE35977"),
  list(name = "GSE35974", gemma = 5939, cov = "GSE35974"),
  list(name = "GSE53987_BA46", gemma = 8359, cov = "GSE53987", region = "Pre-frontal cortex (BA46)"),
  list(name = "GSE53987_hippocampus", gemma = 8359, cov = "GSE53987", region = "hippocampus"),
  list(name = "GSE53987_striatum", gemma = 8359, cov = "GSE53987", region = "Associative striatum"),
  list(name = "GSE17612", gemma = 1381, cov = "GSE17612"),
  list(name = "GSE21935", gemma = 4150, cov = "GSE21935"),
  list(name = "GSE92538_GPL570", gemma = 18742, cov = "GSE92538_GPL10526"),
  list(name = "GSE92538_GPL96", gemma = 18741, cov = "GSE92538_GPL17027"),
  list(name = "GSE80655_DLPFC", gemma = 13014, cov = "GSE80655", region = "DLPFC"),
  list(name = "GSE80655_AnCg", gemma = 13014, cov = "GSE80655", region = "AnCg"),
  list(name = "GSE80655_nAcc", gemma = 13014, cov = "GSE80655", region = "nAcc"),
  list(name = "GSE208338", gemma = 26880, cov = "GSE208338"),
  # GSE208338 without the donors it shares with GSE21138 (same diagnosis, sex, pH, PMI and age within 1 year)
  list(name = "GSE208338_noshared", gemma = 26880, cov = "GSE208338", drop_shared_with = "GSE21138")
)
args <- commandArgs(TRUE)
if (length(args)) DATASETS <- Filter(function(d) d$name %in% args, DATASETS)

.cache <- new.env()
read_gemma <- function(id) {
  key <- as.character(id)
  if (!is.null(.cache[[key]])) return(.cache[[key]])
  x <- fread(cmd = paste("gzip -dc", shQuote(sprintf("data/gemma/%s_processed.tsv.gz", id)), "| grep -v '^#'"), sep = "\t", quote = "")
  x <- x[!is.na(NCBIid) & !grepl("\\|", as.character(NCBIid)) & as.character(NCBIid) != ""]
  sym <- ent2sym[as.character(x$NCBIid)]
  keep <- !is.na(sym)
  sc <- colnames(x)[-(1:6)]
  m <- as.matrix(x[keep, ..sc]); rownames(m) <- sym[keep]
  if (all(grepl("GSM[0-9]+", sc))) {
    colnames(m) <- regmatches(sc, regexpr("GSM[0-9]+", sc))
  } else {  # older export: map BioAssayId to GEO accession via /datasets/{id}/samples
    sj <- rjson::fromJSON(file = sprintf("data/gemma/%s_samples.json", id))$data
    acc <- setNames(sapply(sj, function(s) s$accession$accession), sapply(sj, function(s) as.character(s$id)))
    colnames(m) <- acc[sub(".*BioAssayId=([0-9]+).*", "\\1", sc)]
  }
  .cache[[key]] <- m
  m
}

gene_level <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  o <- order(rownames(m), -mu)
  m <- m[o, , drop = FALSE]
  m[!duplicated(rownames(m)), , drop = FALSE]
}

fit_factor <- function(E, cv, fac) {
  ok <- !is.na(cv[[fac]]) & !is.na(cv$age) & !is.na(cv$sex)
  if (sum(ok) < 10 || length(unique(cv[[fac]][ok])) < 2 || (length(unique(cv[[fac]][ok])) < 3 && fac != "agonal")) return(NULL)
  # a factor with few distinct values (agonal) needs at least 5 samples away from its most common value
  v <- cv[[fac]][ok]
  if (length(unique(v)) <= 4 && sum(v != as.numeric(names(which.max(table(v))))) < 5) return(NULL)
  d <- cv[ok]
  des <- ctl_design(d, fac)
  f <- eBayes(lmFit(E[, d$gsm], des))
  tt <- topTable(f, coef = 2, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), coef = tt$logFC, t = tt$t, p = tt$P.Value, n = nrow(d))
}

fit_dx <- function(E, cv, extra) {
  vars <- c("age", "sex", extra)
  ok <- complete.cases(cv[, ..vars])
  d <- cv[ok]
  if (sum(d$dx == "control") < 5 || sum(d$dx == "schizophrenia") < 5) return(NULL)
  d[, sz := as.integer(dx == "schizophrenia")]
  fo <- paste("~ sz + age", if (length(unique(d$sex)) > 1) "+ sex" else "", if (length(extra)) paste("+", paste(extra, collapse = " + ")) else "")
  des <- model.matrix(as.formula(fo), data = d)
  f <- eBayes(lmFit(E[, d$gsm], des))
  tt <- topTable(f, coef = "sz", number = Inf, sort.by = "none")
  list(tab = data.table(gene = rownames(tt), logfc = tt$logFC, t = tt$t, p = tt$P.Value), n_sz = sum(d$sz), n_ctl = sum(1 - d$sz))
}

matched_null <- function(stat, bins, idx, fun) {
  # stat: named vector over all genes; idx: recurrent gene indices; fun(values, which) -> scalar
  by_bin <- split(seq_along(stat), bins)
  need <- table(bins[idx])
  replicate(NPERM, {
    s <- unlist(lapply(names(need), function(b) sample(by_bin[[b]], need[[b]])))
    fun(s)
  })
}

# Self-contained test that keeps gene-gene correlation: shuffle the factor among controls (age and sex stay with
# their samples), refit ordinary least squares for the set's genes, recompute the signed score.
NSAMPPERM <- 2000
ols_t <- function(Y, X, j) {  # Y samples x genes; t for column j of X
  q <- qr(X); b <- qr.coef(q, Y); r <- qr.resid(q, Y)
  s2 <- colSums(r^2) / (nrow(X) - q$rank)
  v <- chol2inv(qr.R(q))[j, j]
  b[j, ] / sqrt(s2 * v)
}
perm_score <- function(E, cv, fac, genes_idx, sg) {
  ok <- !is.na(cv[[fac]]) & !is.na(cv$age) & !is.na(cv$sex)
  d <- cv[ok]
  X <- ctl_design(d, fac)
  Y <- t(E[genes_idx, d$gsm, drop = FALSE])
  obs <- mean(-sg * ols_t(Y, X, 2))
  f <- d[[fac]]
  null <- replicate(NSAMPPERM, { X[, 2] <- sample(f); mean(-sg * ols_t(Y, X, 2)) })
  c(score_ols = obs, perm_p_greater = (1 + sum(null >= obs)) / (NSAMPPERM + 1), perm_p_less = (1 + sum(null <= obs)) / (NSAMPPERM + 1),
    perm_z = (obs - mean(null)) / sd(null))
}

summ <- list(); grp <- list(); med <- list()
for (D in DATASETS) {
  m <- read_gemma(D$gemma)
  cv <- fread(sprintf("data/covariates/%s.tsv", D$cov), skip = "gsm\t", colClasses = list(character = c("gsm", "subject", "sex", "dx", "region", "cohort")))
  if (!is.null(D$region)) cv <- cv[region == D$region]
  if (!is.null(D$drop_shared_with)) {
    o <- fread(sprintf("data/covariates/%s.tsv", D$drop_shared_with), skip = "gsm\t")
    key <- function(d) paste(d$dx, d$sex, sprintf("%.2f", d$pH), sprintf("%.1f", d$PMI))
    ko <- key(o)
    hit <- vapply(seq_len(nrow(cv)), function(i) any(ko == key(cv[i]) & abs(o$age - cv$age[i]) < 1), logical(1))
    fwrite(data.table(dataset = D$name, gsm = cv$gsm[hit], dx = cv$dx[hit], sex = cv$sex[hit], age = cv$age[hit], pH = cv$pH[hit], PMI = cv$PMI[hit]),
           file.path(OUT, sprintf("shared_donors_%s.tsv", D$name)), sep = "\t")
    cat(D$name, ": dropping", sum(hit), "samples matching a", D$drop_shared_with, "donor (", paste(names(table(cv$dx[hit])), table(cv$dx[hit]), collapse = ", "), ")\n")
    cv <- cv[!hit]
  }
  cv <- cv[gsm %in% colnames(m) & dx %in% c("control", "schizophrenia")]
  E <- m[, cv$gsm, drop = FALSE]
  # samples Gemma left without values (outliers removed in Gemma) are dropped
  okS <- colMeans(is.na(E)) < 0.5
  if (any(!okS)) cat(D$name, ": dropping", sum(!okS), "samples with no Gemma values\n")
  E <- E[, okS, drop = FALSE]; cv <- cv[okS]
  # average arrays from the same subject
  if (any(!is.na(cv$subject) & cv$subject != "") && any(duplicated(cv$subject))) {
    E <- t(apply(E, 1, function(v) tapply(v, cv$subject, mean, na.rm = TRUE)))
    cv <- cv[!duplicated(subject)]
    E <- E[, cv$subject, drop = FALSE]; colnames(E) <- cv$gsm
  }
  E <- gene_level(E)
  E <- E[rowSums(is.na(E)) == 0, , drop = FALSE]
  mu <- rowMeans(E)
  E <- E[mu > quantile(mu, 0.10), , drop = FALSE]
  mu <- rowMeans(E)
  bins <- cut(rank(mu), 10, labels = FALSE)
  genes <- rownames(E)
  cat(sprintf("%s: %d genes, %d controls, %d SZ\n", D$name, nrow(E), sum(cv$dx == "control"), sum(cv$dx == "schizophrenia")))

  # group differences
  for (v in c("pH", "PMI", "RIN", "agonal", "age")) {
    a <- cv[dx == "schizophrenia"][[v]]; b <- cv[dx == "control"][[v]]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) >= 3 && length(b) >= 3 && length(unique(c(a, b))) > 1)
      grp[[length(grp) + 1]] <- data.table(dataset = D$name, variable = v, n_sz = length(a), n_ctl = length(b), mean_sz = mean(a), mean_ctl = mean(b),
                                           diff = mean(a) - mean(b), t_p = t.test(a, b)$p.value, wilcox_p = wilcox.test(a, b, exact = FALSE)$p.value)
  }

  # within controls
  ctl <- cv[dx == "control"]
  wc <- list()
  for (fac in FACTORS) {
    r <- fit_factor(E, ctl, fac)
    if (is.null(r)) next
    r[, factor := fac]
    r[, expr_percentile := round(rank(mu)[match(gene, genes)] / length(mu), 4)]
    wc[[fac]] <- r
    for (S in SETS) {
      rr <- rs[set == S]
      idx <- match(rr$gene, genes); keep <- !is.na(idx); idx <- idx[keep]; sg <- rr$sz_sign[keep]
      if (!length(idx)) next
      sv <- own_stream(S, D$name, fac)
      sig <- r$p < 0.05
      pi0 <- mean(sig)
      obs <- sum(sig[idx])
      ft <- fisher.test(matrix(c(obs, length(idx) - obs, sum(sig) - obs, sum(!sig) - (length(idx) - obs)), 2), alternative = "greater")
      nullc <- matched_null(r$t, bins, idx, function(s) sum(sig[s]))
      score <- mean(-sg * r$t[idx])
      nulls <- matched_null(r$t, bins, idx, function(s) mean(-sg * r$t[s]))
      cons <- sum(sig[idx] & sign(r$t[idx]) == -sg)
      ps <- perm_score(E, ctl, fac, idx, sg)
      summ[[length(summ) + 1]] <- data.table(
        dataset = D$name, factor = fac, set = S, n_ctl = r$n[1], n_genes = length(genes), n_set_measured = length(idx),
        genome_rate_p05 = round(pi0, 4), n_set_p05 = obs, expected_p05 = round(pi0 * length(idx), 2), fisher_p = signif(ft$p.value, 3),
        matched_p = signif(mean(nullc >= obs), 3), matched_mean = round(mean(nullc), 2),
        # same sign counted directly: obs - cons also counted genes with mixed SZ direction (sz_sign 0)
        n_p05_opposite_sign_to_sz = cons, n_p05_same_sign_as_sz = sum(sig[idx] & sg != 0 & sign(r$t[idx]) == sg),
        signed_score = round(score, 3), signed_score_matched_p = signif(mean(nulls >= score), 3),
        score_ols = round(ps[["score_ols"]], 3), sample_perm_p_greater = signif(ps[["perm_p_greater"]], 3), sample_perm_p_less = signif(ps[["perm_p_less"]], 3),
        sample_perm_z = round(ps[["perm_z"]], 3),
        genes_p05 = paste(sprintf("%s(%s,t=%.1f)", genes[idx][sig[idx]], ifelse(sg[sig[idx]] > 0, "SZup", ifelse(sg[sig[idx]] < 0, "SZdown", "SZmixed")), r$t[idx][sig[idx]]), collapse = ";"))
      restore_stream(sv)
    }
  }
  if (length(wc)) fwrite(rbindlist(wc), file.path(OUT, sprintf("within_controls_%s.tsv", D$name)), sep = "\t")

  # medication: within SZ, expression ~ chlorpromazine-equivalent dose + age + sex (datasets that record dose)
  szc <- cv[dx == "schizophrenia" & !is.na(cpz_mg)]
  if (nrow(szc) >= 10) {
    r <- fit_factor(E, szc, "cpz_mg")
    fwrite(r, file.path(OUT, sprintf("within_sz_cpz_%s.tsv", D$name)), sep = "\t")
    for (S in SETS) {
      rr <- rs[set == S]; idx <- match(rr$gene, genes); keep <- !is.na(idx); idx <- idx[keep]; sg <- rr$sz_sign[keep]
      sig <- r$p < 0.05; obs <- sum(sig[idx])
      sv <- own_stream(S, D$name, "cpz_mg")
      ps <- perm_score(E, szc, "cpz_mg", idx, sg)
      restore_stream(sv)
      med[[length(med) + 1]] <- data.table(dataset = D$name, factor = "cpz_mg (within SZ)", set = S, n_sz = r$n[1], n_set_measured = length(idx),
        genome_rate_p05 = round(mean(sig), 4), n_set_p05 = obs, expected_p05 = round(mean(sig) * length(idx), 2),
        fisher_p = signif(fisher.test(matrix(c(obs, length(idx) - obs, sum(sig) - obs, sum(!sig) - (length(idx) - obs)), 2), alternative = "greater")$p.value, 3),
        n_p05_same_sign_as_sz = sum(sig[idx] & sg != 0 & sign(r$t[idx]) == sg),
        score_dose_follows_sz = round(-ps[["score_ols"]], 3), sample_perm_p_dose_follows_sz = signif(ps[["perm_p_less"]], 3),
        genes_p05 = paste(sprintf("%s(%s,t=%.1f)", genes[idx][sig[idx]], ifelse(sg[sig[idx]] > 0, "SZup", ifelse(sg[sig[idx]] < 0, "SZdown", "SZmixed")), r$t[idx][sig[idx]]), collapse = ";"))
    }
  }

  # cases vs controls, with and without quality covariates
  qual <- intersect(c("pH", "PMI", "RIN", "agonal"), names(which(sapply(cv[, .(pH, PMI, RIN, agonal)], function(x) sum(!is.na(x)) >= 0.8 * nrow(cv)))))
  base_vars <- if ("pH" %in% qual) "pH" else character(0)
  # same samples in every model: those with every covariate used
  vars <- c("age", "sex", qual)
  cc <- cv[complete.cases(cv[, ..vars])]
  A <- fit_dx(E, cc, character(0)); B <- if (length(base_vars)) fit_dx(E, cc, "pH") else NULL; C <- if (length(qual)) fit_dx(E, cc, qual) else NULL
  Ag <- if ("agonal" %in% qual) fit_dx(E, cc, "agonal") else NULL
  if (!is.null(A)) {
    out <- A$tab[, .(gene, logfc_base = logfc, t_base = t, p_base = p)]
    if (!is.null(B)) out <- cbind(out, B$tab[, .(logfc_pH = logfc, t_pH = t, p_pH = p)])
    if (!is.null(C)) out <- cbind(out, C$tab[, .(logfc_allq = logfc, t_allq = t, p_allq = p)])
    if (!is.null(Ag)) out <- cbind(out, Ag$tab[, .(logfc_agonal = logfc, t_agonal = t, p_agonal = p)])
    fwrite(out, file.path(OUT, sprintf("sz_adjust_%s.tsv", D$name)), sep = "\t")
    for (S in SETS) {
      rr <- rs[set == S]
      idx <- match(rr$gene, out$gene); keep <- !is.na(idx); idx <- idx[keep]; sg <- rr$sz_sign[keep]
      if (!length(idx)) next
      row <- data.table(dataset = D$name, factor = "SZ_effect", set = S, n_ctl = A$n_ctl, n_sz = A$n_sz, n_set_measured = length(idx),
                        adj_covariates = paste(qual, collapse = "+"),
                        median_signed_logfc_base = round(median(sg * out$logfc_base[idx]), 3),
                        n_consistent_p05_base = sum(out$p_base[idx] < 0.05 & sign(out$t_base[idx]) == sg),
                        genome_rate_p05_base = round(mean(out$p_base < 0.05), 4))
      row[, signed_sum_base := round(sum(sg * out$logfc_base[idx]), 4)]
      if (!is.null(B)) row[, `:=`(signed_sum_pH = round(sum(sg * out$logfc_pH[idx]), 4),
                                   signed_ratio_pH = round(sum(sg * out$logfc_pH[idx]) / sum(sg * out$logfc_base[idx]), 3))]
      if (!is.null(C)) row[, `:=`(signed_sum_allq = round(sum(sg * out$logfc_allq[idx]), 4),
                                   signed_ratio_allq = round(sum(sg * out$logfc_allq[idx]) / sum(sg * out$logfc_base[idx]), 3))]
      if (!is.null(Ag)) row[, signed_ratio_agonal := round(sum(sg * out$logfc_agonal[idx]) / sum(sg * out$logfc_base[idx]), 3)]
      if (!is.null(B)) row[, `:=`(median_signed_logfc_pH = round(median(sg * out$logfc_pH[idx]), 3),
                                   n_consistent_p05_pH = sum(out$p_pH[idx] < 0.05 & sign(out$t_pH[idx]) == sg),
                                   median_ratio_pH = round(median(out$logfc_pH[idx] / out$logfc_base[idx]), 3),
                                   genome_rate_p05_pH = round(mean(out$p_pH < 0.05), 4))]
      if (!is.null(Ag)) row[, `:=`(median_signed_logfc_agonal = round(median(sg * out$logfc_agonal[idx]), 3),
                                    n_consistent_p05_agonal = sum(out$p_agonal[idx] < 0.05 & sign(out$t_agonal[idx]) == sg))]
      row[, n_samples_models := nrow(cc)]
      if (!is.null(C)) row[, `:=`(median_signed_logfc_allq = round(median(sg * out$logfc_allq[idx]), 3),
                                   n_consistent_p05_allq = sum(out$p_allq[idx] < 0.05 & sign(out$t_allq[idx]) == sg),
                                   genome_rate_p05_allq = round(mean(out$p_allq < 0.05), 4))]
      summ[[length(summ) + 1]] <- row
    }
  }
}
S <- rbindlist(summ, fill = TRUE)
fh <- file(file.path(OUT, "within_dataset_summary.tsv"), "w")
writeLines(c(paste("# scripts/within_dataset_confounds.R run", Sys.time()),
             paste("# input:", SETS_FILE, "sha256", substr(system(paste("shasum -a 256", shQuote(SETS_FILE)), intern = TRUE), 1, 12)),
             paste0("# matched null: ", NPERM, " random gene sets matched on mean-expression decile; genes: HGNC via NCBI id, top probe per gene, bottom tenth by mean dropped"),
             "# rows with factor pH/PMI/RIN/agonal: controls only, expr ~ factor + age + sex (limma). n_p05_opposite_sign_to_sz: genes at p < 0.05 whose factor coefficient has the opposite sign to their SZ direction",
             "# rows with factor SZ_effect: SZ + controls, expr ~ dx + age + sex (base), + pH, + all quality covariates (allq); same samples in all three models; signed = multiplied by the recurrent SZ direction"), fh)
close(fh)
fwrite(S, file.path(OUT, "within_dataset_summary.tsv"), sep = "\t", append = TRUE, col.names = TRUE)
fwrite(rbindlist(grp), file.path(OUT, "group_differences.tsv"), sep = "\t")
if (length(med)) fwrite(rbindlist(med), file.path(OUT, "medication_dose.tsv"), sep = "\t")
print(S[factor != "SZ_effect", .(dataset, factor, set, n_ctl, n_set_measured, n_set_p05, expected_p05, fisher_p, matched_p, n_p05_opposite_sign_to_sz, signed_score, signed_score_matched_p)])
print(S[factor == "SZ_effect"])
print(rbindlist(grp))
