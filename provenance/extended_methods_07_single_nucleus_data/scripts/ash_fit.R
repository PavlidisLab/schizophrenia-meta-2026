#!/usr/bin/env Rscript
# Views B and C per source x group (cell type, region or "all") from effect_sizes/data/full/*.tsv.gz,
# written by harmonize.py.
#
#   Rscript effect_sizes/ash_fit.R   -> effect_sizes/data/by_group.tsv
#
# View B: |log2FC| of every tested gene (median, 90th percentile), no SE needed.
# View C: ashr (Stephens 2017) with mixcompdist = "normal" and a point mass at 0, fitted to
#   (log2FC, SE) of every gene in one group. SE = the table's SE; else log2FC / t; else
#   |log2FC| / qt(1 - p/2, df) (qnorm when df is missing). Reported: SD of the fitted prior g
#   (the estimated distribution of true effects, mean 0), its null weight pi0, the mean |posterior
#   mean|, and the fraction of genes with lfsr < 0.05.
# Likelihood (analysis review E3): where the SE is backed out from a t statistic or a p value with a finite df and
#   the group's median df is below DF_T, ash uses a t likelihood (lik_t(df)), which needs a uniform mixture
#   (mixcompdist = "uniform", point mass at 0, null-biased prior). df: the table's df column where given (Student
#   t tests in our reproductions); else, when the table has both t and p, the df at which 2 * pt(-|t|, df) = p,
#   solved per gene over genes with 0.001 <= p <= 0.5 (limma moderated df for Gemma; a DESeq2 Wald statistic
#   solves to a very large df and keeps the normal likelihood). Otherwise the normal likelihood and normal mixture
#   as before. SD of a uniform component [a, b]: variance (b - a)^2 / 12 around (a + b) / 2.
# Groups with fewer than MIN_GENES genes with a usable SE are not fitted.
suppressMessages({ library(ashr); library(parallel) })
set.seed(1)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
full_dir <- file.path(here, "data", "full")
MIN_GENES <- 500
DF_T <- 30

se_of <- function(d) {
  se <- suppressWarnings(as.numeric(d$se))
  t <- suppressWarnings(as.numeric(d$t))
  p <- suppressWarnings(as.numeric(d$p))
  df <- suppressWarnings(as.numeric(d$df))
  lfc <- d$log2fc
  method <- rep(NA_character_, nrow(d))
  have <- is.finite(se) & se > 0
  method[have] <- "se"
  use_t <- !have & is.finite(t) & t != 0
  se[use_t] <- abs(lfc[use_t] / t[use_t]); method[use_t] <- "t"
  use_p <- !have & !use_t & is.finite(p) & p > 0 & p < 1 & lfc != 0
  q <- ifelse(is.finite(df) & df > 0, qt(1 - p / 2, pmax(df, 1)), qnorm(1 - p / 2))
  se[use_p] <- abs(lfc[use_p]) / q[use_p]; method[use_p] <- ifelse(is.finite(df[use_p]), "p_t", "p_norm")
  se[!(have | use_t | use_p)] <- NA
  list(se = se, method = method)
}

solve_df <- function(t, p) {
  ok <- is.finite(t) & is.finite(p) & p >= 0.001 & p <= 0.5 & t != 0
  if (sum(ok) < 50) return(NA_real_)
  i <- which(ok)
  if (length(i) > 400) i <- i[round(seq(1, length(i), length.out = 400))]
  v <- vapply(i, function(k) {
    f <- function(lg) log(2 * pt(-abs(t[k]), exp(lg))) - log(p[k])
    lo <- f(log(0.5)); hi <- f(log(1e5))
    if (!is.finite(lo) || !is.finite(hi)) return(NA_real_)
    if (hi > 0) return(Inf)       # p at least as large as the normal p: effectively normal
    if (lo < 0) return(0.5)
    exp(uniroot(f, c(log(0.5), log(1e5)))$root)
  }, numeric(1))
  median(v, na.rm = TRUE)
}

fit_group <- function(src, g, d) {
  lfc <- d$log2fc
  out <- data.frame(source = src, group = g, n_genes_B = sum(is.finite(lfc)),
                    median_abs_B = median(abs(lfc), na.rm = TRUE),
                    p90_abs_B = unname(quantile(abs(lfc), 0.9, na.rm = TRUE)),
                    n_genes_C = NA, se_method = NA, median_se = NA, sd_g = NA, pi0 = NA,
                    mean_abs_pm = NA, frac_lfsr05 = NA, loglik = NA, df_source = NA, df_lik = NA, lik = NA,
                    stringsAsFactors = FALSE)
  s <- se_of(d)
  ok <- is.finite(lfc) & is.finite(s$se) & s$se > 0
  out$n_genes_C <- sum(ok)
  if (sum(ok) < MIN_GENES) return(out)
  tab <- table(s$method[ok])
  out$se_method <- paste(names(tab), tab, sep = ":", collapse = ",")
  out$median_se <- median(s$se[ok])
  # likelihood df
  dfc <- suppressWarnings(as.numeric(d$df))[ok]
  m <- s$method[ok]
  df_used <- NA_real_
  if (any(m == "p_t") && any(is.finite(dfc))) {
    df_used <- median(dfc[m == "p_t" & is.finite(dfc)])
    out$df_source <- "table df"
  } else if (any(m == "t")) {
    df_used <- solve_df(suppressWarnings(as.numeric(d$t))[ok], suppressWarnings(as.numeric(d$p))[ok])
    out$df_source <- if (is.finite(df_used)) "solved from t and p" else "no p to solve df"
  } else out$df_source <- if (all(m == "se")) "SE given" else "p (normal)"
  out$df_lik <- df_used
  use_t <- is.finite(df_used) && df_used < DF_T
  out$lik <- if (use_t) sprintf("t(df=%.1f)", df_used) else "normal"
  a <- tryCatch(if (use_t) ash(lfc[ok], s$se[ok], lik = lik_t(df_used), mixcompdist = "uniform", pointmass = TRUE, prior = "nullbiased")
                else ash(lfc[ok], s$se[ok], mixcompdist = "normal", pointmass = TRUE, prior = "nullbiased"),
                error = function(e) NULL)
  if (is.null(a)) return(out)
  fg <- get_fitted_g(a)
  if (use_t) {
    mu <- (fg$a + fg$b) / 2; v <- (fg$b - fg$a)^2 / 12
    out$sd_g <- sqrt(sum(fg$pi * (v + mu^2)))
    out$pi0 <- sum(fg$pi[fg$a == 0 & fg$b == 0])
  } else {
    out$sd_g <- sqrt(sum(fg$pi * (fg$sd^2 + fg$mean^2)))
    out$pi0 <- sum(fg$pi[fg$sd == 0])
  }
  pm <- get_pm(a)
  out$mean_abs_pm <- mean(abs(pm))
  out$frac_lfsr05 <- mean(get_lfsr(a) < 0.05)
  out$loglik <- get_loglik(a)
  out
}

files <- sort(list.files(full_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE))
jobs <- list()
for (f in files) {
  src <- sub("\\.tsv\\.gz$", "", basename(f))
  d <- read.delim(gzfile(f), colClasses = c(group = "character", gene = "character"), quote = "")
  d$log2fc <- as.numeric(d$log2fc)
  for (g in unique(d$group)) jobs[[length(jobs) + 1]] <- list(src = src, g = g, d = d[d$group == g, ])
}
cat(length(files), "sources,", length(jobs), "groups\n")
res <- mclapply(jobs, function(j) fit_group(j$src, j$g, j$d), mc.cores = max(1, detectCores() - 2))
bad <- vapply(res, function(r) inherits(r, "try-error") || !is.data.frame(r), logical(1))
if (any(bad)) stop("fits failed: ", sum(bad))
res <- do.call(rbind, res)
out <- file.path(here, "data", "by_group.tsv")
con <- file(out, "w")
writeLines(c(sprintf("# built %s by effect_sizes/ash_fit.R (ashr %s, R %s); pointmass TRUE, prior nullbiased; MIN_GENES %d; normal likelihood and normal mixture, or t likelihood (lik_t, df < %d) and uniform mixture (column lik)",
                     Sys.Date(), packageVersion("ashr"), getRversion(), MIN_GENES, DF_T),
             paste("# inputs: data/full/*.tsv.gz from harmonize.py; data/harmonize_meta.json sha256",
                   substr(system(paste("shasum -a 256", shQuote(file.path(here, "data", "harmonize_meta.json"))), intern = TRUE), 1, 12))), con)
write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat("wrote", out, nrow(res), "rows;", sum(is.finite(res$sd_g)), "fitted\n")
