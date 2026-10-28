# S2: per-dataset limma fits. Usage: Rscript s2_fit.R <mode: primary|v1first|v1repro>   (run from ph_axis_v2/)
# Models per dataset (assigned persons only; replicates of one person averaged):
#   ctrl : expr ~ pH + age + sex   (controls)           case : same in cases      both : expr ~ pH + age + sex + sz (cases and controls)
#   sz_base : expr ~ sz + age + sex                      sz_pH : expr ~ sz + age + sex + pH
# sex is dropped from a model when the less frequent sex has fewer than 2 persons. Bottom tenth of genes by mean expression (all assigned persons of the dataset) dropped.
suppressMessages({library(limma); library(data.table); library(edgeR)})
set.seed(20261014)
args <- commandArgs(TRUE); mode <- args[1]
OUT <- file.path("out", mode); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- fread(file.path("input", sprintf("samples_v2_%s.tsv", mode)))
S <- S[assigned == TRUE]
reg <- fread(file.path("input", ifelse(mode == "v1repro", "dataset_registry_v2_primary.tsv", sprintf("dataset_registry_v2_%s.tsv", mode))))
MIN1 <- if (mode == 'v1repro') 8 else 10   # persons for a single-group model (v1repro: 8, as the stored within-group script)
MINB <- 14   # persons for the both-groups model (and at least 5 per group)
read_expr <- function(d) {
  if (d == "NSW_OFC_snRNA") {
    x <- fread(cmd = "gzip -dc nserc_pilots_2026-10-14/epigenetics_h2az/frz_nsw/pb/pseudobulk_all.tsv.gz", sep = "\t")
    don <- x[[1]]; m <- t(as.matrix(x[, -1])); colnames(m) <- don
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
    y <- DGEList(m); keep <- rowSums(cpm(y) >= 1) >= 0.2 * ncol(m); y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
    return(cpm(y, log = TRUE, prior.count = 2))
  }
  f <- if (grepl("^GSE93987|^GSE87610", d)) file.path("input", sprintf("lcm_%s.expr.tsv.gz", d)) else file.path("input", sprintf("%s.expr.tsv.gz", d))
  x <- fread(f, sep = "\t"); g <- x[[1]]; m <- as.matrix(x[, -1]); rownames(m) <- g; m[!duplicated(rownames(m)), , drop = FALSE]
}
fit1 <- function(Y, cov, terms, coef) {
  # terms: character vector of covariate columns in cov; coef: the term reported
  X <- model.matrix(as.formula(paste("~", paste(terms, collapse = "+"))), data = cov)
  fit <- lmFit(Y, X); fe <- eBayes(fit); j <- which(colnames(X) == coef)
  se <- sqrt(fe$s2.post) * fe$stdev.unscaled[, j]; ols_t <- fit$coefficients[, j] / (fit$sigma * fit$stdev.unscaled[, j])
  data.table(gene = rownames(Y), beta = fit$coefficients[, j], se = se, t = fe$t[, j], p = fe$p.value[, j], t_ols = ols_t, sd = apply(Y, 1, sd), n = nrow(cov), df_total = fe$df.total[1])
}
for (d in unique(S$dataset)) {
  s <- S[dataset == d]; E <- read_expr(d)
  sid <- if (d == "NSW_OFC_snRNA") s$sample else s$sample
  ok <- sid %in% colnames(E); if (!all(ok)) message(d, ": ", sum(!ok), " samples without expression, dropped"); s <- s[ok]
  # collapse replicates of one person
  pers <- unique(s$person); Y <- sapply(pers, function(p) rowMeans(E[, s$sample[s$person == p], drop = FALSE]))
  cov <- s[match(pers, s$person), .(person, sz, age, sex, pH)]; cov[, sexM := as.numeric(sex == "M")]
  # the 18 Gemma-processed dataset-regions arrive with the bottom tenth by mean expression already removed (v1 filter, over all samples of the dataset); LCM and pseudobulk get the same filter here
  if (grepl('^GSE93987|^GSE87610|^NSW_OFC', d)) { ym <- rowMeans(Y); Y <- Y[ym > quantile(ym, 0.10), , drop = FALSE] }
  res <- list(); nn <- c(ctrl = sum(cov$sz == 0), case = sum(cov$sz == 1))
  dir.create(file.path(OUT,'matrices'),showWarnings=FALSE); fwrite(data.table(gene=rownames(Y),Y),file.path(OUT,'matrices',sprintf('%s.expr.tsv.gz',d)),sep='\t'); fwrite(cov[,.(person,sz,age,sex,pH)],file.path(OUT,'matrices',sprintf('%s.cov.tsv',d)),sep='\t')
  mk <- function(sub, terms, coef, name, minn) {
    c2 <- cov[sub]; if (nrow(c2) < minn) return(invisible())
    tm <- terms; if ("sexM" %in% tm && min(table(factor(c2$sexM, 0:1))) < 2) tm <- setdiff(tm, "sexM")
    if (sd(c2$pH) == 0) return(invisible())
    r <- fit1(Y[, sub, drop = FALSE], c2, tm, coef); r[, model := name]; r[, dataset := d]; res[[name]] <<- r
  }
  mk(cov$sz == 0, c("pH", "age", "sexM"), "pH", "ctrl", MIN1)
  mk(cov$sz == 1, c("pH", "age", "sexM"), "pH", "case", MIN1)
  if (mode != "v1repro") {
    if (nrow(cov) >= MINB && min(nn) >= 5) {
      mk(rep(TRUE, nrow(cov)), c("pH", "age", "sexM", "sz"), "pH", "both", MINB)
      mk(rep(TRUE, nrow(cov)), c("sz", "age", "sexM"), "sz", "sz_base", MINB)
      mk(rep(TRUE, nrow(cov)), c("sz", "age", "sexM", "pH"), "sz", "sz_pH", MINB)
    }
  }
  if (length(res)) {
    R <- rbindlist(res); fwrite(R, file.path(OUT, sprintf("%s.tsv.gz", d)), sep = "\t")
    cat(sprintf("%s: persons %d (ctrl %d, case %d), genes %d, models %s\n", d, nrow(cov), nn[1], nn[2], nrow(Y), paste(names(res), collapse = ",")))
  } else cat(sprintf("%s: persons %d (ctrl %d, case %d): no model run\n", d, nrow(cov), nn[1], nn[2]))
}
