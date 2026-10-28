# SZ vs control, sex and age, in the same datasets, with one pipeline (limma).
#
# Run from controls_age_sex/:
#   Rscript scripts/de_contrasts.R full                 -> de/<dataset>.tsv, results/de_n.tsv, results/sex_check*.tsv
#   Rscript scripts/de_contrasts.R matched <workdir>    -> <workdir>/rep_<r>/<dataset>.bin (float32; read by pairs.py),
#                                                          results/matched_design.tsv
# Inputs (read only): data/datasets.tsv (scripts/build_covariates.py); Gemma processed expression and covariate files named
# there; ../../data/reference/hgnc_complete_set.txt (NCBI gene id -> HGNC symbol, as in ../confounds).
#
# Preprocessing is that of ../confounds/scripts/within_dataset_confounds.R (read_gemma, gene_level and the filters below are
# copied from it): probes annotated to one NCBI gene, mapped to HGNC symbols; per gene the probe with the highest mean;
# samples Gemma left without values (>= 50% missing) dropped; arrays of one subject averaged; genes with any missing value
# dropped; genes in the bottom tenth of mean expression dropped. Samples: SZ and controls with age and sex recorded.
#
# Full-data contrasts (limma lmFit + eBayes; one model per contrast; same covariates except the variable of interest):
#   sz       SZ vs control, all subjects:        ~ sz + age + sex
#   sex      male vs female, all subjects:        ~ male + sz + age
#   sex_ctl  male vs female, controls only:       ~ male + age
#   age      age in years, all subjects:          ~ age + sz + sex
#   age_ctl  age in years, controls only:         ~ age + sex
#   age2     older vs younger half, all subjects:  ~ old + sz + sex   (median split; ties at the median broken by sample order)
# A covariate with one level in the samples at hand is left out of that model.
#
# Matched subsampling (R = 50 repeats per dataset): every contrast uses the same number of samples and the same two group
# sizes (a, b) in a dataset. a = min(smaller SZ/control group, smaller sex group, younger half), b = min(larger SZ/control
# group, larger sex group, older half). Per repeat: sz_m draws a from the smaller diagnosis group (SZ on a tie) and b from the
# other; sex_m draws a from the smaller sex and b from the other; age2_m draws a from the younger and b from the older half;
# age_m draws a + b subjects at random and fits age as a continuous variable. Same models as the full data. Datasets with
# a < 4 are left out of the matched analysis.
suppressMessages({library(limma); library(data.table); library(parallel)})
args <- commandArgs(TRUE)
MODE <- if (length(args)) args[1] else "full"
WORK <- if (length(args) > 1) args[2] else NA
R_REPS <- 50
NCORES <- 6
SEED <- 20260925
MARKERS <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")

hg <- fread("../../data/reference/hgnc_complete_set.txt", select = c("symbol", "entrez_id"), quote = "")
hg <- hg[!is.na(entrez_id)]
ent2sym <- setNames(hg$symbol, as.character(hg$entrez_id))
REG <- fread("data/datasets.tsv", skip = "dataset\t", colClasses = list(character = c("region", "cluster")))

read_gemma <- function(path, samples_json) {  # from within_dataset_confounds.R
  x <- fread(cmd = paste("gzip -dc", shQuote(path), "| grep -v '^#'"), sep = "\t", quote = "", colClasses = list(character = c("Probe", "NCBIid")))
  if (all(is.na(x$NCBIid) | x$NCBIid == "")) {
    # gemma2's processed file for some array datasets carries no gene columns; map probes with the platform's Gemma
    # annotation file (scripts/fetch_gemma.py, data/gemma/platform_<GPL>_annotations.tsv.gz)
    gpl <- unique(vapply(rjson::fromJSON(file = samples_json)$data, function(s) s$arrayDesign$shortName, ""))
    an <- rbindlist(lapply(gpl, function(g) fread(cmd = paste("gzip -dc", shQuote(sprintf("data/gemma/platform_%s_annotations.tsv.gz", g)), "| grep -v '^#'"),
                                                    sep = "\t", quote = "", select = c("ElementName", "NCBIids"), colClasses = "character")))
    x[, NCBIid := an$NCBIids[match(Probe, an$ElementName)]]
    cat(basename(path), ": NCBI ids from platform annotation", paste(gpl, collapse = ","), "for", sum(!is.na(x$NCBIid) & x$NCBIid != ""), "probes\n")
  }
  x <- x[!is.na(NCBIid) & !grepl("\\|", as.character(NCBIid)) & as.character(NCBIid) != ""]
  sym <- ent2sym[as.character(x$NCBIid)]
  keep <- !is.na(sym)
  sc <- colnames(x)[-(1:6)]
  m <- as.matrix(x[keep, ..sc]); rownames(m) <- sym[keep]
  attr(m, "probe") <- x$Probe[keep]
  if (all(grepl("GSM[0-9]+", sc))) {
    colnames(m) <- regmatches(sc, regexpr("GSM[0-9]+", sc))
  } else {
    sj <- rjson::fromJSON(file = samples_json)$data
    acc <- setNames(sapply(sj, function(s) s$accession$accession), sapply(sj, function(s) as.character(s$id)))
    colnames(m) <- acc[sub(".*BioAssayId=([0-9]+).*", "\\1", sc)]
  }
  m
}
gene_level <- function(m) {  # from within_dataset_confounds.R
  mu <- rowMeans(m, na.rm = TRUE)
  o <- order(rownames(m), -mu)
  m <- m[o, , drop = FALSE]
  m[!duplicated(rownames(m)), , drop = FALSE]
}

load_dataset <- function(D) {
  m <- read_gemma(D$expr, D$samples_json)
  cv <- fread(D$cov, skip = "gsm\t", colClasses = list(character = intersect(c("gsm", "subject", "sex", "dx", "region", "cohort", "bank", "gemma_outlier"),
                                                                             names(fread(D$cov, skip = "gsm\t", nrows = 1)))))
  if (!"subject" %in% names(cv)) cv[, subject := ""]
  cv[is.na(subject), subject := ""]
  if (nzchar(D$region)) cv <- cv[region == D$region]
  if ("gemma_outlier" %in% names(cv)) cv <- cv[is.na(gemma_outlier) | gemma_outlier != "yes"]
  cv <- cv[gsm %in% colnames(m) & dx %in% c("control", "schizophrenia") & !is.na(age) & sex %in% c("M", "F")]
  ismk <- rownames(m) %in% MARKERS
  Mk <- m[ismk, cv$gsm, drop = FALSE]; rownames(Mk) <- paste(rownames(m)[ismk], attr(m, "probe")[ismk])  # probe level, sex-label check
  E <- m[, cv$gsm, drop = FALSE]
  okS <- colMeans(is.na(E)) < 0.5
  E <- E[, okS, drop = FALSE]; Mk <- Mk[, okS, drop = FALSE]; cv <- cv[okS]
  if (any(cv$subject != "") && any(duplicated(cv$subject[cv$subject != ""]))) {
    avg <- function(X) { X <- t(apply(X, 1, function(v) tapply(v, cv$subject, mean, na.rm = TRUE))); X }
    E <- avg(E); if (nrow(Mk)) Mk <- avg(Mk)
    cv <- cv[!duplicated(subject)]
    E <- E[, cv$subject, drop = FALSE]; colnames(E) <- cv$gsm
    if (nrow(Mk)) { Mk <- Mk[, cv$subject, drop = FALSE]; colnames(Mk) <- cv$gsm }
  }
  E <- gene_level(E)
  E <- E[rowSums(is.na(E)) == 0, , drop = FALSE]
  Eall <- E
  mu <- rowMeans(E)
  E <- E[mu > quantile(mu, 0.10), , drop = FALSE]
  cv[, sz := as.integer(dx == "schizophrenia")]
  cv[, male := as.integer(sex == "M")]
  o <- order(cv$age, seq_len(nrow(cv)))
  cv[, old := 0L]; cv$old[o[(floor(nrow(cv) / 2) + 1):nrow(cv)]] <- 1L
  list(E = E, Eall = Eall, Mk = Mk[rowSums(is.na(Mk)) == 0, , drop = FALSE], cv = cv)
}

fit <- function(E, d, coef, covs) {
  covs <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
  des <- model.matrix(as.formula(paste("~", paste(c(coef, covs), collapse = " + "))), data = d)
  f <- eBayes(lmFit(E[, d$gsm, drop = FALSE], des))
  tt <- topTable(f, coef = coef, number = Inf, sort.by = "none")
  list(logfc = tt$logFC, t = tt$t, p = tt$P.Value, covs = covs)
}
MODELS <- list(sz = c("sz", "age", "male"), sex = c("male", "sz", "age"), sex_ctl = c("male", "age"),
               age = c("age", "sz", "male"), age_ctl = c("age", "male"), age2 = c("old", "sz", "male"))

# Sex-label check, label-free, at probe level (the gene-level value keeps the probe with the highest mean, which for XIST
# on GPL570 is 243712_at, a probe with no sex difference): each marker probe (XIST, RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY)
# is split into two clusters (kmeans, k = 2); per marker gene the probe with the widest gap between cluster means is kept,
# and the marker is informative when that gap is at least 1.5 log2 units. Each informative marker votes: XIST high = F,
# a Y gene high = M. A sample is flagged when the majority of informative markers disagrees with its label.
sex_check <- function(Mk, cv, name) {
  genes <- sub(" .*", "", rownames(Mk))
  votes <- list(); info <- character(0)
  for (g in intersect(MARKERS, genes)) {
    best <- NULL
    for (pr in rownames(Mk)[genes == g]) {
      v <- Mk[pr, cv$gsm]
      if (length(unique(v)) < 2) next
      km <- kmeans(v, centers = c(min(v), max(v)))
      gap <- diff(range(km$centers))
      if (is.null(best) || gap > best$gap) best <- list(gap = gap, hi = km$cluster == which.max(km$centers), probe = sub("^\\S+ ", "", pr))
    }
    if (is.null(best) || best$gap < 1.5) next
    votes[[g]] <- if (g == "XIST") ifelse(best$hi, "F", "M") else ifelse(best$hi, "M", "F")
    info <- c(info, paste0(g, "(", best$probe, ")"))
  }
  meas <- paste(intersect(MARKERS, genes), collapse = ",")
  if (!length(votes)) return(data.table(dataset = name, gsm = cv$gsm, dx = cv$dx, sex_label = cv$sex, n_informative = 0L,
                                        votes_M = NA_integer_, sex_predicted = NA_character_, mismatch = FALSE,
                                        markers_measured = meas, markers_informative = ""))
  V <- do.call(cbind, votes)
  nM <- rowSums(V == "M")
  pred <- ifelse(nM > ncol(V) / 2, "M", ifelse(nM < ncol(V) / 2, "F", NA))
  data.table(dataset = name, gsm = cv$gsm, dx = cv$dx, sex_label = cv$sex, n_informative = ncol(V), votes_M = nM,
             sex_predicted = pred, mismatch = !is.na(pred) & pred != cv$sex,
             markers_measured = meas, markers_informative = paste(info, collapse = ","))
}

if (MODE == "full") {
  dir.create("de", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
  nrows <- list(); sc <- list()
  for (i in seq_len(nrow(REG))) {
    D <- as.list(REG[i]); L <- load_dataset(D); E <- L$E; cv <- L$cv
    out <- data.table(gene = rownames(E), mean_expr = round(rowMeans(E), 4))
    for (cn in names(MODELS)) {
      d <- if (grepl("_ctl$", cn)) cv[dx == "control"] else cv
      v <- MODELS[[cn]]
      if (nrow(d) < 8 || length(unique(d[[v[1]]])) < 2) next
      r <- fit(E, d, v[1], v[-1])
      out[, (paste0("logfc_", cn)) := signif(r$logfc, 5)]
      out[, (paste0("t_", cn)) := signif(r$t, 5)]
      out[, (paste0("p_", cn)) := signif(r$p, 5)]
      out[, (paste0("q_", cn)) := signif(p.adjust(r$p, "BH"), 5)]
      g1 <- d[[v[1]]]
      nrows[[length(nrows) + 1]] <- data.table(dataset = D$dataset, contrast = cn, n = nrow(d),
        n_group1 = if (v[1] == "age") NA_integer_ else sum(g1 == 1), n_group0 = if (v[1] == "age") NA_integer_ else sum(g1 == 0),
        group1 = c(sz = "SZ", male = "male", old = "older half", age = "per year")[[v[1]]],
        n_sz = sum(d$sz), n_ctl = sum(1 - d$sz), n_male = sum(d$male), n_female = sum(1 - d$male),
        age_min = min(d$age), age_median = median(d$age), age_max = max(d$age),
        model = paste("~", paste(c(v[1], r$covs), collapse = " + ")), n_genes = nrow(E),
        n_q05 = sum(p.adjust(r$p, "BH") < 0.05))
    }
    fwrite(out, sprintf("de/%s.tsv", D$dataset), sep = "\t")
    s <- sex_check(L$Mk, cv, D$dataset); if (!is.null(s)) sc[[length(sc) + 1]] <- s
    cat(D$dataset, nrow(E), "genes,", nrow(cv), "samples\n")
  }
  fwrite(rbindlist(nrows), "results/de_n.tsv", sep = "\t")
  S <- rbindlist(sc)
  fwrite(S, "results/sex_check_samples.tsv", sep = "\t")
  fwrite(S[, .(n = .N, n_label_F = sum(sex_label == "F"), markers_measured = markers_measured[1], markers_informative = markers_informative[1], n_mismatch = sum(mismatch), n_undecided = sum(is.na(sex_predicted)),
               mismatched = paste(gsm[mismatch], sex_label[mismatch], sep = ":labelled ", collapse = "; ")), by = dataset],
         "results/sex_check.tsv", sep = "\t")
} else if (MODE == "matched") {
  stopifnot(!is.na(WORK))
  design <- mclapply(seq_len(nrow(REG)), function(i) {
    D <- as.list(REG[i]); L <- load_dataset(D); E <- L$E; cv <- L$cv
    nsz <- sum(cv$sz); nct <- sum(1 - cv$sz); nm <- sum(cv$male); nf <- sum(1 - cv$male); ny <- sum(1 - cv$old); no <- sum(cv$old)
    a <- min(min(nsz, nct), min(nm, nf), ny); b <- min(max(nsz, nct), max(nm, nf), no)
    row <- data.table(dataset = D$dataset, n_sz = nsz, n_ctl = nct, n_male = nm, n_female = nf, n_younger = ny, n_older = no,
                      a = a, b = b, n_matched = a + b, used = a >= 4, reps = if (a >= 4) R_REPS else 0L)
    if (a < 4) return(row)
    small_dx <- if (nsz <= nct) 1L else 0L
    small_sex <- if (nf <= nm) 0L else 1L
    draw <- function(v, small) { i1 <- which(v == small); i2 <- which(v != small); c(i1[sample.int(length(i1), a)], i2[sample.int(length(i2), b)]) }
    for (r in seq_len(R_REPS)) {
      set.seed(SEED + 1000 * r + i)
      dir.create(file.path(WORK, sprintf("rep_%02d", r)), showWarnings = FALSE, recursive = TRUE)
      res <- list(
        sz_m = fit(E, cv[draw(cv$sz, small_dx)], "sz", c("age", "male")),
        sex_m = fit(E, cv[draw(cv$male, small_sex)], "male", c("sz", "age")),
        age2_m = fit(E, cv[draw(cv$old, 0L)], "old", c("sz", "male")),
        age_m = fit(E, cv[sample.int(nrow(cv), a + b)], "age", c("sz", "male")))
      con <- file(file.path(WORK, sprintf("rep_%02d", r), paste0(D$dataset, ".bin")), "wb")
      for (x in res) { writeBin(as.numeric(x$t), con, size = 4); writeBin(as.numeric(x$p), con, size = 4) }
      close(con)
    }
    writeLines(rownames(E), file.path(WORK, paste0(D$dataset, "_genes.txt")))
    row
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  errs <- vapply(design, function(x) inherits(x, "try-error"), logical(1))
  if (any(errs)) { print(design[errs]); stop("matched fits failed") }
  fwrite(rbindlist(design), "results/matched_design.tsv", sep = "\t")
}
