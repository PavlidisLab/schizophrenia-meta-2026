# Shared loaders and agreement metrics for fairness_sample_level/ (sourced by the a_, b_, c_, d_ scripts).
# Run the scripts from fairness_sample_level/.
#
# Datasets: the registry ../controls_age_sex/data/datasets.tsv (read only; built by the controls_age_sex analysis,
# scripts/build_covariates.py there). Its expression, sample and covariate paths are relative to ../controls_age_sex/.
# read_gemma(), gene_level() and load_dataset() are copied from ../controls_age_sex/scripts/de_contrasts.R, which copied
# them from ../confounds/scripts/within_dataset_confounds.R: probes annotated to one NCBI gene mapped to HGNC symbols; per
# gene the probe with the highest mean; samples Gemma left without values (>= 50% missing) dropped; arrays of one subject
# averaged; genes with any missing value dropped; genes in the bottom tenth of mean expression dropped; SZ and controls
# with age and sex recorded. Additions here: every covariate column is kept; Gemma's batch factor (category 'block') and
# GSE138082's antipsychotic on/off field are attached per GSM from data/gemma_batch/<id>_samples.json
# (scripts/fetch_gemma_batch.py).
suppressMessages({library(limma); library(data.table); library(parallel)})

CAS <- "../controls_age_sex"
hg <- fread("../../data/reference/hgnc_complete_set.txt", select = c("symbol", "entrez_id"), quote = "")
hg <- hg[!is.na(entrez_id)]
ent2sym <- setNames(hg$symbol, as.character(hg$entrez_id))
REG <- fread(file.path(CAS, "data/datasets.tsv"), skip = "dataset\t", colClasses = list(character = c("region", "cluster", "gemma_id")))
cas_path <- function(p) if (startsWith(p, "/")) p else file.path(CAS, p)

sha256 <- function(f) sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE))

read_gemma <- function(path, samples_json) {
  x <- fread(cmd = paste("gzip -dc", shQuote(path), "| grep -v '^#'"), sep = "\t", quote = "", colClasses = list(character = c("Probe", "NCBIid")))
  if (all(is.na(x$NCBIid) | x$NCBIid == "")) {
    gpl <- unique(vapply(rjson::fromJSON(file = samples_json)$data, function(s) s$arrayDesign$shortName, ""))
    an <- rbindlist(lapply(gpl, function(g) fread(cmd = paste("gzip -dc", shQuote(file.path(CAS, sprintf("data/gemma/platform_%s_annotations.tsv.gz", g))), "| grep -v '^#'"),
                                                    sep = "\t", quote = "", select = c("ElementName", "NCBIids"), colClasses = "character")))
    x[, NCBIid := an$NCBIids[match(Probe, an$ElementName)]]
  }
  x <- x[!is.na(NCBIid) & !grepl("\\|", as.character(NCBIid)) & as.character(NCBIid) != ""]
  sym <- ent2sym[as.character(x$NCBIid)]
  keep <- !is.na(sym)
  sc <- colnames(x)[-(1:6)]
  m <- as.matrix(x[keep, ..sc]); rownames(m) <- sym[keep]
  if (all(grepl("GSM[0-9]+", sc))) {
    colnames(m) <- regmatches(sc, regexpr("GSM[0-9]+", sc))
  } else {
    sj <- rjson::fromJSON(file = samples_json)$data
    acc <- setNames(sapply(sj, function(s) s$accession$accession), sapply(sj, function(s) as.character(s$id)))
    colnames(m) <- acc[sub(".*BioAssayId=([0-9]+).*", "\\1", sc)]
  }
  m
}
gene_level <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  o <- order(rownames(m), -mu)
  m <- m[o, , drop = FALSE]
  m[!duplicated(rownames(m)), , drop = FALSE]
}

# Gemma batch ('block' factor value) and, for GSE138082, the 'drug' characteristic, per GSM
gemma_sample_fields <- function(gid) {
  f <- sprintf("data/gemma_batch/%s_samples.json", gid)
  sj <- jsonlite::fromJSON(f, simplifyVector = FALSE)$data
  rbindlist(lapply(sj, function(s) {
    fv <- s$sample$factorValues
    blk <- NA_character_
    for (v in fv) if (identical(v$experimentalFactorCategory$category, "block")) blk <- if (!is.null(v$value)) v$value else v$summary
    drug <- NA_character_
    for (ch in s$sample$characteristics) if (identical(ch$category, "drug")) drug <- ch$value
    data.table(gsm = s$accession$accession, batch = blk, drug = drug)
  }))
}

load_dataset <- function(D, bottom_filter = TRUE) {
  m <- read_gemma(cas_path(D$expr), cas_path(D$samples_json))
  cvf <- cas_path(D$cov)
  hdr <- names(fread(cvf, skip = "gsm\t", nrows = 1))
  cv <- fread(cvf, skip = "gsm\t", colClasses = list(character = intersect(c("gsm", "subject", "sex", "dx", "region", "cohort", "bank", "gemma_outlier"), hdr)))
  if (!"subject" %in% names(cv)) cv[, subject := ""]
  cv[is.na(subject), subject := ""]
  for (v in c("pH", "PMI", "RIN", "agonal", "cpz_mg")) if (!v %in% names(cv)) cv[, (v) := NA_real_]
  if (nzchar(D$region)) cv <- cv[region == D$region]
  if ("gemma_outlier" %in% names(cv)) cv <- cv[is.na(gemma_outlier) | gemma_outlier != "yes"]
  cv <- cv[gsm %in% colnames(m) & dx %in% c("control", "schizophrenia") & !is.na(age) & sex %in% c("M", "F")]
  gf <- gemma_sample_fields(D$gemma_id)
  cv <- merge(cv, gf, by = "gsm", all.x = TRUE, sort = FALSE)
  E <- m[, cv$gsm, drop = FALSE]
  okS <- colMeans(is.na(E)) < 0.5
  E <- E[, okS, drop = FALSE]; cv <- cv[okS]
  n_arrays <- nrow(cv)
  if (any(cv$subject != "") && any(duplicated(cv$subject[cv$subject != ""]))) {
    E <- t(apply(E, 1, function(v) tapply(v, cv$subject, mean, na.rm = TRUE)))
    # batch of an averaged subject: kept only when all its arrays share one batch
    bb <- cv[, .(b = if (uniqueN(batch) == 1) batch[1] else NA_character_), by = subject]
    cv <- cv[!duplicated(subject)]
    cv[, batch := bb$b[match(subject, bb$subject)]]
    E <- E[, cv$subject, drop = FALSE]; colnames(E) <- cv$gsm
  }
  E <- gene_level(E)
  E <- E[rowSums(is.na(E)) == 0, , drop = FALSE]
  if (bottom_filter) {
    mu <- rowMeans(E)
    E <- E[mu > quantile(mu, 0.10), , drop = FALSE]
  }
  cv[, sz := as.integer(dx == "schizophrenia")]
  cv[, male := as.integer(sex == "M")]
  cv[, on_drug := ifelse(!is.na(drug), as.integer(drug == "ON"), NA_integer_)]
  attr(E, "n_arrays") <- n_arrays
  list(E = E, cv = cv, name = D$dataset, compartment = D$compartment, cluster = D$cluster, gse = D$gse, gemma_id = D$gemma_id)
}

load_all <- function(which = NULL) {
  R <- if (is.null(which)) REG else REG[dataset %in% which]
  L <- lapply(seq_len(nrow(R)), function(i) load_dataset(as.list(R[i])))
  names(L) <- R$dataset
  L
}

# limma SZ vs control; covs dropped when they have one level in the samples at hand
fit_sz <- function(E, d, covs = c("age", "male"), extra_design = NULL) {
  covs <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
  des <- model.matrix(as.formula(paste("~", paste(c("sz", covs), collapse = " + "))), data = d)
  if (!is.null(extra_design)) des <- cbind(des, extra_design)
  if (qr(des)$rank < ncol(des)) return(NULL)
  f <- eBayes(lmFit(E[, d$gsm, drop = FALSE], des))
  j <- which(colnames(des) == "sz")
  list(gene = rownames(E), t = f$t[, j], logfc = f$coefficients[, j], p = f$p.value[, j],
       se = f$stdev.unscaled[, j] * sqrt(f$s2.post), sigma = f$sigma, df = f$df.total[1], covs = covs, n_sz = sum(d$sz), n_ctl = sum(1 - d$sz))
}

pi0_safe <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) < 10) return(NA_real_)
  r <- tryCatch(qvalue::pi0est(p)$pi0, error = function(e) NA_real_)
  if (is.na(r)) r <- min(1, mean(p > 0.5) / 0.5)
  r
}

TOPN <- c(20, 50, 100, 200)
PI1_NULL_DRAWS <- 20
PI1_NULL_SEED <- 20260926
# agreement between two SZ results a, b (lists from fit_sz); metrics as in ../h1.py where they exist
agree <- function(a, b, pi1 = TRUE) {
  g <- intersect(a$gene, b$gene)
  ia <- match(g, a$gene); ib <- match(g, b$gene)
  ta <- a$t[ia]; tb <- b$t[ib]; la <- a$logfc[ia]; lb <- b$logfc[ib]; pa <- a$p[ia]; pb <- b$p[ib]
  U <- length(g)
  oa <- order(pa, -abs(ta)); ob <- order(pb, -abs(tb))
  out <- list(n_common = U, rho_t = cor(ta, tb, method = "spearman"), rho_logfc = cor(la, lb, method = "spearman"),
              r_logfc = cor(la, lb))
  for (N in TOPN) {
    sa <- oa[1:N]; sb <- ob[1:N]
    ov <- length(intersect(sa, sb))
    out[[paste0("top", N, "_overlap")]] <- ov
    out[[paste0("top", N, "_fold")]] <- ov / (N^2 / U)
  }
  sa <- oa[1:100]; sb <- ob[1:100]; sh <- intersect(sa, sb)
  out$top100_n_dir <- length(sh)
  out$top100_concordance <- if (length(sh)) mean(sign(ta[sh]) == sign(tb[sh])) else NA_real_
  out$top100_sign_in_other <- mean(c(mean(sign(ta[sa]) == sign(tb[sa])), mean(sign(ta[sb]) == sign(tb[sb]))))
  qa <- p.adjust(pa, "BH"); qb <- p.adjust(pb, "BH")
  na <- sum(qa < 0.05); nb <- sum(qb < 0.05); ov <- sum(qa < 0.05 & qb < 0.05)
  out$fdr_n_a <- na; out$fdr_n_b <- nb; out$fdr_overlap <- ov
  out$fdr_fold <- if (na > 0 && nb > 0) ov / (na * nb / U) else NA_real_
  out$fdr_concordance <- if (ov > 0) mean(sign(ta[qa < 0.05 & qb < 0.05]) == sign(tb[qa < 0.05 & qb < 0.05])) else NA_real_
  if (pi1) {
    out$pi1_ab <- 1 - pi0_safe(pb[sa]); out$pi1_ba <- 1 - pi0_safe(pa[sb])
    out$pi1 <- mean(c(out$pi1_ab, out$pi1_ba))
    out$pi1_base <- mean(c(1 - pi0_safe(pb), 1 - pi0_safe(pa)))
    # null reference: pi1 of the other result's p at 100 random genes (PI1_NULL_DRAWS draws each way). pi1 from 100 p
    # values is floored at 0, so with no shared signal it averages about 0.11, not 0 (tests/sample_level/r/
    # test_agree_splithalf.R); read pi1 against pi1_null, or use pi1_excess = pi1 - pi1_null
    # fixed seed, and the caller's RNG state restored: pi1_null is the same on every run, and the split-half and
    # bootstrap draws of the calling scripts are the ones they were before pi1_null existed
    had <- exists(".Random.seed", envir = .GlobalEnv); if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(PI1_NULL_SEED)
    rnd <- function(p) mean(vapply(seq_len(PI1_NULL_DRAWS), function(i) 1 - pi0_safe(p[sample.int(U, 100)]), numeric(1)))
    out$pi1_null <- mean(c(rnd(pb), rnd(pa)))
    if (had) assign(".Random.seed", saved, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv)
    out$pi1_excess <- out$pi1 - out$pi1_null
  }
  as.data.table(out)
}

# dataset pair relations: the project's canonical table ../subject_relations/subject_relations.tsv (built by
# ../subject_relations/subject_relations.py; tracker, Vocabulary). Column `relation` here = its `relation_group` (one reporting
# group per relation, never pooled; 'subject-separate' is the primary comparison set).
SUBJECT_RELATIONS <- "../subject_relations/subject_relations.tsv"
read_pairs <- function() {
  P <- fread(SUBJECT_RELATIONS, skip = "dataset_a\t")
  P <- P[dataset_a %in% REG$dataset & dataset_b %in% REG$dataset]
  P[, relation := relation_group]
  P[, both_bulk := compartment_a == "bulk" & compartment_b == "bulk"]
  P
}
