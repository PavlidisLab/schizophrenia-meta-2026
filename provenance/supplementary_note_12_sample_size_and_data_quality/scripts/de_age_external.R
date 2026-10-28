# Age (and sex, as a label check) contrasts in the external datasets, with the functions of the SZ-study analysis.
#
# Run from controls_age_external/:
#   Rscript scripts/de_age_external.R external       -> de/<dataset>.tsv, results/de_n.tsv, results/sex_check*.tsv
#
# ../controls_age_sex/scripts/de_contrasts.R is sourced unchanged: with the argument 'external' it only defines its
# functions (read_gemma, gene_level, load_dataset, fit, sex_check) and reads data/datasets.tsv here as its registry.
# Preprocessing is therefore identical to the SZ-study datasets: probes annotated to one NCBI gene, HGNC symbols, highest-
# mean probe per gene, samples Gemma left without values dropped, arrays of one subject averaged, genes with any missing
# value dropped, bottom tenth of mean expression dropped.
#
# Samples: the covariate files (scripts/build_covariates.py) mark the included subjects dx = 'control': non-disease,
# age >= 18, not a Gemma outlier. load_dataset keeps dx control or schizophrenia, so only those enter.
#
# Contrasts (limma lmFit + eBayes; de_contrasts.R fit(), which drops a covariate with one level):
#   age       age in years:                 ~ age + sex      (the SZ-study 'age' model ~ age + sz + sex, with the sz term
#                                                             dropped: no SZ subjects; identical to its 'age_ctl' model)
#   age_le80  age in years, subjects aged 18 to 80 (the SZ-study bulk datasets reach 80 or less in 17 of 22):  ~ age + sex
#   sex       male vs female:               ~ male + age     (label check; X/Y genes, results/sex_sanity.tsv)
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1, args[1] == "external")
source("../controls_age_sex/scripts/de_contrasts.R")

EXT_MODELS <- list(age = c("age", "male"), age_le80 = c("age", "male"), sex = c("male", "age"))
dir.create("de", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
nrows <- list(); sc <- list()
for (i in seq_len(nrow(REG))) {
  D <- as.list(REG[i]); L <- load_dataset(D); E <- L$E; cv <- L$cv
  out <- data.table(gene = rownames(E), mean_expr = round(rowMeans(E), 4))
  for (cn in names(EXT_MODELS)) {
    d <- if (cn == "age_le80") cv[age <= 80] else cv
    v <- EXT_MODELS[[cn]]
    if (nrow(d) < 8 || length(unique(d[[v[1]]])) < 2) next
    r <- fit(E, d, v[1], v[-1])
    out[, (paste0("logfc_", cn)) := signif(r$logfc, 5)]
    out[, (paste0("t_", cn)) := signif(r$t, 5)]
    out[, (paste0("p_", cn)) := signif(r$p, 5)]
    out[, (paste0("q_", cn)) := signif(p.adjust(r$p, "BH"), 5)]
    nrows[[length(nrows) + 1]] <- data.table(dataset = D$dataset, contrast = cn, n = nrow(d),
      n_male = sum(d$male), n_female = sum(1 - d$male), age_min = min(d$age), age_median = median(d$age), age_max = max(d$age),
      model = paste("~", paste(c(v[1], r$covs), collapse = " + ")), n_genes = nrow(E), n_q05 = sum(p.adjust(r$p, "BH") < 0.05))
  }
  fwrite(out, sprintf("de/%s.tsv", D$dataset), sep = "\t")
  s <- sex_check(L$Mk, cv, D$dataset); if (!is.null(s)) sc[[length(sc) + 1]] <- s
  cat(D$dataset, nrow(E), "genes,", nrow(cv), "subjects\n")
}
fwrite(rbindlist(nrows), "results/de_n.tsv", sep = "\t")
S <- rbindlist(sc)
fwrite(S, "results/sex_check_samples.tsv", sep = "\t")
fwrite(S[, .(n = .N, n_label_F = sum(sex_label == "F"), markers_measured = markers_measured[1], markers_informative = markers_informative[1],
             n_mismatch = sum(mismatch), n_undecided = sum(is.na(sex_predicted)),
             mismatched = paste(gsm[mismatch], sex_label[mismatch], sep = ":labelled ", collapse = "; ")), by = dataset],
       "results/sex_check.tsv", sep = "\t")
