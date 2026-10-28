# Combine the per-dataset results over datasets with non-overlapping subjects.
# Run from confounds/ after scripts/within_dataset_confounds.R:  Rscript scripts/combine_datasets.R
#
# Independent set (one entry per brain collection): GSE21138 (NSW, BA46), GSE35977 (Stanley, parietal),
# GSE53987_BA46 (Pittsburgh), GSE17612 (BA10), GSE92538_GPL96 (Pritzker, the platform with more controls).
# Left out as the same subjects: GSE35974 (the GSE35977 subjects, cerebellum), GSE53987 hippocampus and striatum,
# GSE21935 (the GSE17612 subjects, BA22, per its GEO summary), GSE92538_GPL570 and GSE80655 (Pritzker subjects).
# GSE208338 shares donors with GSE21138 (same diagnosis, sex, pH, PMI and age within 1 year; results/shared_donors_*),
# so the primary combination leaves it out (analysis review C2). Sensitivities: + GSE208338 without the shared donors;
# + all of GSE208338 (the earlier six-dataset combination).
#
# Outputs:
#   results/combined_controls_pH.tsv   per gene, factor pH (and PMI): Stouffer Z over the control-only t statistics
#   results/combined_summary.tsv       set-level combined tests (sample-permutation z, SZ effect before / after pH)
#   results/recurrent_genes_pH.tsv     per recurrent gene: combined pH Z in controls and SZ Z before / after pH
suppressMessages(library(data.table))
IND <- c("GSE21138", "GSE35977", "GSE53987_BA46", "GSE17612", "GSE92538_GPL96")
VARIANTS <- list(primary_without_GSE208338 = IND, plus_GSE208338_shared_donors_dropped = c(IND, "GSE208338_noshared"),
                 plus_GSE208338_all = c(IND, "GSE208338"))
# CONFOUNDS_SETS / CONFOUNDS_OUT as in within_dataset_confounds.R (unset: recurrent_sets.tsv, results/).
# CONFOUNDS_EXTRA_VARIANTS: more dataset combinations for the set-level tests, "name=ds1,ds2;name2=ds3" (unset: none).
SETS_FILE <- Sys.getenv("CONFOUNDS_SETS", "recurrent_sets.tsv")
OUT <- Sys.getenv("CONFOUNDS_OUT", "results")
for (x in strsplit(Sys.getenv("CONFOUNDS_EXTRA_VARIANTS", ""), ";")[[1]]) if (nzchar(x)) {
  kv <- strsplit(x, "=", fixed = TRUE)[[1]]
  VARIANTS[[kv[1]]] <- strsplit(kv[2], ",", fixed = TRUE)[[1]]
}

zt <- function(t, p) sign(t) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
stouffer <- function(z, w) sum(w * z) / sqrt(sum(w^2))

# per-gene combined control associations
wc <- rbindlist(lapply(IND, function(d) { x <- fread(file.path(OUT, sprintf("within_controls_%s.tsv", d))); x[, dataset := d] }))
cp <- wc[factor %in% c("pH", "PMI"), .(n_datasets = .N, Z = stouffer(zt(t, p), sqrt(n)), datasets = paste(dataset, collapse = ",")), by = .(factor, gene)]
ep <- wc[factor == "pH", .(expr_percentile = mean(expr_percentile)), by = gene]
cp <- merge(cp, ep, by = "gene", all.x = TRUE)
cp <- cp[n_datasets >= 4]  # of 5
cp[, p := 2 * pnorm(-abs(Z))]
cp[, q := p.adjust(p, "BH"), by = factor]
fh <- file(file.path(OUT, "combined_controls_pH.tsv"), "w")
writeLines(c(paste("# scripts/combine_datasets.R run", Sys.Date(), "; input", SETS_FILE, "sha256", substr(system(paste("shasum -a 256", shQuote(SETS_FILE)), intern = TRUE), 1, 12)),
             "# Stouffer Z (weights sqrt(n controls)) over control-only limma t for the factor (expr ~ factor + age + sex),",
             paste0("# datasets: ", paste(IND, collapse = ", "), "; genes measured in >= 4 of them; q = BH within factor"),
             "# Z > 0: expression higher with higher factor value"), fh)
close(fh)
fwrite(cp, file.path(OUT, "combined_controls_pH.tsv"), sep = "\t", append = TRUE, col.names = TRUE)

rs <- fread(SETS_FILE, skip = "set\t")
rs[, sz_sign := ifelse(sz_direction == "up", 1, ifelse(sz_direction == "down", -1, 0))]
S <- fread(file.path(OUT, "within_dataset_summary.tsv"), skip = "dataset\t")

out <- list()
for (v in names(VARIANTS)) {
  ds <- VARIANTS[[v]]
  for (st in unique(rs$set)) {
    for (fac in c("pH", "PMI")) {
      x <- S[dataset %in% ds & set == st & factor == fac]
      out[[length(out) + 1]] <- data.table(variant = v, set = st, test = paste0("control ", fac, " association, signed (sample permutation)"),
        n_datasets = nrow(x), datasets_positive = sum(x$score_ols > 0),
        combined_z = round(sum(x$sample_perm_z) / sqrt(nrow(x)), 3),
        combined_p_one_sided = signif(pnorm(sum(x$sample_perm_z) / sqrt(nrow(x)), lower.tail = FALSE), 3),
        note = paste("z > 0: SZ-up genes higher at lower factor values, SZ-down genes lower. Per-dataset z = (observed - null mean) / null SD from shuffling the",
                     "factor among controls; combined z = unweighted Stouffer of those z; p = normal approximation, not a combined permutation p"))
    }
    # competitive count: recurrent genes at p < 0.05 for pH vs genome-wide rate, summed over datasets
    x <- S[dataset %in% ds & set == st & factor == "pH"]
    out[[length(out) + 1]] <- data.table(variant = v, set = st, test = "control pH, recurrent genes at p < 0.05 (summed over datasets)",
      n_datasets = nrow(x), observed = sum(x$n_set_p05), expected = round(sum(x$expected_p05), 2),
      n_opposite_sign_to_sz = sum(x$n_p05_opposite_sign_to_sz), n_same_sign_as_sz = sum(x$n_p05_same_sign_as_sz),
      combined_p_one_sided = signif(ppois(sum(x$n_set_p05) - 1, sum(x$expected_p05), lower.tail = FALSE), 3),
      note = "Poisson upper tail on the summed expected count; ignores gene-gene correlation")
    y <- S[dataset %in% ds & set == st & factor == "SZ_effect"]
    out[[length(out) + 1]] <- data.table(variant = v, set = st, test = "SZ effect of recurrent genes, before vs after pH adjustment",
      n_datasets = nrow(y), median_signed_logfc_base = round(median(y$median_signed_logfc_base), 3),
      median_signed_logfc_pH = round(median(y$median_signed_logfc_pH), 3),
      median_of_dataset_median_ratio_pH = round(median(y$median_ratio_pH), 3),
      n_consistent_p05_base = sum(y$n_consistent_p05_base), n_consistent_p05_pH = sum(y$n_consistent_p05_pH),
      n_consistent_p05_allq = sum(y$n_consistent_p05_allq, na.rm = TRUE),
      pooled_signed_ratio_pH = round(sum(y$signed_sum_pH) / sum(y$signed_sum_base), 3),
      pooled_signed_ratio_allq = round(sum(y$signed_sum_allq, na.rm = TRUE) / sum(y$signed_sum_base[!is.na(y$signed_sum_allq)]), 3),
      datasets_signed_ratio_pH_below_075 = paste(y$dataset[y$signed_ratio_pH < 0.75], collapse = ","),
      note = paste("signed = multiplied by the recurrent SZ direction; median_of_dataset_median_ratio_pH = median over datasets of the median per-gene",
                   "logFC ratio (ignores sign); pooled_signed_ratio = sum over datasets and genes of signed logFC after / before (1 - ratio = attenuation)"))
  }
}
# expression-matched competitive test on the combined pH association: recurrent genes at combined q < 0.05 with
# the confound-consistent sign, vs 10,000 random gene sets matched on mean expression percentile (deciles) that
# carry the same SZ directions
set.seed(20260924)
cpH <- cp[factor == "pH"]
cpH[, bin := cut(expr_percentile, quantile(expr_percentile, 0:10 / 10, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
by_bin <- split(seq_len(nrow(cpH)), cpH$bin)
for (st in unique(rs$set)) {
  r <- rs[set == st]
  idx <- match(r$gene, cpH$gene); keep <- !is.na(idx); idx <- idx[keep]; sg <- r$sz_sign[keep]
  stat <- function(i) sum(cpH$q[i] < 0.05 & sign(cpH$Z[i]) == -sg)
  stat_any <- function(i) sum(cpH$q[i] < 0.05)
  need <- table(factor(cpH$bin[idx], levels = names(by_bin)))
  null <- replicate(10000, { s <- unlist(lapply(names(need), function(b) if (need[[b]]) sample(by_bin[[b]], need[[b]]))); c(stat(s), stat_any(s)) })
  out[[length(out) + 1]] <- data.table(variant = "primary_without_GSE208338", set = st, test = "combined control pH q < 0.05, expression-matched null (10,000 sets)",
    n_datasets = length(IND), observed = stat(idx), expected = round(mean(null[1, ]), 2), combined_p_one_sided = signif((1 + sum(null[1, ] >= stat(idx))) / 10001, 3),
    note = sprintf("observed = recurrent genes at q < 0.05 with the confound-consistent sign (of %d measured); any sign: observed %d, matched mean %.2f, p %.3g",
                   length(idx), stat_any(idx), mean(null[2, ]), (1 + sum(null[2, ] >= stat_any(idx))) / 10001))
}
O <- rbindlist(out, fill = TRUE)
fwrite(O, file.path(OUT, "combined_summary.tsv"), sep = "\t")

# per recurrent gene: combined pH association in controls; combined SZ effect before and after pH
sz <- rbindlist(lapply(IND, function(d) {
  x <- fread(file.path(OUT, sprintf("sz_adjust_%s.tsv", d))); x[, dataset := d]
  x[, .(dataset, gene, z_base = zt(t_base, p_base), z_pH = zt(t_pH, p_pH))]
}))
nn <- unique(S[factor == "SZ_effect", .(dataset, n = n_ctl + n_sz)], by = "dataset")
sz <- merge(sz, nn, by = "dataset")
szc <- sz[, .(n_datasets_sz = .N, Z_sz_base = stouffer(z_base, sqrt(n)), Z_sz_pH = stouffer(z_pH, sqrt(n))), by = gene]
g <- unique(rs[, .(gene, sz_direction, sz_sign)], by = "gene")
g[, sets := sapply(gene, function(x) paste(rs[gene == x, set], collapse = ","))]
g <- merge(g, cp[factor == "pH", .(gene, n_datasets_pH = n_datasets, Z_pH_controls = round(Z, 2), q_pH_controls = signif(q, 3))], by = "gene", all.x = TRUE)
g <- merge(g, cp[factor == "PMI", .(gene, Z_PMI_controls = round(Z, 2), q_PMI_controls = signif(q, 3))], by = "gene", all.x = TRUE)
g <- merge(g, szc[, .(gene, n_datasets_sz, Z_sz_base = round(Z_sz_base, 2), Z_sz_pH = round(Z_sz_pH, 2))], by = "gene", all.x = TRUE)
g[, pH_confound_consistent := sign(Z_pH_controls) == -sz_sign]
setorder(g, q_pH_controls)
fwrite(g, file.path(OUT, "recurrent_genes_pH.tsv"), sep = "\t")
options(width = 250)
print(O)
print(g)
cat("genes with combined pH q < 0.05 in controls:", cp[factor == "pH" & q < 0.05, .N], "of", cp[factor == "pH", .N], "\n")
