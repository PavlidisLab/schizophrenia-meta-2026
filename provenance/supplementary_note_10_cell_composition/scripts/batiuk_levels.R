# NEW 2026-10-18, analysis 1 on Batiuk 2022 (Gemma 57261; adapted from nsw_levels.R: neurons plus a few glial types as labelled in Gemma, 18 cell types, 3 classes).
# Only the all-donor population is analysable (controls with pH: 5). MIN_DON is lowered to 10 because 11 donors have pH; exploratory.
# (original header) analysis 1. pH association per gene (limma t of per-donor pH) in the Batiuk snRNA-seq (Gemma 57261)
# at each pseudobulk level built from the same nuclei by nsw_build_pb.py: total (all nuclei), total with class-fraction covariates, class, subclass.
# Populations: controls only (n from Control donors) and 'all' = Control + Schizophrenia + Schizoaffective (the 79 persons of the axis v2 inventory) with diagnosis as a covariate.
# Model: ~ ph + age + sex + pmi + log(nuclei) [+ dx]; pseudobulks >= 20 nuclei, >= 10 donors per cell type; filterByExpr(design), TMM, voom, limma.
source(file.path(Sys.getenv("RUN_DIR", normalizePath("ph_levels_2026-10-18")), "scripts/levels_common.R"))
AX <- load_axes()
PB <- file.path(RUN, "out", "batiuk_pb"); OUTD <- file.path(RUN, "out")
MIN_NUC <- 20; MIN_DON <- 10
cols <- fread(file.path(PB, "pb_cols.tsv"))
cov <- fread(file.path(RUN, "out", "batiuk_donor_covariates.tsv"))
cov[, sex01 := as.integer(sex == "Male")]
cov[, dxf := factor(dx)]
mats <- lapply(c(total = "total", class = "class", subclass = "subclass"), function(l) {
  m <- fread(cmd = paste("zcat", shQuote(file.path(PB, paste0("pb_", l, ".tsv.gz")))), header = TRUE)
  ids <- m[[1]]; m <- as.matrix(m[, -1]); rownames(m) <- ids; t(m)   # genes x samples
})
tabs <- list(); gene_t <- list(); unit_tabs <- list(); bootrows <- list()
for (pop in c("all")) {
  cv <- cov
  cv <- cv[is.finite(ph) & is.finite(age) & is.finite(pmi) & !is.na(sex01)]
  rhs <- "ph + age + sex01 + pmi + lognuc"
  if (pop == "all") { cv[, dxf := droplevels(relevel(factor(dx), "CTRL"))]; rhs <- paste(rhs, "+ dxf") }
  units <- list()
  for (lev in c("total", "class", "subclass")) {
    cl <- cols[level == lev & n_nuclei >= MIN_NUC]
    for (ct in unique(cl$cell_type)) {
      d <- merge(cl[cell_type == ct], cv, by = "donor"); d <- as.data.frame(d[order(donor)])
      if (nrow(d) < MIN_DON) next
      d$lognuc <- log(d$n_nuclei)
      u <- fit_unit(mats[[lev]][, d$pb_id, drop = FALSE], d, rhs, "ph")
      if (!is.null(u)) units[[paste(lev, ct, sep = "|")]] <- u
    }
  }
  dt <- as.data.frame(merge(cols[level == "total" & n_nuclei >= MIN_NUC], cv, by = "donor")); dt <- dt[order(dt$donor), ]; dt$lognuc <- log(dt$n_nuclei)
  cc <- cols[level == "class" & donor %in% dt$donor, .(donor, cell_type, n_nuclei)]
  wide <- dcast(cc, donor ~ cell_type, value.var = "n_nuclei", fill = 0); mat <- as.matrix(wide[, -1]); rownames(mat) <- wide$donor
  frac <- mat / rowSums(mat); smallest <- names(which.min(colSums(mat))); keepc <- setdiff(colnames(mat), smallest); fn <- paste0("frac_", make.names(keepc))
  for (i in seq_along(keepc)) dt[[fn[i]]] <- frac[dt$donor, keepc[i]]
  uc <- fit_unit(mats$total[, dt$pb_id, drop = FALSE], dt, paste(rhs, "+", paste(fn, collapse = " + ")), "ph")
  units[["total_composition|all"]] <- uc
  cat(sprintf("[Batiuk %s] composition covariates %s; dropped %s; total fit n=%d\n", pop, paste(keepc, collapse = ","), smallest, nrow(dt)))
  an <- analyse_units(units, dt$donor, NULL, paste("Batiuk", pop))
  tab <- an$table[, `:=`(dataset = "Batiuk_2022", analysis = "1_pH_association_vs_axes",
         population = "all donors, diagnosis adjusted (control + SZ)")]
  tabs[[pop]] <- tab
  gene_t[[pop]] <- rbindlist(lapply(names(units), function(k) { u <- units[[k]]
    data.table(population = pop, level = sub("\\|.*", "", k), cell_type = sub("^[^|]*\\|", "", k), gene = names(u$t), t = round(u$t, 4), logfc = round(u$logfc, 5), p = signif(u$p, 4)) }))
  unit_tabs[[pop]] <- data.table(population = pop, unit = names(units), n_donors = vapply(units, function(u) u$n, 1), n_nuclei = vapply(units, function(u) u$nuclei, 1),
                                 n_genes = vapply(units, function(u) length(u$t), 1), df = vapply(units, function(u) u$df, 1), ph_sd = vapply(units, function(u) sd(u$d$ph), 1))
  B <- an$boot; dimnames(B) <- list(unit = an$names, stat = STATS, boot = seq_len(dim(B)[3]))
  bootrows[[pop]] <- as.data.table(as.data.frame.table(B, responseName = "rho"))[, population := pop]
  rm(units, an, B); gc()
}
hdr <- c(sprintf("NEW 2026-10-18 analysis 1 (Batiuk 2022 pH association by level) by ph_levels_2026-10-18/scripts/batiuk_levels.R; seed %d; %d donor resamples", SEED, NBOOT),
  sprintf("inputs sha256(16): batiuk_donor_covariates.tsv %s; axis_v2_controls.tsv.gz %s; hardy_t_by_region.tsv.gz %s; recurrent_sets.tsv %s; pseudobulk from Gemma 57261 aggregated log2cpm (input_batiuk/manifest.tsv)", substr(sha256(file.path(RUN, "out", "batiuk_donor_covariates.tsv")), 1, 16), SHA("axis_v2_controls.tsv.gz"), SHA("hardy_t_by_region.tsv.gz"), SHA("recurrent_sets.tsv")),
  "t = limma moderated t of per-donor pH from ~ ph + age + sex + pmi + lognuc (+ dx; + class fractions for total_composition); filterByExpr(design), TMM, voom; >= 20 nuclei per pseudobulk, >= 10 donors per cell type",
  "rho = Spearman across genes of that t with the bulk pH axis v2 (positive = higher expression at higher pH) or the agonal axis (Hardy t); gene_set all = all shared genes, rec383 = study_brain, rec107 = study_brain_ge10",
  "interval: rho +- 1.96 SD of Fisher-z of rho over donor bootstrap resamples (donors drawn with replacement); moderated t recomputed with lmFit + eBayes on the resample with voom weights held at full-data values; pct_lo/pct_hi = 2.5/97.5 percentiles",
  "MEDIAN_over_cell_types_weighted_nuclei = median of per-cell-type rho weighted by nuclei summed over donors; unweighted = plain median")
write_tsv_prov(rbindlist(tabs), file.path(OUTD, "batiuk_levels_rho.tsv"), hdr)
write_tsv_prov(rbindlist(gene_t), file.path(OUTD, "batiuk_levels_gene_t.tsv.gz"), hdr[1:2])
write_tsv_prov(rbindlist(unit_tabs), file.path(OUTD, "batiuk_levels_units.tsv"), hdr[1])
write_tsv_prov(rbindlist(bootrows), file.path(OUTD, "batiuk_levels_bootstrap_rho.tsv.gz"), hdr[1])
cat("done\n")
