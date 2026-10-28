# NEW 2026-10-18, analysis 2. Schizophrenia t per gene at each pseudobulk level in PsychAD (HBCC_main, MSSM_main donor sets),
# correlated with the bulk pH axis v2 and the agonal axis. Reuses the pipeline of rerun-full/analyses/psychad_cellxgene
# (scripts/pb_common.R: load_donors, donor_set, load_level, cell_samples, eligible; model ~ sz + male + age + pmi + log(nuclei);
# pseudobulks >= 20 nuclei; >= 20 donors per group; filterByExpr(group = sz), TMM, voom, limma). Adds the 'total with composition
# covariates' fit (class fractions of nuclei as covariates, smallest class dropped) and donor-bootstrap intervals.
# Run: RUN_DIR=ph_levels_2026-10-18 NCORES=8 nice -n 19 Rscript psychad_levels.R
PCDIR <- "psychad_pseudobulk/gl/psychad_cellxgene"
setwd(PCDIR)
SCRIPT <- "ph_levels_2026-10-18/scripts/psychad_levels.R"
source("scripts/pb_common.R")
source(file.path(Sys.getenv("RUN_DIR", normalizePath("ph_levels_2026-10-18")), "scripts/levels_common.R"))  # resets SEED, NCORES
OUTD <- file.path(RUN, "out"); dir.create(OUTD, showWarnings = FALSE)
AX <- load_axes()
D <- load_donors()
DEDIR <- "rerun-full/analyses/psychad_cellxgene/results/de"
RHS <- "sz + male + age + pmi + lognuc"
tabs <- list(); gene_t <- list(); chk <- list(); donor_tabs <- list(); bootrows <- list()

for (bank in c("HBCC", "MSSM")) {
  set <- paste0(bank, "_main")
  Ls <- list(total = load_level(bank, "total"), class = load_level(bank, "class"), subclass = load_level(bank, "subclass"))
  gmap <- setNames(Ls$total$genes$gene_name, Ls$total$genes$ensembl)
  old <- fread(cmd = paste("zcat", shQuote(file.path(DEDIR, paste0(set, ".tsv.gz"))), "| grep -v '^#'"))
  units <- list(); ud <- list()
  for (lev in names(Ls)) {
    L <- Ls[[lev]]
    for (ct in unique(L$cols$cell_type)) {
      d <- cell_samples(L, D, set, ct)
      if (!eligible(d)) next
      d <- as.data.frame(d); d$donor <- d$donor_id
      u <- fit_unit(as.matrix(L$M[, d$pb_id, drop = FALSE]), d, RHS, "sz", group = d$sz, gmap = gmap)
      if (is.null(u)) next
      units[[paste(lev, ct, sep = "|")]] <- u
      o <- old[level == lev & cell_type == ct]
      if (nrow(o)) { g <- intersect(names(u$t), o$gene); chk[[length(chk) + 1]] <- data.table(bank = bank, level = lev, cell_type = ct,
        n_genes_mine = length(u$t), n_genes_stored = nrow(o), n_shared = length(g),
        pearson_t = cor(u$t[g], o$t[match(g, o$gene)]), max_abs_diff_t = max(abs(u$t[g] - o$t[match(g, o$gene)]))) }
    }
  }
  # total with composition covariates: fractions of each donor's nuclei in each class (class-level pseudobulk n_nuclei)
  dt <- as.data.frame(cell_samples(Ls$total, D, set, "all")); dt$donor <- dt$donor_id
  cc <- as.data.frame(Ls$class$cols[donor_id %in% dt$donor_id, c("donor_id", "cell_type", "n_nuclei")])
  wide <- reshape(cc, idvar = "donor_id", timevar = "cell_type", direction = "wide"); wide[is.na(wide)] <- 0
  cn <- sub("^n_nuclei\\.", "", colnames(wide)[-1]); mat <- as.matrix(wide[, -1]); colnames(mat) <- cn; rownames(mat) <- wide$donor_id
  frac <- mat / rowSums(mat)
  smallest <- names(which.min(colSums(mat)))
  keepc <- setdiff(cn, smallest); fn <- paste0("frac_", make.names(keepc))
  for (i in seq_along(keepc)) dt[[fn[i]]] <- frac[match(dt$donor_id, rownames(frac)), keepc[i]]
  stopifnot(!anyNA(dt[, fn]))
  uc <- fit_unit(as.matrix(Ls$total$M[, dt$pb_id, drop = FALSE]), dt, paste(RHS, "+", paste(fn, collapse = " + ")), "sz", group = dt$sz, gmap = gmap)
  units[["total_composition|all"]] <- uc
  cat(sprintf("[%s] composition covariates: %s; dropped smallest class %s\n", bank, paste(keepc, collapse = ","), smallest))
  donors <- dt$donor_id; strata <- dt$sz[match(donors, dt$donor_id)]
  ident <- paste(bank, "SZ")
  an <- analyse_units(units, donors, strata, ident)
  tab <- an$table[, `:=`(dataset = paste0("PsychAD_", bank), analysis = "2_SZ_effect_vs_axes", population = "SZ vs control (donor set main)")]
  tab[, n_donors_total_fit := length(donors)]
  tabs[[bank]] <- tab
  gene_t[[bank]] <- rbindlist(lapply(names(units), function(k) { u <- units[[k]]
    data.table(dataset = paste0("PsychAD_", bank), level = sub("\\|.*", "", k), cell_type = sub("^[^|]*\\|", "", k),
               gene = names(u$t), t = round(u$t, 4), logfc = round(u$logfc, 5), p = signif(u$p, 4)) }))
  donor_tabs[[bank]] <- data.table(dataset = paste0("PsychAD_", bank), unit = names(units), n_samples = vapply(units, function(u) u$n, 1),
                                   n_sz = vapply(units, function(u) sum(u$d$sz == 1), 1), n_ctl = vapply(units, function(u) sum(u$d$sz == 0), 1),
                                   n_nuclei = vapply(units, function(u) u$nuclei, 1), n_genes = vapply(units, function(u) length(u$t), 1),
                                   df = vapply(units, function(u) u$df, 1))
  B <- an$boot; dimnames(B) <- list(unit = an$names, stat = STATS, boot = seq_len(dim(B)[3]))
  bootrows[[bank]] <- as.data.table(as.data.frame.table(B, responseName = "rho"))[, dataset := paste0("PsychAD_", bank)]
  rm(Ls, units, an, B); gc()
}
hdr <- c(sprintf("NEW 2026-10-18 analysis 2 (PsychAD SZ effect by level) by ph_levels_2026-10-18/scripts/psychad_levels.R; seed %d; %d donor resamples (within SZ and control)", SEED, NBOOT),
         sprintf("inputs sha256(16): axis_v2_controls.tsv.gz %s; hardy_t_by_region.tsv.gz %s; recurrent_sets.tsv %s", SHA("axis_v2_controls.tsv.gz"), SHA("hardy_t_by_region.tsv.gz"), SHA("recurrent_sets.tsv")),
         sprintf("PsychAD pseudobulk counts: %s", paste(sprintf("%s/%s %s", rep(c("HBCC","MSSM"), each = 3), rep(c("total","class","subclass"), 2),
           sapply(c("HBCC","MSSM"), function(b) sapply(c("total","class","subclass"), function(l) substr(sha256(file.path("psychad_pseudobulk/out", paste0(b, "_Cohort"), sprintf("counts_%s.npz", l))), 1, 12)))), collapse = "; ")),
         "t = limma moderated t of sz (sz=1 schizophrenia) from ~ sz + male + age + pmi + lognuc (+ class fractions for total_composition); filterByExpr(group=sz), TMM, voom; >= 20 nuclei per pseudobulk, >= 20 donors per group",
         "rho = Spearman across genes of that t with the bulk pH axis v2 (z_stouffer, positive = higher at higher pH) or the agonal axis (Hardy t slow vs fast death, mean of BA9, BA24, Cortex); gene_set all = all shared genes, rec383 = study_brain, rec107 = study_brain_ge10",
         "interval: rho +- 1.96 SD of Fisher-z of rho over donor bootstrap resamples; moderated t recomputed with limma lmFit + eBayes on the resample with voom weights held at the full-data values (voom not re-estimated); pct_lo/pct_hi = 2.5/97.5 percentiles",
         "MEDIAN_over_cell_types_weighted_nuclei = weighted median of per-cell-type rho (weight = nuclei summed over donors in that cell type's fit); unweighted = plain median; the bootstrap repeats the median on every resample")
write_tsv_prov(rbindlist(tabs), file.path(OUTD, "psychad_levels_rho.tsv"), hdr)
write_tsv_prov(rbindlist(gene_t), file.path(OUTD, "psychad_levels_gene_t.tsv.gz"), hdr[1:2])
write_tsv_prov(rbindlist(donor_tabs), file.path(OUTD, "psychad_levels_units.tsv"), hdr[1])
write_tsv_prov(rbindlist(bootrows), file.path(OUTD, "psychad_levels_bootstrap_rho.tsv.gz"), hdr[1])
write_tsv_prov(rbindlist(chk), file.path(OUTD, "psychad_levels_check_vs_stored_DE.tsv"),
               c("agreement of this script's fit with results/de/<set>.tsv.gz (stored limma t) per unit", hdr[1]))
cat("done\n")
