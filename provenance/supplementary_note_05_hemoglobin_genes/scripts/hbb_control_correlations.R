# hbb_control_correlations.R
# Producer for the Note 3 HBB sentences "median r 0.91 and 0.85 in 10 RNA-sequencing datasets" and
# "correlates with the erythrocyte score in 17 of 22 datasets (median r 0.47 in controls)".
# Base statement: the 28 sample-level datasets of analyses/controls_age_sex/data/datasets.tsv (22 bulk), as in analyses/hemoglobin_2026-10-11/02_sample_level_blood.R;
#   no gene list, no pH axis and no study base enters, so core-80 does not change it.
# Loader: analyses/fairness_sample_level/scripts/common.R (read only), bottom-tenth filter OFF (same as 02_sample_level_blood.R). Seed: none used (no randomness).
# Run: cd analyses/fairness_sample_level && Rscript <this file> <outdir>      Date: 2026-10-19
args <- commandArgs(TRUE); OUT <- normalizePath(args[1])
source("scripts/common.R")
HB <- c("HBB","HBA1","HBA2","HBD","HBG1","HBG2","HBM","HBQ1","HBE1","HBZ")
ERY <- c("HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA","GYPB","AHSP","CA1","BPGM","SLC25A37","EPB42")
NONHB <- setdiff(ERY, HB)
zrow <- function(E, g) { g <- intersect(g, rownames(E)); if (length(g) < 1) return(NULL); t(scale(t(E[g, , drop = FALSE]))) }
smean <- function(E, g, min_g = 3) { Z <- zrow(E, g); if (is.null(Z) || nrow(Z) < min_g) return(NULL); colMeans(Z, na.rm = TRUE) }
spc1 <- function(E, g, min_g = 3) { Z <- zrow(E, g); if (is.null(Z) || nrow(Z) < min_g) return(NULL)
  Z <- Z[apply(Z, 1, function(v) all(is.finite(v))), , drop = FALSE]; if (nrow(Z) < min_g) return(NULL)
  s <- prcomp(t(Z), center = TRUE, scale. = FALSE)$x[, 1]; if (cor(s, colMeans(Z)) < 0) s <- -s; s }
rows <- list()
for (nm in REG$dataset) {
  D <- as.list(REG[dataset == nm]); L <- load_dataset(D, bottom_filter = FALSE); E <- L$E; cv <- L$cv
  tech <- tryCatch({ sj <- rjson::fromJSON(file = cas_path(D$samples_json))$data
    paste(unique(vapply(sj, function(s) s$arrayDesign$technologyType, "")), collapse = "|") }, error = function(e) NA)
  sc <- list(ery_mean_z_incl_HBB = smean(E, ERY), ery_mean_z_excl_HBB = smean(E, setdiff(ERY, "HBB")), nonHb_PC1 = spc1(E, NONHB))
  for (scope in c("controls", "all")) {
    k <- if (scope == "controls") cv$sz == 0 else rep(TRUE, nrow(cv))
    r <- data.frame(dataset = nm, compartment = L$compartment, technology = tech, scope = scope, n = sum(k))
    cr <- function(a, b) if (a %in% rownames(E) && b %in% rownames(E)) cor(E[a, k], E[b, k]) else NA
    r$r_HBB_HBA1 <- cr("HBB", "HBA1"); r$r_HBB_HBA2 <- cr("HBB", "HBA2")
    for (s in names(sc)) r[[paste0("r_HBB_", s)]] <- if (!is.null(sc[[s]]) && "HBB" %in% rownames(E)) cor(E["HBB", k], sc[[s]][k]) else NA
    rows[[length(rows) + 1]] <- r
  }
  cat(nm, ncol(E), "\n")
}
R <- do.call(rbind, rows)
con <- file(file.path(OUT, "hbb_control_correlations.tsv"), "w")
writeLines(c("# script: CD/scripts/hbb_control_correlations.R (run 2026-10-19, laptop, R base + data.table/limma, no randomness)",
  "# input: analyses/controls_age_sex/data/datasets.tsv sha256 (see REPORT.md); loader analyses/fairness_sample_level/scripts/common.R",
  "# base: 28 sample-level datasets (22 bulk); Pearson r over samples (controls = sz == 0); log2 scale as delivered by Gemma; bottom-tenth filter off",
  "# technology = Gemma arrayDesign technologyType of the experiment"), con)
write.table(R, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
