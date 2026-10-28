# G. Same Pittsburgh donors: can each donor be picked out from the other dataset by expression profile alone?
# Paul, 2026-09-28: "can you predict which person this was based on their similarity to an external profile, relative to
# other samples in the same data? Like Gillis' MetaNeighbor method sort of."
#
# Run from analyses/sample_level_reanalysis_of_datasets/ (the directory common.R expects), after f_fullstats.R:
#   Rscript same_donors_pitt/g_individuals.R > same_donors_pitt/results/g_individuals.log 2>&1
#
# The test is the identification step of ../psychad_cellxgene/scripts/e_individuals.R (PsychAD HBCC x LIBD DLPFC) and
# ../libd_brainseq/analysis/scripts/e_individuals_regions.R (LIBD DLPFC x hippocampus); its functions are sourced from
# e_individuals.R (main() does not run): resid_mat, identify, strat_perm, identify_auroc.
# Units: the eight per-sample Pittsburgh units of f_fullstats.R that carry a donor identity per sample: GSE53987 BA46 /
# hippocampus / striatum (Lanz 2019), GSE202537 NAc / caudate / putamen (Ketchesin 2023), GSE93987 L3 / L5 (Arion 2015).
# Loaded exactly as f_fullstats.R loads them (common.R load_dataset: SZ and controls with age and sex, replicate arrays
# of one donor averaged, bottom tenth of mean expression dropped); GSE93987's covariate file is the one f_fullstats.R
# wrote (data/cov/GSE93987.tsv), read here, not rewritten. GSE87610 (Arion 2017) has no donor identity per sample and is
# not used.
# Donor identity: across studies, the Pittsburgh case number (results/donor_codes.tsv by GSM; GSE202537's subject =
# HU#); within one study, the series' own donor label (as f_fullstats.R's shared(): GSE53987 GEO subject, GSE202537 HU#,
# GSE93987 case number). Candidates: donors in both units with the same diagnosis (SZ or control) in both.
# Residuals: per unit, every gene regressed on age and sex (the covariates of the f_fullstats.R fits) over all of the
# unit's donors, without diagnosis ('nodx') or with it ('dx'). Genes: those in both units, non-constant over the candidates.
# Identification: identify() (Pearson over genes of the residual profiles z-scored within the candidates; rank of the true
# partner), and identify_auroc() (AUROC = 1 - (rank - 1) / (n - 1), both directions; null: pairing shuffled within
# dx x sex, NPERM shuffles).
# Outputs (same_donors_pitt/results/): g_individuals_identify.tsv (per donor), g_individuals_summary.tsv (per pair, mode and
# direction), g_individuals_meta.json; g_individuals_C.rds (each pair's donor x donor similarity matrix, with donor ids,
# diagnosis and sex; read by i_dx_neighbours.R)
source("scripts/common.R")
suppressMessages(library(jsonlite))
PE <- new.env()
local({ owd <- setwd("../psychad_cellxgene"); on.exit(setwd(owd)); sys.source("scripts/e_individuals.R", envir = PE) })
for (f in c("resid_mat", "identify", "strat_perm", "auroc_of_rank", "identify_ranks", "identify_auroc")) assign(f, get(f, envir = PE))
HERE <- "same_donors_pitt"
RES <- file.path(HERE, "results")
SCRIPT <- "fairness_sample_level/same_donors_pitt/g_individuals.R"
NPERM <- 1000
SEED_G <- get("SEED", envir = PE) + 13
COVS <- c("age", "male")
INS <- list()
note_in <- function(f) { INS[[f]] <<- sha256(f); f }

# sample fields of the two series outside data/gemma_batch: the wrapper of f_fullstats.R
gemma_sample_fields_common <- gemma_sample_fields
gemma_sample_fields <- function(gid) {
  f <- sprintf("data/gemma_batch/%s_samples.json", gid)
  if (file.exists(f)) return(gemma_sample_fields_common(gid))
  sj <- jsonlite::fromJSON(note_in(file.path(HERE, "data/gemma", sprintf("%s_samples.json", gid))), simplifyVector = FALSE)$data
  data.table(gsm = vapply(sj, function(s) s$accession$accession, ""), batch = NA_character_, drug = NA_character_)
}
codes <- fread(note_in(file.path(RES, "donor_codes.tsv")), skip = "dataset\t", colClasses = "character")
cov93 <- normalizePath(note_in(file.path(HERE, "data/cov/GSE93987.tsv")))
e93 <- normalizePath(note_in(file.path(HERE, "data/gemma/10879_processed.tsv.gz")))
s93 <- normalizePath(note_in(file.path(HERE, "data/gemma/10879_samples.json")))
NEW <- list(GSE93987_L3 = list(dataset = "GSE93987_L3", expr = e93, samples_json = s93, cov = cov93, region = "L3", compartment = "LCM layer 3 pyramidal neurons",
                               cluster = "pitt", gse = "GSE93987", gemma_id = "10879"),
            GSE93987_L5 = list(dataset = "GSE93987_L5", expr = e93, samples_json = s93, cov = cov93, region = "L5", compartment = "LCM layer 5 pyramidal neurons",
                               cluster = "pitt", gse = "GSE93987", gemma_id = "10879"))
UA <- rbindlist(list(
  list("GSE53987_BA46", "Lanz 2019", "DLPFC (BA46)", "DLPFC", "bulk", "Affymetrix U133 Plus 2.0"),
  list("GSE53987_hippocampus", "Lanz 2019", "hippocampus", "hippocampus", "bulk", "Affymetrix U133 Plus 2.0"),
  list("GSE53987_striatum", "Lanz 2019", "associative striatum", "striatum:associative", "bulk", "Affymetrix U133 Plus 2.0"),
  list("GSE202537_NAc", "Ketchesin 2023", "nucleus accumbens", "striatum:NAc", "bulk", "RNA-seq (total RNA)"),
  list("GSE202537_caudate", "Ketchesin 2023", "caudate", "striatum:caudate", "bulk", "RNA-seq (total RNA)"),
  list("GSE202537_putamen", "Ketchesin 2023", "putamen", "striatum:putamen", "bulk", "RNA-seq (total RNA)"),
  list("GSE93987_L3", "Arion 2015", "DLPFC", "DLPFC", "LCM layer 3 pyramidal neurons", "Affymetrix HT HG-U133+ PM"),
  list("GSE93987_L5", "Arion 2015", "DLPFC", "DLPFC", "LCM layer 5 pyramidal neurons", "Affymetrix HT HG-U133+ PM")))
setnames(UA, c("unit", "study", "region", "region_key", "compartment", "platform"))
for (u in UA$unit) if (!u %in% names(NEW)) for (k in c("expr", "samples_json", "cov")) invisible(note_in(cas_path(REG[dataset == u][[k]])))

X <- list()
for (u in UA$unit) {
  D <- if (u %in% names(NEW)) NEW[[u]] else as.list(REG[dataset == u])
  if (u %in% names(NEW)) { CAS_common <- CAS; CAS <- HERE }
  x <- load_dataset(D)
  if (u %in% names(NEW)) CAS <- CAS_common
  cv <- x$cv
  cv[, code := if (startsWith(u, "GSE202537")) subject else codes$code[match(gsm, codes$gsm)]]
  cv[, donor := if (startsWith(u, "GSE53987")) codes$subject[match(gsm, codes$gsm)] else subject]
  cv[, dx := ifelse(sz == 1, "SZ", "CTRL")]; cv[, sex := ifelse(male == 1, "M", "F")]
  stopifnot(!anyDuplicated(cv$donor), all(cv$donor != "" & !is.na(cv$donor)))
  R <- lapply(c(nodx = FALSE, dx = TRUE), function(dx) resid_mat(x$E[, cv$gsm], cv, COVS, dx = dx))
  X[[u]] <- list(cv = cv, R = R, E = x$E[, cv$gsm, drop = FALSE])   # E: raw (pre-residual) matrix, for UNITS_RDS reuse only; unread below
  cat(sprintf("%s: %d SZ, %d CTRL, %d genes; case numbers %d\n", u, sum(cv$sz), sum(1 - cv$sz), nrow(x$E), sum(!is.na(cv$code) & cv$code != "")))
}

# Opt-in, backward-compatible hook (default: both env vars unset, no effect): dump the per-unit cv (with pH/PMI/RIN/
# agonal columns, when load_dataset found them) and residual matrices, for reanalysis under different gene sets or
# covariates without duplicating this loading loop (identification_ph_2026-10-01/scripts/load_relations.R).
# UNITS_ONLY also stops before the pair loop, so a units-only run is fast.
if (nzchar(Sys.getenv("UNITS_RDS"))) {
  saveRDS(X, Sys.getenv("UNITS_RDS"))
  if (nzchar(Sys.getenv("UNITS_ONLY"))) { cat("UNITS_ONLY: wrote", Sys.getenv("UNITS_RDS"), "and stopped before the pair loop\n"); quit(save = "no") }
}

region_rel <- function(a, b) {   # f_fullstats.R's rule
  if (a == b) return("same region")
  if (startsWith(a, "striatum") && startsWith(b, "striatum")) return("overlapping region (striatal subregions)")
  "different region"
}
PD <- list(); PS <- list(); CM <- list()
us <- UA$unit
for (i in seq_along(us)) for (j in seq_along(us)) if (i < j) {
  a <- us[i]; b <- us[j]; ua <- UA[unit == a]; ub <- UA[unit == b]
  same_study <- ua$study == ub$study
  key <- if (same_study) "donor" else "code"
  A <- X[[a]]$cv[!is.na(get(key)) & get(key) != ""]; B <- X[[b]]$cv[!is.na(get(key)) & get(key) != ""]
  m <- merge(A[, .(id = get(key), gsm_a = gsm, dx_a = dx, sex_a = sex)], B[, .(id = get(key), gsm_b = gsm, dx_b = dx, sex_b = sex)], by = "id")
  n_dx_mismatch <- sum(m$dx_a != m$dx_b); n_sex_mismatch <- sum(m$dx_a == m$dx_b & m$sex_a != m$sex_b)
  m <- m[dx_a == dx_b][order(id)]
  if (nrow(m) < 3) next
  strata <- paste(m$dx_a, m$sex_a)
  genes <- intersect(rownames(X[[a]]$R$nodx), rownames(X[[b]]$R$nodx))
  for (mode in c("nodx", "dx")) {
    Am <- X[[a]]$R[[mode]][genes, m$gsm_a, drop = FALSE]; Bm <- X[[b]]$R[[mode]][genes, m$gsm_b, drop = FALSE]
    ok <- apply(Am, 1, sd) > 0 & apply(Bm, 1, sd) > 0
    idn <- identify(Am[ok, , drop = FALSE], Bm[ok, , drop = FALSE])
    ia <- identify_auroc(idn$C, strata, NPERM, SEED_G)
    stopifnot(all(ia$per_donor$rank_ab == idn$rank))
    CM[[length(CM) + 1]] <- list(unit_a = a, unit_b = b, mode = mode, C = idn$C, id = m$id, dx = m$dx_a, sex = m$sex_a)
    base <- data.table(unit_a = a, unit_b = b, study_a = ua$study, study_b = ub$study, series = if (same_study) "same study" else "different studies",
                       region_a = ua$region, region_b = ub$region, region_relation = region_rel(ua$region_key, ub$region_key),
                       compartment_relation = if (ua$compartment == ub$compartment) "same compartment" else "different compartment",
                       platform_a = ua$platform, platform_b = ub$platform, platform_relation = if (ua$platform == ub$platform) "same platform" else "different platform",
                       donor_key = if (same_study) "series donor label" else "Pittsburgh case number", mode = mode)
    PD[[length(PD) + 1]] <- data.table(base[, .(unit_a, unit_b, mode)], donor = m$id, dx = m$dx_a, sex = m$sex_a, n_candidates = nrow(m), ia$per_donor)
    PS[[length(PS) + 1]] <- data.table(base, n_SZ = sum(m$dx_a == "SZ"), n_CTRL = sum(m$dx_a == "CTRL"), n_genes = sum(ok),
                                       n_dx_mismatch_excluded = n_dx_mismatch, n_sex_mismatch = n_sex_mismatch, ia$summary)
    if (mode == "nodx") cat(sprintf("%s x %s: %d candidates, %d genes; top-1 %d / %d, mean AUROC %.3f / %.3f (null %.3f), p %.3g\n", a, b, nrow(m), sum(ok),
                                    ia$summary$top1[1], ia$summary$top1[2], ia$summary$mean_auroc[1], ia$summary$mean_auroc[2], ia$summary$null_mean_auroc[1], ia$summary$p_perm[1]))
  }
}
PD <- rbindlist(PD); PS <- rbindlist(PS)
saveRDS(CM, file.path(RES, "g_individuals_C.rds"))   # similarity matrices per pair and mode, for i_dx_neighbours.R
hdr <- function(f, lines) { writeLines(paste("#", lines), f); f }
stamp <- sprintf("built %s by %s; identification of ../psychad_cellxgene/scripts/e_individuals.R (identify, identify_auroc); residuals on %s over each unit's donors ('nodx'; 'dx' adds diagnosis); NPERM %d; seed %d",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT, paste(COVS, collapse = " + "), NPERM, SEED_G)
f1 <- hdr(file.path(RES, "g_individuals_identify.tsv"), c(stamp, "per candidate donor: rank of the true partner, a -> b (ab: among unit b's candidates) and b -> a (ba); auroc = 1 - (rank - 1) / (n_candidates - 1); donor = case number (different studies) or the series' donor label (same study)"))
fwrite(PD, f1, sep = "\t", append = TRUE, col.names = TRUE)
f2 <- hdr(file.path(RES, "g_individuals_summary.tsv"), c(stamp, "per unit pair x residual mode x direction; candidates: donors in both units with the same diagnosis (SZ or control); null = pairing shuffled within dx x sex; p = (1 + #null mean AUROC >= observed) / (NPERM + 1); n_strata = expected number of donors a shuffle leaves with their true partner"))
fwrite(PS, f2, sep = "\t", append = TRUE, col.names = TRUE)
write_json(list(script = SCRIPT, built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), R = R.version.string, nperm = NPERM, seed = SEED_G, covs = COVS,
                method = "MetaNeighbor-style identification AUROC, Crow et al. 2018 Nat Commun doi 10.1038/s41467-018-03282-0",
                functions_from = "../psychad_cellxgene/scripts/e_individuals.R", functions_from_sha256 = sha256("../psychad_cellxgene/scripts/e_individuals.R"),
                inputs_sha256 = INS, outputs_sha256 = list(g_individuals_identify.tsv = sha256(f1), g_individuals_summary.tsv = sha256(f2))),
           file.path(RES, "g_individuals_meta.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
