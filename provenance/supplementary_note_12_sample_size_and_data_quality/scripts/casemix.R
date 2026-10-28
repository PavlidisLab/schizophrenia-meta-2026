# Does case mix explain disagreement between datasets? (paper outline section 6; cohorts/cohorts_report.md)
#
# Run from analyses/case_mix_agreement_test/:  Rscript scripts/casemix.R        (about 30-60 minutes with NCORES=6)
#
# Units: the bulk tissue datasets of the fairness_sample_level registry (../controls_age_sex/data/datasets.tsv, compartment
# 'bulk'), LIBD BrainSeq DLPFC and hippocampus (../libd_brainseq/analysis/data/libd_<region>.rds, adults as prep.R built
# them) and PsychAD HBCC and Mount Sinai whole tissue (../psychad_cellxgene/pseudobulk, level 'total', donor sets
# HBCC_main and MSSM_main). Relations between units: ../subject_relations/subject_relations.tsv (relation_group);
# 'subject-separate' pairs are the primary set, other groups are reported apart and never pooled.
#
# SZ vs control per unit, one model family per data type, age and sex as covariates in all (fairness level L1):
#   Gemma arrays / RNA-seq  limma ~ sz + age + male (common.R fit_sz)
#   LIBD                    limma-trend on logCPM ~ sz + age + male
#   PsychAD pseudobulk      voom-limma ~ sz + male + age + PMI + log(nuclei) (pb_common.R fit_pb)
# Agreement between two units: common.R agree() (rho_t = Spearman of t over shared genes; top100_fold).
#
# Case mix per unit (the SZ group of the samples fitted): median age at death, SZ minus control median age, median pH,
# median PMI, suicide share, pre-mortem assessment. Suicide per sample where the data record it (GSE92538 'suicide (1=yes)',
# GSE208338 'suicide' Y/N, GSE202537 'manner of death', LIBD Manner_Of_Death); otherwise the bank-level share from
# ../cohorts/cohort_diagnostics.tsv (source column suicide_source says which). Pre-mortem assessment: banks mssm and penn
# (cohorts_report.md); illness duration is recorded in 2 units only and is described, not tested.
#
# Test 1 (pairs) and its power: scripts/casemix_pairs.R, from the tables written here.
# Test 2 (intervention): each unit refitted on case-mix-comparable subsets (age at death <= 65; suicides excluded where
# recorded; both), agreement recomputed; null = random subsets of the same per-unit SZ and control counts (NREP), so the
# smaller n is not read as an effect. Statistics: median rho_t and top100_fold over the subject-separate pairs, the same
# over the pairs with a unit that lost >= 20% of a group (affected), and single pairs PsychAD HBCC vs Mount Sinai and
# LIBD DLPFC vs Mount Sinai.
suppressMessages({library(data.table); library(limma); library(edgeR); library(parallel); library(Matrix)})
CM <- normalizePath(".")
GL <- normalizePath("..")
FSL <- file.path(GL, "fairness_sample_level")
PSY <- file.path(GL, "psychad_cellxgene")
local({ owd <- setwd(FSL); on.exit(setwd(owd)); source("scripts/common.R", local = FALSE) })
local({ owd <- setwd(PSY); on.exit(setwd(owd)); source("scripts/pb_common.R", local = FALSE) })
setwd(CM)
SCRIPT <- "scripts/casemix.R"
NCORES <- as.integer(Sys.getenv("NCORES", "6"))
SEED <- 20260926
NREP <- as.integer(Sys.getenv("NREP", "100"))
AGE_MAX <- 65
MIN_GROUP <- 5
dir.create("results", showWarnings = FALSE)

writeT <- function(x, path, header = character(0)) {
  writeLines(paste("#", c(sprintf("built %s by casemix/%s", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT), header)), path)
  fwrite(x, path, sep = "\t", append = TRUE, col.names = TRUE)
}

# ---------------------------------------------------------------- units
units_tsv <- fread(file.path(GL, "subject_relations/subject_relations_units.tsv"), skip = "unit\t")
bank_of <- setNames(units_tsv$bank, units_tsv$unit)

per_sample_suicide <- function(gid, gsm) {
  f <- file.path(FSL, sprintf("data/gemma_batch/%s_samples.json", gid))
  if (!file.exists(f)) return(rep(NA_integer_, length(gsm)))
  sj <- jsonlite::fromJSON(f, simplifyVector = FALSE)$data
  v <- vapply(sj, function(s) {
    x <- NA_character_
    for (ch in s$sample$characteristics) {
      cat_ <- ch$category
      if (identical(cat_, "suicide (1=yes)")) x <- ifelse(ch$value == "1", "1", "0")
      if (identical(cat_, "suicide")) x <- ifelse(ch$value == "Y", "1", "0")
      if (identical(cat_, "manner of death")) x <- ifelse(ch$value == "Suicide", "1", "0")
    }
    x
  }, "")
  acc <- vapply(sj, function(s) s$accession$accession, "")
  as.integer(v[match(gsm, acc)])
}

message("loading Gemma bulk datasets")
bulk <- REG[compartment == "bulk"]
G <- local({ owd <- setwd(FSL); on.exit(setwd(owd)); lapply(seq_len(nrow(bulk)), function(i) {
  x <- load_dataset(as.list(bulk[i]))
  x$cv[, suicide := per_sample_suicide(bulk$gemma_id[i], gsm)]
  x$kind <- "gemma"; x$tech <- ifelse(grepl("GSE(80655|202537|224683|138082)", x$name), "seq", "array")
  x
}) })   # common.R's loaders read paths relative to fairness_sample_level/
names(G) <- bulk$dataset

message("loading LIBD")
LB <- lapply(c(DLPFC = "DLPFC", HIPPO = "HIPPO"), function(r) {
  x <- readRDS(file.path(GL, sprintf("libd_brainseq/analysis/data/libd_%s.rds", r)))
  cv <- copy(x$cv)
  cv[, suicide := fifelse(is.na(manner), NA_integer_, as.integer(manner == "Suicide"))]
  cv[, pH := NA_real_]
  cv[, PMI := as.numeric(PMI)]
  list(E = x$E, cv = cv, name = paste0("LIBD_", r), kind = "libd", tech = "seq")
})
names(LB) <- c("LIBD_DLPFC", "LIBD_HIPPO")

message("loading PsychAD whole tissue")
PA <- local({
  owd <- setwd(PSY); on.exit(setwd(owd))
  D <- load_donors()
  out <- list()
  for (b in c("HBCC", "MSSM")) {
    L <- load_level(b, "total")
    d <- cell_samples(L, D, paste0(b, "_main"), "all")
    d[, gsm := pb_id]; d[, PMI := pmi]; d[, pH := NA_real_]; d[, suicide := NA_integer_]
    out[[paste0("PsychAD_", b)]] <- list(M = L$M, genes = L$genes, cv = d, name = paste0("PsychAD_", b), kind = "psychad", tech = "sn")
  }
  out
})

U <- c(G, LB, PA)
message(length(U), " units: ", paste(names(U), collapse = ", "))

# ---------------------------------------------------------------- fits
fit_unit <- function(u, d) {
  if (sum(d$sz == 1) < MIN_GROUP || sum(d$sz == 0) < MIN_GROUP) return(NULL)
  if (u$kind == "gemma") return(fit_sz(u$E, d, covs = c("age", "male")))
  if (u$kind == "libd") {
    covs <- c("age", "male")[vapply(c("age", "male"), function(v) length(unique(d[[v]])) > 1, logical(1))]
    des <- model.matrix(as.formula(paste("~ sz +", paste(covs, collapse = " + "))), data = d)
    f <- eBayes(lmFit(u$E[, d$gsm, drop = FALSE], des), trend = TRUE)
    return(list(gene = rownames(u$E), t = f$t[, 2], logfc = f$coefficients[, 2], p = f$p.value[, 2],
                n_sz = sum(d$sz), n_ctl = sum(1 - d$sz)))
  }
  r <- fit_pb(u$M, u$genes, d)
  r
}

# ---------------------------------------------------------------- case mix per unit
bank_suicide <- c(stanley = 14 / 50, pitt = 27 / 99, mssm = 0, aus = 16 / 67, pritzker = 3 / 17)  # cohort_diagnostics.tsv suicide_share_SZ rows (aus: NSW 8/37 + Victoria 8/30)
premortem_banks <- c("mssm", "penn")
casemix_of <- function(nm, d) {
  s <- d[sz == 1]; c_ <- d[sz == 0]
  b <- bank_of[[nm]]; if (is.null(b)) b <- NA_character_
  ps <- if (any(!is.na(s$suicide))) mean(s$suicide, na.rm = TRUE) else NA_real_
  data.table(unit = nm, bank = b, n_sz = nrow(s), n_ctl = nrow(c_),
             age_sz = median(s$age), age_ctl = median(c_$age), age_gap = median(s$age) - median(c_$age),
             pH_sz = if (all(is.na(s$pH))) NA_real_ else median(s$pH, na.rm = TRUE),
             PMI_sz = if (all(is.na(s$PMI))) NA_real_ else median(as.numeric(s$PMI), na.rm = TRUE),
             suicide_sz = if (!is.na(ps)) ps else if (!is.na(b) && b %in% names(bank_suicide)) bank_suicide[[b]] else NA_real_,
             suicide_source = if (!is.na(ps)) "per sample" else if (!is.na(b) && b %in% names(bank_suicide)) "bank level (cohorts)" else "not recorded",
             premortem = as.integer(!is.na(b) && b %in% premortem_banks),
             neff = 1 / (1 / nrow(s) + 1 / nrow(c_)))
}

full_cv <- lapply(U, function(u) u$cv)
CMX <- rbindlist(lapply(names(U), function(nm) casemix_of(nm, full_cv[[nm]])))
CMX[, tech := vapply(unit, function(nm) U[[nm]]$tech, "")]
writeT(CMX, "results/casemix_units.tsv", c("case mix of each unit's SZ group in the samples fitted (full data)",
       "suicide_sz: share of SZ donors with suicide as manner of death; suicide_source says per sample or bank level"))

# ---------------------------------------------------------------- pairs and relations
rel <- fread(file.path(GL, "subject_relations/subject_relations.tsv"), skip = "dataset_a\t")
rel <- rel[dataset_a %in% names(U) & dataset_b %in% names(U)]
pairs <- rel[, .(a = dataset_a, b = dataset_b, relation = relation_group)]

agree_all <- function(R, P) {
  rbindlist(lapply(seq_len(nrow(P)), function(i) {
    ra <- R[[P$a[i]]]; rb <- R[[P$b[i]]]
    if (is.null(ra) || is.null(rb)) return(NULL)
    x <- agree(ra, rb, pi1 = FALSE)
    data.table(a = P$a[i], b = P$b[i], relation = P$relation[i], n_common = x$n_common, rho_t = x$rho_t, top100_fold = x$top100_fold)
  }))
}

message("full fits")
R0 <- mclapply(names(U), function(nm) fit_unit(U[[nm]], full_cv[[nm]]), mc.cores = NCORES)
names(R0) <- names(U)
A0 <- agree_all(R0, pairs)
A0 <- merge(A0, CMX[, .(a = unit, neff_a = neff, tech_a = tech)], by = "a")
A0 <- merge(A0, CMX[, .(b = unit, neff_b = neff, tech_b = tech)], by = "b")
setcolorder(A0, c("a", "b"))
writeT(A0, "results/casemix_pairs_full.tsv", c("agreement of the full-data SZ vs control fits, every unit pair with its relation_group"))

# test 1 (pairs) and power: scripts/casemix_pairs.R, from results/casemix_pairs_full.tsv and casemix_units.tsv

# ---------------------------------------------------------------- test 2: comparable subsets vs random subsets of the same size
CONDS <- list(
  age65 = function(d) d[age <= AGE_MAX],
  no_suicide = function(d) d[is.na(suicide) | suicide == 0],
  age65_no_suicide = function(d) d[age <= AGE_MAX & (is.na(suicide) | suicide == 0)])
primary <- pairs[relation == "subject-separate"]
key_pairs <- data.table(a = c("PsychAD_HBCC", "LIBD_DLPFC"), b = c("PsychAD_MSSM", "PsychAD_MSSM"))
AFFECTED <- character(0)   # units that lost >= 20% of SZ or control samples in the current condition (set per condition)
summ <- function(A) {
  p <- A[relation == "subject-separate"]
  pa <- p[a %in% AFFECTED | b %in% AFFECTED]
  kp <- merge(A, key_pairs, by = c("a", "b"))
  kp2 <- merge(A, key_pairs[, .(a = b, b = a)], by = c("a", "b"))
  kp <- rbind(kp, kp2)
  c(median_rho_t = median(p$rho_t), median_top100_fold = median(p$top100_fold), n_pairs = nrow(p),
    median_rho_t_affected = if (nrow(pa)) median(pa$rho_t) else NA_real_,
    median_top100_fold_affected = if (nrow(pa)) median(pa$top100_fold) else NA_real_, n_pairs_affected = nrow(pa),
    setNames(kp$rho_t, paste0("rho_t:", pmin(kp$a, kp$b), "|", pmax(kp$a, kp$b))))
}
T2 <- list(); T2n <- list()
for (cn in names(CONDS)) {
  message("condition ", cn)
  sub_cv <- lapply(full_cv, CONDS[[cn]])
  Rc <- mclapply(names(U), function(nm) fit_unit(U[[nm]], sub_cv[[nm]]), mc.cores = NCORES); names(Rc) <- names(U)
  kept <- names(U)[!vapply(Rc, is.null, TRUE)]
  AFFECTED <- kept[vapply(kept, function(nm) {
    f <- full_cv[[nm]]; s_ <- sub_cv[[nm]]
    sum(s_$sz == 1) < 0.8 * sum(f$sz == 1) || sum(s_$sz == 0) < 0.8 * sum(f$sz == 0) }, TRUE)]
  message("  affected units: ", paste(AFFECTED, collapse = ", "))
  full_c <- summ(A0[a %in% kept & b %in% kept])
  obs <- summ(agree_all(Rc, pairs[a %in% kept & b %in% kept]))
  counts <- lapply(kept, function(nm) c(sz = sum(sub_cv[[nm]]$sz == 1), ctl = sum(sub_cv[[nm]]$sz == 0))); names(counts) <- kept
  reps <- mclapply(seq_len(NREP), function(r) {
    set.seed(SEED + 1000 * match(cn, names(CONDS)) + r)
    Rr <- lapply(kept, function(nm) {
      d <- full_cv[[nm]]
      d <- rbind(d[sz == 1][sample.int(sum(d$sz == 1), counts[[nm]][["sz"]])], d[sz == 0][sample.int(sum(d$sz == 0), counts[[nm]][["ctl"]])])
      fit_unit(U[[nm]], d)
    }); names(Rr) <- kept
    summ(agree_all(Rr, pairs[a %in% kept & b %in% kept]))
  }, mc.cores = NCORES)
  reps <- do.call(rbind, lapply(reps, function(x) x[names(obs)]))
  T2[[cn]] <- rbindlist(lapply(names(obs), function(k) data.table(condition = cn, statistic = k, observed = obs[[k]], affected_units = paste(AFFECTED, collapse = ", "),
           random_median = median(reps[, k], na.rm = TRUE), random_q05 = quantile(reps[, k], 0.05, na.rm = TRUE),
           random_q95 = quantile(reps[, k], 0.95, na.rm = TRUE),
           p_greater = (1 + sum(reps[, k] >= obs[[k]], na.rm = TRUE)) / (1 + sum(!is.na(reps[, k]))), n_rep = sum(!is.na(reps[, k])),
           full_data = full_c[k])))
  T2n[[cn]] <- data.table(condition = cn, unit = kept, n_sz = vapply(counts, `[[`, 0, "sz"), n_ctl = vapply(counts, `[[`, 0, "ctl"),
                          n_sz_full = vapply(kept, function(nm) sum(full_cv[[nm]]$sz == 1), 0), n_ctl_full = vapply(kept, function(nm) sum(full_cv[[nm]]$sz == 0), 0))
}
T2 <- rbindlist(T2)
writeT(T2, "results/casemix_test2.tsv", c(sprintf("comparable subsets (age at death <= %d; suicides excluded where recorded; both) against %d random subsets with the same SZ and control counts per unit", AGE_MAX, NREP),
       "p_greater: share of random subsets at or above the observed value (comparable subsets raising agreement); full_data: the full-data value"))
writeT(rbindlist(T2n), "results/casemix_test2_counts.tsv", c("SZ and control samples per unit in each condition (units with fewer than 5 per group dropped)"))

meta <- list(script = SCRIPT, built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), seed = SEED, nrep = NREP,
             inputs = lapply(c(subject_relations = file.path(GL, "subject_relations/subject_relations.tsv"),
                               registry = file.path(GL, "controls_age_sex/data/datasets.tsv"),
                               cohorts = file.path(GL, "cohorts/cohort_diagnostics.tsv"),
                               libd_DLPFC = file.path(GL, "libd_brainseq/analysis/data/libd_DLPFC.rds"),
                               libd_HIPPO = file.path(GL, "libd_brainseq/analysis/data/libd_HIPPO.rds"),
                               psychad_donors = file.path(PSY, "donors.tsv")), sha256))
jsonlite::write_json(meta, "results/casemix_meta.json", auto_unbox = TRUE, pretty = TRUE)
message("done")
