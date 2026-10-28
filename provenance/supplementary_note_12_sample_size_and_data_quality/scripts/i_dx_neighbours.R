# I. When a donor's nearest profile in the other dataset is not their own, is it someone with the same diagnosis?
# Paul, 2026-09-28: "if the first choice isn't the right person, is there a tendency for them to share the diagnosis? Or
# is that random?"
#
# Run from analyses/sample_level_reanalysis_of_datasets/, after g_individuals.R and h_stanley_individuals.R:
#   Rscript same_donors_pitt/i_dx_neighbours.R > same_donors_pitt/results/i_dx_neighbours.log 2>&1
#
# Sources (SZ and controls only; the similarity matrix is identify()'s: Pearson over genes of residual profiles z-scored
# within the candidates, row = dataset a donor, column = dataset b donor, true partner on the diagonal):
#   LIBD regions   DLPFC x hippocampus, recomputed here as ../libd_brainseq/analysis/scripts/e_individuals_regions.R
#                  computes it (its functions and libd_common.R sourced; that script is not run or edited).
#   PsychAD x LIBD whole tissue ('all'), recomputed here as ../psychad_cellxgene/scripts/e_individuals.R main() computes it
#                  (functions sourced; not run or edited). Cell classes not done.
#   Pittsburgh     results/g_individuals_C.rds (g_individuals.R), all 28 unit pairs.
#   Stanley        results/h_stanley_C.rds (h_stanley_individuals.R), subject set 'sz_ctrl'; pairs where one diagnosis has
#                  fewer than 2 candidates are skipped (chen has no SZ subjects).
# Residual modes: 'nodx' (covariates only: the diagnosis signal stays in) and 'dx' (diagnosis also regressed out), where the
# source's script has them: LIBD regions, PsychAD x LIBD and Pittsburgh have both; h_stanley_individuals.R has no 'dx'
# residuals, so Stanley is 'nodx' only. The 'dx' rows are not a valid test here: removing each diagnosis group's mean over
# the same donors makes same-group donors slightly anti-correlated (about -1/(n-1)); their rows carry a note saying so.
# Per donor and direction (a -> b: row i against every b donor; b -> a: column i), the true partner excluded:
#   nearest-wrong same diagnosis: whether the most similar non-self candidate has the donor's diagnosis; expected = share of
#     same-diagnosis donors among the donor's non-self candidates. Summarized over all donors and over donors whose true
#     partner is not top-1 (rank from identify_ranks(), unchanged by the diagnosis labels).
#   diagnosis AUROC: over the non-self candidates, P(similarity to a same-diagnosis candidate > to a different-diagnosis one)
#     (ties half); mean over donors.
# Null: diagnosis labels shuffled among the relation's donors within sex (strat_perm over sex), NPERM times, similarity
# matrix fixed; p = (1 + #null >= observed) / (NPERM + 1). Pooled per source: donor-direction entries summed over the
# source's relations, with each relation's k-th shuffle summed into the pooled k-th shuffle; the relations share donors,
# so pooled values are not independent across relations and the pooled p treats them as if they were.
# Outputs (same_donors_pitt/results/): i_dx_neighbours_summary.tsv, i_dx_neighbours_donors.tsv.gz, i_dx_neighbours_meta.json
source("scripts/common.R")
suppressMessages(library(jsonlite))
PE <- new.env()
local({ owd <- setwd("../psychad_cellxgene"); on.exit(setwd(owd)); sys.source("scripts/e_individuals.R", envir = PE) })
LE <- new.env()
local({ owd <- setwd("../libd_brainseq/analysis"); on.exit(setwd(owd)); sys.source("scripts/libd_common.R", envir = LE) })
source("scripts/common.R")   # libd_common.R and pb_common.R re-source it globally; this keeps its definitions in force here
for (f in c("resid_mat", "identify", "strat_perm", "identify_ranks")) assign(f, get(f, envir = PE))
RES <- "same_donors_pitt/results"
SCRIPT <- "fairness_sample_level/same_donors_pitt/i_dx_neighbours.R"
NPERM <- 1000
SEED_I <- 20260928
INS <- list()
note_in <- function(f) { INS[[f]] <<- sha256(f); f }

# ---------------------------------------------------------------- similarity matrices
REL <- list()   # each: source, relation, mode, C, dx, sex, id
libd_regions <- function() {
  owd <- setwd("../libd_brainseq/analysis"); on.exit(setwd(owd))
  for (f in c("data/libd_DLPFC.rds", "data/libd_HIPPO.rds")) INS[[file.path("../libd_brainseq/analysis", f)]] <<- sha256(f)
  e <- new.env(parent = LE)
  evalq({
    D <- list(DLPFC = load_region("DLPFC"), HIPPO = load_region("HIPPO"))
    for (r in names(D)) { D[[r]]$cv <- as.data.table(D[[r]]$cv); D[[r]]$cv[, mito := as.numeric(mito)] }
    B <- both_region_donors(D)
    gD <- D$DLPFC$cv[match(B$BrNum, BrNum), gsm]; gH <- D$HIPPO$cv[match(B$BrNum, BrNum), gsm]
    genes <- intersect(rownames(D$DLPFC$E), rownames(D$HIPPO$E))
  }, e)
  lapply(c(nodx = FALSE, dx = TRUE), function(dx) {
    R <- lapply(e$D, function(x) resid_mat(x$E, x$cv, get("BASE_COVS", envir = LE), dx = dx))
    idn <- identify(R$DLPFC[e$genes, e$gD], R$HIPPO[e$genes, e$gH])
    list(C = idn$C, rank = idn$rank, id = e$B$BrNum, dx = ifelse(e$B$sz == 1, "SZ", "CTRL"), sex = ifelse(e$B$male == 1, "M", "F"))
  })
}
psychad_libd <- function() {
  owd <- setwd("../psychad_cellxgene"); on.exit(setwd(owd))
  for (f in c("results/fingerprint_matches.tsv", "donors.tsv", "pseudobulk/HBCC_Cohort/counts_total.npz", "../libd_brainseq/analysis/data/libd_DLPFC.rds"))
    INS[[file.path("../psychad_cellxgene", f)]] <<- sha256(f)
  e <- new.env(parent = PE)
  evalq({
    P <- read_pairs_fp(); pr <- P$unamb
    Dn <- load_donors()
    lib <- readRDS("../libd_brainseq/analysis/data/libd_DLPFC.rds")
    lcv <- as.data.table(lib$cv); lcv[, mito := as.numeric(mito)]
    pr <- pr[br %in% lcv$BrNum]
    lg <- lcv[match(pr$br, BrNum), gsm]
    H <- hbcc_logcpm(load_level("HBCC", "total"), Dn, "all")
    ok <- pr$hbcc %in% colnames(H$E)
    genes <- intersect(rownames(H$E), rownames(lib$E))
  }, e)
  lapply(c(nodx = FALSE, dx = TRUE), function(dx) {
    HR <- resid_mat(e$H$E, e$H$cv, get("HBCC_COVS", envir = PE), dx = dx)
    LR <- resid_mat(e$lib$E, e$lcv, get("LIBD_COVS", envir = PE), dx = dx)
    idn <- identify(HR[e$genes, e$pr$hbcc[e$ok]], LR[e$genes, e$lg[e$ok]])
    list(C = idn$C, rank = idn$rank, id = e$pr$hbcc[e$ok], dx = e$pr$dx[e$ok], sex = e$pr$sex[e$ok])
  })
}
# checks against the scripts' own identification output (nodx ranks)
lr <- libd_regions()
chk <- fread(cmd = paste("grep -v '^#'", shQuote(note_in("../libd_brainseq/analysis/results/e_individuals_regions_identify.tsv"))))
stopifnot(identical(as.numeric(chk$rank_true), as.numeric(lr$nodx$rank)))
for (m in names(lr)) REL[[length(REL) + 1]] <- c(list(source = "LIBD regions", relation = "LIBD DLPFC x hippocampus", mode = m), lr[[m]][c("C", "id", "dx", "sex")])
pl <- psychad_libd()
chk <- fread(cmd = paste("grep -v '^#'", shQuote(note_in("../psychad_cellxgene/results/e_individuals_identify.tsv"))))[level == "all" & mode == "nodx"]
stopifnot(identical(as.numeric(chk$rank_true), as.numeric(pl$nodx$rank)))
for (m in names(pl)) REL[[length(REL) + 1]] <- c(list(source = "PsychAD x LIBD", relation = "HBCC snRNA-seq pseudobulk x LIBD DLPFC", mode = m), pl[[m]][c("C", "id", "dx", "sex")])
for (x in readRDS(note_in(file.path(RES, "g_individuals_C.rds"))))
  REL[[length(REL) + 1]] <- list(source = "Pittsburgh", relation = paste(x$unit_a, "x", x$unit_b), mode = x$mode, C = x$C, id = x$id, dx = x$dx, sex = x$sex)
SKIP <- character(0)
for (x in readRDS(note_in(file.path(RES, "h_stanley_C.rds")))) if (x$subjects == "sz_ctrl") {
  rel <- paste(x$study_a, "x", x$study_b)
  if (min(sum(x$dx == "SZ"), sum(x$dx == "CTRL")) < 2) { SKIP <- c(SKIP, sprintf("%s (%d SZ, %d CTRL)", rel, sum(x$dx == "SZ"), sum(x$dx == "CTRL"))); next }
  REL[[length(REL) + 1]] <- list(source = paste("Stanley", sub("Neuropathology ", "", x$collection)), relation = rel, mode = "nodx", C = x$C, id = x$id, dx = x$dx, sex = x$sex)
}
stopifnot(all(vapply(REL, function(r) all(r$dx %in% c("SZ", "CTRL")), TRUE)))

# ---------------------------------------------------------------- statistics
auroc_row <- function(x, s) {   # x: similarities, s: TRUE same diagnosis; NA where either group is empty
  a <- x[s]; b <- x[!s]
  if (!length(a) || !length(b)) return(NA_real_)
  rk <- rank(c(a, b))
  (sum(rk[seq_along(a)]) - length(a) * (length(a) + 1) / 2) / (length(a) * length(b))
}
per_donor <- function(M, d) {   # M: focal donors x candidates, self on the diagonal
  n <- nrow(M); diag(M) <- -Inf
  nn <- max.col(M, ties.method = "first")
  S <- outer(d, d, "==")
  au <- vapply(seq_len(n), function(i) auroc_row(M[i, -i], S[i, -i]), 0)
  list(same_nn = d[nn] == d, exp_same = (rowSums(S) - 1) / (n - 1), auroc = au, nn = nn)
}
agg <- function(p, idx) c(n = length(idx), same = sum(p$same_nn[idx]), exp = sum(p$exp_same[idx]),
                         au_sum = sum(p$auroc[idx], na.rm = TRUE), au_n = sum(!is.na(p$auroc[idx])))
DON <- list(); SUMS <- list(); NULLS <- list()
set.seed(SEED_I)
for (k in seq_along(REL)) {
  r <- REL[[k]]; n <- nrow(r$C)
  rk <- identify_ranks(r$C)
  for (dir in c("ab", "ba")) {
    M <- if (dir == "ab") r$C else t(r$C)
    nt <- which(rk[[dir]] > 1)
    p <- per_donor(M, r$dx)
    obs <- rbind(all = agg(p, seq_len(n)), not_top1 = agg(p, nt))
    set.seed(SEED_I + k * 2 + (dir == "ba"))
    nul <- replicate(NPERM, { q <- per_donor(M, r$dx[strat_perm(r$sex)]); rbind(all = agg(q, seq_len(n)), not_top1 = agg(q, nt)) }, simplify = "array")
    key <- paste(r$source, r$relation, r$mode, dir, sep = "|")
    SUMS[[key]] <- list(source = r$source, relation = r$relation, mode = r$mode, direction = dir, obs = obs, n_cand = n,
                        n_SZ = sum(r$dx == "SZ"), n_CTRL = sum(r$dx == "CTRL"))
    NULLS[[key]] <- nul
    DON[[length(DON) + 1]] <- data.table(source = r$source, relation = r$relation, mode = r$mode, direction = dir, donor = r$id, dx = r$dx, sex = r$sex,
                                         rank_true = rk[[dir]], nearest_wrong = r$id[p$nn], nearest_wrong_dx = r$dx[p$nn], same_dx = p$same_nn,
                                         expected_same = round(p$exp_same, 4), dx_auroc = round(p$auroc, 4))
  }
  cat(sprintf("%s | %s | %s done\n", r$source, r$relation, r$mode))
}
DON <- rbindlist(DON)

row_of <- function(src, rel, mode, dir, subset, obs, nul, n_rel, n_cand = NA, n_SZ = NA, n_CTRL = NA) {
  o <- obs[subset, ]; nl <- nul[subset, , ]
  fr <- if (o["n"] > 0) o["same"] / o["n"] else NA; ex <- if (o["n"] > 0) o["exp"] / o["n"] else NA
  au <- if (o["au_n"] > 0) o["au_sum"] / o["au_n"] else NA
  nfr <- ifelse(nl["n", ] > 0, nl["same", ] / nl["n", ], NA); nau <- ifelse(nl["au_n", ] > 0, nl["au_sum", ] / nl["au_n", ], NA)
  data.table(source = src, relation = rel, mode = mode, direction = dir, donors = subset, n_relations = n_rel, n_candidates = n_cand, n_SZ = n_SZ, n_CTRL = n_CTRL,
             n_donor_directions = unname(o["n"]), obs_same_frac = round(unname(fr), 4), expected_same_frac = round(unname(ex), 4),
             null_same_frac = round(mean(nfr, na.rm = TRUE), 4), p_same = if (is.na(fr)) NA else signif((1 + sum(nfr >= fr, na.rm = TRUE)) / (NPERM + 1), 3),
             dx_auroc = round(unname(au), 4), null_dx_auroc = round(mean(nau, na.rm = TRUE), 4),
             p_auroc = if (is.na(au)) NA else signif((1 + sum(nau >= au, na.rm = TRUE)) / (NPERM + 1), 3))
}
OUT <- list()
for (key in names(SUMS)) { s <- SUMS[[key]]; for (sub in c("all", "not_top1"))
  OUT[[length(OUT) + 1]] <- row_of(s$source, s$relation, s$mode, s$direction, sub, s$obs, NULLS[[key]], 1, s$n_cand, s$n_SZ, s$n_CTRL) }
# pooled per source, mode and direction
grp <- unique(rbindlist(lapply(SUMS, function(s) data.table(source = s$source, mode = s$mode, direction = s$direction))))
for (g in seq_len(nrow(grp))) {
  ks <- names(SUMS)[vapply(SUMS, function(s) s$source == grp$source[g] && s$mode == grp$mode[g] && s$direction == grp$direction[g], TRUE)]
  obs <- Reduce(`+`, lapply(ks, function(k) SUMS[[k]]$obs)); nul <- Reduce(`+`, NULLS[ks])
  for (sub in c("all", "not_top1")) OUT[[length(OUT) + 1]] <- row_of(grp$source[g], "POOLED", grp$mode[g], grp$direction[g], sub, obs, nul, length(ks))
}
OUT <- rbindlist(OUT)
# 'dx' residuals remove each diagnosis group's mean over the same donors, which makes members of one group slightly
# anti-correlated (about -1/(n-1)): same-diagnosis similarity is pushed below chance by construction
OUT[, note := ifelse(mode == "dx", "dx residuals: biased toward different-diagnosis neighbours by construction (group means removed over the same donors); not a valid test of this question", "")]
hdr <- function(f, lines) { writeLines(paste("#", lines), f); f }
stamp <- sprintf("built %s by %s; SZ and controls only; NPERM %d; seed %d", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT, NPERM, SEED_I)
f1 <- hdr(file.path(RES, "i_dx_neighbours_summary.tsv"), c(stamp,
  "per source x relation x residual mode x direction x donors (all; not_top1 = true partner not ranked first); relation POOLED = all of the source's relations summed (relations share donors: not independent)",
  "obs_same_frac: share of donors whose nearest non-self candidate has their diagnosis; expected_same_frac: mean share of same-diagnosis donors among each donor's non-self candidates; dx_auroc: mean over donors of P(similarity to same-diagnosis > different-diagnosis candidate)",
  "null: diagnosis labels shuffled within sex among the relation's donors, similarity fixed; p = (1 + #null >= observed) / (NPERM + 1)",
  if (length(SKIP)) paste("Stanley pairs skipped (a diagnosis with < 2 candidates):", paste(SKIP, collapse = "; ")) else "no pairs skipped"))
fwrite(OUT, f1, sep = "\t", append = TRUE, col.names = TRUE)
f2 <- file.path(RES, "i_dx_neighbours_donors.tsv")
hdr(f2, c(stamp, "per donor and direction: nearest non-self candidate (nearest_wrong), its diagnosis, whether it matches; expected_same; per-donor diagnosis AUROC"))
fwrite(DON, f2, sep = "\t", append = TRUE, col.names = TRUE)
system(paste("gzip -f", shQuote(f2))); f2 <- paste0(f2, ".gz")
write_json(list(script = SCRIPT, built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), R = R.version.string, nperm = NPERM, seed = SEED_I,
                modes = list(`LIBD regions` = c("nodx", "dx"), `PsychAD x LIBD` = c("nodx", "dx"), Pittsburgh = c("nodx", "dx"), Stanley = "nodx (h_stanley_individuals.R has no dx residuals)"),
                stanley_skipped = SKIP, functions_from = c("../psychad_cellxgene/scripts/e_individuals.R", "../libd_brainseq/analysis/scripts/libd_common.R"),
                functions_from_sha256 = list(e_individuals.R = sha256("../psychad_cellxgene/scripts/e_individuals.R"), libd_common.R = sha256("../libd_brainseq/analysis/scripts/libd_common.R")),
                inputs_sha256 = INS, outputs_sha256 = list(i_dx_neighbours_summary.tsv = sha256(f1), i_dx_neighbours_donors.tsv.gz = sha256(f2))),
           file.path(RES, "i_dx_neighbours_meta.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
