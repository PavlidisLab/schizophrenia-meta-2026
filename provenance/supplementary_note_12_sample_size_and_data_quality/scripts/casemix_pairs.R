# Test 1 of the case-mix analysis (pairs) and its power, from the tables scripts/casemix.R writes.
# Run from analyses/case_mix_agreement_test/ after casemix.R:  Rscript scripts/casemix_pairs.R   (a few minutes)
#
# Inputs: results/casemix_pairs_full.tsv (agreement of every unit pair, full-data fits), results/casemix_units.tsv (case mix
# of each unit's SZ group).
# Test: agreement (rho_t, top100_fold) residualized on log effective n of both units and on same technology; Spearman of
# the residual with the absolute case-mix difference of the pair. Pairs share units, so the null permutes the case-mix
# values across units (NPERM), not pairs. p_less: one-sided (larger difference, lower agreement); p_two: two-sided around
# the null mean. Sets: subject-separate bulk pairs; subject-separate pairs with PsychAD; the same without PsychAD Mount
# Sinai (the one unit with pre-mortem assessment and the oldest SZ group).
# Power: an age assignment drawn at random across units, rho_t residuals lowered by beta x |delta age| under that
# assignment, tested as above (so the residuals keep their observed dependence and the planted effect is the only link).
suppressMessages({library(data.table)})
SCRIPT <- "scripts/casemix_pairs.R"
SEED <- 20260926
NPERM <- as.integer(Sys.getenv("NPERM", "5000"))
NSIM <- as.integer(Sys.getenv("NSIM", "500"))
NPERM_SIM <- 500
writeT <- function(x, path, header = character(0)) {
  writeLines(paste("#", c(sprintf("built %s by casemix/%s", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT), header)), path)
  fwrite(x, path, sep = "\t", append = TRUE, col.names = TRUE)
}
rd <- function(f) fread(cmd = paste("grep -v '^#'", shQuote(f)))
A0 <- rd("results/casemix_pairs_full.tsv")
CMX <- rd("results/casemix_units.tsv")

residualize <- function(A, measure) {
  lneff <- log(A$neff_a) + log(A$neff_b); same_tech <- as.integer(A$tech_a == A$tech_b)
  resid(lm(A[[measure]] ~ lneff + same_tech))
}
DESC <- c("age_sz", "age_gap", "pH_sz", "PMI_sz", "suicide_sz", "premortem")
pair_test <- function(A, measure, set_name) {
  A <- copy(A); A[, resid := residualize(A, measure)]
  set.seed(SEED)
  rbindlist(lapply(DESC, function(v) {
    x <- setNames(CMX[[v]], CMX$unit)
    ok <- !is.na(x[A$a]) & !is.na(x[A$b])
    Ak <- A[ok]; units_k <- unique(c(Ak$a, Ak$b)); xk <- x[units_k]
    out <- data.table(set = set_name, measure = measure, descriptor = v, n_pairs = nrow(Ak), n_units = length(units_k),
                      rho = NA_real_, p_less = NA_real_, p_two = NA_real_)
    if (nrow(Ak) < 8 || length(unique(xk)) < 2) return(out)
    stat <- function(xx) suppressWarnings(cor(Ak$resid, abs(xx[Ak$a] - xx[Ak$b]), method = "spearman"))
    obs <- stat(xk)
    nul <- replicate(NPERM, stat(setNames(sample(xk), names(xk))))
    nul <- nul[!is.na(nul)]
    out[, `:=`(rho = obs, p_less = (1 + sum(nul <= obs)) / (length(nul) + 1),
               p_two = (1 + sum(abs(nul - mean(nul)) >= abs(obs - mean(nul)))) / (length(nul) + 1))]
    out
  }))
}
bulk_units <- CMX[tech != "sn", unit]
ss <- A0[relation == "subject-separate"]
SETS <- list("subject-separate, bulk" = ss[a %in% bulk_units & b %in% bulk_units],
             "subject-separate, with PsychAD" = ss,
             "subject-separate, with PsychAD, without Mount Sinai" = ss[a != "PsychAD_MSSM" & b != "PsychAD_MSSM"])
T1 <- rbindlist(lapply(names(SETS), function(s) rbind(pair_test(SETS[[s]], "rho_t", s), pair_test(SETS[[s]], "top100_fold", s))))
writeT(T1, "results/casemix_test1.tsv", c(sprintf("agreement residualized on log neff (both units) and same technology; Spearman with |case-mix difference|; null permutes case-mix values across units (%d permutations)", NPERM),
       "p_less: one-sided, larger difference with lower agreement; p_two: two-sided around the null mean",
       "premortem is 1 for PsychAD Mount Sinai only, so its test is the Mount Sinai pairs against the rest (smallest attainable p = 1/n_units)"))

power_one <- function(A, beta, nsim, seed) {
  set.seed(seed)
  res <- residualize(A, "rho_t"); units_k <- unique(c(A$a, A$b))
  ages <- CMX$age_sz[match(units_k, CMX$unit)]
  pv <- replicate(nsim, {
    xs <- setNames(sample(ages), units_k)
    d <- abs(xs[A$a] - xs[A$b])
    y <- res - beta * d
    obs <- cor(y, d, method = "spearman")
    nul <- replicate(NPERM_SIM, { xp <- setNames(sample(ages), units_k); cor(y, abs(xp[A$a] - xp[A$b]), method = "spearman") })
    (1 + sum(nul <= obs)) / (NPERM_SIM + 1)
  })
  c(power = mean(pv <= 0.05), median_d = median(abs(CMX$age_sz[match(A$a, CMX$unit)] - CMX$age_sz[match(A$b, CMX$unit)])))
}
betas <- c(0, 0.0025, 0.005, 0.01, 0.02)
PW <- rbindlist(lapply(c("subject-separate, bulk", "subject-separate, with PsychAD"), function(s) rbindlist(lapply(betas, function(b) {
  r <- power_one(SETS[[s]], b, NSIM, SEED + round(b * 1e4))
  data.table(set = s, beta_rho_per_year = b, rho_lost_per_10y = 10 * b, power_p05 = r[["power"]], n_sim = NSIM,
             median_abs_delta_age_observed = r[["median_d"]], n_pairs = nrow(SETS[[s]]), sd_resid_rho_t = sd(residualize(SETS[[s]], "rho_t")))
}))))
writeT(PW, "results/casemix_power.tsv", c("planted effect: rho_t residuals lowered by beta x |delta SZ median age| under a random assignment of the units' ages; one-sided unit-permutation test as in test 1",
       sprintf("%d simulations per row, %d permutations each; beta = 0 row is the null (share below 0.05 should be about 0.05)", NSIM, NPERM_SIM)))
sha256 <- function(f) sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE))
meta <- jsonlite::read_json("results/casemix_meta.json")
meta[c("nperm", "nsim")] <- NULL   # written by earlier versions of casemix.R, which ran test 1 itself
meta$pairs <- list(script = SCRIPT, built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), seed = SEED, nperm = NPERM, nsim = NSIM, nperm_sim = NPERM_SIM,
                   inputs = lapply(c(pairs_full = "results/casemix_pairs_full.tsv", units = "results/casemix_units.tsv"), sha256))
jsonlite::write_json(meta, "results/casemix_meta.json", auto_unbox = TRUE, pretty = TRUE)
message("done")
