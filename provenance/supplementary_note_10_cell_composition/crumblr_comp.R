# NEW 2026-10-18 (crumblr addendum). Analyses 1 and 2 of composition_ph_2026-10-18 redone with the real crumblr + variancePartition pipeline:
#   cobj <- crumblr(counts)  (counts: donors x cell types; pseudocount 0.5, tau = 1, max.ratio = 5, quant = 0.05: package defaults, as in the vignette)
#   fit <- eBayes(dream(cobj, ~ formula, info));  topTable(fit, coef, confint = TRUE)
# Model skeleton (Lee et al. 2026): ~ scale(age) + sex + effect. Fixed effects only (one observation per donor, no repeated measures, so there is no random term to fit).
# dream() with a fixed-effects-only formula returns the same coefficients and t as limma::lmFit with the crumblr weights (validation C); bootstrap and permutation loops therefore use lmFit.
# Point estimates, intervals and P values: dream + eBayes + topTable (moderated t, df.total, 95% CI from the moderated t). Usage: Rscript crumblr_comp.R <crumblr_dir> <input_dir>
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressMessages({library(crumblr); library(variancePartition); library(limma); library(data.table); library(parallel)})
args <- commandArgs(TRUE); CR <- args[1]; I <- args[2]; O <- file.path(CR, "out")
SEED <- 20261018; NBOOT <- 500; NPERM <- 10000; NC <- 6
sha <- function(f) substr(sub(" .*", "", system(paste("sha256sum", shQuote(f)), intern = TRUE)), 1, 16)
vers <- sprintf("crumblr %s; variancePartition %s; limma %s; R %s", packageVersion("crumblr"), packageVersion("variancePartition"), packageVersion("limma"), getRversion())
hdr <- c(sprintf("NEW 2026-10-18 by composition_ph_2026-10-18/crumblr/scripts/crumblr_comp.R; seed %d; %d donor bootstrap resamples (stratified by diagnosis); %d label permutations; %s", SEED, NBOOT, NPERM, vers),
  sprintf("inputs sha256(16): nsw_pb_cols.tsv %s; nsw_donor_covariates_used.tsv %s; batiuk_counts.tsv %s; batiuk_donors.tsv %s", sha(file.path(I, "nsw_pb_cols.tsv")), sha(file.path(I, "nsw_donor_covariates_used.tsv")), sha(file.path(I, "batiuk_counts.tsv")), sha(file.path(I, "batiuk_donors.tsv"))),
  "method: crumblr (CLR of counts + 0.5, Dirichlet-multinomial asymptotic precision weights, tau = 1, ratio cap 5 at the 5% quantile) then dream (fixed effects) + eBayes; est = logFC in CLR units; scale(age), scale(ph), scale(pmi) standardized within the fitted donors")
wr <- function(x, f, extra = NULL) { p <- file.path(O, f); writeLines(paste("#", c(hdr, extra)), p); fwrite(x, p, sep = "\t", append = TRUE, col.names = TRUE) }
WARN <- new.env(); WARN$log <- list()
note <- function(tag, msg) WARN$log[[length(WARN$log) + 1]] <- data.frame(call = tag, message = msg)

# ---------------- data (same selection as comp.R) ----------------
cv <- fread(file.path(I, "nsw_donor_covariates_used.tsv")); cv[, sex01 := as.integer(sex == "Male")]
cv <- cv[dx %in% c("Control", "Schizophrenia", "SchizoaffectiveDisorder") & is.finite(ph) & is.finite(age) & is.finite(pmi) & !is.na(sex01)]
cols <- fread(file.path(I, "nsw_pb_cols.tsv")); tot <- cols[level == "total", .(donor, n_tot = n_nuclei)]; cv <- merge(cv, tot, by = "donor")[n_tot >= 20]
cm <- function(lev) { w <- dcast(cols[level == lev & donor %in% cv$donor], donor ~ cell_type, value.var = "n_nuclei", fill = 0); m <- as.matrix(w[, -1]); rownames(m) <- w$donor; m[cv$donor, , drop = FALSE] }
nsw_dat <- as.data.frame(cv); rownames(nsw_dat) <- cv$donor
bd <- fread(file.path(I, "batiuk_donors.tsv")); bd[, `:=`(sex01 = as.integer(sex == "M"), ph = suppressWarnings(as.numeric(ph)), age = as.numeric(age), pmi = as.numeric(pmi))]
bc <- fread(file.path(I, "batiuk_counts.tsv"))
bm <- function(col) { w <- dcast(bc, donor ~ get(col), value.var = "n_nuclei", fill = 0, fun.aggregate = sum); m <- as.matrix(w[, -1]); rownames(m) <- w$donor; m[bd$donor, , drop = FALSE] }
bat_dat <- as.data.frame(bd); rownames(bat_dat) <- bd$donor
DS <- list(NSW_OFC = list(dat = nsw_dat, mats = list(class = cm("class"), subclass = cm("subclass"))),
           Batiuk2022 = list(dat = bat_dat, mats = list(group = bm("cls"), celltype = bm("cell_type"))))

# ---------------- crumblr + dream wrappers ----------------
mkform <- function(rhs) as.formula(paste("~", rhs))
design_ok <- function(d, rhs) { des <- model.matrix(mkform(rhs), d); if (qr(des)$rank < ncol(des) || nrow(des) - ncol(des) < 1) NULL else des }
lm_coef <- function(M, d, rhs, term) {      # fast path: lmFit with crumblr weights (identical to fixed-effects dream; validation C)
  des <- design_ok(d, rhs); if (is.null(des) || !term %in% colnames(des)) return(NULL)
  rownames(M) <- make.unique(rownames(M)); cobj <- crumblr(M)
  lmFit(cobj, des)$coefficients[, term]
}
dream_tab <- function(M, d, rhs, terms, tag) {   # point estimates, intervals, P: dream + eBayes + topTable(confint)
  des <- design_ok(d, rhs); if (is.null(des)) { note(tag, "design rank deficient or no residual df: not estimable"); return(NULL) }
  rownames(M) <- make.unique(rownames(M)); rownames(d) <- rownames(M)
  cobj <- crumblr(M)
  res <- withCallingHandlers(tryCatch({ f <- dream(cobj, mkform(rhs), d); eBayes(f) }, error = function(e) { note(tag, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { note(tag, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  if (is.null(res)) return(NULL)
  rbindlist(lapply(terms, function(tm) { tb <- topTable(res, coef = tm, number = Inf, confint = TRUE, sort.by = "none")
    data.table(category = rownames(tb), term = tm, est = tb$logFC, se = tb$logFC / tb$t, lo = tb$CI.L, hi = tb$CI.R, t = tb$t, p = tb$P.Value, df = res$df.total[1], df_residual = res$df.residual[1], prior_df = if (is.null(res$df.prior)) NA else res$df.prior[1]) }))
}
strat_draw <- function(dx) unlist(lapply(split(seq_along(dx), dx), function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]), use.names = FALSE)

res1 <- list(); res2 <- list(); resS <- list(); resC <- list()
# ---------------- analysis 1 ----------------
for (dn in names(DS)) {
  dat <- DS[[dn]]$dat; ctl <- dat[dat$dx == "Control" & is.finite(dat$ph), ]
  for (lev in names(DS[[dn]]$mats)) {
    M <- DS[[dn]]$mats[[lev]]; fr <- colMeans(M / rowSums(M)); sdp <- sd(ctl$ph)
    for (rhs in c("scale(age) + sex01 + scale(ph)", "scale(age) + sex01 + scale(ph) + scale(pmi)")) {
      tb <- dream_tab(M[rownames(ctl), , drop = FALSE], ctl, rhs, "scale(ph)", paste(dn, lev, "A1", rhs)); if (is.null(tb)) next
      tb[, `:=`(dataset = dn, level = lev, population = "controls", model = paste("crumblr ~", rhs), n = nrow(ctl), mean_fraction = fr[category], ph_sd = sdp, est_per_sd = est, lo_per_sd = lo, hi_per_sd = hi)]
      tb[, `:=`(est = est_per_sd / sdp, se = se / sdp, lo = lo_per_sd / sdp, hi = hi_per_sd / sdp)]   # est, se, lo, hi per 1 pH unit; _per_sd per SD
      res1[[length(res1) + 1]] <- tb
    }
  }
}
# ---------------- analysis 2 ----------------
set.seed(SEED)
for (dn in names(DS)) {
  dat <- DS[[dn]]$dat; dat$dx <- factor(dat$dx, levels = c("Control", setdiff(unique(dat$dx), "Control")))
  alld <- rownames(dat); phd <- rownames(dat)[is.finite(dat$ph)]
  dxterms <- setdiff(colnames(model.matrix(~ dx, dat)), "(Intercept)"); szterm <- "dxSchizophrenia"
  sets <- if (length(phd) == length(alld)) list(all = alld) else list(all = alld, phset = phd)
  mods <- list(M0 = "scale(age) + sex01 + dx", M1 = "scale(age) + sex01 + dx + scale(ph)", M2 = "scale(age) + sex01 + dx + scale(pmi)", M3 = "scale(age) + sex01 + dx + scale(ph) + scale(pmi)")
  for (lev in names(DS[[dn]]$mats)) {
    M <- DS[[dn]]$mats[[lev]]; fr <- colMeans(M / rowSums(M)); cats <- colnames(M); fitsB <- list()
    for (sn in names(sets)) for (mn in names(mods)) {
      if (sn == "all" && length(sets) > 1 && mn %in% c("M1", "M3")) next
      r <- sets[[sn]]; d <- dat[r, ]
      tb <- dream_tab(M[r, , drop = FALSE], d, mods[[mn]], dxterms, paste(dn, lev, sn, mn)); if (is.null(tb)) next
      tb[, `:=`(dataset = dn, level = lev, donor_set = sn, model = paste("crumblr ~", mods[[mn]]), model_id = mn, n = nrow(d), n_sz = sum(d$dx == "Schizophrenia"), n_control = sum(d$dx == "Control"), mean_fraction = fr[category])]
      res2[[length(res2) + 1]] <- tb; fitsB[[paste(sn, mn)]] <- tb[term == szterm, setNames(est, category)]
    }
    refset <- if (length(sets) > 1) "phset" else "all"; rd <- sets[[refset]]; dxr <- as.character(dat[rd, "dx"])
    pr <- list(c("M0", "M1"), c("M0", "M3"), c("M2", "M3"), c("M0", "M2")); slope <- function(b0, b1) sum(b0 * b1) / sum(b0 * b0)
    draws <- lapply(seq_len(NBOOT), function(b) strat_draw(dxr))
    bcoef <- function(ix, mn) lm_coef(M[rd, , drop = FALSE][ix, , drop = FALSE], dat[rd, ][ix, ], mods[[mn]], szterm)
    BT <- mclapply(draws, function(ix) { o <- list(); for (p in pr) { a <- tryCatch(bcoef(ix, p[1]), error = function(e) NULL); b <- tryCatch(bcoef(ix, p[2]), error = function(e) NULL); o[[paste(p, collapse = "->")]] <- if (is.null(a) || is.null(b)) NULL else rbind(a, b) }; o }, mc.cores = NC, mc.preschedule = FALSE)
    for (p in pr) {
      nm <- paste(p, collapse = "->"); b0 <- fitsB[[paste(refset, p[1])]]; b1 <- fitsB[[paste(refset, p[2])]]; if (is.null(b0) || is.null(b1)) next
      ok <- Filter(is.matrix, lapply(BT, function(x) if (is.list(x)) x[[nm]] else NULL))
      sh_b <- t(sapply(ok, function(m) 1 - m[2, ] / m[1, ])); sl_b <- sapply(ok, function(m) 1 - slope(m[1, ], m[2, ]))
      rows <- data.table(dataset = dn, level = lev, donor_set = refset, comparison = nm, category = cats, mean_fraction = fr[cats], b_before = b0[cats], b_after = b1[cats], share_removed = 1 - b1[cats] / b0[cats],
        share_lo = apply(sh_b, 2, quantile, 0.025, na.rm = TRUE), share_hi = apply(sh_b, 2, quantile, 0.975, na.rm = TRUE), n_boot_ok = length(ok))
      rows <- rbind(rows, data.table(dataset = dn, level = lev, donor_set = refset, comparison = nm, category = "ALL_TYPES_1_minus_slope", mean_fraction = NA, b_before = NA, b_after = NA, share_removed = 1 - slope(b0[cats], b1[cats]),
        share_lo = quantile(sl_b, 0.025, na.rm = TRUE), share_hi = quantile(sl_b, 0.975, na.rm = TRUE), n_boot_ok = length(ok)))
      resS[[length(resS) + 1]] <- rows
    }
    # ---- across-type correlation ----
    cor_fun <- function(a, b) c(pearson = cor(a, b), spearman = cor(a, b, method = "spearman"))
    Ma <- M[alld, , drop = FALSE]; da <- dat[alld, ]
    pheff <- function(ix, variant) { d <- da[ix, ]; Mx <- Ma[ix, , drop = FALSE]; hp <- is.finite(d$ph)
      s <- if (variant == "controls") which(d$dx == "Control" & hp) else which(hp); if (length(s) < 5) return(NULL)
      lm_coef(Mx[s, , drop = FALSE], d[s, ], if (variant == "controls") "scale(age) + sex01 + scale(ph)" else "scale(age) + sex01 + dx + scale(ph)", "scale(ph)") }
    sz0 <- function(ix) lm_coef(Ma[ix, , drop = FALSE], da[ix, ], mods$M0, szterm)
    szph <- function(ix) { d <- da[ix, ]; hp <- which(is.finite(d$ph)); lm_coef(Ma[ix, , drop = FALSE][hp, , drop = FALSE], d[hp, ], mods$M1, szterm) }
    full <- seq_along(alld); dxa <- as.character(da$dx); cdraws <- lapply(seq_len(NBOOT), function(b) strat_draw(dxa))
    for (variant in c("controls", "dx_adjusted_all_with_ph")) {
      p0 <- tryCatch(pheff(full, variant), error = function(e) NULL); if (is.null(p0)) next
      for (szv in c("sz_unadjusted_all", "sz_after_ph")) {
        s0 <- tryCatch(if (szv == "sz_unadjusted_all") sz0(full) else szph(full), error = function(e) NULL); if (is.null(s0)) next
        r0 <- cor_fun(p0, s0)
        bb <- mclapply(cdraws, function(ix) tryCatch({ p <- pheff(ix, variant); s <- if (szv == "sz_unadjusted_all") sz0(ix) else szph(ix); if (is.null(p) || is.null(s)) NULL else cor_fun(p, s) }, error = function(e) NULL), mc.cores = NC, mc.preschedule = FALSE)
        bb <- do.call(rbind, Filter(is.numeric, bb)); bb <- bb[is.finite(rowSums(bb)), , drop = FALSE]; if (nrow(bb) < 20) next
        set.seed(SEED + 1); perm <- t(replicate(NPERM, cor_fun(p0, sample(s0))))
        for (cn in c("pearson", "spearman")) { z <- atanh(pmin(pmax(bb[, cn], -0.999999), 0.999999)); sdz <- sd(z)
          resC[[length(resC) + 1]] <- data.table(dataset = dn, level = lev, ph_effect = variant, sz_effect = szv, method = cn, n_types = length(cats), r = unname(r0[cn]),
            boot_pct_lo = unname(quantile(bb[, cn], .025)), boot_pct_hi = unname(quantile(bb[, cn], .975)), boot_fisher_lo = tanh(atanh(r0[cn]) - 1.96 * sdz), boot_fisher_hi = tanh(atanh(r0[cn]) + 1.96 * sdz),
            n_boot_ok = nrow(bb), perm_p_two_sided = (1 + sum(abs(perm[, cn]) >= abs(r0[cn]))) / (NPERM + 1), perm_null_sd = sd(perm[, cn])) }
      }
    }
    cat("done", dn, lev, "\n"); flush.console()
  }
}
wr(rbindlist(res1, fill = TRUE), "comp_ph_effect_controls.tsv", "analysis 1 (crumblr): est = change in CLR per 1 pH unit; est_per_sd per SD of pH in the fitted donors; interval/P from dream + eBayes + topTable")
wr(rbindlist(res2, fill = TRUE), "comp_dx_effect.tsv", "analysis 2 (crumblr): est = difference in CLR vs Control; donor_set phset = donors with pH (Batiuk)")
wr(rbindlist(resS, fill = TRUE), "comp_dx_share_removed.tsv", "share_removed = 1 - b_after/b_before on the same donors (dream point estimates); interval = percentile over donor bootstrap (crumblr + lmFit refit in each resample)")
wr(rbindlist(resC, fill = TRUE), "comp_corr_ph_vs_sz.tsv", "correlation across cell types of per-type pH effect and per-type SZ effect (crumblr)")
wl <- if (length(WARN$log)) rbindlist(WARN$log) else data.table(call = character(), message = character()); wl <- unique(wl); wr(wl, "dream_messages.tsv", "warnings and errors raised by dream/eBayes (unique)")
cat("all done\n")
