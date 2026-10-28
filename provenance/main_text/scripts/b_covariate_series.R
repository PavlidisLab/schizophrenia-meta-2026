# B. Covariate adjustment series: SZ vs control per dataset under increasing adjustment, then cross-dataset agreement
# per level.
#
# Run from fairness_sample_level/:  Rscript scripts/b_covariate_series.R
# Inputs: see common.R; ../subject_relations/subject_relations.tsv (read_pairs()).
# Levels (limma; each level adds to the previous one what the dataset records for at least 80% of its samples):
#   L0  ~ sz
#   L1  + age + sex
#   L2  + PMI
#   L3  + pH, RIN, agonal factor (whichever are recorded)
#   L4  + medication: GSE21138 chlorpromazine-equivalent dose (controls 0), GSE138082 antipsychotic ON/OFF (Gemma 'drug'
#       characteristic; controls OFF); no other dataset records medication, so L4 = L3 elsewhere
#   L5a L1 + surrogate variables (sva; mod ~ sz + age + sex, mod0 ~ age + sex; number by num.sv method 'be')
#   L5b L4 + surrogate variables (sva; mod = L4 design, mod0 = L4 design without sz; number by num.sv 'be')
# Samples: complete cases for every covariate used at any level of that dataset, so all levels of a dataset use the same
# samples. In GSE21138 that drops SZ donors with no recorded dose.
# Outputs:
#   results/b_de/<dataset>.tsv.gz        per gene: logfc_, t_, p_ for each level
#   results/b_levels.tsv                 per dataset: samples, covariates per level, number of SVs, genes at q < 0.05 per level
#   results/b_pairs.tsv                  per level x dataset pair: agreement metrics (common.R agree())
#   results/b_independent_set.tsv        bulk datasets that are pairwise donor-separate (greedy by effective n)
#   results/b_recurrent_genes.tsv        per recurrent gene x level: Stouffer z over the independent set, signed by the
#                                        recurrent SZ direction, and the number of datasets agreeing in sign
#   results/b_gene_consistency.tsv       per gene x level: Stouffer z (unsigned) over the independent set
suppressMessages(library(sva))
source("scripts/common.R")
set.seed(20260925)
NCORES <- as.integer(Sys.getenv("NCORES", "4"))
LEVELS <- c("L0", "L1", "L2", "L3", "L4", "L5a", "L5b")
dir.create("results/b_de", showWarnings = FALSE)
L <- load_all()
L <- L[vapply(L, function(x) min(sum(x$cv$sz), sum(1 - x$cv$sz)) >= 5, logical(1))]

series <- function(x) {
  cv <- copy(x$cv); E <- x$E
  cover <- function(v) mean(!is.na(cv[[v]])) >= 0.8 && length(unique(na.omit(cv[[v]]))) > 1
  l2 <- if (cover("PMI")) "PMI" else character(0)
  l3 <- Filter(cover, c("pH", "RIN", "agonal"))
  l4 <- character(0)
  if (sum(!is.na(cv$cpz_mg[cv$sz == 1])) >= 5) { cv[sz == 0, cpz_mg := 0]; l4 <- "cpz_mg" }
  if (sum(!is.na(cv$on_drug)) >= 0.8 * nrow(cv) && length(unique(na.omit(cv$on_drug[cv$sz == 1]))) > 1) l4 <- "on_drug"
  used <- c("age", "male", l2, l3, l4)
  cv <- cv[complete.cases(cv[, ..used])]
  covs <- list(L0 = character(0), L1 = c("age", "male"), L2 = c("age", "male", l2), L3 = c("age", "male", l2, l3), L4 = c("age", "male", l2, l3, l4))
  res <- lapply(covs, function(cc) fit_sz(E, cv, cc))
  # surrogate variables
  sv_fit <- function(cc) {
    cc <- cc[vapply(cc, function(v) length(unique(cv[[v]])) > 1, logical(1))]
    mod <- model.matrix(as.formula(paste("~", paste(c("sz", cc), collapse = " + "))), data = cv)
    mod0 <- mod[, colnames(mod) != "sz", drop = FALSE]
    Y <- E[, cv$gsm]
    nsv <- tryCatch(num.sv(Y, mod, method = "be", seed = 20260925), error = function(e) NA_integer_)
    if (is.na(nsv) || nsv == 0) return(list(fit = fit_sz(E, cv, cc), nsv = ifelse(is.na(nsv), NA_integer_, 0L)))
    s <- tryCatch({ capture.output(s <- sva(Y, mod, mod0, n.sv = nsv)); s }, error = function(e) NULL)
    if (is.null(s)) return(list(fit = NULL, nsv = nsv))
    S <- as.matrix(s$sv); colnames(S) <- paste0("SV", seq_len(ncol(S)))
    list(fit = fit_sz(E, cv, cc, extra_design = S), nsv = nsv)
  }
  a <- sv_fit(covs$L1); b <- sv_fit(covs$L4)
  res$L5a <- a$fit; res$L5b <- b$fit
  out <- data.table(gene = rownames(E))
  for (lv in LEVELS) if (!is.null(res[[lv]])) {
    out[, (paste0("logfc_", lv)) := signif(res[[lv]]$logfc, 5)]
    out[, (paste0("t_", lv)) := signif(res[[lv]]$t, 5)]
    out[, (paste0("p_", lv)) := signif(res[[lv]]$p, 5)]
  }
  fwrite(out, sprintf("results/b_de/%s.tsv.gz", x$name), sep = "\t")
  lv <- data.table(dataset = x$name, compartment = x$compartment, n_sz = sum(cv$sz), n_ctl = sum(1 - cv$sz), n_genes = nrow(E),
                   L2_adds = paste(l2, collapse = "+"), L3_adds = paste(l3, collapse = "+"), L4_adds = paste(l4, collapse = "+"),
                   nsv_L5a = a$nsv, nsv_L5b = b$nsv)
  for (l in LEVELS) lv[, (paste0("nq05_", l)) := if (is.null(res[[l]])) NA_integer_ else sum(p.adjust(res[[l]]$p, "BH") < 0.05)]
  list(levels = lv, res = lapply(res, function(r) if (is.null(r)) NULL else r[c("gene", "t", "logfc", "p")]))
}
S <- mclapply(L, function(x) tryCatch(series(x), error = function(e) { message(x$name, ": ", conditionMessage(e)); NULL }),
              mc.cores = NCORES, mc.preschedule = FALSE)
names(S) <- names(L)
S <- S[!vapply(S, is.null, logical(1))]
LV <- rbindlist(lapply(S, `[[`, "levels"))
fwrite(LV, "results/b_levels.tsv", sep = "\t")
print(LV)

# cross-dataset agreement per level
P <- read_pairs()
P <- P[dataset_a %in% names(S) & dataset_b %in% names(S)]
rows <- mclapply(seq_len(nrow(P)), function(i) {
  a <- S[[P$dataset_a[i]]]$res; b <- S[[P$dataset_b[i]]]$res
  rbindlist(lapply(LEVELS, function(lv) if (!is.null(a[[lv]]) && !is.null(b[[lv]]))
    cbind(P[i, .(dataset_a, dataset_b, relation, both_bulk)], data.table(level = lv), suppressWarnings(agree(a[[lv]], b[[lv]])))))
}, mc.cores = NCORES)
BP <- rbindlist(rows)
fwrite(BP, "results/b_pairs.tsv", sep = "\t")

# independent set: bulk datasets pairwise donor-separate, greedy by effective n = 4 / (1/n_sz + 1/n_ctl)
LV[, eff_n := 4 / (1 / n_sz + 1 / n_ctl)]
cand <- LV[compartment == "bulk"][order(-eff_n)]
sep <- function(a, b) { r <- P[(dataset_a == a & dataset_b == b) | (dataset_a == b & dataset_b == a), relation]; length(r) == 1 && r == "subject-separate" }
ind <- character(0)
for (d in cand$dataset) if (all(vapply(ind, function(o) sep(d, o), logical(1)))) ind <- c(ind, d)
fwrite(cand[, .(dataset, n_sz, n_ctl, eff_n = round(eff_n, 1), in_independent_set = dataset %in% ind)], "results/b_independent_set.tsv", sep = "\t")
cat("independent set:", ind, "\n")

zof <- function(r) { z <- qnorm(pmax(r$p, 1e-300) / 2, lower.tail = FALSE) * sign(r$t); setNames(z, r$gene) }
rs <- fread("../confounds/recurrent_sets.tsv", skip = "set\t")[set == "study_brain"]
rs[, sz_sign := ifelse(sz_direction == "up", 1, ifelse(sz_direction == "down", -1, 0))]
RG <- list(); GC <- list()
for (lv in LEVELS) {
  Z <- sapply(ind, function(d) { r <- S[[d]]$res[[lv]]; if (is.null(r)) return(rep(NA_real_, 0)); zof(r) }, simplify = FALSE)
  genes <- Reduce(union, lapply(Z, names))
  M <- sapply(Z, function(z) z[genes]); rownames(M) <- genes
  k <- rowSums(!is.na(M))
  st <- rowSums(M, na.rm = TRUE) / sqrt(k)
  GC[[lv]] <- data.table(level = lv, gene = genes, n_datasets = k, stouffer_z = round(st, 4), n_pos = rowSums(M > 0, na.rm = TRUE), n_neg = rowSums(M < 0, na.rm = TRUE))
  mm <- M[match(rs$gene, genes), , drop = FALSE]
  sg <- rs$sz_sign
  RG[[lv]] <- data.table(level = lv, gene = rs$gene, sz_direction = rs$sz_direction, n_datasets = rowSums(!is.na(mm)),
                         signed_stouffer_z = round(sg * rowSums(mm, na.rm = TRUE) / sqrt(rowSums(!is.na(mm))), 4),
                         n_agree_sign = rowSums(sign(mm) == sg, na.rm = TRUE),
                         n_agree_p05 = rowSums(sign(mm) == sg & abs(mm) > qnorm(0.975), na.rm = TRUE))
}
fwrite(rbindlist(RG), "results/b_recurrent_genes.tsv", sep = "\t")
fwrite(rbindlist(GC)[n_datasets >= 3], "results/b_gene_consistency.tsv.gz", sep = "\t")
