#!/usr/bin/env Rscript
# Donor-level cell-type composition, SCZ vs control, PsychAD HBCC and Mount Sinai (MSSM), across donor-age rules.
# Run: Rscript scripts/composition_age_cut.R <rerun_full_root>   (default rerun-full)
# Reads (read-only): <root>/analyses/psychad_pseudobulk/pseudobulk/{HBCC_Cohort,MSSM_Cohort}/cols_subclass.tsv  (donor x subclass nuclei counts, 27 PsychAD subclasses)
#                    <root>/data/raw/42778699/publisher/moesm3/2024-07-15352A-s3/PsychAD_SupplementaryTable1.csv (Lee donor table: Cohort, Age, Sex, PMI, Tier1_crossDis, Tier1_crossDis_dx)
#                    <root>/analyses/psychad_pseudobulk/results/donors_used.tsv (this project's HBCC_main / MSSM_main donor sets)
# Models (limma on CLR-transformed counts with crumblr precision weights; empirical-Bayes moderated t):
#   lee      all 27 subclasses as one composition; Composition ~ Dx + scale(Age) + Sex [+ Source when pooled]   (Lee 2026 Crumblr methods)
#   endresz  neuronal (EN_*, IN_*) and non-neuronal subclasses as two separate compositions; Composition ~ Dx + scale(Age) + Sex + scale(PMI) [+ Source]
#            (Endresz 2026 Methods, eq. 2; supertype resolution there, subclass here)
# Output: results/composition_age_cut.tsv
suppressMessages({library(data.table); library(limma)})
# crumblr (DiseaseNeurogenomics/crumblr 1.4.5, GitHub main tarball sha256 4b8d9d9b4a3b54cab944a4207ead9e7d0f6a03e4ff2540a49a3e179fba1d2d37) is not installable on remote
# (its Bioconductor dependency chain fails to compile); its three small functions are copied verbatim from R/crumblr.R (scripts/crumblr_upstream_crumblr.R).
clr <- function(counts, pseudocount = 0.5) { if (!is.data.frame(counts) & !is.matrix(counts)) counts <- matrix(counts, nrow = 1); log(counts + pseudocount) - rowMeans(log(counts + pseudocount)) }
.cap_ratio <- function(W, max.ratio = 5, quant = 0.05) { t(apply(W, 1, function(x) { x <- x / quantile(x, quant); pmin(x, max.ratio) })) }
crumblr <- function(counts, pseudocount = 0.5, tau = 1, max.ratio = 5, quant = 0.05) {
  D <- ncol(counts)
  var_asymp <- apply(counts + pseudocount, 1, function(x) { C <- sum(x); p <- x / C; tau * (1 / p - 2 / (p * D) + sum(1 / p) / D^2) / C })
  weights <- .cap_ratio(1 / var_asymp, max.ratio, quant)
  new("EList", list(E = t(clr(counts, pseudocount)), weights = weights))
}
args <- commandArgs(TRUE); ROOT <- if (length(args)) args[1] else "rerun-full"
PB <- file.path(ROOT, "analyses/psychad_pseudobulk/pseudobulk")
T1 <- file.path(ROOT, "data/raw/42778699/publisher/moesm3/2024-07-15352A-s3/PsychAD_SupplementaryTable1.csv")
DU <- file.path(ROOT, "analyses/psychad_pseudobulk/results/donors_used.tsv")
cnt <- rbindlist(lapply(c(HBCC = "HBCC_Cohort", MSSM = "MSSM_Cohort"), function(b) fread(file.path(PB, b, "cols_subclass.tsv"))[, .(donor_id, cell_type, n_nuclei)]), idcol = "Source")
W <- dcast(cnt, Source + donor_id ~ cell_type, value.var = "n_nuclei", fill = 0)
types <- setdiff(names(W), c("Source", "donor_id"))
neuronal <- grepl("^(EN|IN)_", types)
lee <- fread(T1); lee[, Age := as.numeric(sub("\\+", "", Age))]
lee[, Source := fifelse(Cohort == "HBCC", "HBCC", fifelse(Cohort == "MSSM", "MSSM", "RADC"))]
duse <- fread(DU, skip = "set\t")[set %in% c("HBCC_main", "MSSM_main")]
duse[, Source := sub("_main", "", set)]
mk <- function(base) {  # base: data.table donor_id, Source, sz, age, male, pmi
  merge(base, W, by = c("Source", "donor_id"))
}
setA <- mk(lee[Tier1_crossDis == "Y" & Source %in% c("HBCC", "MSSM") & Tier1_crossDis_dx %in% c("CTRL", "SCZ"),
               .(donor_id = DonorID, Source, sz = as.integer(Tier1_crossDis_dx == "SCZ"), dxf = Tier1_crossDis_dx, age = Age, male = as.integer(Sex == "M"), pmi = PMI)])
setB <- mk(duse[, .(donor_id, Source, sz, dxf = fifelse(sz == 1, "SCZ", "CTRL"), age, male, pmi)])
# all Tier1 cross-disorder donors of HBCC and MSSM with every diagnosis in the model (CTRL is the reference), as Lee's shared-control contrasts do
setC <- mk(lee[Tier1_crossDis == "Y" & Source %in% c("HBCC", "MSSM"),
               .(donor_id = DonorID, Source, sz = as.integer(Tier1_crossDis_dx == "SCZ"), dxf = Tier1_crossDis_dx, age = Age, male = as.integer(Sex == "M"), pmi = PMI)])
SETS <- list(lee_crossdisorder = setA, lee_crossdisorder_alldx = setC, project_main = setB)
fit_one <- function(d, model) {
  if (sum(d$sz == 1) < 3 || sum(d$dxf == "CTRL") < 3) return(NULL)
  idx <- if (model == "lee") list(all = seq_along(types)) else list(neuronal = which(neuronal), nonneuronal = which(!neuronal))
  rbindlist(lapply(names(idx), function(g) {
    cm <- as.matrix(d[, ..types])[, idx[[g]], drop = FALSE]; rownames(cm) <- d$donor_id
    ok <- rowSums(cm) > 0; d2 <- d[ok]; cm <- cm[ok, , drop = FALSE]
    keep <- colSums(cm) > 0; cm <- cm[, keep, drop = FALSE]
    ec <- suppressWarnings(crumblr(cm))
    df <- data.frame(dx = relevel(factor(d2$dxf), "CTRL"), age = as.numeric(scale(d2$age)), sex = d2$male)
    if (model == "endresz") df$pmi <- as.numeric(scale(d2$pmi))
    if (length(unique(d2$Source)) > 1) df$src <- factor(d2$Source)
    if (model == "endresz") { df <- df[, !sapply(df, function(x) all(is.na(x))), drop = FALSE]; keepd <- complete.cases(df); df <- df[keepd, , drop = FALSE]; ec <- ec[, keepd] }
    des <- model.matrix(~ ., df)
    f <- eBayes(lmFit(ec, des))
    tt <- topTable(f, coef = "dxSCZ", number = Inf, sort.by = "none")
    se <- sqrt(f$s2.post) * f$stdev.unscaled[, "dxSCZ"]
    data.table(cell_type = rownames(tt), estimate = tt$logFC, SE = se, p = tt$P.Value, FDR = p.adjust(tt$P.Value, "BH"),
               n_sz = sum(df$dx == "SCZ"), n_ctl = sum(df$dx == "CTRL"), composition = g)
  }))
}
cut_sets <- function(d) {
  med <- d[, .(m = median(age)), by = Source]; d <- merge(d, med, by = "Source")
  list(none = d,
       age_lt70 = d[age < 70],
       endresz_rule = d[age < 70 & (Source != "HBCC" | age > 20)],
       age_le_cohort_median = d[age <= m],
       age_gt_cohort_median = d[age > m])
}
res <- list()
for (sn in names(SETS)) for (cn in names(cut_sets(SETS[[sn]]))) for (mod in c("lee", "endresz")) {
  d <- cut_sets(SETS[[sn]])[[cn]]
  for (co in c("HBCC", "MSSM", "pooled")) {
    dd <- if (co == "pooled") d else d[Source == co]
    r <- fit_one(dd, mod); if (is.null(r)) next
    res[[length(res) + 1]] <- cbind(donor_set = sn, cut = cn, model = mod, cohort = co, r)
  }
}
R <- rbindlist(res)
R[, cohort_median_age := NA_real_]
write_hdr <- c(paste0("# built ", format(Sys.Date()), " by psychad_cell_proportions_by_age/scripts/composition_age_cut.R; inputs under ", ROOT),
  paste0("# sha256 ", paste(c(basename(T1), basename(DU)), system2("sha256sum", c(T1, DU), stdout = TRUE) |> sub(" .*", "", x = _), collapse = "; ")),
  "# estimate = SCZ coefficient on the CLR composition scale (crumblr); SE from limma eBayes; FDR = BH across subclasses within one (donor_set, cut, model, cohort, composition) fit")
out <- "results/composition_age_cut.tsv"; writeLines(write_hdr, out)
fwrite(R[, !"cohort_median_age"], out, sep = "\t", append = TRUE, col.names = TRUE)
cat("rows", nrow(R), "\n")
