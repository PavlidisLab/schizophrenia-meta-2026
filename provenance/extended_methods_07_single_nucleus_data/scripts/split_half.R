# (b) Split halves within HBCC and within MSSM at matched n, and HBCC vs MSSM at the same n, per level x cell type.
#
# Run from analyses/psychad_pseudobulk/:  Rscript scripts/split_half.R      (after scripts/de.R)
# For each cell type analysed in both HBCC_main and MSSM_main: m_sz = floor(min(SZ donors HBCC, MSSM) / 2), m_ctl likewise
# for controls. Each replicate draws, in each bank, 2 m_sz SZ and 2 m_ctl control donors at random and splits them into
# halves H1, H2 (HBCC) and M1, M2 (MSSM), each m_sz vs m_ctl. Four fits (pipeline of scripts/pb_common.R) and four
# agreements: H1 vs H2 (split halves, HBCC), M1 vs M2 (split halves, MSSM), H1 vs M1 and H2 vs M2 (HBCC vs MSSM at the
# same n). Also MSSM halves at MSSM's own n (floor(n/2) per group): split_own_MSSM. NREP replicates.
# Outputs: results/b_reps.tsv.gz (one row per replicate x comparison), results/b_summary.tsv (medians per cell type x comparison)
SCRIPT <- "scripts/split_half.R"
source("scripts/pb_common.R")
NREP <- as.integer(Sys.getenv("NREP", "20"))
D <- load_donors()
S <- fread("results/de_summary.tsv", skip = "set\t")
LEVELS <- c("total", "class", "subclass")
out <- list()
for (lv in LEVELS) {
  LH <- load_level("HBCC", lv); LM <- load_level("MSSM", lv)
  cts <- intersect(S[set == "HBCC_main" & level == lv, cell_type], S[set == "MSSM_main" & level == lv, cell_type])
  for (ct in cts) {
    dH <- cell_samples(LH, D, "HBCC_main", ct); dM <- cell_samples(LM, D, "MSSM_main", ct)
    m_sz <- floor(min(sum(dH$sz), sum(dM$sz)) / 2); m_ctl <- floor(min(sum(1 - dH$sz), sum(1 - dM$sz)) / 2)
    o_sz <- floor(sum(dM$sz) / 2); o_ctl <- floor(sum(1 - dM$sz) / 2)
    cat(lv, ct, "matched n", m_sz, "vs", m_ctl, "; MSSM own", o_sz, "vs", o_ctl, "\n")
    rows <- mclapply(seq_len(NREP), function(r) {
      set.seed(SEED + r * 1000 + match(ct, cts) + 100000 * match(lv, LEVELS))
      halves <- function(d, ns, nc) {
        s <- sample(d[sz == 1, pb_id], 2 * ns); c <- sample(d[sz == 0, pb_id], 2 * nc)
        list(d[pb_id %in% c(s[1:ns], c[1:nc])], d[pb_id %in% c(s[(ns + 1):(2 * ns)], c[(nc + 1):(2 * nc)])])
      }
      h <- halves(dH, m_sz, m_ctl); m <- halves(dM, m_sz, m_ctl); mo <- halves(dM, o_sz, o_ctl)
      fH <- lapply(h, function(d) fit_pb(LH$M, LH$genes, d)); fM <- lapply(m, function(d) fit_pb(LM$M, LM$genes, d))
      fO <- lapply(mo, function(d) fit_pb(LM$M, LM$genes, d))
      mk <- function(cmp, a, b, ns, nc) if (is.null(a) || is.null(b)) NULL else cbind(data.table(level = lv, cell_type = ct, rep = r, comparison = cmp, n_sz = ns, n_ctl = nc), agree(a, b))
      rbind(mk("split_HBCC", fH[[1]], fH[[2]], m_sz, m_ctl), mk("split_MSSM", fM[[1]], fM[[2]], m_sz, m_ctl),
            mk("HBCC_vs_MSSM", fH[[1]], fM[[1]], m_sz, m_ctl), mk("HBCC_vs_MSSM", fH[[2]], fM[[2]], m_sz, m_ctl),
            mk("split_own_MSSM", fO[[1]], fO[[2]], o_sz, o_ctl))
    }, mc.cores = NCORES)
    out[[paste(lv, ct)]] <- rbindlist(rows)
  }
  rm(LH, LM); gc()
}
R <- rbindlist(out)
write_tsv(R, "results/b_reps.tsv.gz", c(sprintf("%d replicates per cell type; design in the header of scripts/split_half.R", NREP),
                                        "agree() from ../fairness_sample_level/scripts/common.R"))
num <- c("rho_t", "rho_logfc", "top20_fold", "top100_fold", "top200_fold", "top100_sign_in_other", "fdr_fold", "pi1", "pi1_base")
B <- R[, c(list(n_sz = n_sz[1], n_ctl = n_ctl[1], n_reps = .N), lapply(.SD, function(v) round(median(v, na.rm = TRUE), 4))),
       by = .(level, cell_type, comparison), .SDcols = num]
write_tsv(B, "results/b_summary.tsv", c("medians over replicates of results/b_reps.tsv.gz (HBCC_vs_MSSM: 2 pairs per replicate)"))
print(dcast(B, level + cell_type ~ comparison, value.var = "rho_t"))
cat("done\n")
