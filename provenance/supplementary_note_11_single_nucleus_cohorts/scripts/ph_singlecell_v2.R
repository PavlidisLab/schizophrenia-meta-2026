# Does the bulk pH-association signature also track the SZ effect in single-nucleus data, the way it tracks the
# bulk SZ effect (r = -0.68, n = 1,111; ../confounds_report.md, Per-gene relationship; report/report.md #confounds)?
# Extension of the pH-confound analysis to single-cell data (Paul, 2026-10-02): "the pH study needs to be extended
# to that [the single-cell section]."
#
# Z_pH_controls: ../results/recurrent_genes_pH.tsv (built by ../scripts/combine_datasets.R from 5 bulk dataset-
# regions with recorded pH: GSE21138, GSE35977, GSE53987_BA46, GSE17612, GSE92538_GPL96); read only, not rerun here.
#
# Single-nucleus SZ effect: Endresz 2026's uniform voom-limma reprocessing of 7 donor cohorts (same source as
# ../../single_cell/single_cell_report.md and ../../effect_sizes/harmonize.py ENDRESZ_N), per raw Endresz subclass
# (not coarsened to Gemma's 14 classes: that mapping, ../../single_cell/celltype_map.tsv, is built for cross-study
# label reconciliation and buys nothing here, where every subclass already comes from one pipeline). Read from
# ../../effect_sizes/data/endresz_DE_subclass_all_cohorts.tsv.gz, the plain-TSV conversion of the Zenodo parquet
# that ../../single_cell/sc_cross_prep.py and ../../effect_sizes/harmonize.py already use (no R/Python parquet
# reader needed here or on remote: this repo has no R `arrow` package installed on remote, and the project's own
# convention for this file is to go through pandas once, not read the parquet directly in every consumer).
# zt() and stouffer() are copied verbatim from ../scripts/combine_datasets.R (same formulas; that script also does
# GEO downloads this one does not need).
#
# Patched copy of confounds/ph_singlecell_2026-10-02/scripts/ph_singlecell.R: absolute input and output paths only (inputs from confounds64_2026-10-13/out/results_H; outputs in followups_confounds64_2026-10-13/fig09/)
suppressMessages(library(data.table))

RG_PATH <- "fdv2_2026-10-15/work/axis_v2_tables/recurrent_genes_pH.tsv"
ENDRESZ_PATH <- "sztr-rebuild-infra/repo/analyses/effect_sizes/data/endresz_DE_subclass_all_cohorts.tsv.gz"

zt <- function(t, p) sign(t) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
stouffer <- function(z, w) sum(w * z) / sqrt(sum(w^2))

# Endresz 2026 Table 1 (donors < 70), SCZ + CON; effect_sizes/harmonize.py ENDRESZ_N
ENDRESZ_N <- c(Batiuk = 15, McLean = 32, MSSM1 = 35, Multiome = 11, Frohlich = 61, HBCC = 130, MSSM2 = 185)
ENDRESZ_OF <- c(Batiuk = "Batiuk 2022", McLean = "Ruzicka 2024 (McLean)", MSSM1 = "Ruzicka 2024 (MSSM)",
                Multiome = "Emani 2024 (MultiomeBrain)", Frohlich = "Frohlich 2024", HBCC = "PsychAD (HBCC)",
                MSSM2 = "PsychAD (MSSM)")

sha256 <- function(f) sub(" .*", "", system(paste("shasum -a 256", shQuote(f)), intern = TRUE))

rg <- fread(RG_PATH)
rg <- rg[!is.na(Z_pH_controls), .(gene, Z_pH_controls)]
cat("recurrent genes with Z_pH_controls:", nrow(rg), "\n")

d_all <- fread(ENDRESZ_PATH)
d_all <- d_all[cohort %in% names(ENDRESZ_N)]
d_all[, z := zt(logFC / SE, PValue)]
d_all[, w := sqrt(ENDRESZ_N[cohort])]
cat("rows (cohort x subclass x gene), 7 real cohorts, before any gene restriction:", nrow(d_all), "\n")

d <- d_all[gene %in% rg$gene]
cat("rows after restricting to the pH gene set:", nrow(d), "\n")

# per gene x subclass: Stouffer-combined z across the Endresz cohorts holding that subclass (weight sqrt(cohort N),
# same weighting rule as the bulk Z_sz_base / Z_pH_controls combination)
per_class <- d[, .(n_cohorts = .N, Z_sc = stouffer(z, w)), by = .(cell_type, gene)]
per_class <- merge(per_class, rg, by = "gene")

class_summary <- per_class[, .(n_genes = .N, median_n_cohorts = as.double(median(n_cohorts)),
                                r = round(cor(Z_pH_controls, Z_sc), 3)), by = cell_type]
setorder(class_summary, cell_type)
cat("\nPer-subclass correlation (Z_pH_controls vs. Endresz Stouffer-combined single-nucleus SZ Z):\n")
print(class_summary)
cat("\nMedian r over", nrow(class_summary), "subclasses:", median(class_summary$r), "\n")

# sign check: direction of the published bulk result (genes that respond to pH in controls tend to move the
# opposite way in SZ) -- how many subclasses show a negative r, and how extreme is each
cat("Subclasses with r < 0:", sum(class_summary$r < 0), "of", nrow(class_summary), "\n")

per_class_gene <- per_class[, .(cell_type, gene, n_cohorts, Z_sc, Z_pH_controls)]
setorder(per_class_gene, cell_type, gene)

# ---- robustness: genome-wide (not restricted to the recurrent gene list), same single-nucleus side ------------
GW_PATH <- "fdv2_2026-10-15/work/axis_v2_tables/combined_controls_pH.tsv"
gw <- fread(GW_PATH)[factor == "pH", .(gene, Z_pH_controls_genomewide = Z)]
cat("\ngenome-wide pH-controls genes:", nrow(gw), "\n")
per_class_gw <- d_all[, .(Z_sc = stouffer(z, w)), by = .(cell_type, gene)]
per_class_gw <- merge(per_class_gw, gw, by = "gene")
class_summary_gw <- per_class_gw[, .(n_genes = .N, r = round(cor(Z_pH_controls_genomewide, Z_sc), 3)), by = cell_type]
setorder(class_summary_gw, cell_type)
cat("\nGenome-wide per-subclass correlation (not restricted to recurrent genes):\n")
print(class_summary_gw)
cat("\nMedian r genome-wide over", nrow(class_summary_gw), "subclasses:", median(class_summary_gw$r), "\n")
cat("Subclasses with r < 0 (genome-wide):", sum(class_summary_gw$r < 0), "of", nrow(class_summary_gw), "\n")

# ---- added 2026-10-15: study_brain_ge10 population (genes in 10 or more lists with a combined control-pH z) ----
ge10_genes <- fread(RG_PATH)[grepl("study_brain_ge10", sets) & !is.na(Z_pH_controls), gene]
cat("\nge10 genes with Z_pH_controls:", length(ge10_genes), "\n")
per_class_ge10 <- per_class[gene %in% ge10_genes]
class_summary_ge10 <- per_class_ge10[, .(n_genes = .N, median_n_cohorts = as.double(median(n_cohorts)),
                                r = round(cor(Z_pH_controls, Z_sc), 3)), by = cell_type]
setorder(class_summary_ge10, cell_type)
print(class_summary_ge10)
cat("\nMedian r over", nrow(class_summary_ge10), "subclasses (ge10):", median(class_summary_ge10$r, na.rm = TRUE), "\n")
cat("Subclasses with r < 0 (ge10):", sum(class_summary_ge10$r < 0, na.rm = TRUE), "of", nrow(class_summary_ge10), "\n")

prov <- c(
  sprintf("# built %s by confounds/ph_singlecell_2026-10-02/scripts/ph_singlecell.R", Sys.Date()),
  sprintf("# input: confounds/results/recurrent_genes_pH.tsv sha256 %s", sha256(RG_PATH)),
  sprintf("# input: confounds/results/combined_controls_pH.tsv sha256 %s", sha256(GW_PATH)),
  sprintf("# input: effect_sizes/data/endresz_DE_subclass_all_cohorts.tsv.gz sha256 %s", sha256(ENDRESZ_PATH)),
  "# zt(t,p) = sign(t)*qnorm(p/2, upper tail); stouffer(z,w) = sum(w*z)/sqrt(sum(w^2)); weight w = sqrt(Endresz cohort N) (SCZ+CON, Table 1, donors < 70)"
)
write_tsv <- function(x, path) {
  writeLines(prov, path)
  fwrite(x, path, sep = "\t", append = TRUE, col.names = TRUE)
}
write_tsv(class_summary, "fdv2_2026-10-15/work/fig09/v2/ph_singlecell_by_class.tsv")
write_tsv(per_class_gene, "fdv2_2026-10-15/work/fig09/v2/ph_singlecell_per_gene.tsv")
write_tsv(class_summary_ge10, "fdv2_2026-10-15/work/fig09/v2/ph_singlecell_by_class_ge10.tsv")
write_tsv(class_summary_gw, "fdv2_2026-10-15/work/fig09/v2/ph_singlecell_by_class_genomewide.tsv")

cat("\nWrote ../results/ph_singlecell_by_class.tsv, ph_singlecell_per_gene.tsv, ph_singlecell_by_class_genomewide.tsv\n")
