suppressPackageStartupMessages(library(limma))

# GSE130404 covariate adjusted sensitivity analysis.
# Reruns the EMR failure versus EMR achieved comparison with age and gender
# added to the design matrix, alongside the original unadjusted contrast, so
# the two can be compared directly. This addresses the referee's request for
# an adjusted sensitivity model given the 13 versus 83 group imbalance.
#
# Usage:
#   Rscript run_limma_GSE130404_covariate_adjusted.R \
#       GSE130404_series_matrix.txt \
#       GSE130404_sample_metadata_derived.csv \
#       output_dir
#
# GSE130404_series_matrix.txt is the plain text series matrix downloaded from
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130404/matrix/
# GSE130404_sample_metadata_derived.csv is produced alongside this script and
# supplies group, age, and gender per sample in the same row order as the
# expression matrix columns.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_limma_GSE130404_covariate_adjusted.R series_matrix.txt sample_metadata.csv output_dir")
}
matrix_file <- args[1]
meta_file <- args[2]
outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Read the series matrix expression table -------------------------------
# The series matrix file is metadata lines prefixed with "!", followed by the
# expression table between !series_matrix_table_begin and _end. Probes are
# rows (ILMN_ identifiers), samples are columns (GSM accessions), values are
# neqc-normalized log2 intensities as deposited by the submitters.
lines <- readLines(matrix_file)
start <- grep("^!series_matrix_table_begin", lines)
end <- grep("^!series_matrix_table_end", lines)
tbl <- read.delim(matrix_file, skip = start, nrows = end - start - 1,
                   header = TRUE, check.names = FALSE, row.names = 1,
                   quote = "\"")
expr <- as.matrix(tbl)
storage.mode(expr) <- "double"

# --- Read sample metadata and align sample order ----------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
rownames(meta) <- meta$geo_accession
meta <- meta[colnames(expr), ]
stopifnot(all(rownames(meta) == colnames(expr)))

group <- factor(meta$group, levels = c("EMR_achieved", "EMR_failure"))
age_c <- meta$age - mean(meta$age)
gender <- factor(meta$gender)

# --- Unadjusted model, group only -------------------------------------------
design_unadj <- model.matrix(~group)
fit_unadj <- lmFit(expr, design_unadj)
fit_unadj <- eBayes(fit_unadj)
res_unadj <- topTable(fit_unadj, coef = "groupEMR_failure", number = Inf, sort.by = "none")
res_unadj$ID <- rownames(res_unadj)

# --- Adjusted model, group plus age plus gender -----------------------------
design_adj <- model.matrix(~group + age_c + gender)
fit_adj <- lmFit(expr, design_adj)
fit_adj <- eBayes(fit_adj)
res_adj <- topTable(fit_adj, coef = "groupEMR_failure", number = Inf, sort.by = "none")
res_adj$ID <- rownames(res_adj)

write.csv(res_unadj, file.path(outdir, "GSE130404_limma_unadjusted_full.csv"), row.names = FALSE)
write.csv(res_adj, file.path(outdir, "GSE130404_limma_age_gender_adjusted_full.csv"), row.names = FALSE)

# --- B56 subset and side by side comparison ---------------------------------
b56_probes <- c(
  ILMN_1738784 = "PPP2R5A",
  ILMN_2124082 = "PPP2R5B",
  ILMN_3279757 = "PPP2R5C", ILMN_1795846 = "PPP2R5C", ILMN_1780913 = "PPP2R5C",
  ILMN_1780940 = "PPP2R5D", ILMN_1699384 = "PPP2R5D", ILMN_2359887 = "PPP2R5D",
  ILMN_1666761 = "PPP2R5E"
)

merge_b56 <- function(res, suffix) {
  x <- res[res$ID %in% names(b56_probes), c("ID", "logFC", "P.Value", "adj.P.Val")]
  x$gene <- b56_probes[x$ID]
  colnames(x)[2:4] <- paste0(colnames(x)[2:4], "_", suffix)
  x
}

b56_unadj <- merge_b56(res_unadj, "unadj")
b56_adj <- merge_b56(res_adj, "adj")
b56_compare <- merge(b56_unadj, b56_adj, by = c("ID", "gene"))
b56_compare <- b56_compare[order(b56_compare$gene, b56_compare$ID), ]
write.csv(b56_compare, file.path(outdir, "GSE130404_B56_unadjusted_vs_adjusted_comparison.csv"),
          row.names = FALSE)

cat("Done. Wrote full unadjusted and age/gender adjusted limma results,\n",
    "plus a B56 side by side comparison, to:", outdir, "\n")
