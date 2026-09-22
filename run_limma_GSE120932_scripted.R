suppressPackageStartupMessages(library(limma))

# Reproduce the GSE120932 GEO2R drug-withdrawn K562-IR vs parental contrast
# from the public, preprocessed GEO series matrix (GPL10558; 12 arrays).
# Usage: Rscript run_limma_GSE120932_scripted.R \
#     GSE120932_series_matrix.txt.gz results_gse120932 GSE120932.top.table.tsv
# The GEO2R table is optional; when supplied it is used only to annotate and
# cross-check the newly fitted limma output, not as a model input.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2, 3))) {
  stop("Usage: Rscript run_limma_GSE120932_scripted.R series_matrix.txt.gz output_dir [GSE120932.top.table.tsv]")
}
source_file <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lines <- readLines(if (grepl("\\.gz$", source_file)) gzfile(source_file) else source_file,
                   warn = FALSE)
begin <- grep("^!series_matrix_table_begin$", lines)
end <- grep("^!series_matrix_table_end$", lines)
stopifnot(length(begin) == 1L, length(end) == 1L, end > begin + 1L)

sample_field <- function(prefix) {
  hits <- lines[startsWith(lines, paste0(prefix, "\t"))]
  stopifnot(length(hits) == 1L)
  fields <- strsplit(hits, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", fields)
}
all_gsm <- sample_field("!Sample_geo_accession")
all_titles <- sample_field("!Sample_title")
stopifnot(length(all_gsm) == 12L, length(all_titles) == 12L)
metadata <- data.frame(sample = all_gsm, title = all_titles, stringsAsFactors = FALSE)
metadata$group <- ifelse(grepl("^K562, replicate ", all_titles), "Parental",
                         ifelse(grepl("^K562-IR w/o imatinib, replicate ", all_titles),
                                "Resistant_without_imatinib", "Excluded_other_group"))
write.csv(metadata, file.path(outdir, "GSE120932_all_array_metadata.csv"), row.names = FALSE)

parental <- c("GSM3421746", "GSM3421747", "GSM3421748")
resistant <- c("GSM3421752", "GSM3421753", "GSM3421754")
selected <- c(parental, resistant)
stopifnot(identical(metadata$sample[metadata$group == "Parental"], parental))
stopifnot(identical(metadata$sample[metadata$group == "Resistant_without_imatinib"], resistant))

con <- if (grepl("\\.gz$", source_file)) gzfile(source_file) else file(source_file)
table <- read.delim(con, skip = begin, nrows = end - begin - 2L,
                    check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
stopifnot(identical(names(table)[-1], all_gsm), !anyDuplicated(table[[1]]))
expression <- as.matrix(table[, selected])
rownames(expression) <- table[[1]]
storage.mode(expression) <- "double"
stopifnot(nrow(expression) == 47223L, all(is.finite(expression)))

# The submitted series matrix is already log2, variance stabilized and
# quantile normalized. No second log transform or normalization is applied.
stopifnot(all(expression > 0), max(expression) < 30)
group <- factor(c(rep("Parental", 3), rep("Resistant", 3)),
                levels = c("Parental", "Resistant"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(Resistant_vs_Parental = Resistant - Parental, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expression, design), contrast))
result <- topTable(fit, coef = "Resistant_vs_Parental", number = Inf,
                   adjust.method = "BH", sort.by = "none")
result$ID <- rownames(result)
result <- result[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(result, file.path(outdir, "GSE120932_limma_scripted_full.csv"), row.names = FALSE)

# Array-level QC for the six included, already normalized arrays.
qc <- data.frame(sample = selected, group = as.character(group),
                 median_log2 = apply(expression, 2, median),
                 IQR_log2 = apply(expression, 2, IQR),
                 stringsAsFactors = FALSE)
sample_cor <- cor(expression, method = "pearson")
qc$mean_correlation_to_other_five <- vapply(seq_along(selected), function(i) {
  mean(sample_cor[i, -i])
}, numeric(1))
qc$mean_correlation_within_group <- vapply(seq_along(selected), function(i) {
  others <- which(group == group[i] & seq_along(selected) != i)
  mean(sample_cor[i, others])
}, numeric(1))
write.csv(qc, file.path(outdir, "GSE120932_array_QC_summary.csv"), row.names = FALSE)
write.csv(sample_cor, file.path(outdir, "GSE120932_array_correlations.csv"))
pca <- prcomp(t(expression), center = TRUE, scale. = FALSE)
coordinates <- data.frame(sample = selected, group = as.character(group),
                          PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                          PC1_percent = 100 * summary(pca)$importance[2, 1],
                          PC2_percent = 100 * summary(pca)$importance[2, 2],
                          row.names = NULL)
write.csv(coordinates, file.path(outdir, "GSE120932_array_PCA_coordinates.csv"), row.names = FALSE)
pdf(file.path(outdir, "GSE120932_array_QC_plots.pdf"), width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
boxplot(as.data.frame(expression), las = 2, cex.axis = 0.7,
        col = rep(c("#4477AA", "#EE8855"), each = 3),
        main = "GSE120932 six arrays", ylab = "Source log2 signal")
plot(coordinates$PC1, coordinates$PC2, pch = 19, cex = 1.4,
     col = rep(c("#4477AA", "#EE8855"), each = 3),
     xlim = range(coordinates$PC1) + c(-16, 16),
     xlab = sprintf("PC1 (%.1f%%)", coordinates$PC1_percent[1]),
     ylab = sprintf("PC2 (%.1f%%)", coordinates$PC2_percent[1]),
     main = "Unsupervised PCA")
text(coordinates$PC1, coordinates$PC2,
     labels = c("WT1", "WT2", "WT3", "IR1", "IR2", "IR3"),
     pos = c(1, 3, 4, 2, 2, 2), cex = 0.75)
dev.off()

if (length(args) == 3) {
  geo2r <- read.delim(args[3], check.names = FALSE, stringsAsFactors = FALSE)
  matched <- match(result$ID, geo2r$ID)
  stopifnot(!anyNA(matched))
  result$Gene.symbol <- geo2r$Gene.symbol[matched]
  b56 <- c("PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E")
  write.csv(result[result$Gene.symbol %in% b56, ],
            file.path(outdir, "GSE120932_B56_scripted_limma_results.csv"), row.names = FALSE)
  comparison <- data.frame(ID = result$ID,
                           GEO2R_logFC = geo2r$logFC[matched],
                           scripted_logFC = result$logFC,
                           GEO2R_adjP = geo2r$adj.P.Val[matched],
                           scripted_adjP = result$adj.P.Val)
  write.csv(comparison, file.path(outdir, "GSE120932_GEO2R_vs_scripted_comparison.csv"), row.names = FALSE)
  cat("Maximum absolute GEO2R vs refitted logFC difference:",
      max(abs(comparison$GEO2R_logFC - comparison$scripted_logFC)), "\n")
}
cat("GSE120932: six selected arrays;", nrow(expression), "probes.\n")
cat("B56 and full limma results, sample QC, correlations and PCA written to", outdir, "\n")
