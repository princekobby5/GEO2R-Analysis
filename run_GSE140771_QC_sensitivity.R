suppressPackageStartupMessages(library(DESeq2))

# GSE140771 quality control sensitivity analysis.
# Addresses two specific referee requests that the primary
# run_featurecounts_deseq2_full.R analysis did not cover:
#   1. A blind=TRUE PCA, as a sensitivity check on the blind=FALSE PCA
#      already reported (Figure 1 / GSE140771_PCA_coordinates.csv).
#   2. A sensitivity analysis rerunning the IR2 versus WT contrast with
#      SRR10507806 (IR2_2), the library with weaker alignment and mapping
#      QC metrics, excluded.
#      This leaves one IR2 library and two WT libraries. The resulting
#      fitted SE and adjusted P value are descriptive model outputs,
#      not evidence of biological replication in IR2 after exclusion.
#
# Usage:
#   Rscript run_GSE140771_QC_sensitivity.R featureCounts_filtered_count_matrix.csv output_dir
#
# featureCounts_filtered_count_matrix.csv is the file already included in
# Supplementary File S1 (6 samples: WT_1, WT_2, IR1_1, IR1_2, IR2_1, IR2_2).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript run_GSE140771_QC_sensitivity.R featureCounts_filtered_count_matrix.csv output_dir")
}
counts_file <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_df <- read.csv(counts_file, check.names = FALSE)
counts <- as.matrix(counts_df[, -1])
rownames(counts) <- counts_df[[1]]
storage.mode(counts) <- "integer"

coldata <- data.frame(
  row.names = colnames(counts),
  group = factor(c("WT", "WT", "IR1", "IR1", "IR2", "IR2"), levels = c("WT", "IR1", "IR2"))
)

# --- 1. blind=TRUE PCA sensitivity check ------------------------------------
dds_full <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~group)
dds_full <- DESeq(dds_full)

vsd_blind <- vst(dds_full, blind = TRUE)
pca_blind <- plotPCA(vsd_blind, intgroup = "group", returnData = TRUE)
percentVar_blind <- round(100 * attr(pca_blind, "percentVar"), 1)
pca_blind_out <- data.frame(sample = pca_blind$name, group = pca_blind$group,
                             PC1 = pca_blind$PC1, PC2 = pca_blind$PC2,
                             PC1_percent = percentVar_blind[1], PC2_percent = percentVar_blind[2])
write.csv(pca_blind_out, file.path(outdir, "GSE140771_PCA_coordinates_blindTRUE.csv"), row.names = FALSE)

cat("blind=TRUE PCA coordinates:\n")
print(pca_blind_out)
cat("\nCompare against the existing blind=FALSE coordinates in\n",
    "GSE140771_PCA_coordinates.csv to see whether WT/IR1/IR2 separation is\n",
    "sensitive to the blind parameter.\n\n")

# --- 2. SRR10507806 (IR2_2) exclusion sensitivity ---------------------------
counts_excl <- counts[, colnames(counts) != "IR2_2"]
stopifnot(identical(colnames(counts_excl), c("WT_1", "WT_2", "IR1_1", "IR1_2", "IR2_1")))
cat("IR2_2 excluded: one IR2 library remains. Interpret the IR2-vs-WT ",
    "fold change as a sensitivity check; the fitted P value and SE ",
    "do not establish biological replication within IR2.\n")
coldata_excl <- data.frame(
  row.names = colnames(counts_excl),
  group = factor(c("WT", "WT", "IR1", "IR1", "IR2"), levels = c("WT", "IR1", "IR2"))
)

dds_excl <- DESeqDataSetFromMatrix(countData = counts_excl, colData = coldata_excl, design = ~group)
dds_excl <- DESeq(dds_excl)
res_ir2_excl <- results(dds_excl, contrast = c("group", "IR2", "WT"), alpha = 0.05)
res_ir2_excl_df <- as.data.frame(res_ir2_excl)
res_ir2_excl_df$gene_id <- rownames(res_ir2_excl_df)
write.csv(res_ir2_excl_df, file.path(outdir, "GSE140771_IR2_vs_WT_excl_SRR10507806_full.csv"), row.names = FALSE)

# B56 + broader PP2A panel subset for direct comparison against the primary result
targets <- c("PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E",
             "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B", "SET")

# gene_name is not in this matrix (Ensembl IDs only); if a gene_id-to-symbol
# map is available (e.g. from the GENCODE GTF used originally), merge it here.
# As a minimal fallback, PPP2R5A's Ensembl ID is provided directly so the
# primary gene of interest can always be checked even without the GTF:
ppp2r5a_id <- rownames(res_ir2_excl_df)[grepl("^ENSG00000066027", rownames(res_ir2_excl_df))]
if (length(ppp2r5a_id) == 1) {
  cat("PPP2R5A (ENSG00000066027), IR2 vs WT, SRR10507806 excluded:\n")
  print(res_ir2_excl_df[ppp2r5a_id, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")])
} else {
  cat("PPP2R5A Ensembl ID not found by exact match after filtering; check gene_id column directly.\n")
}

cat("\nDone. Wrote blind=TRUE PCA coordinates and the SRR10507806-excluded\n",
    "IR2 vs WT full result table to:", outdir, "\n",
    "Compare GSE140771_IR2_vs_WT_excl_SRR10507806_full.csv against\n",
    "DESeq2_featureCounts_IR2_vs_WT_full.csv (the primary, all-sample result)\n",
    "for PPP2R5A and the rest of the B56/PP2A panel.\n")
