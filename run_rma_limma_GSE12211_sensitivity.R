suppressPackageStartupMessages({
  library(affy)
  library(limma)
})

# GSE12211 RMA plus limma sensitivity analysis.
# The manuscript's primary GSE12211 result uses MAS5 signal values extracted
# from the public CHP files, paired within patient with a two sided t test
# at each of 22,277 probe positions. This script instead reprocesses the raw
# CEL files with RMA, a modern background correction, normalization, and
# summarization method, and fits a paired limma model with patient as a
# blocking factor and empirical Bayes moderation, which is more stable than
# unmoderated probe wise t tests at n = 6. The two approaches are compared
# side by side for the ten B56 probes.
#
# Usage:
#   Rscript run_rma_limma_GSE12211_sensitivity.R CEL_DIR output_dir
#
# CEL_DIR should contain the 12 CEL files from GSE12211_RAW.tar, downloaded
# directly from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12211
# (Download family > GSE12211_RAW.tar), extracted so the .CEL files sit
# directly in CEL_DIR.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript run_rma_limma_GSE12211_sensitivity.R CEL_DIR output_dir")
}
cel_dir <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Sample pairing, same six patients as the MAS5 paired analysis ---------
sample_pairing <- data.frame(
  patient = paste0("Patient_", 1:6),
  pretreatment_GSM = c("GSM307000", "GSM307012", "GSM307013", "GSM307018", "GSM307020", "GSM307023"),
  day7_GSM = c("GSM307028", "GSM307030", "GSM307031", "GSM307032", "GSM307053", "GSM307072"),
  stringsAsFactors = FALSE
)

pheno <- rbind(
  data.frame(GSM = sample_pairing$pretreatment_GSM, patient = sample_pairing$patient, timepoint = "pretreatment"),
  data.frame(GSM = sample_pairing$day7_GSM, patient = sample_pairing$patient, timepoint = "day7")
)

cel_files <- list.files(cel_dir, pattern = "(\\.|_)CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
if (length(cel_files) == 0) stop("No .CEL files found in ", cel_dir)

# Match each CEL file to its GSM accession, however the filename is decorated
# (GEO CEL files are typically named GSMxxxxxxx_something.CEL.gz).
matched_gsm <- sapply(basename(cel_files), function(fn) {
  hit <- pheno$GSM[sapply(pheno$GSM, function(g) grepl(g, fn, fixed = TRUE))]
  if (length(hit) == 1) hit else NA
})
keep <- !is.na(matched_gsm)
cel_files <- cel_files[keep]
matched_gsm <- matched_gsm[keep]
if (length(cel_files) != 12) {
  warning(sprintf("Expected 12 matched CEL files for the six patient pairs, found %d. Check file naming.", length(cel_files)))
}

pheno <- pheno[match(matched_gsm, pheno$GSM), ]
rownames(pheno) <- basename(cel_files)

# --- RMA background correction, normalization, summarization ---------------
raw <- ReadAffy(filenames = cel_files, phenoData = AnnotatedDataFrame(pheno))
eset <- rma(raw)
expr <- exprs(eset)
colnames(expr) <- pheno$GSM[match(colnames(expr), rownames(pheno))]

write.csv(expr, file.path(outdir, "GSE12211_RMA_expression_matrix.csv"))

# --- Paired limma model, patient blocked ------------------------------------
patient <- factor(pheno$patient)
timepoint <- factor(pheno$timepoint, levels = c("pretreatment", "day7"))
design <- model.matrix(~ patient + timepoint)

fit <- lmFit(expr, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "timepointday7", number = Inf, sort.by = "none")
res$Probe_ID <- rownames(res)
write.csv(res, file.path(outdir, "GSE12211_RMA_limma_paired_full.csv"), row.names = FALSE)

# --- B56 probe subset and comparison against the MAS5 paired result --------
b56_probes <- c(
  "202186_x_at" = "PPP2R5A", "202187_s_at" = "PPP2R5A",
  "204611_s_at" = "PPP2R5B", "635_s_at" = "PPP2R5B",
  "201877_s_at" = "PPP2R5C", "213305_s_at" = "PPP2R5C", "214083_at" = "PPP2R5C",
  "202513_s_at" = "PPP2R5D", "211159_s_at" = "PPP2R5D",
  "203338_at" = "PPP2R5E"
)

b56_res <- res[res$Probe_ID %in% names(b56_probes), c("Probe_ID", "logFC", "P.Value", "adj.P.Val")]
b56_res$gene <- b56_probes[b56_res$Probe_ID]
colnames(b56_res)[2:4] <- paste0(colnames(b56_res)[2:4], "_RMA_limma")
b56_res <- b56_res[order(b56_res$gene, b56_res$Probe_ID), ]
write.csv(b56_res, file.path(outdir, "GSE12211_B56_RMA_limma_results.csv"), row.names = FALSE)

cat("Done. Wrote RMA expression matrix, full paired limma results, and the B56\n",
    "subset to:", outdir, "\n",
    "Compare GSE12211_B56_RMA_limma_results.csv against the existing\n",
    "GSE12211_B56_paired_gene_summary.csv (MAS5, unmoderated paired t test) to\n",
    "see whether RMA processing and empirical Bayes moderation change the\n",
    "B56 conclusion.\n")
