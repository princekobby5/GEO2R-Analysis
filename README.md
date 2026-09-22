# PP2A B56 transcriptomic analysis workflows

This repository contains analysis code and input/output tables for GSE12211, GSE130404 and GSE140771, plus the deposited GSE120932 GEO2R top table and a probe-count-matched gene-set benchmark. The GSE120932 table supports the benchmark but does not contain sample-level expression values. It cannot be used to reproduce an expression-and-variance-matched benchmark or an array-level QC analysis.

### Project Structure & Datasets
Includes **GSE130404** covariate-adjusted modeling, **GSE12211** paired RMA/limma sensitivity analysis from CEL files, **GSE140771** RNA-seq PCA and exclusion analysis, and the **GSE120932** probe-count benchmark.

### Prerequisites
* **Git LFS:** Required for tracking large files like `GSE130404_series_matrix.txt`.
* **R packages:** `limma`, `affy`, `Biobase`, `GEOquery`, and `DESeq2` as required by the individual scripts.
* **Python packages for the GSE120932 benchmark:** `numpy`, `pandas`, and `scipy`.

### Usage
Run the analysis scripts via Terminal:

```bash
Rscript run_limma_GSE130404_covariate_adjusted.R GSE130404_series_matrix.txt GSE130404_sample_metadata_derived.csv results
Rscript run_rma_limma_GSE12211_sensitivity.R GSE12211_CEL results_gse12211
Rscript run_GSE140771_QC_sensitivity.R featureCounts_filtered_count_matrix.csv results_gse140771
python3 random_gene_benchmark_GSE120932.py
```

The benchmark script reads `GSE120932.top.table.tsv` and writes `GSE120932_B56_probe_results.csv`, `GSE120932_random_gene_benchmark_summary.csv`, and `GSE120932_random_gene_benchmark_distribution.csv` in the working directory. With seed 20260910 and 100,000 draws, the distributed input produces an empirical P value of approximately 0.0221. This is exploratory and matches probe count only; it does not account for baseline expression or variance.

For the GSE140771 exclusion check, removing SRR10507806 leaves only one IR2 sample. The similar PPP2R5A fold change shows that the retained IR2 sample supports the direction of the primary estimate. The fitted adjusted P value and standard error after exclusion cannot establish replicate-level robustness for IR2.

## License
MIT License
