# PP2A B56 transcriptomic analysis workflows

This repository contains analysis code and input/output tables for GSE12211, GSE120932, GSE130404, and GSE140771. For GSE120932 it includes the complete public, preprocessed 12-array [series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120932/matrix/GSE120932_series_matrix.txt.gz), a six-array scripted limma model, array quality checks, and probe-count-only and expression/variance-matched exploratory gene benchmarks.

### Project Structure & Datasets
Includes **GSE130404** covariate-adjusted modeling, **GSE12211** paired RMA/limma sensitivity analysis from CEL files, **GSE140771** RNA-seq PCA and exclusion analysis, and **GSE120932** array-level reanalysis and exploratory benchmarks.

### Prerequisites
* **Git LFS:** Required for tracking large files like `GSE130404_series_matrix.txt`.
* **R packages:** `limma`, `affy`, `Biobase`, `GEOquery`, and `DESeq2` as required by the individual scripts.
* **Python packages for the GSE120932 benchmarks:** `numpy`, `pandas`, and `scipy` (the probe-count script uses SciPy).

### Usage
Run the analysis scripts via Terminal:

```bash
Rscript run_limma_GSE130404_covariate_adjusted.R GSE130404_series_matrix.txt GSE130404_sample_metadata_derived.csv results
Rscript run_rma_limma_GSE12211_sensitivity.R GSE12211_CEL results_gse12211
Rscript run_GSE140771_QC_sensitivity.R featureCounts_filtered_count_matrix.csv results_gse140771
Rscript run_limma_GSE120932_scripted.R GSE120932_series_matrix.txt.gz results_gse120932 GSE120932.top.table.tsv
python3 random_gene_benchmark_GSE120932.py
python3 run_GSE120932_matched_benchmark.py GSE120932_series_matrix.txt.gz GSE120932.top.table.tsv results_gse120932
python3 run_GSE120932_matched_benchmark.py GSE120932_series_matrix.txt.gz GSE120932.top.table.tsv results_gse120932 --variance-mode within-group
```

The benchmark script reads `GSE120932.top.table.tsv` and writes `GSE120932_B56_probe_results.csv`, `GSE120932_random_gene_benchmark_summary.csv`, and `GSE120932_random_gene_benchmark_distribution.csv` in the working directory. With seed 20260910 and 100,000 draws, the distributed input produces an empirical P value of approximately 0.0221. This is exploratory and matches probe count only; it does not account for baseline expression or variance.

The series matrix contains all 12 samples; the main contrast uses parental GSM3421746–GSM3421748 and K562-IR without imatinib GSM3421752–GSM3421754. The R script fits this contrast with limma, writes full results, sample medians and IQRs, correlations, PCA, and a QC plot. With R 4.3.3 and limma 3.58.1 the refitted 47,223 probe log fold changes match the GEO2R table to within 5 × 10⁻⁹; B56 adjusted P values match to the printed precision. The matrix was already log2 transformed and normalized by the source investigators; the script does not normalize it a second time.

The expression/variance benchmark matches each B56 gene to 200 non-B56 candidates with the same probe count using parental baseline log2 expression and log2 variance across the six arrays, standardized within probe-count stratum. Gene expression is the mean of its mapped probes in each array. With 100,000 five-gene sets and seed 20260922, the empirical P value is **0.04584**. Pool sizes of 100 and 500 give **0.06754** and **0.04055**. As a variance-definition sensitivity check, pooled residual variance within each group excludes the difference between parental and resistant means: its pool sizes of 100, 200, and 500 give **0.04561**, **0.04124**, and **0.04301**. The all-array variance results cross the conventional threshold as pool size changes; this unprespecified benchmark cannot establish B56-specific enrichment. Full candidate pools, distributions, and sensitivity results are in `results_gse120932/`. To reproduce a different pool size, pass `--pool-size 100` or `--pool-size 500` to the benchmark script with a separate output directory.

For the GSE140771 exclusion check, removing SRR10507806 leaves only one IR2 sample. The similar PPP2R5A fold change shows that the retained IR2 sample supports the direction of the primary estimate. The fitted adjusted P value and standard error after exclusion cannot establish replicate-level robustness for IR2.

The public [GSE140771 sample metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140771/matrix/GSE140771_series_matrix.txt.gz) list original and "Repeated" library entries. Distinct BioSample accessions alone cannot establish whether these were made from separate cultures and RNA extractions. DESeq2 statistics assume the six libraries can be treated as independent observations; see the manuscript limitation.

## License
MIT License
