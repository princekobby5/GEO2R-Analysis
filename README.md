[![DOI](https://shields.io)](https://doi.org)

# GEO2R Expression Analysis: GSE130404 & GSE12211

This repository contains the R programming workflows, scripts, and processed transcriptomic datasets used to perform differential gene expression analysis using the limma package across two GEO datasets.

### Project Structure & Datasets
Includes workflows for **GSE130404** (covariate-adjusted linear modeling) and **GSE12211** (RMA normalization and sensitivity analysis using raw CEL files).

### Prerequisites
* **Git LFS:** Required for tracking large files like `GSE130404_series_matrix.txt`.
* **R Packages:** `limma`, `affy`, `Biobase`, `tidyverse`, and `data.table`.

### Usage
Run the analysis scripts via Terminal:
```bash
Rscript run_limma_GSE130404_covariate_adjusted.R
Rscript run_rma_limma_GSE12211_sensitivity.R
```

## License
MIT License
