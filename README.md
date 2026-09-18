[![DOI](https://shields.io)](https://doi.org)

# GEO2R Expression Analysis: GSE130404 & GSE12211

This repository contains the R programming workflows, scripts, and processed transcriptomic datasets used to perform differential gene expression analysis using the limma package. The analysis includes covariate adjustments and sensitivity testing across two Gene Expression Omnibus (GEO) datasets. 

### Project Structure

text

├── GSE12211_CEL/                     # Raw Affymetrix CEL files for GSE12211
├── results/                          # Output directory for primary analysis results
├── results_gse12211/                 # Output directory for sensitivity analysis results
├── GSE12211_CEL_files.zip            # Archived raw CEL data
├── GSE130404_sample_metadata_derived.csv # Phenotypic and covariate metadata for GSE130404
├── GSE130404_series_matrix.txt       # Processed expression data matrix (Tracked via Git LFS)
├── run_limma_GSE130404_covariate_adjusted.R # Main differential expression script with covariate handling
└── run_rma_limma_GSE12211_sensitivity.R     # Robust Multi-array Average (RMA) & sensitivity analysis script

Use code with caution.

### Datasets

1. **GSE130404** 

  * **Type:** Expression profiling by array.
  * **Data Input:** Processed series matrix file (GSE130404_series_matrix.txt).
  * **Analysis:** Linear modeling with covariate adjustment to evaluate differential expression while controlling for confounding variables.
2. **GSE12211** 

  * **Type:** Expression profiling by array (Affymetrix).
  * **Data Input:** Raw .CEL files.
  * **Analysis:** Background correction, normalization, and summarization via the RMA algorithm followed by a downstream sensitivity assessment.

### Prerequisites & Installation

### Git Large File Storage (LFS)

The dataset GSE130404_series_matrix.txt is larger than 50MB and is tracked using Git LFS. To clone this repository with the complete data tracking history, ensure you have Git LFS installed on your system before cloning: 

bash

# Install Git LFS (macOS example via Homebrew)
brew install git-lfs

# Clone the repository
git clone https://github.com/princekobby5/GEO2R-Analysis.git
cd GEO2R-Analysis

# Pull LFS files
git lfs pull

Use code with caution.

### R Dependencies

To execute the analysis scripts, ensure you have R installed along with the following Bioconductor and CRAN packages: 

R

if (!requireNamespace("BiocManager", quiet = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "affy", "Biobase"))
install.packages(c("tidyverse", "data.table"))

Use code with caution.

### Usage

### 1. Main Covariate-Adjusted Analysis (GSE130404)

To execute the linear modeling and extract the differentially expressed genes while accounting for target covariates: 

bash

Rscript run_limma_GSE130404_covariate_adjusted.R

Use code with caution.

### 2. Sensitivity and RMA Normalization (GSE12211)

To process the raw Affymetrix CEL microarrays from scratch and evaluate statistical sensitivity: 

bash

Rscript run_rma_limma_GSE12211_sensitivity.R

Use code with caution.

### Reproducibility and Citation

If you utilize this code or resource in your academic work, please cite the associated manuscript: 

*[Insert your publication title, authors, journal, and year here once available]* 

For archival permanence, a specific snapshot of this release is registered with Zenodo under DOI: **[Insert Zenodo DOI here if applicable]**. 

### License

This project is licensed under the MIT License - see the [LICENSE](/url?sa=i&source=web&rct=j&url=LICENSE&ved=2ahUKEwjX-pu2l_eWAxWfNIYAHQxYA6YQg5wRegYIAAgREF0&opi=89978449&cd&psig=AOvVaw1MDIUMGMGk9qDh33mPlfr_&ust=1789788380146000) file for details.
