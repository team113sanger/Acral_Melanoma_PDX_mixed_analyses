# Acral Melanoma PDX - Mixed Analyses

This repository contains analysis scripts, metadata, and results for the Acral Melanoma PDX study, including code used to generate figures for the article.

## Repository Overview

### Key Figures Generated
- **Figure 1b**: Overview of samples and biorepository information
- **Figure S3b-c**: Comparison of BCS/FCS copy number burden across tumor subtypes
- **Figure S3a**: Correlation of copy number gains and gene expression in hailstorms
- **Figure 4a**: SNV/indel and BCS/FCS correlations between AM-PDXs and patient tumors
- **Figure S4c**: Relative CNA and mutational signature plots
- **Figure 5a**: Heatmap showing an overview of the pharmacological screen performed in AM-PDX-derived cell lines

## Setup & Reproducibility

### R Environment
If you're interested in [reproducing the R environment using `renv`](https://rstudio.github.io/renv/reference/index.html) please follow the official documentation.

```R
# Install renv if needed
if (!require("renv")) install.packages("renv")

# Restore the project environment
renv::restore()
```

**Environment Details:**
- R version: 4.4.1
- Dependencies: See `renv.lock` in each analysis folder

### File Path Management
We used the `here` package for robust path handling:
- Automatically detects project root
- Creates platform-independent paths
- Seamless RStudio integration

## Data Sources

1. **Sequenza Output Data**:  
   Available at [Figshare](https://doi.org/10.6084/m9.figshare.29088173)

2. **Copy Number Scores (FCS/BCS)**:  
   Generated using [CNApp](https://tools.idibaps.org/CNApp/) with:
   - SEQUENZA segmentation data as input
   - Default parameters

## Additional Notes
- All scripts are contained in the `scripts/` directory
- Each analysis folder contains its specific dependencies
- Analysis from 06 script was created by Aguiar F.C., Carvalho D.G and Sousa Squiavinato A.C.M.

## ✉️ Contact
For questions about this repository or analyses:  
**Annie Cristhine Moraes Sousa Squiavinato**  
📧 as81@sanger.ac.uk