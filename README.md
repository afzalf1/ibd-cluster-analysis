# IBD Cluster Analysis

This repository contains R scripts for analyzing associations between patient clusters, diagnoses, complications, and surgeries in an IBD (Inflammatory Bowel Disease) dataset. The analysis includes:

- Generating multi-way contingency tables.
- Calculating observed vs. expected ratios.
- Performing chi-square tests for independence.
- Computing adjusted residuals and FDR-corrected p-values for each cell.

## ⚠ Note on Data

The actual dataset used in this analysis is not included in this repository due to confidential information. To run the scripts, you must provide your own data in the same format as the original:

- Tab-separated values (`.tsv` file)
- Columns required:
  - `Cluster`
  - `Diagnosis`
  - `Crohn s disease phenotype`
  - `IBD surgery final`

## Usage

1. Clone this repository.
2. Install required packages:

```R
install.packages(c("readr", "MASS"))
