# IOBR 2.2.0

# IOBR 2.1.1

## New Features

* Added Docker support with automated image builds and publishing to GitHub Container Registry (`ghcr.io/iobr/iobr`). The image is based on `rocker/tidyverse` and includes IOBR with all dependencies pre-installed (#Docker).
* Added Conda/Mamba installation instructions to README for easier environment setup.
* Expanded `iobr_cor_plot()` to support multi-group (3+) comparisons using Kruskal-Wallis test, in addition to the existing two-group Wilcoxon test (#28).

## Bug Fixes

* Fixed `generateRef_DEseq2()` to handle sparse single-cell RNA-seq data by using `type = "poscounts"` for size factor estimation, preventing errors when all genes contain at least one zero.
* Fixed `batch_kruskal()` output format to correctly display mean-centered values for all groups.
* Fixed `count2tpm()` computation results to align with previous correct versions.
* Updated `count2tpm()` documentation to clarify that gene identifiers are converted to gene symbols in the output regardless of input ID type.
* Wrapped some long-running examples in `\dontrun{}` and `\donttest{}` to comply with CRAN check time limits.

## Improvements

* Refactored internal helper functions in `count2tpm.R` for better code organization and maintainability.
* Updated `load_data()` documentation to include additional available datasets.

# IOBR 2.1.0

* Initial CRAN submission.
