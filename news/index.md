# Changelog

## IOBR 2.1.1

### New Features

- Added Docker support with automated image builds and publishing to
  GitHub Container Registry (`ghcr.io/iobr/iobr`). The image is based on
  `rocker/tidyverse` and includes IOBR with all dependencies
  pre-installed (#Docker).
- Added Conda/Mamba installation instructions to README for easier
  environment setup.
- Expanded
  [`iobr_cor_plot()`](https://iobr.github.io/IOBR/reference/iobr_cor_plot.md)
  to support multi-group (3+) comparisons using Kruskal-Wallis test, in
  addition to the existing two-group Wilcoxon test
  ([\#28](https://github.com/IOBR/IOBR/issues/28)).

### Bug Fixes

- Fixed
  [`generateRef_DEseq2()`](https://iobr.github.io/IOBR/reference/generateRef_DEseq2.md)
  to handle sparse single-cell RNA-seq data by using
  `type = "poscounts"` for size factor estimation, preventing errors
  when all genes contain at least one zero.
- Fixed
  [`batch_kruskal()`](https://iobr.github.io/IOBR/reference/batch_kruskal.md)
  output format to correctly display mean-centered values for all
  groups.
- Fixed
  [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md)
  computation results to align with previous correct versions.
- Updated
  [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md)
  documentation to clarify that gene identifiers are converted to gene
  symbols in the output regardless of input ID type.
- Wrapped some long-running examples in `\dontrun{}` and `\donttest{}`
  to comply with CRAN check time limits.

### Improvements

- Refactored internal helper functions in `count2tpm.R` for better code
  organization and maintainability.
- Updated
  [`load_data()`](https://iobr.github.io/IOBR/reference/load_data.md)
  documentation to include additional available datasets.

## IOBR 2.1.0

- Initial CRAN submission.
