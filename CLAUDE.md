# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

IOBR (Immune Oncology Biological Research) is an R package for tumor microenvironment (TME) analysis. It provides:

- **TME Deconvolution**: 11 methods (CIBERSORT, TIMER, xCell, MCPcounter, ESTIMATE, EPIC, IPS, quanTIseq, SVR, lsei)
- **Signature Scoring**: PCA, z-score, ssGSEA methods for 322+ curated gene signatures
- **Visualization**: Survival plots, correlation heatmaps, box plots, forest plots
- **Data Preprocessing**: Batch correction, normalization, format conversion

## Development Commands

### Building

```bash
# Quick build (without vignettes) - for development
R CMD build . --no-build-vignettes

# Full build (with vignettes) - for release
R CMD build .
```

### Checking

```bash
# Quick check (skip suggested packages)
_R_CHECK_FORCE_SUGGESTS_=false R CMD check IOBR_*.tar.gz --no-manual

# Using devtools (recommended)
R -e "devtools::check()"

# Full check
R CMD check IOBR_*.tar.gz --no-manual
```

### Documentation

```bash
# Update roxygen2 documentation (use devtools)
R -e "devtools::document()"
```

### Code Quality

```bash
# Format all R code
R -e "styler::style_pkg()"

# Check code style
R -e "lintr::lint_package()"
```

### Installation

```bash
# Install dependencies
R --slave -e "pak::local_install_deps()"

# Install locally
R --slave -e "pak::local_install('.')"
# or
R -e "devtools::install()"
```

**Note**: This package does not have a `tests/` directory (it's excluded in `.Rbuildignore`). Testing is done via examples in roxygen documentation and vignettes.

## Architecture

### Core Modules

| Module | Key Files | Purpose |
|--------|-----------|---------|
| TME Deconvolution | `R/deconvo_tme.R`, `R/CIBERSORT.R`, `R/timer.R`, `R/xCell.R`, `R/mcpcounter.R`, `R/EPIC_helper.R`, `R/quantiseq.R`, `R/estimate_helper.R`, `R/IPS_helper.R` | Cell fraction estimation from expression data |
| Signature Scoring | `R/sigScore.R`, `R/calculate_sig_score.R` | Calculate signature scores (PCA, z-score, ssGSEA) |
| Visualization | `R/sig_heatmap.R`, `R/sig_pheatmap.R`, `R/iobr_cor_plot.R`, `R/sig_surv_plot.R`, `R/batch_surv.R`, `R/sig_box.R`, `R/sig_forest.R` | Plotting functions for TME results |
| Data Preprocessing | `R/remove_batcheffect.R`, `R/count2tpm.R`, `R/check_eset.R`, `R/anno_eset.R`, `R/transform_data.R` | Expression data preparation |
| Statistical Analysis | `R/batch_cor.R`, `R/batch_wilcoxon.R`, `R/batch_kruskal.R`, `R/get_cor_matrix.R` | Statistical testing utilities |
| Prognostic Models | `R/PrognosticModel.R`, `R/BinomialModel.R`, `R/add_riskscore.R`, `R/roc_time.R` | Survival analysis and risk modeling |

### Key Data Structures

- **Expression sets (`eset`)**: Matrix/data.frame with genes as rownames, samples as columns. Expects HGNC gene symbols (not Ensembl IDs).
- **Signature collection**: Lists of gene vectors (e.g., `signature_collection` with 323 signatures)
- **Deconvolution results**: Data frames with `ID` column + cell fraction columns (suffixed with method name)

### Important Internal Data

Stored in `R/sysdata.rda` (internal) and `data/` (exported):

- `lm22`: CIBERSORT reference matrix (22 immune cell types)
- `signature_collection`: 323 curated gene signatures with citations
- `xCell.data`, `mcp_genes`, `quantiseq_data`: Method-specific reference data
- `TRef`/`BRef`: EPIC reference matrices

### Package Initialization

- `R/zzz.R`: `.onLoad()` checks required packages; `.onAttach()` prints citation info
- `R/globalVariables.R`: Suppresses R CMD check notes for tidyverse variables
- `R/imports.R`: Centralized `@importFrom` declarations

## Coding Standards

- **Style**: tidyverse style guide, 80-character line limit
- **Documentation**: roxygen2 with `@examples` for all exported functions
- **Dependencies**: Use specific packages (not tidyverse meta-package), declare in DESCRIPTION
- **Error handling**: Use `rlang::check_installed()` for optional dependencies
- **Data validation**: Input checks for Ensembl IDs, log-transformed data in deconvolution functions

## CI/CD

GitHub Actions workflow `.github/workflows/R-CMD-check.yaml` tests on:
- macOS (R release)
- Windows (R release)
- Ubuntu (R devel, R release)

## Common Patterns

### Adding a New Deconvolution Method

1. Create function `deconvo_<method>()` in `R/deconvo_tme.R` or new file
2. Follow naming convention: output columns suffixed with `_<METHOD>`
3. Add to `tme_deconvolution_methods` vector
4. Add switch case in `deconvo_tme()` dispatcher
5. Include `ID` column as first column in output

### Adding a New Signature

1. Add to `signature_collection` list in `data-raw/`
2. Update `signature_collection_citation` with reference
3. Run `R -e "devtools::document()"` to update documentation

### Data Access Pattern

**All data loading should use `load_data()` function:**

```r
# Load internal data (stored in R/sysdata.rda)
lm22 <- load_data("lm22")
xCell.data <- load_data("xCell.data")

# Load exported data (stored in data/)
eset_stad <- load_data("eset_stad")
signature_collection <- load_data("signature_collection")
```

The `load_data()` function handles both internal and exported datasets uniformly and provides helpful error messages if the dataset name is not found.

## Important Notes

- **No test suite**: Package relies on examples and vignettes for testing
- **Vignette build**: Requires pandoc and takes 5-10 minutes; skip during development with `--no-build-vignettes`
- **Bioconductor dependencies**: Some features require Bioconductor packages (GSVA, ComplexHeatmap)
- **CIBERSORT**: Requires separate license from Stanford for the reference matrix

## Resources

- Documentation: https://iobr.github.io/book/
- Issues: https://github.com/IOBR/IOBR/issues
- Paper: Cell Reports Methods (2024) <https://doi.org/10.1016/j.crmeth.2024.100910>
