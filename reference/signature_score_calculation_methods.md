# Signature Score Calculation Methods

A named vector of supported methods for calculating signature scores.

## Usage

``` r
signature_score_calculation_methods
```

## Format

Named character vector:

- PCA:

  Principal Component Analysis method ("pca")

- ssGSEA:

  Single-sample Gene Set Enrichment Analysis ("ssgsea")

- z-score:

  Z-score transformation method ("zscore")

- Integration:

  Integration of multiple methods ("integration")

## Examples

``` r
signature_score_calculation_methods
#>           PCA        ssGSEA       z-score   Integration 
#>         "pca"      "ssgsea"      "zscore" "integration" 
signature_score_calculation_methods["PCA"]
#>   PCA 
#> "pca" 
```
