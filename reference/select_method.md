# Select a Signature Scoring Method Subset

Filters an integrated signature score matrix to retain results from a
specified method (PCA, ssGSEA, or zscore) and strips method suffixes
from column names.

## Usage

``` r
select_method(data, method = c("ssGSEA", "PCA", "zscore"))
```

## Arguments

- data:

  Data frame or matrix. Integrated signature score matrix.

- method:

  Character. One of "PCA", "ssGSEA", or "zscore" (case-insensitive).
  Default is "ssGSEA".

## Value

Matrix or data frame containing only the selected method's scores.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate data with multiple method columns
set.seed(123)
sim_res <- data.frame(
  ID = paste0("Sample", 1:10),
  Sig1_PCA = rnorm(10),
  Sig1_ssGSEA = rnorm(10),
  Sig1_zscore = rnorm(10),
  Sig2_PCA = rnorm(10),
  Sig2_ssGSEA = rnorm(10),
  Sig2_zscore = rnorm(10)
)
# Select PCA method columns only
pca_res <- select_method(sim_res, method = "PCA")
print(colnames(pca_res))
#> [1] "ID"   "Sig1" "Sig2"
```
