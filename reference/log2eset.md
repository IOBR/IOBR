# Log2 Transformation of Gene Expression Matrix

Determines whether a gene expression matrix requires log2 transformation
based on the distribution of values, and applies it if necessary. This
is useful for automatically detecting raw counts or linear-scale data
that should be log-transformed for downstream analysis.

## Usage

``` r
log2eset(eset)
```

## Arguments

- eset:

  Numeric matrix. Gene expression data with genes as rows and samples in
  columns.

## Value

Numeric matrix. Log2-transformed gene expression data (if transformation
was needed), or the original data otherwise.

## Examples

``` r
set.seed(123)
eset <- matrix(rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 10)
rownames(eset) <- paste0("Gene", 1:100)
colnames(eset) <- paste0("Sample", 1:10)
eset_transformed <- log2eset(eset)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
```
