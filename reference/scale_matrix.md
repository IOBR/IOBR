# Scale and Manipulate a Matrix

Scales a gene expression matrix, optionally applies logarithmic
transformation, and performs feature manipulation to remove problematic
variables (e.g., genes with zero variance or missing values).

## Usage

``` r
scale_matrix(matrix, log2matrix = TRUE, manipulate = TRUE)
```

## Arguments

- matrix:

  Numeric matrix with genes as rows and samples as columns.

- log2matrix:

  Logical indicating whether to apply log2 transformation using
  \[log2eset()\]. Default is \`TRUE\`.

- manipulate:

  Logical indicating whether to perform feature manipulation to remove
  problematic features. Default is \`TRUE\`.

## Value

A scaled matrix (genes as rows, samples as columns).

## Author

Dongqiang Zeng

## Examples

``` r
eset_gse62254 <- load_data("eset_gse62254")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "eset_gse62254"
eset2 <- scale_matrix(eset_gse62254, log2matrix = FALSE, manipulate = TRUE)
#> ✔ Retained 54675 of 54675 features
#> ℹ Retained 54675 features after manipulation
```
