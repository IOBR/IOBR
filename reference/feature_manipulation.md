# Feature Quality Control and Filtering

Filters features (variables) in a matrix or data frame by removing those
with missing values, non-numeric types, infinite values, or zero
variance. This is useful for preparing data for downstream statistical
analyses.

## Usage

``` r
feature_manipulation(
  data,
  feature = NULL,
  is_matrix = FALSE,
  print_result = FALSE
)
```

## Arguments

- data:

  A matrix or data frame containing features to filter.

- feature:

  Character vector of feature names to check. If \`is_matrix = TRUE\`,
  features are extracted from row names of the matrix.

- is_matrix:

  Logical indicating whether \`data\` is a gene expression matrix
  (features as rows, samples as columns). If \`TRUE\`, the matrix is
  transposed for processing. Default is \`FALSE\`.

- print_result:

  Logical indicating whether to print filtering statistics. Default is
  \`FALSE\`.

## Value

Character vector of feature names that pass all quality checks.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
feas <- feature_manipulation(
  data = eset_stad,
  feature = rownames(eset_stad),
  is_matrix = TRUE,
  print_result = TRUE
)
#> ℹ Removing 9060 features with zero variance
#> ✔ Retained 51423 of 60483 features
# }
```
