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
set.seed(123)
test_data <- data.frame(
  A = c(1, 2, 3),
  B = c(1, 1, 1), # zero variance
  C = c(1, NA, 3), # missing value
  D = c("a", "b", "c") # non-numeric
)
feas <- feature_manipulation(data = test_data, feature = colnames(test_data), print_result = TRUE)
#> ℹ Removing 1 feature with NA values
#> ℹ Removing 1 non-numeric feature
#> ℹ Removing 1 feature with zero variance
#> ✔ Retained 1 of 4 features
print(feas)
#> [1] "A"
```
