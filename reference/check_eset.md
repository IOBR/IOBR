# Check Integrity and Outliers of Expression Set

Performs quality checks on an expression matrix to identify missing
values, infinite values, and features with zero variance. Issues
warnings when potential problems are detected that may affect downstream
analyses.

## Usage

``` r
check_eset(eset, print_result = FALSE, estimate_sd = FALSE)
```

## Arguments

- eset:

  Expression matrix or data frame with genes/features in rows and
  samples in columns.

- print_result:

  Logical indicating whether to print detailed check results to the
  console. Default is \`FALSE\`.

- estimate_sd:

  Logical indicating whether to check for features with zero standard
  deviation. Default is \`FALSE\`.

## Value

Invisibly returns \`NULL\`. Function is called for its side effects
(printing messages and issuing warnings).

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(eset) <- paste0("Gene", 1:100)
colnames(eset) <- paste0("Sample", 1:10)
check_eset(eset, print_result = TRUE)
#> ℹ Checking for NA values: 0 found
#> ℹ Checking for -Inf values: 0 found
#> ℹ Checking for +Inf values: 0 found
```
