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
# Load TCGA-STAD expression data
eset_stad <- load_data("eset_stad")

# Convert counts to TPM
eset <- count2tpm(eset_stad, idType = "ensembl")
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> Warning: longer object length is not a multiple of shorter object length
#> ℹ No duplicate gene symbols found.

# Check expression set integrity
check_eset(eset)

# Check with detailed output
check_eset(eset, print_result = TRUE, estimate_sd = TRUE)
#> ℹ Checking for NA values: 0 found
#> ℹ Checking for -Inf values: 0 found
#> ℹ Checking for +Inf values: 0 found
#> ℹ Features with zero variance: 7768
#> Warning: 7768 features have zero variance.
#> ℹ Zero-variance features may affect score calculation.
#> • Set `adjust_eset = TRUE` to remove them automatically.
```
