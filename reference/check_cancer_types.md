# Process Batch Table and Validate Cancer Types

Processes input data containing cancer types and validates each category
against a predefined list of supported cancer types
(\`timer_available_cancers\`).

## Usage

``` r
check_cancer_types(args)
```

## Arguments

- args:

  A list containing input parameters:

  batch

  :   Character. Path to a CSV file (optional).

  expression

  :   Character vector of expression identifiers (used if \`batch\` is
      NULL).

  category

  :   Character vector of cancer types (used if \`batch\` is NULL).

## Value

A character matrix with two columns:

- Column 1:

  Expression identifiers

- Column 2:

  Cancer categories

## Examples

``` r
# \donttest{
# Batch mode
tf <- tempfile(fileext = ".csv")
write.table(data.frame("exp1", "luad", "exp2", "brca"),
  file = tf, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
)
args <- list(batch = tf)
result <- check_cancer_types(args)
#> ℹ Enter batch mode

# Direct input mode
args <- list(
  expression = c("exp1", "exp2"),
  category = c("luad", "brca"),
  batch = NULL
)
result <- check_cancer_types(args)
# }
```
