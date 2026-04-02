# Remove Patterns from Column Names or Variables

Modifies column names or specified variables in a data frame by
replacing specified patterns with empty strings or spaces.

## Usage

``` r
remove_names(
  input_df,
  variable = "colnames",
  patterns_to_na = patterns_to_na,
  patterns_space = NULL
)
```

## Arguments

- input_df:

  Data frame. Input data to modify.

- variable:

  Character. Column to modify: "colnames" for column names, or a
  specific column name. Default is "colnames".

- patterns_to_na:

  Character vector. Patterns to replace with empty string. Default uses
  \[patterns_to_na\].

- patterns_space:

  Character vector or \`NULL\`. Patterns to replace with spaces. Default
  is \`NULL\`.

## Value

Modified data frame with patterns replaced.

## Author

Dongqiang Zeng

## Examples

``` r
df <- data.frame(
  "CellA_cibersort" = 1:5,
  "CellB_xCell" = 6:10,
  "CellC_TIMER" = 11:15
)
result <- remove_names(df, variable = "colnames", patterns_to_na = patterns_to_na)
colnames(result)
#> [1] "CellA"  "CellB_" "CellC" 
```
