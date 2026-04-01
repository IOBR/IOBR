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
imvigor210_sig <- load_data("imvigor210_sig")
input <- remove_names(
  imvigor210_sig,
  variable = "colnames",
  patterns_to_na = patterns_to_na,
  patterns_space = NULL
)
```
