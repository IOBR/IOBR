# Calculate Exact P-Value for Correlation

Computes the exact p-value for the correlation between two numeric
variables using a specified correlation method.

## Usage

``` r
exact_pvalue(x, y, method)
```

## Arguments

- x:

  Numeric vector representing the first variable.

- y:

  Numeric vector representing the second variable.

- method:

  Character string specifying the correlation method: \`"spearman"\`,
  \`"pearson"\`, or \`"kendall"\`.

## Value

Numeric value representing the exact p-value.

## Author

Dongqiang Zeng

## Examples

``` r
sig_stad <- load_data("sig_stad")
#> ℹ Trying mirror 1/11: <https://github.com>
#> ✔ Download complete: "sig_stad"
p_val <- exact_pvalue(
  x = sig_stad$CD8.T.cells,
  y = sig_stad$CD_8_T_effector,
  method = "spearman"
)
print(p_val)
#> [1] 7.285891e-48
```
