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
# Simulate data
set.seed(123)
x <- rnorm(100)
y <- rnorm(100)
p_val <- exact_pvalue(x = x, y = y, method = "spearman")
print(p_val)
#> [1] 0.9692781
```
