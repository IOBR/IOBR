# Permutation Test for CIBERSORT

Performs permutation-based sampling to generate an empirical null
distribution of correlation coefficients for p-value calculation in
CIBERSORT analysis. Randomly samples from the mixture data to create
null distributions.

## Usage

``` r
doPerm(perm, X, Y, absolute, abs_method, seed = NULL)
```

## Arguments

- perm:

  Integer. Number of permutations to perform (≥100 recommended for
  reliable p-value estimation).

- X:

  Matrix or data frame containing signature matrix (predictor
  variables).

- Y:

  Numeric vector containing the mixture sample expression.

- absolute:

  Logical indicating whether to use absolute space for weights.

- abs_method:

  String specifying the method for absolute space weights:
  \`"sig.score"\` or \`"no.sumto1"\`.

- seed:

  Integer. Random seed for reproducibility. If NULL (default), uses
  current random state.

## Value

List containing:

- dist:

  Numeric vector of correlation coefficients from permutations

## Examples

``` r
X <- matrix(rnorm(100), nrow = 10)
Y <- rnorm(10)
result <- doPerm(100, X, Y, absolute = FALSE, abs_method = "sig.score")
```
