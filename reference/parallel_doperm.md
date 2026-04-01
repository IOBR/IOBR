# Parallel Permutation Test for CIBERSORT

Parallel version of doPerm. Performs permutation-based sampling and runs
the CoreAlg function iteratively using multiple CPU cores to accelerate
computation. This function generates an empirical null distribution of
correlation coefficients for p-value calculation in CIBERSORT analysis.

## Usage

``` r
parallel_doperm(
  perm1,
  X1,
  Y1,
  absolute1,
  abs_method1,
  num_cores1 = 2,
  seed = NULL
)
```

## Arguments

- perm1:

  Integer. Number of permutations to perform (≥100 recommended).

- X1:

  Matrix or data frame. Signature matrix (cell type GEP barcode).

- Y1:

  Matrix. Mixture file containing gene expression profiles.

- absolute1:

  Logical. Whether to run in absolute mode (default: FALSE).

- abs_method1:

  String. Method for absolute mode: 'sig.score' or 'no.sumto1'.

- num_cores1:

  Integer. Number of CPU cores for parallel computation (default: 2).

- seed:

  Integer. Random seed for reproducibility. If NULL (default), uses
  current random state. Set to a specific value (e.g., 123) for
  reproducible results across runs.

## Value

List containing:

- dist:

  Numeric vector of correlation coefficients from permutations,
  representing the empirical null distribution.

## Details

This function utilizes the `foreach` and `doParallel` packages to
distribute permutation iterations across multiple cores. It
automatically handles cluster setup/teardown via
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) to prevent resource
leaks.

Note: Windows users may experience slower performance due to
socket-based parallelization (PSOCK) versus forking on Unix systems.

## See also

[`doPerm`](https://iobr.github.io/IOBR/reference/doPerm.md) for the
sequential version,
[`CoreAlg`](https://iobr.github.io/IOBR/reference/CoreAlg.md),
[`CIBERSORT`](https://iobr.github.io/IOBR/reference/CIBERSORT.md)

## Examples

``` r
# \donttest{
X <- matrix(rnorm(1000), nrow = 100)
Y <- matrix(rnorm(500), nrow = 100)
rownames(X) <- rownames(Y) <- paste0("Gene", 1:100)

result <- parallel_doperm(
  perm1 = 100, X1 = X, Y1 = Y,
  absolute1 = FALSE, abs_method1 = "sig.score", num_cores1 = 2
)
#> Warning: already exporting variable(s): perm_indices, Y1, X1, absolute1, abs_method1
str(result$dist)
#>  num [1:100, 1] 0.1036 0.2187 0.1911 0.2319 0.0303 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:100] "result.1" "result.2" "result.3" "result.4" ...
#>   ..$ : NULL
# }
```
