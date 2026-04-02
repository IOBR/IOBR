# Core Algorithm for CIBERSORT Deconvolution

Performs nu-regression using support vector machines (SVM) to estimate
cell type proportions. This is the core computational engine of
CIBERSORT, using nu-SVR with linear kernel to decompose mixed gene
expression signals.

## Usage

``` r
CoreAlg(X, y, absolute, abs_method)
```

## Arguments

- X:

  Matrix or data frame containing signature matrix (predictor
  variables). Rows are genes, columns are cell types.

- y:

  Numeric vector containing the mixture sample expression (response
  variable).

- absolute:

  Logical indicating whether to use absolute space for weights. Default
  is FALSE (relative proportions).

- abs_method:

  String specifying the method for absolute space weights:
  \`"sig.score"\` or \`"no.sumto1"\`.

## Value

List containing:

- w:

  Estimated cell type weights/proportions

- mix_rmse:

  Root mean squared error of the deconvolution

- mix_r:

  Correlation coefficient between observed and predicted mixture

## Examples

``` r
X <- matrix(rnorm(100), nrow = 10)
y <- rnorm(10)
result <- CoreAlg(X, y, absolute = FALSE, abs_method = "sig.score")
```
