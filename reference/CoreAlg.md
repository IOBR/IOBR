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
# Simulate data
set.seed(123)
X <- matrix(rnorm(100), nrow = 10)
rownames(X) <- paste0("Gene", 1:10)
colnames(X) <- paste0("Cell", 1:10)
y <- rnorm(10)
names(y) <- paste0("Gene", 1:10)

# Run core algorithm
result <- CoreAlg(X, y, absolute = FALSE, abs_method = "sig.score")
if (!is.null(result)) str(result)
#> List of 3
#>  $ w       : num [1, 1:10] 0.373 0.2 0.301 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:10] "Cell1" "Cell2" "Cell3" "Cell4" ...
#>  $ mix_rmse: num 0.651
#>  $ mix_r   : num 0.582
```
