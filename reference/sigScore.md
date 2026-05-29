# Calculate Signature Score Using PCA, Mean, or Z-score Methods

Computes signature scores from gene expression data using Principal
Component Analysis (PCA), mean-based, or z-score approaches.

## Usage

``` r
sigScore(eset, methods = c("PCA", "mean", "zscore"))
```

## Arguments

- eset:

  Normalized expression matrix with genes (signature) as rows and
  samples as columns.

- methods:

  Scoring method: \`"PCA"\` (default) for principal component 1,
  \`"mean"\` for mean expression, or \`"zscore"\` for z-score normalized
  mean.

## Value

Numeric vector of length \`ncol(eset)\`; a score summarizing the rows of
\`eset\`.

## Author

Dorothee Nickles, Dongqiang Zeng

## Examples

``` r
# Create small example expression matrix
set.seed(123)
test_eset <- matrix(rnorm(1000), nrow = 10, ncol = 100)
rownames(test_eset) <- paste0("Gene", 1:10)
colnames(test_eset) <- paste0("Sample", 1:100)

# Calculate scores
score_pca <- sigScore(eset = test_eset, methods = "PCA")
score_mean <- sigScore(eset = test_eset, methods = "mean")
score_zscore <- sigScore(eset = test_eset, methods = "zscore")
```
