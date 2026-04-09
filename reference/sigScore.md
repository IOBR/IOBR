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
# Load example data
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
eset <- count2tpm(eset_stad, idType = "ensembl")
#> ℹ Loading cached data: "anno_grch38"
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> ℹ Found 1679 duplicate symbols. Using "mean" for ranking.
#> ✔ Reduced to 54658 unique genes

# Get signature genes
signature_tme <- load_data("signature_tme")
genes <- signature_tme[["CD_8_T_effector"]]
genes <- genes[genes %in% rownames(eset)]

# Calculate scores (only if enough genes are available)
if (length(genes) >= 2) {
  score_pca <- sigScore(eset = eset[genes, ], methods = "PCA")
  score_mean <- sigScore(eset = eset[genes, ], methods = "mean")
  score_zscore <- sigScore(eset = eset[genes, ], methods = "zscore")
}
```
