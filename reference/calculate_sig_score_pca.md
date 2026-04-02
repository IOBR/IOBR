# Calculate Signature Score Using PCA Method

Computes signature scores using Principal Component Analysis. The first
principal component is used as the signature score.

## Usage

``` r
calculate_sig_score_pca(
  pdata = NULL,
  eset,
  signature,
  mini_gene_count = 3,
  column_of_sample = "ID",
  adjust_eset = FALSE
)
```

## Arguments

- pdata:

  Data frame with phenotype data. If \`NULL\`, created from \`eset\`
  column names.

- eset:

  Expression matrix (genes as rows, samples as columns).

- signature:

  List of gene signatures.

- mini_gene_count:

  Minimum genes required per signature. Default is 3.

- column_of_sample:

  Column in \`pdata\` with sample IDs. Default is \`"ID"\`.

- adjust_eset:

  Logical: remove problematic features. Default is \`FALSE\`.

## Value

Tibble with signature scores.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(eset) <- paste0("Gene", 1:100)
colnames(eset) <- paste0("Sample", 1:10)
signature <- list(
  Signature1 = paste0("Gene", 1:10),
  Signature2 = paste0("Gene", 11:20)
)
result <- calculate_sig_score_pca(eset = eset, signature = signature)
#> ℹ Calculating signature scores using PCA method
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Calculating scores for 2 signature(s)
```
