# Calculate Signature Score Using Integration Method

Computes signature scores using PCA, z-score, and ssGSEA methods
combined.

## Usage

``` r
calculate_sig_score_integration(
  pdata = NULL,
  eset,
  signature,
  mini_gene_count = 2,
  column_of_sample = "ID",
  adjust_eset = FALSE,
  parallel.size = 1L
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

- parallel.size:

  Number of parallel workers. Default is 1.

## Value

Tibble with signature scores from all three methods.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(eset) <- paste0("Gene", 1:100)
colnames(eset) <- paste0("Sample", 1:10)
signature <- list(
  Signature1 = paste0("Gene", 1:15),
  Signature2 = paste0("Gene", 16:30)
)
result <- calculate_sig_score_integration(eset = eset, signature = signature)
#> ℹ Calculating signature scores using PCA, z-score, and ssGSEA methods
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Step 1/3: PCA method
#> ℹ Step 2/3: z-score method
#> ℹ Step 3/3: ssGSEA method
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘HDF5Array’
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
#> ℹ GSVA version 2.4.8
#> ℹ Searching for rows with constant values
#> ℹ Calculating ssGSEA scores for 2 gene sets
#> ℹ Calculating ranks
#> ℹ Calculating rank weights
#> ℹ Normalizing ssGSEA scores
#> ✔ Calculations finished
```
