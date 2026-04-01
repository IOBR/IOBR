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
# \donttest{
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
signature_collection <- load_data("signature_collection")
result <- calculate_sig_score_integration(eset = eset, signature = signature_collection[1:4])
#> ℹ Calculating signature scores using PCA, z-score, and ssGSEA methods
#> ✔ Applied log2 transformation
#> ℹ Step 1/3: PCA method
#> ℹ Step 2/3: z-score method
#> ℹ Step 3/3: ssGSEA method
#> ℹ GSVA version 2.4.8
#> ℹ Searching for rows with constant values
#> ℹ Calculating GSVA ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores
#> ✔ Calculations finished
# }
```
