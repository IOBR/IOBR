# Calculate Signature Score Using ssGSEA Method

Computes signature scores using single-sample Gene Set Enrichment
Analysis.

## Usage

``` r
calculate_sig_score_ssgsea(
  pdata = NULL,
  eset,
  signature,
  mini_gene_count = 5,
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

Tibble with signature scores.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
eset <- count2tpm(eset_stad, idType = "ensembl")
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> Warning: longer object length is not a multiple of shorter object length
#> ℹ No duplicate gene symbols found.
signature_tme <- load_data("signature_tme")
result <- calculate_sig_score_ssgsea(eset = eset, signature = signature_tme)
#> ℹ Calculating signature scores using ssGSEA method
#> Warning: No signatures have enough genes (min: 5). Returning pdata unchanged.
# }
```
