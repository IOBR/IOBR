# Calculate Signature Score

Main interface for calculating signature scores from gene expression
data. Supports PCA, z-score, ssGSEA, and integration methods.

## Usage

``` r
calculate_sig_score(
  pdata = NULL,
  eset,
  signature = NULL,
  method = c("pca", "ssgsea", "zscore", "integration"),
  mini_gene_count = 3,
  column_of_sample = "ID",
  print_gene_proportion = FALSE,
  print_filtered_signatures = FALSE,
  adjust_eset = FALSE,
  parallel.size = 1L,
  ...
)
```

## Arguments

- pdata:

  Phenotype data. If \`NULL\`, created from \`eset\` column names.

- eset:

  Expression matrix (CPM, TPM, RPKM, FPKM, etc.).

- signature:

  List of gene signatures. Can also be a character string naming a
  built-in signature collection (e.g., \`"signature_collection"\`,
  \`"signature_tme"\`, \`"go_bp"\`, \`"kegg"\`, \`"hallmark"\`).

- method:

  Scoring method: \`"pca"\`, \`"ssgsea"\`, \`"zscore"\`, or
  \`"integration"\`. Default is \`"pca"\`.

- mini_gene_count:

  Minimum genes required per signature. Default is 3 (or 5 for ssGSEA).

- column_of_sample:

  Column with sample IDs in \`pdata\`. Default is \`"ID"\`.

- print_gene_proportion:

  Logical: print gene coverage. Default is \`FALSE\`.

- print_filtered_signatures:

  Logical: print filtered signatures. Default is \`FALSE\`.

- adjust_eset:

  Logical: clean problematic features. Default is \`FALSE\`.

- parallel.size:

  Parallel workers for ssGSEA. Default is 1.

- ...:

  Additional arguments passed to specific methods.

## Value

Tibble with phenotype data and signature scores.

## References

1.  Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation
    analysis. BMC Bioinformatics. 2013;14:7.

2.  Mariathasan S, et al. TGF\\\beta\\ attenuates tumour response to
    PD-L1 blockade. Nature. 2018;554:544-548.

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
# \donttest{
result <- calculate_sig_score(eset = eset, signature = signature, method = "pca")
#> ℹ Calculating signature scores using PCA method
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Calculating scores for 2 signature(s)
# }
```
