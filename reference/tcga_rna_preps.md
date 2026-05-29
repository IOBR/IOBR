# Preprocess TCGA RNA-seq Data

Preprocesses TCGA RNA-seq data by modifying sample types, transforming
data, and annotating genes based on specified parameters.

## Usage

``` r
tcga_rna_preps(
  eset,
  id_type = c("ensembl", "symbol"),
  input_type = c("log2count", "count"),
  output = c("tumor", "tumor_normal"),
  output_type = c("tpm", "log2tpm", "count"),
  annotation = TRUE
)
```

## Arguments

- eset:

  Matrix or data frame. RNA-seq gene expression matrix from TCGA.

- id_type:

  Character. Gene identifier type: "ensembl" or "symbol". Default is
  "ensembl".

- input_type:

  Character. Input data type: "log2count" or "count". Default is
  "log2count".

- output:

  Character. Sample type: "tumor" or "tumor_normal". Default is "tumor".

- output_type:

  Character. Output data type: "tpm", "log2tpm", or "count". Default is
  "tpm".

- annotation:

  Logical. Whether to perform gene annotation. Default is TRUE.

## Value

Preprocessed gene expression matrix.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate TCGA-style data with Ensembl IDs
set.seed(123)
sim_eset <- matrix(
  abs(rnorm(500)),
  nrow = 100,
  ncol = 5,
  dimnames = list(
    paste0("ENSG00000", 101:200),
    paste0("TCGA-AB-123", 4:8, "-01A")
  )
)
# Process tumor samples only (output as count for offline testing)
eset <- tcga_rna_preps(
  eset = sim_eset, id_type = "ensembl", input_type = "count",
  output = "tumor", output_type = "count", annotation = FALSE
)
#> ℹ Sample type distribution:
#> 
#> 01A 
#>   5 
#> ℹ Adjacent normal samples removed
#> ℹ TCGA barcode reduced to 12 digits
#> ℹ Duplicate barcode check:
#>    Mode   FALSE 
#> logical       5 
print(dim(eset))
#> [1] 100   5
```
