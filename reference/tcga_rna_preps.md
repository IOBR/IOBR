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
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
eset <- tcga_rna_preps(
  eset = eset_stad, id_type = "ensembl", input_type = "count",
  output = "tumor", output_type = "tpm", annotation = TRUE
)
#> ℹ Sample type distribution:
#> 
#>    
#> 10 
#> ℹ Adjacent normal samples removed
#> ℹ TCGA barcode reduced to 12 digits
#> ℹ Duplicate barcode check:
#>    Mode   FALSE 
#> logical      10 
#> ℹ Converting count to TPM
#> ℹ Loading cached data: "anno_grch38"
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> ℹ Found 1679 duplicate symbols. Using "mean" for ranking.
#> ✔ Reduced to 54658 unique genes
# }
```
