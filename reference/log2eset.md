# Log2 Transformation of Gene Expression Matrix

Determines whether a gene expression matrix requires log2 transformation
based on the distribution of values, and applies it if necessary. This
is useful for automatically detecting raw counts or linear-scale data
that should be log-transformed for downstream analysis.

## Usage

``` r
log2eset(eset)
```

## Arguments

- eset:

  Numeric matrix. Gene expression data with genes as rows and samples in
  columns.

## Value

Numeric matrix. Log2-transformed gene expression data (if transformation
was needed), or the original data otherwise.

## Examples

``` r
# \donttest{
# Load TCGA-STAD expression data (raw count matrix)
eset_stad <- load_data("eset_stad")

# Transform count data to TPM
eset <- count2tpm(eset_stad, idType = "ensembl")
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> Warning: longer object length is not a multiple of shorter object length
#> ℹ No duplicate gene symbols found.

# Apply log2 transformation if needed
eset <- log2eset(eset)
#> ✔ Applied log2 transformation
# }
```
