# Combine Phenotype Data and Expression Set

Merges phenotype data with an expression matrix by matching sample IDs.
Optionally filters features, applies feature manipulation, and scales
expression data before combining.

## Usage

``` r
combine_pd_eset(
  eset,
  pdata,
  id_pdata = "ID",
  feas = NULL,
  feature_manipulation = TRUE,
  scale = TRUE,
  choose_who_when_duplicate = c("eset", "pdata")
)
```

## Arguments

- eset:

  Expression matrix with genes/features in rows and samples in columns.

- pdata:

  Data frame containing phenotype/clinical data.

- id_pdata:

  Character string specifying the column name in \`pdata\` containing
  sample identifiers. Default is \`"ID"\`.

- feas:

  Character vector specifying features to include from \`eset\`. If
  \`NULL\`, all features are used. Default is \`NULL\`.

- feature_manipulation:

  Logical indicating whether to apply feature manipulation to filter
  valid features. Default is \`TRUE\`.

- scale:

  Logical indicating whether to scale (standardize) expression data.
  Default is \`TRUE\`.

- choose_who_when_duplicate:

  Character string specifying which data to prefer when duplicate
  columns exist. Options are \`"eset"\` or \`"pdata"\`. Default is
  \`"eset"\`.

## Value

Data frame combining phenotype data and (transposed) expression data,
with samples in rows and features/phenotypes in columns.

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
colnames(eset) <- substring(colnames(eset), 1, 12)

sig_stad <- load_data("sig_stad")

# Example 1: Combine all features (no filtering)
input <- combine_pd_eset(
  eset = eset,
  pdata = sig_stad
)
#> ℹ Retained 48730 features after QC
#> ✔ Combined data: 10 samples x 49053 variables

# Example 2: Combine with specific features
# Note: features must match rownames of eset after ID conversion
input2 <- combine_pd_eset(
  eset = eset,
  pdata = sig_stad,
  feas = rownames(eset)[1:100] # Use first 100 genes as example
)
#> ℹ Filtered to 100 features
#> ℹ Retained 100 features after QC
#> ✔ Combined data: 10 samples x 423 variables
# }
```
