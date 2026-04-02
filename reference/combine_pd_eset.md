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
set.seed(123)
eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(eset) <- paste0("Gene", 1:100)
colnames(eset) <- paste0("Sample", 1:10)
pdata <- data.frame(
  ID = colnames(eset),
  group = rep(c("A", "B"), each = 5),
  age = rnorm(10, 50, 10)
)
result <- combine_pd_eset(eset = eset, pdata = pdata, scale = FALSE)
#> ℹ Retained 100 features after QC
#> ✔ Combined data: 10 samples x 103 variables
dim(result)
#> [1]  10 103
```
