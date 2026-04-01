# Identify Variable Genes in Expression Data

Identifies variable genes from a gene expression dataset using specified
selection criteria. Supports multiple methods, including expression
thresholding and variability estimation via median absolute deviation
(MAD).

Identifies variable genes from a gene expression dataset using specified
selection criteria. Supports multiple methods, including expression
thresholding and variability estimation via median absolute deviation
(MAD).

## Usage

``` r
find_variable_genes(
  eset,
  data_type = c("count", "normalized"),
  methods = c("low", "mad"),
  prop = 0.7,
  quantile = c(0.75, 0.5, 0.25),
  min.mad = 0.1,
  feas = NULL
)

find_variable_genes(
  eset,
  data_type = c("count", "normalized"),
  methods = c("low", "mad"),
  prop = 0.7,
  quantile = c(0.75, 0.5, 0.25),
  min.mad = 0.1,
  feas = NULL
)
```

## Arguments

- eset:

  Numeric matrix. Gene expression data (genes as rows, samples as
  columns).

- data_type:

  Character. Type of data: \`"count"\` or \`"normalized"\`. Default is
  \`"count"\`.

- methods:

  Character vector. Methods for gene selection: \`"low"\`, \`"mad"\`.
  Default is \`c("low", "mad")\`.

- prop:

  Numeric. Proportion of samples in which a gene must be expressed.
  Default is 0.7.

- quantile:

  Numeric. Quantile threshold for minimum MAD (0.25, 0.5, 0.75). Default
  is 0.75.

- min.mad:

  Numeric. Minimum allowable MAD value. Default is 0.1.

- feas:

  Character vector or \`NULL\`. Additional features to include. Default
  is \`NULL\`.

## Value

Matrix subset of \`eset\` containing variable genes.

Matrix subset of \`eset\` containing variable genes.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_tme_stad <- load_data("eset_tme_stad")
eset <- find_variable_genes(
  eset = eset_tme_stad,
  data_type = "normalized",
  methods = "mad",
  quantile = 0.25
)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ min.mad = 0.1
#> ℹ Range of MAD: 0.1 to 2.22
#> ℹ 25% of variables will be filtered out...
# }
# \donttest{
eset_tme_stad <- load_data("eset_tme_stad")
eset <- find_variable_genes(
  eset = eset_tme_stad,
  data_type = "normalized",
  methods = "mad",
  quantile = 0.25
)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ min.mad = 0.1
#> ℹ Range of MAD: 0.1 to 2.22
#> ℹ 25% of variables will be filtered out...
# }
```
