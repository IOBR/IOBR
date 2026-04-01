# Visualize Expression Set Distribution

Generates boxplots and density plots to analyze the distribution of
expression values in an expression set. Useful for quality control and
assessing data normalization.

## Usage

``` r
eset_distribution(eset, quantile = 3, log = TRUE, project = NULL)
```

## Arguments

- eset:

  Expression matrix or data frame with genes in rows and samples in
  columns.

- quantile:

  Integer specifying the divisor for sampling columns. Default is 3
  (samples 1/3 of columns).

- log:

  Logical indicating whether to perform log2 transformation. Default is
  \`TRUE\`.

- project:

  Optional project name for output files. Default is \`NULL\` (uses
  "ESET").

## Value

Invisibly returns \`NULL\`. Side effect: saves PNG files to disk.

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
anno_rnaseq <- load_data("anno_rnaseq")
eset <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2098 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50139
eset_distribution(eset, project = file.path(tempdir(), "ESET"))
#> ✔ Applied log2 transformation
# }
```
