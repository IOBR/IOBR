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

  Optional output directory path for saving files. If \`NULL\`, no files
  are saved. Default is \`NULL\`.

## Value

Invisibly returns \`NULL\`. If \`project\` is provided, saves PNG files
to disk.

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "eset_stad"
eset_distribution(eset_stad[1:1000, ])
#> ✔ Applied log2 transformation
#anno_rnaseq <- load_data("anno_rnaseq")
#eset <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#eset_distribution(eset)
#eset_distribution(eset, project = file.path(tempdir(), "ESET"))
# }
```
