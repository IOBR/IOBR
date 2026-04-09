# Deconvolve Immune Microenvironment Using xCell

Estimates immune cell fractions using the xCell algorithm. xCell
provides cell type enrichment scores for 64 immune and stromal cell
types from gene expression data.

## Usage

``` r
deconvo_xcell(eset, project = NULL, arrays = FALSE)
```

## Arguments

- eset:

  Gene expression matrix with HGNC gene symbols as row names and samples
  as columns.

- project:

  Optional project name to add as \`ProjectID\` column. Default is
  \`NULL\`.

- arrays:

  Logical indicating microarray data (\`TRUE\`) or RNA-seq (\`FALSE\`).
  Default is \`FALSE\`.

## Value

Data frame with xCell enrichment scores. Cell type columns are suffixed
with \`\_xCell\`.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_grch38 <- load_data("anno_grch38")
#> ℹ Loading cached data: "anno_grch38"
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
xcell_result <- deconvo_xcell(eset = eset[, 1:3], project = "TCGA-STAD")
#> ℹ Running xCell deconvolution
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "xCell.data"
#> ℹ Number of genes: 10783
#> ℹ GSVA version 2.4.9
#> ℹ Searching for rows with constant values
#> ℹ Calculating GSVA ranks
#> ℹ kcdf='auto' (default)
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores
#> ✔ Calculations finished
head(xcell_result)[, 1:5]
#>             ID ProjectID Adipocytes_xCell Astrocytes_xCell B-cells_xCell
#> 1 TCGA-BR-6455 TCGA-STAD     3.755491e-23     0.000000e+00  0.000000e+00
#> 2 TCGA-BR-7196 TCGA-STAD     1.794635e-23     3.702040e-22  1.311341e-23
#> 3 TCGA-BR-8371 TCGA-STAD     0.000000e+00     1.881182e-21  3.149663e-08
# }
```
