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
if (FALSE) { # \dontrun{
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
xcell_result <- deconvo_xcell(eset = eset, project = "TCGA-STAD")
head(xcell_result)
} # }
```
