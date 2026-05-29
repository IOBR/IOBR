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
xcell <- load_data("xCell.data")
if (!is.null(xcell)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(length(xcell$genes) * 2), length(xcell$genes), 2)
  rownames(sim_eset) <- xcell$genes
  colnames(sim_eset) <- paste0("Sample", 1:2)
  result <- deconvo_xcell(eset = sim_eset, project = "TCGA-STAD")
  if (!is.null(result)) head(result)
}
} # }
```
