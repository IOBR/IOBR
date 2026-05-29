# Tumor Microenvironment (TME) Deconvolution Pipeline

Executes an integrated TME analysis on a gene expression matrix:
performs immune/stromal cell deconvolution using multiple algorithms,
computes signature scores, and aggregates results. Designed for
exploratory immunogenomic profiling.

## Usage

``` r
iobr_deconvo_pipeline(
  eset,
  project,
  array,
  tumor_type,
  path = NULL,
  permutation = 1000
)
```

## Arguments

- eset:

  Numeric matrix. Gene expression (TPM/log scale) with genes in rows.

- project:

  Character. Project name (used in output naming).

- array:

  Logical. Whether data originated from an array platform. Affects
  deconvolution choices.

- tumor_type:

  Character. Tumor type code (e.g., "stad") used by certain methods.

- path:

  Character. Output directory. Default is NULL (uses tempdir()).

- permutation:

  Integer. Number of permutations for CIBERSORT (and similar). Default
  is 1000.

## Value

Data frame integrating cell fractions and signature scores (also writes
intermediate outputs to disk).

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
lm22 <- load_data("lm22")
cancer_genes <- load_data("cancer_type_genes")
if (!is.null(lm22) && !is.null(cancer_genes)) {
  set.seed(123)
  genes <- rownames(lm22)
  xcell <- load_data("xCell.data")
  if (!is.null(xcell)) genes <- unique(c(genes, xcell$genes))
  genes <- unique(c(genes, cancer_genes[["stad"]]))
  eset <- matrix(runif(length(genes) * 2), nrow = length(genes), ncol = 2)
  rownames(eset) <- genes
  colnames(eset) <- paste0("Sample", 1:2)
  res <- iobr_deconvo_pipeline(
    eset = eset, project = "TEST",
    array = FALSE, tumor_type = "stad",
    path = tempdir(), permutation = 2
  )
  if (!is.null(res)) head(res)
}
} # }
```
