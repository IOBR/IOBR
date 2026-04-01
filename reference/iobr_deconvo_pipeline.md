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
  path = "1-TME",
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

  Character. Output directory. Default is "1-TME".

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
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
eset <- eset[1:500, 1:5]
res <- iobr_deconvo_pipeline(
  eset = eset, project = "STAD",
  array = FALSE, tumor_type = "stad",
  path = tempdir(), permutation = 10
)
} # }
```
