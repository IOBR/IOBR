# Deconvolve Using CIBERSORT

CIBERSORT is freely available to academic users. License and binary can
be obtained from https://cibersortx.stanford.edu.

## Usage

``` r
deconvo_cibersort(
  eset,
  project = NULL,
  arrays = FALSE,
  perm = 1000,
  absolute = FALSE,
  abs_method = "sig.score",
  parallel = FALSE,
  num_cores = 2,
  seed = NULL
)
```

## Arguments

- eset:

  Expression matrix with gene symbols as row names.

- project:

  Optional project name. Default is \`NULL\`.

- arrays:

  Logical: optimized for microarray data. Default is \`FALSE\`.

- perm:

  Permutations for statistical analysis. Default is 1000.

- absolute:

  Logical: run in absolute mode. Default is \`FALSE\`.

- abs_method:

  Method for absolute mode: \`"sig.score"\` or \`"no.sumto1"\`. Default
  is \`"sig.score"\`.

- parallel:

  Enable parallel execution. Default is \`FALSE\`.

- num_cores:

  Number of cores for parallel mode. Default is 2.

- seed:

  Random seed for reproducibility. Default is \`NULL\`.

## Value

Data frame with CIBERSORT cell fractions. Columns suffixed with
\`\_CIBERSORT\`.

## Author

Dongqiang Zeng

## Examples

``` r
eset_tme_stad <- load_data("eset_tme_stad")
lm22 <- load_data("lm22")
cibersort_result <- deconvo_cibersort(
  eset = eset_tme_stad,
  project = "TCGA-STAD",
  perm = 100
)
#> ℹ Running CIBERSORT
```
