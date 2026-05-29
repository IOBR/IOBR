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
if (FALSE) { # \dontrun{
lm22 <- load_data("lm22")
if (!is.null(lm22)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(lm22) * 2), nrow(lm22), 2)
  rownames(sim_eset) <- rownames(lm22)
  colnames(sim_eset) <- paste0("Sample", 1:2)
  result <- deconvo_cibersort(eset = sim_eset, project = "TCGA-STAD", perm = 10)
  if (!is.null(result)) head(result)
}
} # }
```
