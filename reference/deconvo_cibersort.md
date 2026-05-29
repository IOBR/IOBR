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
lm22 <- load_data("lm22")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "lm22"
if (!is.null(lm22)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(lm22) * 3), nrow(lm22), 3)
  rownames(sim_eset) <- rownames(lm22)
  colnames(sim_eset) <- paste0("Sample", 1:3)

  # Run deconvolution
  result <- deconvo_cibersort(
    eset = sim_eset,
    project = "TCGA-STAD",
    perm = 10
  )
  if (!is.null(result)) head(result)
}
#> ℹ Running CIBERSORT
#> ℹ Loading cached data: "lm22"
#>        ID ProjectID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT
#> 1 Sample1 TCGA-STAD              0.00000000               0.08619764
#> 2 Sample2 TCGA-STAD              0.06359520               0.00000000
#> 3 Sample3 TCGA-STAD              0.03145937               0.00000000
#>   Plasma_cells_CIBERSORT T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1             0.03406151            0.00000000                   0.0000000
#> 2             0.00000000            0.03341067                   0.2006443
#> 3             0.01943398            0.00000000                   0.0000000
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                            0.1686361                             0.00000000
#> 2                            0.0000000                             0.04465699
#> 3                            0.3518399                             0.00000000
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                         0.277934421                          0.000000000
#> 2                         0.000000000                          0.002049631
#> 3                         0.005228186                          0.084194087
#>   T_cells_gamma_delta_CIBERSORT NK_cells_resting_CIBERSORT
#> 1                    0.00000000                          0
#> 2                    0.08455261                          0
#> 3                    0.06415089                          0
#>   NK_cells_activated_CIBERSORT Monocytes_CIBERSORT Macrophages_M0_CIBERSORT
#> 1                   0.07127415          0.06590437               0.01297899
#> 2                   0.00000000          0.40590254               0.09312735
#> 3                   0.06329834          0.03628107               0.00000000
#>   Macrophages_M1_CIBERSORT Macrophages_M2_CIBERSORT
#> 1               0.05736787               0.00000000
#> 2               0.00000000               0.00000000
#> 3               0.06291034               0.08127811
#>   Dendritic_cells_resting_CIBERSORT Dendritic_cells_activated_CIBERSORT
#> 1                        0.00000000                                   0
#> 2                        0.01000038                                   0
#> 3                        0.00000000                                   0
#>   Mast_cells_resting_CIBERSORT Mast_cells_activated_CIBERSORT
#> 1                   0.19312749                              0
#> 2                   0.06206029                              0
#> 3                   0.14263721                              0
#>   Eosinophils_CIBERSORT Neutrophils_CIBERSORT P-value_CIBERSORT
#> 1            0.03251746            0.00000000               0.1
#> 2            0.00000000            0.00000000               0.3
#> 3            0.00000000            0.05728852               0.0
#>   Correlation_CIBERSORT RMSE_CIBERSORT
#> 1            0.10553554       1.053760
#> 2            0.08675359       1.059597
#> 3            0.12032892       1.044903
```
