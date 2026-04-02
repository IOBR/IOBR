# Deconvolve Using Custom Reference

Cell fraction estimation using SVR or lsei methods with custom
reference.

## Usage

``` r
deconvo_ref(
  eset,
  project = NULL,
  arrays = TRUE,
  method = c("svr", "lsei"),
  perm = 100,
  reference,
  scale_reference = TRUE,
  absolute.mode = FALSE,
  abs.method = "sig.score"
)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- arrays:

  Logical: use quantile normalization. Default is \`TRUE\`.

- method:

  Method: \`"svr"\` or \`"lsei"\`. Default is \`"svr"\`.

- perm:

  Permutations for SVR. Default is 100.

- reference:

  Custom reference matrix (e.g., lm22, lm6).

- scale_reference:

  Logical: scale reference. Default is \`TRUE\`.

- absolute.mode:

  Logical: absolute mode for SVR. Default is \`FALSE\`.

- abs.method:

  Method for absolute mode. Default is \`"sig.score"\`.

## Value

Data frame with cell fractions. Columns suffixed with \`\_CIBERSORT\`.

## Author

Dongqiang Zeng, Rongfang Shen

## Examples

``` r
lm22 <- load_data("lm22")
common_genes <- rownames(lm22)[1:500]
sim_eset <- as.data.frame(matrix(
  rnorm(length(common_genes) * 5, mean = 5, sd = 2),
  nrow = length(common_genes), ncol = 5
))
rownames(sim_eset) <- common_genes
colnames(sim_eset) <- paste0("Sample", 1:5)
deconvo_ref(eset = sim_eset, reference = lm22, method = "lsei")
#> ℹ Found 500 common genes
#> ℹ Running lsei deconvolution
#>   ID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT Plasma_cells_CIBERSORT
#> 1  1              0.04545455               0.04545455             0.04545455
#> 2  2              0.04545455               0.04545455             0.04545455
#> 3  3              0.04545455               0.04545455             0.04545455
#> 4  4              0.04545455               0.04545455             0.04545455
#> 5  5              0.04545455               0.04545455             0.04545455
#>   T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1            0.04545455                  0.04545455
#> 2            0.04545455                  0.04545455
#> 3            0.04545455                  0.04545455
#> 4            0.04545455                  0.04545455
#> 5            0.04545455                  0.04545455
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                           0.04545455                             0.04545455
#> 2                           0.04545455                             0.04545455
#> 3                           0.04545455                             0.04545455
#> 4                           0.04545455                             0.04545455
#> 5                           0.04545455                             0.04545455
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                          0.04545455                           0.04545455
#> 2                          0.04545455                           0.04545455
#> 3                          0.04545455                           0.04545455
#> 4                          0.04545455                           0.04545455
#> 5                          0.04545455                           0.04545455
#>   T_cells_gamma_delta_CIBERSORT NK_cells_resting_CIBERSORT
#> 1                    0.04545455                 0.04545455
#> 2                    0.04545455                 0.04545455
#> 3                    0.04545455                 0.04545455
#> 4                    0.04545455                 0.04545455
#> 5                    0.04545455                 0.04545455
#>   NK_cells_activated_CIBERSORT Monocytes_CIBERSORT Macrophages_M0_CIBERSORT
#> 1                   0.04545455          0.04545455               0.04545455
#> 2                   0.04545455          0.04545455               0.04545455
#> 3                   0.04545455          0.04545455               0.04545455
#> 4                   0.04545455          0.04545455               0.04545455
#> 5                   0.04545455          0.04545455               0.04545455
#>   Macrophages_M1_CIBERSORT Macrophages_M2_CIBERSORT
#> 1               0.04545455               0.04545455
#> 2               0.04545455               0.04545455
#> 3               0.04545455               0.04545455
#> 4               0.04545455               0.04545455
#> 5               0.04545455               0.04545455
#>   Dendritic_cells_resting_CIBERSORT Dendritic_cells_activated_CIBERSORT
#> 1                        0.04545455                          0.04545455
#> 2                        0.04545455                          0.04545455
#> 3                        0.04545455                          0.04545455
#> 4                        0.04545455                          0.04545455
#> 5                        0.04545455                          0.04545455
#>   Mast_cells_resting_CIBERSORT Mast_cells_activated_CIBERSORT
#> 1                   0.04545455                     0.04545455
#> 2                   0.04545455                     0.04545455
#> 3                   0.04545455                     0.04545455
#> 4                   0.04545455                     0.04545455
#> 5                   0.04545455                     0.04545455
#>   Eosinophils_CIBERSORT Neutrophils_CIBERSORT
#> 1            0.04545455            0.04545455
#> 2            0.04545455            0.04545455
#> 3            0.04545455            0.04545455
#> 4            0.04545455            0.04545455
#> 5            0.04545455            0.04545455
```
