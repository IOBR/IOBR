# CIBERSORT Deconvolution Algorithm

An analytical tool to estimate cell type abundances in mixed cell
populations using gene expression data.

## Usage

``` r
CIBERSORT(
  sig_matrix = NULL,
  mixture_file,
  perm,
  QN = TRUE,
  absolute = FALSE,
  abs_method = "sig.score",
  parallel = FALSE,
  num_cores = 2,
  seed = NULL
)
```

## Arguments

- sig_matrix:

  Cell type GEP barcode matrix: row 1 = sample labels; column 1 = gene
  symbols; no missing values; default = LM22.txt download from CIBERSORT
  (https://cibersortx.stanford.edu/runcibersort.php)

- mixture_file:

  GEP matrix: row 1 = sample labels; column 1 = gene symbols; no missing
  values

- perm:

  Set permutations for statistical analysis (\>=100 permutations
  recommended).

- QN:

  Quantile normalization of input mixture (default = TRUE)

- absolute:

  Run CIBERSORT in absolute mode (default = FALSE) - note that cell
  subsets will be scaled by their absolute levels and will not be
  represented as fractions (to derive the default output, normalize
  absolute levels such that they sum to 1 for each mixture sample) - the
  sum of all cell subsets in each mixture sample will be added to the
  output ('Absolute score'). If LM22 is used, this score will capture
  total immune content.

- abs_method:

  If absolute is set to TRUE, choose method: 'no.sumto1' or
  'sig.score' - sig.score = for each mixture sample, define S as the
  median expression level of all genes in the signature matrix divided
  by the median expression level of all genes in the mixture. Multiple
  cell subset fractions by S. - no.sumto1 = remove sum to 1 constraint

- parallel:

  Logical. Enable parallel execution? (default = FALSE)

- num_cores:

  Integer. Number of cores to use when `parallel = TRUE` (default = 2)

- seed:

  Integer. Random seed for reproducible permutation testing. If `NULL`
  (default), uses current random state. Set to a specific value
  (e.g., 123) for reproducible results across runs. Applies to both
  parallel and serial permutation.

## Value

A matrix object containing the estimated cibersort-cell fractions,
p-values, correlation coefficients, and RMSE values.

## Author

Aaron M. Newman, Stanford University (amnewman@stanford.edu)

## Examples

``` r
# \donttest{
lm22 <- load_data("lm22")
#> ℹ Trying mirror 1/11: <https://github.com>
#> ✔ Download complete: "lm22"
common_genes <- rownames(lm22)[1:500]
sim_mixture <- as.data.frame(matrix(
  rnorm(length(common_genes) * 10, mean = 5, sd = 2),
  nrow = length(common_genes), ncol = 10
))
rownames(sim_mixture) <- common_genes
colnames(sim_mixture) <- paste0("Sample", 1:10)
result <- CIBERSORT(
  sig_matrix = lm22,
  mixture_file = sim_mixture,
  perm = 10, QN = FALSE, absolute = FALSE,
  parallel = FALSE
)
head(result)
#>         B cells naive B cells memory Plasma cells T cells CD8 T cells CD4 naive
#> Sample1    0.02412905     0.00000000  0.037892551 0.019998246         0.0000000
#> Sample2    0.00000000     0.07172372  0.002375459 0.044081770         0.1042969
#> Sample3    0.11411301     0.00000000  0.000000000 0.000000000         0.1783976
#> Sample4    0.11715171     0.00000000  0.037357552 0.000000000         0.0265595
#> Sample5    0.26525360     0.00000000  0.000000000 0.004074691         0.0000000
#> Sample6    0.00000000     0.20875624  0.000000000 0.000000000         0.2307710
#>         T cells CD4 memory resting T cells CD4 memory activated
#> Sample1                 0.40885968                    0.0000000
#> Sample2                 0.31821463                    0.0000000
#> Sample3                 0.00000000                    0.0000000
#> Sample4                 0.03621121                    0.2325042
#> Sample5                 0.37245162                    0.0000000
#> Sample6                 0.00000000                    0.2167616
#>         T cells follicular helper T cells regulatory (Tregs)
#> Sample1               0.000000000                0.013494054
#> Sample2               0.000000000                0.000000000
#> Sample3               0.009737929                0.003535199
#> Sample4               0.000000000                0.000000000
#> Sample5               0.000000000                0.086789363
#> Sample6               0.000000000                0.000000000
#>         T cells gamma delta NK cells resting NK cells activated  Monocytes
#> Sample1           0.0000000       0.09884303        0.000000000 0.19894746
#> Sample2           0.0000000       0.00000000        0.182426896 0.00000000
#> Sample3           0.1649016       0.05216157        0.000000000 0.04518743
#> Sample4           0.1515656       0.00000000        0.181914131 0.00000000
#> Sample5           0.0000000       0.02657193        0.009252027 0.04745452
#> Sample6           0.0000000       0.00000000        0.121123571 0.04279619
#>         Macrophages M0 Macrophages M1 Macrophages M2 Dendritic cells resting
#> Sample1    0.000000000    0.000000000              0             0.039739550
#> Sample2    0.009086296    0.000000000              0             0.000000000
#> Sample3    0.055247741    0.013230023              0             0.214292844
#> Sample4    0.038951497    0.000000000              0             0.001887852
#> Sample5    0.022309752    0.006959781              0             0.081581053
#> Sample6    0.070107782    0.000000000              0             0.010842482
#>         Dendritic cells activated Mast cells resting Mast cells activated
#> Sample1              0.0512166001         0.07712769           0.02975208
#> Sample2              0.0309039896         0.17390011           0.01904856
#> Sample3              0.0000000000         0.07987597           0.00000000
#> Sample4              0.0175130310         0.09108622           0.00000000
#> Sample5              0.0000000000         0.06087293           0.00000000
#> Sample6              0.0008915717         0.09794957           0.00000000
#>         Eosinophils Neutrophils P-value   Correlation     RMSE
#> Sample1  0.00000000  0.00000000     0.0  0.0812753089 1.051259
#> Sample2  0.00000000  0.04394169     0.5 -0.0017905115 1.122588
#> Sample3  0.06931906  0.00000000     0.0  0.1174831102 1.031201
#> Sample4  0.06729750  0.00000000     0.5  0.0003266841 1.129563
#> Sample5  0.01642874  0.00000000     0.2  0.0340252078 1.091239
#> Sample6  0.00000000  0.00000000     0.4  0.0150830890 1.087891
# }
```
