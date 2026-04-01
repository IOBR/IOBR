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
  (https://cibersort.stanford.edu/runcibersort.php)

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
# Create simulated data matching LM22 signature matrix gene names
data(lm22)
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
#> Sample1             0     0.18225155  0.000000000  0.00000000         0.0000000
#> Sample2             0     0.00000000  0.009197916  0.09708868         0.1543896
#> Sample3             0     0.00000000  0.031932727  0.00000000         0.1179203
#> Sample4             0     0.02848117  0.000000000  0.00000000         0.0000000
#> Sample5             0     0.13461164  0.031886333  0.00000000         0.1400873
#> Sample6             0     0.04968648  0.000000000  0.30487427         0.0000000
#>         T cells CD4 memory resting T cells CD4 memory activated
#> Sample1                0.000000000                            0
#> Sample2                0.000000000                            0
#> Sample3                0.000000000                            0
#> Sample4                0.141609900                            0
#> Sample5                0.000000000                            0
#> Sample6                0.006337446                            0
#>         T cells follicular helper T cells regulatory (Tregs)
#> Sample1                0.06818487                 0.06879419
#> Sample2                0.18578157                 0.00000000
#> Sample3                0.08071338                 0.11463869
#> Sample4                0.03146895                 0.00000000
#> Sample5                0.16872688                 0.00000000
#> Sample6                0.00000000                 0.04591456
#>         T cells gamma delta NK cells resting NK cells activated  Monocytes
#> Sample1         0.007867152       0.03349354        0.000000000 0.00000000
#> Sample2         0.032333507       0.04837349        0.000000000 0.04976289
#> Sample3         0.017656394       0.00000000        0.112270547 0.25979544
#> Sample4         0.000000000       0.02212775        0.067613042 0.00000000
#> Sample5         0.157583874       0.00000000        0.005653372 0.00000000
#> Sample6         0.000000000       0.08006571        0.009674645 0.13915701
#>         Macrophages M0 Macrophages M1 Macrophages M2 Dendritic cells resting
#> Sample1      0.1305253   0.000000e+00    0.220600655              0.00000000
#> Sample2      0.0000000   0.000000e+00    0.221272980              0.00000000
#> Sample3      0.0000000   0.000000e+00    0.023165264              0.00000000
#> Sample4      0.1477710   4.818640e-02    0.008736663              0.09074849
#> Sample5      0.0000000   1.954162e-05    0.046899915              0.02008421
#> Sample6      0.0000000   2.501610e-02    0.000000000              0.09413481
#>         Dendritic cells activated Mast cells resting Mast cells activated
#> Sample1                0.00000000         0.00000000           0.12404322
#> Sample2                0.11746095         0.05457561           0.00000000
#> Sample3                0.03702223         0.17010062           0.00000000
#> Sample4                0.05891092         0.35434569           0.00000000
#> Sample5                0.00000000         0.14752717           0.03462643
#> Sample6                0.00000000         0.00000000           0.24513896
#>         Eosinophils Neutrophils P-value  Correlation     RMSE
#> Sample1  0.16423953  0.00000000     0.9 -0.028504120 1.101778
#> Sample2  0.02976284  0.00000000     0.9 -0.024426465 1.088628
#> Sample3  0.03478440  0.00000000     0.8 -0.003341701 1.081645
#> Sample4  0.00000000  0.00000000     0.9 -0.049245746 1.103864
#> Sample5  0.08494039  0.02735297     0.3  0.051526620 1.062632
#> Sample6  0.00000000  0.00000000     0.8 -0.007026086 1.109070
# }
```
