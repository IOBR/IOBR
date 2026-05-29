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
# Simulate data
set.seed(123)
sim_sig <- matrix(rnorm(100 * 5), 100, 5)
rownames(sim_sig) <- paste0("Gene", 1:100)
colnames(sim_sig) <- paste0("Cell", 1:5)

sim_mixture <- as.data.frame(matrix(
  rnorm(100 * 3), 100, 3
))
rownames(sim_mixture) <- paste0("Gene", 1:100)
colnames(sim_mixture) <- paste0("Sample", 1:3)

# Run deconvolution
result <- CIBERSORT(
  sig_matrix = sim_sig,
  mixture_file = sim_mixture,
  perm = 10, QN = FALSE, absolute = FALSE,
  parallel = FALSE
)
if (!is.null(result)) head(result)
#>             Cell1     Cell2      Cell3     Cell4     Cell5 P-value Correlation
#> Sample1 0.2661124 0.3560011 0.09902726 0.0000000 0.2788592     0.2  0.14713328
#> Sample2 0.1987843 0.3010143 0.50020143 0.0000000 0.0000000     0.8  0.04060941
#> Sample3 0.4232823 0.1873284 0.00000000 0.3893893 0.0000000     0.1  0.17900877
#>             RMSE
#> Sample1 1.031594
#> Sample2 1.132430
#> Sample3 1.062058
```
