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
# Simulate data
set.seed(123)
sim_ref <- matrix(rnorm(100 * 5), 100, 5)
rownames(sim_ref) <- paste0("Gene", 1:100)
colnames(sim_ref) <- paste0("CellType", 1:5)

sim_eset <- matrix(rnorm(100 * 3), 100, 3)
rownames(sim_eset) <- paste0("Gene", 1:100)
colnames(sim_eset) <- paste0("Sample", 1:3)

# Run deconvolution
result <- deconvo_ref(eset = sim_eset, reference = sim_ref, method = "lsei")
#> ℹ Found 100 common genes
#> ℹ Running lsei deconvolution
if (!is.null(result)) head(result)
#>   ID CellType1_CIBERSORT CellType2_CIBERSORT CellType3_CIBERSORT
#> 1  1                 0.2                 0.2                 0.2
#> 2  2                 0.2                 0.2                 0.2
#> 3  3                 0.2                 0.2                 0.2
#>   CellType4_CIBERSORT CellType5_CIBERSORT
#> 1                 0.2                 0.2
#> 2                 0.2                 0.2
#> 3                 0.2                 0.2
```
