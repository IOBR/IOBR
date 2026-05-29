# Visualize Expression Set Distribution

Generates boxplots and density plots to analyze the distribution of
expression values in an expression set. Useful for quality control and
assessing data normalization.

## Usage

``` r
eset_distribution(eset, quantile = 3, log = TRUE, project = NULL)
```

## Arguments

- eset:

  Expression matrix or data frame with genes in rows and samples in
  columns.

- quantile:

  Integer specifying the divisor for sampling columns. Default is 3
  (samples 1/3 of columns).

- log:

  Logical indicating whether to perform log2 transformation. Default is
  \`TRUE\`.

- project:

  Optional output directory path for saving files. If \`NULL\`, no files
  are saved. Default is \`NULL\`.

## Value

Invisibly returns \`NULL\`. If \`project\` is provided, saves PNG files
to disk.

## Examples

``` r
# Simulate data
set.seed(123)
sim_eset <- matrix(rnorm(1000 * 10, mean = 5, sd = 2), 1000, 10)
rownames(sim_eset) <- paste0("Gene", 1:1000)
colnames(sim_eset) <- paste0("Sample", 1:10)

# Run distribution plot
result <- eset_distribution(sim_eset)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
```
