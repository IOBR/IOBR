# Test for Cell Population Infiltration

Returns p-values for the null hypothesis that samples are not
infiltrated by the corresponding cell population.

## Usage

``` r
test_for_infiltration(MCPcounterMatrix, platform = c("133P2", "133A", "HG1"))
```

## Arguments

- MCPcounterMatrix:

  Matrix, usually output from MCPcounter.estimate.

- platform:

  Expression platform: "133P2", "133A", or "HG1". Default is "133P2".

## Value

Matrix with samples in rows and cell populations in columns. Elements
are p-values.

## Author

Etienne Becht

## Examples

``` r
# This function requires null_models data which is loaded internally
# Create example data
scores <- matrix(runif(30), nrow = 3, ncol = 10)
rownames(scores) <- c("T cells", "B cells", "NK cells")
pvals <- test_for_infiltration(scores, platform = "133P2")
```
