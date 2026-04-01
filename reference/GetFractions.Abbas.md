# Constrained Regression Method (Abbas et al., 2009)

Implements a constrained regression approach described by Abbas et al.
(2009). Estimates proportions of immune cell types within mixed cancer
tissue samples based on gene expression data. Iteratively adjusts
regression coefficients to ensure non-negative values.

## Usage

``` r
GetFractions.Abbas(XX, YY, w = NA)
```

## Arguments

- XX:

  Matrix. Immune expression data with genes as rows and cell types as
  columns.

- YY:

  Vector. Cancer expression data with gene expression levels.

- w:

  Vector or NA. Weights for regression. Default is NA (no weights).

## Value

Vector with non-negative coefficients representing proportions of each
cell type.

## Examples

``` r
XX <- matrix(runif(100), nrow = 10, ncol = 10)
colnames(XX) <- paste("CellType", 1:10, sep = "")
YY <- runif(10)
results <- GetFractions.Abbas(XX, YY)
print(results)
#>  CellType1  CellType2  CellType3  CellType4  CellType5  CellType6  CellType7 
#>  0.1303706  0.1698972  0.0000000  0.0000000  0.0000000  0.4695632  0.0000000 
#>  CellType8  CellType9 CellType10 
#>  0.0000000  0.3164912  0.0000000 
```
