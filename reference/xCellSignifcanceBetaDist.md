# Calculate Significance P-values Using Beta Distribution

Calculates FDR-adjusted p-values for the null hypothesis that a cell
type is not present in the mixture.

## Usage

``` r
xCellSignifcanceBetaDist(
  scores,
  beta_params = NULL,
  rnaseq = TRUE,
  file.name = NULL
)
```

## Arguments

- scores:

  xCell scores matrix.

- beta_params:

  Pre-calculated beta distribution parameters.

- rnaseq:

  Logical for RNA-seq vs array parameters. Default is \`TRUE\`.

- file.name:

  Character string for saving p-values. Default is \`NULL\`.

## Value

Matrix of p-values.
