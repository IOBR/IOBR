# Calculate Significance Using Random Matrix

Calculates FDR-adjusted p-values using a random shuffled matrix.

## Usage

``` r
xCellSignifcanceRandomMatrix(
  scores,
  expr,
  spill,
  alpha = 0.5,
  nperm = 250,
  file.name = NULL
)
```

## Arguments

- scores:

  xCell scores matrix.

- expr:

  Input expression matrix.

- spill:

  Spillover object.

- alpha:

  Spillover alpha parameter. Default is \`0.5\`.

- nperm:

  Number of permutations. Default is \`250\`.

- file.name:

  Character string for saving p-values. Default is \`NULL\`.

## Value

List containing p-values, shuffled xCell scores, shuffled expression,
and beta distributions.
