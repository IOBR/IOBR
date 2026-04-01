# The xCell Analysis Pipeline

Returns the xCell cell types enrichment scores for tumor
microenvironment deconvolution. Uses ssGSEA-based enrichment analysis
with spillover compensation to estimate cell type proportions from gene
expression data.

## Usage

``` r
xCellAnalysis(
  expr,
  signatures = NULL,
  genes = NULL,
  spill = NULL,
  rnaseq = TRUE,
  file.name = NULL,
  scale = TRUE,
  alpha = 0.5,
  save.raw = FALSE,
  cell.types.use = NULL
)
```

## Arguments

- expr:

  Gene expression data matrix with row names as gene symbols and columns
  as samples.

- signatures:

  GMT object of signatures. If \`NULL\`, uses xCell defaults.

- genes:

  Character vector of genes to use in the analysis. If \`NULL\`, uses
  xCell defaults.

- spill:

  Spillover object for adjusting scores. If \`NULL\`, uses xCell
  defaults.

- rnaseq:

  Logical indicating whether to use RNA-seq (TRUE) or array (FALSE)
  spillover parameters. Default is \`TRUE\`.

- file.name:

  Character string for saving scores. Default is \`NULL\`.

- scale:

  Logical indicating whether to use scaling. Default is \`TRUE\`.

- alpha:

  Numeric value to override spillover alpha parameter. Default is
  \`0.5\`.

- save.raw:

  Logical indicating whether to save raw scores. Default is \`FALSE\`.

- cell.types.use:

  Character vector of cell types to use. If \`NULL\`, uses all available
  cell types. Default is \`NULL\`.

## Value

A matrix of adjusted xCell scores.
