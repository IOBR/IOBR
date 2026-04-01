# Principal Component Analysis (PCA) Visualization

This function performs Principal Component Analysis (PCA) on gene
expression data, reduces dimensionality while preserving variance, and
generates a scatter plot visualization.

## Usage

``` r
iobr_pca(
  data,
  is.matrix = TRUE,
  scale = TRUE,
  is.log = FALSE,
  pdata,
  id_pdata = "ID",
  group = NULL,
  geom.ind = "point",
  cols = "normal",
  palette = "jama",
  repel = FALSE,
  ncp = 5,
  axes = c(1, 2),
  addEllipses = TRUE
)
```

## Arguments

- data:

  Input data for PCA: matrix or data frame.

- is.matrix:

  Logical indicating if input is a matrix. Default is TRUE.

- scale:

  Logical indicating whether to scale the data. Default is TRUE.

- is.log:

  Logical indicating whether to log-transform the data. Default is
  FALSE.

- pdata:

  Data frame with sample IDs and grouping information.

- id_pdata:

  Column name in \`pdata\` for sample IDs. Default is "ID".

- group:

  Column name in \`pdata\` for grouping variable. Default is NULL.

- geom.ind:

  Type of geometric representation for points. Default is "point".

- cols:

  Color scheme for groups. Default is "normal".

- palette:

  Color palette for groups. Default is "jama".

- repel:

  Logical indicating whether to repel overlapping points. Default is
  FALSE.

- ncp:

  Number of principal components to retain. Default is 5.

- axes:

  Principal components to plot (e.g., c(1, 2)). Default is c(1, 2).

- addEllipses:

  Logical indicating whether to add concentration ellipses. Default is
  TRUE.

## Value

A ggplot object of the PCA plot.

## Author

Dongqiang Zeng

## Examples

``` r
eset_stad <- load_data("eset_stad")
eset <- count2tpm(eset_stad)
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> Warning: longer object length is not a multiple of shorter object length
#> ℹ No duplicate gene symbols found.
stad_group <- load_data("stad_group")
iobr_pca(eset,
  is.matrix = TRUE, scale = TRUE,
  is.log = TRUE, pdata = stad_group,
  id_pdata = "ID", group = "subtype"
)
#> ✔ Applied log2 transformation
#> 
#> EBV  GS 
#>   5   5 
#> 55
#> [1] ">>== colors for group: "
#> >>== #374E55FF>>== #DF8F44FF
#> Ignoring unknown labels:
#> • linetype : "subtype"

```
