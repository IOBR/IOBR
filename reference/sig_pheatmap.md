# Generate Heatmap for Signature Data

Creates a heatmap from signature data with grouping variables, offering
flexible options for colors, clustering, and output formats using
ComplexHeatmap.

## Usage

``` r
sig_pheatmap(
  input,
  feas,
  group,
  group2 = NULL,
  group3 = NULL,
  ID = "ID",
  path = NULL,
  cols1 = "random",
  cols2 = "random",
  cols3 = "random",
  seed = 54321,
  show_col = FALSE,
  palette1 = 1,
  palette2 = 2,
  palette3 = 3,
  cluster_cols = TRUE,
  palette_for_heatmape = 6,
  scale.matrix = TRUE,
  cellwidth = 1,
  cellheight = 9,
  show_colnames = FALSE,
  fig.type = "pdf",
  width = 6,
  height = NULL,
  file_name_prefix = 1
)
```

## Arguments

- input:

  Data frame with variables in columns.

- feas:

  Character vector. Feature names (columns) to include in heatmap.

- group:

  Character string. Column name for primary grouping variable.

- group2:

  Character string or \`NULL\`. Optional secondary grouping variable.

- group3:

  Character string or \`NULL\`. Optional tertiary grouping variable.

- ID:

  Character string. Column name for sample identifiers. Default is
  \`"ID"\`.

- path:

  Character string or \`NULL\`. Directory to save output files. Default
  creates \`"Marker-heatmap-average"\`.

- cols1:

  Character vector or \`"random"\` or \`"normal"\`. Colors for primary
  group. Default is \`"random"\`.

- cols2:

  Character vector or \`"random"\` or \`"normal"\`. Colors for secondary
  group. Default is \`"random"\`.

- cols3:

  Character vector or \`"random"\` or \`"normal"\`. Colors for tertiary
  group. Default is \`"random"\`.

- seed:

  Integer. Random seed for color generation. Default is \`54321\`.

- show_col:

  Logical indicating whether to display colors. Default is \`FALSE\`.

- palette1:

  Integer. Palette for primary group. Default is \`1\`.

- palette2:

  Integer. Palette for secondary group. Default is \`2\`.

- palette3:

  Integer. Palette for tertiary group. Default is \`3\`.

- cluster_cols:

  Logical indicating whether to cluster columns. Default is \`TRUE\`.

- palette_for_heatmape:

  Integer. Palette number for heatmap. Default is \`6\`.

- scale.matrix:

  Logical indicating whether to scale the matrix. Default is \`TRUE\`.

- cellwidth:

  Numeric. Width of each cell in points. Default is \`1\`.

- cellheight:

  Numeric. Height of each cell in points. Default is \`9\`.

- show_colnames:

  Logical indicating whether to show column names. Default is \`FALSE\`.

- fig.type:

  Character string. File format for saving. Default is \`"pdf"\`.

- width:

  Numeric. Width of saved figure in inches. Default is \`6\`.

- height:

  Numeric or \`NULL\`. Height of saved figure in inches. Calculated if
  \`NULL\`.

- file_name_prefix:

  Character or numeric. Prefix for saved file name. Default is \`1\`.

## Value

A list containing:

- p_anno:

  Annotation data frame

- p_cols:

  List of cluster colors

- plot:

  ComplexHeatmap object

- eset:

  Transformed expression matrix

## Author

Dongqiang Zeng

## Examples

``` r
tcga_stad_sig <- load_data("tcga_stad_sig")
tcga_stad_pdata <- load_data("tcga_stad_pdata")
input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
feas <- grep("MCPcounter", colnames(input), value = TRUE)
sig_pheatmap(
  input = input, feas = feas, group = "subtype",
  scale.matrix = TRUE, path = tempdir()
)
#> ℹ Heatmap palettes: 1 (pheatmap), 2 (peach), 3 (blues), 4 (virids), 5 (reds), 6 (RdBu), 7 (navy_firebrick), 8 (magma)
#> ℹ Random palettes: 1 (palette1), 2 (palette2), 3 (palette3), 4 (palette4)
#> ℹ Using random seed: 54321
#> $subtype
#>         EBV          GS         CIN         MSI 
#> "#8F7700FF"   "#386CB0" "#008B45FF"   "#e31a1c" 
#> 
#> ✔ Heatmap saved to: /tmp/Rtmpupw0yP/1-pheatmap-subtype.pdf
```
