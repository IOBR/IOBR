# Signature Heatmap with Optional Annotations

Generates a heatmap of selected features grouped by a categorical
variable, with optional conditional (annotation) bars. Supports palette
customization, scaling, size controls, and output saving.

## Usage

``` r
sig_heatmap(
  input,
  id = "ID",
  features,
  group,
  condition = NULL,
  id_condition = "vars",
  col_condition = "condition",
  cols_condition = NULL,
  scale = FALSE,
  palette = 2,
  cols_heatmap = NULL,
  palette_group = "jama",
  show_col = FALSE,
  show_palettes = FALSE,
  cols_group = NULL,
  show_plot = TRUE,
  width = 8,
  height = NULL,
  size_col = 10,
  size_row = 8,
  angle_col = 90,
  column_title = NULL,
  row_title = NULL,
  show_heatmap_col_name = FALSE,
  path = NULL,
  index = NULL
)
```

## Arguments

- input:

  Data frame containing ID, grouping variable, and feature columns.

- id:

  Character string. Column name for sample identifier. Default is
  \`"ID"\`.

- features:

  Character vector. Feature (column) names to include in the heatmap.

- group:

  Character string. Grouping variable column name.

- condition:

  Data frame or \`NULL\`. Optional annotation table with
  variable-condition mapping. Default is \`NULL\`.

- id_condition:

  Character string. Column name in \`condition\` for feature IDs.
  Default is \`"vars"\`.

- col_condition:

  Character string. Column name in \`condition\` for condition labels.
  Default is \`"condition"\`.

- cols_condition:

  Character vector. Colors for conditions.

- scale:

  Logical indicating whether to scale values by row. Default is
  \`FALSE\`.

- palette:

  Integer or character. Palette index/name for heatmap colors. Default
  is \`2\`.

- cols_heatmap:

  Character vector. Custom colors for heatmap gradient.

- palette_group:

  Character string. Palette name for group colors. Default is
  \`"jama"\`.

- show_col:

  Logical indicating whether to display the color vector. Default is
  \`FALSE\`.

- show_palettes:

  Logical indicating whether to print palette options. Default is
  \`FALSE\`.

- cols_group:

  Character vector. Custom colors for groups.

- show_plot:

  Logical indicating whether to print the heatmap. Default is \`TRUE\`.

- width:

  Numeric. Plot width in inches. Default is \`8\`.

- height:

  Numeric or \`NULL\`. Plot height in inches. Auto-calculated if
  \`NULL\`.

- size_col:

  Numeric. Font size for column labels. Default is \`10\`.

- size_row:

  Numeric. Font size for row labels. Default is \`8\`.

- angle_col:

  Numeric. Rotation angle for column labels in degrees. Default is
  \`90\`.

- column_title:

  Character string or \`NULL\`. Title for column annotation.

- row_title:

  Character string or \`NULL\`. Title for row annotation.

- show_heatmap_col_name:

  Logical indicating whether to show column names. Default is \`FALSE\`.

- path:

  Character string or \`NULL\`. Output directory for saving the heatmap.

- index:

  Integer or \`NULL\`. Index appended to filename. Default is \`NULL\`.

## Value

A tidyHeatmap object. Saves PDF only when \`path\` is provided.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
tcga_stad_sig <- load_data("tcga_stad_sig")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "tcga_stad_sig"
tcga_stad_pdata <- load_data("tcga_stad_pdata")
input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
feas <- grep("MCPcounter", colnames(input), value = TRUE)
sig_heatmap(input = input, features = feas, group = "subtype", scale = TRUE)
#> ℹ Creating heatmap with 10 features
# }
```
