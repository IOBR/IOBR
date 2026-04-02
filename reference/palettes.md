# Select Color Palettes for Visualization

Provides curated qualitative, sequential, and diverging palettes for
multiple plot types. Supports intensity adjustment and preview.

## Usage

``` r
palettes(
  category = "box",
  palette = "nrc",
  alpha = 1,
  counts = 50,
  show_col = TRUE,
  show_message = FALSE
)
```

## Arguments

- category:

  Character. Plot/palette category: one of \`box\`, \`continue2\`,
  \`continue\`, \`random\`, \`heatmap\`, \`heatmap3\`, \`tidyheatmap\`.

- palette:

  Character or numeric. Palette name or index (varies by category).

- alpha:

  Numeric. Alpha (transparency) scaling factor. Default is 1.

- counts:

  Integer. Number of colors (for continuous palettes). Default is 50.

- show_col:

  Logical. If TRUE, prints the palette. Default is TRUE.

- show_message:

  Logical. If TRUE, prints available options. Default is FALSE.

## Value

Character vector of hex color codes.

## Author

Dongqiang Zeng

## Examples

``` r
colors <- palettes(category = "box", palette = "nrc", show_col = FALSE)
heatmap_colors <- palettes(
  category = "heatmap", palette = 1, counts = 10, show_col = FALSE
)
#> ℹ Heatmap palettes: 1 (pheatmap), 2 (peach), 3 (blues), 4 (virids), 5 (reds), 6 (RdBu), 7 (navy_firebrick), 8 (magma)
```
