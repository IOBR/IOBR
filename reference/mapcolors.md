# Map Score to Color

Maps a numeric input value to a color from a blue-white-red gradient
palette. Values are mapped to a 1001-color palette where -3 maps to
blue, 0 maps to white, and +3 maps to red.

## Usage

``` r
mapcolors(x, my_palette = NULL)
```

## Arguments

- x:

  Numeric value to be mapped to a color (typically between -3 and 3).

- my_palette:

  Color palette vector (should have 1001 colors). Default uses
  blue-white-red gradient.

## Value

A color from the palette as a hex code.

## Examples

``` r
my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(1001)
color <- mapcolors(2, my_palette)
color <- mapcolors(-2, my_palette)
```
