# Map Score to Black and White Color

Maps a numeric input value to a color from a black-white gradient
palette. Values are mapped to a 1001-color palette where -2 maps to
black and +2 maps to white.

## Usage

``` r
mapbw(x, my_palette2 = NULL)
```

## Arguments

- x:

  Numeric value to be mapped to a color (typically between -2 and 2).

- my_palette2:

  Color palette vector (should have 1001 colors). Default uses
  black-white gradient.

## Value

A color from the black-white palette as a hex code.

## Examples

``` r
# \donttest{
my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(1001)
color <- mapbw(1.5, my_palette2)
color <- mapbw(-1, my_palette2)
# }
```
