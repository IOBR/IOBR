# Design Custom Theme for ggplot2 Plots

Creates a customized ggplot2 theme based on user-specified parameters
for plot elements such as title size, axis sizes, legend settings, and
theme style. Supports various base themes and allows fine-tuning of
visual aspects.

## Usage

``` r
design_mytheme(
  theme = c("light", "bw", "classic", "classic2"),
  plot_title_size = 2,
  axis_title_size = 2,
  axis_text_size = 12,
  axis_angle = 60,
  hjust = 1,
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.size = 0.25,
  legend.key.height = 0.5,
  legend.key.width = 0.5,
  legend.size.text = 10,
  legend.box = "horizontal"
)
```

## Arguments

- theme:

  Base theme: "light", "bw", "classic", "classic2". Default is "light".

- plot_title_size:

  Relative size of plot title. Default is 2.

- axis_title_size:

  Relative size of axis titles. Default is 2.

- axis_text_size:

  Size of axis tick labels. Default is 12.

- axis_angle:

  Angle of x-axis tick labels. Default is 60.

- hjust:

  Horizontal justification for x-axis text. Default is 1.

- legend.position:

  Legend position: "none", "left", "right", "bottom", "top". Default is
  "bottom".

- legend.direction:

  Direction of legend items: "horizontal" or "vertical". Default is
  "horizontal".

- legend.size:

  Size of legend key. Default is 0.25.

- legend.key.height:

  Height of legend key in cm. Default is 0.5.

- legend.key.width:

  Width of legend key in cm. Default is 0.5.

- legend.size.text:

  Size of legend text labels. Default is 10.

- legend.box:

  Orientation of legend box: "horizontal" or "vertical". Default is
  "horizontal".

## Value

A ggplot2 theme object.

## Author

Dongqiang Zeng

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point()
mytheme <- design_mytheme(theme = "bw", plot_title_size = 1.5, axis_text_size = 14)
p + mytheme + ggtitle("Example Plot")
```
