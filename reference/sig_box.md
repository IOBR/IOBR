# Signature Box Plot with Statistical Comparisons

Creates box plots to visualize signature distributions across groups
with optional statistical pairwise comparisons. Supports both data
frames and Seurat objects for single-cell data visualization.

## Usage

``` r
sig_box(
  data,
  signature,
  variable,
  palette = "nrc",
  cols = NULL,
  jitter = FALSE,
  point_size = 5,
  angle_x_text = 0,
  hjust = 0.5,
  show_pairwise_p = TRUE,
  show_overall_p = FALSE,
  return_stat_res = FALSE,
  size_of_pvalue = 6,
  size_of_font = 10,
  assay = NULL,
  slot = "scale.data",
  scale = FALSE
)
```

## Arguments

- data:

  Data frame or Seurat object containing the signature and grouping
  variable.

- signature:

  Character string specifying the column name (or feature name in
  Seurat) for the signature values to plot on the y-axis.

- variable:

  Character string specifying the grouping variable column name for the
  x-axis.

- palette:

  Character string specifying the color palette name. Default is
  \`"nrc"\`.

- cols:

  Character vector of custom fill colors. If \`NULL\`, palette is used.
  Default is \`NULL\`.

- jitter:

  Logical indicating whether to add jittered points to the box plot.
  Default is \`FALSE\`.

- point_size:

  Numeric value specifying the size of jittered points. Default is 5.

- angle_x_text:

  Numeric value specifying the rotation angle for x-axis labels (in
  degrees). Default is 0.

- hjust:

  Numeric value specifying the horizontal justification of x-axis
  labels. Default is 0.5.

- show_pairwise_p:

  Logical indicating whether to display pairwise comparison p-values
  between groups. Default is \`TRUE\`.

- show_overall_p:

  Logical indicating whether to display the overall group difference
  p-value. Default is \`FALSE\`.

- return_stat_res:

  Logical indicating whether to return statistical test results instead
  of the plot. Default is \`FALSE\`.

- size_of_pvalue:

  Numeric value specifying the font size for p-values. Default is 6.

- size_of_font:

  Numeric value specifying the base font size. Default is 10.

- assay:

  Character string specifying the assay name (for Seurat objects).
  Default is \`NULL\`.

- slot:

  Character string specifying the slot name (for Seurat objects).
  Default is \`"scale.data"\`.

- scale:

  Logical indicating whether to scale signature values (z-score
  transformation). Default is \`FALSE\`.

## Value

If \`return_stat_res = FALSE\`, returns a ggplot2 object. If
\`return_stat_res = TRUE\`, returns a data frame containing statistical
test results.

## Author

Dongqiang Zeng

## Examples

``` r
tcga_stad_pdata <- load_data("tcga_stad_pdata")
sig_box(
  data = tcga_stad_pdata,
  signature = "TMEscore_plus",
  variable = "subtype",
  jitter = TRUE,
  palette = "jco"
)
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the IOBR package.
#>   Please report the issue at <https://github.com/IOBR/IOBR/issues>.
#> # A tibble: 6 × 8
#>   .y.       group1 group2        p    p.adj p.format p.signif method  
#>   <chr>     <chr>  <chr>     <dbl>    <dbl> <chr>    <chr>    <chr>   
#> 1 signature CIN    EBV    1.81e-10 7.30e-10 1.8e-10  ****     Wilcoxon
#> 2 signature CIN    GS     3.67e- 6 7.3 e- 6 3.7e-06  ****     Wilcoxon
#> 3 signature EBV    GS     2.17e-14 1.10e-13 2.2e-14  ****     Wilcoxon
#> 4 signature CIN    MSI    3.57e-10 1.10e- 9 3.6e-10  ****     Wilcoxon
#> 5 signature EBV    MSI    1.10e- 2 1.1 e- 2 0.011    *        Wilcoxon
#> 6 signature GS     MSI    1.83e-15 1.10e-14 1.8e-15  ****     Wilcoxon
```
