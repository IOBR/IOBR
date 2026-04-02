# Calculate and Visualize Correlation Matrix Between Two Variable Sets

Calculates and visualizes the correlation matrix between two sets of
variables. Supports Pearson, Spearman, and Kendall correlation methods.
The function generates a customizable heatmap with significance stars.

## Usage

``` r
get_cor_matrix(
  data,
  feas1,
  feas2,
  method = c("pearson", "spearman", "kendall"),
  path = NULL,
  index = 1,
  fig.type = "pdf",
  width = NULL,
  height = NULL,
  project = NULL,
  is.matrix = FALSE,
  scale = TRUE,
  font.size = 15,
  fill_by_cor = FALSE,
  round.num = 1,
  font.size.star = 8,
  cols = NULL
)
```

## Arguments

- data:

  Input data frame or matrix. Variables should be in columns.

- feas1:

  Character vector of variable names for the first set.

- feas2:

  Character vector of variable names for the second set.

- method:

  Correlation method: \`"pearson"\`, \`"spearman"\`, or \`"kendall"\`.
  Default is \`"pearson"\`.

- path:

  Directory to save the plot. If \`NULL\`, plot is not saved. Default is
  \`NULL\`.

- index:

  Numeric prefix for output filename. Default is 1.

- fig.type:

  File format: \`"pdf"\`, \`"png"\`, etc. Default is \`"pdf"\`.

- width:

  Plot width in inches. Auto-calculated if \`NULL\`.

- height:

  Plot height in inches. Auto-calculated if \`NULL\`.

- project:

  Project name for plot title. Default is \`NULL\`.

- is.matrix:

  Logical: if \`TRUE\`, data is transposed. Default is \`FALSE\`.

- scale:

  Logical: scale variables before correlation. Default is \`TRUE\`.

- font.size:

  Font size for axis labels. Default is 15.

- fill_by_cor:

  Logical: show correlation values instead of stars. Default is
  \`FALSE\`.

- round.num:

  Decimal places for correlation values. Default is 1.

- font.size.star:

  Font size for significance stars. Default is 8.

- cols:

  Custom colors for gradient (low, mid, high). If \`NULL\`, uses
  blue-white-red. Default is \`NULL\`.

## Value

ggplot object displaying the correlation matrix heatmap.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
data <- as.data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
colnames(data) <- paste0("Gene_", 1:10)

feas1 <- c("Gene_1", "Gene_2", "Gene_3")
feas2 <- c("Gene_4", "Gene_5", "Gene_6")

cor_plot <- get_cor_matrix(
  data = data,
  feas1 = feas1,
  feas2 = feas2,
  method = "spearman",
  project = "Example Correlation"
)
#> ℹ Calculating spearman correlation: 3 x 3
```
