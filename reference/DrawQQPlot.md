# Draw QQ Plot Comparing Cancer and Immune Expression

Creates a quantile-quantile (QQ) plot to compare gene expression
distributions between cancer and immune samples. Points along the
diagonal indicate similar distributions.

## Usage

``` r
DrawQQPlot(cancer.exp, immune.exp, name = "")
```

## Arguments

- cancer.exp:

  Vector. Gene expression data for cancer samples.

- immune.exp:

  Vector. Gene expression data for immune samples.

- name:

  Character. Optional subtitle with additional information.

## Value

Generates a QQ plot.

## Examples

``` r
cancer_exp <- rnorm(100, mean = 5, sd = 1.5)
immune_exp <- rnorm(100, mean = 5, sd = 1.5)
DrawQQPlot(
  cancer.exp = cancer_exp,
  immune.exp = immune_exp,
  name = "Comparison of Gene Expression"
)
```
