# Batch Calculation of Partial Correlation Coefficients

Computes partial correlation coefficients between multiple features and
a target variable while controlling for an interference (confounding)
variable. Adjusts p-values for multiple testing using the
Benjamini-Hochberg method.

## Usage

``` r
batch_pcc(
  input,
  interferenceid,
  target,
  features,
  method = c("pearson", "spearman", "kendall")
)
```

## Arguments

- input:

  Data frame containing feature variables, target variable, and
  interference variable.

- interferenceid:

  Character string specifying the column name of the interference
  (confounding) variable to control for.

- target:

  Character string specifying the column name of the target variable.

- features:

  Character vector specifying the column names of feature variables to
  correlate with the target.

- method:

  Character string specifying the correlation method. Options are
  \`"pearson"\`, \`"spearman"\`, or \`"kendall"\`. Default is
  \`"pearson"\`.

## Value

Tibble containing the following columns for each feature:

- sig_names:

  Feature name

- p.value:

  Raw p-value

- statistic:

  Partial correlation coefficient

- p.adj:

  Adjusted p-value (Benjamini-Hochberg method)

- log10pvalue:

  Negative log10-transformed p-value

- stars:

  Significance stars: \*\*\*\* p.adj\<0.0001, \*\*\* p.adj\<0.001, \*\*
  p.adj\<0.01, \* p.adj\<0.05, + p.adj\<0.5

## Author

Rongfang Shen

## Examples

``` r
# Create small example data
set.seed(123)
test_data <- data.frame(
  TumorPurity = runif(100),
  TargetVar = rnorm(100),
  Signature1 = rnorm(100),
  Signature2 = rnorm(100)
)
# Calculate partial correlations controlling for tumor purity
res <- batch_pcc(
  input = test_data,
  interferenceid = "TumorPurity",
  target = "TargetVar",
  method = "pearson",
  features = c("Signature1", "Signature2")
)
#> ℹ Computing pearson partial correlation for
#> ✔ Partial correlation analysis complete
head(res)
#> # A tibble: 2 × 6
#>   sig_names  p.value statistic p.adj log10pvalue stars
#>   <chr>        <dbl>     <dbl> <dbl>       <dbl> <fct>
#> 1 Signature1   0.387   -0.0878 0.775      0.412  ""   
#> 2 Signature2   0.877   -0.0158 0.877      0.0570 ""   
```
