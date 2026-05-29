# Batch Kruskal-Wallis Test

Performs Kruskal-Wallis rank sum tests on multiple continuous features
across different groups. Computes p-values, adjusts for multiple
testing, and ranks features by significance.

## Usage

``` r
batch_kruskal(data, group, feature = NULL, feature_manipulation = FALSE)
```

## Arguments

- data:

  Data frame containing the dataset for analysis.

- group:

  Character string specifying the name of the grouping variable.

- feature:

  Character vector specifying the names of feature variables to test. If
  \`NULL\`, the user is prompted to select features (interactive mode
  only). Default is \`NULL\`.

- feature_manipulation:

  Logical indicating whether to apply feature manipulation to filter
  valid features. Default is \`FALSE\`.

## Value

Tibble containing:

- sig_names:

  Feature name

- p.value:

  Raw p-value from Kruskal-Wallis test

- statistic:

  Test statistic (chi-squared)

- p.adj:

  Adjusted p-value (Benjamini-Hochberg)

- log10pvalue:

  Negative log10-transformed p-value

- stars:

  Significance stars: \*\*\*\* p\<0.0001, \*\*\* p\<0.001, \*\* p\<0.01,
  \* p\<0.05, + p\<0.5

- group columns:

  Mean-centered values for each group

## Author

Dongqiang Zeng

## Examples

``` r
# Create small example data
set.seed(123)
test_data <- data.frame(
  Gender = rep(c("Male", "Female"), each = 50),
  Signature1 = rnorm(100),
  Signature2 = rnorm(100)
)
# Test features by gender
res <- batch_kruskal(
  data = test_data,
  group = "Gender",
  feature = c("Signature1", "Signature2")
)
#> ℹ Groups: 2 ("Male" and "Female")
#> ℹ Features: 2
#> ✔ Kruskal-Wallis test complete
head(res)
#> # A tibble: 2 × 8
#>   sig_names  p.value statistic Female    Male p.adj log10pvalue stars
#>   <chr>        <dbl>     <dbl>  <dbl>   <dbl> <dbl>       <dbl> <fct>
#> 1 Signature2   0.128     2.32  0.146  -0.146  0.255       0.894 +    
#> 2 Signature1   0.482     0.494 0.0560 -0.0560 0.482       0.317 +    
```
