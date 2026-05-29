# Batch Wilcoxon Rank-Sum Test Between Two Groups

Performs Wilcoxon rank-sum tests (Mann-Whitney U tests) to compare the
distribution of specified features between two groups. Computes
p-values, adjusts for multiple testing, and ranks features by
significance.

## Usage

``` r
batch_wilcoxon(
  data,
  target = "group",
  feature = NULL,
  feature_manipulation = FALSE
)
```

## Arguments

- data:

  Data frame containing the dataset for analysis.

- target:

  Character string specifying the column name of the grouping variable.
  Default is \`"group"\`.

- feature:

  Character vector specifying feature names to analyze. If \`NULL\`,
  prompts for selection (interactive mode only). Default is \`NULL\`.

- feature_manipulation:

  Logical indicating whether to apply feature manipulation filtering.
  Default is \`FALSE\`.

## Value

Tibble with columns:

- sig_names:

  Feature name

- p.value:

  Raw p-value

- statistic:

  Difference in means between groups

- p.adj:

  Adjusted p-value (Benjamini-Hochberg)

- log10pvalue:

  Negative log10-transformed p-value

- stars:

  Significance stars: \*\*\*\* p\<0.0001, \*\*\* p\<0.001, \*\* p\<0.01,
  \* p\<0.05, + p\<0.5

- group1, group2:

  Mean values for each group

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
# Compare features by gender
res <- batch_wilcoxon(
  data = test_data,
  target = "Gender",
  feature = c("Signature1", "Signature2")
)
#> ℹ Groups: "Male" and "Female"
#> ℹ Features: 2
#> ✔ Wilcoxon test complete
head(res)
#> # A tibble: 2 × 8
#>   sig_names  p.value Female    Male statistic p.adj log10pvalue stars
#>   <chr>        <dbl>  <dbl>   <dbl>     <dbl> <dbl>       <dbl> <fct>
#> 1 Signature2   0.128 0.0388 -0.254      0.293 0.257       0.891 +    
#> 2 Signature1   0.484 0.146   0.0344     0.112 0.484       0.315 +    
```
