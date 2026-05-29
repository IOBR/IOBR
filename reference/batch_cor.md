# Batch Correlation Analysis

Performs correlation analysis between a target variable and multiple
feature variables. Computes correlation coefficients, p-values, and
adjusts for multiple testing using the Benjamini-Hochberg method.

## Usage

``` r
batch_cor(data, target, feature, method = c("spearman", "pearson", "kendall"))
```

## Arguments

- data:

  Data frame containing the target and feature variables.

- target:

  Character string specifying the name of the target variable.

- feature:

  Character vector specifying the names of feature variables to
  correlate with the target.

- method:

  Character string specifying the correlation method. Options are
  \`"spearman"\`, \`"pearson"\`, or \`"kendall"\`. Default is
  \`"spearman"\`.

## Value

Tibble containing the following columns for each feature:

- sig_names:

  Feature name

- p.value:

  Raw p-value

- statistic:

  Correlation coefficient

- p.adj:

  Adjusted p-value (Benjamini-Hochberg)

- log10pvalue:

  Negative log10-transformed p-value

- stars:

  Significance stars: \*\*\*\* p\<0.0001, \*\*\* p\<0.001, \*\* p\<0.01,
  \* p\<0.05, + p\<0.5

## Author

Dongqiang Zeng

## Examples

``` r
# Create a small example dataset
set.seed(123)
data_df <- as.data.frame(matrix(runif(100 * 10), 100, 10))
colnames(data_df) <- paste0("Signature", 1:10)

# Perform batch correlation
results <- batch_cor(
  data = data_df,
  target = "Signature1",
  feature = colnames(data_df)[2:10]
)
#> ℹ Computing spearman correlation for 9 features
#> ✔ Correlation analysis complete
head(results)
#> # A tibble: 6 × 6
#>   sig_names  p.value statistic p.adj log10pvalue stars
#>   <chr>        <dbl>     <dbl> <dbl>       <dbl> <fct>
#> 1 Signature9  0.0514   -0.195  0.462       1.29  "+"  
#> 2 Signature5  0.171    -0.138  0.698       0.766 "+"  
#> 3 Signature2  0.380    -0.0887 0.698       0.420 "+"  
#> 4 Signature6  0.406    -0.0841 0.698       0.392 "+"  
#> 5 Signature4  0.438    -0.0784 0.698       0.358 "+"  
#> 6 Signature3  0.532    -0.0633 0.698       0.274 ""   
```
