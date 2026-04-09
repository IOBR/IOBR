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
# Load TCGA-STAD signature data
sig_stad <- load_data("sig_stad")
#> ℹ Loading cached data: "sig_stad"

# Test features by gender (if available in your dataset)
if ("Gender" %in% colnames(sig_stad)) {
  res <- batch_kruskal(
    data = sig_stad,
    group = "Gender",
    feature = colnames(sig_stad)[69:ncol(sig_stad)]
  )
  head(res)
}
#> ℹ Groups: 2 ("M" and "F")
#> ℹ Features: 255
#> ✔ Kruskal-Wallis test complete
#> # A tibble: 6 × 8
#>   sig_names               p.value statistic      F     M p.adj log10pvalue stars
#>   <chr>                     <dbl>     <dbl>  <dbl> <dbl> <dbl>       <dbl> <fct>
#> 1 Steroid_Hormone_Metabo… 0.00652      7.40 -0.409 0.409 0.715        2.19 **   
#> 2 Steroid_Hormone_Biosyn… 0.00693      7.29 -0.400 0.400 0.715        2.16 **   
#> 3 Linoleic_Acid_Metaboli… 0.00989      6.65 -0.268 0.268 0.715        2.00 **   
#> 4 Ascorbate_and_Aldrate_… 0.0112       6.43 -0.288 0.288 0.715        1.95 *    
#> 5 Pentose_and_Glucuronat… 0.0189       5.51 -0.294 0.294 0.853        1.72 *    
#> 6 Drug_Metabolism_by_Cyt… 0.0201       5.41 -0.422 0.422 0.853        1.70 *    
```
