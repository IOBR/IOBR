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
# \donttest{
# Load TCGA-STAD signature data
sig_stad <- load_data("sig_stad")
#> ℹ Loading cached data: "sig_stad"

# Compare features by gender
res <- batch_wilcoxon(
  data = sig_stad,
  target = "Gender",
  feature = colnames(sig_stad)[69:ncol(sig_stad)]
)
#> ℹ Groups: "M" and "F"
#> ℹ Features: 255
#> ✔ Wilcoxon test complete
head(res)
#> # A tibble: 6 × 8
#>   sig_names               p.value      F     M statistic p.adj log10pvalue stars
#>   <chr>                     <dbl>  <dbl> <dbl>     <dbl> <dbl>       <dbl> <fct>
#> 1 Steroid_Hormone_Metabo… 0.00653 -0.526 0.292    -0.818 0.717        2.18 **   
#> 2 Steroid_Hormone_Biosyn… 0.00694 -0.516 0.284    -0.800 0.717        2.16 **   
#> 3 Linoleic_Acid_Metaboli… 0.00991 -0.345 0.192    -0.537 0.717        2.00 **   
#> 4 Ascorbate_and_Aldrate_… 0.0112  -0.371 0.204    -0.575 0.717        1.95 *    
#> 5 Pentose_and_Glucuronat… 0.0189  -0.378 0.210    -0.588 0.854        1.72 *    
#> 6 Drug_Metabolism_by_Cyt… 0.0201  -0.548 0.295    -0.843 0.854        1.70 *    
# }
```
