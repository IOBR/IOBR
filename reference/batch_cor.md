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
# \donttest{
# Load TCGA-STAD signature data
sig_stad <- load_data("sig_stad")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "sig_stad"

# Perform batch correlation
results <- batch_cor(
  data = sig_stad,
  target = "CD_8_T_effector",
  feature = colnames(sig_stad)[69:ncol(sig_stad)]
)
#> ℹ Computing spearman correlation for 254 features
#> ✔ Correlation analysis complete
head(results)
#> # A tibble: 6 × 6
#>   sig_names                        p.value statistic     p.adj log10pvalue stars
#>   <chr>                              <dbl>     <dbl>     <dbl>       <dbl> <fct>
#> 1 CD8_Rooney_et_al               4.10e-234     0.971 1.04e-231        233. **** 
#> 2 TIP_Killing_of_cancer_cells_1  1.35e-232     0.971 1.72e-230        232. **** 
#> 3 Cytotoxic_cells_Danaher_et_al  5.18e-212     0.962 4.38e-210        211. **** 
#> 4 T_cell_inflamed_GEP_Ayers_et_… 5.95e-204     0.958 3.78e-202        203. **** 
#> 5 IFNG_signature_Ayers_et_al     8.63e-173     0.938 4.38e-171        172. **** 
#> 6 TMEscoreA_plus                 7.35e-154     0.920 3.11e-152        153. **** 
# }
```
