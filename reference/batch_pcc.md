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
# \donttest{
# Load TCGA-STAD signature data
sig_stad <- load_data("sig_stad")

# Calculate partial correlations controlling for tumor purity
res <- batch_pcc(
  input = sig_stad,
  interferenceid = "TumorPurity_estimate",
  target = "Pan_F_TBRs",
  method = "pearson",
  features = colnames(sig_stad)[70:ncol(sig_stad)]
)
#> ℹ Computing pearson partial correlation for
#> ✔ Partial correlation analysis complete
head(res)
#> # A tibble: 6 × 6
#>   sig_names                    p.value statistic     p.adj log10pvalue stars
#>   <chr>                          <dbl>     <dbl>     <dbl>       <dbl> <fct>
#> 1 EMT2                       1.37e-147     0.914 3.45e-145       147.  **** 
#> 2 Normal_mucosa_Bindea_et_al 1.41e-134     0.898 1.77e-132       134.  **** 
#> 3 TGFb.myCAF                 1.07e-105     0.851 9.00e-104       105.  **** 
#> 4 CAF.S1                     1.58e- 99     0.838 9.94e- 98        98.8 **** 
#> 5 TMEscoreB_plus             3.05e- 98     0.835 1.54e- 96        97.5 **** 
#> 6 CAF_Peng_et_al             4.95e- 98     0.834 2.08e- 96        97.3 **** 
# }
```
