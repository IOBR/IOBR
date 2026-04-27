# Batch Survival Analysis

Performs Cox proportional hazards regression analysis on multiple
variables. Optionally determines optimal cutoffs to dichotomize
continuous predictors before modeling. Returns hazard ratios, confidence
intervals, and p-values for each variable.

## Usage

``` r
batch_surv(
  pdata,
  variable,
  time = "time",
  status = "status",
  best_cutoff = FALSE
)
```

## Arguments

- pdata:

  Data frame containing survival time, event status, and predictor
  variables.

- variable:

  Character vector specifying the names of predictor variables to
  analyze.

- time:

  Character string specifying the column name containing follow-up time.
  Default is \`"time"\`.

- status:

  Character string specifying the column name containing event status (1
  = event occurred, 0 = censored). Default is \`"status"\`.

- best_cutoff:

  Logical indicating whether to compute optimal cutoffs for continuous
  variables and analyze dichotomized versions. Default is \`FALSE\`.

## Value

Data frame containing hazard ratios (HR), 95 and p-values for each
variable, sorted by p-value.

## Author

Dongqiang Zeng

## Examples

``` r
sig_stad <- load_data("sig_stad")
#> ℹ Loading cached data: "sig_stad"
batch_surv(
  pdata = sig_stad,
  variable = colnames(sig_stad)[69:ncol(sig_stad)],
  time = "OS_time",
  status = "OS_status"
)
#> # A tibble: 255 × 5
#>    ID                                        P    HR CI_low_0.95 CI_up_0.95
#>    <chr>                                 <dbl> <dbl>       <dbl>      <dbl>
#>  1 Taurine_and_Hypotaurine_Metabolism 0.000420  1.21        1.09       1.35
#>  2 TGFb_Family_Member_Li_et_al        0.00263   1.10        1.03       1.17
#>  3 Cytokines_Li_et_al                 0.00273   1.03        1.01       1.06
#>  4 Retinoic_Acid_Metabolism           0.00273   1.13        1.04       1.23
#>  5 CAF.S1                             0.00278   1.04        1.01       1.06
#>  6 Arachidonic_Acid_Metabolism        0.00305   1.09        1.03       1.16
#>  7 GPAGs                              0.00353   1.09        1.03       1.15
#>  8 PPARgama_target_genes              0.00354   1.21        1.06       1.37
#>  9 EMT2                               0.00358   1.13        1.04       1.23
#> 10 MDSC_Peng_et_al                    0.00393   1.09        1.03       1.16
#> # ℹ 245 more rows
```
