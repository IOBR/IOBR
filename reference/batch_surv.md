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
# Create small example data
set.seed(123)
test_data <- data.frame(
  OS_time = runif(100, 0, 100),
  OS_status = sample(c(0, 1), 100, replace = TRUE),
  Signature1 = rnorm(100),
  Signature2 = rnorm(100)
)
batch_surv(
  pdata = test_data,
  variable = c("Signature1", "Signature2"),
  time = "OS_time",
  status = "OS_status"
)
#> # A tibble: 2 × 5
#>   ID             P    HR CI_low_0.95 CI_up_0.95
#>   <chr>      <dbl> <dbl>       <dbl>      <dbl>
#> 1 Signature2 0.317 0.862       0.644       1.15
#> 2 Signature1 0.365 0.863       0.628       1.19
```
