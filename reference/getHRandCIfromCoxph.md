# Extract Hazard Ratio and Confidence Intervals from Cox Model

Extracts hazard ratio (HR) and 95 Cox proportional hazards model.

## Usage

``` r
getHRandCIfromCoxph(coxphData)
```

## Arguments

- coxphData:

  Fitted Cox model object from \`survival::coxph()\`.

## Value

Data frame with p-values, HR, and confidence intervals.

## Author

Dorothee Nickles

Dongqiang Zeng

## Examples

``` r
library(survival)
set.seed(123)
df <- data.frame(
  TTE = rexp(200, rate = 0.1),
  Cens = rbinom(200, size = 1, prob = 0.7),
  group = sample(c("Treatment", "Control"), 200, replace = TRUE)
)
coxphData <- survival::coxph(survival::Surv(TTE, Cens) ~ group, data = df)
results <- getHRandCIfromCoxph(coxphData)
print(results)
#>                       P    HR CI_low_0.95 CI_up_0.95
#> groupTreatment 0.400853 0.868      0.6238     1.2077
```
