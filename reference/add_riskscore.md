# Add Risk Score to Dataset

Computes a risk score for each observation based on Cox proportional
hazards regression or binary logistic regression. The function fits the
specified model and returns the dataset with an added risk score column.

## Usage

``` r
add_riskscore(
  input,
  family = c("cox", "binary"),
  target = NULL,
  time = NULL,
  status = NULL,
  vars,
  new_var_name = "riskscore"
)
```

## Arguments

- input:

  Data frame containing the variables for analysis.

- family:

  Character string specifying the model family: \`"cox"\` for Cox
  proportional hazards regression or \`"binary"\` for logistic
  regression. Default is \`"cox"\`.

- target:

  Character string specifying the target variable name. Required when
  \`family = "binary"\`.

- time:

  Character string specifying the time-to-event variable name. Required
  when \`family = "cox"\`.

- status:

  Character string specifying the event status variable name. Required
  when \`family = "cox"\`.

- vars:

  Character vector of variable names to include in the model.

- new_var_name:

  Character string specifying the name for the new risk score column.
  Default is \`"riskscore"\`.

## Value

Data frame identical to \`input\` with an additional column containing
risk scores (linear predictors for Cox models or predicted probabilities
for logistic models).

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
if (requireNamespace("survival", quietly = TRUE)) {
  lung <- survival::lung
  result <- add_riskscore(
    lung,
    time = "time", status = "status",
    vars = c("age", "sex")
  )
  head(result)
}
#>   inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
#> 1    3  306      2  74   1       1       90       100     1175      NA
#> 2    3  455      2  68   1       0       90        90     1225      15
#> 3    3 1010      1  56   1       0       90        90       NA      15
#> 4    5  210      2  57   1       1       90        60     1150      11
#> 5    1  883      2  60   1       0      100        90       NA       0
#> 6   12 1022      1  74   1       1       50        80      513       0
#>    riskscore
#> 1 0.39950470
#> 2 0.29723270
#> 3 0.09268872
#> 4 0.10973405
#> 5 0.16087005
#> 6 0.39950470
# }
```
