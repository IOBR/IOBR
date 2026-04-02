# Extract Best Cutoff and Add Binary Variable to Data Frame

Finds the best cutoff point for a continuous variable in survival
analysis. Takes input data containing a continuous variable and survival
information (time and status). Returns the modified input data with a
new binary variable created based on the best cutoff point.

## Usage

``` r
best_cutoff2(
  pdata,
  variable,
  time = "time",
  status = "status",
  print_result = TRUE
)
```

## Arguments

- pdata:

  Data frame containing survival information and features.

- variable:

  Character string specifying the continuous variable name.

- time:

  Character string specifying the time-to-event column name. Default is
  \`"time"\`.

- status:

  Character string specifying the event status column name. Default is
  \`"status"\`.

- print_result:

  Logical indicating whether to print results. Default is \`TRUE\`.

## Value

List containing:

- pdata:

  Data frame with binary variable added

- best_cutoff:

  Numeric cutoff value

- cox_continuous_object:

  Cox model summary for continuous variable

- summary_binary_variable:

  Summary of binary variable

- cox_binary_object:

  Cox model summary for binary variable

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
pdata <- data.frame(
  time = rexp(100),
  status = rbinom(100, 1, 0.5),
  score = rnorm(100, mean = 50, sd = 10)
)
result <- best_cutoff2(pdata, variable = "score", print_result = FALSE)
#> ✔ Best cutoff for "score": 46.503
result$best_cutoff
#> [1] 46.5035
```
