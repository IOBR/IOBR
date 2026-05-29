# Subgroup Survival Analysis Using Cox Proportional Hazards Models

Extracts hazard ratios (HR) and 95 proportional hazards models across
specified subgroups.

## Usage

``` r
subgroup_survival(
  pdata,
  time_name = "time",
  status_name = "status",
  variables,
  object
)
```

## Arguments

- pdata:

  Data frame containing variables, follow-up time, and outcome.

- time_name:

  Character. Column name of follow-up time. Default is \`"time"\`.

- status_name:

  Character. Column name of event status (1/0). Default is \`"status"\`.

- variables:

  Character vector. Subgrouping variables (each processed
  independently).

- object:

  Character. Variable of interest used in Cox model. If it has levels
  "High" and "Low", recode "High" to 0 and "Low" to 1 before calling.

## Value

Data frame summarizing subgroup Cox results (HR, CI, p-value).

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
test_pdata <- data.frame(
  time = runif(100, 1, 100),
  status = sample(c(0, 1), 100, replace = TRUE),
  subgroup = sample(c("A", "B"), 100, replace = TRUE),
  score = rnorm(100)
)
res <- subgroup_survival(
  pdata = test_pdata,
  time_name = "time",
  status_name = "status",
  variables = "subgroup",
  object = "score"
)
print(res)
#>                    P     HR CI_low_0.95 CI_up_0.95
#> subgroup_A 0.8258723 1.0460      0.7009     1.5608
#> subgroup_B 0.5137623 0.8618      0.5515     1.3467
```
