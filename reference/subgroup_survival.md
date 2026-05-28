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
# \donttest{
subgroup_data <- load_data("subgroup_data")
input <- subset(subgroup_data, time > 0 & !is.na(status) & !is.na(AJCC_stage))

# Binary variable example
res_bin <- subgroup_survival(
  pdata = input, time_name = "time", status_name = "status",
  variables = c("ProjectID", "AJCC_stage"), object = "score_binary"
)

# Continuous variable example
res_cont <- subgroup_survival(
  pdata = input, time_name = "time", status_name = "status",
  variables = c("ProjectID", "AJCC_stage"), object = "score"
)
# }
```
