# Batch Signature Survival Plot

Generates Kaplan-Meier survival plots for multiple projects or cohorts
based on signature scores. Automatically determines optimal cutpoints
for signature stratification and creates publication-ready survival
curves.

## Usage

``` r
batch_sig_surv_plot(
  input_pdata,
  signature,
  id = "ID",
  column_of_project = "ProjectID",
  project = NULL,
  time = "time",
  status = "status",
  time_type = "day",
  break_month = "auto",
  palette = "jama",
  cols = NULL,
  mini_sig = "score",
  save_path = NULL,
  show_col = TRUE,
  fig_type = "pdf"
)
```

## Arguments

- input_pdata:

  Data frame containing survival data and signature scores.

- signature:

  Character string specifying the column name of the target signature
  for survival analysis.

- id:

  Character string specifying the column name containing unique
  identifiers. Default is \`"ID"\`.

- column_of_project:

  Character string specifying the column name containing project
  identifiers. Default is \`"ProjectID"\`.

- project:

  Character string or vector specifying project name(s) to analyze. If
  \`NULL\`, all projects are analyzed. Default is \`NULL\`.

- time:

  Character string specifying the column name containing time-to-event
  data. Default is \`"time"\`.

- status:

  Character string specifying the column name containing event status.
  Default is \`"status"\`.

- time_type:

  Character string specifying the time unit. Options are \`"day"\` or
  \`"month"\`. Default is \`"day"\`.

- break_month:

  Numeric value or \`"auto"\` specifying the interval for time axis
  breaks in months. Default is \`"auto"\`.

- palette:

  Character string specifying the color palette. Default is \`"jama"\`.

- cols:

  Character vector of custom colors. If \`NULL\`, palette is used.
  Default is \`NULL\`.

- mini_sig:

  Character string for the signature label in the legend. Default is
  \`"score"\`.

- save_path:

  Character string or \`NULL\`. Directory path for saving plots. If
  \`NULL\`, plots are not saved. Default is \`NULL\`.

- show_col:

  Logical indicating whether to display color information. Default is
  \`TRUE\`.

- fig_type:

  Character string specifying the output file format (\`"pdf"\`,
  \`"png"\`, etc.). Default is \`"pdf"\`.

## Value

Data frame containing combined survival analysis results from all
projects.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
test_pdata <- data.frame(
  ID = paste0("S", 1:30),
  ProjectID = rep(c("P1", "P2"), each = 15),
  OS_time = runif(30, 1, 60),
  OS_status = sample(c(0, 1), 30, replace = TRUE),
  Marker = rnorm(30)
)
result <- batch_sig_surv_plot(
  input_pdata = test_pdata, signature = "Marker",
  id = "ID", column_of_project = "ProjectID",
  time = "OS_time", status = "OS_status", time_type = "month"
)
#> ℹ Processing project: "P1"
#> ℹ Survival follow-up time range: 3.69 to 57.45 months
#> ℹ Best cutoff for "Marker": -0.21
#> ✔ Best cutoff for "Marker": -0.208
#> ℹ High Marker: 9
#> ℹ Low Marker: 6
#> ℹ Maximum follow-up time is 57.5 months; divided into 6 sections
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ggpubr package.
#>   Please report the issue at <https://github.com/kassambara/ggpubr/issues>.
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> ℹ Processing project: "P2"
#> ℹ Survival follow-up time range: 3.48 to 59.66 months
#> ℹ Best cutoff for "Marker": -1.55
#> ✔ Best cutoff for "Marker": -1.549
#> ℹ High Marker: 14
#> ℹ Low Marker: 1
#> ℹ Maximum follow-up time is 59.7 months; divided into 6 sections
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
#> Ignoring unknown labels:
#> • colour : "Strata"
```
