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
  ID = paste0("S", 1:20),
  ProjectID = rep("P1", 20),
  OS_time = runif(20, 1, 60),
  OS_status = sample(c(0, 1), 20, replace = TRUE),
  Marker = rnorm(20)
)
result <- batch_sig_surv_plot(
  input_pdata = test_pdata, signature = "Marker",
  id = "ID", column_of_project = "ProjectID",
  time = "OS_time", status = "OS_status", time_type = "month"
)
#> ℹ Processing project: "P1"
#> ℹ Survival follow-up time range: 3.48 to 57.45 months
#> ℹ Best cutoff for "Marker": 0.69
#> ✔ Best cutoff for "Marker": 0.689
#> ℹ High Marker: 5
#> ℹ Low Marker: 15
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
if (!is.null(result)) head(result)
#>   ID      time status     Marker group3 group2 bestcutoff
#> 1 S1 17.967074      0 -1.0678237    Low    Low        Low
#> 2 S2 47.510003      1 -0.2179749 Middle    Low        Low
#> 3 S3 25.129638      0 -1.0260044    Low    Low        Low
#> 4 S4 53.098027      0 -0.7288912    Low    Low        Low
#> 5 S5 56.487570      0 -0.6250393    Low    Low        Low
#> 6 S6  3.687833      0 -1.6866933    Low    Low        Low
```
