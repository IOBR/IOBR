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
  save_path = file.path(tempdir(), "Multiple-KM-plot"),
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

  Character string specifying the directory path for saving plots.

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
if (FALSE) { # \dontrun{
sig_stad <- load_data("sig_stad")
result <- batch_sig_surv_plot(
  input_pdata = sig_stad,
  signature = "T.cells.CD8",
  id = "ID",
  column_of_project = "ProjectID",
  project = NULL,
  time = "OS_time",
  status = "OS_status",
  time_type = "month",
  break_month = "auto",
  palette = "jama",
  cols = NULL,
  mini_sig = "score",
  show_col = TRUE,
  fig_type = "pdf"
)
} # }
```
