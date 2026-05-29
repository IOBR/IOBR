# Time-dependent ROC Curve for Survival Analysis

Generates a time-dependent Receiver Operating Characteristic (ROC) plot
to evaluate the predictive performance of one or more variables in
survival analysis. Calculates the Area Under the Curve (AUC) for each
specified time point and variable, and creates a multi-line ROC plot
with annotated AUC values.

## Usage

``` r
roc_time(
  input,
  vars,
  time = "time",
  status = "status",
  time_point = 12,
  time_type = "month",
  palette = "jama",
  cols = "normal",
  seed = 1234,
  show_col = FALSE,
  path = NULL,
  main = "PFS",
  index = 1,
  fig.type = "pdf",
  width = 5,
  height = 5.2
)
```

## Arguments

- input:

  Data frame containing variables for analysis.

- vars:

  Character vector. Variable(s) to be evaluated.

- time:

  Character string. Name of the time variable. Default is \`"time"\`.

- status:

  Character string. Name of the status variable. Default is
  \`"status"\`.

- time_point:

  Integer or vector. Time point(s) for ROC analysis. Default is \`12\`.

- time_type:

  Character string. Time unit (\`"day"\` or \`"month"\`). Default is
  \`"month"\`.

- palette:

  Character string. Color palette for the plot. Default is \`"jama"\`.

- cols:

  Character vector or string. Color specification: \`"normal"\`,
  \`"random"\`, or custom color vector.

- seed:

  Integer. Random seed for reproducibility. Default is \`1234\`.

- show_col:

  Logical indicating whether to display the color palette. Default is
  \`FALSE\`.

- path:

  Character string or \`NULL\`. Path to save the plot. Default is
  \`NULL\`.

- main:

  Character string. Main title of the plot. Default is \`"PFS"\`.

- index:

  Integer. Index for plot filename. Default is \`1\`.

- fig.type:

  Character string. Output file type (e.g., \`"pdf"\`, \`"png"\`).
  Default is \`"pdf"\`.

- width:

  Numeric. Width of the plot. Default is \`5\`.

- height:

  Numeric. Height of the plot. Default is \`5.2\`.

## Value

A ggplot object representing the time-dependent ROC plot.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate data for offline testing
set.seed(123)
sim_input <- data.frame(
  ID = paste0("Sample", 1:100),
  time = runif(100, 1, 60),
  status = sample(0:1, 100, replace = TRUE),
  Marker1 = rnorm(100),
  Marker2 = rnorm(100),
  Marker3 = rnorm(100)
)
if (interactive()) {
  roc_time(
    input = sim_input, vars = c("Marker1", "Marker2", "Marker3"),
    time = "time", status = "status", time_point = 12, path = NULL, main = "OS"
  )
}
```
