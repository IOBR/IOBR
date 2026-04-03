# Generate Kaplan-Meier Survival Plot for Signature

Creates Kaplan-Meier survival plots for a given signature or gene, with
automatic cutoff determination. Generates three types of plots: optimal
cutoff (best cutoff), tertile-based (3 groups), and median split (2
groups).

## Usage

``` r
sig_surv_plot(
  input_pdata,
  signature,
  project = "KM",
  ID = "ID",
  time = "time",
  status = "status",
  time_type = "month",
  break_month = "auto",
  cols = NULL,
  palette = "jama",
  show_col = TRUE,
  mini_sig = "score",
  fig.type = "png",
  save_path = NULL,
  index = 1
)
```

## Arguments

- input_pdata:

  Data frame with survival data and signature scores.

- signature:

  Character string. Column name of the target signature.

- project:

  Character string. Project name for output. Default is \`"KM"\`.

- ID:

  Character string. Column name for sample IDs. Default is \`"ID"\`.

- time:

  Character string. Column name for survival time. Default is
  \`"time"\`.

- status:

  Character string. Column name for survival status. Default is
  \`"status"\`.

- time_type:

  Character string. Time unit (\`"month"\` or \`"day"\`). Default is
  \`"month"\`.

- break_month:

  Numeric or \`"auto"\`. Time axis breaks. Default is \`"auto"\`.

- cols:

  Character vector. Optional custom colors.

- palette:

  Character string. Color palette if \`cols\` not provided. Default is
  \`"jama"\`.

- show_col:

  Logical indicating whether to show colors. Default is \`TRUE\`.

- mini_sig:

  Character string. Label for low score group. Default is \`"score"\`.

- fig.type:

  Character string. File format. Default is \`"png"\`.

- save_path:

  Character string or \`NULL\`. Directory for saving plots. If \`NULL\`,
  plots are not saved. Default is \`NULL\`.

- index:

  Integer. Index for multiple plots. Default is \`1\`.

## Value

A list containing:

- data:

  Processed input data with group assignments

- plots:

  Combined survival plots

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
tcga_stad_pdata <- load_data("tcga_stad_pdata")
sig_surv_plot(
  input_pdata = tcga_stad_pdata,
  signature = "TMEscore_plus",
  time = "time",
  status = "OS_status"
)
} # }
```
