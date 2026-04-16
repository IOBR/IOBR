# Analyze Mutations Related to Signature Scores

This function identifies mutations associated with a specific signature
score, performs statistical tests for significance, and generates
oncoprints and box plots to visualize relationships.

## Usage

``` r
find_mutations(
  mutation_matrix,
  signature_matrix,
  id_signature_matrix = "ID",
  signature,
  min_mut_freq = 0.05,
  plot = TRUE,
  method = "multi",
  point_alpha = 0.1,
  save_path = NULL,
  palette = "jco",
  cols = NULL,
  show_plot = TRUE,
  show_col = FALSE,
  width = 8,
  height = 4,
  oncoprint_group_by = "mean",
  oncoprint_col = "#224444",
  gene_counts = 10,
  jitter = FALSE,
  genes = NULL,
  point_size = 4.5
)
```

## Arguments

- mutation_matrix:

  A matrix of mutation data with samples in rows and genes in columns.

- signature_matrix:

  A data frame with sample identifiers and signature scores.

- id_signature_matrix:

  Column name in \`signature_matrix\` for sample identifiers.

- signature:

  Name of the target signature for analysis.

- min_mut_freq:

  Minimum mutation frequency required for gene inclusion. Default is
  0.05.

- plot:

  Logical indicating whether to generate and save plots. Default is
  TRUE.

- method:

  Statistical test method: "multi" for both Cuzick and Wilcoxon, or
  "Wilcoxon" only. Default is "multi".

- point_alpha:

  Transparency of points in box plot. Default is 0.1.

- save_path:

  Directory to save plots and results. If NULL, uses signature name.

- palette:

  Color palette for box plots(used when cols is NULL). Default is "jco".

- cols:

  Character vector. Custom colors for box plots. If NULL, uses palette.
  Default is NULL.

- show_plot:

  Logical indicating whether to display plots. Default is TRUE.

- show_col:

  Logical indicating whether to show color codes. Default is FALSE.

- width:

  Width of oncoprint plot. Default is 8.

- height:

  Height of oncoprint plot. Default is 4.

- oncoprint_group_by:

  Grouping method for oncoprint: "mean" or "quantile". Default is
  "mean".

- oncoprint_col:

  Color for mutations in oncoprint. Default is "#224444".

- gene_counts:

  Number of genes to display in oncoprint. Default is 10.

- jitter:

  Logical indicating whether to add jitter to box plot points. Default
  is FALSE.

- genes:

  Optional vector of gene names; if NULL, selects based on frequency.

- point_size:

  Size of points in box plot. Default is 4.5.

## Value

A list containing statistical test results, oncoprint plots, and box
plots.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
# This example requires a MAF file
mut_list <- make_mut_matrix(
  maf = "path_to_maf_file", isTCGA = TRUE,
  category = "multi"
)
#> -Reading
#> Error in data.table::fread(file = maf, sep = "\t", stringsAsFactors = FALSE,     verbose = FALSE, data.table = TRUE, showProgress = TRUE,     header = TRUE, fill = TRUE, skip = "Hugo_Symbol", quote = ""): File 'path_to_maf_file' does not exist or is non-readable. getwd()=='/home/runner/work/IOBR/IOBR/docs/reference'
mut <- mut_list$snp
#> Error: object 'mut_list' not found
results <- find_mutations(
  mutation_matrix = mut, signature_matrix = signature_data,
  id_signature_matrix = "ID", signature = "CD_8_T_effector",
  min_mut_freq = 0.01, plot = TRUE, method = "multi"
)
#> Error: object 'signature_data' not found
# }
```
