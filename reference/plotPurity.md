# plotPurity

This function generates scatterplots of tumor purity based on ESTIMATE
scores for given samples.

## Usage

``` r
plotPurity(
  scores,
  samples = "all_samples",
  platform = c("affymetrix", "agilent", "illumina"),
  output.dir = "estimated_purity_plots"
)
```

## Arguments

- scores:

  A character string specifying the path to the input file containing
  ESTIMATE scores. The file should be a tab-separated table with
  appropriate headers.

- samples:

  A character vector specifying the sample names to plot. The default is
  "all_samples", which plots all samples in the input file.

- platform:

  A character string specifying the platform used for data collection.
  Can be "affymetrix", "agilent", or "illumina". Currently, only
  "affymetrix" is implemented.

- output.dir:

  A character string specifying the directory to save the output plots.
  The default is "estimated_purity_plots".

## Value

No return value. The function generates and saves scatterplots in the
specified output directory.

## Examples

``` r
# \donttest{
# Create a sample ESTIMATE score matrix
scores_data <- data.frame(
  Sample1 = c(100, 200, 500, 0.80),
  Sample2 = c(120, 220, 450, 0.70),
  Sample3 = c(140, 240, 600, 0.90),
  row.names = c(
    "StromalScore", "ImmuneScore", "ESTIMATEScore",
    "TumorPurity"
  ),
  check.names = FALSE
)

# Write to a temporary GCT file
scores_file <- tempfile(fileext = ".gct")
outputGCT(scores_data, scores_file)
# }
```
