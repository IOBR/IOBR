# estimateScore

This function reads a gene expression dataset in GCT format, calculates
enrichment scores for specific gene sets, and writes the computed scores
to an output file. It supports multiple platform types and performs
platform-specific calculations if necessary.

## Usage

``` r
estimateScore(
  input.ds,
  output.ds,
  platform = c("affymetrix", "agilent", "illumina")
)
```

## Arguments

- input.ds:

  A character string specifying the path to the input dataset file in
  GCT format. The file should have gene expression data with appropriate
  headers.

- output.ds:

  A character string specifying the path to the output dataset file,
  where the calculated scores will be written.

- platform:

  A character vector indicating the platform type. Must be one of
  "affymetrix", "agilent", or "illumina". Platform-specific calculations
  are performed based on this parameter.

## Value

This function does not return a value but writes the computed scores to
the specified output file in GCT format.

## Examples

``` r
if (FALSE) { # \dontrun{
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
eset <- tibble::rownames_to_column(eset, var = "symbol")
input_file <- tempfile(pattern = "estimate_", fileext = ".gct")
output_file <- tempfile(pattern = "estimate_score_", fileext = ".gct")
writeLines(c("#1.2", paste(nrow(eset), ncol(eset) - 1, sep = "\t")), input_file)
utils::write.table(
eset, 
input_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
estimateScore(input.ds = input_file, output.ds = output_file, platform = "affymetrix")
} # }
```
