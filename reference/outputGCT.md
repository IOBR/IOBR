# outputGCT

This function converts a gene expression dataset to a GCT format file.

## Usage

``` r
outputGCT(input.f, output.f)
```

## Arguments

- input.f:

  A data frame or a character string specifying the path to the input
  file. If a character string, the file should be a tab-separated table
  with row names.

- output.f:

  A character string specifying the path to the output file.

## Value

No return value. The function writes the dataset to the specified output
file in GCT format.

## Examples

``` r
# Create a sample input data frame
sample_data <- data.frame(
  Gene = c("BRCA1", "TP53", "EGFR"),
  Sample1 = c(10, 15, 8),
  Sample2 = c(12, 18, 7),
  stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_data$Gene
sample_data <- sample_data[, -1]

# Convert the input data frame to GCT format and save to temporary file
output_file <- tempfile(fileext = ".gct")
outputGCT(sample_data, output_file)
```
