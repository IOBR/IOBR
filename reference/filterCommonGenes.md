# filterCommonGenes

This function filters and merges a dataset with a set of common genes.

## Usage

``` r
filterCommonGenes(input.f, output.f, id = c("GeneSymbol", "EntrezID"))
```

## Arguments

- input.f:

  A character string specifying the path to the input file or a
  connection object. The file should be a tab-separated table with row
  names.

- output.f:

  A character string specifying the path to the output file.

- id:

  A character string indicating the type of gene identifier to use. Can
  be either "GeneSymbol" or "EntrezID".

## Value

No return value. The function writes the merged dataset to the specified
output file.

## Examples

``` r
# \donttest{
# Create a sample input dataframe
input_data <- data.frame(
  GeneSymbol = c("BRCA1", "TP53", "EGFR", "NOTCH1"),
  Value = c(10, 15, 8, 12),
  stringsAsFactors = FALSE
)

# Write the input data to temporary file
input_file <- tempfile(fileext = ".txt")
output_file <- tempfile(fileext = ".txt")
write.table(input_data,
  file = input_file, sep = "\t", row.names = TRUE,
  quote = FALSE
)

# Call the filterCommonGenes function
filterCommonGenes(input_file, output_file, id = "GeneSymbol")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "common_genes"
#> Merged dataset includes 0 genes (10412 mismatched).
# }
```
