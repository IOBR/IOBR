# Transform Signature Data into List Format

Converts signature data from a data frame (with signatures as columns
and genes as rows) into a list format suitable for IOBR functions.
Handles NA values appropriately.

## Usage

``` r
format_signatures(sig_data, save_signature = FALSE, output_name = "signatures")
```

## Arguments

- sig_data:

  Data frame with signature names as columns and genes in rows. Use
  \`NA\` for missing values.

- save_signature:

  Logical. Whether to save the signature list as RData. Default is
  \`FALSE\`.

- output_name:

  Character. Output RData file name. Default is \`"signatures"\`.

## Value

List of signatures.

## Author

Dongqiang Zeng

## Examples

``` r
sig_data <- data.frame(
  Signature1 = c("Gene1", "Gene2", "Gene3", NA),
  Signature2 = c("Gene4", "Gene5", NA, NA)
)
format_signatures(sig_data)
#> ℹ There are 2 signatures
#> $Signature1
#> [1] "Gene1" "Gene2" "Gene3"
#> 
#> $Signature2
#> [1] "Gene4" "Gene5"
#> 
```
