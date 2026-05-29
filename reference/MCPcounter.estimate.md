# MCP-counter Cell Population Abundance Estimation

Estimates the abundance of different immune and stromal cell populations
using the MCP-counter method. Works with various gene identifiers
including Affymetrix probesets, HUGO gene symbols, Entrez IDs, and
Ensembl IDs.

## Usage

``` r
MCPcounter.estimate(
  expression,
  featuresType = c("affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID"),
  probesets = read.table(system.file("extdata/probesets.txt", package = "IOBR"), sep =
    "\t", stringsAsFactors = FALSE, colClasses = "character"),
  genes = read.table(system.file("extdata/genes.txt", package = "IOBR"), sep = "\t",
    stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names =
    FALSE)
)
```

## Arguments

- expression:

  Matrix or data.frame with features in rows and samples in columns.

- featuresType:

  Type of identifiers for expression features. Options:
  "affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID".
  Default is "affy133P2_probesets".

- probesets:

  Probesets data table. Default loads from GitHub.

- genes:

  Genes data table. Default loads from GitHub.

## Value

Matrix with cell populations in rows and samples in columns.

## Author

Etienne Becht

## Examples

``` r
expr <- matrix(runif(1000), nrow = 100, ncol = 10)
rownames(expr) <- paste0("Gene", 1:100)
estimates <- MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
```
