# Parse Input Gene Expression Data

Reads gene expression data from a tab-delimited text file, using the
first column as row names. Converts data into a numeric matrix for
analysis.

## Usage

``` r
ParseInputExpression(path)
```

## Arguments

- path:

  Character. Path to a tab-delimited gene expression file. First column
  should contain gene identifiers.

## Value

Numeric matrix of gene expression values with genes as rows and samples
as columns.

## Examples

``` r
tf <- tempfile(fileext = ".tsv")
expr <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC"),
  Sample1 = c(10, 20, 30),
  Sample2 = c(15, 25, 35)
)
write.table(expr, tf, sep = "\t", row.names = FALSE, quote = FALSE)
gene_expression_data <- ParseInputExpression(tf)
print(gene_expression_data)
#>       Sample1 Sample2
#> GeneA      10      15
#> GeneB      20      25
#> GeneC      30      35
```
