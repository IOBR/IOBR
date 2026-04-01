# Convert Rowname To Loci

Processes a gene expression data matrix by modifying its row names.
Extracts the gene identifier from row names formatted as 'GENE\|ID',
simplifying them to 'GENE'.

## Usage

``` r
ConvertRownameToLoci(cancerGeneExpression)
```

## Arguments

- cancerGeneExpression:

  Matrix or data frame. Gene expression data with row names in the
  format 'GENE\|ID'.

## Value

Matrix with modified gene expression data with updated row names. Rows
without a valid identifier are removed.

## Examples

``` r
example_data <- matrix(runif(20), ncol = 5)
rownames(example_data) <- c("LOC101", "LOC102", "LOC103", "LOC104")
processed_data <- ConvertRownameToLoci(example_data)
print(processed_data)
#>             [,1]       [,2]      [,3]      [,4]      [,5]
#> LOC101 0.9674823 0.04420631 0.3817684 0.6047118 0.1886970
#> LOC102 0.7401439 0.92745678 0.5211051 0.2385416 0.7027989
#> LOC103 0.3484356 0.85054762 0.4146174 0.9692808 0.9780016
#> LOC104 0.3843045 0.68256832 0.6845463 0.8410889 0.9306995
```
