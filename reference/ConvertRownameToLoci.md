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
#>             [,1]      [,2]       [,3]      [,4]       [,5]
#> LOC101 0.7455072 0.6227921 0.98834162 0.7216793 0.08006179
#> LOC102 0.3432129 0.1607459 0.53566298 0.5830229 0.12839282
#> LOC103 0.3912402 0.2464355 0.04530769 0.5264734 0.36587539
#> LOC104 0.3521896 0.9476689 0.60226645 0.3537667 0.60014650
```
