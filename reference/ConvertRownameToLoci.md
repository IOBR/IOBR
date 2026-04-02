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
#>              [,1]      [,2]      [,3]        [,4]      [,5]
#> LOC101 0.46114846 0.9781086 0.8539091 0.005264029 0.9156719
#> LOC102 0.54836452 0.0257430 0.7604734 0.955281245 0.0211167
#> LOC103 0.06104055 0.3999507 0.1678785 0.303305713 0.9013085
#> LOC104 0.33051549 0.4053564 0.3870172 0.849672277 0.9738908
```
