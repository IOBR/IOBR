# Construct Contrast Matrix

Creates a contrast matrix for differential analysis, where each
phenotype level is contrasted against all other levels combined.

## Usage

``` r
Construct_con(pheno)
```

## Arguments

- pheno:

  Factor with different levels representing groups to contrast.

## Value

Square matrix with dimensions equal to the number of levels in
\`pheno\`. Each row represents a contrast where the corresponding level
is compared against the average of others.

## Examples

``` r
pheno <- factor(c("A", "B", "C", "D"))
contrast_matrix <- Construct_con(pheno)
print(contrast_matrix)
#>   A - others B - others C - others D - others
#> A          1         -1         -1         -1
#> B         -1          1         -1         -1
#> C         -1         -1          1         -1
#> D         -1         -1         -1          1
```
