# Remove Batch Effect of Expression Set

Removes batch effects between two gene expression datasets, typically
representing different sample types such as cancer cells and immune
cells. Uses ComBat from the sva package for batch correction.

## Usage

``` r
RemoveBatchEffect(cancer.exp, immune.exp, immune.cellType)
```

## Arguments

- cancer.exp:

  Matrix or data frame. Cancer cell expression data with genes as rows
  and samples as columns.

- immune.exp:

  Matrix or data frame. Immune cell expression data with genes as rows
  and samples as columns.

- immune.cellType:

  Vector. Cell type for each column in \`immune.exp\`.

## Value

A list containing:

- 1:

  Batch effect corrected cancer expression data

- 2:

  Batch effect corrected immune expression data

- 3:

  Aggregated immune expression data (median per cell type)

## Author

Bo Li

## Examples

``` r
set.seed(123)
gene_names <- paste0("Gene", 1:100)
sample_names_cancer <- paste0("CancerSample", 1:10)
cancer.exp <- matrix(runif(1000, 1, 1000),
  nrow = 100, ncol = 10,
  dimnames = list(gene_names, sample_names_cancer)
)

sample_names_immune <- paste0("ImmuneSample", 1:5)
immune.exp <- matrix(runif(500, 1, 1000),
  nrow = 100, ncol = 5,
  dimnames = list(gene_names, sample_names_immune)
)

immune.cellType <- c("T-cell", "B-cell", "T-cell", "NK-cell", "B-cell")
names(immune.cellType) <- sample_names_immune

result <- RemoveBatchEffect(cancer.exp, immune.exp, immune.cellType)
if (!is.null(result)) str(result)
#> List of 3
#>  $ : num [1:100, 1:10] 290 777 395 872 927 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:100] "Gene1" "Gene2" "Gene3" "Gene4" ...
#>   .. ..$ : chr [1:10] "CancerSample1" "CancerSample2" "CancerSample3" "CancerSample4" ...
#>  $ : num [1:100, 1:5] 329 554 213 811 783 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:100] "Gene1" "Gene2" "Gene3" "Gene4" ...
#>   .. ..$ : chr [1:5] "ImmuneSample1" "ImmuneSample2" "ImmuneSample3" "ImmuneSample4" ...
#>  $ : num [1:100, 1:5] 329 554 213 811 783 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:100] "Gene1" "Gene2" "Gene3" "Gene4" ...
#>   .. ..$ : chr [1:5] "ImmuneSample1" "ImmuneSample2" "ImmuneSample3" "ImmuneSample4" ...
```
