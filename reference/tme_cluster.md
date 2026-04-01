# Identification of TME Cluster

Performs TME (Tumor Microenvironment) clustering analysis using various
clustering methods. Supports feature selection, scaling, and automatic
determination of optimal cluster number.

## Usage

``` r
tme_cluster(
  input,
  features = NULL,
  pattern = NULL,
  id = NULL,
  scale = TRUE,
  method = "kmeans",
  min_nc = 2,
  max.nc = 6
)
```

## Arguments

- input:

  Data frame containing the input dataset.

- features:

  Vector of features to use for clustering. Default is NULL (uses all
  columns or pattern-selected columns).

- pattern:

  Regular expression pattern for selecting features. Default is NULL.

- id:

  Column name for identifiers. Default is NULL (uses row names).

- scale:

  Logical indicating whether to scale features. Default is TRUE.

- method:

  Clustering method. Default is "kmeans".

- min_nc:

  Minimum number of clusters to evaluate. Default is 2.

- max.nc:

  Maximum number of clusters to evaluate. Default is 6.

## Value

Data frame with cluster assignments appended.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
tcga_stad_sig <- load_data("tcga_stad_sig")
res <- tme_cluster(
  input = tcga_stad_sig,
  pattern = "xCell",
  id = "ID",
  method = "kmeans"
)
#> ℹ Best number of TME clusters: 3
#> ℹ Cluster distribution:
#>   1   2   3 
#> 264  13  98 
# }
```
