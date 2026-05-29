# Identify Marker Features in Bulk Expression Data

Identifies informative marker features across groups from bulk gene
expression or signature score matrices using Seurat workflows. Performs
feature selection, scaling, PCA, clustering, and marker discovery.

## Usage

``` r
find_markers_in_bulk(
  pdata,
  eset,
  group,
  id_pdata = "ID",
  nfeatures = 2000,
  top_n = 20,
  thresh.use = 0.25,
  only.pos = TRUE,
  min.pct = 0.25,
  npcs = 30
)
```

## Arguments

- pdata:

  Data frame. Sample metadata.

- eset:

  Matrix. Gene expression or signature score matrix.

- group:

  Character. Column name in pdata specifying grouping variable.

- id_pdata:

  Character. Column name for sample IDs. Default is "ID".

- nfeatures:

  Integer. Number of top variable features to select. Default is 2000.

- top_n:

  Integer. Number of top markers to retain per cluster. Default is 20.

- thresh.use:

  Numeric. Threshold for marker selection. Default is 0.25.

- only.pos:

  Logical. Whether to retain only positive markers. Default is TRUE.

- min.pct:

  Numeric. Minimum expression percentage threshold. Default is 0.25.

- npcs:

  Integer. Number of principal components to use. Default is 30.

## Value

List with components: \`sce\` (Seurat object), \`markers\` (all
markers), \`top_markers\` (top markers per group).

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
  # Simulate data
  set.seed(123)
  sim_eset <- matrix(abs(rnorm(100 * 30)), 100, 30)
  rownames(sim_eset) <- paste0("Gene", 1:100)
  colnames(sim_eset) <- paste0("Sample", 1:30)
  
  sim_pdata <- data.frame(
    ID = paste0("Sample", 1:30),
    TMEcluster = rep(c("A", "B", "C"), each = 10)
  )
  
  res <- find_markers_in_bulk(
    pdata = sim_pdata, eset = sim_eset,
    group = "TMEcluster", npcs = 5
  )
  if (!is.null(res)) head(res$top_markers)
}
#> Using Seurat v5+ workflow
#> Normalizing layer: counts
#> Centering and scaling data matrix
#> Finding variable features for layer counts
#> PC_ 1 
#> Positive:  Gene81, Gene98, Gene73, Gene29, Gene99, Gene97, Gene44, Gene52, Gene14, Gene15 
#>     Gene32, Gene100, Gene36, Gene62, Gene65, Gene75, Gene91, Gene3, Gene25, Gene17 
#>     Gene83, Gene46, Gene1, Gene35, Gene13, Gene74, Gene7, Gene59, Gene82, Gene19 
#> Negative:  Gene53, Gene51, Gene78, Gene77, Gene34, Gene70, Gene22, Gene28, Gene42, Gene33 
#>     Gene88, Gene57, Gene38, Gene61, Gene87, Gene76, Gene93, Gene37, Gene66, Gene63 
#>     Gene47, Gene94, Gene30, Gene50, Gene24, Gene21, Gene5, Gene89, Gene55, Gene45 
#> PC_ 2 
#> Positive:  Gene18, Gene17, Gene85, Gene11, Gene14, Gene68, Gene5, Gene48, Gene6, Gene89 
#>     Gene84, Gene28, Gene94, Gene88, Gene91, Gene12, Gene73, Gene70, Gene61, Gene75 
#>     Gene62, Gene97, Gene40, Gene92, Gene23, Gene81, Gene1, Gene8, Gene27, Gene24 
#> Negative:  Gene96, Gene15, Gene69, Gene80, Gene39, Gene30, Gene71, Gene74, Gene64, Gene79 
#>     Gene83, Gene10, Gene13, Gene99, Gene22, Gene63, Gene82, Gene100, Gene93, Gene19 
#>     Gene7, Gene52, Gene26, Gene50, Gene38, Gene25, Gene95, Gene78, Gene33, Gene51 
#> PC_ 3 
#> Positive:  Gene54, Gene35, Gene58, Gene45, Gene3, Gene44, Gene76, Gene63, Gene49, Gene75 
#>     Gene37, Gene87, Gene8, Gene31, Gene86, Gene53, Gene25, Gene89, Gene98, Gene57 
#>     Gene79, Gene100, Gene94, Gene72, Gene47, Gene24, Gene12, Gene51, Gene5, Gene91 
#> Negative:  Gene68, Gene41, Gene50, Gene92, Gene20, Gene88, Gene36, Gene11, Gene2, Gene71 
#>     Gene21, Gene39, Gene34, Gene16, Gene29, Gene4, Gene85, Gene62, Gene55, Gene80 
#>     Gene10, Gene60, Gene59, Gene67, Gene40, Gene84, Gene82, Gene43, Gene7, Gene66 
#> PC_ 4 
#> Positive:  Gene49, Gene64, Gene8, Gene5, Gene47, Gene31, Gene27, Gene67, Gene2, Gene40 
#>     Gene74, Gene84, Gene41, Gene83, Gene81, Gene88, Gene15, Gene33, Gene71, Gene57 
#>     Gene10, Gene44, Gene91, Gene46, Gene87, Gene43, Gene100, Gene99, Gene51, Gene25 
#> Negative:  Gene26, Gene38, Gene86, Gene19, Gene4, Gene77, Gene72, Gene23, Gene9, Gene13 
#>     Gene90, Gene48, Gene58, Gene30, Gene69, Gene79, Gene93, Gene20, Gene68, Gene65 
#>     Gene97, Gene39, Gene1, Gene3, Gene17, Gene60, Gene95, Gene61, Gene45, Gene37 
#> PC_ 5 
#> Positive:  Gene63, Gene78, Gene50, Gene40, Gene62, Gene25, Gene96, Gene21, Gene95, Gene54 
#>     Gene90, Gene29, Gene87, Gene41, Gene100, Gene6, Gene17, Gene8, Gene16, Gene35 
#>     Gene57, Gene61, Gene43, Gene45, Gene11, Gene81, Gene38, Gene18, Gene39, Gene20 
#> Negative:  Gene31, Gene2, Gene19, Gene71, Gene3, Gene59, Gene37, Gene60, Gene89, Gene9 
#>     Gene46, Gene66, Gene92, Gene55, Gene22, Gene93, Gene28, Gene80, Gene24, Gene94 
#>     Gene33, Gene48, Gene75, Gene52, Gene7, Gene88, Gene73, Gene47, Gene13, Gene64 
#> Calculating cluster A
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> Calculating cluster B
#> Calculating cluster C
#> # A tibble: 3 × 7
#> # Groups:   cluster [2]
#>     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
#>     <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
#> 1 0.00597       1.51     1     1     0.597 B       Gene28
#> 2 0.00778       1.26     1     1     0.778 B       Gene61
#> 3 0.00521       1.36     1     1     0.521 C       Gene54
```
