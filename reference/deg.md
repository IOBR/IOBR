# Single-cell RNA-seq Differential Expression Analysis Results

Example dataset containing differential expression analysis results from
single-cell RNA sequencing (scRNA-seq) analysis. This dataset serves as
input for signature analysis functions in the IOBR package, particularly
[`get_sig_sc`](https://iobr.github.io/IOBR/reference/get_sig_sc.md) and
[`sig_gsea`](https://iobr.github.io/IOBR/reference/sig_gsea.md).

## Usage

``` r
data(deg)
```

## Format

A data frame with 3,212 rows (genes) and 7 columns (statistics):

- p_val:

  Raw p-value for differential expression test

- avg_log2FC:

  Average log2 fold-change of gene expression

- pct.1:

  Percentage of cells expressing the gene in cluster 1

- pct.2:

  Percentage of cells expressing the gene in cluster 2

- p_val_adj:

  Adjusted p-value (e.g., Bonferroni, BH FDR correction)

- cluster:

  Cell cluster or cell type identifier (e.g., "Epithelial cells 2")

- gene:

  Gene symbol (e.g., "IGFBP3", "PCDH7")

## Examples

``` r
data(deg)
head(deg)
#> # A tibble: 6 × 7
#>       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster            gene  
#>       <dbl>      <dbl> <dbl> <dbl>     <dbl> <chr>              <chr> 
#> 1 0               2.21 0.994 0.169 0         Epithelial cells 2 IGFBP3
#> 2 0               1.64 0.925 0.124 0         Epithelial cells 2 PCDH7 
#> 3 0               1.23 0.701 0.045 0         Epithelial cells 2 PTPRN2
#> 4 0               1.08 0.786 0.119 0         Epithelial cells 2 ADAM28
#> 5 6.33e-307       1.29 0.95  0.199 1.90e-303 Epithelial cells 2 CASP1 
#> 6 1.21e-280       1.32 0.897 0.159 3.62e-277 Epithelial cells 2 TRAF2 
```
