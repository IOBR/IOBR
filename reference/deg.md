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
