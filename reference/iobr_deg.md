# Differential Expression Analysis

Performs differential expression analysis on gene expression data using
either DESeq2 or limma. Includes pre-processing steps like filtering low
count data, and calculates fold changes and adjusted p-values.
Optionally generates volcano plots and heatmaps.

## Usage

``` r
iobr_deg(
  eset,
  annotation = NULL,
  id_anno = NULL,
  pdata,
  group_id = "group",
  pdata_id = "ID",
  array = FALSE,
  method = c("DESeq2", "limma"),
  contrast = c("High", "Low"),
  path = NULL,
  padj_cutoff = 0.01,
  logfc_cutoff = 0.5,
  volcano_plot = FALSE,
  col_volcano = 1,
  heatmap = TRUE,
  col_heatmap = 1,
  parallel = FALSE
)
```

## Arguments

- eset:

  A matrix of gene expression data where rows represent genes and
  columns represent samples.

- annotation:

  Optional data frame for mapping gene IDs to gene names. Default is
  \`NULL\`.

- id_anno:

  Character string specifying the identifier column in annotation.
  Default is \`NULL\`.

- pdata:

  A data frame containing sample information and grouping labels.

- group_id:

  Character string specifying the column name in \`pdata\` containing
  grouping labels. Default is \`"group"\`.

- pdata_id:

  Character string specifying the column name in \`pdata\` for sample
  IDs. Default is \`"ID"\`.

- array:

  Logical indicating whether to perform quantile normalization. Default
  is \`FALSE\`.

- method:

  Character string specifying the method: \`"DESeq2"\` or \`"limma"\`.
  Default is \`"DESeq2"\`.

- contrast:

  Character vector of length 2 specifying contrast groups. Default is
  \`c("High", "Low")\`.

- path:

  Character string for output directory. Default is \`NULL\`.

- padj_cutoff:

  Numeric cutoff for adjusted p-values. Default is \`0.01\`.

- logfc_cutoff:

  Numeric log2 fold change cutoff. Default is \`0.5\`.

- volcano_plot:

  Logical indicating whether to generate a volcano plot. Default is
  \`FALSE\`.

- col_volcano:

  Integer specifying color index for volcano plot. Default is \`1\`.

- heatmap:

  Logical indicating whether to generate a heatmap. Default is \`TRUE\`.

- col_heatmap:

  Integer specifying color index for heatmap. Default is \`1\`.

- parallel:

  Logical indicating whether to run in parallel. Default is \`FALSE\`.

## Value

Data frame containing differentially expressed genes with statistics
including log2 fold changes and adjusted p-values.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate data
set.seed(123)
sim_eset <- matrix(abs(rnorm(100 * 20)), 100, 20)
rownames(sim_eset) <- paste0("Gene", 1:100)
colnames(sim_eset) <- paste0("Sample", 1:20)

sim_pdata <- data.frame(
  ID = paste0("Sample", 1:20),
  group = rep(c("High", "Low"), each = 10)
)

# Run DEG analysis
deg <- iobr_deg(
  eset = sim_eset, pdata = sim_pdata,
  group_id = "group", pdata_id = "ID",
  method = "limma", contrast = c("High", "Low"),
  heatmap = FALSE
)
#> ℹ Matching grouping information and expression matrix
#> ℹ Using limma for array differential analysis
#> ℹ Group 1 = High
#> ℹ Group 2 = Low
if (!is.null(deg)) head(deg)
#> # A tibble: 6 × 11
#>   symbol log2FoldChange AveExpr     t  pvalue  padj     B sigORnot label    High
#>   <chr>           <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <chr>    <chr>   <dbl>
#> 1 Gene29          0.813   0.808  3.27 0.00170 0.170 -1.25 NOT      log2FC… 1.21 
#> 2 Gene16          0.697   1.00   2.56 0.0129  0.430 -2.94 NOT      log2FC… 1.35 
#> 3 Gene96          0.686   0.732  2.57 0.0123  0.430 -2.90 NOT      log2FC… 1.07 
#> 4 Gene43          0.605   0.903  2.23 0.0294  0.734 -3.60 NOT      log2FC… 1.21 
#> 5 Gene28         -0.544   0.775 -2.10 0.0400  0.800 -3.84 NOT      log2FC… 0.503
#> 6 Gene12         -0.528   0.927 -1.89 0.0632  0.851 -4.20 NOT      log2FC… 0.663
#> # ℹ 1 more variable: Low <dbl>
```
