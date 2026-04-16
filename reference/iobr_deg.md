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
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
stad_group <- load_data("stad_group")
deg <- iobr_deg(
  eset = eset_stad, pdata = stad_group,
  group_id = "subtype", pdata_id = "ID", array = FALSE,
  method = "DESeq2", contrast = c("EBV", "GS"),
  path = file.path(tempdir(), "STAD")
)
#> ℹ Matching grouping information and expression matrix
#> ℹ Using DESeq2 for RNA-seq differential analysis
#> ! Ensure `eset` is a count expression matrix
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> ℹ Differential analysis summary:
#> ℹ   Adj.p < 0.001: 2509
#> ℹ   Adj.p < 0.05: 7446
#> ℹ   Adj.p < 0.1: 9565
#> ℹ   Adj.p < 0.25: 13933
#> ℹ Using built-in anno_grch38 for annotation
#> ℹ Loading cached data: "anno_grch38"
#> ℹ Group 1 = EBV: 5 samples
#> ℹ Group 2 = GS: 5 samples
#> ℹ Group 1 samples: TCGA-BR-6455, TCGA-BR-7196, TCGA-BR-8686, TCGA-BR-A4J4, TCGA-FP-7916
#> ℹ Group 2 samples: TCGA-BR-8371, TCGA-BR-8380, TCGA-BR-8592, TCGA-BR-A4IV, TCGA-BR-A4J9
#> ✔ DEG results written to: /tmp/Rtmpq9LRNq/STAD/2-DEGs.csv
head(deg)
#> # A tibble: 6 × 21
#>   row     baseMean log2FoldChange lfcSE  stat   pvalue     padj eff_length    gc
#>   <chr>      <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>      <dbl> <dbl>
#> 1 ENSG00…    3256.           3.94 0.303  13.0 1.41e-38 4.93e-34       3156 0.519
#> 2 ENSG00…     312.           5.88 0.492  12.0 5.79e-33 1.01e-28       2494 0.510
#> 3 ENSG00…  128738.          -4.27 0.373 -11.5 2.15e-30 2.50e-26       9251 0.623
#> 4 ENSG00…  392400.          -4.95 0.436 -11.3 8.41e-30 7.33e-26       3811 0.498
#> 5 ENSG00…    3248.           2.55 0.231  11.1 1.78e-28 1.24e-24       9980 0.430
#> 6 ENSG00…    3225.           3.12 0.302  10.3 4.60e-25 2.67e-21       3545 0.521
#> # ℹ 12 more variables: entrez <dbl>, symbol <chr>, chr <chr>, start <dbl>,
#> #   end <dbl>, strand <dbl>, biotype <chr>, description <chr>, sigORnot <chr>,
#> #   label <chr>, EBV <dbl>, GS <dbl>
# }
```
