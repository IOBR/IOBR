# Identify High-Variance Features from Statistical Results

Selects top variable (up- and down-regulated) features based on adjusted
p-value and log fold-change thresholds.

## Usage

``` r
high_var_fea(
  result,
  target,
  name_padj = "padj",
  padj_cutoff = 1,
  name_logfc,
  logfc_cutoff = 0,
  n = 10,
  data_type = NULL
)
```

## Arguments

- result:

  Data frame or tibble. Statistical results containing feature, adjusted
  p-value, and logFC columns.

- target:

  Character. Column name of feature identifiers.

- name_padj:

  Character. Adjusted p-value column name. Default is \`"padj"\`.

- padj_cutoff:

  Numeric. Adjusted p-value threshold. Default is 1.

- name_logfc:

  Character. log2 fold-change column name.

- logfc_cutoff:

  Numeric. Absolute log2 fold-change threshold. Default is 0.

- n:

  Integer. Number of top up and top down features to select. Default is
  10.

- data_type:

  Character or \`NULL\`. If \`"survival"\`, adjusts logFC
  interpretation. Default is \`NULL\`.

## Value

Character vector of selected feature names (combined up and down sets).

## Author

Dongqiang Zeng

## Examples

``` r
result_data <- data.frame(
  gene = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
  padj = c(0.01, 0.02, 0.05, 0.001, 0.03),
  logfc = c(-2, 1.5, -3, 2.5, 0.5)
)
high_var_fea(
  result = result_data,
  target = "gene",
  name_padj = "padj",
  name_logfc = "logfc",
  n = 2,
  padj_cutoff = 0.05,
  logfc_cutoff = 1.5
)
#> ! Cutoff too strict for down-regulated features, only 1 found
#> ! Cutoff too strict for up-regulated features, only 1 found
#> [1] "Gene1" "Gene4"
```
