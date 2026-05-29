# Deconvolve Using TIMER

TIMER deconvolution for cancer-specific immune estimation.

## Usage

``` r
deconvo_timer(eset, project = NULL, indications = NULL)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- indications:

  Cancer type for each sample (e.g., \`"brca"\`, \`"stad"\`). Must match
  number of columns in \`eset\`.

## Value

Data frame with TIMER cell fractions. Columns suffixed with \`\_TIMER\`.

## Author

Dongqiang Zeng

## Examples

``` r
immune <- load_data("immuneCuratedData")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "immuneCuratedData"
cancer_genes <- load_data("cancer_type_genes")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "cancer_type_genes"
if (!is.null(immune) && !is.null(cancer_genes)) {
  set.seed(123)
  # Extract a subset of genes to avoid "NA" gene names breaking read.table
  genes <- unique(c(head(rownames(immune$genes), 100), cancer_genes[["stad"]]))
  sim_eset <- matrix(rnorm(length(genes) * 3), length(genes), 3)
  rownames(sim_eset) <- genes
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_timer(
    eset = sim_eset, project = "TCGA-STAD",
    indications = rep("stad", 3)
  )
  if (!is.null(result)) head(result)
}
#> ℹ Running TIMER deconvolution
#> ℹ Enter batch mode
#> ℹ Loading immune gene expression
#> ℹ Loading cached data: "immuneCuratedData"
#> ℹ Outlier genes: AATK ABCA8 CD40 FGFBP2 ICAM2 IL3RA LST1 LTA PDCD1LG2 PLXNC1 PREX1 RASGRP2 RASGRP3 RORA SLAMF8
#> ℹ Removing batch effects for stad
#> ℹ Loading cached data: "cancer_type_genes"
#>        ID ProjectID B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER
#> 1 Sample1 TCGA-STAD    0.1044913        0.1291149        0.1947369
#> 2 Sample2 TCGA-STAD    0.1047552        0.1306022        0.1924542
#> 3 Sample3 TCGA-STAD    0.1028564        0.1271369        0.1994449
#>   Neutrophil_TIMER Macrophage_TIMER  DC_TIMER
#> 1        0.1066019       0.04875218 0.4950243
#> 2        0.1081234       0.04625540 0.4952995
#> 3        0.1083485       0.04674189 0.4949487
```
