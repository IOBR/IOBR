# Deconvolve Using quanTIseq

quanTIseq deconvolution for RNA-seq immune cell fractions.

## Usage

``` r
deconvo_quantiseq(eset, project = NULL, tumor, arrays, scale_mrna)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- tumor:

  Logical: tumor samples. Must be specified.

- arrays:

  Logical: microarray data. Must be specified.

- scale_mrna:

  Logical: correct for mRNA content. Must be specified.

## Value

Data frame with quanTIseq cell fractions. Columns suffixed with
\`\_quantiseq\`.

## Author

Dongqiang Zeng

## Examples

``` r
quantiseq_data <- load_data("quantiseq_data")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "quantiseq_data"
if (!is.null(quantiseq_data)) {
  set.seed(123)
  n_sig <- nrow(quantiseq_data$TIL10_signature)
  sim_eset <- matrix(rnorm(n_sig * 3), n_sig, 3)
  rownames(sim_eset) <- rownames(quantiseq_data$TIL10_signature)
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_quantiseq(eset = sim_eset, project = "Example", tumor = TRUE,
                              arrays = FALSE, scale_mrna = FALSE)
  if (!is.null(result)) head(result)
}
#> ℹ Running quanTIseq deconvolution
#> ℹ Running quanTIseq deconvolution module
#> ℹ Loading cached data: "quantiseq_data"
#> ℹ Gene expression normalization and re-annotation (arrays: FALSE)
#> ℹ Loading cached data: "quantiseq_data"
#> ℹ Removing 17 noisy genes
#> ℹ Removing 15 genes with high expression in tumors
#> ℹ Signature genes found in data set: 138/138 (100%)
#> ℹ Mixture deconvolution (method: lsei)
#> ✔ Deconvolution successful!
#>        ID ProjectID B_cells_quantiseq Macrophages_M1_quantiseq
#> 1 Sample1   Example         0.1431669                        0
#> 2 Sample2   Example         0.0000000                        0
#> 3 Sample3   Example         0.2047434                        0
#>   Macrophages_M2_quantiseq Monocytes_quantiseq Neutrophils_quantiseq
#> 1                        0                   0             0.8568331
#> 2                        0                   0             0.0000000
#> 3                        0                   0             0.0000000
#>   NK_cells_quantiseq T_cells_CD4_quantiseq T_cells_CD8_quantiseq
#> 1          0.0000000                     0                     0
#> 2          1.0000000                     0                     0
#> 3          0.7952566                     0                     0
#>   Tregs_quantiseq Dendritic_cells_quantiseq Other_quantiseq
#> 1               0                         0               0
#> 2               0                         0               0
#> 3               0                         0               0
```
