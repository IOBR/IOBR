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
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
eset <- eset[1:500, 1:3]
# \donttest{
res <- deconvo_quantiseq(
  eset = eset, project = "stad", tumor = TRUE,
  arrays = FALSE, scale_mrna = FALSE
)
#> ℹ Running quanTIseq deconvolution
#> ℹ Running quanTIseq deconvolution module
#> ℹ Gene expression normalization and re-annotation (arrays: FALSE)
#> ℹ Removing 17 noisy genes
#> ℹ Removing 15 genes with high expression in tumors
#> ℹ Signature genes found in data set: 0/138 (0%)
#> ℹ Mixture deconvolution (method: lsei)
#> ✔ Deconvolution successful!
head(res)
#>             ID ProjectID B_cells_quantiseq Macrophages_M1_quantiseq
#> 1 TCGA-BR-6455      stad               NaN                      NaN
#> 2 TCGA-BR-7196      stad               NaN                      NaN
#> 3 TCGA-BR-8371      stad               NaN                      NaN
#>   Macrophages_M2_quantiseq Monocytes_quantiseq Neutrophils_quantiseq
#> 1                      NaN                 NaN                   NaN
#> 2                      NaN                 NaN                   NaN
#> 3                      NaN                 NaN                   NaN
#>   NK_cells_quantiseq T_cells_CD4_quantiseq T_cells_CD8_quantiseq
#> 1                NaN                   NaN                   NaN
#> 2                NaN                   NaN                   NaN
#> 3                NaN                   NaN                   NaN
#>   Tregs_quantiseq Dendritic_cells_quantiseq Other_quantiseq
#> 1             NaN                       NaN             NaN
#> 2             NaN                       NaN             NaN
#> 3             NaN                       NaN             NaN
# }
```
