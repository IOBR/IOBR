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
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
# \donttest{
res <- deconvo_timer(
  eset = eset, project = "stad",
  indications = rep("stad", ncol(eset))
)
#> ℹ Running TIMER deconvolution
#> ℹ Enter batch mode
#> ℹ Loading immune gene expression
#> ℹ Outlier genes: ACTB ACTG1 CD74 COL1A1 EEF1A1 ERBB2 FLNA IGHG1 IGKC MT-CO1 MT-CO2 MT-ND4 MT-RNR2 MYH11
#> ℹ Removing batch effects for stad
head(res)
#>             ID ProjectID B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER
#> 1 TCGA-BR-6455      stad   0.06826546       0.21842268       0.15152999
#> 2 TCGA-BR-7196      stad   0.23999602       0.02734689       0.31100535
#> 3 TCGA-BR-8371      stad   0.09497758       0.23729463       0.00000000
#> 4 TCGA-BR-8380      stad   0.01639131       0.01427079       0.24176014
#> 5 TCGA-BR-8592      stad   0.09192814       0.27745725       0.08620296
#> 6 TCGA-BR-8686      stad   0.11028687       0.12045346       0.20750275
#>   Neutrophil_TIMER Macrophage_TIMER  DC_TIMER
#> 1       0.03064815        0.0000000 0.3631485
#> 2       0.03863240        0.2485108 0.4523019
#> 3       0.13239111        0.0000000 0.4644136
#> 4       0.18324600        0.1265858 0.4282349
#> 5       0.10411675        0.2958384 0.3811259
#> 6       0.01544118        0.0000000 0.2927147
# }
```
