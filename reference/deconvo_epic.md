# Deconvolve Immune Microenvironment Using EPIC

Estimates immune cell fractions using EPIC algorithm.

## Usage

``` r
deconvo_epic(eset, project = NULL, tumor = TRUE)
```

## Arguments

- eset:

  Gene expression matrix with genes as row names.

- project:

  Optional project name. Default is \`NULL\`.

- tumor:

  Logical indicating tumor (\`TRUE\`) or normal (\`FALSE\`) samples.

## Value

Data frame with EPIC cell fraction estimates. Columns suffixed with
\`\_EPIC\`.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_grch38 <- load_data("anno_grch38")
#> ℹ Loading cached data: "anno_grch38"
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
eset <- eset[1:500, 1:5]
epic_result <- deconvo_epic(eset = eset, project = "Example", tumor = TRUE)
#> ℹ Running EPIC deconvolution
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "TRef"
#> Warning: there are few genes in common between the bulk samples and reference cells:471, so the data scaling might be an issue
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "mRNA_cell_default"
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
head(epic_result)
#>             ID ProjectID  Bcells_EPIC  CAFs_EPIC CD4_Tcells_EPIC
#> 1 TCGA-BR-6455   Example 1.718530e-01 0.21892728    2.202905e-01
#> 2 TCGA-BR-7196   Example 4.592485e-05 0.56234195    2.505652e-08
#> 3 TCGA-BR-8371   Example 4.505341e-07 0.06686809    2.429061e-03
#> 4 TCGA-BR-8380   Example 2.564102e-02 0.33933005    9.011961e-02
#> 5 TCGA-BR-8592   Example 1.202835e-04 0.37508784    1.278933e-02
#>   CD8_Tcells_EPIC Endothelial_EPIC Macrophages_EPIC NKcells_EPIC
#> 1    4.321134e-03        0.2925362     1.362744e-06 9.206739e-02
#> 2    9.172491e-06        0.4375160     8.671181e-05 1.773362e-07
#> 3    8.258464e-01        0.1020327     2.330986e-04 2.590220e-03
#> 4    6.812041e-02        0.4477680     2.515113e-05 2.899495e-02
#> 5    2.265081e-01        0.3822764     9.192460e-04 2.296996e-03
#>   otherCells_EPIC
#> 1    3.197966e-06
#> 2    1.850589e-09
#> 3    1.058650e-09
#> 4    7.642237e-07
#> 5    1.789082e-06
# }
```
