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
TRef <- load_data("TRef")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "TRef"
if (!is.null(TRef)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(TRef$refProfiles) * 3), nrow(TRef$refProfiles), 3)
  rownames(sim_eset) <- rownames(TRef$refProfiles)
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_epic(eset = sim_eset, project = "Example", tumor = TRUE)
  if (!is.null(result)) head(result)
}
#> ℹ Running EPIC deconvolution
#> ℹ Loading cached data: "TRef"
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "mRNA_cell_default"
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
#>        ID ProjectID Bcells_EPIC    CAFs_EPIC CD4_Tcells_EPIC CD8_Tcells_EPIC
#> 1 Sample1   Example   0.6241581 7.681545e-05    3.384637e-06    3.721190e-06
#> 2 Sample2   Example   0.4979907 1.000478e-05    3.054977e-01    3.408085e-06
#> 3 Sample3   Example   0.2806504 6.059457e-01    1.566452e-04    3.720188e-03
#>   Endothelial_EPIC Macrophages_EPIC NKcells_EPIC otherCells_EPIC
#> 1     1.524474e-05       0.07651994 2.986568e-01    5.660059e-04
#> 2     8.786633e-03       0.17879307 8.915900e-03    2.599778e-06
#> 3     6.816456e-04       0.10884494 1.251044e-07    3.095910e-07
```
