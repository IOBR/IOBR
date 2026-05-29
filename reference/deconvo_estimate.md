# Calculate ESTIMATE Scores

Calculates stromal, immune, and ESTIMATE scores from gene expression.

## Usage

``` r
deconvo_estimate(eset, project = NULL, platform = "affymetrix")
```

## Arguments

- eset:

  Gene expression matrix with gene symbols.

- project:

  Optional project name. Default is \`NULL\`.

- platform:

  Platform type: \`"affymetrix"\` or \`"illumina"\`. Default is
  \`"affymetrix"\`.

## Value

Data frame with ESTIMATE scores. Columns suffixed with \`\_estimate\`.

## Author

Dongqiang Zeng

## Examples

``` r
common_genes <- load_data("common_genes")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "common_genes"
if (!is.null(common_genes)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(common_genes) * 3), nrow(common_genes), 3)
  rownames(sim_eset) <- common_genes$GeneSymbol
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run calculation
  result <- deconvo_estimate(sim_eset, project = "TCGA-STAD")
  if (!is.null(result)) head(result)
}
#> ℹ Running ESTIMATE
#> ℹ Loading cached data: "common_genes"
#> Merged dataset includes 10412 genes (0 mismatched).
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "SI_geneset"
#> 1 gene set: StromalSignature overlap=141
#> 2 gene set: ImmuneSignature overlap=141
#>        ID ProjectID StromalScore_estimate ImmuneScore_estimate
#> 1 Sample1 TCGA-STAD              337.3758             11.49439
#> 2 Sample2 TCGA-STAD              532.8902            814.93003
#> 3 Sample3 TCGA-STAD              418.7180            625.30846
#>   ESTIMATEScore_estimate TumorPurity_estimate
#> 1               348.8702            0.7923180
#> 2              1347.8202            0.6946727
#> 3              1044.0265            0.7260486
```
