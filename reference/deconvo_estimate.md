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
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_grch38 <- load_data("anno_grch38")
#> ℹ Loading cached data: "anno_grch38"
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
# \donttest{
estimate_result <- deconvo_estimate(eset, project = "TCGA-STAD")
#> ℹ Running ESTIMATE
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "common_genes"
#> Merged dataset includes 10148 genes (264 mismatched).
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "SI_geneset"
#> 1 gene set: StromalSignature overlap=138
#> 2 gene set: ImmuneSignature overlap=140
# }
```
