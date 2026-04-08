# Calculate Immunophenoscore (IPS)

Calculates immune phenotype scores from gene expression data.

## Usage

``` r
deconvo_ips(eset, project = NULL, plot = FALSE)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- plot:

  Logical: generate visualization. Default is \`FALSE\`.

## Value

Data frame with IPS scores. Columns suffixed with \`\_IPS\`.

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
ips_result <- deconvo_ips(eset = eset, project = "TCGA-STAD")
#> ℹ Running IPS calculation
# }
```
