# Deconvolve Immune Microenvironment Using MCP-Counter

Estimates immune cell abundances using MCP-counter.

## Usage

``` r
deconvo_mcpcounter(eset, project = NULL)
```

## Arguments

- eset:

  Gene expression matrix with HGNC symbols as row names.

- project:

  Optional project name. Default is \`NULL\`.

## Value

Data frame with MCP-counter scores. Columns suffixed with
\`\_MCPcounter\`.

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
eset <- eset[1:500, 1:3]
# \donttest{
mcp_result <- deconvo_mcpcounter(eset = eset, project = "TCGA-STAD")
#> ℹ Running MCP-counter deconvolution
head(mcp_result)
#>             ID ProjectID B_lineage_MCPcounter Endothelial_cells_MCPcounter
#> 1 TCGA-BR-6455 TCGA-STAD               113670                        15734
#> 2 TCGA-BR-7196 TCGA-STAD               952322                        46828
#> 3 TCGA-BR-8371 TCGA-STAD                 9390                         6607
#>   Fibroblasts_MCPcounter
#> 1               51296.86
#> 2              143209.86
#> 3               99398.14
# }
```
