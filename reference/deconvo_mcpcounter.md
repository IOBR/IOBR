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
if (FALSE) { # \dontrun{
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
mcp_result <- deconvo_mcpcounter(eset = eset, project = "TCGA-STAD")
} # }
```
