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
mcp_genes <- load_data("mcp_genes")
if (!is.null(mcp_genes)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(mcp_genes) * 3), nrow(mcp_genes), 3)
  rownames(sim_eset) <- mcp_genes$`HUGO symbols`
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_mcpcounter(eset = sim_eset, project = "TCGA-STAD")
  if (!is.null(result)) head(result)
}
#> ℹ Running MCP-counter deconvolution
#>        ID ProjectID T_cells_MCPcounter CD8_T_cells_MCPcounter
#> 1 Sample1 TCGA-STAD         0.25454239             0.49785048
#> 2 Sample2 TCGA-STAD        -0.17518511             0.07796085
#> 3 Sample3 TCGA-STAD        -0.02209937             0.32430434
#>   Cytotoxic_lymphocytes_MCPcounter B_lineage_MCPcounter NK_cells_MCPcounter
#> 1                     -0.682678137          -0.01981956          0.14347809
#> 2                     -0.224517277          -0.19681888         -0.02269104
#> 3                      0.001387157          -0.13694623          0.03475688
#>   Monocytic_lineage_MCPcounter Myeloid_dendritic_cells_MCPcounter
#> 1                   0.12840542                          0.2068939
#> 2                  -0.08789732                          0.3109649
#> 3                   0.35663172                          0.1382419
#>   Neutrophils_MCPcounter Endothelial_cells_MCPcounter Fibroblasts_MCPcounter
#> 1             0.14153095                   0.09613642             -0.4792015
#> 2             0.02724544                   0.07829873             -0.3697520
#> 3             0.26682783                   0.11812885              0.2132765
```
