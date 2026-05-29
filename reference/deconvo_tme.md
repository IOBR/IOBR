# Main TME Deconvolution Function

Unified interface for multiple TME deconvolution methods.

## Usage

``` r
deconvo_tme(
  eset,
  project = NULL,
  method = tme_deconvolution_methods,
  arrays = FALSE,
  tumor = TRUE,
  perm = 1000,
  reference,
  scale_reference = TRUE,
  plot = FALSE,
  scale_mrna = TRUE,
  group_list = NULL,
  platform = "affymetrix",
  absolute.mode = FALSE,
  abs.method = "sig.score",
  ...
)
```

## Arguments

- eset:

  Gene expression matrix with HGNC symbols as row names.

- project:

  Optional project name. Default is \`NULL\`.

- method:

  Deconvolution method. See \[tme_deconvolution_methods\].

- arrays:

  Logical: microarray-optimized mode. Default is \`FALSE\`.

- tumor:

  Logical: tumor-optimized mode (EPIC). Default is \`TRUE\`.

- perm:

  Permutations (CIBERSORT/SVR). Default is 1000.

- reference:

  Custom reference matrix (SVR/lsei).

- scale_reference:

  Logical: scale reference (SVR/lsei).

- plot:

  Logical: generate plots (IPS). Default is \`FALSE\`.

- scale_mrna:

  Logical: mRNA correction (quanTIseq/EPIC).

- group_list:

  Cancer types for TIMER (vector).

- platform:

  Platform for ESTIMATE. Default is \`"affymetrix"\`.

- absolute.mode:

  Logical: absolute mode (CIBERSORT/SVR). Default is \`FALSE\`.

- abs.method:

  Absolute mode method. Default is \`"sig.score"\`.

- ...:

  Additional arguments passed to method.

## Value

Tibble with cell fractions and \`ID\` column.

## References

1.  Newman et al. (2015). Robust enumeration of cell subsets from tissue
    expression profiles. Nature Methods.

2.  Vegesna et al. (2013). Inferring tumour purity and stromal/immune
    cell admixture. Nature Communications.

3.  Finotello et al. (2019). Molecular and pharmacological modulators of
    the tumor immune contexture. Genome Medicine.

4.  Li et al. (2016). Comprehensive analyses of tumor immunity. Genome
    Biology.

5.  Charoentong et al. (2017). Pan-cancer Immunogenomic Analyses. Cell
    Reports.

6.  Becht et al. (2016). Estimating population abundance of
    tissue-infiltrating immune cells. Genome Biology.

7.  Aran et al. (2017). xCell: digitally portraying tissue cellular
    heterogeneity. Genome Biology.

8.  Racle et al. (2017). Simultaneous enumeration of cancer and immune
    cell types. ELife.

## Author

Dongqiang Zeng, Rongfang Shen

## Examples

``` r
mcp_genes <- load_data("mcp_genes")
if (!is.null(mcp_genes)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(mcp_genes) * 3), nrow(mcp_genes), 3)
  rownames(sim_eset) <- mcp_genes$`HUGO symbols`
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_tme(eset = sim_eset, method = "mcpcounter")
  if (!is.null(result)) head(result)
}
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running MCP-counter deconvolution
#> # A tibble: 3 × 11
#>   ID      T_cells_MCPcounter CD8_T_cells_MCPcounter Cytotoxic_lymphocytes_MCPc…¹
#>   <chr>                <dbl>                  <dbl>                        <dbl>
#> 1 Sample1             0.255                  0.498                      -0.683  
#> 2 Sample2            -0.175                  0.0780                     -0.225  
#> 3 Sample3            -0.0221                 0.324                       0.00139
#> # ℹ abbreviated name: ¹​Cytotoxic_lymphocytes_MCPcounter
#> # ℹ 7 more variables: B_lineage_MCPcounter <dbl>, NK_cells_MCPcounter <dbl>,
#> #   Monocytic_lineage_MCPcounter <dbl>,
#> #   Myeloid_dendritic_cells_MCPcounter <dbl>, Neutrophils_MCPcounter <dbl>,
#> #   Endothelial_cells_MCPcounter <dbl>, Fibroblasts_MCPcounter <dbl>
```
