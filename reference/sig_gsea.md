# Perform Gene Set Enrichment Analysis (GSEA)

Conducts Gene Set Enrichment Analysis to identify significantly enriched
gene sets from differential gene expression data. Supports MSigDB gene
sets or custom gene signatures, and generates comprehensive
visualizations and statistical results.

## Usage

``` r
sig_gsea(
  deg,
  genesets = NULL,
  path = NULL,
  gene_symbol = "symbol",
  logfc = "log2FoldChange",
  org = c("hsa", "mus"),
  msigdb = TRUE,
  category = "H",
  subcategory = NULL,
  palette_bar = "jama",
  palette_gsea = 2,
  cols_gsea = NULL,
  cols_bar = NULL,
  show_bar = 10,
  show_col = FALSE,
  show_plot = FALSE,
  show_gsea = 8,
  show_path_n = 20,
  plot_single_sig = FALSE,
  project = "custom_sig",
  minGSSize = 10,
  maxGSSize = 500,
  verbose = TRUE,
  seed = FALSE,
  fig.type = "pdf",
  print_bar = TRUE
)
```

## Arguments

- deg:

  Data frame containing differential expression results with gene
  symbols and log fold changes.

- genesets:

  List of custom gene sets for enrichment analysis. If \`NULL\`, MSigDB
  gene sets are used based on \`org\` and \`category\`. Default is
  \`NULL\`.

- path:

  Character string specifying the directory path for saving results.
  Default is \`NULL\`.

- gene_symbol:

  Character string specifying the column name in \`deg\` containing gene
  symbols. Default is \`"symbol"\`.

- logfc:

  Character string specifying the column name in \`deg\` containing log
  fold change values. Default is \`"log2FoldChange"\`.

- org:

  Character string specifying the organism. Options are \`"hsa"\` (Homo
  sapiens) or \`"mus"\` (Mus musculus). Default is \`"hsa"\`.

- msigdb:

  Logical indicating whether to use MSigDB gene sets. Default is
  \`TRUE\`.

- category:

  Character string specifying the MSigDB category (e.g., \`"H"\` for
  Hallmark, \`"C2"\` for curated gene sets). Default is \`"H"\`.

- subcategory:

  Character string specifying the MSigDB subcategory to filter gene
  sets. Default is \`NULL\`.

- palette_bar:

  Character string or integer specifying the color palette for bar
  plots. Default is \`"jama"\`.

- palette_gsea:

  Integer specifying the color palette for GSEA plots. Default is \`2\`.

- cols_gsea:

  Character vector specifying custom colors for GSEA enrichment plots.
  If \`NULL\`, colors are automatically generated. Default is \`NULL\`.

- cols_bar:

  Character vector specifying custom colors for the enrichment bar plot.
  If \`NULL\`, colors are automatically generated. Default is \`NULL\`.

- show_bar:

  Integer specifying the number of top enriched gene sets to display in
  the bar plot. Default is \`10\`.

- show_col:

  Logical indicating whether to display color names in the bar plot.
  Default is \`FALSE\`.

- show_plot:

  Logical indicating whether to display GSEA enrichment plots. Default
  is \`FALSE\`.

- show_gsea:

  Integer specifying the number of top significant gene sets for which
  to generate GSEA plots. Default is \`8\`.

- show_path_n:

  Integer specifying the number of pathways to display in GSEA plots.
  Default is \`20\`.

- plot_single_sig:

  Logical indicating whether to generate separate plots for each
  significant gene set. Default is \`TRUE\`.

- project:

  Character string specifying the project name for output files. Default
  is \`"custom_sig"\`.

- minGSSize:

  Integer specifying the minimum gene set size for analysis. Default is
  \`10\`.

- maxGSSize:

  Integer specifying the maximum gene set size for analysis. Default is
  \`500\`.

- verbose:

  Logical indicating whether to display progress messages. Default is
  \`TRUE\`.

- seed:

  Logical indicating whether to set a random seed for reproducibility.
  Default is \`FALSE\`.

- fig.type:

  Character string specifying the file format for saving plots (e.g.,
  \`"pdf"\`, \`"png"\`). Default is \`"pdf"\`.

- print_bar:

  Logical indicating whether to save and print the bar plot. Default is
  \`TRUE\`.

## Value

List containing:

- up:

  Data frame of up-regulated enriched gene sets

- down:

  Data frame of down-regulated enriched gene sets

- all:

  Complete GSEA results

- plot_top:

  GSEA enrichment plot for top gene sets

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
stad_group <- load_data("stad_group")
deg <- iobr_deg(
  eset = eset_stad, pdata = stad_group, group_id = "subtype",
  pdata_id = "ID", array = FALSE, method = "DESeq2",
  contrast = c("EBV", "GS"), path = file.path(tempdir(), "STAD")
)
#> ℹ Matching grouping information and expression matrix
#> ℹ Using DESeq2 for RNA-seq differential analysis
#> ! Ensure `eset` is a count expression matrix
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> ℹ Differential analysis summary:
#> ℹ   Adj.p < 0.001: 2509
#> ℹ   Adj.p < 0.05: 7446
#> ℹ   Adj.p < 0.1: 9565
#> ℹ   Adj.p < 0.25: 13933
#> ℹ Using built-in anno_grch38 for annotation
#> ℹ Group 1 = EBV: 5 samples
#> ℹ Group 2 = GS: 5 samples
#> ℹ Group 1 samples: TCGA-BR-6455, TCGA-BR-7196, TCGA-BR-8686, TCGA-BR-A4J4, TCGA-FP-7916
#> ℹ Group 2 samples: TCGA-BR-8371, TCGA-BR-8380, TCGA-BR-8592, TCGA-BR-A4IV, TCGA-BR-A4J9
#> ✔ DEG results written to: /tmp/Rtmp8wyU3P/STAD/2-DEGs.csv
signature_tme <- load_data("signature_tme")
res <- sig_gsea(deg = deg, genesets = signature_tme, path = tempdir())
#> ℹ Mapping gene symbols to Entrez IDs for "hsa"
#> 
#> 'select()' returned 1:many mapping between keys and columns
#> Warning: 39.39% of input gene IDs are fail to map...
#> Warning: number of columns of result is not a multiple of vector length (arg 1)
#> ✔ Signature data saved to sig.csv
#> 'select()' returned 1:many mapping between keys and columns
#> Warning: 5.59% of input gene IDs are fail to map...
#> ℹ Running GSEA analysis...
#> using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
#> preparing geneSet collections...
#> GSEA analysis...
#> Warning: There are ties in the preranked stats (3.41% of the list).
#> The order of those tied genes will be arbitrary, which may produce unexpected results.
#> leading edge analysis...
#> done...
#> ✔ GSEA results written to: /tmp/Rtmp8wyU3P/1-H_GSEA_significant_results.csv
#> ℹ Most significant gene sets: TMEscoreB_CIR, HLA_signature_gene, TMEscoreA_CIR, TMEscoreA_plus, Antigen_Processing_and_Presentation_Li_et_al, Natural_Killer_Cell_Cytotoxicity_Li_et_al, Normal_Fibroblast, CD8_c5_Tisg
#> Warning: Removed 27564 rows containing missing values or values outside the scale range
#> (`geom_line()`).
# }
```
