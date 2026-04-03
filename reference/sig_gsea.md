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
set.seed(123)
genes <- c(
  "TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PTEN", "APC", "RB1",
  "CDKN2A", "VHL", "ATM", "ATR", "CHEK2", "PALB2", "RAD51", "MDM2",
  "CDK4", "CDK6", "CCND1", "CCNE1", "CDK2", "E2F1", "E2F2", "E2F3",
  "ARF1", "ARF3", "ARF4", "ARF5", "ARF6", "GSK3B", "AKT1", "AKT2",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3"
)
deg <- data.frame(
  symbol = genes,
  log2FoldChange = rnorm(length(genes), mean = 0, sd = 2),
  padj = runif(length(genes), 0, 0.1)
)
signature <- list(
  DNA_Repair = c(
    "TP53", "BRCA1", "ATM",
    "ATR", "CHEK2", "PALB2", "RAD51"
  ),
  Cell_Cycle = c(
    "TP53", "MYC",
    "RB1", "CDKN2A", "CDK4",
    "CDK6", "CCND1", "CCNE1",
    "CDK2", "E2F1", "E2F2", "E2F3"
  ),
  PI3K_AKT = c(
    "AKT1", "AKT2",
    "PIK3CA", "PIK3CB", "PIK3CD",
    "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3"
  )
)
res <- sig_gsea(
  deg = deg,
  genesets = signature,
  path = tempdir(),
  show_plot = FALSE,
  print_bar = FALSE
)
#> ℹ Mapping gene symbols to Entrez IDs for "hsa"
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: number of columns of result is not a multiple of vector length (arg 1)
#> ✔ Signature data saved to sig.csv
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: 3.57% of input gene IDs are fail to map...
#> ℹ Running GSEA analysis...
#> using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
#> preparing geneSet collections...
#> GSEA analysis...
#> leading edge analysis...
#> done...
#> ✔ GSEA results written to: /tmp/RtmpEjY5c4/1-H_GSEA_significant_results.csv
#> ℹ Most significant gene sets: Cell_Cycle
print(names(res))
#> [1] "up"       "down"     "all"      "plot_top"
```
