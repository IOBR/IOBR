# Identify Marker Features in Bulk Expression Data

Identifies informative marker features across groups from bulk gene
expression or signature score matrices using Seurat workflows. Performs
feature selection, scaling, PCA, clustering, and marker discovery.

## Usage

``` r
find_markers_in_bulk(
  pdata,
  eset,
  group,
  id_pdata = "ID",
  nfeatures = 2000,
  top_n = 20,
  thresh.use = 0.25,
  only.pos = TRUE,
  min.pct = 0.25,
  npcs = 30
)
```

## Arguments

- pdata:

  Data frame. Sample metadata.

- eset:

  Matrix. Gene expression or signature score matrix.

- group:

  Character. Column name in pdata specifying grouping variable.

- id_pdata:

  Character. Column name for sample IDs. Default is "ID".

- nfeatures:

  Integer. Number of top variable features to select. Default is 2000.

- top_n:

  Integer. Number of top markers to retain per cluster. Default is 20.

- thresh.use:

  Numeric. Threshold for marker selection. Default is 0.25.

- only.pos:

  Logical. Whether to retain only positive markers. Default is TRUE.

- min.pct:

  Numeric. Minimum expression percentage threshold. Default is 0.25.

- npcs:

  Integer. Number of principal components to use. Default is 30.

## Value

List with components: \`sce\` (Seurat object), \`markers\` (all
markers), \`top_markers\` (top markers per group).

## Examples

``` r
eset_tme_stad <- load_data("eset_tme_stad")
colnames(eset_tme_stad) <- substring(colnames(eset_tme_stad), 1, 12)
pdata_sig_tme <- load_data("pdata_sig_tme")
# \donttest{
res <- find_markers_in_bulk(
  pdata = pdata_sig_tme, eset = eset_tme_stad,
  group = "TMEcluster"
)
#> Using Seurat v5+ workflow
#> Normalizing layer: counts
#> Normalizing layer: counts
#> Centering and scaling data matrix
#> Finding variable features for layer counts
#> PC_ 1 
#> Positive:  ACTA2, TAGLN, VIM, TNS1, AXL, CNN1, ZEB1, ACTG2, CTGF, TWIST2 
#>     TPM1, ZEB2, FOXF1, ROR2, ADAM19, TGFBI, PER1, FSTL3, FAM101B, HSPB1 
#>     PXDC1, IGFBP3, COL4A1, PDCD1LG2, B2M, HAVCR2, SH3PXD2A, TWIST1, CD8A, FAP 
#> Negative:  FANCI, POLQ, FANCD2, EME1, XRCC2, EXO1, RAD54L, BRIP1, GEN1, BRCA2 
#>     BLM, FANCA, RAD54B, CHEK1, RECQL4, RAD51, MSH2, FANCB, TOPBP1, BRCA1 
#>     POLE, CHEK2, FEN1, NEIL3, LIG1, RDM1, TDP1, FANCM, FANCG, RAD18 
#> PC_ 2 
#> Positive:  ATM, REV3L, SHPRH, TP53BP1, ERCC4, RECQL, FAN1, POLK, ERCC6, ATR 
#>     RIF1, ZEB2, WRN, POLI, FANCM, MLH3, PDCD1LG2, UVSSA, ADAM19, ERCC8 
#>     ATRIP, ZEB1, XPC, RAD52, HELQ, MSH6, POLH, LIG3, CTLA4, FAP 
#> Negative:  H2AFX, APEX1, PRPF19, RAD23A, NUDT1, PCNA, HLA-A, RPA3, APEX2, NTHL1 
#>     HLA-C, UBE2N, HLA-B, XRCC6, UBE2A, RAD23B, HSPB1, SHFM1, MPG, XRCC5 
#>     SMUG1, XRCC1, UNG, B2M, CETN2, ERCC1, CLDN7, FEN1, MAD2L2, ALKBH2 
#> PC_ 3 
#> Positive:  LAG3, GZMB, PRF1, CXCL10, IFNG, TBX21, CXCL9, PDCD1, TIGIT, PDCD1LG2 
#>     CD8A, HAVCR2, CD274, GZMA, CTLA4, TAP1, RAD54L, TAP2, ADAM12, FANCA 
#>     RAD51, LIG1, POLD1, CHAF1A, TWIST1, EXO1, HIST1H2BL, ADAM19, FAP, SEMA7A 
#> Negative:  XRCC5, RAD17, LIG4, RAD50, HUS1, NEIL1, GTF2H1, ERCC5, XPA, CLK2 
#>     FAN1, UBE2V2, CETN2, MSH3, RAD23B, DDB1, UBE2B, POLK, MMS19, POLM 
#>     POLI, MLH3, PMS1, XRCC6, TDP2, ERCC6, MRE11A, TPM1, ERCC3, SPRTN 
#> PC_ 4 
#> Positive:  ERCC2, POLG, MUS81, MDC1, RECQL5, POLD1, LIG3, PNKP, XRCC3, ATRIP 
#>     SH3PXD2A, POLM, HIST2H2BF, FANCE, SEMA7A, XAB2, FOXF1, ROR2, RAD9A, CLDN3 
#>     TNS1, OGG1, FSTL3, POLL, DDB1, NTHL1, RECQL4, HIST1H2BL, PRPF19, MMS19 
#> Negative:  CCNH, UBE2B, GZMA, B2M, CDK7, IFNG, RRM2B, RAD17, POLB, GZMB 
#>     NBN, FANCL, CD274, CXCL10, UBE2N, DUT, CXCL9, TAP2, PMS1, TDP2 
#>     MDM2, XRCC4, ERCC5, RB1, GTF2H3, PCNA, UBE2V2, MNAT1, GTF2H5, ERCC8 
#> PC_ 5 
#> Positive:  XPC, POLL, RNF4, XAB2, PARP3, DDB1, POLG, MMS19, DDB2, RPA1 
#>     NHEJ1, RAD23A, MPG, ERCC2, PDCD1, XRCC1, XRCC6, DCLRE1B, MDM2, PRF1 
#>     TP53, PARP1, POLH, LIG1, OGG1, CD8A, TAP1, TDP1, TP53BP1, POLD1 
#> Negative:  RPA4, CDKN2A, SHFM1, FAP, MSH5, CCNE1, RAD54B, RAD1, TWIST1, MRE11A 
#>     ADAM12, CLDN3, HLTF, GTF2H5, E2F3, CETN2, GTF2H4, NUDT1, LOXL2, ZEB1 
#>     WNT5A, COL4A1, RPA3, TWIST2, ACTA2, CLK2, XRCC2, UBE2V2, FSTL3, BLM 
#> Calculating cluster IA
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> Calculating cluster IE
#> Calculating cluster IS
# Extract top 15 markers per cluster
top15 <- res$top_markers |>
  dplyr::group_by(cluster) |>
  dplyr::top_n(15, avg_log2FC)
# }
```
