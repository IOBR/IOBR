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
#> ℹ Loading cached data: "eset_tme_stad"
colnames(eset_tme_stad) <- substring(colnames(eset_tme_stad), 1, 12)
pdata_sig_tme <- load_data("pdata_sig_tme")
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "pdata_sig_tme"
# \donttest{
res <- find_markers_in_bulk(
  pdata = pdata_sig_tme, eset = eset_tme_stad,
  group = "TMEcluster"
)
#> Using Seurat v5+ workflow
#> Normalizing layer: counts
#> Centering and scaling data matrix
#> Finding variable features for layer counts
#> PC_ 1 
#> Positive:  HSPB1, CNN1, ACTG2, TAGLN, ACTA2, TWIST2, RPA4, FOXF1, TNS1, TPM1 
#>     ZEB1, ROR2, VIM, TWIST1, CLDN3, PXDC1, CDKN2A, GTF2H5, CTGF, AXL 
#>     FSTL3, PER1, NEIL1, ZEB2, HLA-A, MGMT, SHFM1, FAM101B, GATA6, CDKN1A 
#> Negative:  TOPBP1, MSH2, POLE, MSH6, PRKDC, BRCA1, TDP1, RAD54L, BRIP1, DCLRE1A 
#>     FANCD2, FANCI, PALB2, GTF2H3, CTPS1, FANCM, FANCA, RIF1, E2F3, POLQ 
#>     ERCC3, XRCC5, DCLRE1B, FANCB, EXO1, RPA1, DDB1, XRCC2, CHEK1, BRCA2 
#> PC_ 2 
#> Positive:  RECQL4, RPA3, EME1, NUDT1, RAD54B, RDM1, PCNA, H2AFX, RAD51, EXO1 
#>     FEN1, CHEK2, SHFM1, RAD54L, CCNE1, XRCC2, POLQ, CLDN7, HIST1H2BL, CHEK1 
#>     FANCI, NEIL3, FANCG, BLM, CHAF1A, MUTYH, CDH1, FANCD2, FANCA, CLDN4 
#> Negative:  ZEB1, ZEB2, TNS1, ACTA2, TAGLN, AXL, VIM, ADAM19, ROR2, FOXF1 
#>     CNN1, TWIST2, TPM1, ACTG2, CTGF, PDCD1LG2, PER1, FAP, SH3PXD2A, POLK 
#>     TGFBI, FAM101B, REV3L, HAVCR2, COL4A1, ATM, ERCC4, HELQ, XPC, FSTL3 
#> PC_ 3 
#> Positive:  GZMB, LAG3, CXCL10, GZMA, PRF1, IFNG, CXCL9, CD8A, TAP1, PDCD1 
#>     TBX21, HAVCR2, HLA-C, HLA-B, TIGIT, CD274, HLA-A, B2M, PDCD1LG2, H2AFX 
#>     TAP2, MPG, CTLA4, DUT, NUDT1, NTHL1, ALKBH2, PCNA, RAD23A, RPA3 
#> Negative:  NEIL1, LIG4, MSH5, SHPRH, FAN1, POLI, ERCC6, HLTF, ATM, MRE11A 
#>     RAD52, MLH3, HUS1, TP53BP1, ERCC5, UVSSA, CLK2, REV3L, POLK, WRN 
#>     TREX1, RPA4, RAD50, GTF2H4, PMS1, MSH3, LIG3, SH3PXD2A, SPRTN, ATR 
#> PC_ 4 
#> Positive:  IFNG, B2M, GZMA, ERCC5, CD274, RAD17, CXCL9, RRM2B, CXCL10, CCNH 
#>     UBE2B, TIGIT, GZMB, MDM2, FBXW7, CTLA4, FANCL, PMS1, TBX21, LIG4 
#>     TAP2, RB1, NBN, ATM, CD8A, GEN1, POLK, ATR, WRN, PDCD1LG2 
#> Negative:  NTHL1, MPG, ERCC1, PRPF19, HSPB1, FSTL3, RAD23A, MUS81, MAD2L2, TWIST1 
#>     ERCC2, CLDN3, PNKP, ROR2, H2AFX, HIST1H2BL, SEMA7A, XAB2, SMUG1, HIST2H2BF 
#>     NUDT1, TAGLN, OGG1, POLL, LOXL2, SH3PXD2A, XRCC1, TWIST2, POLM, FOXF1 
#> PC_ 5 
#> Positive:  UBE2B, UBE2V2, UBE2N, CDK7, APEX1, MNAT1, FAP, RAD51C, POLB, TGFBI 
#>     CCNH, ADAM12, CETN2, RAD1, TDP2, ALKBH2, LOXL2, GTF2H5, XRCC4, TWIST1 
#>     XRCC6, PARP2, XRCC5, UBE2A, ERCC8, RNF8, RAD23B, FSTL3, RPA2, CHEK1 
#> Negative:  EME2, RECQL5, TBX21, PDCD1, PNKP, TIGIT, POLM, NEIL1, UVSSA, XRCC3 
#>     POLG, ATRIP, CD8A, MSH5, FBXW7, RAD9A, GATA6, PER1, ENDOV, MUTYH 
#>     CTLA4, DCLRE1C, POLD1, ATM, POLE, HLA-A, RAD52, PRF1, XAB2, DDB2 
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
# Extract top markers per cluster using base R
top_markers <- res$top_markers
# }
```
