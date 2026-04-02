# Gene signature collection for pathway and immune analysis

A named list of gene signatures used in the IOBR package for immune
deconvolution, pathway scoring, functional annotation, and tumour
microenvironment (TME) characterization. Each element corresponds to a
predefined biological signature and contains a character vector of HGNC
gene symbols.

## Usage

``` r
data(signature_collection)
```

## Format

A named list of length 323. Each element is a character vector of gene
symbols. Representative entries include:

- CD_8_T_effector:

  Markers of CD8\\^{+}\\ effector T cells.

- DDR:

  DNA damage response and repair genes.

- Immune_Checkpoint:

  Immune checkpoint molecules.

- CellCycle_Reg:

  Core regulators of cell-cycle progression.

- Mismatch_Repair:

  Mismatch-repair pathway genes.

- TMEsocreA_CIR:

  TME-related signature used in TMEscore analysis.

- ...:

  Additional signatures are included in the list but are not
  individually listed here; all follow the same structure.

## Examples

``` r
data(signature_collection)
head(signature_collection)
#> $CD_8_T_effector
#> [1] "CD8A"   "GZMA"   "GZMB"   "IFNG"   "CXCL9"  "CXCL10" "PRF1"   "TBX21" 
#> 
#> $DDR
#>   [1] "UNG"     "SMUG1"   "MBD4"    "OGG1"    "MUTYH"   "NTHL1"   "MPG"    
#>   [8] "NEIL1"   "NEIL2"   "NEIL3"   "APEX1"   "APEX2"   "LIG3"    "XRCC1"  
#>  [15] "PNKP"    "APLF"    "PARP1"   "PARP2"   "PARP3"   "MGMT"    "ALKBH2" 
#>  [22] "ALKBH3"  "TDP1"    "TDP2"    "MSH2"    "MSH3"    "MSH6"    "MLH1"   
#>  [29] "PMS2"    "MSH4"    "MSH5"    "MLH3"    "PMS1"    "XPC"     "RAD23B" 
#>  [36] "CETN2"   "RAD23A"  "XPA"     "DDB1"    "DDB2"    "RPA1"    "RPA2"   
#>  [43] "RPA3"    "ERCC3"   "ERCC2"   "GTF2H1"  "GTF2H2"  "GTF2H3"  "GTF2H4" 
#>  [50] "GTF2H5"  "CDK7"    "CCNH"    "MNAT1"   "ERCC5"   "ERCC1"   "ERCC4"  
#>  [57] "LIG1"    "ERCC8"   "ERCC6"   "UVSSA"   "XAB2"    "MMS19"   "RAD51"  
#>  [64] "RAD51B"  "RAD51D"  "DMC1"    "XRCC2"   "XRCC3"   "RAD52"   "RAD54L" 
#>  [71] "RAD54B"  "BRCA1"   "SHFM1"   "RAD50"   "MRE11A"  "NBN"     "RBBP8"  
#>  [78] "MUS81"   "EME1"    "EME2"    "GEN1"    "FANCA"   "FANCB"   "FANCC"  
#>  [85] "BRCA2"   "FANCD2"  "FANCE"   "FANCF"   "FANCG"   "FANCI"   "BRIP1"  
#>  [92] "FANCL"   "FANCM"   "PALB2"   "RAD51C"  "XRCC6"   "XRCC5"   "PRKDC"  
#>  [99] "LIG4"    "XRCC4"   "DCLRE1C" "NHEJ1"   "NUDT1"   "DUT"     "RRM2B"  
#> [106] "POLB"    "POLG"    "POLD1"   "POLE"    "PCNA"    "REV3L"   "MAD2L2" 
#> [113] "POLH"    "POLI"    "POLQ"    "POLK"    "POLL"    "POLM"    "POLN"   
#> [120] "FEN1"    "FAN1"    "TREX1"   "EXO1"    "APTX"    "ENDOV"   "UBE2A"  
#> [127] "UBE2B"   "RAD18"   "SHPRH"   "HLTF"    "RNF168"  "SPRTN"   "RNF8"   
#> [134] "RNF4"    "UBE2V2"  "UBE2N"   "H2AFX"   "CHAF1A"  "SETMAR"  "BLM"    
#> [141] "WRN"     "RECQL4"  "ATM"     "DCLRE1A" "DCLRE1B" "RPA4"    "PRPF19" 
#> [148] "RECQL"   "RECQL5"  "HELQ"    "RDM1"    "ATR"     "ATRIP"   "MDC1"   
#> [155] "RAD1"    "RAD9A"   "HUS1"    "RAD17"   "CHEK1"   "CHEK2"   "TP53"   
#> [162] "TP53BP1" "RIF1"    "TOPBP1"  "CLK2"    "PER1"   
#> 
#> $APM
#> [1] "B2M"   "HLA-A" "HLA-B" "HLA-C" "TAP1"  "TAP2" 
#> 
#> $Immune_Checkpoint
#> [1] "CD274"    "PDCD1LG2" "CTLA4"    "PDCD1"    "LAG3"     "HAVCR2"   "TIGIT"   
#> 
#> $CellCycle_Reg
#>  [1] "ATM"    "CDKN1A" "CDKN2A" "MDM2"   "TP53"   "CCND1"  "RB1"    "CCNE1" 
#>  [9] "FBXW7"  "E2F3"  
#> 
#> $Pan_F_TBRs
#>  [1] "ACTA2"    "ACTG2"    "ADAM12"   "ADAM19"   "CNN1"     "COL4A1"  
#>  [7] "CTGF"     "CTPS1"    "FAM101B"  "FSTL3"    "HSPB1"    "IGFBP3"  
#> [13] "PXDC1"    "SEMA7A"   "SH3PXD2A" "TAGLN"    "TGFBI"    "TNS1"    
#> [19] "TPM1"    
#> 
```
