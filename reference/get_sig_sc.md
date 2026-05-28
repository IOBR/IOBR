# Extract Top Marker Genes from Single-Cell Differential Results

Selects the top N marker genes per cluster from a ranked differential
expression result table.

## Usage

``` r
get_sig_sc(
  deg,
  cluster = "cluster",
  gene = "gene",
  avg_log2FC = "avg_log2FC",
  n = 100
)
```

## Arguments

- deg:

  Data frame or matrix. Ranked marker statistics.

- cluster:

  Character. Column name containing cluster identifiers. Default is
  \`"cluster"\`.

- gene:

  Character. Column name containing gene identifiers. Default is
  \`"gene"\`.

- avg_log2FC:

  Character. Column name for average log2 fold change. Default is
  \`"avg_log2FC"\`.

- n:

  Integer. Number of top markers per cluster. Default is 100.

## Value

List of character vectors; each element contains the top N genes for a
cluster.

## Examples

``` r
# \donttest{
deg <- load_data("deg")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "deg"
get_sig_sc(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100)
#> $`Epithelial cells 11`
#>   [1] "AQP4"        "BCAT1"       "RGS16"       "LRRN4"       "TIMP3"      
#>   [6] "SUSD2"       "AQP1"        "CYP4B1"      "AGER"        "FAM46B"     
#>  [11] "SCGB3A1"     "PEG10"       "AK1"         "DLX3"        "GSTM3"      
#>  [16] "BTG2"        "GDF15"       "CTSE"        "CLDN18"      "EFNA1"      
#>  [21] "DEGS2"       "CYR61"       "CKB"         "KLHL35"      "ARL4D"      
#>  [26] "C4BPA"       "SRRM2"       "HMGN2"       "SFTPB"       "RNASE1"     
#>  [31] "EGR1"        "C16orf89"    "IFIT1"       "FBLN5"       "SFTA2"      
#>  [36] "IRX5"        "GGH"         "FABP3"       "SULT1A2"     "SLC34A2"    
#>  [41] "CDKN2A"      "HMGN3"       "RSRP1"       "C12orf65"    "SFTPD"      
#>  [46] "ARHGAP24"    "TMEM37"      "NAPSA"       "IGFBP2"      "TUBA4A"     
#>  [51] "CA2"         "PIFO"        "MLF1"        "THUMPD3-AS1" "CLU"        
#>  [56] "TRA2A"       "IFI27"       "ATHL1"       "ARHGEF2"     "CYP51A1"    
#>  [61] "STMN1"       "PGC"         "SFPQ"        "CDKN1C"      "ACBD3"      
#>  [66] "DGKD"        "C19orf53"    "FAM177A1"    "HIST1H2AM"   "SCGB3A2"    
#>  [71] "ZNF593"      "ARHGEF17"    "ATP13A4-AS1" "ISG15"       "MMP15"      
#>  [76] "GADD45GIP1"  "PEBP4"       "ZFAND2A"     "PARP1"       "TPPP3"      
#>  [81] "GADD45G"     "CD55"        "SERTAD1"     "RDH10"       "OSR2"       
#>  [86] "CLIC3"       "TPPP"        "SNRNP25"     "C4BPB"       "RANBP10"    
#>  [91] "TNNC1"       "CRYM"        "PHLDA2"      "TFAP2C"      "EPHX1"      
#>  [96] "SFN"         "MRPL2"       "ZNF493"      "SNX22"       "ALKBH4"     
#> 
#> $`Epithelial cells 15`
#>   [1] "SFTPC"       "SFTPA1"      "SFTPA2"      "FABP5"       "SFTPD"      
#>   [6] "SFTPB"       "PEBP4"       "ABCA3"       "CSF3R"       "NAPSA"      
#>  [11] "PGC"         "WIF1"        "MACROD2"     "CLDN18"      "SLPI"       
#>  [16] "HP"          "DUOX1"       "C16orf89"    "AQP1"        "AK1"        
#>  [21] "SCGB3A1"     "CA2"         "SCGB3A2"     "SCD"         "AFF3"       
#>  [26] "NECAB1"      "FGG"         "SLC34A2"     "MFSD2A"      "GADD45G"    
#>  [31] "SFTA2"       "PTPN13"      "RASGRF1"     "AQP4"        "C1R"        
#>  [36] "ACOXL"       "NR0B2"       "PLA2G4F"     "DCXR"        "SUSD2"      
#>  [41] "ADIRF"       "FILIP1"      "RGS16"       "MYLK"        "HLA-DQB2"   
#>  [46] "ORM1"        "ETV5"        "TSPAN7"      "CXCL17"      "PHYHD1"     
#>  [51] "RNASE1"      "TTN"         "SNX30"       "C4BPA"       "CHI3L1"     
#>  [56] "IFITM1"      "FAM46C"      "SFTA1P"      "CHI3L2"      "GUCY1A3"    
#>  [61] "CLIC3"       "GOLIM4"      "EGR1"        "TPPP3"       "SULT1A2"    
#>  [66] "CYP4B1"      "ATP13A4-AS1" "KCNJ2"       "SLC5A2"      "TMEM37"     
#>  [71] "PLLP"        "PTTG1"       "ARL4D"       "PLIN5"       "HMGN2"      
#>  [76] "PMM1"        "RSRP1"       "CLU"         "PIGR"        "QDPR"       
#>  [81] "AGER"        "DMBT1"       "PEG10"       "RSBN1"       "CKB"        
#>  [86] "RMDN2"       "CYR61"       "NUDT14"      "F3"          "CTSE"       
#>  [91] "C3"          "TUBA4A"      "IRX5"        "C12orf65"    "CCND3"      
#>  [96] "TMEM116"     "TUBA1B"      "PLEKHB1"     "IVD"         "ACP5"       
#> 
#> $`Epithelial cells 2`
#>   [1] "IGFBP3"   "PCDH7"    "PTPRN2"   "ADAM28"   "CASP1"    "TRAF2"   
#>   [7] "CYBA"     "SERPINE1" "IL32"     "RFK"      "TRIM16"   "FXYD5"   
#>  [13] "RAMP1"    "EGFL7"    "SNAP23"   "ATP6V1E1" "HSD11B2"  "BATF"    
#>  [19] "UACA"     "TCTN3"    "HSPA9"    "BIK"      "ACSL1"    "POLR3K"  
#>  [25] "CCL28"    "AUH"      "CMC4"     "PRPS1"    "EFHD2"    "RAB24"   
#>  [31] "RER1"     "HN1"      "MTCH2"    "NTMT1"    "COMTD1"   "KRT16"   
#>  [37] "STX7"     "PDCD6IP"  "UBTD2"    "CHMP1B"   "DSC2"     "UBA7"    
#>  [43] "VTA1"     "P4HA1"    "UFM1"     "PTP4A1"   "CEBPG"    "AMZ2"    
#>  [49] "PLIN2"    "PRPS2"    "RNASEH1"  "KATNBL1"  "PSMD5"    "MECR"    
#>  [55] "DCUN1D1"  "PEX19"    "TIMM17A"  "SGMS1"    "C7orf73"  "ORC5"    
#>  [61] "UBQLN1"   "VAPA"     "PAFAH1B3" "PEX11B"   "DRG1"     "ADSL"    
#>  [67] "UBN1"     "CLNS1A"   "SUMF1"    "NSMCE1"   "COG7"     "DGAT1"   
#>  [73] "NDUFS3"   "YIPF5"    "CCDC78"   "RAB6A"    "LPP"      "CHMP2B"  
#>  [79] "G3BP1"    "AAAS"     "MOV10"    "C3"       "CRIP1"    "LSG1"    
#>  [85] "RSPH3"    "SNRPA1"   "DNAJA3"   "HIGD1A"   "ARHGDIB"  "PIGG"    
#>  [91] "CCT4"     "TES"      "CCNJL"    "LEPROT"   "ATP6V1D"  "AFMID"   
#>  [97] "TALDO1"   "NGRN"     "PSG5"     "CD55"    
#> 
#> $`Epithelial cells 23`
#>   [1] "NECAB1"      "SFTPC"       "SFTPA1"      "PGC"         "SFTPA2"     
#>   [6] "CA2"         "PEBP4"       "ABCA3"       "ETV5"        "HNRNPH1"    
#>  [11] "CHI3L2"      "SAR1A"       "MFSD2A"      "C4BPA"       "C1R"        
#>  [16] "LTF"         "F3"          "CXCL2"       "NR0B2"       "PIGR"       
#>  [21] "FABP5"       "SFTPD"       "NAPSA"       "EGR1"        "PTPN13"     
#>  [26] "SFTPB"       "CLU"         "RRAD"        "WIF1"        "DUOX1"      
#>  [31] "AQP1"        "NFKBIZ"      "TOB1"        "SLC34A2"     "TRA2A"      
#>  [36] "RASGRF1"     "UGCG"        "HACD1"       "PHYHD1"      "AFF3"       
#>  [41] "CHI3L1"      "GOLIM4"      "CD274"       "KCNJ2"       "SLC38A2"    
#>  [46] "PPP3CA"      "C16orf89"    "ACOXL"       "SCD"         "SKIL"       
#>  [51] "RDH10"       "WTAP"        "TSC22D2"     "SCGB3A1"     "ARGLU1"     
#>  [56] "DMBT1"       "AZGP1"       "C8orf4"      "NNMT"        "RNASE1"     
#>  [61] "PLA2G4F"     "NUCKS1"      "FGG"         "LRRC75A-AS1" "UBE2B"      
#>  [66] "SERTAD1"     "HMGCS1"      "MBNL1"       "PLIN5"       "ADAM17"     
#>  [71] "DDX3X"       "CTNNB1"      "FAM46B"      "PPP1CB"      "CLDN18"     
#>  [76] "ATP13A4-AS1" "SNX30"       "CYP51A1"     "PTP4A1"      "SUSD2"      
#>  [81] "TFAP2C"      "HP"          "STOM"        "STAM"        "RAB21"      
#>  [86] "MPZL3"       "CRIM1"       "MACROD2"     "STK17A"      "UGDH"       
#>  [91] "KCTD9"       "ALDH6A1"     "FNIP1"       "ZC3H4"       "DUSP14"     
#>  [96] "CXCL17"      "DCBLD2"      "PPP6R3"      "ARL8B"       "SLMO2"      
#> 
#> $`Epithelial cells 25`
#>   [1] "LYPD2"       "REG1A"       "SERPINB3"    "SERPINB5"    "CLDN10"     
#>   [6] "ARL14"       "SERPINB4"    "MUC5AC"      "AKR1C2"      "KRT17"      
#>  [11] "IL1RN"       "ALDH3A1"     "C12orf75"    "SERPINA3"    "KRT6A"      
#>  [16] "BASP1"       "BPIFB1"      "PRKCDBP"     "CDKN2A"      "HNRNPH1"    
#>  [21] "SFN"         "MUC5B"       "CXCL1"       "RARRES1"     "F3"         
#>  [26] "IFITM1"      "NUCKS1"      "MLF1"        "ST6GAL1"     "EREG"       
#>  [31] "OXCT1"       "MAFIP"       "IFI27"       "KLK11"       "PHLDA2"     
#>  [36] "TIMP1"       "TOB1"        "LRRC75A-AS1" "S100A2"      "EEF1B2"     
#>  [41] "PADI2"       "AEN"         "BCL10"       "EGR1"        "HMGCS1"     
#>  [46] "WTAP"        "TRA2A"       "HACD1"       "BTG2"        "AKR1C1"     
#>  [51] "SKIL"        "OLFM1"       "MPRIP"       "SNRPB"       "WEE1"       
#>  [56] "JAG1"        "SLC7A5P2"    "ATHL1"       "EFHD2"       "GMDS"       
#>  [61] "C12orf57"    "IL17RB"      "GGH"         "UGCG"        "MIR205HG"   
#>  [66] "STAT2"       "SERTAD1"     "R3HDM2"      "CTNNB1"      "NCOR1"      
#>  [71] "ARGLU1"      "PDCD4"       "CLTB"        "GMNN"        "RIN2"       
#>  [76] "PPP3CA"      "UBE2B"       "SLC5A2"      "PPP1CB"      "STON2"      
#>  [81] "DST"         "KIAA1217"    "ARFIP2"      "FAM208B"     "RASSF6"     
#>  [86] "HIST1H2AM"   "SMURF2"      "FKBP9"       "ARL8B"       "RREB1"      
#>  [91] "UGDH"        "HIST1H4C"    "HMGN3"       "BIK"         "B3GALNT2"   
#>  [96] "C4BPB"       "TRIM16"      "PDGFA"       "ZC3H4"       "HMGCR"      
#> 
#> $`Epithelial cells 26`
#>   [1] "ITLN2"    "FMO2"     "RTKN2"    "AGER"     "HEG1"     "CRYAB"   
#>   [7] "IGFBP7"   "CLDN18"   "SPOCK2"   "UPK3B"    "CYP4B1"   "TNNC1"   
#>  [13] "PDPN"     "ADIRF"    "SFTA1P"   "TIMP3"    "KLK11"    "COL4A2"  
#>  [19] "AQP4"     "ABCA1"    "PLLP"     "CLIC3"    "RGCC"     "CCND2"   
#>  [25] "TAGLN"    "SUSD2"    "SBSPON"   "RNASE1"   "WFS1"     "MGLL"    
#>  [31] "MAP2"     "SULT1A2"  "TUBA1A"   "UNC13D"   "IFI27"    "CD55"    
#>  [37] "FKBP1B"   "AKR1C1"   "RIN2"     "F3"       "FILIP1"   "FN1"     
#>  [43] "PEBP4"    "FAM46B"   "SFTA2"    "PRKCDBP"  "CRIP1"    "TMEM98"  
#>  [49] "CEACAM6"  "CKB"      "FBLN5"    "PHACTR2"  "EPB41L5"  "DST"     
#>  [55] "ARHGAP24" "PEG10"    "ABI3BP"   "LRRN4"    "TSPAN7"   "CTGF"    
#>  [61] "FBXL15"   "STOM"     "HMGN2"    "SDC1"     "STON2"    "QKI"     
#>  [67] "RFC1"     "PDGFA"    "FAS"      "IFT43"    "TXNRD2"   "DUOX1"   
#>  [73] "ARL8B"    "NHLRC3"   "TERF1"    "SNX22"    "CYR61"    "RBM17"   
#>  [79] "TJP1"     "SAP30BP"  "SGMS1"    "ARAP2"    "HSD17B8"  "MICA"    
#>  [85] "STX7"     "TRIM5"    "VPS4B"    "PDS5B"    "NUCKS1"   "PCMTD1"  
#>  [91] "SFTPB"    "MDM4"     "FAM134A"  "ACAA2"    "CTSE"     "ABCA7"   
#>  [97] "SNRPB"    "CIAPIN1"  "HNRNPA0"  "PHF10"   
#> 
#> $`Epithelial cells 27`
#>   [1] "APOBEC3H"         "ALDH3A1"          "LDHD"            
#>   [4] "AZGP1"            "SCGB3A2"          "MT1G"            
#>   [7] "KLK11"            "TPPP3"            "TMEM37"          
#>  [10] "RASL11A"          "CTSE"             "CEACAM5"         
#>  [13] "DEGS2"            "ITPR2"            "NHS"             
#>  [16] "LL22NC03-75H12.2" "TMEM98"           "C16orf89"        
#>  [19] "C4BPA"            "HLA-DQB2"         "SLC1A7"          
#>  [22] "GSTA1"            "AKR7A2"           "BTG2"            
#>  [25] "CXCL17"           "RARRES2"          "IGFBP2"          
#>  [28] "TMEM230"          "SUSD2"            "CYP4B1"          
#>  [31] "PPIE"             "SYT2"             "LINC00578"       
#>  [34] "TERF2IP"          "SRRM2"            "SFTA2"           
#>  [37] "PMM1"             "DHRS4-AS1"        "GDF15"           
#>  [40] "SFTPB"            "NME3"             "CXCL14"          
#>  [43] "PLLP"             "CLU"              "ELOVL1"          
#>  [46] "SFTPD"            "PEMT"             "MLF1"            
#>  [49] "IFI27"            "MRPL57"           "NGDN"            
#>  [52] "SNRPB"            "SNRNP25"          "MTSS1L"          
#>  [55] "UBE2I"            "MAFIP"            "AQP4"            
#>  [58] "SCGB3A1"          "CHMP4A"           "PIGR"            
#>  [61] "ITPA"             "IRX5"             "CEACAM6"         
#>  [64] "SZT2"             "TMEM205"          "RGS16"           
#>  [67] "DNAJA3"           "NAPSA"            "TMEM52"          
#>  [70] "POLR3K"           "GFRA3"            "MSRB1"           
#>  [73] "CCDC101"          "IVD"              "PCNA"            
#>  [76] "NUDT16L1"         "HLA-DOB"          "BCL6"            
#>  [79] "FTO"              "SFTA1P"           "C10orf32"        
#>  [82] "MLYCD"            "AK1"              "CCDC115"         
#>  [85] "DMBT1"            "RNASE1"           "CIAPIN1"         
#>  [88] "FAM71E1"          "COQ9"             "CPAMD8"          
#>  [91] "ALDH6A1"          "MMP15"            "PLEKHB1"         
#>  [94] "RMDN2"            "LYRM1"            "HSD17B8"         
#>  [97] "YRDC"             "ERAL1"            "SLC5A2"          
#> [100] "MMP7"            
#> 
#> $`Epithelial cells 28`
#>   [1] "UPK1B"      "BPIFB1"     "KLK6"       "PSG3"       "PPBP"      
#>   [6] "LDLRAD1"    "PADI2"      "MSMB"       "ADGRE2"     "PSG5"      
#>  [11] "IGFBP3"     "SLC30A2"    "APOD"       "TALDO1"     "KRT6A"     
#>  [16] "KRT16"      "TMEM140"    "IL3RA"      "KRT4"       "CHMP1B"    
#>  [21] "ZNF793-AS1" "TRIM16"     "P4HA1"      "APOL2"      "WFS1"      
#>  [26] "KRT17"      "MARCKS"     "KRT13"      "MOV10"      "HDAC9"     
#>  [31] "CLN5"       "RNF13"      "CLP1"       "DAZAP2"     "WEE1"      
#>  [36] "BIK"        "UACA"       "RNF114"     "DYNLL1"     "S100A9"    
#>  [41] "NDRG3"      "TMEM106C"   "QRICH1"     "SUMF1"      "CSTF1"     
#>  [46] "PPAP2B"     "P4HTM"      "LOC283788"  "CHID1"      "UBE2D1"    
#>  [51] "ARFIP2"     "RARRES1"    "UNG"        "TBC1D1"     "SERPINE1"  
#>  [56] "PEX11B"     "TMTC3"      "TOMM34"     "NARS"       "STOM"      
#>  [61] "UBQLN1"     "PIGG"       "COMMD7"     "RFK"        "BLOC1S4"   
#>  [66] "ADAM10"     "TFAP2C"     "CTNNB1"     "FUBP1"      "HAT1"      
#>  [71] "CRIP1"      "C3"         "HNRNPC"     "VTA1"       "CCT4"      
#>  [76] "PDCD6IP"    "LSG1"       "SNAP23"     "HIGD1A"     "PAFAH1B3"  
#>  [81] "AAAS"       "PIP4K2C"    "CLIC3"      "RER1"       "UBA7"      
#>  [86] "RSPH3"      "UFM1"       "AUH"        "CEACAM5"    "VPS41"     
#>  [91] "TCTN3"      "TJP1"       "IFIT2"      "CCL28"      "CHMP2B"    
#>  [96] "COG7"       "ACTR10"     "FXYD5"      "PTP4A1"     "HSPA9"     
#> 
#> $`Epithelial cells 29`
#>   [1] "PRSS2"      "PRSS1"      "PRSS3"      "ATG9B"      "PLA1A"     
#>   [6] "PAEP"       "AZGP1"      "UBE2C"      "TNFRSF18"   "G0S2"      
#>  [11] "CDA"        "ALDOC"      "NME1-NME2"  "MIR205HG"   "CDKN2A"    
#>  [16] "AKR1B1"     "CHI3L1"     "C4orf48"    "YEATS4"     "MAP1B"     
#>  [21] "TNNC2"      "RGS17"      "PPAT"       "MSRB1"      "ZNF593"    
#>  [26] "LAGE3"      "RAN"        "ICT1"       "NOP16"      "RASSF3"    
#>  [31] "EEF1B2"     "SLC27A5"    "MYEOV2"     "GEMIN6"     "TIMP1"     
#>  [36] "CEACAM5"    "CD320"      "HN1"        "MRPL13"     "APOD"      
#>  [41] "S100A2"     "MRPL32"     "MRPL21"     "PAFAH1B3"   "GDF15"     
#>  [46] "GADD45GIP1" "PHLDA2"     "QTRTD1"     "ANGPTL4"    "C19orf24"  
#>  [51] "POLE4"      "SUPT4H1"    "MTIF2"      "MRPL57"     "DDX18"     
#>  [56] "ACOX2"      "SDC1"       "CCDC59"     "DST"        "ATP5D"     
#>  [61] "SNHG17"     "SLC34A2"    "DCBLD2"     "SNRNP25"    "NIFK"      
#>  [66] "MPHOSPH10"  "NAA10"      "TBCB"       "S100A9"     "PLIN2"     
#>  [71] "MRPL40"     "ITPA"       "TFAM"       "TMEM61"     "TMEM98"    
#>  [76] "GPATCH4"    "TRMT10C"    "SLMO2"      "RPS19BP1"   "RPP21"     
#>  [81] "NUFIP2"     "NAPSA"      "RAMP1"      "MEA1"       "NSMCE1"    
#>  [86] "CEBPZOS"    "APOL2"      "TMEM14A"    "COMTD1"     "PSMD9"     
#>  [91] "DTYMK"      "SNHG11"     "HINT2"      "NT5C"       "TMEM91"    
#>  [96] "SFTA2"      "AP1AR"      "POLR3K"     "GOLT1A"     "SMIM6"     
#> 
#> $`Epithelial cells 6`
#>   [1] "PTN"        "CHGB"       "NR2F1-AS1"  "CRLF1"      "PGF"       
#>   [6] "MUC5B"      "MGP"        "COL1A1"     "CCDC80"     "KRT81"     
#>  [11] "MIR205HG"   "CEACAM6"    "IFITM1"     "NHS"        "EREG"      
#>  [16] "RGS17"      "CDKN2A"     "CCL2"       "CAMK1D"     "S100A9"    
#>  [21] "IFI27"      "PHGDH"      "SNHG18"     "DMBT1"      "CEMIP"     
#>  [26] "NNMT"       "ARL4A"      "CXCL14"     "OXCT1"      "CEACAM5"   
#>  [31] "TNFRSF18"   "TGM2"       "THBS1"      "CYSLTR1"    "COL17A1"   
#>  [36] "HDGFRP3"    "LOC648987"  "EMB"        "TNFSF10"    "SPON2"     
#>  [41] "ARHGDIB"    "PIGR"       "NR1D1"      "CLEC2D"     "RN7SK"     
#>  [46] "RNASE1"     "TIMP1"      "SLC38A1"    "CXCL1"      "CDKN1C"    
#>  [51] "ATHL1"      "TGFBI"      "MARCKS"     "SLC12A2"    "GUCY1A3"   
#>  [56] "PHLDA2"     "TRA2A"      "EXOC3-AS1"  "LRRC37A4P"  "C4orf48"   
#>  [61] "C21orf2"    "WEE1"       "PIP4K2C"    "CD320"      "C12orf57"  
#>  [66] "IGFBP2"     "C8orf4"     "TMEM205"    "XRRA1"      "IRF2BPL"   
#>  [71] "STON2"      "RRAD"       "SNHG1"      "EPHX1"      "ISG15"     
#>  [76] "TMEM238"    "MRPL39"     "IRF7"       "CHI3L1"     "NFKBIZ"    
#>  [81] "CCDC85B"    "YBEY"       "LINC00649"  "C4BPB"      "ZNF195"    
#>  [86] "KRT17"      "ATP5D"      "PDK3"       "CAPS"       "PMAIP1"    
#>  [91] "IVD"        "DST"        "HSPA6"      "SAC3D1"     "UBE2S"     
#>  [96] "ABCA7"      "SEPT11"     "C19orf60"   "HPGD"       "GADD45GIP1"
#> 
# }
```
