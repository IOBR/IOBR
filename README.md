
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IOBR: Immuno-Oncology Biological Research

[![GitHub
release](https://img.shields.io/github/release/IOBR/IOBR.svg)](https://GitHub.com/IOBR/IOBR/releases/)
[![GitHub
stars](https://img.shields.io/github/stars/IOBR/IOBR.svg)](https://GitHub.com/IOBR/IOBR/stargazers/)
[![GitHub
issues](https://img.shields.io/github/issues/IOBR/IOBR.svg)](https://GitHub.com/IOBR/IOBR/issues/)
[![License](https://img.shields.io/badge/license-GPL%203.0-blue.svg)](https://opensource.org/licenses/GPL-3.0)

IOBR is a comprehensive R package designed for immuno-oncology research,
providing a one-stop solution for tumor microenvironment (TME)
deconvolution, signature analysis, and integrated visualization. It
integrates multiple state-of-the-art algorithms and curated gene sets to
facilitate in-depth analysis of tumor immunity.

## 📊 Package Workflow

<figure>
<img src="./man/figures/IOBR-Workflow.png" alt="IOBR Workflow" />
<figcaption aria-hidden="true">IOBR Workflow</figcaption>
</figure>

<figure>
<img src="./man/figures/IOBR-Package.png" alt="IOBR Package" />
<figcaption aria-hidden="true">IOBR Package</figcaption>
</figure>

## 🎯 Key Features

### 1. Extensive Signature Collection

- **322+ published signature gene sets** covering TME, metabolism, m6A,
  exosomes, microsatellite instability, tertiary lymphoid structure, and
  more
- Easy access to signature genes and source citations
- Flexible signature management and customization

### 2. Multi-algorithm TME Deconvolution

Integrates 8 cutting-edge TME decoding methodologies:

- `CIBERSORT` - Cell-type identification by estimating relative subsets
  of RNA transcripts
- `TIMER` - Tumor Immune Estimation Resource
- `xCell` - Digital portrayal of tissue cellular heterogeneity
- `MCPcounter` - Estimation of immune and stromal cell populations
- `ESTIMATE` - Inference of tumor purity and stromal/immune cell
  admixture
- `EPIC` - Enumeration of cancer and immune cell types
- `IPS` - Immunophenoscore calculation
- `quanTIseq` - Quantification of tumor-infiltrating immune cells

### 3. Signature Score Calculation

Three robust computational methods for signature scoring:

- `PCA` - Principal Component Analysis
- `z-score` - Standardized expression scoring
- `ssGSEA` - Single-sample Gene Set Enrichment Analysis

### 4. Integrated Analysis Pipeline

- Data preprocessing and normalization
- Batch effect correction
- Feature selection and dimensionality reduction
- Survival analysis and risk modeling
- Statistical testing and visualization

### 5. Rich Visualization Tools

- Batch survival analysis plots
- Subgroup characteristic visualization
- Correlation heatmaps and scatter plots
- Forest plots for biomarker validation
- Customizable themes and color palettes

## 📦 Installation

### Prerequisites

- R version 3.6.0 or higher
- Bioconductor version 3.10 or higher

### Install from GitHub

#### Standard Installation

``` r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install IOBR from GitHub
BiocManager::install("IOBR/IOBR")
```

#### For Chinese Users (Faster Download)

``` r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install IOBR using a mirror for faster download
remotes::install_git("https://ghfast.top/https://github.com/IOBR/IOBR")
```

### Load the Package

``` r
library(IOBR)
```

## 🚀 Quick Start

### TME Deconvolution

``` r
# List available TME deconvolution methods
tme_deconvolution_methods

# Perform TME deconvolution using multiple methods
# Assuming you have an expression set object 'eset'
tme_result <- deconvo_tme(eset, 
                          methods = c("cibersort", "timer", "xcell"),
                          output_format = "data.frame")
```

### Signature Score Calculation

``` r
# List available signature score calculation methods
signature_score_calculation_methods

# Calculate signature scores using ssGSEA
# Assuming you have an expression set object 'eset' and signature list 'sig_list'
sig_scores <- sigScore(eset, 
                       signature = sig_list,
                       method = "ssgsea")
```

## 📖 Documentation

### IOBR Book

For detailed tutorials and case studies, please refer to the [IOBR
Book](https://iobr.github.io/book/), which provides comprehensive
guidance on:

- Installation and setup
- Data preprocessing
- TME deconvolution
- Signature analysis
- Survival analysis
- Visualization techniques

### Package Vignettes

Vignettes are available within the package and can be accessed using:

``` r
browseVignettes("IOBR")
```

## 🛠️ Available Methods

### TME Deconvolution Methods

``` r
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVR 
#>    "cibersort_abs"              "ips"         "estimate"              "svr" 
#>               lsei              TIMER          quanTIseq 
#>             "lsei"            "timer"        "quantiseq"
```

### Signature Score Calculation Methods

``` r
signature_score_calculation_methods
#>           PCA        ssGSEA       z-score   Integration 
#>         "pca"      "ssgsea"      "zscore" "integration"
```

### Signature Collection

``` r
data("signature_collection")
# Number of available signatures
length(signature_collection)
#> [1] 323

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

data("signature_collection_citation")
head(signature_collection_citation)
#> # A tibble: 6 × 6
#>   Signatures        `Published year` Journal Title                   PMID  DOI  
#>   <chr>                        <dbl> <chr>   <chr>                   <chr> <chr>
#> 1 CD_8_T_effector               2018 Nature  TGFβ attenuates tumour… 2944… 10.1…
#> 2 DDR                           2018 Nature  TGFβ attenuates tumour… 2944… 10.1…
#> 3 APM                           2018 Nature  TGFβ attenuates tumour… 2944… 10.1…
#> 4 Immune_Checkpoint             2018 Nature  TGFβ attenuates tumour… 2944… 10.1…
#> 5 CellCycle_Reg                 2018 Nature  TGFβ attenuates tumour… 2944… 10.1…
#> 6 Pan_F_TBRs                    2018 Nature  TGFβ attenuates tumour… 2944… 10.1…

data("sig_group")
sig_group[1:3]
#> $tumor_signature
#>  [1] "CellCycle_Reg"                            
#>  [2] "Cell_cycle"                               
#>  [3] "DDR"                                      
#>  [4] "Mismatch_Repair"                          
#>  [5] "Histones"                                 
#>  [6] "Homologous_recombination"                 
#>  [7] "Nature_metabolism_Hypoxia"                
#>  [8] "Molecular_Cancer_m6A"                     
#>  [9] "MT_exosome"                               
#> [10] "Positive_regulation_of_exosomal_secretion"
#> [11] "Ferroptosis"                              
#> [12] "EV_Cell_2020"                             
#> 
#> $EMT
#> [1] "Pan_F_TBRs" "EMT1"       "EMT2"       "EMT3"       "WNT_target"
#> 
#> $io_biomarkers
#>  [1] "TMEscore_CIR"                    "TMEscoreA_CIR"                  
#>  [3] "TMEscoreB_CIR"                   "T_cell_inflamed_GEP_Ayers_et_al"
#>  [5] "CD_8_T_effector"                 "IPS_IPS"                        
#>  [7] "Immune_Checkpoint"               "Exhausted_CD8_Danaher_et_al"    
#>  [9] "Pan_F_TBRs"                      "Mismatch_Repair"                
#> [11] "APM"
```

## 📊 Licenses and Citations

### TME Deconvolution Methods

| Method | License | Citation |
|----|----|----|
| [CIBERSORT](https://cibersort.stanford.edu/) | Free for non-commercial use only | Newman, A. M., et al. (2015). Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337> |
| [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) | Free ([GPL2.0](https://bioinformatics.mdanderson.org/estimate/)) | Vegesna R, et al. (2013). Nature Communications, 4, 2612. <http://doi.org/10.1038/ncomms3612> |
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) | Free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)) | Finotello, F., et al. (2019). Genome Medicine, 11(1), 34. <https://doi.org/10.1186/s13073-019-0638-6> |
| [TIMER](http://cistrome.org/TIMER/) | Free ([GPL 2.0](http://cistrome.org/TIMER/download.html)) | Li, B., et al. (2016). Genome Biology, 17(1), 174. <https://doi.org/10.1186/s13059-016-1028-7> |
| [IPS](https://github.com/icbi-lab/Immunophenogram) | Free ([BSD](https://github.com/icbi-lab/Immunophenogram/blob/master/LICENSE)) | Charoentong P, et al. (2017). Cell Reports, 18, 248-262. <https://doi.org/10.1016/j.celrep.2016.12.019> |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | Free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License)) | Becht, E., et al. (2016). Genome Biology, 17(1), 218. <https://doi.org/10.1186/s13059-016-1070-5> |
| [xCell](http://xcell.ucsf.edu/) | Free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION)) | Aran, D., et al. (2017). Genome Biology, 18(1), 220. <https://doi.org/10.1186/s13059-017-1349-1> |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | Free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., et al. (2017). eLife, 6, e26476. <https://doi.org/10.7554/eLife.26476> |

### Signature Estimation Methods

| Method | License | Citation |
|----|----|----|
| [GSVA](http://www.bioconductor.org/packages/release/bioc/html/GSVA.html) | Free ([GPL (\>= 2)](https://github.com/rcastelo/GSVA)) | Hänzelmann S, et al. (2013). BMC Bioinformatics, 14, 7. <https://doi.org/10.1186/1471-2105-14-7> |

## 📝 Citing IOBR

If you use IOBR in your research, please cite both the IOBR package and
the specific methods you employ.

1.  Zeng DQ, Fang YR, …, Liao WJ. Enhancing Immuno-Oncology
    Investigations Through Multidimensional Decoding of Tumour
    Microenvironment with IOBR 2.0, **Cell Reports Methods**, 2024
    <https://doi.org/10.1016/j.crmeth.2024.100910>

2.  Fang YR, …, Liao WJ, Zeng DQ, Systematic Investigation of Tumor
    Microenvironment and Antitumor Immunity With IOBR, **Med Research**,
    2025 <https://onlinelibrary.wiley.com/doi/epdf/10.1002/mdr2.70001>

## 🤝 Contributing

We welcome contributions to IOBR! If you’re interested in contributing,
please:

1.  Fork the repository
2.  Create a feature branch
3.  Make your changes
4.  Submit a pull request

Please ensure your code follows the project’s coding standards and
includes appropriate documentation and tests.

## 🐛 Reporting Bugs

If you encounter any bugs or issues, please report them to the [GitHub
issues page](https://github.com/IOBR/IOBR/issues). When reporting bugs,
please include:

- A clear and descriptive title
- A detailed description of the issue
- Minimal reproducible example (if possible)
- Your R and IOBR versions
- Any error messages or warnings

## 📧 Contact

For questions or inquiries, please contact:

- Dr. Deqiang Zeng: <interlaken@smu.edu.cn>
- Dr. Yunru Fang: <fyr_nate@163.com>

## 📄 License

IOBR is released under the [GNU General Public License v3.0](LICENSE).
