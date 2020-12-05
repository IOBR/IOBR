
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IOBR: Immuno-Oncology Biological Research

IOBR is an R package to perform comprehensive analysis of tumor
microenvironment and signatures for immuno-oncology.

### 1.Introduction

  - 1.  IOBR integrates 8 published methodologies decoding tumor
        microenvironment (TME) contexture: `CIBERSORT`, `TIMER`,
        `xCell`, `MCPcounter`, `ESITMATE`, `EPIC`, `IPS`, `quanTIseq`;

  - 2.  IOBR collects 255 published signature gene sets, involving tumor
        microenvironment, tumor metabolism, m6A, exosomes,
        microsatellite instability, and tertiary lymphoid structure.
        Running the function `signature_collection_citation` to attain
        the source papers. The function `signature_collection` returns
        the detail signature genes of all given signatures.

  - 3.  IOBR adopts three computational methods to calculate the
        signature score, comprising `PCA`,`z-score`, and `ssGSEA`;

  - 4.  IOBR integrates multiple approaches for variable transition,
        visualization, batch survival analysis, feature selection, and
        statistical analysis.

  - 5.  IOBR also integrates methods for batch visualization of subgroup
        characteristics.

IOBR package workflow ![IOBR logo](./man/figures/IOBR-Workflow.png)

### 2.Installation

Before installing IOBR, please install all dependencies by executing the
following command in R console:

The dependencies includs `tibble`, `survival`, `survminer`, `limma`,
`limSolve`, `GSVA`, `e1071`, `preprocessCore`, `ggplot2` and `ggpubr`.

``` r
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor","timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen)
}

if (!requireNamespace("EPIC", quietly = TRUE))
  devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
if (!requireNamespace("MCPcounter", quietly = TRUE))
  devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
if (!requireNamespace("estimate", quietly = TRUE)){
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
```

The package is not yet on CRAN or Bioconductor. You can install it from
Github:

``` r
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR",ref="master")
```

Library R packages

``` r
library(IOBR) 
```

### 3.Manual

IOBR pipeline diagram below outlines the data processing flow of this
package, and a detailed guidance of how to use IOBR could be find in the
[HTML vignette](https://iobr.github.io/IOBR/IOBR-VIGNETTE.html), or by
typing the following code to R session.

``` r
browseVignettes("IOBR") 
```

![IOBR logo](./man/figures/IOBR-Package.png)

## 3.Availabie methods to decode TME contexture

``` r
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVM 
#>    "cibersort_abs"              "ips"         "estimate"              "svm" 
#>               lsei              TIMER          quanTIseq 
#>             "lsei"            "timer"        "quantiseq"
# Return available parameter options of TME deconvolution.
```

If you use this package in your work, please cite both our package and
the method(s) you are using.

#### Licenses of the deconvolution methods

| method                                                                      | license                                                                                                       | citation                                                                                                                                                                                                                                                                                                 |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [CIBERSORT](https://cibersort.stanford.edu/)                                | free for non-commerical use only                                                                              | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>                                                            |
| [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) | free ([GPL2.0](https://bioinformatics.mdanderson.org/estimate/))                                              | Vegesna R, Kim H, Torres-Garcia W, …, Verhaak R. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4, 2612. <http://doi.org/10.1038/ncomms3612>                                                                                          |
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)               | free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md))                                 | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., …, Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. <https://doi.org/10.1186/s13073-019-0638-6>           |
| [TIMER](http://cistrome.org/TIMER/)                                         | free ([GPL 2.0](http://cistrome.org/TIMER/download.html))                                                     | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174. <https://doi.org/10.1186/s13059-016-1028-7>                                                          |
| [IPS](https://github.com/icbi-lab/Immunophenogram)                          | free ([BSD](https://github.com/icbi-lab/Immunophenogram/blob/master/LICENSE))                                 | P. Charoentong et al., Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Reports 18, 248-262 (2017). <https://doi.org/10.1016/j.celrep.2016.12.019>                                                                |
| [MCPCounter](https://github.com/ebecht/MCPcounter)                          | free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License))                             | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. <https://doi.org/10.1186/s13059-016-1070-5> |
| [xCell](http://xcell.ucsf.edu/)                                             | free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))                                   | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. <https://doi.org/10.1186/s13059-017-1349-1>                                                                                                                |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)                           | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. <https://doi.org/10.7554/eLife.26476>                                                                  |

## 4.Availabie methods to estimate signatures

``` r
signature_score_calculation_methods
#>           PCA        ssGSEA       z-score   Integration 
#>         "pca"      "ssgsea"      "zscore" "integration"
# Return available parameter options of signature estimation.
```

#### Licenses of the signature-esitmation method

| method                                                                   | license                                                | citation                                                                                                                                                                                                              |
| ------------------------------------------------------------------------ | ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [GSVA](http://www.bioconductor.org/packages/release/bioc/html/GSVA.html) | free ([GPL (\>= 2)](https://github.com/rcastelo/GSVA)) | Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, 7. doi: 10.1186/1471-2105-14-7, <http://www.biomedcentral.com/1471-2105/14/7> |

## 5.Signature collection

``` r
#References of collected signatures
signature_collection_citation[!duplicated(signature_collection_citation$Journal),]
#> # A tibble: 19 x 6
#>    Signatures       `Published year` Journal      Title           PMID  DOI     
#>    <chr>                       <dbl> <chr>        <chr>           <chr> <chr>   
#>  1 CD_8_T_effector              2018 Nature       TGFβ attenuate~ 2944~ 10.1038~
#>  2 TMEscoreA_CIR                2019 Cancer Immu~ Tumor Microenv~ 3084~ 10.1158~
#>  3 CD8_Rooney_et_al             2015 Cell         Molecular and ~ 2559~ 10.1016~
#>  4 T_cell_inflamed~             2017 The Journal~ IFN-γ–related ~ 2865~ 10.1172~
#>  5 MDSC_Wang_et_al              2016 Canccer Dis~ Targeting YAP-~ 2670~ 10.1158~
#>  6 B_cells_Danaher~             2017 Journal for~ Gene expressio~ 2823~ 10.1186~
#>  7 Nature_metaboli~             2019 Nature Meta~ Characterizati~ 3198~ 10.1038~
#>  8 Winter_hypoxia_~             2007 Cancer Rese~ Relation of a ~ 1740~ 10.1158~
#>  9 Hu_hypoxia_sign~             2019 Molecular B~ The Genome Lan~ 3044~ 10.1093~
#> 10 MT_exosome                   2019 Molecular T~ An EV-Associat~ 3147~ 10.1016~
#> 11 SR_exosome                   2017 Scientific ~ Genetic Mutati~ 2838~ 10.1038~
#> 12 MC_Review_Exoso~             2016 Molcular Ca~ Diagnostic, Pr~ 2718~ 10.1186~
#> 13 CMLS_Review_Exo~             2018 Cellular an~ Current knowle~ 2873~ 10.1007~
#> 14 Positive_regula~             2020 Gene Ontolo~ http://geneont~ <NA>  <NA>    
#> 15 Molecular_Cance~             2020 Molecular C~ m6A regulator-~ <NA>  10.1186~
#> 16 Ferroptosis                  2020 IOBR         Constructed by~ <NA>  <NA>    
#> 17 T_cell_accumula~             2018 Nature Medi~ Signatures of ~ 3012~ 10.1038~
#> 18 Antigen_Process~             2020 Nature Comm~ Pan-cancer Cha~ 3208~ 10.1038~
#> 19 CD8_T_cells_Bin~             2013 Immunity     Spatio-tempora~ 2413~ 10.1016~

#signature groups
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

## References

DQ Zeng, ZL Ye, RF Shen, Y Xiong, JN Wu, WJ Qiu, WJ Liao\*.

IOBR: Multi-omics Immuno-Oncology Biological Research to decode tumor
microenvironment and signatures.

## Reporting bugs

Please report bugs to the [Github issues
page](https://github.com/IOBR/IOBR/issues)

E-mail any questions to <dongqiangzeng0808@gmail.com>
