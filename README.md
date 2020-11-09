
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IOBR: Immuno-Oncology Bioinformatics Research

IOBR is an R package to perform comprehensive analysis of tumor
microenvironment and signatures for immuno-oncology.

## 1 Introduction and Installation

### 1.1 Introduction

  - 1.IOBR integrates 8 published methodologies decoding tumor
    microenvironment (TME) contexture: `CIBERSORT`, `TIMER`, `xCell`,
    `MCPcounter`, `ESITMATE`, `EPIC`, `IPS`, `quanTIseq`;
  - 2.IOBR collects 255 published signature gene sets, involving tumor
    microenvironment, tumor metabolism, m6A, exosomes, microsatellite
    instability, and tertiary lymphoid structure. Running the function
    `signature_collection_citation` to attain the source papers. The
    function `signature_collection` returns the detail signature genes
    of all given signatures.
  - 3.IOBR adopts three computational methods to calculate the signature
    score, comprising `PCA`,`z-score`, and `ssGSEA`;
  - 4.IOBR integrates multiple approaches for variable transition,
    visulization, batch survival analysis, and statistical analysis.
  - 5.IOBR also integrates methods for batch visualization of subgroup
    characteristics.

![IOBR logo](./man/figures/IOBR-Package.png)

### 1.2 Installation

Before installing IOBR, please install all dependencies by executing the
following command in R console:

The dependencies includs `tibble`, `survival`, `survminer`, `limma`,
`limSolve`, `GSVA`, `e1071`, `preprocessCore`, `ggplot2` and `ggpubr`.

``` r
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', 
          'ggpubr',"devtools")
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

The package is not yet on CRAN. You can install it from Github:

``` r
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("DongqiangZeng0808/IOBR",ref="master")
#> Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading
#> 'IOBR'
#> Warning: replacing previous import 'caret::cluster' by 'survival::cluster' when
#> loading 'IOBR'
```

## 1.3 Library R packages

``` r
library(IOBR) 
library(EPIC)
library(estimate) 
library(MCPcounter)
Sys.setenv(LANGUAGE = "en") # The error message would be present in English.
options(stringsAsFactors = FALSE) #To avoid transiting character into factor.
```

## 2 TME deconvolution

### 2.1 Availabie methods to decode TME contexture

``` r
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVM 
#>    "cibersort_abs"              "ips"         "estimate"          "svm_ref" 
#>               lsei              TIMER          quanTIseq 
#>         "lsei_ref"            "timer"        "quantiseq"
# Return available parameter options of deconvolution methods
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

#### Licenses of the signature-esitmation method

| method                                                                   | license                                                | citation                                                                                                                                                                                                              |
| ------------------------------------------------------------------------ | ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [GSVA](http://www.bioconductor.org/packages/release/bioc/html/GSVA.html) | free ([GPL (\>= 2)](https://github.com/rcastelo/GSVA)) | Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, 7. doi: 10.1186/1471-2105-14-7, <http://www.biomedcentral.com/1471-2105/14/7> |

The input data is a matrix (log2(TMP+1) transformed) containing 98
TCGA-STAD samples, with genes in rows and samples in columns. The
rownames must be HGNC symbols and the colnames must be sample names.

``` r
# Check data 
eset_stad[1:5,1:5]
#>         TCGA-B7-5818 TCGA-BR-4187 TCGA-BR-4201 TCGA-BR-4253 TCGA-BR-4256
#> MT-CO1      15.18012     15.55806     14.60960     14.63728     15.23528
#> MT-CO3      14.75536     15.19199     14.55337     13.54925     14.30425
#> MT-ND4      14.19637     15.61564     15.80262     14.98329     14.83764
#> MT-CO2      15.10790     15.50514     15.43261     14.52009     14.65806
#> MT-RNR2     14.22690     14.71157     13.48096     13.32553     13.55689
```

Check detail parameters of the function

``` r
help(deconvo_tme)
#> starting httpd help server ... done
```

#### 2.1.1 Method 1: Use CIBERSORT to decode the TME contexture.

``` r
cibersort<-deconvo_tme(eset = eset_stad,method = "cibersort",arrays = FALSE,perm = 200 )
#> 
#> >>> Running CIBERSORT
head(cibersort)
#>             ID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT
#> 1 TCGA-B7-5818              0.03229758                        0
#> 2 TCGA-BR-4187              0.08664611                        0
#> 3 TCGA-BR-4201              0.04743638                        0
#> 4 TCGA-BR-4253              0.01247806                        0
#> 5 TCGA-BR-4256              0.05441618                        0
#> 6 TCGA-BR-4257              0.02463128                        0
#>   Plasma_cells_CIBERSORT T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1            0.000000000            0.19281089                           0
#> 2            0.000000000            0.08724967                           0
#> 3            0.006442484            0.02859921                           0
#> 4            0.002572915            0.22438042                           0
#> 5            0.009230031            0.09355182                           0
#> 6            0.001619011            0.12445594                           0
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                            0.1049710                             0.07551750
#> 2                            0.1430036                             0.07453247
#> 3                            0.2153594                             0.05864149
#> 4                            0.0000000                             0.30129769
#> 5                            0.1756031                             0.07596749
#> 6                            0.1104545                             0.09309550
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                          0.01745022                          0.095616200
#> 2                          0.00000000                          0.000000000
#> 3                          0.00000000                          0.000000000
#> 4                          0.01190655                          0.003320162
#> 5                          0.03274849                          0.005481215
#> 6                          0.04139932                          0.017168846
#>   T_cells_gamma_delta_CIBERSORT NK_cells_resting_CIBERSORT
#> 1                   0.000000000                 0.05889901
#> 2                   0.008804990                 0.02322682
#> 3                   0.006425457                 0.04398323
#> 4                   0.033623962                 0.06326207
#> 5                   0.000000000                 0.04669955
#> 6                   0.000000000                 0.04331891
#>   NK_cells_activated_CIBERSORT Monocytes_CIBERSORT Macrophages_M0_CIBERSORT
#> 1                            0          0.01116026               0.15037855
#> 2                            0          0.01107405               0.00000000
#> 3                            0          0.00000000               0.11864698
#> 4                            0          0.00000000               0.04664432
#> 5                            0          0.07281264               0.06173783
#> 6                            0          0.00000000               0.20723276
#>   Macrophages_M1_CIBERSORT Macrophages_M2_CIBERSORT
#> 1               0.09342079                0.1365691
#> 2               0.02774632                0.3037347
#> 3               0.06025746                0.1987880
#> 4               0.11419626                0.1433415
#> 5               0.04691225                0.1280897
#> 6               0.04149607                0.1869276
#>   Dendritic_cells_resting_CIBERSORT Dendritic_cells_activated_CIBERSORT
#> 1                       0.013276957                          0.00000000
#> 2                       0.003343398                          0.00000000
#> 3                       0.003612247                          0.00000000
#> 4                       0.004120229                          0.01071592
#> 5                       0.000000000                          0.02940856
#> 6                       0.023247600                          0.00000000
#>   Mast_cells_resting_CIBERSORT Mast_cells_activated_CIBERSORT
#> 1                   0.00000000                   1.380008e-02
#> 2                   0.18609633                   3.292483e-05
#> 3                   0.00000000                   9.116860e-02
#> 4                   0.01744220                   0.000000e+00
#> 5                   0.07624986                   2.927279e-02
#> 6                   0.00000000                   3.515756e-02
#>   Eosinophils_CIBERSORT Neutrophils_CIBERSORT P-value_CIBERSORT
#> 1           0.000000000            0.00383182             0.000
#> 2           0.005580183            0.03892844             0.000
#> 3           0.033511519            0.08712756             0.070
#> 4           0.000000000            0.01069773             0.000
#> 5           0.000000000            0.06181853             0.045
#> 6           0.000000000            0.04979512             0.000
#>   Correlation_CIBERSORT RMSE_CIBERSORT
#> 1             0.2538008      0.9851012
#> 2             0.4387743      0.8997892
#> 3             0.1138376      1.0366015
#> 4             0.2838719      0.9730242
#> 5             0.1365922      1.0197311
#> 6             0.2417254      0.9874249
res<-cell_bar_plot(input = cibersort[1:12,], title = "CIBERSORT Cell Fraction")
#> There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap
```

<img src="man/figuresunnamed-chunk-8-1.png" width="100%" />

#### 2.1.2 Method 2: Use EPIC to decode the TME contexture.

``` r
help(deconvo_epic)
epic<-deconvo_tme(eset = eset_stad,method = "epic",arrays = FALSE)
#> 
#> >>> Running EPIC
#> Warning in EPIC::EPIC(bulk = eset, reference = ref, mRNA_cell = NULL, scaleExprs = TRUE): The optimization didn't fully converge for some samples:
#> TCGA-BR-6705; TCGA-BR-7958; TCGA-BR-8366; TCGA-BR-A4J4; TCGA-BR-A4J9; TCGA-CD-A4MG; TCGA-FP-8209; TCGA-HU-A4G9
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning in EPIC::EPIC(bulk = eset, reference = ref, mRNA_cell = NULL, scaleExprs
#> = TRUE): mRNA_cell value unknown for some cell types: CAFs, Endothelial - using
#> the default value of 0.4 for these but this might bias the true cell proportions
#> from all cell types.
head(epic)
#>             ID Bcells_EPIC  CAFs_EPIC CD4_Tcells_EPIC CD8_Tcells_EPIC
#> 1 TCGA-B7-5818  0.02260847 0.02056021       0.1631332      0.10990559
#> 2 TCGA-BR-4187  0.05108222 0.02541769       0.1892978      0.10720878
#> 3 TCGA-BR-4201  0.04348886 0.02447877       0.2423876      0.07255348
#> 4 TCGA-BR-4253  0.03527550 0.01627839       0.2119763      0.16930595
#> 5 TCGA-BR-4256  0.03301111 0.02255326       0.2076062      0.11087587
#> 6 TCGA-BR-4257  0.01447665 0.02230610       0.2052179      0.10247566
#>   Endothelial_EPIC Macrophages_EPIC NKcells_EPIC otherCells_EPIC
#> 1        0.1236952       0.01459545 6.221647e-09       0.5455018
#> 2        0.1819173       0.01529498 3.320319e-10       0.4297812
#> 3        0.1456124       0.01319948 1.854999e-09       0.4582795
#> 4        0.1276442       0.01558499 3.993855e-10       0.4239347
#> 5        0.1768126       0.01491919 9.890730e-11       0.4342218
#> 6        0.1404188       0.01474483 2.891369e-10       0.5003601
```

#### 2.1.3 Method 3: Use MCPcounter to decode the TME contexture.

``` r
mcp<-deconvo_tme(eset = eset_stad,method = "mcpcounter")
#> 
#> >>> Running MCP-counter
head(mcp)
#>             ID T_cells_MCPcounter CD8_T_cells_MCPcounter
#> 1 TCGA-B7-5818           2.260030               2.809712
#> 2 TCGA-BR-4187           2.073911               3.136468
#> 3 TCGA-BR-4201           2.278328               1.108326
#> 4 TCGA-BR-4253           3.694865               4.664873
#> 5 TCGA-BR-4256           2.713699               3.397456
#> 6 TCGA-BR-4257           2.063659               2.537676
#>   Cytotoxic_lymphocytes_MCPcounter NK_cells_MCPcounter B_lineage_MCPcounter
#> 1                         1.594603           0.3867827             1.622539
#> 2                         1.533179           0.4870056             3.248121
#> 3                         1.959712           1.1467281             3.457685
#> 4                         3.833974           1.1215646             2.979089
#> 5                         2.255906           0.8003135             3.012940
#> 6                         1.761330           0.4889650             1.681222
#>   Monocytic_lineage_MCPcounter Myeloid_dendritic_cells_MCPcounter
#> 1                     3.733155                           1.960880
#> 2                     4.494092                           2.270000
#> 3                     4.233543                           1.976727
#> 4                     5.478763                           2.451218
#> 5                     4.880174                           1.932627
#> 6                     4.111373                           2.616387
#>   Neutrophils_MCPcounter Endothelial_cells_MCPcounter Fibroblasts_MCPcounter
#> 1               1.887732                     2.244634               6.447366
#> 2               2.753930                     3.750963               8.809173
#> 3               2.875811                     3.240503               8.369350
#> 4               2.265412                     2.820912               5.819634
#> 5               3.178209                     3.792050               8.504134
#> 6               2.301324                     2.674796               7.466128
```

#### 2.1.4 Method 4: Use xCELL to decode the TME contexture.

``` r
xcell<-deconvo_tme(eset = eset_stad,method = "xcell",arrays = FALSE)
#> [1] "Num. of genes: 10761"
#> Estimating ssGSEA scores for 489 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 4 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
head(xcell)
#>             ID  aDC_xCell Adipocytes_xCell Astrocytes_xCell B-cells_xCell
#> 1 TCGA-B7-5818 0.18857132      0.000000000      0.001588524    0.06640974
#> 2 TCGA-BR-4187 0.06039922      0.008961759      0.132348373    0.05993867
#> 3 TCGA-BR-4201 0.09927959      0.001419346      0.115689899    0.05037067
#> 4 TCGA-BR-4253 0.31057277      0.000000000      0.000000000    0.20181972
#> 5 TCGA-BR-4256 0.31536771      0.007678063      0.126574848    0.04577051
#> 6 TCGA-BR-4257 0.08046210      0.000000000      0.050764977    0.04985408
#>   Basophils_xCell CD4+_memory_T-cells_xCell CD4+_naive_T-cells_xCell
#> 1    1.197904e-01               0.005506425             8.794994e-03
#> 2    5.330722e-18               0.021017929             2.520513e-18
#> 3    0.000000e+00               0.031860833             0.000000e+00
#> 4    7.630865e-02               0.139723991             5.240607e-19
#> 5    2.002366e-01               0.030674768             0.000000e+00
#> 6    4.212687e-02               0.041222540             1.108219e-18
#>   CD4+_T-cells_xCell CD4+_Tcm_xCell CD4+_Tem_xCell CD8+_naive_T-cells_xCell
#> 1       8.219759e-19   3.166981e-18   4.807029e-02              0.007065725
#> 2       1.750760e-18   8.015496e-18   0.000000e+00              0.001486785
#> 3       1.896213e-19   0.000000e+00   0.000000e+00              0.013295314
#> 4       3.739354e-02   0.000000e+00   1.879084e-18              0.017925164
#> 5       8.327629e-19   0.000000e+00   1.951922e-19              0.003720685
#> 6       9.194384e-04   1.344560e-19   1.284826e-04              0.001285986
#>   CD8+_T-cells_xCell CD8+_Tcm_xCell CD8+_Tem_xCell   cDC_xCell
#> 1       2.705866e-02     0.05780763   3.398977e-02 0.007902785
#> 2       6.125270e-03     0.08096435   0.000000e+00 0.397304984
#> 3       4.382314e-20     0.03356524   3.148284e-19 0.115152220
#> 4       1.968640e-01     0.34239476   1.319459e-01 0.186482631
#> 5       0.000000e+00     0.07237217   3.694286e-18 0.075318397
#> 6       8.030273e-03     0.04681635   1.406805e-02 0.178145482
#>   Chondrocytes_xCell Class-switched_memory_B-cells_xCell    CLP_xCell
#> 1       1.338926e-18                          0.02677159 1.866824e-02
#> 2       1.377808e-01                          0.02107003 1.250650e-18
#> 3       1.086657e-01                          0.01405732 1.597001e-02
#> 4       3.747651e-19                          0.10452611 7.853919e-02
#> 5       7.394959e-02                          0.02676998 0.000000e+00
#> 6       0.000000e+00                          0.03164731 2.633398e-02
#>      CMP_xCell   DC_xCell Endothelial_cells_xCell Eosinophils_xCell
#> 1 7.601006e-19 0.01419908              0.01704328      0.000000e+00
#> 2 1.232317e-02 0.06698100              0.15339967      1.788159e-02
#> 3 0.000000e+00 0.03056513              0.03602349      7.931558e-22
#> 4 2.482432e-18 0.04808759              0.03066050      0.000000e+00
#> 5 1.816014e-01 0.07165982              0.13407367      7.368637e-03
#> 6 6.678673e-19 0.03538559              0.03120105      5.710163e-19
#>   Epithelial_cells_xCell Erythrocytes_xCell Fibroblasts_xCell    GMP_xCell
#> 1             0.09724818       0.000000e+00      3.122828e-20 2.428381e-19
#> 2             0.04298145       0.000000e+00      3.951265e-01 2.229407e-01
#> 3             0.12377682       0.000000e+00      1.163427e-01 1.649112e-01
#> 4             0.09196166       4.314514e-04      0.000000e+00 2.571311e-02
#> 5             0.07902924       3.296836e-20      2.422493e-01 3.771340e-01
#> 6             0.13366633       1.709959e-03      1.105295e-02 1.988928e-02
#>   Hepatocytes_xCell   HSC_xCell  iDC_xCell Keratinocytes_xCell
#> 1      3.844542e-19 0.027161695 0.04782129          0.02428391
#> 2      1.939939e-05 0.458945899 1.05023777          0.01307245
#> 3      0.000000e+00 0.173637749 0.34379207          0.06295765
#> 4      0.000000e+00 0.027527034 0.16396619          0.01693531
#> 5      2.165702e-03 0.187160362 0.52171365          0.02455099
#> 6      0.000000e+00 0.009484859 0.30197690          0.03711747
#>   ly_Endothelial_cells_xCell Macrophages_xCell Macrophages_M1_xCell
#> 1                0.007462473        0.08182209           0.04201815
#> 2                0.054556958        0.05017048           0.04886551
#> 3                0.011455881        0.03802175           0.02375797
#> 4                0.003675218        0.13331553           0.07957648
#> 5                0.062787726        0.09967781           0.07835602
#> 6                0.006822980        0.09802533           0.04788241
#>   Macrophages_M2_xCell Mast_cells_xCell Megakaryocytes_xCell Melanocytes_xCell
#> 1           0.05454049     1.526811e-23         0.0000000000       0.006662637
#> 2           0.03419120     1.203509e-02         0.0065153345       0.000000000
#> 3           0.01132362     3.641574e-03         0.0006099907       0.001937422
#> 4           0.04780682     1.292472e-02         0.0005100685       0.002673386
#> 5           0.03678692     1.401345e-02         0.0035463550       0.011340987
#> 6           0.05439590     2.612071e-03         0.0000000000       0.005075628
#>   Memory_B-cells_xCell  MEP_xCell Mesangial_cells_xCell Monocytes_xCell
#> 1         3.378453e-03 0.06138737          0.0219200604     0.002916951
#> 2         3.016879e-18 0.03722859          0.0312161855     0.088081545
#> 3         1.022182e-19 0.05058237          0.0127222688     0.066281553
#> 4         7.796270e-03 0.06025132          0.0004618847     0.067217815
#> 5         2.270415e-18 0.00000000          0.0301026547     0.196124036
#> 6         0.000000e+00 0.05350684          0.0249753484     0.053365578
#>      MPP_xCell  MSC_xCell mv_Endothelial_cells_xCell Myocytes_xCell
#> 1 0.000000e+00 0.28883315               0.0163575193   1.586153e-20
#> 2 1.542792e-18 0.00000000               0.0559953012   0.000000e+00
#> 3 0.000000e+00 0.11363803               0.0061850890   3.068153e-03
#> 4 0.000000e+00 0.00000000               0.0023791004   2.015474e-19
#> 5 2.615186e-02 0.08890403               0.0544846147   2.001219e-03
#> 6 1.058351e-18 0.18110478               0.0008244167   4.160355e-19
#>   naive_B-cells_xCell Neurons_xCell Neutrophils_xCell NK_cells_xCell
#> 1        8.065775e-20  0.000000e+00        0.00000000   8.222944e-19
#> 2        0.000000e+00  6.785817e-03        0.01315471   9.766268e-19
#> 3        8.273345e-19  2.859851e-04        0.01976858   2.780922e-03
#> 4        1.777302e-04  0.000000e+00        0.01858329   6.092509e-02
#> 5        3.489858e-18  7.818667e-04        0.05355418   0.000000e+00
#> 6        0.000000e+00  3.665537e-05        0.01382641   1.139179e-02
#>      NKT_xCell Osteoblast_xCell  pDC_xCell Pericytes_xCell Plasma_cells_xCell
#> 1 6.752586e-02     1.444113e-02 0.05312166    3.557793e-02       1.571993e-02
#> 2 0.000000e+00     5.679401e-18 0.01084697    0.000000e+00       1.494602e-02
#> 3 0.000000e+00     0.000000e+00 0.03676934    1.069611e-20       1.678725e-02
#> 4 1.580197e-02     1.637449e-20 0.14599535    0.000000e+00       3.015153e-02
#> 5 1.543747e-18     7.923329e-20 0.02715485    0.000000e+00       5.561710e-20
#> 6 6.506464e-03     0.000000e+00 0.04458637    1.328906e-18       9.484202e-03
#>   Platelets_xCell Preadipocytes_xCell pro_B-cells_xCell Sebocytes_xCell
#> 1    0.000000e+00        2.007880e-18      1.320989e-02     0.011902720
#> 2    0.000000e+00        1.466821e-01      0.000000e+00     0.008430347
#> 3    0.000000e+00        0.000000e+00      1.739112e-18     0.019636957
#> 4    5.948610e-19        1.614552e-02      2.111642e-02     0.007138945
#> 5    1.988042e-18        3.228353e-19      2.179121e-02     0.013391420
#> 6    0.000000e+00        1.935703e-02      7.385504e-03     0.015259198
#>   Skeletal_muscle_xCell Smooth_muscle_xCell Tgd_cells_xCell Th1_cells_xCell
#> 1          2.390410e-18           0.0000000    0.000000e+00      0.22218272
#> 2          0.000000e+00           0.3461178    0.000000e+00      0.00000000
#> 3          1.477955e-19           0.3273036    1.547124e-18      0.03022514
#> 4          0.000000e+00           0.1288147    3.576465e-02      0.26037555
#> 5          3.060294e-03           0.2718761    1.622825e-18      0.04087900
#> 6          2.102653e-03           0.1912370    1.080367e-02      0.07137173
#>   Th2_cells_xCell  Tregs_xCell ImmuneScore_xCell StromaScore_xCell
#> 1      0.05082409 7.631932e-03         0.1282710       0.008521642
#> 2      0.05042155 1.750307e-20         0.2095789       0.278743983
#> 3      0.16424222 2.771978e-19         0.1409535       0.076892768
#> 4      0.34765629 6.709947e-02         0.5180875       0.015330251
#> 5      0.19474253 6.054240e-20         0.3254456       0.192000530
#> 6      0.15053676 1.604214e-02         0.1822737       0.021126998
#>   MicroenvironmentScore_xCell
#> 1                   0.1367927
#> 2                   0.4883229
#> 3                   0.2178462
#> 4                   0.5334178
#> 5                   0.5174462
#> 6                   0.2034007
```

#### 2.1.5 Method 5: Use ESTIMATE to estimate tumor purity, and generate immune score and stromal score.

``` r
estimate<-deconvo_tme(eset = eset_stad,method = "estimate")
#> [1] "Merged dataset includes 10156 genes (256 mismatched)."
#> [1] "1 gene set: StromalSignature  overlap= 139"
#> [1] "2 gene set: ImmuneSignature  overlap= 140"
head(estimate)
#>             ID StromalScore_estimate ImmuneScore_estimate
#> 1 TCGA-B7-5818             -110.5759             1761.803
#> 2 TCGA-BR-4187             2054.3902             2208.033
#> 3 TCGA-BR-4201             1410.9584             1769.551
#> 4 TCGA-BR-4253              483.0597             2905.491
#> 5 TCGA-BR-4256             1659.2327             2540.646
#> 6 TCGA-BR-4257              831.1475             1721.680
#>   ESTIMATEScore_estimate TumorPurity_estimate
#> 1               1651.227            0.6619581
#> 2               4262.423            0.3336142
#> 3               3180.509            0.4785014
#> 4               3388.551            0.4514674
#> 5               4199.878            0.3422549
#> 6               2552.827            0.5572612
```

#### 2.1.6 Method 6: Use TIMER to decode the TME contexture.

``` r
timer<-deconvo_tme(eset = eset_stad,method = "timer",group_list = rep("stad",dim(eset_stad)[2]))
#> [1] "Outlier genes: FTL IGF2 IGHA1 IGHM IGKC IGKV4-1 LYZ MT-ATP6 MT-CO1 MT-CO2 MT-CO3 MT-CYB MT-ND1 MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-RNR1 MT-RNR2 MT-TP PGC"
#> Standardizing Data across genes
head(timer)
#>             ID B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER Neutrophil_TIMER
#> 1 TCGA-B7-5818   0.09544149        0.1304849        0.1911229        0.1135122
#> 2 TCGA-BR-4187   0.09827376        0.1305894        0.2020704        0.1252187
#> 3 TCGA-BR-4201   0.09959006        0.1221931        0.1997947        0.1274679
#> 4 TCGA-BR-4253   0.10110096        0.1310875        0.2404586        0.1294138
#> 5 TCGA-BR-4256   0.09451320        0.1325623        0.2126199        0.1372099
#> 6 TCGA-BR-4257   0.09071870        0.1256934        0.1987924        0.1206287
#>   Macrophage_TIMER  DC_TIMER
#> 1       0.03885340 0.5145410
#> 2       0.06265559 0.5079503
#> 3       0.05506044 0.5097286
#> 4       0.04978486 0.5305424
#> 5       0.05810445 0.5197208
#> 6       0.05125203 0.5153890
```

#### 2.1.7 Method 7: Use quanTIseq to decode the TME contexture.

``` r
quantiseq<-deconvo_tme(eset = eset_stad, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE,method = "quantiseq")
#> 
#> Running quanTIseq deconvolution module
#> Gene expression normalization and re-annotation (arrays: FALSE)
#> Removing 17 noisy genes
#> Removing 15 genes with high expression in tumors
#> Signature genes found in data set: 138/138 (100%)
#> Mixture deconvolution (method: lsei)
#> Deconvolution sucessful!
head(quantiseq)
#>             ID B_cells_quantiseq Macrophages_M1_quantiseq
#> 1 TCGA-B7-5818       0.008971498               0.09430627
#> 2 TCGA-BR-4187       0.026915913               0.01551815
#> 3 TCGA-BR-4201       0.016043863               0.05700276
#> 4 TCGA-BR-4253       0.024036081               0.17147320
#> 5 TCGA-BR-4256       0.012881564               0.12627623
#> 6 TCGA-BR-4257       0.014768360               0.06853805
#>   Macrophages_M2_quantiseq Monocytes_quantiseq Neutrophils_quantiseq
#> 1               0.02986523                   0           0.003866517
#> 2               0.16912650                   0           0.030072763
#> 3               0.06610113                   0           0.047839086
#> 4               0.12777493                   0           0.068495811
#> 5               0.09222656                   0           0.066754502
#> 6               0.06729849                   0           0.027692418
#>   NK_cells_quantiseq T_cells_CD4_quantiseq T_cells_CD8_quantiseq
#> 1        0.021595936            0.00000000           0.018980102
#> 2        0.011939039            0.01051120           0.016382968
#> 3        0.009433017            0.00000000           0.004489431
#> 4        0.005478200            0.00000000           0.103014149
#> 5        0.012183634            0.00000000           0.031481991
#> 6        0.012469349            0.02979917           0.015648317
#>   Tregs_quantiseq Dendritic_cells_quantiseq Other_quantiseq
#> 1      0.01757416                         0       0.8048403
#> 2      0.02165007                         0       0.6978834
#> 3      0.03983729                         0       0.7592534
#> 4      0.06604435                         0       0.4336833
#> 5      0.05510340                         0       0.6030921
#> 6      0.02262833                         0       0.7411575
res<-cell_bar_plot(input = quantiseq[1:12,], title = "quanTIseq Cell Fraction")
#> There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap
```

<img src="man/figuresunnamed-chunk-14-1.png" width="100%" />

#### 2.1.8 Method 8: Use IPS to estimate immune phenotype.

``` r
ips<-deconvo_tme(eset = eset_stad,method = "ips",plot= FALSE)
#>    Mode    TRUE 
#> logical     161 
#> [1] GENE   NAME   CLASS  WEIGHT
#> <0 rows> (or 0-length row.names)
#> [1] GENE   NAME   CLASS  WEIGHT
#> <0 rows> (or 0-length row.names)
head(ips)
#>             ID  MHC_IPS   EC_IPS    SC_IPS     CP_IPS   AZ_IPS IPS_IPS
#> 1 TCGA-B7-5818 3.557041 1.088319 -1.365492 -0.5758197 2.704049       9
#> 2 TCGA-BR-4187 3.367949 1.182982 -1.709257 -0.3082528 2.533421       8
#> 3 TCGA-BR-4201 3.037724 1.223938 -1.621719 -0.4394294 2.200514       7
#> 4 TCGA-BR-4253 3.787724 1.626261 -1.738435 -1.2510594 2.424492       8
#> 5 TCGA-BR-4256 3.654959 1.401826 -1.954562 -0.8474067 2.254817       8
#> 6 TCGA-BR-4257 3.111207 1.230521 -1.741116 -0.5290984 2.071514       7
```

#### 2.1.9 Combine the above deconvolution results for subsequent analyses.

``` r
tme_combine<-cibersort %>% 
  inner_join(.,mcp,by = "ID") %>% 
  inner_join(.,xcell,by = "ID") %>%
  inner_join(.,epic,by = "ID") %>% 
  inner_join(.,estimate,by = "ID") %>% 
  inner_join(.,timer,by = "ID") %>% 
  inner_join(.,quantiseq,by = "ID") %>% 
  inner_join(.,ips,by = "ID")
dim(tme_combine)
#> [1]  98 138
colnames(tme_combine)
#>   [1] "ID"                                    
#>   [2] "B_cells_naive_CIBERSORT"               
#>   [3] "B_cells_memory_CIBERSORT"              
#>   [4] "Plasma_cells_CIBERSORT"                
#>   [5] "T_cells_CD8_CIBERSORT"                 
#>   [6] "T_cells_CD4_naive_CIBERSORT"           
#>   [7] "T_cells_CD4_memory_resting_CIBERSORT"  
#>   [8] "T_cells_CD4_memory_activated_CIBERSORT"
#>   [9] "T_cells_follicular_helper_CIBERSORT"   
#>  [10] "T_cells_regulatory_(Tregs)_CIBERSORT"  
#>  [11] "T_cells_gamma_delta_CIBERSORT"         
#>  [12] "NK_cells_resting_CIBERSORT"            
#>  [13] "NK_cells_activated_CIBERSORT"          
#>  [14] "Monocytes_CIBERSORT"                   
#>  [15] "Macrophages_M0_CIBERSORT"              
#>  [16] "Macrophages_M1_CIBERSORT"              
#>  [17] "Macrophages_M2_CIBERSORT"              
#>  [18] "Dendritic_cells_resting_CIBERSORT"     
#>  [19] "Dendritic_cells_activated_CIBERSORT"   
#>  [20] "Mast_cells_resting_CIBERSORT"          
#>  [21] "Mast_cells_activated_CIBERSORT"        
#>  [22] "Eosinophils_CIBERSORT"                 
#>  [23] "Neutrophils_CIBERSORT"                 
#>  [24] "P-value_CIBERSORT"                     
#>  [25] "Correlation_CIBERSORT"                 
#>  [26] "RMSE_CIBERSORT"                        
#>  [27] "T_cells_MCPcounter"                    
#>  [28] "CD8_T_cells_MCPcounter"                
#>  [29] "Cytotoxic_lymphocytes_MCPcounter"      
#>  [30] "NK_cells_MCPcounter"                   
#>  [31] "B_lineage_MCPcounter"                  
#>  [32] "Monocytic_lineage_MCPcounter"          
#>  [33] "Myeloid_dendritic_cells_MCPcounter"    
#>  [34] "Neutrophils_MCPcounter"                
#>  [35] "Endothelial_cells_MCPcounter"          
#>  [36] "Fibroblasts_MCPcounter"                
#>  [37] "aDC_xCell"                             
#>  [38] "Adipocytes_xCell"                      
#>  [39] "Astrocytes_xCell"                      
#>  [40] "B-cells_xCell"                         
#>  [41] "Basophils_xCell"                       
#>  [42] "CD4+_memory_T-cells_xCell"             
#>  [43] "CD4+_naive_T-cells_xCell"              
#>  [44] "CD4+_T-cells_xCell"                    
#>  [45] "CD4+_Tcm_xCell"                        
#>  [46] "CD4+_Tem_xCell"                        
#>  [47] "CD8+_naive_T-cells_xCell"              
#>  [48] "CD8+_T-cells_xCell"                    
#>  [49] "CD8+_Tcm_xCell"                        
#>  [50] "CD8+_Tem_xCell"                        
#>  [51] "cDC_xCell"                             
#>  [52] "Chondrocytes_xCell"                    
#>  [53] "Class-switched_memory_B-cells_xCell"   
#>  [54] "CLP_xCell"                             
#>  [55] "CMP_xCell"                             
#>  [56] "DC_xCell"                              
#>  [57] "Endothelial_cells_xCell"               
#>  [58] "Eosinophils_xCell"                     
#>  [59] "Epithelial_cells_xCell"                
#>  [60] "Erythrocytes_xCell"                    
#>  [61] "Fibroblasts_xCell"                     
#>  [62] "GMP_xCell"                             
#>  [63] "Hepatocytes_xCell"                     
#>  [64] "HSC_xCell"                             
#>  [65] "iDC_xCell"                             
#>  [66] "Keratinocytes_xCell"                   
#>  [67] "ly_Endothelial_cells_xCell"            
#>  [68] "Macrophages_xCell"                     
#>  [69] "Macrophages_M1_xCell"                  
#>  [70] "Macrophages_M2_xCell"                  
#>  [71] "Mast_cells_xCell"                      
#>  [72] "Megakaryocytes_xCell"                  
#>  [73] "Melanocytes_xCell"                     
#>  [74] "Memory_B-cells_xCell"                  
#>  [75] "MEP_xCell"                             
#>  [76] "Mesangial_cells_xCell"                 
#>  [77] "Monocytes_xCell"                       
#>  [78] "MPP_xCell"                             
#>  [79] "MSC_xCell"                             
#>  [80] "mv_Endothelial_cells_xCell"            
#>  [81] "Myocytes_xCell"                        
#>  [82] "naive_B-cells_xCell"                   
#>  [83] "Neurons_xCell"                         
#>  [84] "Neutrophils_xCell"                     
#>  [85] "NK_cells_xCell"                        
#>  [86] "NKT_xCell"                             
#>  [87] "Osteoblast_xCell"                      
#>  [88] "pDC_xCell"                             
#>  [89] "Pericytes_xCell"                       
#>  [90] "Plasma_cells_xCell"                    
#>  [91] "Platelets_xCell"                       
#>  [92] "Preadipocytes_xCell"                   
#>  [93] "pro_B-cells_xCell"                     
#>  [94] "Sebocytes_xCell"                       
#>  [95] "Skeletal_muscle_xCell"                 
#>  [96] "Smooth_muscle_xCell"                   
#>  [97] "Tgd_cells_xCell"                       
#>  [98] "Th1_cells_xCell"                       
#>  [99] "Th2_cells_xCell"                       
#> [100] "Tregs_xCell"                           
#> [101] "ImmuneScore_xCell"                     
#> [102] "StromaScore_xCell"                     
#> [103] "MicroenvironmentScore_xCell"           
#> [104] "Bcells_EPIC"                           
#> [105] "CAFs_EPIC"                             
#> [106] "CD4_Tcells_EPIC"                       
#> [107] "CD8_Tcells_EPIC"                       
#> [108] "Endothelial_EPIC"                      
#> [109] "Macrophages_EPIC"                      
#> [110] "NKcells_EPIC"                          
#> [111] "otherCells_EPIC"                       
#> [112] "StromalScore_estimate"                 
#> [113] "ImmuneScore_estimate"                  
#> [114] "ESTIMATEScore_estimate"                
#> [115] "TumorPurity_estimate"                  
#> [116] "B_cell_TIMER"                          
#> [117] "T_cell_CD4_TIMER"                      
#> [118] "T_cell_CD8_TIMER"                      
#> [119] "Neutrophil_TIMER"                      
#> [120] "Macrophage_TIMER"                      
#> [121] "DC_TIMER"                              
#> [122] "B_cells_quantiseq"                     
#> [123] "Macrophages_M1_quantiseq"              
#> [124] "Macrophages_M2_quantiseq"              
#> [125] "Monocytes_quantiseq"                   
#> [126] "Neutrophils_quantiseq"                 
#> [127] "NK_cells_quantiseq"                    
#> [128] "T_cells_CD4_quantiseq"                 
#> [129] "T_cells_CD8_quantiseq"                 
#> [130] "Tregs_quantiseq"                       
#> [131] "Dendritic_cells_quantiseq"             
#> [132] "Other_quantiseq"                       
#> [133] "MHC_IPS"                               
#> [134] "EC_IPS"                                
#> [135] "SC_IPS"                                
#> [136] "CP_IPS"                                
#> [137] "AZ_IPS"                                
#> [138] "IPS_IPS"
```

## 3 Signature score estimation

IOBR integrates 255 published signature gene sets, involving tumor
microenvironment, tumor metabolism, m6A, exosomes, microsatellite
instability, and tertiary lymphoid structure. Running the function
`signature_collection_citation` to attain the source papers. The
function `signature_collection` returns the detail signature genes of
all given signatures.

### 3.1 Obtain the included signatures

``` r
library(IOBR)
#TME associated signatures
names(signature_tme)[1:20]
#>  [1] "CD_8_T_effector"            "DDR"                       
#>  [3] "APM"                        "Immune_Checkpoint"         
#>  [5] "CellCycle_Reg"              "Pan_F_TBRs"                
#>  [7] "Histones"                   "EMT1"                      
#>  [9] "EMT2"                       "EMT3"                      
#> [11] "WNT_target"                 "FGFR3_related"             
#> [13] "Cell_cycle"                 "Mismatch_Repair"           
#> [15] "Homologous_recombination"   "Nucleotide_excision_repair"
#> [17] "DNA_replication"            "Base_excision_repair"      
#> [19] "TMEscoreA_CIR"              "TMEscoreB_CIR"
#Metabolism related signatures
names(signature_metabolism)[1:20]
#>  [1] "Cardiolipin_Metabolism"                    
#>  [2] "Cardiolipin_Biosynthesis"                  
#>  [3] "Cholesterol_Biosynthesis"                  
#>  [4] "Citric_Acid_Cycle"                         
#>  [5] "Cyclooxygenase_Arachidonic_Acid_Metabolism"
#>  [6] "Prostaglandin_Biosynthesis"                
#>  [7] "Purine_Biosynthesis"                       
#>  [8] "Pyrimidine_Biosynthesis"                   
#>  [9] "Dopamine_Biosynthesis"                     
#> [10] "Epinephrine_Biosynthesis"                  
#> [11] "Norepinephrine_Biosynthesis"               
#> [12] "Fatty_Acid_Degradation"                    
#> [13] "Fatty_Acid_Elongation"                     
#> [14] "Fatty_Acid_Biosynthesis"                   
#> [15] "Folate_One_Carbon_Metabolism"              
#> [16] "Folate_biosynthesis"                       
#> [17] "Gluconeogenesis"                           
#> [18] "Glycolysis"                                
#> [19] "Glycogen_Biosynthesis"                     
#> [20] "Glycogen_Degradation"
#Signatures associated with biomedical basic research: such as m6A and exosomes
names(signature_tumor)
#>  [1] "Nature_metabolism_Hypoxia"                
#>  [2] "Winter_hypoxia_signature"                 
#>  [3] "Hu_hypoxia_signature"                     
#>  [4] "Molecular_Cancer_m6A"                     
#>  [5] "MT_exosome"                               
#>  [6] "SR_exosome"                               
#>  [7] "Positive_regulation_of_exosomal_secretion"
#>  [8] "Negative_regulation_of_exosomal_secretion"
#>  [9] "Exosomal_secretion"                       
#> [10] "Exosome_assembly"                         
#> [11] "Extracellular_vesicle_biogenesis"         
#> [12] "MC_Review_Exosome1"                       
#> [13] "MC_Review_Exosome2"                       
#> [14] "CMLS_Review_Exosome"                      
#> [15] "Ferroptosis"                              
#> [16] "EV_Cell_2020"
#signature collection including all aforementioned signatures 
names(signature_collection)[1:20]
#>  [1] "CD_8_T_effector"            "DDR"                       
#>  [3] "APM"                        "Immune_Checkpoint"         
#>  [5] "CellCycle_Reg"              "Pan_F_TBRs"                
#>  [7] "Histones"                   "EMT1"                      
#>  [9] "EMT2"                       "EMT3"                      
#> [11] "WNT_target"                 "FGFR3_related"             
#> [13] "Cell_cycle"                 "Mismatch_Repair"           
#> [15] "Homologous_recombination"   "Nucleotide_excision_repair"
#> [17] "DNA_replication"            "Base_excision_repair"      
#> [19] "TMEscoreA_CIR"              "TMEscoreB_CIR"
```

### 3.2 Methods for signature calculation

### 3.2.1 Calculate TME associated signatures-(through PCA method).

``` r
sig_tme<-calculate_sig_score(pdata = NULL,
                             eset = eset_stad,
                             signature = signature_tme,
                             method = "pca",
                             mini_gene_count = 2)
#> 
#> >>> Calculating signature score with PCA method
sig_tme[1:5,1:10]
#>   Index           ID CD_8_T_effector        DDR        APM Immune_Checkpoint
#> 1     1 TCGA-B7-5818       2.4078956  1.0184565  0.7822385        0.83337098
#> 2     2 TCGA-BR-4187      -0.5814857 -1.9366977  0.1720316        0.07939075
#> 3     3 TCGA-BR-4201       1.0856177  0.7347264 -0.4170782        0.67105977
#> 4     4 TCGA-BR-4253       5.6853246  3.1936537  1.3327474        3.09475114
#> 5     5 TCGA-BR-4256       2.2337399  0.3293904  1.2457966        1.79448256
#>   CellCycle_Reg Pan_F_TBRs   Histones        EMT1
#> 1    -0.4629385 -1.6210226  1.4364836 -2.52906550
#> 2    -0.5838805  2.9507234  0.8786425  0.61053612
#> 3    -0.2466607 -0.3590154 -0.1593197 -0.48920322
#> 4    -0.5532432 -3.9517259  1.0731159 -2.70655819
#> 5     0.9852912  1.9124010 -0.2739138  0.07716814
```

### 3.2.2 Calculate TME associated signatures-(through ssGSEA method).

``` r
sig_tme<-calculate_sig_score(pdata = NULL,
                             eset = eset_stad,
                             signature = signature_tme,
                             method = "ssgsea",
                             mini_gene_count = 5)
#> Estimating ssGSEA scores for 98 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_tme[1:5,1:10]
#>             ID Index CD_8_T_effector       DDR       APM Immune_Checkpoint
#> 1 TCGA-B7-5818     1       0.4558828 0.3683733 0.5549803         0.3396046
#> 2 TCGA-BR-4187     2       0.3337869 0.3432917 0.5373613         0.2776199
#> 3 TCGA-BR-4201     3       0.4002791 0.3702784 0.5476764         0.3070661
#> 4 TCGA-BR-4253     4       0.5390946 0.3943726 0.5730135         0.4530388
#> 5 TCGA-BR-4256     5       0.4471666 0.3617588 0.5686462         0.3813583
#>   CellCycle_Reg Pan_F_TBRs      EMT1      EMT2
#> 1     0.3621017  0.4787365 0.4141455 0.3852234
#> 2     0.4186531  0.5127919 0.5155580 0.4720354
#> 3     0.4427704  0.5123139 0.5038751 0.4627777
#> 4     0.4416835  0.4306756 0.4192358 0.3353618
#> 5     0.4515601  0.5288249 0.5026945 0.4548524
```

### 3.2.3 Calculate metabolism related signatures.

``` r
sig_meta<-calculate_sig_score(pdata = NULL,
                              eset = eset_stad,
                              signature = signature_metabolism,
                              method = "pca",
                              mini_gene_count = 2)
#> 
#> >>> Calculating signature score with PCA method
sig_meta[1:5,1:10]
#>   Index           ID Cardiolipin_Metabolism Cardiolipin_Biosynthesis
#> 1     1 TCGA-B7-5818             -0.6344980              -0.04527508
#> 2     2 TCGA-BR-4187              0.3493315              -0.15579130
#> 3     3 TCGA-BR-4201              0.7894225               0.37846339
#> 4     4 TCGA-BR-4253             -0.9625402              -0.04090572
#> 5     5 TCGA-BR-4256              0.3106262              -0.08954902
#>   Cholesterol_Biosynthesis Citric_Acid_Cycle
#> 1                1.9349209       -0.39043546
#> 2               -1.2167907        1.83912921
#> 3                1.4816557       -0.52630368
#> 4               -0.6966216        0.14933490
#> 5                0.3527329        0.06397952
#>   Cyclooxygenase_Arachidonic_Acid_Metabolism Prostaglandin_Biosynthesis
#> 1                                  -1.031111                 -0.1787793
#> 2                                   1.259596                  1.6076195
#> 3                                   1.336709                  0.7230978
#> 4                                   0.403215                  0.3138395
#> 5                                   1.059021                  0.4990814
#>   Purine_Biosynthesis Pyrimidine_Biosynthesis
#> 1          -0.1610375               0.8975888
#> 2          -0.5924434              -0.4154989
#> 3           0.6757635               0.5565299
#> 4           1.8088674               1.5473453
#> 5           0.2081434               0.4303811
```

### 3.2.4 Calculate all collected signature scores (integrating three methods: PCA, ssGSEA and z-score).

``` r
sig_res<-calculate_sig_score(pdata = NULL,
                             eset = eset_stad,
                             signature = signature_collection,
                             method = "integration",
                             mini_gene_count = 2)
#> Estimating ssGSEA scores for 213 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_res[1:5,1:10]
#>             ID Index CD_8_T_effector_PCA    DDR_PCA    APM_PCA
#> 1 TCGA-B7-5818     1           2.4078956  1.0184565  0.7822385
#> 2 TCGA-BR-4187     2          -0.5814857 -1.9366977  0.1720316
#> 3 TCGA-BR-4201     3           1.0856177  0.7347264 -0.4170782
#> 4 TCGA-BR-4253     4           5.6853246  3.1936537  1.3327474
#> 5 TCGA-BR-4256     5           2.2337399  0.3293904  1.2457966
#>   Immune_Checkpoint_PCA CellCycle_Reg_PCA Pan_F_TBRs_PCA Histones_PCA
#> 1            0.83337098        -0.4629385     -1.6210226    1.4364836
#> 2            0.07939075        -0.5838805      2.9507234    0.8786425
#> 3            0.67105977        -0.2466607     -0.3590154   -0.1593197
#> 4            3.09475114        -0.5532432     -3.9517259    1.0731159
#> 5            1.79448256         0.9852912      1.9124010   -0.2739138
#>      EMT1_PCA
#> 1 -2.52906550
#> 2  0.61053612
#> 3 -0.48920322
#> 4 -2.70655819
#> 5  0.07716814
```

### 3.2.5 The signature gene sets derived from GO, KEGG, HALLMARK and REACTOME datasets.

IOBR also enrolls the signature gene sets, containing GO, KEGG, HALLMARK
and REACTOME gene sets obtained from
[MsigDB](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H). The
ssGSEA method is recommended for these signature estimation. This
process may take a while for big datasets or calculating a large number
of signatures.

``` r
sig_hallmark<-calculate_sig_score(pdata = NULL,
                                  eset = eset_stad,
                                  signature = hallmark,
                                  method = "ssgsea",
                                  mini_gene_count = 2)
#> Estimating ssGSEA scores for 50 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_hallmark[1:5,1:10]
#>             ID Index HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_HYPOXIA
#> 1 TCGA-B7-5818     1                        0.9111075        0.7935822
#> 2 TCGA-BR-4187     2                        0.9228259        0.8646110
#> 3 TCGA-BR-4201     3                        0.9797193        0.8695803
#> 4 TCGA-BR-4253     4                        0.9519414        0.7538464
#> 5 TCGA-BR-4256     5                        1.0138328        0.8630841
#>   HALLMARK_CHOLESTEROL_HOMEOSTASIS HALLMARK_MITOTIC_SPINDLE
#> 1                        0.8865270                0.8661656
#> 2                        0.8888053                0.8581112
#> 3                        0.9358607                0.8842899
#> 4                        0.8582068                0.9072948
#> 5                        0.9073095                0.8897920
#>   HALLMARK_WNT_BETA_CATENIN_SIGNALING HALLMARK_TGF_BETA_SIGNALING
#> 1                           0.7279822                   0.9356868
#> 2                           0.7598011                   0.9541815
#> 3                           0.7137900                   0.9665762
#> 4                           0.6828606                   0.9357797
#> 5                           0.7549524                   0.9777782
#>   HALLMARK_IL6_JAK_STAT3_SIGNALING HALLMARK_DNA_REPAIR
#> 1                        0.8300878           0.9419390
#> 2                        0.8809411           0.9070991
#> 3                        0.8805711           0.9296297
#> 4                        0.9070770           0.9356410
#> 5                        0.9304483           0.9265556
```

## 4 Find phenotype relevant signatures

### 4.1 Data of IMvigor210 immuntherapy cohort is used for batch analysis of phenotype-related sigantures.

``` r
## Load the signatures and Cell fractions calculated in previous steps by IOBR.
imvigor210_sig[1:5,1:10]
#>                ID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT
#> 1 SAMf2ce197162ce              0.02900304              0.045657003
#> 2 SAM698d8d76b934              0.08119867              0.001266419
#> 3 SAMc1b27bc16435              0.01237173              0.000000000
#> 4 SAM85e41e7f33f9              0.00000000              0.000000000
#> 5 SAMf275eb859a39              0.00000000              0.000000000
#>   Plasma_cells_CIBERSORT T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1            0.000000000                     0                   0.1316134
#> 2            0.000000000                     0                   0.0000000
#> 3            0.001273053                     0                   0.1693195
#> 4            0.001186289                     0                   0.0000000
#> 5            0.009493367                     0                   0.0000000
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                            0.2390998                             0.04564597
#> 2                            0.1963620                             0.03865429
#> 3                            0.3168480                             0.08919590
#> 4                            0.3188597                             0.04426192
#> 5                            0.3923644                             0.04485934
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                          0.00000000                                    0
#> 2                          0.00822820                                    0
#> 3                          0.00000000                                    0
#> 4                          0.01379417                                    0
#> 5                          0.00000000                                    0

# Check therapeutic response and survival outcome in the phenotype data.
head(imvigor210_pdata)
#> # A tibble: 6 x 12
#>   ID    BOR   BOR_binary OS_days OS_status Mutation_Load Neo_antigen_Load
#>   <chr> <chr> <chr>      <chr>   <chr>     <chr>         <chr>           
#> 1 SAM0~ NA    NA         57.166~ 1         NA            NA              
#> 2 SAM0~ SD    NR         469.15~ 1         18            4.6862745099999~
#> 3 SAM0~ PD    NR         263.16~ 1         1             0.31372549      
#> 4 SAM0~ PD    NR         74.907~ 1         44            6.1960784310000~
#> 5 SAM0~ NA    NA         20.698~ 0         50            NA              
#> 6 SAM0~ SD    NR         136.01~ 1         2             1.4705882349999~
#> # ... with 5 more variables: CD_8_T_effector <dbl>, Immune_Checkpoint <dbl>,
#> #   Pan_F_TBRs <chr>, Mismatch_Repair <chr>, TumorPurity <dbl>
```

### 4.2 sig\_group() function for batch visualization of response associated signatures.

``` r
# The signature group list.
names(sig_group)
#>  [1] "tumor_signature"        "EMT"                    "io_biomarkers"         
#>  [4] "immu_microenvironment"  "immu_suppression"       "immu_exclusion"        
#>  [7] "immu_exhaustion"        "TCR_BCR"                "tme_signatures1"       
#> [10] "tme_signatures2"        "Bcells"                 "Tcells"                
#> [13] "DCs"                    "Macrophages"            "Neutrophils"           
#> [16] "Monocytes"              "CAFs"                   "NK"                    
#> [19] "tme_cell_types"         "CIBERSORT"              "MCPcounter"            
#> [22] "EPIC"                   "xCell"                  "quanTIseq"             
#> [25] "ESTIMATE"               "IPS"                    "TIMER"                 
#> [28] "fatty_acid_metabolism"  "hypoxia_signature"      "cholesterol_metabolism"
#> [31] "Metabolism"             "hallmark"               "hallmark1"             
#> [34] "hallmark2"              "hallmark3"              "Rooney_et_al"          
#> [37] "Bindea_et_al"           "Li_et_al"               "Peng_et_al"

# The tumor-relevant signatures in first group.
names(sig_group)[1]
#> [1] "tumor_signature"
sig_group[[1]]
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

# The signatures associated with immunotherapy biomarkers.
names(sig_group)[3]
#> [1] "io_biomarkers"
sig_group[[3]]
#>  [1] "TMEscore_CIR"                    "TMEscoreA_CIR"                  
#>  [3] "TMEscoreB_CIR"                   "T_cell_inflamed_GEP_Ayers_et_al"
#>  [5] "CD_8_T_effector"                 "IPS_IPS"                        
#>  [7] "Immune_Checkpoint"               "Exhausted_CD8_Danaher_et_al"    
#>  [9] "Pan_F_TBRs"                      "Mismatch_Repair"                
#> [11] "APM"

# The signatures of immunosuppression.
names(sig_group)[5]
#> [1] "immu_suppression"
sig_group[[5]]
#> [1] "Pan_F_TBRs"                          
#> [2] "Fibroblasts_MCPcounter"              
#> [3] "Immune_Checkpoint"                   
#> [4] "Exhausted_CD8_Danaher_et_al"         
#> [5] "MDSC_Wang_et_al"                     
#> [6] "Macrophages_M2_cibersort"            
#> [7] "Tregs_quantiseq"                     
#> [8] "T_cells_regulatory_(Tregs)_CIBERSORT"
```

### 4.3 Construct phenotype group and bath visualize correlation betweeen therapy response and selected signatures.

``` r
pdata_group<-imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR","BOR_binary")]
res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_data = imvigor210_sig,
                   id2 = "ID",
                   target= NULL,
                   group = "BOR_binary",
                   is_target_continuous = FALSE,
                   padj_cutoff = 1,
                   index = 1,
                   category = "signature",
                   signature_group =sig_group[c(1,3,5)],
                   ProjectID = "IMvigor210",
                   palette_box = "paired1",
                   palette_corplot = "pheatmap",
                   palette_heatmap = 2,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>>  Processing signature: tumor_signature"
```

<img src="man/figuresunnamed-chunk-25-1.png" width="100%" /><img src="man/figuresunnamed-chunk-25-2.png" width="100%" />

    #> [1] ">>>  Processing signature: io_biomarkers"

<img src="man/figuresunnamed-chunk-25-3.png" width="100%" /><img src="man/figuresunnamed-chunk-25-4.png" width="100%" />

    #> [1] ">>>  Processing signature: immu_suppression"

<img src="man/figuresunnamed-chunk-25-5.png" width="100%" /><img src="man/figuresunnamed-chunk-25-6.png" width="100%" />

### 4.4 Estemate the expression of response relevant signature genes.

``` r

imvigor210_eset[1:5,1:6]
#>           SAMf2ce197162ce SAM698d8d76b934 SAMc1b27bc16435 SAM85e41e7f33f9
#> LINC00696      -2.8043469      -3.1185484      -3.0885166       -3.660966
#> LINC01126      -2.1876756      -1.0233912      -1.7100050       -3.008889
#> LINC01460      -6.7112375      -4.2560520      -1.0181273       -2.076004
#> LINC00552      -3.9038826       1.4276445      -2.6030898       -6.468321
#> LINC00910       0.1342525       0.7349029       0.6119231       -1.258868
#>           SAMf275eb859a39 SAM7f0d9cc7f001
#> LINC00696      -3.1070971      -2.8845109
#> LINC01126      -2.7200739      -1.6044029
#> LINC01460      -4.4856087      -3.3699377
#> LINC00552      -4.4856087      -1.7849752
#> LINC00910       0.7849203      -0.2000127

res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_data = imvigor210_eset,
                   id2 = "ID",
                   target= NULL,
                   group = "BOR_binary",
                   is_target_continuous = FALSE,
                   padj_cutoff = 1,
                   index = 1,
                   category = "gene",
                   signature_group =signature_collection[c(1:3)],
                   ProjectID = "IMvigor210",
                   palette_box = "paired1",
                   palette_corplot = "pheatmap",
                   palette_heatmap = 3,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>>  Processing signature: CD_8_T_effector"
```

<img src="man/figuresunnamed-chunk-26-1.png" width="100%" /><img src="man/figuresunnamed-chunk-26-2.png" width="100%" />

    #> [1] ">>>  Processing signature: DDR"

<img src="man/figuresunnamed-chunk-26-3.png" width="100%" /><img src="man/figuresunnamed-chunk-26-4.png" width="100%" />

    #> [1] ">>>  Processing signature: APM"

<img src="man/figuresunnamed-chunk-26-5.png" width="100%" /><img src="man/figuresunnamed-chunk-26-6.png" width="100%" />

## 5 Analyze signatures relevant to target signature.

### 5.1 Construct `pdata_group`as a phenotype data frame of patients.

``` r
# The data have been saved in `imvigor210_pdata`.
head(imvigor210_pdata)
#> # A tibble: 6 x 12
#>   ID    BOR   BOR_binary OS_days OS_status Mutation_Load Neo_antigen_Load
#>   <chr> <chr> <chr>      <chr>   <chr>     <chr>         <chr>           
#> 1 SAM0~ NA    NA         57.166~ 1         NA            NA              
#> 2 SAM0~ SD    NR         469.15~ 1         18            4.6862745099999~
#> 3 SAM0~ PD    NR         263.16~ 1         1             0.31372549      
#> 4 SAM0~ PD    NR         74.907~ 1         44            6.1960784310000~
#> 5 SAM0~ NA    NA         20.698~ 0         50            NA              
#> 6 SAM0~ SD    NR         136.01~ 1         2             1.4705882349999~
#> # ... with 5 more variables: CD_8_T_effector <dbl>, Immune_Checkpoint <dbl>,
#> #   Pan_F_TBRs <chr>, Mismatch_Repair <chr>, TumorPurity <dbl>
pdata_group<-as.data.frame(imvigor210_pdata[,c("ID","Pan_F_TBRs")])
pdata_group$Pan_F_TBRs<-scale(as.numeric(pdata_group$Pan_F_TBRs))
head(pdata_group)
#>                ID Pan_F_TBRs
#> 1 SAM00b9e5c52da9  0.2681507
#> 2 SAM0257bbbbd388  0.4618378
#> 3 SAM025b45c27e05 -0.8229645
#> 4 SAM032c642382a7 -0.4156874
#> 5 SAM04c589eb3fb3  9.0208975
#> 6 SAM0571f17f4045 -0.4462966
```

### 5.2 Analyze Pan-F-TBRs relevant signatures.

`Pan-F-TBRs` is an immune-exclusion signature published in Nature in
2018: S. Mariathasan et al., TGFbeta attenuates tumour response to PD-L1
blockade by contributing to exclusion of T cells. Nature 554, 544-548
(2018). the Pan-F-TBRs score obtained could be used to batch analyze and
visualize the interaction between this immune-exclusion signature and
other interested `signatures` , as well as its significant `cell
fractions`, and to further explore other biomarkers may impact
corresponding biological process.

``` r
res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_data = imvigor210_sig,
                   id2 = "ID",
                   target= "Pan_F_TBRs",
                   group = "group3",
                   is_target_continuous = TRUE,
                   padj_cutoff = 1,
                   index = 1,
                   category = "signature",
                   signature_group =sig_group[1:2],
                   ProjectID = "IMvigor210",
                   palette_box = "set2",
                   palette_corplot = "pheatmap",
                   palette_heatmap = 2,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>>  Processing signature: tumor_signature"
```

<img src="man/figuresunnamed-chunk-28-1.png" width="100%" /><img src="man/figuresunnamed-chunk-28-2.png" width="100%" /><img src="man/figuresunnamed-chunk-28-3.png" width="100%" />

    #> [1] ">>>  Processing signature: EMT"

<img src="man/figuresunnamed-chunk-28-4.png" width="100%" /><img src="man/figuresunnamed-chunk-28-5.png" width="100%" /><img src="man/figuresunnamed-chunk-28-6.png" width="100%" />

### 5.3 Analyze the immune cell infiltration lanscape related to Pan-F-TBRs.

``` r
names(sig_group)
#>  [1] "tumor_signature"        "EMT"                    "io_biomarkers"         
#>  [4] "immu_microenvironment"  "immu_suppression"       "immu_exclusion"        
#>  [7] "immu_exhaustion"        "TCR_BCR"                "tme_signatures1"       
#> [10] "tme_signatures2"        "Bcells"                 "Tcells"                
#> [13] "DCs"                    "Macrophages"            "Neutrophils"           
#> [16] "Monocytes"              "CAFs"                   "NK"                    
#> [19] "tme_cell_types"         "CIBERSORT"              "MCPcounter"            
#> [22] "EPIC"                   "xCell"                  "quanTIseq"             
#> [25] "ESTIMATE"               "IPS"                    "TIMER"                 
#> [28] "fatty_acid_metabolism"  "hypoxia_signature"      "cholesterol_metabolism"
#> [31] "Metabolism"             "hallmark"               "hallmark1"             
#> [34] "hallmark2"              "hallmark3"              "Rooney_et_al"          
#> [37] "Bindea_et_al"           "Li_et_al"               "Peng_et_al"
res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_data = imvigor210_sig,
                   id2 = "ID",
                   target= "Pan_F_TBRs",
                   group = "group3",
                   is_target_continuous = TRUE,
                   padj_cutoff = 1,
                   index = 1,
                   category = "signature",
                   signature_group =sig_group[20:24],
                   ProjectID = "IMvigor210",
                   palette_box = "jco",
                   palette_corplot = "pheatmap",
                   palette_heatmap = 3,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>>  Processing signature: CIBERSORT"
```

<img src="man/figuresunnamed-chunk-29-1.png" width="100%" /><img src="man/figuresunnamed-chunk-29-2.png" width="100%" />

    #> [1] ">>>  Processing signature: MCPcounter"

<img src="man/figuresunnamed-chunk-29-3.png" width="100%" /><img src="man/figuresunnamed-chunk-29-4.png" width="100%" /><img src="man/figuresunnamed-chunk-29-5.png" width="100%" />

    #> [1] ">>>  Processing signature: EPIC"

<img src="man/figuresunnamed-chunk-29-6.png" width="100%" /><img src="man/figuresunnamed-chunk-29-7.png" width="100%" /><img src="man/figuresunnamed-chunk-29-8.png" width="100%" />

    #> [1] ">>>  Processing signature: xCell"

<img src="man/figuresunnamed-chunk-29-9.png" width="100%" /><img src="man/figuresunnamed-chunk-29-10.png" width="100%" />

    #> [1] ">>>  Processing signature: quanTIseq"

<img src="man/figuresunnamed-chunk-29-11.png" width="100%" /><img src="man/figuresunnamed-chunk-29-12.png" width="100%" /><img src="man/figuresunnamed-chunk-29-13.png" width="100%" />

## 6 Serach for mutations related to signature

### 6.1 Mutation matrix using MAF data format.

MAF data was obtained from [UCSC Xena
website](https://api.gdc.cancer.gov/data/c06465a3-50e7-46f7-b2dd-7bd654ca206b)

``` r
help("make_mut_matrix")

maf_file<-"TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf"
# Each gene of every samples is categorized as binary variables(mutation or non-mutation) in the MAF data. 
mut_list<-make_mut_matrix(maf = maf_file,
                          isTCGA = T, 
                          category = "multi")
#> -Reading
#> -Validating
#> --Removed 2 duplicated variants
#> --Using only `Somatic` variants from Mutation_Status. Set useAll = TRUE to include everything.-Silent variants: 70966 
#> -Summarizing
#> --Possible FLAGS among top ten genes:
#>   TTN
#>   MUC16
#>   SYNE1
#>   FLG
#> -Processing clinical data
#> --Missing clinical data
#> -Finished in 12.1s elapsed (12.1s cpu) 
#>        Frame_Shift_Del        Frame_Shift_Ins           In_Frame_Del 
#>                  18418                   4461                    692 
#>           In_Frame_Ins      Missense_Mutation      Nonsense_Mutation 
#>                    268                 109668                   6011 
#>       Nonstop_Mutation            Splice_Site Translation_Start_Site 
#>                    107                   2445                    106 
#>    DEL    INS    SNP 
#>  19387   4900 117889
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
#> Aggregation function missing: defaulting to length
# If "multi" is set in above "category" parameter,four data frames will be returned, which evaluate all the mutations of every gene or estimate only SNP, indel and frameshift as follow:
names(mut_list)
#> [1] "all"        "snp"        "indel"      "frameshift"
mut_list$all[1:5,1:10]
#>              A1BG A1CF A2M A2ML1 A3GALT2 A4GALT A4GNT AAAS AACS AADAC
#> TCGA-3M-AB46    1    1   0     0       0      0     0    0    0     0
#> TCGA-3M-AB47    0    0   0     0       0      0     0    0    0     0
#> TCGA-B7-5816    0    0   1     0       0      0     0    1    0     0
#> TCGA-B7-5818    0    0   0     0       0      0     0    0    0     0
#> TCGA-B7-A5TI    0    0   0     0       0      0     0    0    0     0
```

### 6.2 sig\_group() function for batch visualization of response associated signatures.

``` r
# Transform the mutation data into categorical variable matrix.----------
mut<-mut_list$snp
max(mut)
#> [1] 22
# If the maximum of mutation counts of a single gene of a person in is over 4, the mutation data would be standardized according to following principles.
# mut[mut>=3&mut<=5]<-3
# mut[mut>5]<-4
##########################
res<-find_mutations(mutation_matrix = mut, 
                    signature_matrix = tcga_stad_sig,
                    id_signature_matrix = "ID",
                    signature = "CD_8_T_effector",
                    min_mut_freq = 0.01,
                    plot = TRUE,
                    method = "multi",
                    save_path = paste0("CD_8_T_effector-relevant-mutations"),
                    palette = "jama",
                    show_plot = T)
#> [1] ">>>> Result of Cuzick Test"
#>             p.value  names statistic adjust_pvalue
#> PIK3CA 3.148160e-09 PIK3CA  5.923680  1.574080e-06
#> SPEG   9.070928e-05   SPEG  3.914187  2.267732e-02
#> TCHH   4.409469e-04   TCHH  3.514281  5.740100e-02
#> PLXNA4 5.420662e-04 PLXNA4  3.459059  5.740100e-02
#> ARID1A 5.805905e-04 ARID1A  3.440523  5.740100e-02
#> WDFY3  6.888120e-04  WDFY3  3.393994  5.740100e-02
#> GTF3C1 8.120095e-04 GTF3C1  3.348668  5.800068e-02
#> DMD    1.675467e-03    DMD  3.142439  6.972915e-02
#> CR1    1.775997e-03    CR1  3.125340  6.972915e-02
#> EP300  2.042083e-03  EP300  3.084043  6.972915e-02
#> [1] ">>>> Result of Wilcoxon test"
#>             p.value  names statistic adjust_pvalue
#> PIK3CA 1.921035e-10 PIK3CA      4125  9.605174e-08
#> TCHH   1.961642e-05   TCHH      3312  4.904106e-03
#> SPEG   3.532750e-05   SPEG      1947  5.887916e-03
#> LRP1   7.511741e-05   LRP1      2649  9.389676e-03
#> WDFY3  1.257659e-04  WDFY3      2964  1.257659e-02
#> ARID1A 2.468609e-04 ARID1A      4878  2.057174e-02
#> PLXNA4 4.215809e-04 PLXNA4      3638  3.011292e-02
#> ANK3   6.399572e-04   ANK3      4446  3.933979e-02
#> DMD    7.364591e-04    DMD      5311  3.933979e-02
#> PLEC   8.026240e-04   PLEC      5562  3.933979e-02
```

<img src="man/figuresunnamed-chunk-31-1.png" width="100%" /><img src="man/figuresunnamed-chunk-31-2.png" width="100%" />

## 7 Other methods for batch analysis

IOBR provide multiple batch analysis functions for statical analysis,
data operation and transformation. Batch survival analysis, for example,
could batch analyze the best cutoff value for subsequent batch survival
analysis, varied batch statistical tests, etc. \#\#\# 7.1 subgroup
survival analyses

``` r
help("subgroup_survival")
##loading data and filter NA
data(subgroup_data)
input <- subgroup_data %>% 
   filter(time > 0) %>% 
   filter(!is.na(status)) %>% 
   filter(!is.na(AJCC_stage))
dim(input)
#> [1] 1333    7
##for binary variable
data1 <- subgroup_survival(pdata = input,
                           time ="time", 
                           status = "status",
                           variable = c("ProjectID", "AJCC_stage"), 
                           object ="score_binary" )
data1
#>                         P     HR CI_low_0.95 CI_up_0.95
#> ProjectID_Dataset1 0.0000 3.0901      1.8648     5.1205
#> ProjectID_Dataset2 0.0023 2.1049      1.3057     3.3931
#> ProjectID_Dataset3 0.1397 2.4431      0.7467     7.9940
#> ProjectID_Dataset4 0.0024 2.0312      1.2854     3.2097
#> ProjectID_Dataset5 0.0000 4.6375      2.3759     9.0520
#> AJCC_stage_1       0.2533 1.7307      0.6753     4.4354
#> AJCC_stage_2       0.0000 2.9450      1.7858     4.8569
#> AJCC_stage_3       0.0000 2.2821      1.5551     3.3490
#> AJCC_stage_4       0.0510 1.5941      0.9979     2.5466
```

## References

DQ Zeng, ZL Ye, RF Shen, Y Xiong, JN Wu, WJ Qiu WJ Liao.

IOBR: Comprehensive analysis of tumor microenvironment and signatures
for immuno-oncology.

## Contact

E-mail any questions to <dongqiangzeng0808@gmail.com>
