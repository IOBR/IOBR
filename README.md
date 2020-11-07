
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

The dependencies including `tibble`, `survival`, `survminer`, `limma`,
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
#> 1           0.000000000            0.00383182              0.00
#> 2           0.005580183            0.03892844              0.00
#> 3           0.033511519            0.08712756              0.05
#> 4           0.000000000            0.01069773              0.00
#> 5           0.000000000            0.06181853              0.04
#> 6           0.000000000            0.04979512              0.00
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

### 3.2.1 Calculate TME associated signatures-(through PCA method)

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

### 3.2.2 Calculate TME associated signatures-(through ssGSEA method)

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

### 3.2.3 Calculate metabolism related signatures

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

### 3.2.4 Calculate all collected signature scores (integrating three methods: PCA, ssGSEA and z-score)

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
ssGSEA method is recommended for these signature estimation. It may take
a while for the code to run, if input a large sample data or to
calculate a large number of signatures.

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

### 4.1 这里利用IMvigor210的免疫治疗数据来批量分析与phenotype相关的signature

``` r
## 加载IOBR中已经计算好的signature和 Cell fraction
head(imvigor210_sig)
#>                ID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT
#> 1 SAMf2ce197162ce              0.02900304              0.045657003
#> 2 SAM698d8d76b934              0.08119867              0.001266419
#> 3 SAMc1b27bc16435              0.01237173              0.000000000
#> 4 SAM85e41e7f33f9              0.00000000              0.000000000
#> 5 SAMf275eb859a39              0.00000000              0.000000000
#> 6 SAM7f0d9cc7f001              0.05819870              0.128559340
#>   Plasma_cells_CIBERSORT T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1            0.000000000                     0                   0.1316134
#> 2            0.000000000                     0                   0.0000000
#> 3            0.001273053                     0                   0.1693195
#> 4            0.001186289                     0                   0.0000000
#> 5            0.009493367                     0                   0.0000000
#> 6            0.000000000                     0                   0.0000000
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                            0.2390998                             0.04564597
#> 2                            0.1963620                             0.03865429
#> 3                            0.3168480                             0.08919590
#> 4                            0.3188597                             0.04426192
#> 5                            0.3923644                             0.04485934
#> 6                            0.3715513                             0.03490259
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                          0.00000000                                    0
#> 2                          0.00822820                                    0
#> 3                          0.00000000                                    0
#> 4                          0.01379417                                    0
#> 5                          0.00000000                                    0
#> 6                          0.02685719                                    0
#>   T_cells_gamma_delta_CIBERSORT NK_cells_resting_CIBERSORT
#> 1                             0                 0.06136277
#> 2                             0                 0.03325457
#> 3                             0                 0.08718412
#> 4                             0                 0.02929384
#> 5                             0                 0.11512150
#> 6                             0                 0.08005595
#>   NK_cells_activated_CIBERSORT Monocytes_CIBERSORT Macrophages_M0_CIBERSORT
#> 1                   0.00000000         0.031217612               0.29544635
#> 2                   0.00000000         0.043832992               0.18030697
#> 3                   0.00000000         0.013799170               0.15507559
#> 4                   0.03507366         0.002027608               0.16842597
#> 5                   0.00000000         0.036451256               0.09402453
#> 6                   0.00000000         0.036758356               0.07425533
#>   Macrophages_M1_CIBERSORT Macrophages_M2_CIBERSORT
#> 1              0.004193591               0.06828464
#> 2              0.040512883               0.26862850
#> 3              0.008568925               0.04674582
#> 4              0.026754131               0.23602698
#> 5              0.009412639               0.10234745
#> 6              0.002842425               0.07806874
#>   Dendritic_cells_resting_CIBERSORT Dendritic_cells_activated_CIBERSORT
#> 1                        0.00000000                          0.03149558
#> 2                        0.00000000                          0.03869560
#> 3                        0.00000000                          0.02492113
#> 4                        0.00000000                          0.03770440
#> 5                        0.02896006                          0.05515404
#> 6                        0.00000000                          0.01462629
#>   Mast_cells_resting_CIBERSORT Mast_cells_activated_CIBERSORT
#> 1                            0                    0.006246756
#> 2                            0                    0.069058871
#> 3                            0                    0.073027769
#> 4                            0                    0.070689539
#> 5                            0                    0.062317814
#> 6                            0                    0.048727474
#>   Eosinophils_CIBERSORT Neutrophils_CIBERSORT T_cells_MCPcounter
#> 1           0.004733079           0.006000408           1.870733
#> 2           0.000000000           0.000000000           0.712941
#> 3           0.001669244           0.000000000           2.472843
#> 4           0.011702681           0.004199108           1.387958
#> 5           0.000000000           0.049493574           2.290112
#> 6           0.016453814           0.028142501           3.083601
#>   CD8_T_cells_MCPcounter Cytotoxic_lymphocytes_MCPcounter NK_cells_MCPcounter
#> 1             -1.9563500                       -0.0718495          -3.6729485
#> 2             -2.8775404                       -0.1231299          -4.4408517
#> 3             -0.6555572                        1.8123831          -2.1480694
#> 4             -1.3390380                        2.2554275          -0.6324088
#> 5             -1.9495558                        0.6824556          -2.4856234
#> 6             -1.4439383                        1.5713859          -1.5777600
#>   B_lineage_MCPcounter Monocytic_lineage_MCPcounter
#> 1            2.1027820                     3.720460
#> 2            2.3395481                     4.118058
#> 3            1.5324034                     3.074732
#> 4           -1.0506905                     4.886917
#> 5            0.6249503                     4.303551
#> 6            4.4192107                     4.077619
#>   Myeloid_dendritic_cells_MCPcounter Neutrophils_MCPcounter
#> 1                         -0.2275544              2.5475104
#> 2                         -0.8044946              0.9494973
#> 3                          0.4080493              1.8840660
#> 4                         -1.2508100              1.5909430
#> 5                          1.9508177              3.7048758
#> 6                          1.1853367              3.2258626
#>   Endothelial_cells_MCPcounter Fibroblasts_MCPcounter  aDC_xCell
#> 1                     3.167174               8.989274 0.08217075
#> 2                     2.441914               8.055582 0.17152703
#> 3                     3.312655               6.926645 0.06062884
#> 4                     2.417368               8.999718 0.21005269
#> 5                     3.191219               6.048701 0.10050575
#> 6                     4.105552               8.776603 0.09805099
#>   Adipocytes_xCell Astrocytes_xCell B-cells_xCell Basophils_xCell
#> 1     1.795976e-20     8.698935e-02   0.002803556     0.017094527
#> 2     0.000000e+00     6.883639e-02   0.038570594     0.041511626
#> 3     3.854502e-20     0.000000e+00   0.000000000     0.029538839
#> 4     2.426480e-03     8.630578e-02   0.013291929     0.019226765
#> 5     0.000000e+00     2.022866e-18   0.005483038     0.033532301
#> 6     7.681178e-03     7.973259e-02   0.147953568     0.002545997
#>   CD4+_memory_T-cells_xCell CD4+_naive_T-cells_xCell CD4+_T-cells_xCell
#> 1               0.000000000             2.688567e-19       1.245304e-19
#> 2               0.000000000             9.508148e-19       1.890989e-19
#> 3               0.013522259             8.072307e-03       1.001745e-02
#> 4               0.006429428             0.000000e+00       7.361980e-03
#> 5               0.007590415             0.000000e+00       1.720743e-03
#> 6               0.005302120             0.000000e+00       1.171757e-02
#>   CD4+_Tcm_xCell CD4+_Tem_xCell CD8+_naive_T-cells_xCell CD8+_T-cells_xCell
#> 1   6.074539e-03   0.000000e+00             2.911591e-20       0.0000000000
#> 2   0.000000e+00   0.000000e+00             1.046325e-02       0.0002703006
#> 3   2.441868e-02   9.240324e-03             2.143943e-02       0.0451232846
#> 4   4.862074e-18   0.000000e+00             7.379417e-03       0.0000000000
#> 5   2.218192e-02   2.058918e-02             2.968490e-03       0.0128073490
#> 6   2.160600e-02   4.864743e-19             3.104862e-03       0.0127521418
#>   CD8+_Tcm_xCell CD8+_Tem_xCell  cDC_xCell Chondrocytes_xCell
#> 1    0.000000000   7.908344e-19 0.03710992       5.104980e-02
#> 2    0.000000000   2.672083e-19 0.05052769       1.164654e-02
#> 3    0.016235059   0.000000e+00 0.03061615       1.904156e-02
#> 4    0.005758151   4.060382e-03 0.08753767       5.803518e-19
#> 5    0.001609873   0.000000e+00 0.11545946       0.000000e+00
#> 6    0.000000000   2.807802e-18 0.12036488       6.069508e-02
#>   Class-switched_memory_B-cells_xCell  CLP_xCell    CMP_xCell     DC_xCell
#> 1                         0.020229582 0.07265988 0.000000e+00 1.926217e-03
#> 2                         0.031417552 0.14750072 7.418611e-03 1.330796e-21
#> 3                         0.020284717 0.15420987 3.874713e-20 4.076250e-03
#> 4                         0.002536295 0.20228345 0.000000e+00 2.716205e-03
#> 5                         0.012734756 0.14844918 0.000000e+00 8.559884e-03
#> 6                         0.085976326 0.02517814 5.229075e-03 1.179110e-02
#>   Endothelial_cells_xCell Eosinophils_xCell Epithelial_cells_xCell
#> 1             0.021657797      0.000000e+00             0.09582442
#> 2             0.014927749      0.000000e+00             0.01064872
#> 3             0.056115417      3.755652e-03             0.15065290
#> 4             0.002274537      1.880995e-19             0.18344048
#> 5             0.034813361      5.504984e-19             0.15395672
#> 6             0.075649335      8.153991e-03             0.08167594
#>   Erythrocytes_xCell Fibroblasts_xCell    GMP_xCell Hepatocytes_xCell
#> 1       0.000000e+00      0.1194671134 3.154064e-18      1.660225e-03
#> 2       7.523285e-05      0.0378708024 7.781268e-19      1.305920e-19
#> 3       0.000000e+00      0.0093959116 1.822108e-18      0.000000e+00
#> 4       4.429546e-19      0.0006427847 0.000000e+00      0.000000e+00
#> 5       0.000000e+00      0.0311242372 0.000000e+00      0.000000e+00
#> 6       9.862431e-20      0.2685284814 4.947730e-02      2.485180e-03
#>    HSC_xCell  iDC_xCell Keratinocytes_xCell ly_Endothelial_cells_xCell
#> 1 0.12490414 0.09424564         0.040804530               3.910666e-04
#> 2 0.11189960 0.02425641         0.001382772               1.202748e-19
#> 3 0.11057710 0.10020700         0.028566886               1.479536e-02
#> 4 0.03298200 0.09652857         0.130938529               0.000000e+00
#> 5 0.07416533 0.09746552         0.061229272               3.587088e-03
#> 6 0.30364157 0.24917819         0.003136851               2.854756e-02
#>   Macrophages_xCell Macrophages_M1_xCell Macrophages_M2_xCell Mast_cells_xCell
#> 1       0.007806223          0.008192630         0.0073138160     0.0009551974
#> 2       0.033458121          0.029928398         0.0192428116     0.0169771211
#> 3       0.003984817          0.003872689         0.0003743472     0.0020520421
#> 4       0.067760302          0.048185487         0.0155254355     0.0016279814
#> 5       0.009127548          0.010451845         0.0011195475     0.0054138335
#> 6       0.006091868          0.012597714         0.0049569141     0.0067563813
#>   Megakaryocytes_xCell Melanocytes_xCell Memory_B-cells_xCell   MEP_xCell
#> 1         3.519704e-03       0.011764110         2.879417e-19 0.019709206
#> 2         8.908976e-20       0.002048811         4.650534e-03 0.102993908
#> 3         5.915070e-03       0.003331928         3.604285e-19 0.031534279
#> 4         1.435880e-03       0.002984078         3.237338e-03 0.009360597
#> 5         6.760601e-03       0.006678787         4.933745e-18 0.003489972
#> 6         1.881286e-02       0.000000000         4.983961e-02 0.015576429
#>   Mesangial_cells_xCell Monocytes_xCell    MPP_xCell    MSC_xCell
#> 1          2.310541e-02     0.027629612 1.331530e-20 2.300920e-01
#> 2          4.426161e-19     0.018573114 5.877777e-19 2.215625e-01
#> 3          1.358396e-02     0.009496829 1.460458e-17 8.510970e-02
#> 4          1.210097e-02     0.046508684 3.144816e-18 1.902891e-01
#> 5          1.235094e-02     0.038569786 1.307121e-17 3.836456e-02
#> 6          5.140956e-02     0.050122516 1.404913e-17 7.655075e-18
#>   mv_Endothelial_cells_xCell Myocytes_xCell naive_B-cells_xCell Neurons_xCell
#> 1                0.016030783    0.000000000        1.339336e-18   0.004849936
#> 2                0.025087822    0.018489500        7.137479e-03   0.021041531
#> 3                0.045807584    0.014913492        1.006773e-04   0.002848984
#> 4                0.003646426    0.008781728        9.528622e-04   0.001618478
#> 5                0.026957842    0.006106227        4.303641e-04   0.005304488
#> 6                0.053713290    0.025142880        2.064891e-02   0.037681121
#>   Neutrophils_xCell NK_cells_xCell   NKT_xCell Osteoblast_xCell    pDC_xCell
#> 1      0.000000e+00   1.179256e-18 0.001517122     6.681496e-18 6.052418e-20
#> 2      0.000000e+00   2.655550e-19 0.000000000     1.715928e-02 0.000000e+00
#> 3      3.006289e-19   3.507185e-18 0.003992388     0.000000e+00 0.000000e+00
#> 4      2.779902e-18   0.000000e+00 0.014943564     7.394898e-18 1.187306e-02
#> 5      3.094353e-03   0.000000e+00 0.006603545     2.361239e-03 9.326823e-19
#> 6      3.513022e-03   0.000000e+00 0.038717814     0.000000e+00 0.000000e+00
#>   Pericytes_xCell Plasma_cells_xCell Platelets_xCell Preadipocytes_xCell
#> 1    1.109950e-01        0.010374345     0.000000000        8.847290e-02
#> 2    0.000000e+00        0.003122261     0.000000000        0.000000e+00
#> 3    2.795545e-18        0.009727119     0.010063869        3.230092e-19
#> 4    3.104109e-03        0.004337165     0.000000000        1.463230e-19
#> 5    2.975865e-18        0.011277349     0.002131944        9.165871e-02
#> 6    5.411425e-02        0.000000000     0.001336909        1.268837e-01
#>   pro_B-cells_xCell Sebocytes_xCell Skeletal_muscle_xCell Smooth_muscle_xCell
#> 1      3.332585e-18    0.0108548858           0.005419181           0.2159813
#> 2      2.135625e-02    0.0002007755           0.001808531           0.2234351
#> 3      0.000000e+00    0.0200361762           0.002334263           0.2401171
#> 4      4.391544e-03    0.0458140515           0.000000000           0.2880682
#> 5      0.000000e+00    0.0095461735           0.000000000           0.1118417
#> 6      0.000000e+00    0.0034465975           0.004472427           0.1093343
#>   Tgd_cells_xCell Th1_cells_xCell Th2_cells_xCell  Tregs_xCell
#> 1    0.000000e+00     0.016938881     0.081510426 1.770388e-02
#> 2    0.000000e+00     0.104980746     0.112521326 6.730478e-21
#> 3    1.907781e-02     0.013719699     0.109963333 3.089447e-02
#> 4    2.025314e-02     0.060176938     0.101678029 2.681846e-02
#> 5    1.967145e-03     0.007904472     0.035931684 3.191049e-02
#> 6    1.416379e-19     0.001326582     0.004012456 9.576664e-02
#>   ImmuneScore_xCell StromaScore_xCell MicroenvironmentScore_xCell  Bcells_EPIC
#> 1        0.02741387       0.070562455                  0.09797633 3.384960e-02
#> 2        0.07189950       0.026399276                  0.09829878 3.302083e-02
#> 3        0.05233755       0.032755664                  0.08509322 2.195392e-02
#> 4        0.09284472       0.002671901                  0.09551662 9.869922e-09
#> 5        0.05651769       0.032968799                  0.08948649 1.638629e-02
#> 6        0.17256811       0.175929497                  0.34849761 5.236272e-02
#>    CAFs_EPIC CD4_Tcells_EPIC CD8_Tcells_EPIC Endothelial_EPIC Macrophages_EPIC
#> 1 0.02668757       0.2808177    5.658129e-07        0.1382526      0.009638566
#> 2 0.02448600       0.2463193    7.709019e-08        0.1052977      0.010340225
#> 3 0.02161894       0.3399107    3.421037e-02        0.1615781      0.007191021
#> 4 0.02742702       0.3358012    9.645048e-03        0.1055330      0.012526942
#> 5 0.01747081       0.3404707    9.465339e-03        0.1504026      0.008216401
#> 6 0.02350852       0.3232805    2.396414e-02        0.1723092      0.008391381
#>   NKcells_EPIC otherCells_EPIC StromalScore_estimate ImmuneScore_estimate
#> 1 4.518694e-10       0.5107534             602.68615             91.14194
#> 2 3.385607e-10       0.5805358             -55.76716           -127.34739
#> 3 5.277051e-11       0.4135369            -471.64985             93.37235
#> 4 8.268158e-09       0.5090667             311.62047            709.67246
#> 5 6.040282e-11       0.4575878            -791.47587            432.65370
#> 6 1.206618e-09       0.3961836             722.35072            794.50766
#>   ESTIMATEScore_estimate TumorPurity_estimate B_cells_quantiseq
#> 1               693.8281            0.7604224       0.025752032
#> 2              -183.1146            0.8374980       0.023543678
#> 3              -378.2775            0.8528064       0.015612822
#> 4              1021.2929            0.7283393       0.007929979
#> 5              -358.8222            0.8513116       0.021638949
#> 6              1516.8584            0.6766121       0.050482100
#>   Macrophages_M1_quantiseq Macrophages_M2_quantiseq Monocytes_quantiseq
#> 1              0.005987321               0.04358373          0.00000000
#> 2              0.016734895               0.01751459          0.04092005
#> 3              0.004295221               0.05721011          0.00000000
#> 4              0.010570050               0.02185907          0.00000000
#> 5              0.000000000               0.04396653          0.00000000
#> 6              0.010925619               0.05933714          0.00000000
#>   Neutrophils_quantiseq NK_cells_quantiseq T_cells_CD4_quantiseq
#> 1            0.12800909        0.007598058                     0
#> 2            0.08246606        0.011144575                     0
#> 3            0.12017542        0.005196617                     0
#> 4            0.23935679        0.003557837                     0
#> 5            0.18773918        0.006384031                     0
#> 6            0.18331350        0.017138731                     0
#>   T_cells_CD8_quantiseq Tregs_quantiseq Dendritic_cells_quantiseq
#> 1           0.000000000      0.02834501                         0
#> 2           0.001519011      0.02485400                         0
#> 3           0.024641275      0.03598988                         0
#> 4           0.004893350      0.01225338                         0
#> 5           0.002638232      0.02681215                         0
#> 6           0.011865695      0.04773090                         0
#>   Other_quantiseq B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER
#> 1       0.7607248   0.09189314        0.1388286        0.1693827
#> 2       0.7813031   0.08170135        0.1377691        0.1460446
#> 3       0.7368787   0.07628695        0.1447232        0.1899296
#> 4       0.6995795   0.05711371        0.1378148        0.1859518
#> 5       0.7108209   0.07380974        0.1535897        0.1609782
#> 6       0.6192063   0.09261939        0.1575034        0.1844467
#>   Neutrophil_TIMER Macrophage_TIMER  DC_TIMER CD_8_T_effector        DDR
#> 1        0.1305428       0.06254971 0.4897297      -1.0722292  0.3858460
#> 2        0.1259930       0.05645570 0.4854951       0.2243254  1.6671380
#> 3        0.1227739       0.04793405 0.4859183       0.7153629  0.8845539
#> 4        0.1393475       0.06455816 0.5032895       1.3451381  0.5424434
#> 5        0.1394671       0.04389780 0.4999189       0.2959889 -1.0584244
#> 6        0.1409084       0.05182806 0.5059829      -0.4314508 -1.2688993
#>           APM Immune_Checkpoint CellCycle_Reg Pan_F_TBRs    Histones       EMT1
#> 1 -0.31657936       -0.33332402     1.0050427  1.1188464  0.05824973 -0.1713622
#> 2  0.24790540        0.15769947     0.9255255  0.1475944  0.48329306  0.5861912
#> 3 -0.30649968        0.20202295     0.2059704 -0.4788580  0.48145159 -0.1242950
#> 4  0.20945501        0.81441444     1.0586794  0.8800992  0.42068365  0.1521361
#> 5  0.04589724        0.06139678    -0.9224057 -1.3172525 -0.68906843 -0.5485783
#> 6 -0.29499364        0.58860675     0.9237428  1.2768864 -0.79989487  0.4065044
#>         EMT2       EMT3  WNT_target FGFR3_related Cell_cycle Mismatch_Repair
#> 1  1.0017703  0.8066020  0.12438920    0.04381166  0.4905936       0.3741222
#> 2  0.1538176 -0.3632004  0.26228167   -2.07295906  1.7786461       0.9299782
#> 3 -0.7935536 -0.8254196  0.98657468    0.81158327  0.7831659       0.3465392
#> 4  1.0865484  0.3370684 -0.06447579    0.24313357  1.4713395       0.5430490
#> 5 -1.2236291 -1.1235620  0.09679941    1.08881903 -1.6588511      -0.6530710
#> 6  0.6279988  1.0736917 -0.16371050   -1.05425967 -1.1760630      -0.4845738
#>   Homologous_recombination Nucleotide_excision_repair DNA_replication
#> 1              -0.02085097                  0.1960518       0.3875528
#> 2               1.02598400                  0.9330273       1.1988024
#> 3               0.52716662                  0.3094811       0.4161723
#> 4               0.29320358                  0.6266390       0.7754689
#> 5              -0.58918968                 -0.6251043      -0.8703067
#> 6              -0.84655286                 -0.4838426      -0.5397244
#>   Base_excision_repair TMEscoreA_CIR TMEscoreB_CIR CD8_Rooney_et_al
#> 1          -0.20569323    -1.1182793    3.01460322       -0.4386698
#> 2           0.78340178     0.4815436    0.08144524       -0.3654156
#> 3           0.02887199     0.3591707   -2.74600243        0.7533909
#> 4           0.64013909     2.3769708    0.46138655        0.9994595
#> 5          -0.43261489     0.6971289   -3.48772483        0.2246819
#> 6          -0.52015304    -0.1556606    4.96154160        0.3018941
#>   B_cells_Rooney_et_al Treg_Rooney_et_al Macrophages_Rooney_et_al
#> 1            0.5917744        -0.1189047                0.2890008
#> 2            0.1060844        -0.2121649                0.7929304
#> 3            0.5455201         0.8719461               -0.5310392
#> 4           -1.3197796         0.1444663                1.0994335
#> 5           -0.3246520         0.9458732               -0.6216495
#> 6            1.7673803         0.5066192               -0.1783162
#>   Neutrophils_Rooney_et_al NK_cells_Rooney_et_al pDCs_Rooney_et_al
#> 1               0.36238971           -0.05560187       -0.87703316
#> 2              -0.39107400           -0.49866228        0.12232472
#> 3              -0.02897612            0.41407196       -0.82715791
#> 4              -0.22717038            0.37809679       -0.09772715
#> 5               1.08868514            1.08163883        0.25423313
#> 6               1.25656366            0.32839630        0.12555102
#>   MHC_Class_I_Rooney_et_al Co_stimulation_APC_Rooney_et_al
#> 1               -0.2761937                      0.32769386
#> 2                0.2072330                     -0.09885921
#> 3               -0.1869448                      0.53378329
#> 4                0.2497558                      0.69846218
#> 5                0.1125277                     -0.01078629
#> 6               -0.1603536                      0.10115458
#>   Co_stimulation_T_cell_Rooney_et_al Co_inhibition_APC_Rooney_et_al
#> 1                         -0.2568277                    -0.30291222
#> 2                         -0.2474173                     0.02019607
#> 3                          0.5568222                    -0.32529577
#> 4                          0.4764837                     0.33650442
#> 5                          0.6101593                     0.15657956
#> 6                          1.0443919                     0.09103965
#>   Co_inhibition_T_cell_Rooney_et_al Type_I_IFN_Reponse_Rooney_et_al
#> 1                      -0.167345709                      -0.9870562
#> 2                      -0.006800161                       0.8027931
#> 3                       0.315978394                      -0.8672649
#> 4                       0.727850180                       1.4421641
#> 5                       0.309432962                      -0.3169421
#> 6                       0.648588050                      -0.2020362
#>   Type_II_IFN_Reponse_Rooney_et_al Cytolytic_Activity_Rooney_et_al MHC_Class_I
#> 1                        0.1893699                      -0.2330544 -0.31683787
#> 2                       -0.4600845                      -0.2983182  0.25802701
#> 3                        0.6939777                       0.4105303 -0.28587691
#> 4                       -0.6592874                       0.8080796  0.14091204
#> 5                        0.3120635                       0.2152258 -0.03720336
#> 6                        0.6193444                       0.2745450 -0.31556165
#>   T_cell_inflamed_GEP_Ayers_et_al IFNG_signature_Ayers_et_al MDSC_Wang_et_al
#> 1                      -0.6213280                 -1.0999504       0.2706644
#> 2                       0.1215910                  0.7178762       0.4646121
#> 3                       0.7220965                  0.3462288      -0.4322426
#> 4                       1.3642053                  0.8903432       1.5129888
#> 5                       0.5742895                  0.5343390       1.1386674
#> 6                       0.4517527                 -0.9223018       1.2252798
#>        GPAGs       PPAGs HLA_signature_gene B_cells_Danaher_et_al
#> 1  0.8023706  0.02525884         -1.1646488            0.62642044
#> 2  0.5048880 -2.53856304          2.0011168            0.69160461
#> 3 -0.7391467  0.60809847         -0.8819123            0.03374711
#> 4 -0.2901284 -0.15910878          1.6010495           -2.18227444
#> 5 -0.5961889  0.84031186         -0.1447313           -0.50659844
#> 6  2.0598534  0.15263183          5.5179216            1.84347678
#>   CD8_T_cells_Danaher_et_al Cytotoxic_cells_Danaher_et_al DC_Danaher_et_al
#> 1               -0.30898559                    -0.6175153       0.12941089
#> 2               -0.39666534                    -0.6716328       0.43642046
#> 3                0.61860187                     0.8794285       0.08233783
#> 4                0.23013375                     1.5398549       0.53489101
#> 5               -0.01219789                     0.3070409       0.20308363
#> 6                0.16534419                     0.5419669       0.50882325
#>   Exhausted_CD8_Danaher_et_al Macrophages_Danaher_et_al
#> 1                 -0.46450627               -0.06240368
#> 2                 -0.14056322                0.35841869
#> 3                  0.41260250               -0.53250090
#> 4                  0.56766814                0.60176187
#> 5                 -0.04566691               -0.30727550
#> 6                  0.03989999                0.05927750
#>   Mast_cells_Danaher_et_al Neutrophils_Danaher_et_al
#> 1                0.4944537                 0.4154951
#> 2                1.2860528                 0.1162757
#> 3                0.2951682                -0.3825528
#> 4               -0.5989872                 0.6707447
#> 5                0.4684027                 1.2749127
#> 6                1.6917196                 1.1322373
#>   NK_CD56dim_cells_Danaher_et_al NK_cells_Danaher_et_al T_cells_Danaher_et_al
#> 1                    -1.45780056            -0.01693590            0.25095580
#> 2                    -0.71045755             0.56902179           -0.46621663
#> 3                     0.17187515             0.08060774            0.72157856
#> 4                     1.12327537             0.01194487            0.08028545
#> 5                     0.46941724            -0.85729445            0.49355053
#> 6                     0.07932662             0.52064641            0.80886732
#>   TIP_Release_of_cancer_cell_antigens TIP_Cancer_antigen_presentation_1
#> 1                         -0.03180457                       0.006100139
#> 2                          0.12427722                       0.039967684
#> 3                         -0.28634458                      -0.633002114
#> 4                          0.85237066                       0.370092085
#> 5                          0.23436379                       0.803613323
#> 6                         -0.40569443                       0.415559841
#>   TIP_Priming_and_activation_1 TIP_Priming_and_activation_2
#> 1                  -0.06918646                  -0.27343848
#> 2                  -0.15974460                  -0.01533992
#> 3                   0.90725008                   0.51045276
#> 4                   1.03996328                   0.78517416
#> 5                   0.32013580                   0.19174846
#> 6                   1.54191487                   0.84818324
#>   TIP_Trafficking_of_immune_cells_to_tumors
#> 1                                0.45076932
#> 2                                0.46842079
#> 3                                0.42178894
#> 4                                1.60653821
#> 5                                0.08017683
#> 6                                1.69599231
#>   TIP_Infiltration_of_immune_cells_into_tumors_1
#> 1                                     -0.1243353
#> 2                                      0.1925681
#> 3                                     -0.1692498
#> 4                                      0.2188219
#> 5                                     -0.1463984
#> 6                                      0.6623051
#>   TIP_Infiltration_of_immune_cells_into_tumors_2
#> 1                                     0.09630952
#> 2                                     0.08531847
#> 3                                    -0.17611738
#> 4                                     0.26877242
#> 5                                    -0.06434480
#> 6                                     0.71938794
#>   TIP_Recognition_of_cancer_cells_by_T_cells_1
#> 1                                    0.1638849
#> 2                                   -0.5843075
#> 3                                    0.7920848
#> 4                                    0.5226644
#> 5                                    0.7743662
#> 6                                    1.2578433
#>   TIP_Recognition_of_cancer_cells_by_T_cells_2 TIP_Killing_of_cancer_cells_1
#> 1                                   -0.2483571                    -0.8409183
#> 2                                    0.5022748                     0.4224869
#> 3                                    0.2116241                     0.6468164
#> 4                                    0.1118034                     1.1716413
#> 5                                    0.1583868                     0.3121693
#> 6                                    0.7609604                    -0.3782959
#>   TIP_Killing_of_cancer_cells_2  TLS_Nature TMEscoreA_plus TMEscoreB_plus
#> 1                    -0.3580841  0.58651104     -1.4193314     0.61673215
#> 2                     0.7402292 -1.51237451      0.4901072     0.05481347
#> 3                     0.1381383 -0.04055625      0.4543883    -0.74495838
#> 4                     0.5950119  0.81763800      1.9594211     0.86922803
#> 5                     0.5994629  0.65067960      0.6766036    -0.54467809
#> 6                     0.9392708  0.82698403     -0.8219949     1.39139794
#>   B_cells_Bindea_et_al T_cells_Bindea_et_al T_helper_cells_Bindea_et_al
#> 1            0.6566410          -0.09104612                 -0.25076513
#> 2            0.9432744          -0.70655232                  0.05378525
#> 3           -0.1893713           0.82050459                  0.36356151
#> 4           -2.3718216           0.14334590                  0.20156158
#> 5           -0.6573705           0.67957990                  0.33638103
#> 6            2.4102216           1.21167267                  0.92174680
#>   Tcm_Bindea_et_al Tem_Bindea_et_al Th1_cells_Bindea_et_al
#> 1        0.6661151       0.42436434             -0.6630341
#> 2       -2.0103695      -1.59489068              0.5632527
#> 3        0.1000980       0.03757222              0.1428366
#> 4        0.5673032       0.65788186              1.3686936
#> 5        0.6637185       0.78727544              0.3228507
#> 6        1.0593604       0.88609939             -0.1370525
#>   Th2_cells_Bindea_et_al TFH_Bindea_et_al Th17_cells_Bindea_et_al
#> 1              0.6565024     -0.344701470             -0.30671445
#> 2             -0.1908479      0.544545671             -1.02680015
#> 3             -0.9924835     -0.271521165             -1.39445393
#> 4             -0.1606679      0.490120114             -0.28562291
#> 5             -0.6046320      0.009812293              0.33776943
#> 6              0.7704207      0.860022686              0.04093458
#>   CD8_T_cells_Bindea_et_al Tgd_Bindea_et_al Cytotoxic_cells_Bindea_et_al
#> 1              -0.36441112     -0.078457058                   -0.4794958
#> 2              -0.59451922     -1.088129202                   -0.7655065
#> 3               0.71050124     -0.012578397                    0.8289005
#> 4               0.33177846      0.284971755                    1.1731478
#> 5              -0.06178685      0.005159577                    0.3881670
#> 6               0.69906719      0.352896065                    0.7059167
#>   NK_cells_Bindea_et_al NK_CD56dim_cells_Bindea_et_al DC_Bindea_et_al
#> 1            0.21407345                    -1.9019014     -0.13081146
#> 2            0.03573202                    -0.7558300     -0.53493909
#> 3           -0.56383754                     0.5972162      0.71057151
#> 4            0.39073846                     1.7645762      0.04874792
#> 5           -0.71755388                     0.9209950      0.76281063
#> 6            1.01137164                     0.4832673      0.73565256
#>   iDC_Bindea_et_al aDC_Bindea_et_al Eosinophils_Bindea_et_al
#> 1        0.2909590      -0.51530209                0.4698789
#> 2        0.1324264       0.02778804                0.7279321
#> 3       -0.1148392       0.37069536               -0.7337183
#> 4        0.1583992       0.03404668                0.2739067
#> 5        1.5408800       0.27607673               -0.3817719
#> 6        0.4512418       0.13326571                1.4830447
#>   Macrophages_Bindea_et_al Mast_cells_Bindea_et_al Neutrophils_Bindea_et_al
#> 1                0.9695245               1.1536636                0.6671336
#> 2                1.1284447               1.8012561                0.2949121
#> 3               -1.1636915              -0.1278880               -0.7179112
#> 4                1.8010536              -1.1756683                1.1227808
#> 5               -1.5772201               0.5499247                1.6688433
#> 6                0.4929379               2.5792129                1.7344041
#>   SW480_cancer_cells_Bindea_et_al Normal_mucosa_Bindea_et_al
#> 1                       0.1261678                  1.1701735
#> 2                      -3.3347978                  0.7117883
#> 3                      -2.2615391                 -1.1216002
#> 4                       2.2849981                  0.7728001
#> 5                       1.3868872                 -1.7287382
#> 6                      -1.5831931                  1.5539341
#>   Lymph_vessels_Bindea_et_al Antigen_Processing_and_Presentation_Li_et_al
#> 1                 -0.1441374                                   -1.2684699
#> 2                 -0.7521155                                    0.2125215
#> 3                  0.6272827                                    0.1557994
#> 4                 -0.9722146                                    1.4842918
#> 5                 -0.5589525                                    0.7927213
#> 6                  0.5633130                                    1.0898181
#>   Antimicrobials_Li_et_al BCR_Signaling_Pathway_Li_et_al Chemokines_Li_et_al
#> 1                0.488577                      0.2192788           0.6570261
#> 2                1.514465                      0.8346698           1.2570658
#> 3               -1.663750                     -0.2634934          -0.6281273
#> 4                3.592499                     -0.9043938           2.2183807
#> 5               -0.494438                     -0.6094461           0.5369521
#> 6                3.194787                      2.3128994           1.7488803
#>   Chemokine_Receptors_Li_et_al Cytokines_Li_et_al Cytokine_Receptors_Li_et_al
#> 1                   0.42537976          1.2778112                   0.7550736
#> 2                   0.02769955          2.2086744                   0.9457428
#> 3                   0.35822512         -2.1194316                  -1.0276466
#> 4                   0.46662254          3.5266765                   0.9599109
#> 5                   0.42220725         -0.9815684                  -0.4191840
#> 6                   1.49041621          3.1058306                   3.8105694
#>   Interferons_Li_et_al Interferon_Receptor_Li_et_al Interleukins_Li_et_al
#> 1           0.01666841                   0.51574843             1.7509033
#> 2           0.51209405                  -0.07296792             0.2224877
#> 3          -0.83273083                  -0.09374811            -0.4122271
#> 4           0.82303534                   0.02494052             1.8071557
#> 5          -0.95137117                   0.55954041             0.4233798
#> 6           0.86985809                   0.08781163             1.7040556
#>   Interleukins_Receptor_Li_et_al Natural_Killer_Cell_Cytotoxicity_Li_et_al
#> 1                    -0.31065185                               -1.32085302
#> 2                    -0.05321163                                0.01769597
#> 3                    -0.19601653                                0.62638591
#> 4                     1.11367004                                2.09049297
#> 5                     0.68883328                                0.31949593
#> 6                     1.27744590                                1.17547265
#>   TCR_signaling_Pathway_Li_et_al TGFb_Family_Member_Li_et_al
#> 1                     -0.5412752                  -0.3223404
#> 2                     -0.3499472                  -0.7585347
#> 3                      0.8521841                   0.9232442
#> 4                      0.7312735                  -1.6809685
#> 5                      0.2842419                   1.8567536
#> 6                      1.8214163                   0.7483625
#>   TGFb_Family_Member_Receptor_Li_et_al TNF_Family_Members_Li_et_al
#> 1                            0.5843470                   0.3344842
#> 2                            0.2689723                  -0.3060312
#> 3                           -1.2634639                   0.6036527
#> 4                            1.0334668                   0.3915612
#> 5                           -1.9375114                  -0.5982484
#> 6                           -0.2195380                   0.4130962
#>   TNF_Family_Members_Receptors_Li_et_al T_cell_accumulation_Peng_et_al
#> 1                            0.27889496                      0.1334780
#> 2                            0.38396959                      0.1980172
#> 3                           -0.05256704                     -0.3934805
#> 4                           -0.32065482                      0.2498266
#> 5                            0.32515158                     -0.5559951
#> 6                            0.89426776                      0.5747637
#>   T_cell_exhaustion_Peng_et_al T_cell_regulatory_Peng_et_al
#> 1                  -0.15763204                  -0.06001087
#> 2                  -0.07398469                   0.44821448
#> 3                  -0.53673714                  -0.51017961
#> 4                   0.92111368                   1.47432778
#> 5                  -0.70391506                  -0.03857981
#> 6                   1.59590262                   0.77658340
#>   ICB_resistance_Peng_et_al MDSC_Peng_et_al TAM_Peng_et_al CAF_Peng_et_al
#> 1                0.05304805       -1.904373     -0.2861066       2.836353
#> 2                0.80330672        2.653911      1.3898696       0.796241
#> 3               -0.95435908       -1.606399     -2.1096741      -1.511045
#> 4                1.58377902        1.569549      1.4431035       2.277270
#> 5                0.47318549       -2.146400      0.2184291      -3.175397
#> 6                0.54063222       -1.577902      1.3460055       1.724684
#>   T_cell_probe_Jeschke_et_al MeTIL_Jeschke_et_al Cardiolipin_Metabolism
#> 1                 0.67985054         -0.82376104            -0.87639685
#> 2                -0.21793951          0.09036723            -0.20779114
#> 3                -0.35236666          1.19432566             0.15371012
#> 4                 0.52517540         -0.58652926            -0.32299277
#> 5                 0.09668882          0.40838996             0.05636391
#> 6                 0.39007038          0.29317114            -0.80107363
#>   Cardiolipin_Biosynthesis Cholesterol_Biosynthesis Citric_Acid_Cycle
#> 1               0.29774848             -0.079335673       0.589557146
#> 2              -0.09311578             -0.267912955      -0.004152856
#> 3               0.26823938              0.457963824      -0.653136363
#> 4              -0.18661084              0.002827104      -0.510773493
#> 5               0.11919024             -0.337586693      -1.330693043
#> 6              -0.07153350             -0.893562788       0.253742101
#>   Cyclooxygenase_Arachidonic_Acid_Metabolism Prostaglandin_Biosynthesis
#> 1                                  0.7637823                 -0.2569709
#> 2                                 -0.3639180                 -1.3805836
#> 3                                 -0.3921126                  0.1670359
#> 4                                  0.4067062                 -0.1488430
#> 5                                 -0.2714927                  0.7189272
#> 6                                  0.8272591                 -0.4704466
#>   Purine_Biosynthesis Pyrimidine_Biosynthesis Dopamine_Biosynthesis
#> 1           0.5052538             -0.12378559             0.7093124
#> 2           0.8796072              0.81762633            -0.9498221
#> 3          -0.2363282              0.02827347             0.5612763
#> 4           0.6336158              0.30854245             0.1417517
#> 5           0.6857567             -0.13948086             0.2687690
#> 6          -0.4109335             -0.56573783            -0.7681105
#>   Epinephrine_Biosynthesis Norepinephrine_Biosynthesis Fatty_Acid_Degradation
#> 1                0.4648963                   0.5243915             0.05117126
#> 2               -0.8389118                  -0.7076005            -2.40741459
#> 3                0.5717094                   0.4842038            -0.16996610
#> 4               -0.1607804                  -0.1638714            -2.22022438
#> 5                0.3290345                   0.4023478            -0.70136113
#> 6               -0.5958421                  -0.7419037            -0.61097931
#>   Fatty_Acid_Elongation Fatty_Acid_Biosynthesis Folate_One_Carbon_Metabolism
#> 1            0.20236754               0.3152738                    0.7927652
#> 2           -0.42397000              -0.7140774                   -0.3006251
#> 3            0.61684247               1.0408020                   -0.5054835
#> 4           -0.63108793              -1.5631382                    1.1586327
#> 5            0.01660887               0.2310818                   -0.6319715
#> 6            0.09472071              -0.5433039                   -0.8659583
#>   Folate_biosynthesis Gluconeogenesis Glycolysis Glycogen_Biosynthesis
#> 1           0.9772973      -0.1892089 -0.1775139            0.18261773
#> 2          -0.4132189      -1.4686688 -1.4748699           -1.00037902
#> 3           0.9827341       0.1228063  0.1028278            0.74865891
#> 4           1.6176266      -2.1489407 -2.1642463           -0.41326037
#> 5          -0.8104045      -1.3964497 -1.3815553            0.05592086
#> 6          -1.4808126      -0.8002555 -0.7649310           -0.17973633
#>   Glycogen_Degradation Heme_Biosynthesis Hexosamine_Biosynthesis
#> 1           -0.5369703        -0.3986826              0.18253764
#> 2           -0.3026259         0.1028228             -0.01419031
#> 3            0.5255660        -0.6352185             -0.10550973
#> 4           -0.4632253         0.2001857              0.13797009
#> 5           -0.3598158         0.3390117             -0.26201558
#> 6           -0.5817274        -1.0795410             -0.15102736
#>   Homocysteine_Biosynthesis Remethylation Transsulfuration
#> 1                -0.4126426     0.4546720        0.8312224
#> 2                -0.7447188    -2.1889581        0.8272014
#> 3                 0.5554031     1.7497581       -0.8092396
#> 4                -0.1044218    -1.0171845        0.5249471
#> 5                -0.1738129    -0.1852285       -0.7750495
#> 6                 0.5847990    -1.3536200       -0.7078339
#>   Ketone_Biosynthesis_and_Metabolism Kynurenine_Metabolism Methionine_Cycle
#> 1                         -0.5446141          -0.617425364        0.1992633
#> 2                         -2.6950830           0.008480219       -2.3583935
#> 3                          1.6236773           0.238632086        1.6764446
#> 4                         -1.2223798           0.304104716       -1.1073145
#> 5                          0.5662125           0.966748879       -0.3464087
#> 6                         -1.0679453           0.326813353       -1.2808038
#>   Pentose_Phosphate Nicotinamide_Adenine_Metabolism ADP_Ribosylation
#> 1        -0.9845481                       0.3971882     -0.305364598
#> 2        -0.5381216                      -0.2076272     -0.005250418
#> 3         0.8491269                      -0.5027833     -0.353697536
#> 4        -0.8705712                      -0.6650589      0.289479309
#> 5        -0.6825064                       0.1560991      0.022468916
#> 6        -0.7876161                       0.2534945      0.042026031
#>   Nicotinamide_Adenine_Dinucleotide_Biosynthesis
#> 1                                     -0.1288045
#> 2                                     -0.1725586
#> 3                                      0.5482463
#> 4                                      0.5704354
#> 5                                      1.1438242
#> 6                                      0.3243339
#>   Nicotinate_and_Nicotinamide_Metabolism Sirtuin_Nicotinamide_Metabolism
#> 1                              0.5363386                    -0.007188352
#> 2                              0.2176954                     0.435902855
#> 3                             -1.0404520                     0.155025090
#> 4                             -0.2070006                     0.361484432
#> 5                             -0.3270581                    -0.031404959
#> 6                              0.8051460                    -0.237533630
#>   Polyamine_Biosynthesis Prostanoid_Biosynthesis Pyruvate_Metabolism
#> 1             0.17528119              0.87300045           0.6818372
#> 2             0.22471920              0.20858970          -0.7649732
#> 3            -0.10623896             -0.70787250          -0.9764174
#> 4             0.05102969              0.02765451          -0.8377932
#> 5            -0.29074709             -0.97484144          -0.7685666
#> 6            -0.48784151              0.85494999           0.5620090
#>   Retinoic_Acid_Metabolism Retinoid_Metabolism Retinol_Metabolism
#> 1              -0.31517079         -0.97135907          1.1888791
#> 2              -0.60790760         -0.66359266         -6.8641843
#> 3               1.03946589          0.20419205         -0.8528126
#> 4              -1.33788672          0.26000603         -1.4483519
#> 5              -1.55027003         -1.58712468          1.8634256
#> 6              -0.08125257         -0.02733729         -4.2502602
#>   Steroid_Hormone_Metabolism Aldosterone_Biosynthesis Cortisol_Biosynthesis
#> 1                  1.2355678                1.4524042             0.2193730
#> 2                 -6.7779065               -0.1389591            -1.3061684
#> 3                 -0.6617941               -0.3447326            -0.9772724
#> 4                 -0.8700307                0.1233044             0.3633240
#> 5                  2.9992239               -1.0924741            -1.1244378
#> 6                 -4.1892890               -0.7147450            -1.5208844
#>   Estradiol_Biosynthesis Testosterone_Biosynthesis  Urea_Cycle   Vitamin_K
#> 1              0.3356946                 0.4128637 -0.72250884  0.17632752
#> 2             -1.4117722                -1.3626964 -0.32616832  0.07140197
#> 3              0.5545593                 1.2199398 -0.50794879 -0.08160066
#> 4             -0.9161043                -0.9359648 -0.06628019 -0.04152470
#> 5             -0.9442747                -0.9643749  0.14254166 -0.25559160
#> 6             -1.9317825                -1.8826810  0.09774957 -0.11454807
#>   Pentose_and_Glucuronate_Interconversions Fructose_and_Mannose_Metabolism
#> 1                                1.4391253                      -1.0176347
#> 2                               -5.8620651                      -1.2819516
#> 3                               -1.5418449                       0.1813223
#> 4                                0.1024901                      -0.8006759
#> 5                                2.8145888                      -0.2758227
#> 6                               -4.4973083                      -1.1530171
#>   Galactose_Metabolism Ascorbate_and_Aldrate_Metabolism
#> 1           -0.1527452                        1.4399847
#> 2           -1.6107905                       -5.7609872
#> 3           -0.8821752                       -1.5071409
#> 4           -0.3270414                        0.1135092
#> 5            0.8401661                        2.7008086
#> 6           -0.6823908                       -4.4485478
#>   Starch_and_Suctose_Metabolism Amino_Sugar_and_Nucleotide_Sugar_Metabolism
#> 1                    -1.4467038                                 0.002387852
#> 2                    -0.5851452                                -0.478834644
#> 3                    -0.6108423                                -0.763960281
#> 4                    -0.9903755                                 0.374050103
#> 5                    -0.4946986                                -0.943161354
#> 6                    -0.6394448                                 0.931719544
#>   Glyoxylate_and_Dicarboxylate_Metabolism Propanoate_Metabolism
#> 1                              0.05638195            -0.5506041
#> 2                             -0.35938217             1.5039057
#> 3                              0.50586219             0.7577437
#> 4                              0.37847771            -0.2371526
#> 5                             -0.39854608            -0.7407263
#> 6                             -0.45266219             0.4158993
#>   Butanoate_Metabolism Inositol_Phosphate_Metabolism Oxidative_Phosphorylation
#> 1            0.2738003                     0.1501629                -1.0070033
#> 2           -3.1458566                    -0.6722825                 0.4580757
#> 3            2.2111725                     1.5569383                 1.2754299
#> 4           -2.1718632                    -0.8141192                 0.4959665
#> 5            0.2067034                     0.8297760                 0.1384076
#> 6           -0.9710174                    -0.2074525                -1.5029417
#>   Nitrogen_Metabolism Sulfur_Metabolism Steroid_Biosynthesis
#> 1          1.06373355        0.05894934            1.0379938
#> 2         -0.04721592       -0.05088345           -1.3567910
#> 3         -0.20901247        0.58532279           -0.7536134
#> 4          1.01320009       -0.34529310            0.3505268
#> 5         -0.22256162        0.31389413            1.7543847
#> 6         -1.49258563       -0.47507000           -0.2746261
#>   Primary_Bile_Acid_Biosynthesis Steroid_Hormone_Biosynthesis
#> 1                    -0.60383380                    1.5078449
#> 2                    -0.27710831                   -6.8000833
#> 3                    -0.48171476                   -0.4965040
#> 4                    -0.44004735                   -0.9261952
#> 5                     0.01824793                    2.8971364
#> 6                    -1.03629287                   -4.1332957
#>   Glycerolipid_Metabolism Glycerophospholipid_Metabolism Ether_Lipid_Metabolism
#> 1              -0.6236669                      0.8166329              0.4455289
#> 2              -0.8828233                     -2.4737065             -2.9691184
#> 3              -1.5201793                      1.9378934              1.8056923
#> 4               0.9294664                     -0.8794945             -0.4903165
#> 5               0.3261957                      0.8036792              0.8675731
#> 6              -0.5665587                     -1.1309004             -0.8098890
#>   Shingolipid_Metabolism Arachidonic_Acid_Metabolism Linoleic_Acid_Metabolism
#> 1             -0.1552188                  -0.1743386                0.3276277
#> 2             -2.2377518                  -5.0695517               -3.0625327
#> 3              2.1047106                   1.2183305                0.4686198
#> 4             -0.1338881                  -2.1951792               -1.1462117
#> 5              1.0019147                   1.1897163                0.4986423
#> 6             -0.3999281                  -2.0410958               -1.0281125
#>   alpha_Linoleic_Acid_Metabolism Biosynthesis_of_Unsaturated_Fatty_Acids
#> 1                      0.5207145                             -0.17601454
#> 2                     -3.1944368                              0.29922914
#> 3                      1.4348617                             -0.43760950
#> 4                     -0.9539400                             -0.03858443
#> 5                      0.8571238                             -0.20652295
#> 6                     -0.6345843                             -1.00691660
#>   Purine_Metabolism Pyrimidine_Metabolism
#> 1        -0.1689800            0.20531439
#> 2         0.6825340           -0.14658066
#> 3         0.9066876            1.14812974
#> 4         0.6308725            0.03545628
#> 5         0.4675419            0.31003474
#> 6        -1.8919682           -1.44541593
#>   Alanine_Aspartate_and_Glutamate_Metabolism
#> 1                                 -1.0412476
#> 2                                  0.4405471
#> 3                                  0.6841108
#> 4                                 -0.7208254
#> 5                                  0.1240910
#> 6                                 -0.9893993
#>   Glycine_Serine_and_Threonine_Metabolism Cysteine_and_Methionine_Metabolism
#> 1                               0.2079144                         0.05223889
#> 2                              -2.2178190                        -1.46442751
#> 3                               1.5635828                         1.80992805
#> 4                              -0.8702598                        -1.26367495
#> 5                              -0.5865305                        -0.57155785
#> 6                              -1.3592609                        -1.30793007
#>   Valine_Leucine_and_Isoleucine_Degradation
#> 1                                -0.7249696
#> 2                                -2.4790681
#> 3                                 1.7718992
#> 4                                -1.4749999
#> 5                                 0.8892679
#> 6                                -0.9614294
#>   Valine_Leucine_and_Isoleucine_Biosynthesis Lysine_Degradation
#> 1                                 0.63548015         -0.7692155
#> 2                                 0.06238261         -0.1704759
#> 3                                -0.79853260         -0.7171521
#> 4                                 0.51945340          0.5526805
#> 5                                -0.85346468          0.4932688
#> 6                                 0.09954195         -0.7659554
#>   Arginine_Biosynthesis Arginine_and_Proline_Metabolism Histidine_Metabolism
#> 1           -0.70986724                      -0.9591439          -0.32830352
#> 2           -0.26113362                      -0.8168890           1.76318310
#> 3           -0.48471166                       0.6744821           0.04546724
#> 4           -0.03941625                      -0.8595765          -0.40044863
#> 5            0.11001704                       0.9697592          -1.12901165
#> 6           -0.05217103                      -1.5101096           0.37503444
#>   Tyrosine_Metabolism Phenylalanine_Metabolism Tryptophan_Metabolism
#> 1           0.3212506              -1.20785720             3.8488450
#> 2          -2.1722026               0.08959426            -1.0967881
#> 3          -1.2258498              -0.97158727             3.0515478
#> 4          -2.2764718              -0.90655903            -1.4381016
#> 5          -0.9510747              -1.66887499            -0.8714761
#> 6          -0.5150489              -0.90795678             0.6237964
#>   Phenylalanine_Tyrosine_and_Tryptophan_Biosynthesis Beta_Alanine_Metabolism
#> 1                                       -1.326187295              0.08697081
#> 2                                        0.001283902             -1.98139486
#> 3                                       -0.901271155             -0.44491674
#> 4                                       -0.054087147              0.14138015
#> 5                                       -1.298074304              1.09880664
#> 6                                       -1.467254299             -0.84172754
#>   Taurine_and_Hypotaurine_Metabolism Selenocompound_Metabolism
#> 1                         -0.2907417                0.40501244
#> 2                          1.9871809               -0.01314174
#> 3                         -0.6714827               -0.24207624
#> 4                          0.2416510                0.01974161
#> 5                         -0.5523903               -0.49220381
#> 6                         -0.5895731                0.27864082
#>   D_Glutamine_and_D_Glutamate_Metabolism Glutathione_Metabolism
#> 1                            -0.23405155              0.7983287
#> 2                             0.38048713             -2.1868292
#> 3                             0.46949944             -0.1159699
#> 4                            -0.15921354             -2.2746492
#> 5                            -0.07535048              2.2466983
#> 6                            -0.61842922             -1.5814406
#>   N_Glycan_Biosynthesis Mucin_Type_O_Glycan_Biosynthesis
#> 1             0.4353414                      -0.67189634
#> 2             0.9025926                      -1.15063902
#> 3            -0.5812715                      -0.53448991
#> 4             0.5791943                       0.80524401
#> 5            -0.9533280                       0.27591662
#> 6             0.4115728                      -0.09387302
#>   Other_Types_of_O_Glycan_Biosynthesis Glycosaminoglycan_Biosynthesis
#> 1                            0.2850761                     -0.1698166
#> 2                           -1.0237209                      0.8620824
#> 3                            0.5847049                     -0.4159240
#> 4                           -0.7590142                      1.8604599
#> 5                            0.2384218                     -0.6934627
#> 6                            0.4426812                     -0.3896467
#>   Glycosaminoglycan_Degradation Glycosphosphatidylinositol
#> 1                    0.24972922                  0.0188355
#> 2                    0.32581390                 -0.1781743
#> 3                   -0.79593818                  0.7709734
#> 4                    0.03422953                 -0.7666202
#> 5                    0.10205619                  0.3910387
#> 6                    0.88833401                 -0.4880614
#>   Glycosphingolipid_Biosynthesis Other_Glycan_Degradation Thiamine_Metabolism
#> 1                      0.0190998               0.05148039         -0.12312600
#> 2                     -2.4815919              -0.51188657          0.05316807
#> 3                      0.3422093              -0.31922361         -0.06236520
#> 4                     -1.3362474              -0.34217018         -0.09652619
#> 5                      1.5354055              -0.18568811          0.25693894
#> 6                      0.1717414              -0.16514652         -0.07324444
#>   Riboflavin_Metabolism Vitamin_B6_Metabolism Pantothenate_and_CoA_Biosynthesis
#> 1          -0.252822679             0.5702018                       -0.24681844
#> 2          -0.030906427            -0.1082251                       -0.03297744
#> 3          -0.191280570            -0.5346490                       -1.35672719
#> 4           0.075681741            -0.5210649                       -0.49539781
#> 5          -0.009496592            -0.6520421                       -0.02635510
#> 6          -0.438454830             0.7523788                        0.95020236
#>   Biotin_Metabolism Lipoic_Acid_Metabolism Porphyrin_and_Chlorophyll_Metabolism
#> 1        0.10049438          -0.0267403507                            1.4475326
#> 2       -0.08106278          -0.0436920565                           -5.7706135
#> 3        0.24562750           0.2291923937                           -1.4855098
#> 4       -0.23012123          -0.5297972842                            0.1755997
#> 5        0.08460314           0.2690609124                            2.7881254
#> 6       -0.15307965          -0.0006316579                           -4.4018210
#>   Ubiquinone_and_other_Terpenoid_Quinone_Biosynthesis
#> 1                                          0.15059045
#> 2                                         -0.27473979
#> 3                                         -0.95683070
#> 4                                         -0.69596667
#> 5                                         -0.76891324
#> 6                                         -0.08057671
#>   Terpenoid_Backbone_Biosynthesis Caffiene_Metabolism
#> 1                      -0.5741884           2.6606297
#> 2                      -2.6426055          -0.8883417
#> 3                       1.5155121           1.4226773
#> 4                      -1.1317680          -0.7463807
#> 5                       0.6290543          -0.7978596
#> 6                      -1.1223782           0.4785393
#>   Neomycin_Kanamysin_and_Gentamicin_Biosynthesis
#> 1                                     -0.5639853
#> 2                                     -0.8583111
#> 3                                     -0.3700284
#> 4                                     -0.1691826
#> 5                                     -0.2751318
#> 6                                     -0.4848240
#>   Metabolism_of_Xenobiotics_by_Cytochrome_P450
#> 1                                    1.5290186
#> 2                                   -6.8485295
#> 3                                   -0.9518103
#> 4                                   -1.9274310
#> 5                                    2.0059461
#> 6                                   -4.5364679
#>   Drug_Metabolism_by_Cytochrome_P450 Drug_Metabolism_by_other_enzymes
#> 1                          0.8877538                       1.39845628
#> 2                         -6.9435982                      -5.73229702
#> 3                         -1.7341265                      -1.48112359
#> 4                         -1.5019997                       0.03005906
#> 5                          2.3554130                       2.75362093
#> 6                         -4.6460036                      -4.40083458
#>   Nature_metabolism_Hypoxia Winter_hypoxia_signature Hu_hypoxia_signature
#> 1                0.18990643               1.15548820            0.9476420
#> 2               -0.26176795              -0.06364701            0.6452329
#> 3               -0.50383455               0.41000138           -1.0850317
#> 4                1.16503638               0.41535135            1.5826813
#> 5               -0.05814433              -0.22108619           -0.9651191
#> 6               -1.22012497              -1.45210722           -0.7038365
#>   Molecular_Cancer_m6A  MT_exosome Positive_regulation_of_exosomal_secretion
#> 1            0.1567624  0.33220137                                0.14374798
#> 2            2.1141484  0.08341014                               -0.85992231
#> 3            0.8293245 -0.31863225                                0.50299954
#> 4            0.6758115  0.09627660                               -0.01432638
#> 5           -1.3560695 -0.28127124                                0.14234492
#> 6           -1.4728362 -0.01607636                               -0.61502427
#>   Ferroptosis EV_Cell_2020 HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION
#> 1  0.08154693    0.3627038             0.2258825                   0.10236953
#> 2  0.22882352    0.2350307             0.2074258                   0.09849066
#> 3 -0.91975939   -0.5883062             0.2235993                   0.08186991
#> 4  0.74487197    0.5811955             0.2231714                   0.13638849
#> 5 -0.55379867   -0.9646009             0.2223973                   0.10566475
#> 6 -0.70858896    0.5513844             0.2133946                   0.12666113
#>   HALLMARK_ANDROGEN_RESPONSE HALLMARK_ANGIOGENESIS HALLMARK_APICAL_JUNCTION
#> 1                  0.2997372             0.2695479                0.2048530
#> 2                  0.2515314             0.2057329                0.1850909
#> 3                  0.2711597             0.2363239                0.1724121
#> 4                  0.2827975             0.2637562                0.2188514
#> 5                  0.2708298             0.2118968                0.1883964
#> 6                  0.2780936             0.2603515                0.2184456
#>   HALLMARK_APICAL_SURFACE HALLMARK_APOPTOSIS HALLMARK_BILE_ACID_METABOLISM
#> 1              0.11994121          0.2376789                    0.10249791
#> 2              0.02979827          0.1977715                    0.07492853
#> 3              0.09729760          0.1975451                    0.08522369
#> 4              0.13886970          0.2402368                    0.06366086
#> 5              0.10184063          0.2183151                    0.10310188
#> 6              0.11342885          0.2302876                    0.09782265
#>   HALLMARK_CHOLESTEROL_HOMEOSTASIS HALLMARK_COAGULATION HALLMARK_COMPLEMENT
#> 1                        0.2118394           0.15621946           0.1694040
#> 2                        0.1703181           0.12595033           0.1481347
#> 3                        0.2088638           0.09585739           0.1298536
#> 4                        0.2200324           0.17454378           0.2007805
#> 5                        0.1882506           0.13524916           0.1806912
#> 6                        0.2021830           0.15455439           0.1834767
#>   HALLMARK_DNA_REPAIR HALLMARK_E2F_TARGETS
#> 1           0.1740762            0.2515750
#> 2           0.2098816            0.2978351
#> 3           0.1907927            0.2624527
#> 4           0.2002712            0.2753313
#> 5           0.1704336            0.1905408
#> 6           0.1479990            0.1914812
#>   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_ESTROGEN_RESPONSE_EARLY
#> 1                                  0.2972794                       0.16027423
#> 2                                  0.2484435                       0.06731134
#> 3                                  0.1970552                       0.12387049
#> 4                                  0.2960018                       0.14322910
#> 5                                  0.1707949                       0.14584780
#> 6                                  0.2733503                       0.15085643
#>   HALLMARK_ESTROGEN_RESPONSE_LATE HALLMARK_FATTY_ACID_METABOLISM
#> 1                      0.14219629                      0.1895202
#> 2                      0.08401326                      0.1922885
#> 3                      0.12434032                      0.2004007
#> 4                      0.15876028                      0.1974802
#> 5                      0.13429057                      0.1902265
#> 6                      0.13108140                      0.1741947
#>   HALLMARK_G2M_CHECKPOINT HALLMARK_GLYCOLYSIS HALLMARK_HEDGEHOG_SIGNALING
#> 1               0.2609024           0.2006028                  0.14587911
#> 2               0.3056389           0.1856343                  0.15677725
#> 3               0.2713343           0.1803533                  0.09483918
#> 4               0.2902021           0.2110865                  0.13851493
#> 5               0.2062817           0.1783986                  0.13498113
#> 6               0.2090930           0.1639742                  0.12695454
#>   HALLMARK_HEME_METABOLISM HALLMARK_HYPOXIA HALLMARK_IL2_STAT5_SIGNALING
#> 1                0.1689150        0.2235898                    0.1642260
#> 2                0.1505027        0.1702739                    0.1224207
#> 3                0.1493296        0.1607908                    0.1294439
#> 4                0.1667817        0.2224210                    0.1506948
#> 5                0.1794674        0.1783082                    0.1395556
#> 6                0.1610178        0.1961416                    0.1700355
#>   HALLMARK_IL6_JAK_STAT3_SIGNALING HALLMARK_INFLAMMATORY_RESPONSE
#> 1                        0.1556276                     0.07125957
#> 2                        0.1115441                     0.03470739
#> 3                        0.1121064                     0.02150742
#> 4                        0.1820841                     0.11507588
#> 5                        0.1604430                     0.06750874
#> 6                        0.1841053                     0.10018913
#>   HALLMARK_INTERFERON_ALPHA_RESPONSE HALLMARK_INTERFERON_GAMMA_RESPONSE
#> 1                          0.1534134                          0.1642819
#> 2                          0.2256531                          0.1915338
#> 3                          0.1661001                          0.1537477
#> 4                          0.2984049                          0.2445330
#> 5                          0.2415602                          0.2126528
#> 6                          0.1963487                          0.2080084
#>   HALLMARK_KRAS_SIGNALING_DN HALLMARK_KRAS_SIGNALING_UP
#> 1                -0.10479749                 0.13058984
#> 2                -0.11458048                 0.06636304
#> 3                -0.13376052                 0.07167808
#> 4                -0.08531253                 0.11564844
#> 5                -0.07355827                 0.09996630
#> 6                -0.08629210                 0.13484113
#>   HALLMARK_MITOTIC_SPINDLE HALLMARK_MTORC1_SIGNALING HALLMARK_MYC_TARGETS_V1
#> 1                0.3213105                 0.2814548               0.3372598
#> 2                0.3312961                 0.2663423               0.3784446
#> 3                0.3256150                 0.2618298               0.3495329
#> 4                0.3374321                 0.2919227               0.3587917
#> 5                0.2984236                 0.2409340               0.3114124
#> 6                0.3039122                 0.2329167               0.2959519
#>   HALLMARK_MYC_TARGETS_V2 HALLMARK_MYOGENESIS HALLMARK_NOTCH_SIGNALING
#> 1               0.1972419          0.14538731                0.1995705
#> 2               0.2817367          0.12155189                0.1379046
#> 3               0.2305404          0.11017413                0.1955100
#> 4               0.2137894          0.13222590                0.1742360
#> 5               0.1789928          0.09908186                0.2015088
#> 6               0.1590548          0.17536150                0.1793302
#>   HALLMARK_OXIDATIVE_PHOSPHORYLATION HALLMARK_P53_PATHWAY
#> 1                          0.2470695            0.1911455
#> 2                          0.2753947            0.1295974
#> 3                          0.2681116            0.1663144
#> 4                          0.2780856            0.2009326
#> 5                          0.2558388            0.1998776
#> 6                          0.2200683            0.1673710
#>   HALLMARK_PANCREAS_BETA_CELLS HALLMARK_PEROXISOME
#> 1                  -0.03578828           0.1761825
#> 2                   0.01433006           0.1680596
#> 3                  -0.07356179           0.1776193
#> 4                  -0.03063279           0.1676174
#> 5                  -0.08465927           0.1811698
#> 6                  -0.03551709           0.1603892
#>   HALLMARK_PI3K_AKT_MTOR_SIGNALING HALLMARK_PROTEIN_SECRETION
#> 1                        0.2341569                  0.3212248
#> 2                        0.2142833                  0.3137465
#> 3                        0.2447455                  0.3189630
#> 4                        0.2457933                  0.3366214
#> 5                        0.2343474                  0.3144233
#> 6                        0.2234422                  0.3015594
#>   HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY HALLMARK_SPERMATOGENESIS
#> 1                                0.2331076              0.032791906
#> 2                                0.2219856              0.094314972
#> 3                                0.2203466              0.027962738
#> 4                                0.2383843              0.041341386
#> 5                                0.2308216              0.007887468
#> 6                                0.2044808              0.008314541
#>   HALLMARK_TGF_BETA_SIGNALING HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 1                   0.3262102                        0.1738830
#> 2                   0.2888964                        0.0964380
#> 3                   0.3006051                        0.1129213
#> 4                   0.3211058                        0.1787389
#> 5                   0.3159277                        0.1563443
#> 6                   0.3133321                        0.2151073
#>   HALLMARK_UNFOLDED_PROTEIN_RESPONSE HALLMARK_UV_RESPONSE_DN
#> 1                          0.2959271               0.2898912
#> 2                          0.2865260               0.2615268
#> 3                          0.2753978               0.2479981
#> 4                          0.2930936               0.2839940
#> 5                          0.2545039               0.2705874
#> 6                          0.2566789               0.3137656
#>   HALLMARK_UV_RESPONSE_UP HALLMARK_WNT_BETA_CATENIN_SIGNALING
#> 1               0.1682020                           0.1563698
#> 2               0.1472323                           0.1668078
#> 3               0.1368461                           0.1421387
#> 4               0.1600101                           0.1516424
#> 5               0.1277375                           0.1539792
#> 6               0.1432875                           0.1471895
#>   HALLMARK_XENOBIOTIC_METABOLISM GO_CELLULAR_RESPONSE_TO_FATTY_ACID
#> 1                     0.16500800                         0.11666224
#> 2                     0.09189626                         0.03902495
#> 3                     0.13325514                         0.12390176
#> 4                     0.12534810                         0.07461225
#> 5                     0.13777436                         0.10934775
#> 6                     0.13314910                         0.09506540
#>   GO_CHOLESTEROL_EFFLUX GO_FATTY_ACID_BETA_OXIDATION
#> 1            0.12720785                    0.1711063
#> 2            0.13613218                    0.1717612
#> 3            0.08695757                    0.2022022
#> 4            0.12562603                    0.1600479
#> 5            0.10074141                    0.2213414
#> 6            0.15282926                    0.1728636
#>   GO_FATTY_ACID_BIOSYNTHETIC_PROCESS GO_FATTY_ACID_CATABOLIC_PROCESS
#> 1                         0.08347456                       0.1635600
#> 2                         0.05930992                       0.1612430
#> 3                         0.07839497                       0.1832077
#> 4                         0.08548022                       0.1526224
#> 5                         0.06032075                       0.2133246
#> 6                         0.07159021                       0.1652061
#>   GO_FATTY_ACID_DERIVATIVE_TRANSPORT GO_FATTY_ACID_ELONGATION
#> 1                       -0.004329622               0.06873250
#> 2                       -0.080023030               0.02275221
#> 3                       -0.016786792               0.05108986
#> 4                       -0.053621039               0.09440015
#> 5                        0.050768621               0.01421334
#> 6                       -0.015064786               0.03763442
#>   GO_FATTY_ACID_HOMEOSTASIS GO_FATTY_ACID_METABOLIC_PROCESS
#> 1                0.15569328                      0.10896534
#> 2                0.15844036                      0.07299275
#> 3                0.12385967                      0.10397307
#> 4                0.11732684                      0.08477918
#> 5                0.09501693                      0.10797537
#> 6                0.09104664                      0.09062819
#>   GO_FATTY_ACID_TRANSPORT GO_LONG_CHAIN_FATTY_ACID_TRANSPORT
#> 1              0.10637487                         0.09149423
#> 2              0.06427564                         0.07121955
#> 3              0.09170987                         0.07592504
#> 4              0.06571926                         0.06449145
#> 5              0.13011471                         0.09012226
#> 6              0.09117602                         0.09430780
#>   GO_OXIDATIVE_DEMETHYLATION GO_OXIDATIVE_PHOSPHORYLATION
#> 1                0.022581880                    0.2125692
#> 2               -0.052461427                    0.2579551
#> 3               -0.006460393                    0.2560846
#> 4               -0.054149149                    0.2595542
#> 5                0.041689681                    0.2283524
#> 6               -0.066618449                    0.1923685
#>   GO_RESPONSE_TO_CORTICOSTEROID GO_RESPONSE_TO_FATTY_ACID
#> 1                    0.11117154                0.11534092
#> 2                    0.06266375                0.04305225
#> 3                    0.07559518                0.11431561
#> 4                    0.08826689                0.07056640
#> 5                    0.08638579                0.09401233
#> 6                    0.10671688                0.10140394
#>   GO_RESPONSE_TO_OXIDATIVE_STRESS GO_RESPONSE_TO_STEROID_HORMONE
#> 1                       0.2001834                      0.1330262
#> 2                       0.1767738                      0.0979506
#> 3                       0.1814513                      0.1100478
#> 4                       0.1961388                      0.1126234
#> 5                       0.1925789                      0.1060669
#> 6                       0.1896079                      0.1285447
#>   GO_REVERSE_CHOLESTEROL_TRANSPORT GO_STEROID_BIOSYNTHETIC_PROCESS
#> 1                       0.20177971                      0.06462203
#> 2                       0.07953961                      0.01533833
#> 3                       0.09323439                      0.05862961
#> 4                       0.11143678                      0.04650181
#> 5                       0.08956274                      0.04543077
#> 6                       0.11811471                      0.01661113
#>   GO_STEROID_CATABOLIC_PROCESS GO_STEROID_METABOLIC_PROCESS
#> 1                   0.04437358                   0.09388795
#> 2                  -0.01934774                   0.02775147
#> 3                  -0.01458488                   0.07315233
#> 4                   0.02971295                   0.05460669
#> 5                  -0.00344240                   0.08379218
#> 6                  -0.02224991                   0.03847498
#>   GO_CD40_RECEPTOR_COMPLEX GO_FATTY_ACID_BINDING GO_FATTY_ACID_LIGASE_ACTIVITY
#> 1                0.2541449          0.0639066387                   0.064266710
#> 2                0.3175371         -0.0694487930                   0.006081296
#> 3                0.2460457         -0.0003518796                   0.150378743
#> 4                0.2749684          0.0716091379                   0.037022046
#> 5                0.2840244          0.0851138240                   0.102389802
#> 6                0.2664919          0.0083054378                   0.074041566
#>   GO_FATTY_ACID_TRANSPORTER_ACTIVITY GO_LONG_CHAIN_FATTY_ACID_BINDING
#> 1                       -0.016820393                       0.06499884
#> 2                        0.008870766                      -0.08338347
#> 3                        0.041919187                       0.05475594
#> 4                       -0.004064023                       0.09017462
#> 5                        0.005223742                       0.11662270
#> 6                        0.019584007                       0.02567145
#>   GO_STEROID_BINDING GO_STEROID_DEHYDROGENASE_ACTIVITY
#> 1         0.08605715                        0.01540340
#> 2         0.02590768                       -0.05325848
#> 3         0.04986943                       -0.03558623
#> 4         0.03981594                       -0.03430492
#> 5         0.06253328                        0.04175350
#> 6         0.06899577                       -0.06413916
#>   GO_STEROID_HYDROXYLASE_ACTIVITY KEGG_FATTY_ACID_METABOLISM
#> 1                      -0.1158503                  0.1627621
#> 2                      -0.2661372                  0.1592814
#> 3                      -0.2001081                  0.1751884
#> 4                      -0.2224021                  0.1618627
#> 5                      -0.2161116                  0.1749305
#> 6                      -0.2312790                  0.1702057
#>   KEGG_GLYCOLYSIS_GLUCONEOGENESIS KEGG_OXIDATIVE_PHOSPHORYLATION
#> 1                       0.1884304                      0.1702298
#> 2                       0.1709801                      0.2026034
#> 3                       0.1573668                      0.2138705
#> 4                       0.1910164                      0.2189309
#> 5                       0.1882730                      0.1892493
#> 6                       0.1772983                      0.1380269
#>   KEGG_STEROID_BIOSYNTHESIS KEGG_STEROID_HORMONE_BIOSYNTHESIS
#> 1                 0.1757052                        0.01819952
#> 2                 0.1681551                       -0.24753733
#> 3                 0.1754245                       -0.06096648
#> 4                 0.1892304                       -0.10421769
#> 5                 0.1413013                        0.03005918
#> 6                 0.1104809                       -0.18845221
#>   REACTOME_CHOLESTEROL_BIOSYNTHESIS REACTOME_DOWNSTREAM_TCR_SIGNALING
#> 1                         0.2023676                         0.2030581
#> 2                         0.1981701                         0.1534079
#> 3                         0.2488644                         0.2255162
#> 4                         0.2151139                         0.2129382
#> 5                         0.1773750                         0.2235095
#> 6                         0.1498952                         0.2189046
#>   REACTOME_GLYCOLYSIS REACTOME_STEROID_HORMONES REACTOME_TCR_SIGNALING
#> 1           0.2682875               -0.03526465              0.2039928
#> 2           0.2921093               -0.13290789              0.1739174
#> 3           0.2347783               -0.15559565              0.2294578
#> 4           0.2685199               -0.06951794              0.2239982
#> 5           0.2189748               -0.06403291              0.2234833
#> 6           0.2336288               -0.09688094              0.2331010

#患者治疗反应数据和生存结局
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

### 4.2 使用sig\_group来批量可视化与治疗反应相关的signature

``` r
#signature group包含这些大类
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

# 与肿瘤相关的signature
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

# 与免疫治疗标志物相关的signature
names(sig_group)[3]
#> [1] "io_biomarkers"
sig_group[[3]]
#>  [1] "TMEscore_CIR"                    "TMEscoreA_CIR"                  
#>  [3] "TMEscoreB_CIR"                   "T_cell_inflamed_GEP_Ayers_et_al"
#>  [5] "CD_8_T_effector"                 "IPS_IPS"                        
#>  [7] "Immune_Checkpoint"               "Exhausted_CD8_Danaher_et_al"    
#>  [9] "Pan_F_TBRs"                      "Mismatch_Repair"                
#> [11] "APM"

# 与免疫抑制相关的signature
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

### 4.3 构建分组信息批量查看免疫治疗反应与选择的signature group的关系

``` r
pdata_group<-imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR","BOR_binary")]
res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_matrix = imvigor210_sig,
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
                   palette_heatmap = 1,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>> Is NA exist:  0"
#>    Mode    TRUE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#> [1] ">>>  Processing signature: tumor_signature"
```

<img src="man/figuresunnamed-chunk-25-1.png" width="100%" /><img src="man/figuresunnamed-chunk-25-2.png" width="100%" />

    #> [1] ">>>  Processing signature: io_biomarkers"

<img src="man/figuresunnamed-chunk-25-3.png" width="100%" /><img src="man/figuresunnamed-chunk-25-4.png" width="100%" />

    #> [1] ">>>  Processing signature: immu_suppression"

<img src="man/figuresunnamed-chunk-25-5.png" width="100%" /><img src="man/figuresunnamed-chunk-25-6.png" width="100%" />

### 4.4 计算与治疗反应相关的signature genes

``` r

#构建表达矩阵的输入文件
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
eset_input<-rownames_to_column(as.data.frame(t(imvigor210_eset)),var = "ID")

res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_matrix = eset_input,
                   id2 = "ID",
                   target= NULL,
                   group = "BOR_binary",
                   is_target_continuous = FALSE,
                   padj_cutoff = 1,
                   index = 1,
                   category = "signature",
                   signature_group =signature_collection[c(1:3)],
                   ProjectID = "IMvigor210",
                   palette_box = "paired1",
                   palette_corplot = "pheatmap",
                   palette_heatmap = 1,
                   feature_limit = 26,
                   character_limit = 30,
                   show_heatmap_col_name = FALSE,
                   show_col = FALSE,
                   show_plot = TRUE)
#> [1] ">>> Is NA exist:  0"
#>    Mode    TRUE 
#> logical     257 
#>    Mode   FALSE 
#> logical     257 
#>    Mode   FALSE 
#> logical     257 
#>    Mode   FALSE 
#> logical     257 
#> [1] ">>>  Processing signature: CD_8_T_effector"
```

<img src="man/figuresunnamed-chunk-26-1.png" width="100%" /><img src="man/figuresunnamed-chunk-26-2.png" width="100%" />

    #> [1] ">>>  Processing signature: DDR"

<img src="man/figuresunnamed-chunk-26-3.png" width="100%" /><img src="man/figuresunnamed-chunk-26-4.png" width="100%" />

    #> [1] ">>>  Processing signature: APM"

<img src="man/figuresunnamed-chunk-26-5.png" width="100%" /><img src="man/figuresunnamed-chunk-26-6.png" width="100%" />

## 5 分析与目标signature相关的signatures

### 5.1 构建pdata\_group作为患者的表型数据

``` r
#数据已经储存在`imvigor210_pdata`中
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

### 5.2 分析与Pan-F-TBRs相关的signatures

`Pan-F-TBRs`是2018年Nature上一篇研究免疫驱逐的评分: S. Mariathasan et al., TGFbeta
attenuates tumour response to PD-L1 blockade by contributing to
exclusion of T cells. Nature 554, 544-548 (2018).
我们可以使用计算出来的signature来批量分析和可视化这个免疫驱逐评分与其他`signature`和`cell
fraction` 的关系，并探索可能影响该生物过程的其他分子标志物。

``` r
res<-iobr_cor_plot(pdata_group = pdata_group,
                   id1 = "ID",
                   feature_matrix = imvigor210_sig,
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
#> [1] ">>> Is NA exist:  0"
#>    Mode    TRUE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#> [1] ">>>  Processing signature: tumor_signature"
```

<img src="man/figuresunnamed-chunk-28-1.png" width="100%" /><img src="man/figuresunnamed-chunk-28-2.png" width="100%" /><img src="man/figuresunnamed-chunk-28-3.png" width="100%" />

    #> [1] ">>>  Processing signature: EMT"

<img src="man/figuresunnamed-chunk-28-4.png" width="100%" /><img src="man/figuresunnamed-chunk-28-5.png" width="100%" /><img src="man/figuresunnamed-chunk-28-6.png" width="100%" />

### 5.3 分析与Pan-F-TBRs相关的免疫细胞浸润情况

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
                   feature_matrix = imvigor210_sig,
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
#> [1] ">>> Is NA exist:  0"
#>    Mode    TRUE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
#>    Mode   FALSE 
#> logical     455 
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

## 6 Seraching for mutations related to signature

### 6.1 Format mutation matrix using MAF data. 这里使用了TCGA-STAD的数据为例.

MAF data was obtained from [UCSC Xena
website](https://api.gdc.cancer.gov/data/c06465a3-50e7-46f7-b2dd-7bd654ca206b)

``` r
help("make_mut_matrix")

maf_file<-"TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf"
#使用MAF的数据将每个患者的每个基因转化成二分类变量（有或者无突变）
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
#> -Finished in 14.8s elapsed (13.1s cpu) 
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
#如下所示，如果上述方法选择”multi", 结果将产生四个data frame, 分别计算了每个基因的所有突变，或者只计算SNP, indel 和 frameshift,
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

### 6.2 使用sig\_group来批量可视化与治疗反应相关的signature

``` r
#'将要分析的突变数据转化为等级变量矩阵----------
mut<-mut_list$snp
max(mut)
#> [1] 22
#如果突变矩阵中单个基因的最大突变数据大于4，将会把突变数量按照以下规则进行转化
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

## 7 其他批量分析方法

IOBR提供了多个批量统计分析、数据过滤和转化的方法：包括生存分析，批量寻找最佳cutoff,批量生存分析，批量统计检验等功能

### 7.1 subgroup survival analyses

``` r
help("subgroup_survival")
##source data and filter NA
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

Contact: E-mail any questions to <dongqiangzeng0808@gmail.com>
