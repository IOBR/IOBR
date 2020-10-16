
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IOBR: Immune Oncology Bioinformatics Research

IOBR is a R package to perform Tumor microenvironment evaluation,
signature estimation.

## IOBR简介

  - 1.IOBR集合了8种已经发表的肿瘤微环境解析方法：`CIBERSORT`, `TIMER`, `xCell`,
    `MCPcounter`, `ESITMATE`, `EPIC`, `IPS`, `quanTIseq`;
  - 2.IOBR集合了共256个已经发表的Signature gene sets：包括肿瘤微环境相关的，肿瘤代谢相关、m6A,
    外泌体相关的, 微卫星不稳定, 三级淋巴结评分等，可通过函数 `signatures_sci`
    获取到signature的出处；通过 `signature_collection` 可以获取到每个signature
    gene;
  - 3.IOBR集合了三种方法用于上述signature评分的计算，包括`PCA`,`z-score`,`ssGSEA`;
  - 4.IOBR集合了多种方法用于变量转化和批量生存分析和统计学分析的方法；
  - 5.IOBR集合了批量可视化分组特征的方法；

![](man/figures/IOBR-Package.png)

## 安装依赖包

IOBR依赖包较多，包括：tibble, survival, survminer, limma, limSolve, GSVA, e1071,
preprocessCore, ggplot2, ggpubr;

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

The package is not yet on CRAN. You can install from Github:

``` r
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("DongqiangZeng0808/IOBR",ref="master")
```

加载包

``` r
library(IOBR) 
library(EPIC)
library(estimate) 
library(MCPcounter)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
```

## 解析肿瘤微环境

### 可选择的肿瘤微环境解析方法

``` r
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVM 
#>    "cibersort_abs"              "ips"         "estimate"          "svm_ref" 
#>               lsei              TIMER          quanTIseq 
#>         "lsei_ref"            "timer"        "quantiseq"
#每种方法所对应选择的参数
```
###所用方法的许可证（还在收集中）
| method | license | citation |
|--------|---------|----------|
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) | free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)) | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., ..., Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6 |
| [TIMER](http://cistrome.org/TIMER/) | free ([GPL 2.0](http://cistrome.org/TIMER/download.html)) | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174.  https://doi.org/10.1186/s13059-016-1028-7 |
| [CIBERSORT](https://cibersort.stanford.edu/) | free for non-commerical use only | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457.  https://doi.org/10.1038/nmeth.3337 |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License)) | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5 |
| [xCell](http://xcell.ucsf.edu/) | free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION)) | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1 |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476 |


输入文件为100例TCGA-CAOD的RNAseq数据,已采用log2(TPM)进行标准化。

``` r
# 查看数据
eset_crc[1:5,1:5]
#>         TCGA-3L-AA1B TCGA-4N-A93T TCGA-4T-AA8H TCGA-5M-AAT4 TCGA-5M-AAT5
#> MT-CO2      14.68459     15.68270     15.95073     15.10946     15.05872
#> MT-CO3      15.04403     16.07435     15.74394     14.68427     14.51087
#> MT-ND4      14.59023     14.76067     15.12249     14.45566     14.49834
#> MT-CO1      13.89565     14.44195     15.19789     14.00815     14.56622
#> MT-ATP6     13.62783     14.42318     14.74342     13.98928     14.15316
```

查看函数具体的参数

``` r
help(deconvo_tme)
#> starting httpd help server ... done
```

方法1：使用CIBERSORT解析肿瘤微环境

``` r
#默认进行1000次permutation, 这里为了节约时间设成了100
cibersort<-deconvo_tme(eset = eset_crc,method = "cibersort",arrays = FALSE,perm = 100 )
#> 
#> >>> Running CIBERSORT
head(cibersort)
#>             ID B_cells_naive_CIBERSORT B_cells_memory_CIBERSORT
#> 1 TCGA-3L-AA1B              0.11107873              0.000000000
#> 2 TCGA-4N-A93T              0.10590915              0.013085986
#> 3 TCGA-4T-AA8H              0.00000000              0.002118305
#> 4 TCGA-5M-AAT4              0.02329158              0.000000000
#> 5 TCGA-5M-AAT5              0.18155805              0.000000000
#> 6 TCGA-5M-AAT6              0.02229106              0.000000000
#>   Plasma_cells_CIBERSORT T_cells_CD8_CIBERSORT T_cells_CD4_naive_CIBERSORT
#> 1             0.01122323            0.06819909                           0
#> 2             0.06802228            0.16677876                           0
#> 3             0.19397170            0.11767817                           0
#> 4             0.00000000            0.27889238                           0
#> 5             0.20533622            0.03335308                           0
#> 6             0.05495993            0.14670277                           0
#>   T_cells_CD4_memory_resting_CIBERSORT T_cells_CD4_memory_activated_CIBERSORT
#> 1                           0.26092818                             0.02311397
#> 2                           0.17495073                             0.00000000
#> 3                           0.35990743                             0.01241794
#> 4                           0.00000000                             0.00000000
#> 5                           0.19829705                             0.04160562
#> 6                           0.08458969                             0.10138945
#>   T_cells_follicular_helper_CIBERSORT T_cells_regulatory_(Tregs)_CIBERSORT
#> 1                          0.02345378                           0.08726958
#> 2                          0.15759492                           0.07592051
#> 3                          0.00000000                           0.03383531
#> 4                          0.05796164                           0.06535554
#> 5                          0.07783184                           0.02678012
#> 6                          0.04863560                           0.03532881
#>   T_cells_gamma_delta_CIBERSORT NK_cells_resting_CIBERSORT
#> 1                             0                 0.02959256
#> 2                             0                 0.00000000
#> 3                             0                 0.01570515
#> 4                             0                 0.02560453
#> 5                             0                 0.00816428
#> 6                             0                 0.01709894
#>   NK_cells_activated_CIBERSORT Monocytes_CIBERSORT Macrophages_M0_CIBERSORT
#> 1                  0.000000000         0.038833418               0.11488268
#> 2                  0.046516938         0.011469881               0.05281655
#> 3                  0.000000000         0.008099187               0.07829807
#> 4                  0.000000000         0.009758635               0.22571531
#> 5                  0.008876187         0.001817252               0.07675599
#> 6                  0.000000000         0.000000000               0.14926728
#>   Macrophages_M1_CIBERSORT Macrophages_M2_CIBERSORT
#> 1              0.018511758               0.09907659
#> 2              0.029589739               0.05294717
#> 3              0.001952015               0.07286956
#> 4              0.055724721               0.11461051
#> 5              0.028243224               0.07386405
#> 6              0.124314639               0.15155714
#>   Dendritic_cells_resting_CIBERSORT Dendritic_cells_activated_CIBERSORT
#> 1                        0.01987360                         0.000000000
#> 2                        0.04242040                         0.001976983
#> 3                        0.05996742                         0.027616054
#> 4                        0.00000000                         0.001191356
#> 5                        0.00000000                         0.016756071
#> 6                        0.01896034                         0.000000000
#>   Mast_cells_resting_CIBERSORT Mast_cells_activated_CIBERSORT
#> 1                   0.05763218                     0.03633067
#> 2                   0.00000000                     0.00000000
#> 3                   0.00000000                     0.01556369
#> 4                   0.00000000                     0.11895257
#> 5                   0.02076096                     0.00000000
#> 6                   0.04490436                     0.00000000
#>   Eosinophils_CIBERSORT Neutrophils_CIBERSORT P-value_CIBERSORT
#> 1            0.00000000                     0              0.01
#> 2            0.00000000                     0              0.01
#> 3            0.00000000                     0              0.00
#> 4            0.02294122                     0              0.37
#> 5            0.00000000                     0              0.01
#> 6            0.00000000                     0              0.01
#>   Correlation_CIBERSORT RMSE_CIBERSORT
#> 1            0.16057040      1.0185423
#> 2            0.18833970      1.0271235
#> 3            0.25557817      0.9923352
#> 4            0.01919799      1.0924515
#> 5            0.19349231      1.0343501
#> 6            0.17832546      1.0114621
```

方法2：使用EPIC解析肿瘤微环境

``` r
epic<-deconvo_tme(eset = eset_crc,method = "epic",arrays = FALSE)
#> 
#> >>> Running EPIC
#> Warning in EPIC::EPIC(bulk = eset, reference = ref, mRNA_cell = NULL, scaleExprs = TRUE): The optimization didn't fully converge for some samples:
#> TCGA-A6-2677; TCGA-A6-5657; TCGA-A6-5659; TCGA-A6-5664; TCGA-AA-3496; TCGA-AA-3502; TCGA-AA-3511; TCGA-AA-3516; TCGA-AA-3520; TCGA-AA-3527; TCGA-AA-3530; TCGA-AA-3553; TCGA-AA-3554; TCGA-AA-3561
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning in EPIC::EPIC(bulk = eset, reference = ref, mRNA_cell = NULL, scaleExprs
#> = TRUE): mRNA_cell value unknown for some cell types: CAFs, Endothelial - using
#> the default value of 0.4 for these but this might bias the true cell proportions
#> from all cell types.
head(epic)
#>             ID Bcells_EPIC  CAFs_EPIC CD4_Tcells_EPIC CD8_Tcells_EPIC
#> 1 TCGA-3L-AA1B  0.04621723 0.01648197       0.1642720      0.07574850
#> 2 TCGA-4N-A93T  0.03279542 0.01203903       0.1405250      0.07509584
#> 3 TCGA-4T-AA8H  0.01851783 0.01204811       0.1848611      0.05236068
#> 4 TCGA-5M-AAT4  0.01011495 0.01535950       0.1086244      0.08721530
#> 5 TCGA-5M-AAT5  0.01825654 0.01442913       0.1273391      0.07064480
#> 6 TCGA-5M-AAT6  0.04043042 0.02223098       0.1868112      0.10462340
#>   Endothelial_EPIC Macrophages_EPIC NKcells_EPIC otherCells_EPIC
#> 1       0.12215267      0.007972055 9.308906e-10       0.5671556
#> 2       0.08184873      0.006157295 6.019345e-11       0.6515387
#> 3       0.07556058      0.005665098 1.156554e-11       0.6509866
#> 4       0.10465578      0.007621791 1.393203e-09       0.6664083
#> 5       0.08779865      0.007065009 6.025925e-10       0.6744668
#> 6       0.15303721      0.012025094 4.479847e-10       0.4808417
```

方法3：使用EPIC解析肿瘤微环境

``` r
mcp<-deconvo_tme(eset = eset_crc,method = "mcpcounter")
#> 
#> >>> Running MCP-counter
head(mcp)
#>             ID T_cells_MCPcounter CD8_T_cells_MCPcounter
#> 1 TCGA-3L-AA1B           2.276328              2.4716953
#> 2 TCGA-4N-A93T           1.332391              1.7721989
#> 3 TCGA-4T-AA8H           1.202586              0.8480327
#> 4 TCGA-5M-AAT4           1.044777              2.5670796
#> 5 TCGA-5M-AAT5           1.244461              1.8985864
#> 6 TCGA-5M-AAT6           2.459058              3.1384888
#>   Cytotoxic_lymphocytes_MCPcounter NK_cells_MCPcounter B_lineage_MCPcounter
#> 1                        0.7747745          0.12963782             3.375311
#> 2                        0.4311466          0.15266260             2.421203
#> 3                        0.3958956          0.11984715             1.569042
#> 4                        0.5136194          0.11343929             1.327973
#> 5                        0.5019803          0.08874996             1.570996
#> 6                        2.0216901          0.73760870             3.092513
#>   Monocytic_lineage_MCPcounter Myeloid_dendritic_cells_MCPcounter
#> 1                    2.4584809                          1.8543679
#> 2                    1.3996255                          1.0973198
#> 3                    0.9545705                          0.8126323
#> 4                    1.5593173                          0.4168602
#> 5                    1.4396720                          0.5333424
#> 6                    3.8595154                          1.5635735
#>   Neutrophils_MCPcounter Endothelial_cells_MCPcounter Fibroblasts_MCPcounter
#> 1               1.806034                     2.714606               6.753796
#> 2               1.241315                     1.715088               3.958849
#> 3               1.652448                     1.477787               3.999647
#> 4               1.554660                     2.190213               5.219359
#> 5               1.446799                     1.976331               4.835505
#> 6               1.744977                     3.169847               7.293699
```

方法4：使用xCell解析肿瘤微环境

    #> [1] "Num. of genes: 10776"
    #> Estimating ssGSEA scores for 489 gene sets.
    #>   |                                                                              |                                                                      |   0%Using parallel with 4 cores
    #>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
    #>             ID    aDC_xCell Adipocytes_xCell Astrocytes_xCell B-cells_xCell
    #> 1 TCGA-3L-AA1B 1.036618e-01     2.398965e-19     9.536865e-19  1.851477e-02
    #> 2 TCGA-4N-A93T 3.629409e-03     0.000000e+00     4.775877e-18  3.652864e-02
    #> 3 TCGA-4T-AA8H 3.554369e-19     0.000000e+00     0.000000e+00  4.697424e-20
    #> 4 TCGA-5M-AAT4 3.576842e-03     1.007067e-18     5.793605e-18  1.374562e-18
    #> 5 TCGA-5M-AAT5 6.229837e-03     9.269526e-19     5.371013e-18  1.787947e-03
    #> 6 TCGA-5M-AAT6 3.292397e-01     1.224713e-19     3.910092e-02  2.879303e-02
    #>   Basophils_xCell CD4+_memory_T-cells_xCell CD4+_naive_T-cells_xCell
    #> 1     0.000000000               0.004305412             9.597158e-03
    #> 2     0.013071529               0.000000000             0.000000e+00
    #> 3     0.000000000               0.006397612             1.653949e-19
    #> 4     0.000759893               0.004195499             9.793798e-19
    #> 5     0.018385533               0.008618021             0.000000e+00
    #> 6     0.012503349               0.026490181             4.669957e-03
    #>   CD4+_T-cells_xCell CD4+_Tcm_xCell CD4+_Tem_xCell CD8+_naive_T-cells_xCell
    #> 1       3.486596e-21    0.018278142   7.011902e-03              0.007550990
    #> 2       6.527096e-20    0.029826405   1.209102e-02              0.005457778
    #> 3       7.178576e-19    0.002092308   8.731603e-03              0.011192672
    #> 4       0.000000e+00    0.000000000   0.000000e+00              0.010977229
    #> 5       0.000000e+00    0.001809237   1.863407e-18              0.003509766
    #> 6       6.167899e-05    0.000000000   4.837433e-19              0.013147208
    #>   CD8+_T-cells_xCell CD8+_Tcm_xCell CD8+_Tem_xCell    cDC_xCell
    #> 1        0.005348034   7.443031e-19   3.398223e-18 0.0318318478
    #> 2        0.010435656   0.000000e+00   6.252190e-18 0.0006313437
    #> 3        0.008675741   0.000000e+00   2.754946e-19 0.0202504165
    #> 4        0.000000000   0.000000e+00   0.000000e+00 0.0004828999
    #> 5        0.001292990   3.583860e-19   5.461325e-18 0.0000000000
    #> 6        0.060124119   4.932981e-02   1.908784e-02 0.0150120118
    #>   Chondrocytes_xCell Class-switched_memory_B-cells_xCell   CLP_xCell
    #> 1       1.599979e-02                         0.009481405 0.008588866
    #> 2       1.370017e-18                         0.050404207 0.003198056
    #> 3       2.318615e-18                         0.010294194 0.016385956
    #> 4       0.000000e+00                         0.007169692 0.024174180
    #> 5       0.000000e+00                         0.001408652 0.025493733
    #> 6       4.278670e-03                         0.006146947 0.030524926
    #>      CMP_xCell     DC_xCell Endothelial_cells_xCell Eosinophils_xCell
    #> 1 0.000000e+00 3.830072e-03            0.0112829805      0.000000e+00
    #> 2 0.000000e+00 8.968523e-20            0.0008155795      1.515909e-02
    #> 3 1.489998e-18 0.000000e+00            0.0002708751      1.026898e-02
    #> 4 6.440033e-19 8.379382e-19            0.0040688697      0.000000e+00
    #> 5 0.000000e+00 8.452054e-20            0.0002722808      4.590789e-19
    #> 6 1.725539e-19 3.316656e-02            0.0516094130      0.000000e+00
    #>   Epithelial_cells_xCell Erythrocytes_xCell Fibroblasts_xCell    GMP_xCell
    #> 1            0.001170113                  0      1.272092e-02 3.815681e-02
    #> 2            0.010785487                  0      1.109346e-17 0.000000e+00
    #> 3            0.011754538                  0      2.472967e-18 3.477081e-18
    #> 4            0.002827428                  0      0.000000e+00 1.322572e-17
    #> 5            0.001923605                  0      5.148224e-18 0.000000e+00
    #> 6            0.003232257                  0      3.611154e-02 1.177597e-19
    #>   Hepatocytes_xCell  HSC_xCell    iDC_xCell Keratinocytes_xCell
    #> 1      0.000000e+00 0.11145539 0.000000e+00         0.000000000
    #> 2      8.795481e-19 0.02269690 3.036764e-03         0.001595504
    #> 3      0.000000e+00 0.02974023 6.498629e-03         0.012228449
    #> 4      3.433135e-21 0.04569535 1.575143e-03         0.005875808
    #> 5      0.000000e+00 0.01805671 1.711236e-19         0.003640640
    #> 6      1.679354e-18 0.04105797 1.898246e-19         0.009184715
    #>   ly_Endothelial_cells_xCell Macrophages_xCell Macrophages_M1_xCell
    #> 1               3.100422e-03      0.000000e+00         0.0060232585
    #> 2               4.550723e-19      0.000000e+00         0.0000000000
    #> 3               3.355326e-19      1.631160e-19         0.0000000000
    #> 4               9.016826e-04      1.427516e-18         0.0005150186
    #> 5               2.187801e-04      0.000000e+00         0.0000000000
    #> 6               2.693647e-02      4.427136e-02         0.0386547107
    #>   Macrophages_M2_xCell Mast_cells_xCell Megakaryocytes_xCell Melanocytes_xCell
    #> 1         1.103505e-05      0.009736625         1.494856e-04      0.000000e+00
    #> 2         7.103414e-20      0.010274437         8.966161e-20      3.384458e-03
    #> 3         0.000000e+00      0.017707357         0.000000e+00      0.000000e+00
    #> 4         0.000000e+00      0.007134796         9.260045e-04      4.201263e-20
    #> 5         0.000000e+00      0.009380248         0.000000e+00      1.641398e-19
    #> 6         1.277626e-02      0.018162033         3.920425e-03      0.000000e+00
    #>   Memory_B-cells_xCell  MEP_xCell Mesangial_cells_xCell Monocytes_xCell
    #> 1          0.004844646 0.04508331          2.916627e-02    2.231907e-18
    #> 2          0.016435014 0.03213538          7.131246e-19    2.257187e-18
    #> 3          0.002770904 0.03468093          0.000000e+00    0.000000e+00
    #> 4          0.000000000 0.07144444          5.647505e-03    0.000000e+00
    #> 5          0.005826091 0.06630378          6.485127e-03    3.006424e-18
    #> 6          0.010459491 0.05493735          1.753210e-03    1.655159e-02
    #>      MPP_xCell    MSC_xCell mv_Endothelial_cells_xCell Myocytes_xCell
    #> 1 0.000000e+00 0.000000e+00                0.006609335   2.065827e-03
    #> 2 0.000000e+00 2.904667e-02                0.003963793   2.285061e-18
    #> 3 2.820314e-18 2.641848e-17                0.000000000   4.853237e-03
    #> 4 1.412277e-17 1.360384e-17                0.002658978   2.458567e-03
    #> 5 2.703356e-18 1.290166e-17                0.001814830   1.402936e-03
    #> 6 0.000000e+00 1.275626e-01                0.034780693   0.000000e+00
    #>   naive_B-cells_xCell Neurons_xCell Neutrophils_xCell NK_cells_xCell
    #> 1        5.205385e-03  0.0017114228      0.0000000000   0.000000e+00
    #> 2        9.687372e-04  0.0019734895      0.0000000000   4.924421e-19
    #> 3        1.152748e-03  0.0016189748      0.0000000000   0.000000e+00
    #> 4        5.351319e-20  0.0014438709      0.0014802524   1.151411e-18
    #> 5        1.561950e-03  0.0020049632      0.0002781054   4.701879e-18
    #> 6        8.413473e-04  0.0003151171      0.0000000000   1.481603e-03
    #>     NKT_xCell Osteoblast_xCell    pDC_xCell Pericytes_xCell Plasma_cells_xCell
    #> 1 0.040452192     4.532036e-03 4.172042e-19    2.824396e-18        0.016228957
    #> 2 0.065312167     3.704818e-18 2.852102e-05    1.360564e-17        0.023053819
    #> 3 0.059163934     0.000000e+00 0.000000e+00    5.012782e-18        0.011493410
    #> 4 0.057289635     2.534327e-03 5.053096e-19    0.000000e+00        0.005325248
    #> 5 0.063735883     5.030752e-18 1.687703e-03    1.527995e-18        0.009231568
    #> 6 0.009853135     9.776863e-20 1.043905e-02    5.994726e-19        0.014289239
    #>   Platelets_xCell Preadipocytes_xCell pro_B-cells_xCell Sebocytes_xCell
    #> 1    1.116439e-18        0.000000e+00                 0    0.000000e+00
    #> 2    0.000000e+00        2.665500e-02                 0    1.912112e-04
    #> 3    4.761077e-03        1.536697e-02                 0    8.549712e-04
    #> 4    3.993451e-20        1.019290e-17                 0    1.037912e-04
    #> 5    5.083930e-18        9.024405e-03                 0    1.490059e-03
    #> 6    0.000000e+00        1.754079e-03                 0    1.086219e-18
    #>   Skeletal_muscle_xCell Smooth_muscle_xCell Tgd_cells_xCell Th1_cells_xCell
    #> 1          0.000000e+00           0.4103968    0.000000e+00     0.008742359
    #> 2          0.000000e+00           0.1859708    1.896191e-18     0.009083326
    #> 3          5.346708e-19           0.2139791    0.000000e+00     0.031631991
    #> 4          8.453081e-18           0.3209835    1.923142e-17     0.059103746
    #> 5          0.000000e+00           0.2255458    2.172994e-18     0.039190064
    #> 6          1.121273e-03           0.2939993    0.000000e+00     0.038067404
    #>   Th2_cells_xCell  Tregs_xCell ImmuneScore_xCell StromaScore_xCell
    #> 1    6.758450e-19 1.128412e-03       0.024952998      0.0120019509
    #> 2    3.865487e-19 9.021057e-20       0.048265216      0.0004077897
    #> 3    1.874381e-02 0.000000e+00       0.024434721      0.0001354375
    #> 4    6.208965e-02 3.207116e-19       0.005743365      0.0020344348
    #> 5    2.968342e-02 0.000000e+00       0.008492861      0.0001361404
    #> 6    1.484422e-01 9.274516e-03       0.135074649      0.0438604758
    #>   MicroenvironmentScore_xCell
    #> 1                 0.036954949
    #> 2                 0.048673006
    #> 3                 0.024570158
    #> 4                 0.007777800
    #> 5                 0.008629001
    #> 6                 0.178935125

方法5：使用ESTIMATE计算肿瘤纯度和免疫、间质评分

``` r
estimate<-deconvo_tme(eset = eset_crc,method = "estimate")
#> 
#> >>> Running ESTIMATE
#> [1] "Merged dataset includes 10186 genes (226 mismatched)."
#> [1] "1 gene set: StromalSignature  overlap= 139"
#> [1] "2 gene set: ImmuneSignature  overlap= 140"
head(estimate)
#>             ID StromalScore_estimate ImmuneScore_estimate
#> 1 TCGA-3L-AA1B             -644.2579             202.9177
#> 2 TCGA-4N-A93T            -1828.3122            -383.8928
#> 3 TCGA-4T-AA8H            -2030.2323            -574.6390
#> 4 TCGA-5M-AAT4            -1398.8294            -560.1262
#> 5 TCGA-5M-AAT5            -1581.8424            -505.8969
#> 6 TCGA-5M-AAT6              314.1607            1480.4165
#>   ESTIMATEScore_estimate TumorPurity_estimate
#> 1              -441.3403            0.8576040
#> 2             -2212.2050            0.9609832
#> 3             -2604.8712            0.9753218
#> 4             -1958.9555            0.9500391
#> 5             -2087.7393            0.9557695
#> 6              1794.5772            0.6460408
```

方法6：使用TIMER解析肿瘤微环境

    #> [1] "Outlier genes: ACTB ACTG1 EEF1A1 EEF2 FTL IGF2 IGHA1 IGHA2 IGHG4 MIR7641-1 MT-ATP6 MT-ATP8 MT-CO1 MT-CO2 MT-CO3 MT-CYB MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-RNR2 OLFM4 PIGR PPBP REG1A RN7SK RN7SL2 RNA5SP141 RNU4-2 RPL8 RPS12 RPS18 RPS2 RPS21 S100A6 SNORA73B TFF1 TFF3 TMSB10"
    #> Standardizing Data across genes
    #>             ID B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER Neutrophil_TIMER
    #> 1 TCGA-3L-AA1B   0.10954038        0.1450304        0.1893913        0.1144261
    #> 2 TCGA-4N-A93T   0.10625470        0.1304173        0.1717525        0.1061153
    #> 3 TCGA-4T-AA8H   0.09952132        0.1243976        0.1834921        0.1049141
    #> 4 TCGA-5M-AAT4   0.09863300        0.1171330        0.1859230        0.1119424
    #> 5 TCGA-5M-AAT5   0.10048711        0.1166910        0.1906776        0.1092820
    #> 6 TCGA-5M-AAT6   0.10199630        0.1319679        0.2247014        0.1271064
    #>   Macrophage_TIMER  DC_TIMER
    #> 1       0.04611901 0.5059655
    #> 2       0.03864497 0.4798004
    #> 3       0.03836999 0.4786281
    #> 4       0.04661717 0.4798606
    #> 5       0.03799281 0.4888228
    #> 6       0.05763071 0.5219651

方法7：使用quanTIseq解析肿瘤微环境

``` r
quantiseq<-deconvo_tme(eset = eset_crc, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE,method = "quantiseq")
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
#> 1 TCGA-3L-AA1B       0.013180386               0.02732029
#> 2 TCGA-4N-A93T       0.010114895               0.03933797
#> 3 TCGA-4T-AA8H       0.005394022               0.01767286
#> 4 TCGA-5M-AAT4       0.008079243               0.04239306
#> 5 TCGA-5M-AAT5       0.006985994               0.07677666
#> 6 TCGA-5M-AAT6       0.016130582               0.06606741
#>   Macrophages_M2_quantiseq Monocytes_quantiseq Neutrophils_quantiseq
#> 1              0.032106198                   0            0.03833859
#> 2              0.028721164                   0            0.03445506
#> 3              0.015349325                   0            0.04669382
#> 4              0.005832801                   0            0.03544526
#> 5              0.000000000                   0            0.02365587
#> 6              0.036378324                   0            0.02533244
#>   NK_cells_quantiseq T_cells_CD4_quantiseq T_cells_CD8_quantiseq
#> 1         0.03661911          0.0000000000          0.000000e+00
#> 2         0.01731620          0.0099323786          1.400995e-03
#> 3         0.02284371          0.0001493545          0.000000e+00
#> 4         0.03125493          0.0000000000          6.572346e-05
#> 5         0.03251279          0.0000000000          0.000000e+00
#> 6         0.01753874          0.0000000000          2.259896e-02
#>   Tregs_quantiseq Dendritic_cells_quantiseq Other_quantiseq
#> 1     0.033211797                         0       0.8192236
#> 2     0.010739959                         0       0.8479814
#> 3     0.012542265                         0       0.8793546
#> 4     0.006725298                         0       0.8702037
#> 5     0.008851506                         0       0.8512172
#> 6     0.040286144                         0       0.7756674
```

方法8：使用IPS评估肿瘤免疫表型

``` r
ips<-deconvo_tme(eset = eset_crc,method = "ips",plot= FALSE)
#>    Mode    TRUE 
#> logical     161 
#> [1] GENE   NAME   CLASS  WEIGHT
#> <0 rows> (or 0-length row.names)
#> [1] GENE   NAME   CLASS  WEIGHT
#> <0 rows> (or 0-length row.names)
head(ips)
#>             ID  MHC_IPS    EC_IPS     SC_IPS       CP_IPS   AZ_IPS IPS_IPS
#> 1 TCGA-3L-AA1B 2.825104 0.9506996 -1.0253824  0.001208434 2.751630       9
#> 2 TCGA-4N-A93T 3.102602 0.7209778 -0.5090616  0.217528831 3.532047      10
#> 3 TCGA-4T-AA8H 2.904127 0.7514136 -0.5030863  0.223891809 3.376346      10
#> 4 TCGA-5M-AAT4 2.762653 0.8138660 -0.7357093  0.119353981 2.960163      10
#> 5 TCGA-5M-AAT5 2.929669 0.8399299 -0.6818371  0.125965787 3.213727      10
#> 6 TCGA-5M-AAT6 3.739094 1.3668529 -1.5456332 -0.775121868 2.785192       9
```

合并所有的解析结果用于后续的分析

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
#> [1] 100 138
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

## Signature score 的评估

IOBR集合了共256个已经发表的Signature gene sets：包括肿瘤微环境相关的，肿瘤代谢相关、m6A, 外泌体相关的,
微卫星不稳定,
三级淋巴结评分等，可通过函数’signatures\_sci’获取到signature的出处；通过’signature\_collection’可以获取到每个signature
gene;

查看有哪些signature

``` r
#微环境相关的signature
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
#代谢相关的signature
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
#与基础研究相关的signature: 比如m6A, 外泌体
names(signature_star)
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
#signature 集合
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

### 评估肿瘤微环境相关的signature-(使用PCA方法）

``` r
sig_tme<-calculate_sig_score(pdata = NULL,eset = eset_crc,
                             signature = signature_tme,
                             method = "pca",
                             mini_gene_count = 2)
#> 
#> >>> Calculating signature score with PCA method
sig_tme[1:5,1:10]
#>   Index           ID CD_8_T_effector        DDR        APM Immune_Checkpoint
#> 1     1 TCGA-3L-AA1B       -2.189621  0.2507418 -1.5677265        -0.2968556
#> 2     2 TCGA-4N-A93T       -2.845395 -1.6537649 -0.2954285        -1.4342785
#> 3     3 TCGA-4T-AA8H       -3.147606  0.3759185 -0.5825841        -1.4957523
#> 4     4 TCGA-5M-AAT4       -2.228884  1.6015096 -1.0314272        -1.4290535
#> 5     5 TCGA-5M-AAT5       -2.153055  1.3332892 -0.7496979        -1.3052797
#>   CellCycle_Reg Pan_F_TBRs   Histones       EMT1
#> 1     0.2304255  0.4732723 -1.2571610  0.4080749
#> 2     0.5520772 -2.3541168 -0.8635123 -1.2441486
#> 3     0.0478045 -1.8223391 -1.0953553 -1.5203636
#> 4     0.2550369 -0.7726379 -0.8457354 -0.6146912
#> 5     0.4405093 -1.3725215 -0.8528409 -1.0677128
```

### 评估肿瘤微环境相关的signature-（使用ssGSEA方法）

``` r
sig_tme<-calculate_sig_score(pdata = NULL,eset = eset_crc,
                                 signature = signature_tme,
                                 method = "ssgsea",
                                 mini_gene_count = 5)
#> Estimating ssGSEA scores for 98 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_tme[1:5,1:10]
#>             ID Index CD_8_T_effector       DDR       APM Immune_Checkpoint
#> 1 TCGA-3L-AA1B     1       0.2278595 0.4118750 0.5543922         0.2486442
#> 2 TCGA-4N-A93T     2       0.1828740 0.4043511 0.5706539         0.1588728
#> 3 TCGA-4T-AA8H     3       0.1979758 0.4200638 0.5711350         0.1514483
#> 4 TCGA-5M-AAT4     4       0.2399799 0.4371537 0.5703392         0.1574774
#> 5 TCGA-5M-AAT5     5       0.2490083 0.4344176 0.5743728         0.1756249
#>   CellCycle_Reg Pan_F_TBRs      EMT1      EMT2
#> 1     0.4505993  0.4996082 0.4866532 0.3910946
#> 2     0.4216721  0.4328508 0.4431105 0.2346882
#> 3     0.4667955  0.4505855 0.4348335 0.2740186
#> 4     0.4494207  0.4817414 0.4636002 0.3204985
#> 5     0.4791662  0.4639342 0.4571042 0.2981027
```

### 评估代谢相关的signature

``` r
sig_meta<-calculate_sig_score(pdata = NULL,eset = eset_crc,
                                 signature = signature_metabolism,
                                 method = "pca",
                                 mini_gene_count = 2)
#> 
#> >>> Calculating signature score with PCA method
sig_meta[1:5,1:10]
#>   Index           ID Cardiolipin_Metabolism Cardiolipin_Biosynthesis
#> 1     1 TCGA-3L-AA1B             -0.3346619               0.29183707
#> 2     2 TCGA-4N-A93T             -0.3719745               0.28752661
#> 3     3 TCGA-4T-AA8H              0.8551322              -0.20626172
#> 4     4 TCGA-5M-AAT4             -0.1914077              -0.09619904
#> 5     5 TCGA-5M-AAT5              0.8174483              -0.39563770
#>   Cholesterol_Biosynthesis Citric_Acid_Cycle
#> 1              -0.52863584        0.12676471
#> 2               0.07114562       -0.45176171
#> 3               0.21171068        0.01730886
#> 4               0.60531372        0.06331044
#> 5               0.48259640       -0.49013422
#>   Cyclooxygenase_Arachidonic_Acid_Metabolism Prostaglandin_Biosynthesis
#> 1                               -0.195428422                 -0.6746913
#> 2                               -2.381801421                  0.7776066
#> 3                               -0.879067676                 -0.2492614
#> 4                               -1.758651867                 -0.7762122
#> 5                               -0.007635041                 -0.3520407
#>   Purine_Biosynthesis Pyrimidine_Biosynthesis
#> 1          -0.1467908               0.2195573
#> 2          -0.7895191               0.5663729
#> 3          -0.3122503              -0.3758337
#> 4           0.9013135               0.3203994
#> 5           0.4255790               0.1777062
```

### 计算所有收集的signature score（综合三种的方法： PCA, ssGSEA和z-score）

``` r
sig_res<-calculate_sig_score(pdata = NULL,eset = eset_crc,
                                 signature = signature_collection,
                                 method = "integration",
                                 mini_gene_count = 2)
#> Estimating ssGSEA scores for 213 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_res[1:5,1:10]
#>             ID Index CD_8_T_effector_PCA    DDR_PCA    APM_PCA
#> 1 TCGA-3L-AA1B     1           -2.189621  0.2507418 -1.5677265
#> 2 TCGA-4N-A93T     2           -2.845395 -1.6537649 -0.2954285
#> 3 TCGA-4T-AA8H     3           -3.147606  0.3759185 -0.5825841
#> 4 TCGA-5M-AAT4     4           -2.228884  1.6015096 -1.0314272
#> 5 TCGA-5M-AAT5     5           -2.153055  1.3332892 -0.7496979
#>   Immune_Checkpoint_PCA CellCycle_Reg_PCA Pan_F_TBRs_PCA Histones_PCA
#> 1            -0.2968556         0.2304255      0.4732723   -1.2571610
#> 2            -1.4342785         0.5520772     -2.3541168   -0.8635123
#> 3            -1.4957523         0.0478045     -1.8223391   -1.0953553
#> 4            -1.4290535         0.2550369     -0.7726379   -0.8457354
#> 5            -1.3052797         0.4405093     -1.3725215   -0.8528409
#>     EMT1_PCA
#> 1  0.4080749
#> 2 -1.2441486
#> 3 -1.5203636
#> 4 -0.6146912
#> 5 -1.0677128
```

### IOBR还集合了GO, KEGG, HALLMARK, REACTOME的signature gene sets

建议使用ssGSEA的方法进行评估,如果样本量比较大且signature比较多，运行时间会比较长，

``` r
sig_hallmark<-calculate_sig_score(pdata = NULL,
                             eset = eset_crc,
                             signature = hallmark ,
                             method = "ssgsea",
                             mini_gene_count = 2)
#> Estimating ssGSEA scores for 50 gene sets.
#>   |                                                                              |                                                                      |   0%Using parallel with 8 cores
#>   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
sig_hallmark[1:5,1:10]
#>             ID Index HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION
#> 1 TCGA-3L-AA1B     1             0.9800903                    0.7819648
#> 2 TCGA-4N-A93T     2             0.9545930                    0.6852996
#> 3 TCGA-4T-AA8H     3             0.9772683                    0.6985885
#> 4 TCGA-5M-AAT4     4             0.9670854                    0.6980879
#> 5 TCGA-5M-AAT5     5             0.9593519                    0.7205536
#>   HALLMARK_ANDROGEN_RESPONSE HALLMARK_ANGIOGENESIS HALLMARK_APICAL_JUNCTION
#> 1                  0.9677298             0.8744085                0.8457583
#> 2                  0.9112817             0.7191048                0.7835878
#> 3                  0.9517066             0.7172882                0.7919595
#> 4                  0.9379049             0.8464839                0.8166375
#> 5                  0.9496281             0.8208311                0.7920659
#>   HALLMARK_APICAL_SURFACE HALLMARK_APOPTOSIS HALLMARK_BILE_ACID_METABOLISM
#> 1               0.7438635          0.9683096                     0.7398229
#> 2               0.6964792          0.9143470                     0.6924335
#> 3               0.7027327          0.9319835                     0.7160735
#> 4               0.7101088          0.9470103                     0.6656917
#> 5               0.6859652          0.9672865                     0.7117280
```

### IOBR提供了多个批量统计分析、数据过滤和转化的方法：包括生存分析，批量寻找最佳cutoff,批量生存分析，批量统计检验等功能

subgroup survival analyses

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
