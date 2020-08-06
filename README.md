
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IOBR: Immune Oncology Bioinformatics Research

IOBR is a R package to perform Tumor microenvironment evaluation,
signature estimation.

## Main advantages:

  - Provides functionality to estimate signature score (PCA, ssGSEA or
    z-score)
  - Integrated prevalent methodologies to deconvolute microenvironment,
    such as CIBERSORT, EPIC, ESTIMATE, xCell
  - Functions to visulize TME data.
  - Identify TME relevant signatures.
  - Functions to perform batch analyses for prognostic and predictive
    biomarkers

The package is not yet on CRAN. You can install from Github:

``` r
# install.packages("devtools")
devtools::install_github("DongqiangZeng0808/IOBR")
```

``` r
library('IOBR')
#> Loading required package: tidyverse
#> -- Attaching packages ----------------------------------------------------------------- tidyverse 1.3.0 --
#> √ ggplot2 3.3.2     √ purrr   0.3.4
#> √ tibble  3.0.1     √ dplyr   1.0.0
#> √ tidyr   1.1.0     √ stringr 1.4.0
#> √ readr   1.3.1     √ forcats 0.5.0
#> -- Conflicts -------------------------------------------------------------------- tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
#> Loading required package: survival
#> Loading required package: survminer
#> Loading required package: ggpubr
#> Loading required package: limma
#> Loading required package: GSVA
#> Loading required package: e1071
#> Loading required package: preprocessCore
#> Warning: replacing previous import 'ggplot2::resolution' by
#> 'limSolve::resolution' when loading 'IOBR'
#> Warning: replacing previous import 'DESeq2::plotMA' by 'limma::plotMA' when
#> loading 'IOBR'
#> Warning: replacing previous import 'dplyr::count' by 'matrixStats::count' when
#> loading 'IOBR'
#> IOBR v1.3.0  For help: https://github.com/DongqiangZeng0808/IOBR
#> 
#>  If you use IOBR in published research, please cite:
#> ==========================================================================
#>  Tumor microenvironment characterization in gastric cancer identifies prognostic
#>  and imunotherapeutically relevant gene signatures.
#>  Cancer Immunology Research, 2019, 7(5), 737-750
#>  DOI: 10.1158/2326-6066.CIR-18-0436 
#>  PMID: 30842092===========================================================================
```

Main documentation is on the `deconvo_tme` function in the package:

  - support methods

<!-- end list -->

``` r
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVM 
#>    "cibersort_abs"              "ips"         "estimate"          "svm_ref" 
#>               lsei 
#>         "lsei_ref"
```

  - quick start

<!-- end list -->

``` r
# Example
lm22<-deconvo_tme(eset = eset_crc, project = "TCGA-CRC", method = "cibersort", arrays = F)
```

subgroup survival analyses

``` r
help("subgroup_survival")
#> starting httpd help server ... done
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
