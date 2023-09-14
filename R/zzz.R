



.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(library("ComplexHeatmap")))
  invisible(suppressPackageStartupMessages(library("tidyHeatmap")))
  invisible(suppressPackageStartupMessages(library("clusterProfiler")))

  invisible(suppressPackageStartupMessages(
    sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2", "patchwork",
             "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA"),
           requireNamespace, quietly = TRUE)
    ))
}



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                "  ", pkgname, " v", pkgVersion, "  Immuno-Oncology Biological Research ", "\n",

                "  For Tutorial: https://iobr.github.io/book/", "\n",
                "  For Help: https://github.com/IOBR/IOBR/issues", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, ZL Ye, RF Sheng, GC Yu, Y Xiong …, WJ Liao*.\n",
                    " IOBR: Multi-omics Immuno-Oncology Biological Research to decode \n ",
                    " tumor microenvironment and signatures. Frontiers in Immunology. 12:687975,(2021). \n",
                    # " XXXX, 2020", "\n",
                    " DOI: 10.3389/fimmu.2021.687975\n",
                    " Higly Cited Paper and Hot Paper of WOS\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



