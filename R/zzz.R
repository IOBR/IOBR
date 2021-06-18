



.onLoad <- function(libname, pkgname) {

  suppressPackageStartupMessages(library("ComplexHeatmap"))
  suppressPackageStartupMessages(library("tidyHeatmap"))

  invisible(suppressPackageStartupMessages(
    sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2",
             "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA"),
           requireNamespace, quietly = TRUE)
    ))
}



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/IOBR/IOBR/issues", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, ZL Ye, RF Sheng, GC Yu, Y Xiong …, WJ Liao*.\n",
                    " IOBR: Multi-omics Immuno-Oncology Biological Research to decode \n ",
                    " tumor microenvironment and signatures. bioRxiv, 2020.2012.2014.422647 (2020). \n",
                    # " XXXX, 2020", "\n",
                    " DOI: 10.1101/2020.12.14.422647\n",
                    # " PMID:  ","\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



