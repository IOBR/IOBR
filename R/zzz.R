



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/DongqiangZeng0808/", pkgname, "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, ZL Ye, RF Shen, Y Xiong, JN Wu, WJ Qiu  WJ Liao\n",
                    " IOBR: Multi-omics analysis process integrating tumor microenvironment \n",
                    " and signatures. \n",
                    # " XXXX, 2020", "\n",
                    # " DOI:   ","\n" ,
                    # " PMID:  ","\n",
                    "===========================================================================")

  packageStartupMessage(paste0(msg, citation))
}




.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2", "ComplexHeatmap",
             "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA","tidyHeatmap"),requireNamespace, quietly = TRUE)))
}
