



.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(
    sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2", "ComplexHeatmap",
             "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA","tidyHeatmap"),
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
                    " Zeng et al., (2020)*\n",
                    " IOBR: Multi-omics Immuno-Oncology Biological Research to decode tumor microenvironment \n",
                    " and signatures. \n",
                    # " XXXX, 2020", "\n",
                    # " DOI:   ","\n" ,
                    # " PMID:  ","\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



