


#
# .onLoad <- function(libname, pkgname) {
#
#   invisible(suppressPackageStartupMessages(library("ComplexHeatmap")))
#   invisible(suppressPackageStartupMessages(library("tidyHeatmap")))
#   invisible(suppressPackageStartupMessages(library("clusterProfiler")))
#
#   invisible(suppressPackageStartupMessages(
#     sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2", "patchwork",
#              "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA"),
#            requireNamespace, quietly = TRUE)
#     ))
# }
#

.onLoad <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap, quietly = TRUE, warn.conflicts = FALSE)
    library(tidyHeatmap, quietly = TRUE, warn.conflicts = FALSE)
    library(clusterProfiler, quietly = TRUE, warn.conflicts = FALSE)
    library(tibble, quietly = TRUE, warn.conflicts = FALSE)
    library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
    library(survival, quietly = TRUE, warn.conflicts = FALSE)
    library(survminer, quietly = TRUE, warn.conflicts = FALSE)
    library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
    library(patchwork, quietly = TRUE, warn.conflicts = FALSE)
    library(ggpubr, quietly = TRUE, warn.conflicts = FALSE)
    library(limma, quietly = TRUE, warn.conflicts = FALSE)
    library(limSolve, quietly = TRUE, warn.conflicts = FALSE)
    library(preprocessCore, quietly = TRUE, warn.conflicts = FALSE)
    library(e1071, quietly = TRUE, warn.conflicts = FALSE)
    library(GSVA, quietly = TRUE, warn.conflicts = FALSE)
    # Add other packages as necessary
  })
}


##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                "  ", pkgname, " v", pkgVersion, "  Immuno-Oncology Biological Research ", "\n",

                "  For Tutorial: https://iobr.github.io/book/", "\n",
                "  For Help: https://github.com/IOBR/IOBR/issues", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",

                    " DQ Zeng, YR Fang, ..., GC Yu*, WJ Liao*, \n",
                    " Enhancing immuno-oncology investigations through multidimensional decoding \n",
                    " of tumor microenvironment with IOBR 2.0. Cell Rep Methods 4, 100910 (2024). \n",
                    " &  \n",
                    " YR Fang, ..., WJ Liao*, DQ Zeng*, \n",
                    " Systematic Investigation of Tumor Microenvironment and , \n",
                    " Antitumor Immunity With IOBR, Med Research (2025). \n",
                    " https://onlinelibrary.wiley.com/doi/epdf/10.1002/mdr2.70001 \n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



