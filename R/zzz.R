#' Package Startup and Initialization
#'
#' @description Internal startup hooks for the IOBR package.
#' The \code{.onLoad} hook performs non-intrusive checks for optional (suggested)
#' namespaces using \code{requireNamespace()}, avoiding hard dependencies. The
#' \code{.onAttach} hook prints a concise startup message that includes the
#' package version, tutorial URL, issue tracker link, and recommended citation
#' information for published use.
#'
#' @details These functions are for internal package initialization and are not
#' intended for direct invocation by end users. The availability checks executed
#' in \code{.onLoad} are silent and do not attach or load optional packages;
#' they only detect presence so runtime code can adapt accordingly.
#'
#' @keywords internal
#' @noRd
#' @importFrom utils packageDescription

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
  # Check that required packages are available
  invisible(suppressPackageStartupMessages({
    sapply(
      c(
        "ComplexHeatmap", "tidyHeatmap", "clusterProfiler", "tibble",
        "tidyverse", "survival", "survminer", "ggplot2", "patchwork",
        "ggpubr", "limma", "limSolve", "preprocessCore", "e1071", "GSVA"
      ),
      requireNamespace,
      quietly = TRUE
    )
  }))
}


.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  msg <- paste0(
    "==========================================================================\n",
    "  ", pkgname, " v", pkgVersion, "  Immuno-Oncology Biological Research ", "\n",
    "  For Tutorial: https://iobr.github.io/book/", "\n",
    "  For Help: https://github.com/IOBR/IOBR/issues", "\n\n"
  )

  citation <- paste0(
    " If you use ", pkgname, " in published research, please cite:\n",
    " DQ Zeng, YR Fang, ..., GC Yu*, WJ Liao*, \n",
    " Enhancing immuno-oncology investigations through multidimensional decoding \n",
    " of tumor microenvironment with IOBR 2.0. Cell Rep Methods 4, 100910 (2024). \n",
    " &  \n",
    " YR Fang, ..., WJ Liao*, DQ Zeng*, \n",
    " Systematic Investigation of Tumor Microenvironment and \n",
    " Antitumor Immunity With IOBR, Med Research (2025). \n",
    " https://onlinelibrary.wiley.com/doi/epdf/10.1002/mdr2.70001 \n",
    "=========================================================================="
  )

  packageStartupMessage(paste0(msg, citation))
}
