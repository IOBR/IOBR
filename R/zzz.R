#' Package Startup and Initialization
#'
#' @description Internal startup hooks for the IOBR package.
#' The `.onLoad` hook performs non-intrusive checks for optional (suggested)
#' namespaces using `requireNamespace()`, avoiding hard dependencies. The
#' `.onAttach` hook prints a concise startup message that includes the
#' package version, tutorial URL, issue tracker link, and recommended citation
#' information for published use.
#'
#' @details These functions are for internal package initialization and are not
#' intended for direct invocation by end users. The availability checks executed
#' in `.onLoad` are silent and do not attach or load optional packages;
#' they only detect presence so runtime code can adapt accordingly.
#'
#' @param libname Character string giving the library directory where
#'   the package defining the package was found.
#' @param pkgname Character string giving the name of the package.
#'
#' @keywords internal
#' @noRd
#' @importFrom utils packageDescription

.onLoad <- function(libname, pkgname) {
  # # Check that suggested packages are available (silently)
  # # These are optional dependencies used by various functions
  # suggested_pkgs <- c(
  #   "ComplexHeatmap", "tidyHeatmap", "clusterProfiler", "tibble",
  #   "tidyverse", "survival", "survminer", "ggplot2",
  #   "ggpubr", "limma", "limSolve", "preprocessCore", "e1071", "GSVA"
  # )

  # invisible(suppressPackageStartupMessages(
  #   vapply(
  #     suggested_pkgs,
  #     requireNamespace,
  #     logical(1),
  #     quietly = TRUE
  #   )
  # ))
}


.onAttach <- function(libname, pkgname) {
  pkgVersion <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0(
    "==========================================================================\n",
    "  ", pkgname, " v", pkgVersion, "  Immuno-Oncology Biological Research ", "\n",
    "  For Documentation: https://iobr.github.io/IOBR/", "\n",
    "  For Tutorial: https://iobr.github.io/book/", "\n",
    "  For Help: https://github.com/IOBR/IOBR/issues", "\n\n"
  )

  citation <- paste0(
    " If you use ", pkgname, " in published research, please cite:\n",
    " DQ Zeng, YR Fang, ..., GC Yu*, WJ Liao*, \n",
    " Enhancing immuno-oncology investigations through multidimensional decoding \n",
    " of tumor microenvironment with IOBR 2.0. Cell Rep Methods 4, 100910 (2024). \n",
    " \u0026  \n",
    " YR Fang, ..., WJ Liao*, DQ Zeng*, \n",
    " Systematic Investigation of Tumor Microenvironment and \n",
    " Antitumor Immunity With IOBR, Med Research (2025). \n",
    " https://onlinelibrary.wiley.com/doi/epdf/10.1002/mdr2.70001 \n",
    "=========================================================================="
  )

  packageStartupMessage(msg, citation)
}
