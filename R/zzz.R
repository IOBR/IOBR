



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/DongqiangZeng0808/", pkgname, "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, ZL Ye, RF Shen, Y Xiong, JN Wu, WJ Qiu  WJ Liao\n",
                    " IOBR: Multi-omics analysis process integrating tumor microenvironment \n",
                    " and tumor signature. \n",
                    " XXXX, 2020", "\n",
                    " DOI:   ","\n" ,
                    " PMID:  ","\n",
                    "===========================================================================")

  packageStartupMessage(paste0(msg, citation))
}




.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "DongqiagnZeng0808",
    devtools.desc.author = '"Dongqiang Zeng <dognqiangzeng0808@gmail.com> [aut, cre]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  requireNamespace(c('GSVA','limma','DESeq2','tidyverse','MASS',"ggplot2"))

  invisible()
}
