



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",

                "For help: https://github.com/DongqiangZeng0808/", pkgname, "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    "==========================================================================\n",
                    " Tumor microenvironment characterization in gastric cancer identifies prognostic\n",
                    " and imunotherapeutically relevant gene signatures.","\n",
                    " Cancer Immunology Research, 2019, 7(5), 737-750", "\n",
                    " DOI: 10.1158/2326-6066.CIR-18-0436 ","\n" ,
                    " PMID: 30842092",
                    "===========================================================================")

  packageStartupMessage(paste0(msg, citation))
}

