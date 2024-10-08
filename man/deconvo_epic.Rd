% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvo_tme.R
\name{deconvo_epic}
\alias{deconvo_epic}
\title{Estimating immune microenvironment using EPIC: FOR RNAseq mostly}
\usage{
deconvo_epic(eset, project = NULL, tumor)
}
\arguments{
\item{eset}{expression set with genes at row, sample ID at column}

\item{project}{project name used to distinguish different data sets}

\item{tumor}{is input sample tumor or normal}
}
\value{
EPIC with immune cell fractions
}
\description{
Estimating immune microenvironment using EPIC: FOR RNAseq mostly
}
\examples{
# Loading TCGA-STAD expresion data(raw count matrix)
data(eset_stad, package = "IOBR")
eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
epic_result<-deconvo_epic(eset = eset,project = "TCGA-STAD",tumor = TRUE)

}
\author{
Dongqiang Zeng
}
