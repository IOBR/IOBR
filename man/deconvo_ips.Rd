% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvo_tme.R
\name{deconvo_ips}
\alias{deconvo_ips}
\title{Calculating immune phenotype score using IPS}
\usage{
deconvo_ips(eset, project = NULL, plot)
}
\arguments{
\item{eset}{expression set with genes at row, sample ID at column}

\item{project}{project name used to distinguish different data sets}
}
\value{
IPS data frame
}
\description{
Calculating immune phenotype score using IPS
}
\examples{
# Loading TCGA-STAD expresion data(raw count matrix)
data(eset_stad, package = "IOBR")
eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
ips_result<-deconvo_ips(eset = eset, project = "TCGA-STAD")

}
\author{
Dongqiang Zeng
}
