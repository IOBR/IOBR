% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_sig_score.R
\name{calculate_sig_score_zscore}
\alias{calculate_sig_score_zscore}
\title{Calculating signature score using z-score method}
\usage{
calculate_sig_score_zscore(
  pdata = NULL,
  eset,
  signature,
  mini_gene_count = 3,
  column_of_sample = "ID",
  adjust_eset = FALSE
)
}
\arguments{
\item{pdata}{phenotype data of input sample;
if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset}

\item{eset}{normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)}

\item{signature}{List of gene signatures}

\item{mini_gene_count}{filter out signatures with genes less than minimal gene in expression set;}

\item{column_of_sample}{Defines in which column of pdata the sample identifier can be found}

\item{adjust_eset}{remove variables with missing value, sd =0, and Inf value}
}
\value{
data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
}
\description{
Calculating signature score using z-score method
}
\examples{
# Loading TCGA-STAD expresion data(raw count matrix)
data("eset_stad", package = "IOBR")
# transform count data to tpm
eset <- count2tpm(eset_stad, idType = "ensembl")
# signature score estimation using z-score method
calculate_sig_score_zscore(eset = eset, signature = signature_tme)
}
\author{
Dongqiang Zeng
}
