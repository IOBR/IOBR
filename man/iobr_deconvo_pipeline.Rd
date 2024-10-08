% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iobr_deconvo_pipeline.R
\name{iobr_deconvo_pipeline}
\alias{iobr_deconvo_pipeline}
\title{Tumor Microenvironment Analysis Pipeline}
\usage{
iobr_deconvo_pipeline(
  eset,
  project,
  array,
  tumor_type,
  path = "1-TME",
  permutation = 1000
)
}
\arguments{
\item{eset}{A gene expression dataset, typically an expression set object, prepared for TME analysis.}

\item{project}{A character string specifying the project or analysis name, used for output file naming.}

\item{array}{A logical indicating whether the data comes from an array platform; influences deconvolution methods.}

\item{tumor_type}{A character string specifying the tumor type, which can tailor certain analysis aspects.}

\item{path}{A string indicating the output file path; defaults to "1-TME".}

\item{permutation}{An integer specifying the number of permutations to use in TME deconvolution methods; default is 1000.}
}
\value{
A comprehensive dataset combining TME cell fractions, signature scores, and integrated TME-signature analysis.
        This function saves several files to the specified path, documenting various stages of the analysis.
}
\description{
The `iobr_deconvo_pipeline` function executes a comprehensive TME analysis on a gene expression dataset.
This analysis includes TME deconvolution using various computational methods, calculation of signature scores,
and integration of these data into a cohesive output. It is designed for in-depth exploration of the
microenvironment's role in tumor biology based on gene expression profiles.
}
\examples{
data("eset_stad", package = "IOBR")
eset <- count2tpm(eset_stad)  # Prepare data
# Run the pipeline
res <- iobr_deconvo_pipeline(eset = eset, project = "STAD", array = FALSE, tumor_type = "stad", path = "1-TME", permutation = 1000)
}
\author{
Dongqiang Zeng
}
