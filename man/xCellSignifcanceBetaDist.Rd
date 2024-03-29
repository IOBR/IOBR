% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xCell.R
\name{xCellSignifcanceBetaDist}
\alias{xCellSignifcanceBetaDist}
\title{Calculate significance p-values for the null hypothesis that the cell type is not present in the mixture using a random matrix.}
\usage{
xCellSignifcanceBetaDist(
  scores,
  beta_params = NULL,
  rnaseq = T,
  file.name = NULL
)
}
\arguments{
\item{scores}{the xCell scores.}

\item{beta_params}{the pre-calculated beta distribution parameters from random mixtures.}

\item{rnaseq}{if beta_params is null, than uses xCell.data beta_params. If TRUE uses sequencing-based params, else array-based params.}

\item{file.name}{file name to write the p-values table.}
}
\value{
a p-values matrix for each score.
}
\description{
xCellSignifcanceBetaDist Returns the FDR adjusted p-values of the chance that the cell is not present in the mixture.
}
