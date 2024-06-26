% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xCell.R
\name{spillOver}
\alias{spillOver}
\title{Adjust scores using the spillover compensation method}
\usage{
spillOver(transformedScores, K, alpha = 0.5, file.name = NULL)
}
\arguments{
\item{transformedScores}{the trasnformed scores of cell types calculated by transformScores}

\item{K}{the Spillover matrix (spill$K).}

\item{alpha}{a value to override the spillover alpha parameter. Deafult = 0.5}

\item{file.name}{string for the file name for saving the scores. Default is NULL.}
}
\value{
the adjusted xCell scores
}
\description{
spillOver Returns the adjusted xCell scores.
}
