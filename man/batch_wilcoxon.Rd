% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_wilcoxon.R
\name{batch_wilcoxon}
\alias{batch_wilcoxon}
\title{Batch to conduct wilcoxon test between two groups}
\usage{
batch_wilcoxon(
  data,
  target = "group",
  feature = NULL,
  feature_manipulation = FALSE
)
}
\arguments{
\item{data}{A data frame containing the input data for analysis.}

\item{target}{The name of the column in the data frame that represents the group labels. The default value is "group".}

\item{feature}{A character vector specifying the names of the features or variables to be analyzed. If NULL, all continuous features in the data frame will be used. Default value is NULL.}

\item{feature_manipulation}{A logical value indicating whether feature manipulation is required. If TRUE, a custom feature manipulation function will be applied to the data. Default value is FALSE.}
}
\value{
statistical result
}
\description{
The batch_wilcoxon function performs Wilcoxon rank-sum tests on a given data set to compare the distribution of a specified feature between two groups. It computes the p-values and ranks the significant features based on the p-values. The function returns a data frame with the feature names, p-values, adjusted p-values, logarithm of p-values, and a star rating based on the p-value ranges.
}
\examples{
# Loading TCGA-STAD microenvironment signature data
data("sig_stad", package = "IOBR")
# Finding microenvironmental scores associated with Gender
batch_wilcoxon(data = sig_stad, target = "Gender", feature = colnames(sig_stad)[69:ncol(sig_stad)])
}
\author{
Dongqiang Zeng
}
