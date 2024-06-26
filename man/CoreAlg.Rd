% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIBERSORT.R
\name{CoreAlg}
\alias{CoreAlg}
\title{Title CoreAlg}
\usage{
CoreAlg(X, y, absolute, abs_method)
}
\arguments{
\item{X}{Input matrix or data frame containing the predictor variables.}

\item{y}{Numeric vector containing the response variable.}

\item{absolute}{Logical value indicating whether to use absolute space or relative space for the weights.}

\item{abs_method}{String indicating the method to calculate the weights in absolute space. Can be either 'sig.score' or 'no.sumto1'.}
}
\value{
A list containing the weights (`w`), root mean squared error (`mix_rmse`), and correlation coefficient (`mix_r`).
}
\description{
CoreAlg performs nu-regression using support vector machines (SVM) and calculates weights, root mean squared error (RMSE), and correlation coefficient (R).
}
\examples{
X <- matrix(rnorm(100), nrow = 10)
y <- rnorm(10)
result <- CoreAlg(X, y, absolute = FALSE, abs_method = "sig.score")

}
