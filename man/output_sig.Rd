% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_sig.R
\name{output_sig}
\alias{output_sig}
\title{save signature data into a data frame}
\usage{
output_sig(signatures, format = "csv", file.name)
}
\arguments{
\item{signatures}{output signatures, a list or a string}

\item{format}{output format}

\item{file.name}{file name of output data}
}
\value{
Returns the data frame containing the processed signature data. This data frame is also saved to the specified file format.
}
\description{
This function saves the provided signature data into a specified file format. It supports saving as CSV or RData format. 
The function can handle single signature strings or lists of signatures, converting them into a data frame for storage.
}
\examples{
output_sig(signature_collection, file.name = "signature")
}
