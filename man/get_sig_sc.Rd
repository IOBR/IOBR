% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_sc.R
\name{get_sig_sc}
\alias{get_sig_sc}
\title{get_sig_sc}
\usage{
get_sig_sc(
  deg,
  cluster = "cluster",
  gene = "gene",
  avg_log2FC = "avg_log2FC",
  n = 100
)
}
\arguments{
\item{deg}{Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)}

\item{cluster}{Name of the column in which the clusters are located}

\item{gene}{Name of the column in which the markers are located}

\item{avg_log2FC}{Name of the column in which the average log2FC values are located}

\item{n}{Number of selected top ranked markers}
}
\value{
A list containing top n gene markers of each cell types
}
\description{
Get gene signature from single-cell differential analysis
}
\examples{
data('deg', package = "IOBR")
get_sig_sc(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100)
}
