% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_msigdb.R
\name{format_msigdb}
\alias{format_msigdb}
\title{Format input signatures from MsgiDB}
\usage{
format_msigdb(gmt, ont = "ont", gene = "gene")
}
\arguments{
\item{ont}{column name of signatures}

\item{gene}{column name of gene}

\item{sig}{`signature` data frame downloaded from `https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H`}
}
\value{
signature list
}
\description{
Format input signatures from MsgiDB
}
