% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno_eset.R
\name{anno_eset}
\alias{anno_eset}
\title{Annotating gene expression matrix and remove duplicated genes}
\usage{
anno_eset(
  eset,
  annotation,
  symbol = "symbol",
  probe = "probe_id",
  method = "mean"
)
}
\arguments{
\item{eset}{(Required): An ExpressionSet object containing gene expression data.}

\item{annotation}{(Required): A data.frame that contains annotation information for the probes in the expression set.
user can choose `anno_hug133plus2`, `anno_rnaseq` and `anno_illumina` as input}

\item{symbol}{(Optional, defaults to "symbol"): The column name in the annotation data.frame that represents the gene symbol.}

\item{probe}{(Optional, defaults to "probe_id"): The column name in the annotation data.frame that represents the probe identifiers.}

\item{method}{(Optional, defaults to "mean"): The method used to handle duplicate gene symbols; can be either "mean", "sum" or "sd".}
}
\value{
modified gene expression set
}
\description{
The anno_eset function is used to annotate an ExpressionSet object (eset) with gene symbols using the provided annotation data. Within the function, probes with missing symbols or symbols labelled as "NA_NA" are filtered out. The function calculates and prints the percentage of probes in the expression set that were annotated successfully. It then filters the eset to keep only the rows with probes that have matching identifiers in the annotation data. The function performs additional operations such as merging annotation data with the eset, removing unnecessary columns, transforming rows and columns, handing duplicates based on the specified method. In case of duplicates, gene symbols may be retained based on their mean values (if method is set as "mean") or standard deviation values (if method is set as "sd"). Finally, the function removes rows with all zero values, all NA values, or a NA for the first column. The output is the annotated and cleaned expression set.
}
\examples{
# For affymatrix data
data(eset_gse62254, package = "IOBR")
eset <- anno_eset(eset = eset_gse62254, annotation = anno_hug133plus2)

# For RNAseq data with ensembl id
data(eset_stad, package = "IOBR")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
}
