% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/high_var_fea.R
\name{high_var_fea}
\alias{high_var_fea}
\title{Extract most variable features form statistical results}
\usage{
high_var_fea(
  result,
  target,
  name_padj = "padj",
  padj_cutoff = 1,
  name_logfc,
  logfc_cutoff = 0,
  n = 10,
  data_type = NULL
)
}
\arguments{
\item{result}{a tibble or data frame}

\item{target}{The column name of the target variable in the dataset.}

\item{name_padj}{The column name representing the adjusted p-value in the dataset. Default is "padj".}

\item{padj_cutoff}{The significance cutoff for adjusted p-values. Values below this threshold are considered significant. Default is 1.}

\item{name_logfc}{he column name representing the log-fold change in the dataset. Default is NULL.}

\item{logfc_cutoff}{The fold change cutoff for log-fold changes. Values above or below this threshold are considered significant. Default is 0.}

\item{n}{The number of top variable features to select. Default is 10.}

\item{data_type}{The type of data being analyzed. Options include "survival" and other relevant data types. Default is NULL.}
}
\value{
A vector of feature names identified as the most variable, combining both upregulated and downregulated features.
}
\description{
This function is designed to identify high variability features (genes or markers) based on certain criteria from a given dataset. It takes as input the dataset and several parameters, such as the significance cutoff for adjusted p-values, the fold change cutoff for log-fold changes, the number of top variables to select, and the data type. The function then filters the dataset based on these criteria and selects the top variable features that meet the specified conditions.
}
\examples{
# Assume 'result_data' is a data frame with statistical analysis results
result_data <- data.frame(
  gene = c("Gene1", "Gene2", "Gene3"),
  padj = c(0.01, 0.2, 0.05),
  logfc = c(-2, 1.5, -3)
)
high_var_features <- high_var_fea(
  result = result_data,
  target = "gene",
  name_padj = "padj",
  name_logfc = "logfc",
  n = 2,
  padj_cutoff = 0.05,
  logfc_cutoff = 1.5
)
print(high_var_features)
}
\author{
Dongqiang Zeng
}
