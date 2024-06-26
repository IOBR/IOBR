% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_surv_plot.R
\name{sig_surv_plot}
\alias{sig_surv_plot}
\title{Multiple survival KM-plot to a signature or gene}
\usage{
sig_surv_plot(
  input_pdata,
  signature,
  project = "KM",
  ID = "ID",
  time = "time",
  status = "status",
  time_type = "month",
  break_month = "auto",
  cols = NULL,
  palette = "jama",
  show_col = TRUE,
  mini_sig = "score",
  fig.type = "png",
  save_path = "KM-plot",
  index = 1
)
}
\arguments{
\item{input_pdata}{This parameter represents the input data for the survival analysis. It is expected to be a data frame.}

\item{signature}{This parameter specifies the column name in input_pdata that represents the target signature for survival analysis.}

\item{project}{This optional parameter represents the project name. The default value is "KM".}

\item{ID}{This optional parameter represents the column name in input_pdata that contains the unique identifier for each data point. The default value is "ID".}

\item{time}{This optional parameter represents the column name in input_pdata that contains the time variable for survival analysis. The default value is "time".}

\item{status}{This optional parameter represents the column name in input_pdata that contains the survival status variable. The default value is "status".}

\item{time_type}{This optional parameter specifies the time unit used in the analysis. The default value is "month".}

\item{break_month}{This optional parameter specifies the interval at which the time axis breaks should be made in the Kaplan-Meier plot. The default value is "auto".}

\item{cols}{This optional parameter allows the user to specify the color palette for the plot. If not specified, a default palette will be used.}

\item{palette}{This optional parameter allows the user to specify the palette type if cols is not provided. The default value is "jama".}

\item{show_col}{logical variable, if TRUE, color will be print and show in the R studio}

\item{mini_sig}{This optional parameter represents the label for the mini signature in the legend. The default value is "score".}

\item{fig.type}{This optional parameter specifies the file format for saving the plot. The default value is "png".}

\item{save_path}{This optional parameter specifies the directory path for saving the plot. The default value is "KM-plot".}

\item{index}{This parameter represents the index for multiple plots generated in a single run. It is used as the prefix in the filename. The default value is 1.}
}
\value{
survival
}
\description{
Multiple survival KM-plot to a signature or gene
}
\examples{

data("tcga_stad_pdata", package = "IOBR")
sig_surv_plot(input_pdata = tcga_stad_pdata, signature = "TMEscore_plus", time = "time", status = "OS_status")

}
\author{
Dongqiang Zeng
}
