% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_pheatmap.R
\name{sig_pheatmap}
\alias{sig_pheatmap}
\title{Title sig_pheatmap}
\usage{
sig_pheatmap(
  input,
  feas,
  group,
  group2 = NULL,
  group3 = NULL,
  ID = "ID",
  path = NULL,
  cols1 = "random",
  cols2 = "random",
  cols3 = "random",
  seed = 54321,
  show_col = FALSE,
  palette1 = 1,
  palette2 = 2,
  palette3 = 3,
  cluster_cols = TRUE,
  palette_for_heatmape = 6,
  scale.matrix = TRUE,
  cellwidth = 1,
  cellheight = 9,
  show_colnames = FALSE,
  fig.type = "pdf",
  width = 6,
  height = NULL,
  file_name_prefix = 1
)
}
\arguments{
\item{input}{This parameter represents the input data for the heatmap. It should be a data frame with variables in column}

\item{feas}{This parameter specifies the features to include in the heatmap. It should be a vector containing the names of the columns in the "input" data frame.}

\item{group}{This parameter indicates the grouping variable for the heatmap. It should be a column name in the "input" data frame.}

\item{group2}{(optional): This parameter represents an additional grouping variable for the heatmap. It should be a column name in the "input" data frame.}

\item{group3}{(optional): This parameter indicates another additional grouping variable for the heatmap. It should be a column name in the "input" data frame.}

\item{ID}{(default: "ID"): This parameter allows you to specify the column name to be used as the sample identifier. The default value is "ID".}

\item{path}{(optional): This parameter represents the path where the output files will be stored. If not provided, the default path will be used.}

\item{cols1}{(default: "random"): This parameter specifies the colors to be used for the first grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".}

\item{cols2}{(default: "random"): This parameter represents the colors to be used for the second grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".}

\item{cols3}{(default: "random"): This parameter indicates the colors to be used for the third grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".}

\item{seed}{(default: 54321): This parameter specifies the seed for the random number generator used in color selection. The default value is 54321.}

\item{show_col}{(default: FALSE): This parameter determines whether to display the colors used in the heatmap. The default value is FALSE.}

\item{palette1}{(default: 1): This parameter represents the palette number for grouping variable 1. The default value is 1.}

\item{palette2}{(default: 2): This parameter indicates the palette number for grouping variable 2. The default value is 2.}

\item{palette3}{(default: 3): This parameter specifies the palette number for grouping variable 3. The default value is 3.}

\item{cluster_cols}{(default: TRUE): This parameter determines whether to cluster the columns of the heatmap. The default value is TRUE.}

\item{palette_for_heatmape}{(default: 6): This parameter represents the palette number for the heatmap. The default value is 6.}

\item{scale.matrix}{(default: TRUE): This parameter specifies whether to scale the input matrix. The default value is TRUE.}

\item{cellwidth}{(default: 9): This parameter determines the width (in points) of each cell in the heatmap. The default value is 9.}

\item{cellheight}{(default: 9): This parameter specifies the height (in points) of each cell in the heatmap. The default value is 9.}

\item{show_colnames}{(default: FALSE): This parameter determines whether to display colum names}

\item{fig.type}{(default: "pdf"): This parameter indicates the file format for saving the heatmap. The default value is "pdf".}

\item{width}{(default: 9): This parameter represents the width (in inches) of the saved heatmap. The default value is 9.}

\item{height}{(optional): This parameter specifies the height (in inches) of the saved heatmap. If not provided, a default height will be calculated based on the number of features.}

\item{file_name_prefix}{(default: 1): This parameter allows you to specify a prefix for the saved file name. The default value is 1.}
}
\value{
A list containing the annotation data, cluster colors, the plot object, and the transformed input matrix.
}
\description{
The sig_pheatmap function is used to generate heatmaps. It takes a data frame as input and creates a heatmap with grouping variables based on the provided parameters. The parameters include the input data, the features to include in the heatmap, the grouping variable, and optional additional grouping variables. It also offers flexible options to specify colors, adjust the size and layout of the heatmap, and save output files. The function returns a list containing annotation data, cluster colors, the plot object, and the transformed input matrix, allowing users to further analyze and visualize the heatmap data.
}
\examples{
data("tcga_stad_sig", package = "IOBR")
data("tcga_stad_pdata", package = "IOBR")
input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
feas <- colnames(input)[grep(colnames(input), pattern = "MCPcounter")]
sig_pheatmap(input = input, feas = feas , group = "subtype", scale.matrix = TRUE)

}
\author{
Dongqiang Zeng
}
