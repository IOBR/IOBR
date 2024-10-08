% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_heatmap.R
\name{sig_heatmap}
\alias{sig_heatmap}
\title{sig_heatmap}
\usage{
sig_heatmap(
  input,
  ID = "ID",
  features,
  group,
  condiction = NULL,
  id_condiction = "vars",
  col_condiction = "condiction",
  cols_condiction = NULL,
  scale = FALSE,
  palette = 2,
  palette_group = "jama",
  show_col = F,
  show_palettes = F,
  cols_group = NULL,
  show_plot = T,
  width = 8,
  height = NULL,
  size_col = 10,
  size_row = 8,
  angle_col = 90,
  column_title = NULL,
  row_title = NULL,
  show_heatmap_col_name = F,
  path = NULL,
  index = NULL
)
}
\arguments{
\item{input}{data frame with `ID`, `group` and features}

\item{ID}{This parameter specifies the column name in the input data that represents the unique identifier. The default value is "ID".}

\item{features}{This parameter is a vector of column names representing the features to be used in the heatmap.}

\item{group}{This parameter is a column name in the input data that represents the grouping variable.}

\item{condiction}{This parameter is an optional data frame that defines additional conditions to be applied in the heatmap. It should have two columns: one specifying the variables and another specifying the conditions. The default value is NULL.}

\item{id_condiction}{This parameter specifies the column name in the condiction data frame that represents the variables. The default value is "vars".}

\item{col_condiction}{This parameter specifies the column name in the condiction data frame that represents the conditions. The default value is "condiction".}

\item{cols_condiction}{This parameter is a vector of colors to be used for the condiction variable.}

\item{scale}{This parameter specifies whether to scale the heatmap by column. If set to TRUE, the heatmap will be scaled by column. The default value is FALSE.}

\item{palette}{This parameter specifies the palette to be used in the heatmap. It can be a numeric value representing the palette index or a string specifying the palette name. The default value is 2.}

\item{palette_group}{This parameter specifies the palette group to be used in the heatmap. It should be a string representing the palette group name. The default value is "jama".}

\item{show_col}{This parameter specifies whether to show the color boxes in the heatmap. If set to TRUE, the color boxes will be shown. The default value is FALSE.}

\item{show_palettes}{This parameter specifies whether to show the available palettes. If set to TRUE, the available palettes will be shown. The default value is FALSE.}

\item{cols_group}{This parameter is a vector of colors to be used for the grouping variable. If not specified, the colors will be randomly assigned.}

\item{show_plot}{This parameter specifies whether to show the resulting heatmap plot. If set to TRUE, the plot will be shown. The default value is TRUE.}

\item{width}{This parameter specifies the width of the resulting heatmap plot in inches. The default value is 8.}

\item{height}{This parameter specifies the height of the resulting heatmap plot in inches. If not specified, the height will be calculated based on the number of features.}

\item{size_col}{This parameter specifies the font size of the column names in the heatmap plot. The default value is 10.}

\item{size_row}{This parameter specifies the font size of the row names in the heatmap plot. The default value is 8.}

\item{angle_col}{This parameter specifies the rotation angle of the column names in the heatmap plot. The default value is 90.}

\item{column_title}{This parameter specifies the title of the column in the heatmap plot. If not specified, no title will be displayed.}

\item{row_title}{This parameter specifies the title of the row in the heatmap plot. If not specified, no title will be displayed.}

\item{show_heatmap_col_name}{This parameter specifies whether to show the column names in the heatmap plot. If set to TRUE, the column names will be shown. The default value is FALSE.}

\item{path}{This parameter specifies the directory path where the resulting heatmap plot will be saved. If not specified, a default path will be used.}

\item{index}{default is nullThis parameter specifies the index number to be added to the filename of the resulting heatmap plot. If not specified, the index will be set to 1.}
}
\value{
A heatmap plot object.
}
\description{
The sig_heatmap function is used to generate a heatmap plot based on input data, grouping variables, and optional conditions. The function allows customization of various parameters such as palette selection, scaling, color boxes, plot dimensions, and more. It provides flexibility in visualizing relationships between variables and groups in a concise and informative manner.
}
\examples{
data("tcga_stad_sig", package = "IOBR")
data("tcga_stad_pdata", package = "IOBR")
input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
feas <- colnames(input)[grep(colnames(input), pattern = "MCPcounter")]
sig_heatmap(input = input, features = feas , group = "subtype", scale = TRUE)

}
\author{
Dongqiang Zeng
}
