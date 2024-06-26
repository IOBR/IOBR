% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/percent_bar.R
\name{pie_chart}
\alias{pie_chart}
\title{pie_chart}
\usage{
pie_chart(
  input,
  var,
  var2 = NULL,
  type = 2,
  show_freq = FALSE,
  color = NULL,
  palette = "jama",
  title = NULL,
  text_size = 10,
  title_size = 20,
  add_sum = FALSE
)
}
\arguments{
\item{input}{The input dataframe.}

\item{var}{The variable on which the pie chart or donut chart will be based.}

\item{var2}{Optional. A secondary variable for creating a donut chart. Default is NULL.}

\item{type}{The type of chart to be generated. 1 for pie chart, 2 for donut chart, 3 for donut chart based on webr package.}

\item{show_freq}{Boolean indicating whether to show frequencies on the chart. Default is FALSE.}

\item{color}{Optional. The color palette for the chart. If not provided, a default color palette will be used.}

\item{palette}{Optional. The color palette to be used if color is not specified. Default is "jama".}

\item{title}{Optional. The title of the chart. Default is NULL.}

\item{text_size}{The size of the text on the chart. Default is 10.}

\item{title_size}{The size of the title text on the chart. Default is 20.}

\item{add_sum}{Boolean indicating whether to add the sum of frequencies to the chart. Default is FALSE.}
}
\value{
Generate Pie or Donut Charts
}
\description{
This function generates a pie chart or a donut chart from the input data. It allows customization of various visual aspects such as colors, labels, and title.
}
\examples{
pie_chart(input = sig_stad, var = "Subtype",  palette = "jama")
pie_chart(input = sig_stad, var = "Subtype",  palette = "nrc")
}
\author{
Dongqiang Zeng
}
