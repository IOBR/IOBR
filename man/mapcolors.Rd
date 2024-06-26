% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPS_helper.R
\name{mapcolors}
\alias{mapcolors}
\title{Map Colors Based on Input Value}
\usage{
mapcolors(x)
}
\arguments{
\item{x}{A numeric value to be mapped to a color.}
}
\value{
A color from the predefined palette corresponding to the input value.
}
\description{
This function maps a numeric input value to a color from a predefined palette.
}
\examples{
# Define a palette with 1001 colors
my_palette <- colorRampPalette(c("blue", "white", "red"))(1001)
# Map a value to a color
color <- mapcolors(2)
print(color)
}
