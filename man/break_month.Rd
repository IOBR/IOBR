% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv_group.R
\name{break_month}
\alias{break_month}
\title{Break Time Into Blocks}
\usage{
break_month(input, block = 6, time_type = "month")
}
\arguments{
\item{input}{A numeric vector representing time durations that need to be divided, typically in months or days.}

\item{block}{An integer specifying the number of blocks to divide the time into; default is 6.}

\item{time_type}{A character string specifying the units of the input time: "month" for months and "day" for days.
The default is "month". If "day" is specified, the function converts days into months by dividing by 30.}
}
\value{
A numeric vector representing the breakpoints for the time blocks, rounded to the nearest multiple of 5.
}
\description{
Divides a time duration into specified blocks. This function is useful for creating intervals or categories
within a given time period, such as months or days, for further analysis or visualization in studies where
time segmentation might be relevant.
}
\examples{
# Example with time in months
time_data <- c(24, 36, 12, 48)
blocks <- break_month(input = time_data)
print(blocks)

# Example with time in days
day_data <- c(720, 1080, 360, 1440)  # Corresponding to 24, 36, 12, 48 months
blocks_days <- break_month(input = day_data, time_type = "day")
print(blocks_days)
}
\author{
Dongqiang Zeng
}
