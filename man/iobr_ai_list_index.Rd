\name{iobr_ai_list_index}
\alias{iobr_ai_list_index}
\title{List index contents}
\usage{
iobr_ai_list_index(index_path = NULL)
}
\arguments{
\item{index_path}{Path to index RDS file. If NULL, uses default location.}
}
\value{
Data frame with index statistics
}
\description{
Shows information about the current index including number of entries
and sources.
}
\examples{
\dontrun{
index_info <- iobr_ai_list_index()
print(index_info)
}
}
