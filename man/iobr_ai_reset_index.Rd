\name{iobr_ai_reset_index}
\alias{iobr_ai_reset_index}
\title{Reset/delete AI index}
\usage{
iobr_ai_reset_index(index_path = NULL)
}
\arguments{
\item{index_path}{Path to index RDS file. If NULL, uses default location.}
}
\value{
TRUE if successful
}
\description{
Removes the index file from disk.
}
\examples{
\dontrun{
iobr_ai_reset_index()
}
}
