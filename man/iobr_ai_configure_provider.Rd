\name{iobr_ai_configure_provider}
\alias{iobr_ai_configure_provider}
\title{Configure and validate AI provider settings}
\usage{
iobr_ai_configure_provider(provider)
}
\arguments{
\item{provider}{List containing provider configuration}
}
\value{
Validated provider configuration
}
\description{
Validates provider configuration and sets defaults for known providers.
}
\examples{
\dontrun{
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))
}
}
