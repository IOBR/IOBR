\name{iobr_ai_query}
\alias{iobr_ai_query}
\title{Query the AI assistant}
\usage{
iobr_ai_query(
  query,
  index_path = NULL,
  provider = NULL,
  top_k = 5,
  max_tokens = 800,
  temperature = 0.2
)
}
\arguments{
\item{query}{Natural language query}

\item{index_path}{Path to index RDS file. If NULL, uses default location
in inst/ai/iobr_embeddings.rds}

\item{provider}{Provider configuration (same format as iobr_ai_init). If NULL,
must be provided via environment variables.}

\item{top_k}{Number of top relevant chunks to retrieve (default: 5)}

\item{max_tokens}{Maximum tokens to generate in response (default: 800)}

\item{temperature}{Sampling temperature for generation (default: 0.2)}
}
\value{
List containing:
\itemize{
  \item answer: Generated answer text
  \item code: Extracted R code snippets (if any)
  \item sources: List of source documents with similarity scores
  \item raw_response: Raw API response
}
}
\description{
Searches the index for relevant content and generates a response using
the configured AI provider.
}
\examples{
\dontrun{
# Query the assistant
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)
result <- iobr_ai_query(
  "How do I perform TME deconvolution?",
  provider = provider
)
cat(result$answer)
}
}
