\name{iobr_ai_init}
\alias{iobr_ai_init}
\title{Initialize AI index from package content}
\usage{
iobr_ai_init(
  pkg_path = ".",
  provider = list(name = "openai", api_key = Sys.getenv("OPENAI_API_KEY"),
    base_url = NULL, model_embeddings = NULL, model_chat = NULL),
  index_path = NULL,
  chunk_size = 800,
  progress_callback = NULL
)
}
\arguments{
\item{pkg_path}{Path to package root directory (default: current directory)}

\item{provider}{List containing provider configuration:
\itemize{
  \item name: Provider name ('openai', 'anthropic', 'huggingface', or 'custom')
  \item api_key: API key for authentication (can use Sys.getenv())
  \item base_url: Base URL for API (optional, uses defaults for known providers)
  \item model_embeddings: Model name for embeddings (optional, uses defaults)
  \item model_chat: Model name for chat/generation (optional, uses defaults)
  \item headers: Named list of additional HTTP headers (optional)
}}

\item{index_path}{Path to save index RDS file. If NULL, uses
inst/ai/iobr_embeddings.rds in the package.}

\item{chunk_size}{Maximum characters per text chunk (default: 800)}

\item{progress_callback}{Optional function to receive progress messages}
}
\value{
Path to saved index file (invisibly)
}
\description{
Creates an index of package documentation and code with embeddings for
semantic search. The index is saved to disk for later use.
}
\examples{
\dontrun{
# Initialize with OpenAI
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)
iobr_ai_init(provider = provider)

# Initialize with custom provider
provider <- list(
  name = "custom",
  api_key = "your-key",
  base_url = "https://your-api.com/v1/embeddings",
  model_embeddings = "your-model"
)
iobr_ai_init(provider = provider)
}
}
