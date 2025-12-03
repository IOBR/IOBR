#' IOBR AI Assistant Functions
#'
#' Functions for interacting with the IOBR AI assistant that provides
#' Retrieval-Augmented Generation (RAG) capabilities for package documentation
#' and code.
#'
#' @name iobr_ai
NULL

#' Initialize AI index from package content
#'
#' Creates an index of package documentation and code with embeddings for
#' semantic search. The index is saved to disk for later use.
#'
#' @param pkg_path Path to package root directory (default: current directory)
#' @param provider List containing provider configuration:
#'   \itemize{
#'     \item name: Provider name ('openai', 'anthropic', 'huggingface', or 'custom')
#'     \item api_key: API key for authentication (can use Sys.getenv())
#'     \item base_url: Base URL for API (optional, uses defaults for known providers)
#'     \item model_embeddings: Model name for embeddings (optional, uses defaults)
#'     \item model_chat: Model name for chat/generation (optional, uses defaults)
#'     \item headers: Named list of additional HTTP headers (optional)
#'   }
#' @param index_path Path to save index RDS file. If NULL, uses
#'   inst/ai/iobr_embeddings.rds in the package.
#' @param chunk_size Maximum characters per text chunk (default: 800)
#' @param progress_callback Optional function to receive progress messages
#'
#' @return Path to saved index file (invisibly)
#'
#' @examples
#' \dontrun{
#' # Initialize with OpenAI
#' provider <- list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#' iobr_ai_init(provider = provider)
#'
#' # Initialize with custom provider
#' provider <- list(
#'   name = "custom",
#'   api_key = "your-key",
#'   base_url = "https://your-api.com/v1/embeddings",
#'   model_embeddings = "your-model"
#' )
#' iobr_ai_init(provider = provider)
#' }
#'
#' @export
iobr_ai_init <- function(pkg_path = ".", 
                        provider = list(
                          name = "openai",
                          api_key = Sys.getenv("OPENAI_API_KEY"),
                          base_url = NULL,
                          model_embeddings = NULL,
                          model_chat = NULL
                        ),
                        index_path = NULL,
                        chunk_size = 800,
                        progress_callback = NULL) {
  
  # Validate provider
  validated_provider <- iobr_ai_configure_provider(provider)
  
  # Set default index path
  if (is.null(index_path)) {
    index_path <- file.path(pkg_path, "inst", "ai", "iobr_embeddings.rds")
  }
  
  # Collect package content
  if (!is.null(progress_callback)) {
    progress_callback("Collecting package content...")
  }
  content_items <- collect_package_content(pkg_path)
  
  if (length(content_items) == 0) {
    stop("No content found in package")
  }
  
  # Create index with embeddings
  if (!is.null(progress_callback)) {
    progress_callback(sprintf("Creating embeddings for %d content items...", 
                            length(content_items)))
  }
  
  index <- create_index(content_items, validated_provider, chunk_size, progress_callback)
  
  # Save index
  if (!is.null(progress_callback)) {
    progress_callback("Saving index...")
  }
  save_index(index, index_path)
  
  if (!is.null(progress_callback)) {
    progress_callback(sprintf("Index saved to %s", index_path))
  }
  
  message("AI index initialized successfully with ", length(index), " entries")
  invisible(index_path)
}

#' Query the AI assistant
#'
#' Searches the index for relevant content and generates a response using
#' the configured AI provider.
#'
#' @param query Natural language query
#' @param index_path Path to index RDS file. If NULL, uses default location
#'   in inst/ai/iobr_embeddings.rds
#' @param provider Provider configuration (same format as iobr_ai_init). If NULL,
#'   must be provided via environment variables.
#' @param top_k Number of top relevant chunks to retrieve (default: 5)
#' @param max_tokens Maximum tokens to generate in response (default: 800)
#' @param temperature Sampling temperature for generation (default: 0.2)
#'
#' @return List containing:
#'   \itemize{
#'     \item answer: Generated answer text
#'     \item code: Extracted R code snippets (if any)
#'     \item sources: List of source documents with similarity scores
#'     \item raw_response: Raw API response
#'   }
#'
#' @examples
#' \dontrun{
#' # Query the assistant
#' provider <- list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#' result <- iobr_ai_query(
#'   "How do I perform TME deconvolution?",
#'   provider = provider
#' )
#' cat(result$answer)
#' }
#'
#' @export
iobr_ai_query <- function(query, 
                         index_path = NULL,
                         provider = NULL,
                         top_k = 5,
                         max_tokens = 800,
                         temperature = 0.2) {
  
  # Validate provider
  if (is.null(provider)) {
    stop("Provider configuration is required")
  }
  validated_provider <- iobr_ai_configure_provider(provider)
  
  # Set default index path
  if (is.null(index_path)) {
    # Try to find in inst/ai or package installation
    pkg_inst_path <- system.file("ai", "iobr_embeddings.rds", package = "IOBR")
    if (file.exists(pkg_inst_path)) {
      index_path <- pkg_inst_path
    } else {
      index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
    }
  }
  
  # Load index
  if (!file.exists(index_path)) {
    stop("Index file not found. Please run iobr_ai_init() first.")
  }
  index <- load_index(index_path)
  
  # Search for relevant chunks
  search_results <- search_index(query, index, validated_provider, top_k)
  
  # Construct prompt with context
  context_text <- paste(
    sapply(search_results, function(r) {
      sprintf("Source: %s (similarity: %.3f)\n%s", r$source, r$similarity, r$text)
    }),
    collapse = "\n\n---\n\n"
  )
  
  system_prompt <- paste(
    "You are an expert assistant for the IOBR R package for immuno-oncology analysis.",
    "Answer questions based on the provided documentation and code context.",
    "If you provide R code examples, wrap them in ```r code blocks.",
    "Always cite the sources you use in your answer.",
    "If the answer is not in the provided context, say so clearly."
  )
  
  user_prompt <- sprintf(
    "Context from IOBR documentation:\n\n%s\n\n---\n\nQuestion: %s\n\nPlease provide a clear answer with code examples if applicable, and cite your sources.",
    context_text,
    query
  )
  
  # Get response from AI
  response <- send_chat(
    system = system_prompt,
    user = user_prompt,
    provider_config = validated_provider,
    max_tokens = max_tokens,
    temperature = temperature
  )
  
  # Extract code blocks
  code_pattern <- "```r\\s*\n([^`]+)\n```"
  code_matches <- gregexpr(code_pattern, response$content, perl = TRUE)
  code_snippets <- character()
  
  if (code_matches[[1]][1] != -1) {
    for (i in seq_along(code_matches[[1]])) {
      match_start <- code_matches[[1]][i]
      match_length <- attr(code_matches[[1]], "match.length")[i]
      code_block <- substr(response$content, match_start, match_start + match_length - 1)
      
      # Extract just the code
      code_only <- sub("```r\\s*\n", "", code_block)
      code_only <- sub("\n```$", "", code_only)
      code_snippets <- c(code_snippets, code_only)
    }
  }
  
  # Build result
  result <- list(
    answer = response$content,
    code = code_snippets,
    sources = lapply(search_results, function(r) {
      list(
        source = r$source,
        type = r$type,
        similarity = r$similarity,
        text = r$text
      )
    }),
    raw_response = response$raw
  )
  
  return(result)
}

#' Configure and validate AI provider settings
#'
#' Validates provider configuration and sets defaults for known providers.
#'
#' @param provider List containing provider configuration
#'
#' @return Validated provider configuration
#'
#' @examples
#' \dontrun{
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' ))
#' }
#'
#' @export
iobr_ai_configure_provider <- function(provider) {
  if (!is.list(provider)) {
    stop("Provider must be a list")
  }
  
  if (is.null(provider$name)) {
    stop("Provider name is required")
  }
  
  provider_name <- tolower(provider$name)
  supported_providers <- c("openai", "anthropic", "huggingface", "custom")
  
  if (!provider_name %in% supported_providers) {
    stop("Unsupported provider: ", provider_name, 
         ". Supported providers: ", paste(supported_providers, collapse = ", "))
  }
  
  if (is.null(provider$api_key) || provider$api_key == "") {
    stop("API key is required. Set provider$api_key or use Sys.getenv('PROVIDER_API_KEY')")
  }
  
  # Set defaults based on provider
  validated <- provider
  validated$name <- provider_name
  
  if (provider_name == "openai") {
    validated$base_url <- provider$base_url %||% "https://api.openai.com/v1"
    validated$model_embeddings <- provider$model_embeddings %||% "text-embedding-ada-002"
    validated$model_chat <- provider$model_chat %||% "gpt-3.5-turbo"
  } else if (provider_name == "anthropic") {
    validated$base_url <- provider$base_url %||% "https://api.anthropic.com/v1"
    validated$model_embeddings <- provider$model_embeddings %||% "text-embedding-ada-002"  # Note: Anthropic doesn't have embeddings
    validated$model_chat <- provider$model_chat %||% "claude-3-sonnet-20240229"
  } else if (provider_name == "huggingface") {
    validated$base_url <- provider$base_url %||% "https://api-inference.huggingface.co"
    validated$model_embeddings <- provider$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
    validated$model_chat <- provider$model_chat %||% "mistralai/Mistral-7B-Instruct-v0.1"
  } else if (provider_name == "custom") {
    if (is.null(provider$base_url)) {
      stop("base_url is required for custom provider")
    }
  }
  
  return(validated)
}

#' List index contents
#'
#' Shows information about the current index including number of entries
#' and sources.
#'
#' @param index_path Path to index RDS file. If NULL, uses default location.
#'
#' @return Data frame with index statistics
#'
#' @examples
#' \dontrun{
#' index_info <- iobr_ai_list_index()
#' print(index_info)
#' }
#'
#' @export
iobr_ai_list_index <- function(index_path = NULL) {
  # Set default index path
  if (is.null(index_path)) {
    pkg_inst_path <- system.file("ai", "iobr_embeddings.rds", package = "IOBR")
    if (file.exists(pkg_inst_path)) {
      index_path <- pkg_inst_path
    } else {
      index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
    }
  }
  
  if (!file.exists(index_path)) {
    stop("Index file not found: ", index_path)
  }
  
  index <- load_index(index_path)
  
  # Aggregate statistics
  sources <- sapply(index, function(x) x$source)
  types <- sapply(index, function(x) x$type)
  
  stats <- data.frame(
    total_entries = length(index),
    unique_sources = length(unique(sources)),
    index_path = index_path,
    stringsAsFactors = FALSE
  )
  
  # Add type breakdown
  type_counts <- table(types)
  for (type_name in names(type_counts)) {
    stats[[paste0("count_", type_name)]] <- as.integer(type_counts[type_name])
  }
  
  return(stats)
}

#' Reset/delete AI index
#'
#' Removes the index file from disk.
#'
#' @param index_path Path to index RDS file. If NULL, uses default location.
#'
#' @return TRUE if successful
#'
#' @examples
#' \dontrun{
#' iobr_ai_reset_index()
#' }
#'
#' @export
iobr_ai_reset_index <- function(index_path = NULL) {
  # Set default index path
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (file.exists(index_path)) {
    file.remove(index_path)
    message("Index removed: ", index_path)
  } else {
    message("Index file not found: ", index_path)
  }
  
  invisible(TRUE)
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
