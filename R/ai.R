#' Initialize IOBR AI Assistant
#'
#' Build the embedding index from package documentation for the AI assistant.
#'
#' @param pkg_path Character, path to package root (default: current directory)
#' @param provider List with provider configuration (name, api_key, base_url, etc.).
#'   If NULL, defaults to dummy provider.
#' @param index_path Character, path to save index (default: inst/ai/iobr_embeddings.rds)
#' @param chunk_size Integer, approximate character size per chunk (default: 800)
#'
#' @return List with index_path and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize with dummy provider (no API key needed)
#' metadata <- iobr_ai_init()
#' 
#' # Initialize with OpenAI
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' ))
#' metadata <- iobr_ai_init(provider = provider)
#' }
iobr_ai_init <- function(pkg_path = ".", provider = NULL, 
                         index_path = NULL, chunk_size = 800) {
  # Default to dummy provider if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
    message("Using dummy provider for offline testing")
  }
  
  # Validate provider
  provider <- iobr_ai_configure_provider(provider)
  
  # Collect documentation
  message("Collecting package documentation...")
  docs <- collect_package_docs(pkg_path)
  
  # Chunk texts
  message("Chunking texts...")
  chunks <- chunk_texts(docs, chunk_size)
  
  # Create index
  message("Creating embedding index...")
  if (is.null(index_path)) {
    index_path <- file.path(pkg_path, "inst/ai/iobr_embeddings.rds")
  }
  
  create_index(chunks, provider, index_path, batch_size = 16)
  
  metadata <- list(
    index_path = index_path,
    provider = provider$name,
    num_chunks = length(chunks),
    num_docs = length(docs),
    created_at = Sys.time()
  )
  
  message("AI assistant initialized successfully!")
  return(invisible(metadata))
}

#' Query IOBR AI Assistant
#'
#' Ask the AI assistant a question using RAG (Retrieval-Augmented Generation).
#'
#' @param query Character, user question or query
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param provider List with provider configuration. If NULL, defaults to dummy provider.
#' @param top_k Integer, number of context chunks to retrieve (default: 5)
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature (default: 0.2)
#'
#' @return List with fields: answer, code, retrieved, raw
#' @export
#'
#' @examples
#' \dontrun{
#' # Query with dummy provider
#' response <- iobr_ai_query("How do I calculate TME scores?")
#' cat(response$answer)
#' 
#' # Query with OpenAI
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' ))
#' response <- iobr_ai_query("Explain deconvolution methods", provider = provider)
#' }
iobr_ai_query <- function(query, index_path = NULL, provider = NULL,
                          top_k = 5, max_tokens = 800, temperature = 0.2) {
  # Default to dummy provider if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
  }
  
  # Validate provider
  provider <- iobr_ai_configure_provider(provider)
  
  # Default index path
  if (is.null(index_path)) {
    index_path <- "inst/ai/iobr_embeddings.rds"
  }
  
  # Check if index exists
  if (!file.exists(index_path)) {
    stop(sprintf(
      "Index not found at %s. Please run iobr_ai_init() first to create the index.",
      index_path
    ))
  }
  
  # Load index and search
  message("Searching for relevant context...")
  retrieved <- search_index(query, index_path, provider, top_k)
  
  # Build prompt with context
  context_text <- paste(
    sapply(retrieved, function(r) {
      sprintf("Source: %s\n%s", r$source, r$text)
    }),
    collapse = "\n\n---\n\n"
  )
  
  system_msg <- "You are an expert assistant for the IOBR R package, which is used for immune oncology biological research and tumor microenvironment (TME) analysis. Answer questions based on the provided documentation context. If the answer includes R code, format it clearly in a code block. Be concise and accurate. If you cannot answer based on the context, say so."
  
  user_msg <- sprintf(
    "Context from IOBR documentation:\n\n%s\n\n---\n\nUser question: %s\n\nPlease provide a clear answer based on the context above. If relevant, include example R code.",
    context_text,
    query
  )
  
  # Get response from LLM
  message("Generating response...")
  chat_response <- send_chat(system_msg, user_msg, provider, max_tokens, temperature)
  
  # Extract R code if present
  code_blocks <- .extract_code_blocks(chat_response$content)
  
  # Build result
  result <- list(
    answer = chat_response$content,
    code = if (length(code_blocks) > 0) code_blocks else NULL,
    retrieved = retrieved,
    raw = chat_response$raw
  )
  
  return(result)
}

#' Configure AI Provider
#'
#' Validate and normalize provider configuration with sensible defaults.
#'
#' @param provider List with provider configuration. Must include 'name' field.
#'   Optional fields: api_key, base_url, model_embeddings, model_chat, headers
#'
#' @return Validated and normalized provider list
#' @export
#'
#' @examples
#' \dontrun{
#' # Dummy provider
#' provider <- iobr_ai_configure_provider(list(name = "dummy"))
#' 
#' # OpenAI with defaults
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' ))
#' 
#' # Custom provider
#' provider <- iobr_ai_configure_provider(list(
#'   name = "custom",
#'   api_key = "my-key",
#'   base_url = "https://my-api.com",
#'   model_embeddings = "my-embed-model",
#'   model_chat = "my-chat-model"
#' ))
#' }
iobr_ai_configure_provider <- function(provider) {
  if (!is.list(provider)) {
    stop("provider must be a list")
  }
  
  if (is.null(provider$name) || !is.character(provider$name)) {
    stop("provider must have a 'name' field (character)")
  }
  
  provider$name <- tolower(provider$name)
  
  # Apply defaults based on provider name
  if (provider$name == "openai") {
    provider$base_url <- provider$base_url %||% "https://api.openai.com/v1"
    provider$model_embeddings <- provider$model_embeddings %||% "text-embedding-ada-002"
    provider$model_chat <- provider$model_chat %||% "gpt-3.5-turbo"
    
    if (is.null(provider$api_key) || provider$api_key == "") {
      warning("OpenAI provider requires an api_key. Set it or use dummy provider for testing.")
    }
  } else if (provider$name == "huggingface") {
    provider$base_url <- provider$base_url %||% "https://api-inference.huggingface.co"
    provider$model_embeddings <- provider$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
    provider$model_chat <- provider$model_chat %||% "mistralai/Mistral-7B-Instruct-v0.1"
    
    if (is.null(provider$api_key) || provider$api_key == "") {
      warning("Hugging Face provider requires an api_key. Set it or use dummy provider for testing.")
    }
  } else if (provider$name == "anthropic") {
    provider$base_url <- provider$base_url %||% "https://api.anthropic.com/v1"
    provider$model_chat <- provider$model_chat %||% "claude-3-sonnet-20240229"
    
    message("Note: Anthropic provider support is minimal. You may need to customize the implementation.")
    
    if (is.null(provider$api_key) || provider$api_key == "") {
      warning("Anthropic provider requires an api_key.")
    }
  } else if (provider$name == "dummy") {
    # Dummy provider needs no configuration
    message("Dummy provider configured for offline testing")
  } else if (provider$name == "custom") {
    # Custom provider: user must supply all required fields
    if (is.null(provider$base_url)) {
      warning("Custom provider should specify base_url")
    }
  } else {
    warning(sprintf("Unknown provider: %s. Treating as custom.", provider$name))
  }
  
  return(provider)
}

#' List Index Information
#'
#' Display information about the current embedding index.
#'
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#'
#' @return List with index metadata
#' @export
#'
#' @examples
#' \dontrun{
#' info <- iobr_ai_list_index()
#' print(info)
#' }
iobr_ai_list_index <- function(index_path = NULL) {
  if (is.null(index_path)) {
    index_path <- "inst/ai/iobr_embeddings.rds"
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("Index not found at %s", index_path))
    return(invisible(NULL))
  }
  
  index <- readRDS(index_path)
  
  info <- list(
    index_path = index_path,
    num_entries = length(index$entries),
    created_at = index$created_at,
    provider = index$provider_meta$name,
    model = index$provider_meta$model_embeddings,
    file_size = file.info(index_path)$size
  )
  
  message(sprintf("Index: %s", index_path))
  message(sprintf("Entries: %d", info$num_entries))
  message(sprintf("Created: %s", info$created_at))
  message(sprintf("Provider: %s", info$provider))
  message(sprintf("File size: %.2f MB", info$file_size / 1024^2))
  
  return(invisible(info))
}

#' Reset Index
#'
#' Delete the embedding index file.
#'
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#'
#' @return Logical, TRUE if deletion successful
#' @export
#'
#' @examples
#' \dontrun{
#' iobr_ai_reset_index()
#' }
iobr_ai_reset_index <- function(index_path = NULL) {
  if (is.null(index_path)) {
    index_path <- "inst/ai/iobr_embeddings.rds"
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("Index not found at %s (nothing to reset)", index_path))
    return(invisible(FALSE))
  }
  
  success <- file.remove(index_path)
  
  if (success) {
    message(sprintf("Index deleted: %s", index_path))
  } else {
    warning(sprintf("Failed to delete index: %s", index_path))
  }
  
  return(invisible(success))
}

# Internal: Extract code blocks from text
.extract_code_blocks <- function(text) {
  # Extract code blocks marked with ```r or ``` or just code patterns
  patterns <- c(
    "```r\\s*\\n([^`]+)```",
    "```\\s*\\n([^`]+)```"
  )
  
  code_blocks <- list()
  for (pattern in patterns) {
    matches <- gregexpr(pattern, text, perl = TRUE)
    if (matches[[1]][1] != -1) {
      match_starts <- matches[[1]]
      match_lengths <- attr(matches[[1]], "match.length")
      
      for (i in seq_along(match_starts)) {
        code_text <- substr(text, match_starts[i], match_starts[i] + match_lengths[i] - 1)
        # Extract just the code content
        code_content <- gsub("^```r?\\s*\\n|```$", "", code_text, perl = TRUE)
        code_blocks[[length(code_blocks) + 1]] <- trimws(code_content)
      }
    }
  }
  
  return(unique(code_blocks))
}

# Null-coalescing operator helper (in case not defined in provider file)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
