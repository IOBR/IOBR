#' Initialize IOBR AI Assistant
#'
#' @description
#' Initialize the AI assistant by collecting package documentation, chunking it,
#' and creating an embeddings index.
#'
#' @param pkg_path Character, path to package root directory (default: current directory)
#' @param provider List with provider configuration. If NULL, defaults to dummy provider.
#' @param index_path Character, path to save index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param chunk_size Integer, maximum characters per chunk (default: 800)
#'
#' @return List with index metadata including path, number of entries, creation time
#'
#' @examples
#' \dontrun{
#' # Initialize with dummy provider (no API key needed)
#' result <- iobr_ai_init()
#'
#' # Initialize with OpenAI
#' provider <- list(
#'   name = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   model_embeddings = "text-embedding-ada-002"
#' )
#' result <- iobr_ai_init(provider = provider)
#' }
#'
#' @export
iobr_ai_init <- function(pkg_path = ".", provider = NULL, index_path = NULL, chunk_size = 800) {
  message("=== IOBR AI Assistant Initialization ===")
  
  # Default to dummy provider if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
    message("Using dummy provider (offline mode)")
  } else {
    message(sprintf("Using provider: %s", provider$name))
  }
  
  # Validate provider
  provider <- iobr_ai_configure_provider(provider)
  
  # Step 1: Collect documentation
  message("\n1. Collecting package documentation...")
  docs <- collect_package_docs(pkg_path)
  
  # Step 2: Chunk texts
  message("\n2. Chunking texts...")
  chunks <- chunk_texts(docs, chunk_size = chunk_size)
  
  # Step 3: Create index
  message("\n3. Creating embeddings index...")
  index_meta <- create_index(chunks, provider, index_path = index_path)
  
  message("\n=== Initialization Complete ===")
  message(sprintf("Index path: %s", index_meta$path))
  message(sprintf("Total entries: %d", index_meta$num_entries))
  
  return(index_meta)
}


#' Query IOBR AI Assistant
#'
#' @description
#' Query the AI assistant using retrieval-augmented generation (RAG).
#' Retrieves relevant context from the index and generates a response.
#'
#' @param query Character, user query/question
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param provider List with provider configuration. If NULL, defaults to dummy provider.
#' @param top_k Integer, number of context chunks to retrieve (default: 5)
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature 0-1 (default: 0.2)
#'
#' @return List with 'answer', 'code', 'retrieved' (context chunks), and 'raw' (API response)
#'
#' @examples
#' \dontrun{
#' # Query with dummy provider
#' response <- iobr_ai_query("How do I analyze the tumor microenvironment?")
#' cat(response$answer)
#'
#' # Query with OpenAI
#' provider <- list(name = "openai", api_key = Sys.getenv("OPENAI_API_KEY"))
#' response <- iobr_ai_query("What is CIBERSORT?", provider = provider)
#' }
#'
#' @export
iobr_ai_query <- function(query, index_path = NULL, provider = NULL,
                          top_k = 5, max_tokens = 800, temperature = 0.2) {
  if (is.null(query) || nchar(trimws(query)) == 0) {
    stop("Query cannot be empty")
  }
  
  # Default index path
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (!file.exists(index_path)) {
    stop(sprintf("Index not found at: %s. Please run iobr_ai_init() first.", index_path))
  }
  
  # Default to dummy provider if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
  }
  
  # Validate provider
  provider <- iobr_ai_configure_provider(provider)
  
  message("Retrieving relevant context...")
  
  # Step 1: Retrieve relevant chunks
  retrieved <- search_index(query, index_path, provider, top_k = top_k)
  
  message(sprintf("Retrieved %d relevant chunks", length(retrieved)))
  
  # Step 2: Build prompt with context
  context_text <- paste(
    sapply(retrieved, function(r) {
      sprintf("Source: %s\n%s", r$source, r$text)
    }),
    collapse = "\n\n---\n\n"
  )
  
  system_msg <- "You are an AI assistant for the IOBR R package, which provides tools for Immuno-Oncology Biological Research. Your role is to help users understand and use the package based on the provided documentation context. Answer questions clearly and provide R code examples when appropriate. If you provide code, ensure it follows R best practices and uses IOBR functions correctly."
  
  user_msg <- sprintf(
    "Based on the following documentation from the IOBR package, please answer the user's question.\n\nContext:\n%s\n\n---\n\nUser Question: %s\n\nPlease provide a clear answer. If relevant, include R code examples in a code block marked with ```r and ```. Do not execute code, only provide examples.",
    context_text,
    query
  )
  
  message("Generating response...")
  
  # Step 3: Get LLM response
  chat_response <- send_chat(system_msg, user_msg, provider,
                             max_tokens = max_tokens,
                             temperature = temperature)
  
  # Step 4: Extract code if present
  answer_text <- chat_response$content
  code_blocks <- extract_code_blocks(answer_text)
  
  # Build response
  result <- list(
    answer = answer_text,
    code = if (length(code_blocks) > 0) code_blocks else NULL,
    retrieved = retrieved,
    raw = chat_response$raw
  )
  
  return(result)
}


#' Configure and validate provider settings
#'
#' @description
#' Validate and normalize provider configuration, applying sensible defaults.
#'
#' @param provider List with provider configuration
#'
#' @return Validated and normalized provider configuration
#'
#' @examples
#' \dontrun{
#' # Configure OpenAI provider
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = "sk-..."
#' ))
#'
#' # Configure with custom base URL
#' provider <- iobr_ai_configure_provider(list(
#'   name = "custom",
#'   api_key = "...",
#'   base_url = "https://my-api.com/v1",
#'   model_embeddings = "my-embed-model",
#'   model_chat = "my-chat-model"
#' ))
#' }
#'
#' @export
iobr_ai_configure_provider <- function(provider) {
  if (is.null(provider)) {
    stop("Provider configuration cannot be NULL")
  }
  
  if (is.null(provider$name)) {
    stop("Provider must have a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  # Apply defaults based on provider name
  if (provider_name == "openai") {
    provider$base_url <- provider$base_url %||% "https://api.openai.com/v1"
    provider$model_embeddings <- provider$model_embeddings %||% "text-embedding-ada-002"
    provider$model_chat <- provider$model_chat %||% "gpt-3.5-turbo"
    
    if (is.null(provider$api_key) || provider$api_key == "") {
      warning("OpenAI provider requires an API key. Set provider$api_key or use dummy provider for testing.")
    }
  } else if (provider_name == "huggingface") {
    provider$base_url <- provider$base_url %||% "https://api-inference.huggingface.co"
    provider$model_embeddings <- provider$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
    provider$model_chat <- provider$model_chat %||% "microsoft/DialoGPT-medium"
    
    if (is.null(provider$api_key) || provider$api_key == "") {
      warning("Hugging Face provider requires an API key. Set provider$api_key or use dummy provider for testing.")
    }
  } else if (provider_name == "anthropic") {
    # Placeholder for future Anthropic support
    provider$base_url <- provider$base_url %||% "https://api.anthropic.com/v1"
    warning("Anthropic provider is not yet fully implemented. Use OpenAI or Hugging Face instead.")
  } else if (provider_name == "dummy") {
    # Dummy provider needs no configuration
    provider$model_embeddings <- "dummy-embeddings"
    provider$model_chat <- "dummy-chat"
  } else if (provider_name == "custom") {
    # Custom provider requires user to specify all parameters
    if (is.null(provider$base_url)) {
      stop("Custom provider requires 'base_url' to be specified")
    }
    if (is.null(provider$model_embeddings)) {
      warning("Custom provider: 'model_embeddings' not specified")
    }
    if (is.null(provider$model_chat)) {
      warning("Custom provider: 'model_chat' not specified")
    }
  } else {
    warning(sprintf("Unknown provider: %s. Treating as custom provider.", provider$name))
  }
  
  return(provider)
}


#' List index information
#'
#' @description
#' Display information about the current embeddings index.
#'
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#'
#' @return List with index metadata
#'
#' @examples
#' \dontrun{
#' info <- iobr_ai_list_index()
#' print(info)
#' }
#'
#' @export
iobr_ai_list_index <- function(index_path = NULL) {
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("No index found at: %s", index_path))
    return(NULL)
  }
  
  index <- readRDS(index_path)
  
  info <- list(
    path = index_path,
    num_entries = length(index$entries),
    created_at = index$created_at,
    provider = index$provider_meta$name,
    model_embeddings = index$provider_meta$model_embeddings
  )
  
  message("=== Index Information ===")
  message(sprintf("Path: %s", info$path))
  message(sprintf("Entries: %d", info$num_entries))
  message(sprintf("Created: %s", info$created_at))
  message(sprintf("Provider: %s", info$provider))
  if (!is.na(info$model_embeddings)) {
    message(sprintf("Model: %s", info$model_embeddings))
  }
  
  return(invisible(info))
}


#' Reset/delete embeddings index
#'
#' @description
#' Delete the embeddings index file.
#'
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#'
#' @return Logical, TRUE if index was deleted, FALSE if no index existed
#'
#' @examples
#' \dontrun{
#' iobr_ai_reset_index()
#' }
#'
#' @export
iobr_ai_reset_index <- function(index_path = NULL) {
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("No index found at: %s", index_path))
    return(FALSE)
  }
  
  unlink(index_path)
  message(sprintf("Index deleted: %s", index_path))
  return(TRUE)
}


# Internal helper functions -----------------------------------------------

#' Extract code blocks from markdown text
#' @keywords internal
extract_code_blocks <- function(text) {
  # Match code blocks with ```r or ``` markers
  pattern <- "```(?:r)?\\s*\\n([^`]+)```"
  matches <- gregexpr(pattern, text, perl = TRUE)
  
  if (matches[[1]][1] == -1) {
    return(list())
  }
  
  code_blocks <- list()
  match_starts <- matches[[1]]
  match_lengths <- attr(matches[[1]], "match.length")
  
  for (i in seq_along(match_starts)) {
    match_text <- substr(text, match_starts[i], match_starts[i] + match_lengths[i] - 1)
    # Extract content between markers
    code_content <- gsub("```(?:r)?\\s*\\n|```", "", match_text, perl = TRUE)
    code_blocks <- c(code_blocks, list(trimws(code_content)))
  }
  
  return(code_blocks)
}


#' NULL-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
