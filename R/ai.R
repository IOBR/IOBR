#' Initialize IOBR AI assistant and create index
#'
#' Collect package documentation, chunk it, and create an embeddings index.
#'
#' @param pkg_path Character, path to package root (default: current directory)
#' @param provider List, provider configuration (default: list(name = "dummy"))
#' @param index_path Character, path to save index (default: inst/ai/iobr_embeddings.rds)
#' @param chunk_size Integer, chunk size in characters (default: 800)
#'
#' @return List with index_path and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize with dummy provider (no API key needed)
#' result <- iobr_ai_init()
#' 
#' # Initialize with OpenAI
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = "sk-..."
#' ))
#' result <- iobr_ai_init(provider = provider)
#' }
iobr_ai_init <- function(pkg_path = ".", provider = NULL, index_path = NULL, chunk_size = 800) {
  # Default provider to dummy if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
    message("Using dummy provider (local mode, no API key required)")
  }
  
  # Validate and normalize provider
  provider <- iobr_ai_configure_provider(provider)
  
  # Set default index path
  if (is.null(index_path)) {
    index_path <- file.path(pkg_path, "inst", "ai", "iobr_embeddings.rds")
  }
  
  # Collect documents
  message("Collecting package documentation...")
  docs <- collect_package_docs(pkg_path)
  
  if (length(docs) == 0) {
    stop("No documentation found. Please check pkg_path.")
  }
  
  # Chunk documents
  message("Chunking documents...")
  chunks <- chunk_texts(docs, chunk_size = chunk_size)
  
  if (length(chunks) == 0) {
    stop("No chunks created from documents.")
  }
  
  # Create index
  message("Creating embeddings index...")
  create_index(chunks, provider, index_path = index_path)
  
  # Return metadata
  list(
    index_path = index_path,
    num_chunks = length(chunks),
    num_docs = length(docs),
    provider = provider$name,
    created_at = Sys.time()
  )
}

#' Query the IOBR AI assistant
#'
#' Search the index and generate an AI response with retrieved context.
#'
#' @param query Character, user query
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param provider List, provider configuration (default: list(name = "dummy"))
#' @param top_k Integer, number of context chunks to retrieve (default: 5)
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature (default: 0.2)
#'
#' @return List with answer, code (if any), retrieved chunks, and raw response
#' @export
#'
#' @examples
#' \dontrun{
#' # Query with dummy provider
#' result <- iobr_ai_query("How do I use CIBERSORT?")
#' cat(result$answer)
#' 
#' # View retrieved sources
#' print(result$retrieved)
#' }
iobr_ai_query <- function(query, index_path = NULL, provider = NULL, 
                          top_k = 5, max_tokens = 800, temperature = 0.2) {
  # Default provider to dummy if not specified
  if (is.null(provider)) {
    provider <- list(name = "dummy")
  }
  
  # Validate and normalize provider
  provider <- iobr_ai_configure_provider(provider)
  
  # Set default index path
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  # Load index
  if (!file.exists(index_path)) {
    stop(sprintf("Index not found at %s. Please run iobr_ai_init() first.", index_path))
  }
  
  message("Loading index...")
  index <- readRDS(index_path)
  
  # Search index
  message("Searching for relevant context...")
  retrieved <- search_index(query, index, provider, top_k = top_k)
  
  # Build prompt with context
  context_text <- paste(
    sapply(seq_len(nrow(retrieved)), function(i) {
      sprintf("[Source %d: %s (score: %.3f)]\n%s", 
              i, retrieved$source[i], retrieved$score[i], retrieved$text[i])
    }),
    collapse = "\n\n"
  )
  
  system_msg <- "You are a helpful AI assistant for the IOBR R package, specialized in immune-oncology biological research, tumor microenvironment analysis, and bioinformatics. Answer questions based on the provided context from the package documentation. If you provide R code examples, wrap them in ```r ``` code blocks. Be concise and accurate."
  
  user_msg <- sprintf(
    "Context from IOBR package documentation:\n\n%s\n\nUser question: %s\n\nPlease answer based on the context above. If relevant, provide R code examples.",
    context_text,
    query
  )
  
  # Get AI response
  message("Generating response...")
  response <- send_chat(system_msg, user_msg, provider, max_tokens, temperature)
  
  # Extract code blocks if present
  code_blocks <- .extract_code_blocks(response$content)
  
  # Build result
  list(
    answer = response$content,
    code = if (length(code_blocks) > 0) code_blocks else NULL,
    retrieved = retrieved,
    raw = response$raw
  )
}

#' Configure and validate provider settings
#'
#' Normalize provider configuration and apply defaults.
#'
#' @param provider List with provider settings (name, api_key, base_url, models, headers)
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
#'   api_key = "sk-..."
#' ))
#' 
#' # Custom configuration
#' provider <- iobr_ai_configure_provider(list(
#'   name = "openai",
#'   api_key = "sk-...",
#'   base_url = "https://api.openai.com/v1",
#'   model_embeddings = "text-embedding-ada-002",
#'   model_chat = "gpt-4",
#'   headers = list("Organization" = "org-...")
#' ))
#' }
iobr_ai_configure_provider <- function(provider) {
  if (is.null(provider) || !is.list(provider)) {
    stop("Provider must be a list")
  }
  
  if (is.null(provider$name)) {
    stop("Provider must have a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  # Apply defaults based on provider
  if (provider_name == "dummy") {
    # No additional configuration needed
    provider$model_embeddings <- "bag-of-words"
    provider$model_chat <- "template-based"
    
  } else if (provider_name == "openai") {
    # OpenAI defaults
    if (is.null(provider$base_url)) {
      provider$base_url <- "https://api.openai.com/v1"
    }
    if (is.null(provider$model_embeddings)) {
      provider$model_embeddings <- "text-embedding-ada-002"
    }
    if (is.null(provider$model_chat)) {
      provider$model_chat <- "gpt-3.5-turbo"
    }
    if (is.null(provider$api_key) || nchar(provider$api_key) == 0) {
      stop("OpenAI provider requires api_key")
    }
    
  } else if (provider_name == "huggingface") {
    # Hugging Face defaults
    if (is.null(provider$base_url)) {
      provider$base_url <- "https://api-inference.huggingface.co"
    }
    if (is.null(provider$model_embeddings)) {
      provider$model_embeddings <- "sentence-transformers/all-MiniLM-L6-v2"
    }
    if (is.null(provider$model_chat)) {
      provider$model_chat <- "mistralai/Mistral-7B-Instruct-v0.1"
    }
    if (is.null(provider$api_key) || nchar(provider$api_key) == 0) {
      stop("Hugging Face provider requires api_key")
    }
    
  } else if (provider_name == "anthropic") {
    # Anthropic/Claude defaults
    if (is.null(provider$base_url)) {
      provider$base_url <- "https://api.anthropic.com/v1"
    }
    if (is.null(provider$model_chat)) {
      provider$model_chat <- "claude-3-sonnet-20240229"
    }
    if (is.null(provider$api_key) || nchar(provider$api_key) == 0) {
      stop("Anthropic provider requires api_key")
    }
    message("Note: Anthropic provider is configured but may require additional adapter implementation")
    
  } else {
    message(sprintf("Using custom provider: %s (ensure adapters are available)", provider$name))
  }
  
  return(provider)
}

#' List index information
#'
#' Display metadata about the current index.
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
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("No index found at %s", index_path))
    return(invisible(NULL))
  }
  
  index <- readRDS(index_path)
  
  info <- list(
    path = index_path,
    num_entries = length(index$entries),
    created_at = index$created_at,
    provider_name = index$provider_meta$name,
    model_embeddings = index$provider_meta$model_embeddings,
    file_size_mb = file.info(index_path)$size / (1024^2)
  )
  
  message(sprintf("Index: %s", info$path))
  message(sprintf("  Entries: %d", info$num_entries))
  message(sprintf("  Created: %s", info$created_at))
  message(sprintf("  Provider: %s (%s)", info$provider_name, info$model_embeddings))
  message(sprintf("  Size: %.2f MB", info$file_size_mb))
  
  invisible(info)
}

#' Reset (delete) the index
#'
#' Remove the index file from disk.
#'
#' @param index_path Character, path to index RDS file (default: inst/ai/iobr_embeddings.rds)
#'
#' @return Logical, TRUE if deleted successfully
#' @export
#'
#' @examples
#' \dontrun{
#' iobr_ai_reset_index()
#' }
iobr_ai_reset_index <- function(index_path = NULL) {
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  if (!file.exists(index_path)) {
    message(sprintf("No index found at %s", index_path))
    return(invisible(FALSE))
  }
  
  unlink(index_path)
  message(sprintf("Index deleted: %s", index_path))
  invisible(TRUE)
}

# Internal: Extract code blocks from markdown text
.extract_code_blocks <- function(text) {
  # Match code blocks with ```r or ```R
  pattern <- "```[rR]\\s*\n([^`]+)```"
  matches <- gregexpr(pattern, text, perl = TRUE)
  
  if (matches[[1]][1] == -1) {
    return(character(0))
  }
  
  # Extract matched code
  match_data <- regmatches(text, matches)[[1]]
  
  # Remove ``` delimiters
  code_blocks <- gsub("```[rR]\\s*\n", "", match_data)
  code_blocks <- gsub("```$", "", code_blocks)
  code_blocks <- trimws(code_blocks)
  
  return(code_blocks)
}
