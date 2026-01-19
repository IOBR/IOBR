#' Send embeddings request to configured provider
#'
#' Generate embeddings for input texts using the specified provider.
#' Supports 'dummy' (local bag-of-words), 'openai', and 'huggingface' providers.
#'
#' @param texts Character vector of texts to embed
#' @param provider List containing provider configuration with fields:
#'   name, api_key, base_url, model_embeddings, headers
#' @param batch_size Integer, number of texts to process in each batch (default: 16)
#'
#' @return List of numeric vectors (embeddings) for each input text
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' embeddings <- send_embeddings(c("text1", "text2"), provider)
#' }
send_embeddings <- function(texts, provider, batch_size = 16) {
  if (is.null(provider) || is.null(provider$name)) {
    stop("Provider must be specified with at least a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  # Dispatch to appropriate provider
  if (provider_name == "dummy") {
    return(.send_embeddings_dummy(texts))
  } else if (provider_name == "openai") {
    return(.send_embeddings_openai(texts, provider, batch_size))
  } else if (provider_name == "huggingface") {
    return(.send_embeddings_huggingface(texts, provider, batch_size))
  } else {
    stop(sprintf("Unknown provider: %s. Supported: dummy, openai, huggingface", provider$name))
  }
}

#' Send chat completion request to configured provider
#'
#' Generate a chat response using the specified provider.
#' Supports 'dummy' (simple synthesizer), 'openai', and 'huggingface' providers.
#'
#' @param system_msg Character, system message to set context
#' @param user_msg Character, user message/query
#' @param provider List containing provider configuration
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature (default: 0.2)
#'
#' @return List with fields: content (character), raw (response object)
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' response <- send_chat("You are helpful", "What is R?", provider)
#' }
send_chat <- function(system_msg, user_msg, provider, max_tokens = 800, temperature = 0.2) {
  if (is.null(provider) || is.null(provider$name)) {
    stop("Provider must be specified with at least a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  # Dispatch to appropriate provider
  if (provider_name == "dummy") {
    return(.send_chat_dummy(system_msg, user_msg, max_tokens, temperature))
  } else if (provider_name == "openai") {
    return(.send_chat_openai(system_msg, user_msg, provider, max_tokens, temperature))
  } else if (provider_name == "huggingface") {
    return(.send_chat_huggingface(system_msg, user_msg, provider, max_tokens, temperature))
  } else {
    stop(sprintf("Unknown provider: %s. Supported: dummy, openai, huggingface", provider$name))
  }
}

# Internal: Dummy embeddings using bag-of-words
.send_embeddings_dummy <- function(texts) {
  # Create a simple bag-of-words embedding
  # Build vocabulary from all texts
  all_words <- unique(unlist(strsplit(tolower(paste(texts, collapse = " ")), "\\W+")))
  all_words <- all_words[nchar(all_words) > 0]
  
  # Create embeddings for each text
  embeddings <- lapply(texts, function(text) {
    words <- unlist(strsplit(tolower(text), "\\W+"))
    words <- words[nchar(words) > 0]
    
    # Count occurrences of each vocabulary word
    vec <- sapply(all_words, function(w) sum(words == w))
    
    # Normalize to unit vector
    norm <- sqrt(sum(vec^2))
    if (norm > 0) {
      vec <- vec / norm
    }
    
    as.numeric(vec)
  })
  
  return(embeddings)
}

# Internal: OpenAI embeddings via REST API
.send_embeddings_openai <- function(texts, provider, batch_size) {
  if (is.null(provider$api_key)) {
    stop("OpenAI provider requires api_key")
  }
  if (is.null(provider$model_embeddings)) {
    provider$model_embeddings <- "text-embedding-ada-002"
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  endpoint <- paste0(base_url, "/embeddings")
  
  # Process in batches
  embeddings <- list()
  for (i in seq(1, length(texts), by = batch_size)) {
    batch_end <- min(i + batch_size - 1, length(texts))
    batch_texts <- texts[i:batch_end]
    
    body <- list(
      input = batch_texts,
      model = provider$model_embeddings
    )
    
    # Build headers
    headers <- c(
      "Content-Type" = "application/json",
      "Authorization" = paste("Bearer", provider$api_key)
    )
    if (!is.null(provider$headers)) {
      headers <- c(headers, provider$headers)
    }
    
    # Make request with retry
    response <- .http_request_with_retry(
      url = endpoint,
      method = "POST",
      headers = headers,
      body = body
    )
    
    # Parse response
    if (!is.null(response$error)) {
      stop(sprintf("OpenAI embeddings error: %s", response$error))
    }
    
    batch_embeddings <- lapply(response$data, function(item) item$embedding)
    embeddings <- c(embeddings, batch_embeddings)
  }
  
  return(embeddings)
}

# Internal: Hugging Face embeddings via REST API
.send_embeddings_huggingface <- function(texts, provider, batch_size) {
  if (is.null(provider$api_key)) {
    stop("Hugging Face provider requires api_key")
  }
  if (is.null(provider$model_embeddings)) {
    provider$model_embeddings <- "sentence-transformers/all-MiniLM-L6-v2"
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co/pipeline/feature-extraction"
  endpoint <- sprintf("%s/%s", base_url, provider$model_embeddings)
  
  # Process in batches
  embeddings <- list()
  for (i in seq(1, length(texts), by = batch_size)) {
    batch_end <- min(i + batch_size - 1, length(texts))
    batch_texts <- texts[i:batch_end]
    
    body <- list(inputs = batch_texts)
    
    # Build headers
    headers <- c(
      "Content-Type" = "application/json",
      "Authorization" = paste("Bearer", provider$api_key)
    )
    if (!is.null(provider$headers)) {
      headers <- c(headers, provider$headers)
    }
    
    # Make request with retry
    response <- .http_request_with_retry(
      url = endpoint,
      method = "POST",
      headers = headers,
      body = body
    )
    
    # Parse response - HF returns list of embeddings
    if (!is.null(response$error)) {
      stop(sprintf("Hugging Face embeddings error: %s", response$error))
    }
    
    # Response format varies, handle both array and nested formats
    if (is.list(response)) {
      batch_embeddings <- response
    } else {
      stop("Unexpected Hugging Face response format")
    }
    
    embeddings <- c(embeddings, batch_embeddings)
  }
  
  return(embeddings)
}

# Internal: Dummy chat synthesizer
.send_chat_dummy <- function(system_msg, user_msg, max_tokens, temperature) {
  # Simple template-based response for offline testing
  keywords <- c("IOBR", "immune", "TME", "tumor", "microenvironment", 
                "deconvolution", "signature", "analysis")
  
  found_keywords <- keywords[sapply(keywords, function(k) {
    grepl(k, user_msg, ignore.case = TRUE)
  })]
  
  response_text <- if (length(found_keywords) > 0) {
    sprintf(
      "Based on the IOBR package documentation, here's information about %s:\n\n%s\n\nThe IOBR package provides comprehensive tools for immune-oncology biological research, including tumor microenvironment analysis and signature calculations.",
      paste(found_keywords, collapse = ", "),
      "IOBR offers multiple deconvolution methods (CIBERSORT, TIMER, xCell, etc.) and integrates 322 published signature gene sets for TME analysis."
    )
  } else {
    "I'm a dummy AI assistant for the IOBR package. I can help answer questions about immune-oncology analysis, tumor microenvironment deconvolution, and signature scoring methods. Please note this is a test mode without actual AI capabilities."
  }
  
  list(
    content = response_text,
    raw = list(
      model = "dummy",
      system = system_msg,
      user = user_msg,
      temperature = temperature,
      max_tokens = max_tokens
    )
  )
}

# Internal: OpenAI chat completion
.send_chat_openai <- function(system_msg, user_msg, provider, max_tokens, temperature) {
  if (is.null(provider$api_key)) {
    stop("OpenAI provider requires api_key")
  }
  if (is.null(provider$model_chat)) {
    provider$model_chat <- "gpt-3.5-turbo"
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  endpoint <- paste0(base_url, "/chat/completions")
  
  body <- list(
    model = provider$model_chat,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    max_tokens = max_tokens,
    temperature = temperature
  )
  
  # Build headers
  headers <- c(
    "Content-Type" = "application/json",
    "Authorization" = paste("Bearer", provider$api_key)
  )
  if (!is.null(provider$headers)) {
    headers <- c(headers, provider$headers)
  }
  
  # Make request with retry
  response <- .http_request_with_retry(
    url = endpoint,
    method = "POST",
    headers = headers,
    body = body
  )
  
  if (!is.null(response$error)) {
    stop(sprintf("OpenAI chat error: %s", response$error$message %||% response$error))
  }
  
  content <- response$choices[[1]]$message$content
  
  list(
    content = content,
    raw = response
  )
}

# Internal: Hugging Face chat completion
.send_chat_huggingface <- function(system_msg, user_msg, provider, max_tokens, temperature) {
  if (is.null(provider$api_key)) {
    stop("Hugging Face provider requires api_key")
  }
  if (is.null(provider$model_chat)) {
    provider$model_chat <- "mistralai/Mistral-7B-Instruct-v0.1"
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co/models"
  endpoint <- sprintf("%s/%s", base_url, provider$model_chat)
  
  # Combine system and user messages
  prompt <- sprintf("%s\n\nUser: %s\nAssistant:", system_msg, user_msg)
  
  body <- list(
    inputs = prompt,
    parameters = list(
      max_new_tokens = max_tokens,
      temperature = temperature,
      return_full_text = FALSE
    )
  )
  
  # Build headers
  headers <- c(
    "Content-Type" = "application/json",
    "Authorization" = paste("Bearer", provider$api_key)
  )
  if (!is.null(provider$headers)) {
    headers <- c(headers, provider$headers)
  }
  
  # Make request with retry
  response <- .http_request_with_retry(
    url = endpoint,
    method = "POST",
    headers = headers,
    body = body
  )
  
  if (!is.null(response$error)) {
    stop(sprintf("Hugging Face chat error: %s", response$error))
  }
  
  # Extract generated text
  content <- if (is.list(response) && length(response) > 0) {
    response[[1]]$generated_text %||% response[[1]]$text %||% as.character(response)
  } else {
    as.character(response)
  }
  
  list(
    content = content,
    raw = response
  )
}

# Internal: HTTP request with exponential backoff retry
.http_request_with_retry <- function(url, method = "POST", headers = NULL, 
                                      body = NULL, max_retries = 3) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for HTTP requests")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON handling")
  }
  
  attempt <- 1
  last_error <- NULL
  
  while (attempt <= max_retries) {
    tryCatch({
      # Build httr request
      if (method == "POST") {
        response <- httr::POST(
          url = url,
          httr::add_headers(.headers = headers),
          body = jsonlite::toJSON(body, auto_unbox = TRUE),
          encode = "raw"
        )
      } else {
        stop("Only POST method is currently supported")
      }
      
      # Check status
      if (httr::status_code(response) >= 400) {
        error_content <- httr::content(response, as = "text", encoding = "UTF-8")
        error_parsed <- tryCatch(
          jsonlite::fromJSON(error_content),
          error = function(e) list(message = error_content)
        )
        
        # Retry on rate limit or server errors
        status <- httr::status_code(response)
        if (status == 429 || status >= 500) {
          if (attempt < max_retries) {
            wait_time <- 2^attempt  # Exponential backoff
            message(sprintf("Attempt %d failed with status %d, retrying in %d seconds...", 
                          attempt, status, wait_time))
            Sys.sleep(wait_time)
            attempt <- attempt + 1
            next
          }
        }
        
        return(list(error = error_parsed))
      }
      
      # Parse successful response
      content_text <- httr::content(response, as = "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(content_text))
      
    }, error = function(e) {
      last_error <- e$message
      if (attempt < max_retries) {
        wait_time <- 2^attempt
        message(sprintf("Request error on attempt %d: %s. Retrying in %d seconds...", 
                      attempt, e$message, wait_time))
        Sys.sleep(wait_time)
        attempt <<- attempt + 1
      } else {
        stop(sprintf("HTTP request failed after %d attempts. Last error: %s", 
                    max_retries, e$message))
      }
    })
  }
  
  stop(sprintf("HTTP request failed after %d attempts. Last error: %s", 
              max_retries, last_error))
}

# Utility: NULL-coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
