#' Send Embeddings Request
#'
#' Send texts to an embedding provider and retrieve vector embeddings.
#' Supports dummy (local bag-of-words), OpenAI, and Hugging Face providers.
#'
#' @param texts Character vector of texts to embed
#' @param provider List containing provider configuration with fields:
#'   name, api_key, base_url, model_embeddings, headers
#' @param batch_size Integer, number of texts to process per batch (default: 16)
#'
#' @return List of numeric vectors (embeddings)
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' embeddings <- send_embeddings(c("text1", "text2"), provider)
#' }
send_embeddings <- function(texts, provider, batch_size = 16) {
  if (!is.list(provider) || is.null(provider$name)) {
    stop("provider must be a list with at least a 'name' field")
  }
  
  if (provider$name == "dummy") {
    return(.send_embeddings_dummy(texts))
  } else if (provider$name == "openai") {
    return(.send_embeddings_openai(texts, provider, batch_size))
  } else if (provider$name == "huggingface") {
    return(.send_embeddings_huggingface(texts, provider, batch_size))
  } else {
    stop(sprintf("Unsupported provider: %s", provider$name))
  }
}

#' Send Chat Request
#'
#' Send a chat request to a language model provider.
#' Supports dummy (local synthesizer), OpenAI, and Hugging Face providers.
#'
#' @param system_msg Character, system message/prompt
#' @param user_msg Character, user message/query
#' @param provider List containing provider configuration
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature (default: 0.2)
#'
#' @return List with fields: content (character), raw (list of full response)
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' response <- send_chat("You are a helpful assistant", "Hello", provider)
#' }
send_chat <- function(system_msg, user_msg, provider, 
                      max_tokens = 800, temperature = 0.2) {
  if (!is.list(provider) || is.null(provider$name)) {
    stop("provider must be a list with at least a 'name' field")
  }
  
  if (provider$name == "dummy") {
    return(.send_chat_dummy(system_msg, user_msg, max_tokens))
  } else if (provider$name == "openai") {
    return(.send_chat_openai(system_msg, user_msg, provider, max_tokens, temperature))
  } else if (provider$name == "huggingface") {
    return(.send_chat_huggingface(system_msg, user_msg, provider, max_tokens, temperature))
  } else {
    stop(sprintf("Unsupported provider: %s", provider$name))
  }
}

# Internal: Dummy embedding provider (bag-of-words)
.send_embeddings_dummy <- function(texts) {
  # Simple bag-of-words embedding for offline testing
  vocab <- unique(unlist(strsplit(tolower(paste(texts, collapse = " ")), "[^a-z0-9]+")))
  vocab <- vocab[nzchar(vocab)]
  
  embeddings <- lapply(texts, function(text) {
    words <- strsplit(tolower(text), "[^a-z0-9]+")[[1]]
    words <- words[nzchar(words)]
    vec <- sapply(vocab, function(v) sum(words == v))
    # Normalize
    norm <- sqrt(sum(vec^2))
    if (norm > 0) vec <- vec / norm
    as.numeric(vec)
  })
  
  return(embeddings)
}

# Internal: OpenAI embedding provider
.send_embeddings_openai <- function(texts, provider, batch_size) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for OpenAI provider")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for OpenAI provider")
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  model <- provider$model_embeddings %||% "text-embedding-ada-002"
  api_key <- provider$api_key
  
  if (is.null(api_key) || api_key == "") {
    stop("OpenAI provider requires api_key")
  }
  
  embeddings <- list()
  batches <- split(texts, ceiling(seq_along(texts) / batch_size))
  
  for (batch in batches) {
    url <- paste0(base_url, "/embeddings")
    
    headers <- httr::add_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type" = "application/json"
    )
    
    # Add custom headers if provided
    if (!is.null(provider$headers) && length(provider$headers) > 0) {
      for (name in names(provider$headers)) {
        headers <- c(headers, httr::add_headers(.headers = setNames(provider$headers[[name]], name)))
      }
    }
    
    body <- list(
      input = batch,
      model = model
    )
    
    response <- .http_request_with_retry(
      function() {
        httr::POST(url, headers, body = jsonlite::toJSON(body, auto_unbox = TRUE), 
                   encode = "raw")
      }
    )
    
    if (httr::http_error(response)) {
      stop(sprintf("OpenAI API error: %s", httr::content(response, "text", encoding = "UTF-8")))
    }
    
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    batch_embeddings <- lapply(result$data, function(x) x$embedding)
    embeddings <- c(embeddings, batch_embeddings)
  }
  
  return(embeddings)
}

# Internal: Hugging Face embedding provider
.send_embeddings_huggingface <- function(texts, provider, batch_size) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for Hugging Face provider")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for Hugging Face provider")
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co/pipeline/feature-extraction"
  model <- provider$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
  api_key <- provider$api_key
  
  if (is.null(api_key) || api_key == "") {
    stop("Hugging Face provider requires api_key")
  }
  
  # For HF, the model is in the URL path
  url <- paste0(base_url, "/", model)
  
  embeddings <- list()
  batches <- split(texts, ceiling(seq_along(texts) / batch_size))
  
  for (batch in batches) {
    headers <- httr::add_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type" = "application/json"
    )
    
    # Add custom headers if provided
    if (!is.null(provider$headers) && length(provider$headers) > 0) {
      for (name in names(provider$headers)) {
        headers <- c(headers, httr::add_headers(.headers = setNames(provider$headers[[name]], name)))
      }
    }
    
    body <- list(inputs = batch)
    
    response <- .http_request_with_retry(
      function() {
        httr::POST(url, headers, body = jsonlite::toJSON(body, auto_unbox = TRUE),
                   encode = "raw")
      }
    )
    
    if (httr::http_error(response)) {
      stop(sprintf("Hugging Face API error: %s", httr::content(response, "text", encoding = "UTF-8")))
    }
    
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    # HF returns embeddings directly as matrix/list
    if (is.matrix(result)) {
      batch_embeddings <- lapply(1:nrow(result), function(i) as.numeric(result[i, ]))
    } else {
      batch_embeddings <- lapply(result, as.numeric)
    }
    embeddings <- c(embeddings, batch_embeddings)
  }
  
  return(embeddings)
}

# Internal: Dummy chat provider
.send_chat_dummy <- function(system_msg, user_msg, max_tokens) {
  # Simple dummy response for offline testing
  content <- sprintf(
    "This is a dummy response to your query: '%s'.\n\nBased on the provided context, here's what I can tell you:\n\n%s\n\nThis is a local dummy provider for testing purposes. Configure a real provider (openai, huggingface) to get actual AI responses.",
    substr(user_msg, 1, 100),
    paste(rep("Lorem ipsum dolor sit amet.", 5), collapse = " ")
  )
  
  # Limit to approximate max_tokens (rough estimate: 1 token ~ 4 chars)
  max_chars <- max_tokens * 4
  if (nchar(content) > max_chars) {
    content <- substr(content, 1, max_chars)
  }
  
  list(
    content = content,
    raw = list(
      model = "dummy",
      usage = list(total_tokens = floor(nchar(content) / 4))
    )
  )
}

# Internal: OpenAI chat provider
.send_chat_openai <- function(system_msg, user_msg, provider, max_tokens, temperature) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for OpenAI provider")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for OpenAI provider")
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  model <- provider$model_chat %||% "gpt-3.5-turbo"
  api_key <- provider$api_key
  
  if (is.null(api_key) || api_key == "") {
    stop("OpenAI provider requires api_key")
  }
  
  url <- paste0(base_url, "/chat/completions")
  
  headers <- httr::add_headers(
    "Authorization" = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )
  
  # Add custom headers if provided
  if (!is.null(provider$headers) && length(provider$headers) > 0) {
    for (name in names(provider$headers)) {
      headers <- c(headers, httr::add_headers(.headers = setNames(provider$headers[[name]], name)))
    }
  }
  
  body <- list(
    model = model,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    max_tokens = max_tokens,
    temperature = temperature
  )
  
  response <- .http_request_with_retry(
    function() {
      httr::POST(url, headers, body = jsonlite::toJSON(body, auto_unbox = TRUE),
                 encode = "raw")
    }
  )
  
  if (httr::http_error(response)) {
    stop(sprintf("OpenAI API error: %s", httr::content(response, "text", encoding = "UTF-8")))
  }
  
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  list(
    content = result$choices[[1]]$message$content,
    raw = result
  )
}

# Internal: Hugging Face chat provider
.send_chat_huggingface <- function(system_msg, user_msg, provider, max_tokens, temperature) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for Hugging Face provider")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for Hugging Face provider")
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co/models"
  model <- provider$model_chat %||% "mistralai/Mistral-7B-Instruct-v0.1"
  api_key <- provider$api_key
  
  if (is.null(api_key) || api_key == "") {
    stop("Hugging Face provider requires api_key")
  }
  
  url <- paste0(base_url, "/", model)
  
  headers <- httr::add_headers(
    "Authorization" = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )
  
  # Add custom headers if provided
  if (!is.null(provider$headers) && length(provider$headers) > 0) {
    for (name in names(provider$headers)) {
      headers <- c(headers, httr::add_headers(.headers = setNames(provider$headers[[name]], name)))
    }
  }
  
  # Combine system and user messages
  full_prompt <- sprintf("%s\n\nUser: %s\nAssistant:", system_msg, user_msg)
  
  body <- list(
    inputs = full_prompt,
    parameters = list(
      max_new_tokens = max_tokens,
      temperature = temperature,
      return_full_text = FALSE
    )
  )
  
  response <- .http_request_with_retry(
    function() {
      httr::POST(url, headers, body = jsonlite::toJSON(body, auto_unbox = TRUE),
                 encode = "raw")
    }
  )
  
  if (httr::http_error(response)) {
    stop(sprintf("Hugging Face API error: %s", httr::content(response, "text", encoding = "UTF-8")))
  }
  
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  # HF returns array of results
  content <- if (is.list(result) && length(result) > 0) {
    result[[1]]$generated_text
  } else {
    result$generated_text
  }
  
  list(
    content = content,
    raw = result
  )
}

# Internal: HTTP request with exponential backoff retry
.http_request_with_retry <- function(request_fn, max_retries = 3, initial_delay = 1) {
  for (attempt in 1:max_retries) {
    tryCatch({
      response <- request_fn()
      return(response)
    }, error = function(e) {
      if (attempt == max_retries) {
        stop(sprintf("HTTP request failed after %d attempts: %s", max_retries, e$message))
      }
      
      delay <- initial_delay * (2 ^ (attempt - 1))
      message(sprintf("Request failed (attempt %d/%d), retrying in %d seconds...", 
                      attempt, max_retries, delay))
      Sys.sleep(delay)
    })
  }
}

# Null-coalescing operator helper
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
