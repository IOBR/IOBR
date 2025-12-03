#' Provider-agnostic AI client for embeddings and chat
#'
#' Internal functions for sending embedding and chat requests to various AI providers
#'
#' @name ai_provider
#' @keywords internal
NULL

#' Send embedding request to AI provider
#'
#' @param texts Character vector of texts to embed
#' @param provider_config List containing provider configuration:
#'   \itemize{
#'     \item name: Provider name ('openai', 'anthropic', 'huggingface', 'custom')
#'     \item api_key: API key for authentication
#'     \item base_url: Base URL for API (optional, uses defaults for known providers)
#'     \item model_embeddings: Model name for embeddings
#'     \item headers: Named list of additional HTTP headers (optional)
#'   }
#'
#' @return List of numeric vectors (embeddings) or error
#' @keywords internal
send_embedding <- function(texts, provider_config) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for AI features. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for AI features. Please install it.")
  }
  
  # Validate provider config
  if (is.null(provider_config$name)) {
    stop("Provider name is required")
  }
  if (is.null(provider_config$api_key) || provider_config$api_key == "") {
    stop("API key is required")
  }
  
  provider_name <- tolower(provider_config$name)
  
  # Set defaults for known providers
  if (provider_name == "openai") {
    base_url <- provider_config$base_url %||% "https://api.openai.com/v1"
    model <- provider_config$model_embeddings %||% "text-embedding-ada-002"
    endpoint <- paste0(base_url, "/embeddings")
    
    # Build request body
    body <- list(
      input = texts,
      model = model
    )
    
    # Make request with retry logic
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "Authorization" = paste("Bearer", provider_config$api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(60)
      
      # Add custom headers if provided
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      # Extract embeddings
      embeddings <- lapply(content$data, function(x) x$embedding)
      embeddings
    }, error = function(e) {
      stop("Embedding request failed: ", e$message)
    })
    
    return(result)
    
  } else if (provider_name == "huggingface") {
    base_url <- provider_config$base_url %||% "https://api-inference.huggingface.co"
    model <- provider_config$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
    endpoint <- paste0(base_url, "/pipeline/feature-extraction/", model)
    
    # Hugging Face expects different format - send each text separately
    embeddings <- lapply(texts, function(text) {
      tryCatch({
        req <- httr2::request(endpoint) |>
          httr2::req_headers(
            "Authorization" = paste("Bearer", provider_config$api_key),
            "Content-Type" = "application/json"
          ) |>
          httr2::req_body_json(list(inputs = text)) |>
          httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
          httr2::req_timeout(60)
        
        if (!is.null(provider_config$headers)) {
          for (header_name in names(provider_config$headers)) {
            req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
          }
        }
        
        resp <- httr2::req_perform(req)
        content <- httr2::resp_body_json(resp)
        
        # HF returns nested list, flatten to vector
        if (is.list(content) && length(content) > 0) {
          if (is.list(content[[1]])) {
            # Mean pooling if multiple vectors returned
            embedding <- colMeans(do.call(rbind, content))
          } else {
            embedding <- unlist(content)
          }
          return(as.numeric(embedding))
        } else {
          stop("Unexpected response format from Hugging Face")
        }
      }, error = function(e) {
        stop("Embedding request failed for text: ", e$message)
      })
    })
    
    return(embeddings)
    
  } else if (provider_name == "custom") {
    if (is.null(provider_config$base_url)) {
      stop("base_url is required for custom provider")
    }
    endpoint <- provider_config$base_url
    
    # Custom provider - assume OpenAI-compatible format
    body <- list(
      input = texts,
      model = provider_config$model_embeddings %||% "default"
    )
    
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "Authorization" = paste("Bearer", provider_config$api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(60)
      
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      # Try to extract embeddings from response
      if (!is.null(content$data)) {
        embeddings <- lapply(content$data, function(x) x$embedding)
      } else if (!is.null(content$embeddings)) {
        embeddings <- content$embeddings
      } else {
        stop("Could not extract embeddings from response")
      }
      embeddings
    }, error = function(e) {
      stop("Embedding request failed: ", e$message)
    })
    
    return(result)
  } else {
    stop("Unsupported provider: ", provider_name, 
         ". Supported providers: openai, huggingface, custom")
  }
}

#' Send chat/generation request to AI provider
#'
#' @param system System message
#' @param user User message
#' @param provider_config List containing provider configuration
#' @param max_tokens Maximum tokens to generate
#' @param temperature Sampling temperature
#'
#' @return List with 'content' (assistant response text) and 'raw' (full response)
#' @keywords internal
send_chat <- function(system, user, provider_config, max_tokens = 800, temperature = 0.2) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for AI features. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for AI features. Please install it.")
  }
  
  # Validate provider config
  if (is.null(provider_config$name)) {
    stop("Provider name is required")
  }
  if (is.null(provider_config$api_key) || provider_config$api_key == "") {
    stop("API key is required")
  }
  
  provider_name <- tolower(provider_config$name)
  
  if (provider_name == "openai") {
    base_url <- provider_config$base_url %||% "https://api.openai.com/v1"
    model <- provider_config$model_chat %||% "gpt-3.5-turbo"
    endpoint <- paste0(base_url, "/chat/completions")
    
    body <- list(
      model = model,
      messages = list(
        list(role = "system", content = system),
        list(role = "user", content = user)
      ),
      max_tokens = max_tokens,
      temperature = temperature
    )
    
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "Authorization" = paste("Bearer", provider_config$api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(120)
      
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      assistant_content <- content$choices[[1]]$message$content
      list(content = assistant_content, raw = content)
    }, error = function(e) {
      stop("Chat request failed: ", e$message)
    })
    
    return(result)
    
  } else if (provider_name == "anthropic") {
    base_url <- provider_config$base_url %||% "https://api.anthropic.com/v1"
    model <- provider_config$model_chat %||% "claude-3-sonnet-20240229"
    endpoint <- paste0(base_url, "/messages")
    
    body <- list(
      model = model,
      max_tokens = max_tokens,
      temperature = temperature,
      system = system,
      messages = list(
        list(role = "user", content = user)
      )
    )
    
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "x-api-key" = provider_config$api_key,
          "anthropic-version" = "2023-06-01",
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(120)
      
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      assistant_content <- content$content[[1]]$text
      list(content = assistant_content, raw = content)
    }, error = function(e) {
      stop("Chat request failed: ", e$message)
    })
    
    return(result)
    
  } else if (provider_name == "huggingface") {
    base_url <- provider_config$base_url %||% "https://api-inference.huggingface.co"
    model <- provider_config$model_chat %||% "mistralai/Mistral-7B-Instruct-v0.1"
    endpoint <- paste0(base_url, "/models/", model)
    
    # Construct prompt
    prompt <- paste0(system, "\n\nUser: ", user, "\n\nAssistant:")
    
    body <- list(
      inputs = prompt,
      parameters = list(
        max_new_tokens = max_tokens,
        temperature = temperature,
        return_full_text = FALSE
      )
    )
    
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "Authorization" = paste("Bearer", provider_config$api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(120)
      
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      # HF returns array with generated_text
      if (is.list(content) && length(content) > 0) {
        assistant_content <- content[[1]]$generated_text
      } else {
        assistant_content <- content$generated_text
      }
      
      list(content = assistant_content, raw = content)
    }, error = function(e) {
      stop("Chat request failed: ", e$message)
    })
    
    return(result)
    
  } else if (provider_name == "custom") {
    if (is.null(provider_config$base_url)) {
      stop("base_url is required for custom provider")
    }
    endpoint <- provider_config$base_url
    
    # Assume OpenAI-compatible format
    body <- list(
      model = provider_config$model_chat %||% "default",
      messages = list(
        list(role = "system", content = system),
        list(role = "user", content = user)
      ),
      max_tokens = max_tokens,
      temperature = temperature
    )
    
    result <- tryCatch({
      req <- httr2::request(endpoint) |>
        httr2::req_headers(
          "Authorization" = paste("Bearer", provider_config$api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_retry(max_tries = 3, backoff = ~ 2^.x) |>
        httr2::req_timeout(120)
      
      if (!is.null(provider_config$headers)) {
        for (header_name in names(provider_config$headers)) {
          req <- httr2::req_headers(req, !!header_name := provider_config$headers[[header_name]])
        }
      }
      
      resp <- httr2::req_perform(req)
      content <- httr2::resp_body_json(resp)
      
      # Try to extract content
      if (!is.null(content$choices) && length(content$choices) > 0) {
        assistant_content <- content$choices[[1]]$message$content
      } else if (!is.null(content$content)) {
        assistant_content <- content$content
      } else if (!is.null(content$text)) {
        assistant_content <- content$text
      } else {
        stop("Could not extract content from response")
      }
      
      list(content = assistant_content, raw = content)
    }, error = function(e) {
      stop("Chat request failed: ", e$message)
    })
    
    return(result)
  } else {
    stop("Unsupported provider: ", provider_name, 
         ". Supported providers: openai, anthropic, huggingface, custom")
  }
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
