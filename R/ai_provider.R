#' Send embeddings request to provider
#'
#' @description
#' Generate embeddings for a list of texts using the specified provider.
#' Supports 'dummy' (local bag-of-words), 'openai', and 'huggingface' providers.
#'
#' @param texts Character vector of texts to embed
#' @param provider List with provider configuration including 'name', 'api_key',
#'   'base_url', 'model_embeddings', and 'headers'
#' @param batch_size Integer, number of texts to process in each batch (default: 16)
#' @param max_retries Integer, maximum number of retry attempts (default: 3)
#'
#' @return List of numeric vectors (embeddings) for each input text
#'
#' @examples
#' \dontrun{
#' # Dummy provider (no API key needed)
#' provider <- list(name = "dummy")
#' embeddings <- send_embeddings(c("text 1", "text 2"), provider)
#'
#' # OpenAI provider
#' provider <- list(
#'   name = "openai",
#'   api_key = "sk-...",
#'   base_url = "https://api.openai.com/v1",
#'   model_embeddings = "text-embedding-ada-002"
#' )
#' embeddings <- send_embeddings(c("text 1", "text 2"), provider)
#' }
#'
#' @export
send_embeddings <- function(texts, provider, batch_size = 16, max_retries = 3) {
  if (is.null(provider) || is.null(provider$name)) {
    stop("Provider must be specified with a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  if (provider_name == "dummy") {
    return(send_embeddings_dummy(texts))
  } else if (provider_name == "openai") {
    return(send_embeddings_openai(texts, provider, batch_size, max_retries))
  } else if (provider_name == "huggingface") {
    return(send_embeddings_huggingface(texts, provider, batch_size, max_retries))
  } else {
    stop(sprintf("Unsupported provider: %s. Use 'dummy', 'openai', or 'huggingface'.", provider$name))
  }
}


#' Send chat request to provider
#'
#' @description
#' Send a chat completion request with system and user messages to the specified provider.
#'
#' @param system_msg Character, system message to set context
#' @param user_msg Character, user message/query
#' @param provider List with provider configuration
#' @param max_tokens Integer, maximum tokens in response (default: 800)
#' @param temperature Numeric, sampling temperature 0-1 (default: 0.2)
#' @param max_retries Integer, maximum retry attempts (default: 3)
#'
#' @return List with 'content' (character, response text), 'raw' (full API response)
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' response <- send_chat("You are helpful", "What is R?", provider)
#' cat(response$content)
#' }
#'
#' @export
send_chat <- function(system_msg, user_msg, provider,
                      max_tokens = 800, temperature = 0.2, max_retries = 3) {
  if (is.null(provider) || is.null(provider$name)) {
    stop("Provider must be specified with a 'name' field")
  }
  
  provider_name <- tolower(provider$name)
  
  if (provider_name == "dummy") {
    return(send_chat_dummy(system_msg, user_msg))
  } else if (provider_name == "openai") {
    return(send_chat_openai(system_msg, user_msg, provider, max_tokens, temperature, max_retries))
  } else if (provider_name == "huggingface") {
    return(send_chat_huggingface(system_msg, user_msg, provider, max_tokens, temperature, max_retries))
  } else {
    stop(sprintf("Unsupported provider: %s. Use 'dummy', 'openai', or 'huggingface'.", provider$name))
  }
}


# Internal functions -----------------------------------------------------

#' Dummy embeddings using simple bag-of-words
#' @keywords internal
send_embeddings_dummy <- function(texts) {
  # Build vocabulary from all texts
  all_words <- unique(unlist(strsplit(tolower(paste(texts, collapse = " ")), "\\W+")))
  all_words <- all_words[nchar(all_words) > 0]
  
  # Create embeddings as word frequency vectors
  embeddings <- lapply(texts, function(text) {
    words <- unlist(strsplit(tolower(text), "\\W+"))
    words <- words[nchar(words) > 0]
    
    # Count word frequencies
    vec <- sapply(all_words, function(w) sum(words == w))
    
    # Normalize
    norm <- sqrt(sum(vec^2))
    if (norm > 0) {
      vec <- vec / norm
    }
    
    return(as.numeric(vec))
  })
  
  return(embeddings)
}


#' OpenAI embeddings via API
#' @keywords internal
send_embeddings_openai <- function(texts, provider, batch_size, max_retries) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for OpenAI provider. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for OpenAI provider. Please install it.")
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  api_key <- provider$api_key
  model <- provider$model_embeddings %||% "text-embedding-ada-002"
  
  if (is.null(api_key) || api_key == "") {
    stop("OpenAI provider requires 'api_key' in provider configuration")
  }
  
  endpoint <- paste0(base_url, "/embeddings")
  
  # Process in batches
  all_embeddings <- list()
  num_batches <- ceiling(length(texts) / batch_size)
  
  for (i in seq_len(num_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(texts))
    batch_texts <- texts[start_idx:end_idx]
    
    body <- list(
      input = batch_texts,
      model = model
    )
    
    response <- retry_request(
      function() {
        httr::POST(
          endpoint,
          httr::add_headers(
            "Authorization" = paste("Bearer", api_key),
            "Content-Type" = "application/json"
          ),
          body = jsonlite::toJSON(body, auto_unbox = TRUE),
          encode = "raw"
        )
      },
      max_retries = max_retries
    )
    
    if (httr::status_code(response) != 200) {
      stop(sprintf("OpenAI API error: %s - %s",
                   httr::status_code(response),
                   httr::content(response, "text", encoding = "UTF-8")))
    }
    
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    batch_embeddings <- lapply(result$data, function(x) as.numeric(x$embedding))
    all_embeddings <- c(all_embeddings, batch_embeddings)
  }
  
  return(all_embeddings)
}


#' Hugging Face embeddings via API
#' @keywords internal
send_embeddings_huggingface <- function(texts, provider, batch_size, max_retries) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for Hugging Face provider. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for Hugging Face provider. Please install it.")
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co"
  api_key <- provider$api_key
  model <- provider$model_embeddings %||% "sentence-transformers/all-MiniLM-L6-v2"
  
  if (is.null(api_key) || api_key == "") {
    stop("Hugging Face provider requires 'api_key' in provider configuration")
  }
  
  endpoint <- paste0(base_url, "/pipeline/feature-extraction/", model)
  
  # Process in batches
  all_embeddings <- list()
  num_batches <- ceiling(length(texts) / batch_size)
  
  for (i in seq_len(num_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(texts))
    batch_texts <- texts[start_idx:end_idx]
    
    body <- list(
      inputs = batch_texts
    )
    
    headers <- list(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type" = "application/json"
    )
    
    # Add custom headers if provided
    if (!is.null(provider$headers)) {
      headers <- c(headers, provider$headers)
    }
    
    response <- retry_request(
      function() {
        do.call(httr::POST, c(
          list(
            url = endpoint,
            body = jsonlite::toJSON(body, auto_unbox = TRUE),
            encode = "raw"
          ),
          lapply(names(headers), function(nm) httr::add_headers(.headers = setNames(headers[[nm]], nm)))
        ))
      },
      max_retries = max_retries
    )
    
    if (httr::status_code(response) != 200) {
      stop(sprintf("Hugging Face API error: %s - %s",
                   httr::status_code(response),
                   httr::content(response, "text", encoding = "UTF-8")))
    }
    
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    
    # HF returns different formats depending on model
    if (is.list(result) && !is.null(result[[1]])) {
      batch_embeddings <- lapply(result, function(x) {
        if (is.matrix(x)) {
          # Average pooling for token embeddings
          as.numeric(colMeans(x))
        } else {
          as.numeric(x)
        }
      })
    } else {
      batch_embeddings <- list(as.numeric(result))
    }
    
    all_embeddings <- c(all_embeddings, batch_embeddings)
  }
  
  return(all_embeddings)
}


#' Dummy chat response
#' @keywords internal
send_chat_dummy <- function(system_msg, user_msg) {
  # Simple template-based response for testing
  response_text <- sprintf(
    "Based on the IOBR package documentation, here is a response to your query:\n\n%s\n\nThis is a dummy response generated for offline testing. To get real AI-powered responses, configure a provider such as OpenAI or Hugging Face with valid API credentials.\n\nRelevant context was retrieved from the package documentation to answer your question.",
    substr(user_msg, 1, 100)
  )
  
  return(list(
    content = response_text,
    raw = list(
      model = "dummy",
      system = system_msg,
      user = user_msg,
      timestamp = Sys.time()
    )
  ))
}


#' OpenAI chat completion via API
#' @keywords internal
send_chat_openai <- function(system_msg, user_msg, provider, max_tokens, temperature, max_retries) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for OpenAI provider. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for OpenAI provider. Please install it.")
  }
  
  base_url <- provider$base_url %||% "https://api.openai.com/v1"
  api_key <- provider$api_key
  model <- provider$model_chat %||% "gpt-3.5-turbo"
  
  if (is.null(api_key) || api_key == "") {
    stop("OpenAI provider requires 'api_key' in provider configuration")
  }
  
  endpoint <- paste0(base_url, "/chat/completions")
  
  body <- list(
    model = model,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    max_tokens = max_tokens,
    temperature = temperature
  )
  
  response <- retry_request(
    function() {
      httr::POST(
        endpoint,
        httr::add_headers(
          "Authorization" = paste("Bearer", api_key),
          "Content-Type" = "application/json"
        ),
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        encode = "raw"
      )
    },
    max_retries = max_retries
  )
  
  if (httr::status_code(response) != 200) {
    stop(sprintf("OpenAI API error: %s - %s",
                 httr::status_code(response),
                 httr::content(response, "text", encoding = "UTF-8")))
  }
  
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  return(list(
    content = result$choices[[1]]$message$content,
    raw = result
  ))
}


#' Hugging Face chat completion via API
#' @keywords internal
send_chat_huggingface <- function(system_msg, user_msg, provider, max_tokens, temperature, max_retries) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for Hugging Face provider. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for Hugging Face provider. Please install it.")
  }
  
  base_url <- provider$base_url %||% "https://api-inference.huggingface.co"
  api_key <- provider$api_key
  model <- provider$model_chat %||% "microsoft/DialoGPT-medium"
  
  if (is.null(api_key) || api_key == "") {
    stop("Hugging Face provider requires 'api_key' in provider configuration")
  }
  
  endpoint <- paste0(base_url, "/models/", model)
  
  # Combine system and user messages
  prompt <- sprintf("%s\n\nUser: %s\nAssistant:", system_msg, user_msg)
  
  body <- list(
    inputs = prompt,
    parameters = list(
      max_new_tokens = max_tokens,
      temperature = temperature
    )
  )
  
  headers <- list(
    "Authorization" = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )
  
  if (!is.null(provider$headers)) {
    headers <- c(headers, provider$headers)
  }
  
  response <- retry_request(
    function() {
      do.call(httr::POST, c(
        list(
          url = endpoint,
          body = jsonlite::toJSON(body, auto_unbox = TRUE),
          encode = "raw"
        ),
        lapply(names(headers), function(nm) httr::add_headers(.headers = setNames(headers[[nm]], nm)))
      ))
    },
    max_retries = max_retries
  )
  
  if (httr::status_code(response) != 200) {
    stop(sprintf("Hugging Face API error: %s - %s",
                 httr::status_code(response),
                 httr::content(response, "text", encoding = "UTF-8")))
  }
  
  result <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  
  # Extract generated text (format varies by model)
  content <- if (is.list(result) && !is.null(result[[1]]$generated_text)) {
    result[[1]]$generated_text
  } else if (is.character(result)) {
    result
  } else {
    "Unable to parse response"
  }
  
  return(list(
    content = content,
    raw = result
  ))
}


#' Retry HTTP request with exponential backoff
#' @keywords internal
retry_request <- function(request_fn, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    tryCatch({
      response <- request_fn()
      return(response)
    }, error = function(e) {
      if (attempt == max_retries) {
        stop(sprintf("Request failed after %d attempts: %s", max_retries, e$message))
      }
      
      # Exponential backoff: 1s, 2s, 4s, ...
      wait_time <- 2^(attempt - 1)
      message(sprintf("Request failed (attempt %d/%d). Retrying in %ds...",
                      attempt, max_retries, wait_time))
      Sys.sleep(wait_time)
    })
  }
}


#' NULL-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
