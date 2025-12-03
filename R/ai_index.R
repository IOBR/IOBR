#' Collect Package Documentation
#'
#' Gather documentation from README, vignettes, R files, man pages, and examples.
#'
#' @param pkg_path Character, path to package root (default: current directory)
#'
#' @return List with fields: source (file path), content (text), type (doc type)
#' @export
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' }
collect_package_docs <- function(pkg_path = ".") {
  pkg_path <- normalizePath(pkg_path, mustWork = TRUE)
  docs <- list()
  
  # Collect README
  readme_paths <- c("README.md", "README.Rmd", "README")
  for (readme in readme_paths) {
    readme_file <- file.path(pkg_path, readme)
    if (file.exists(readme_file)) {
      content <- paste(readLines(readme_file, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        source = readme,
        content = content,
        type = "readme"
      )
      break
    }
  }
  
  # Collect vignettes
  vignette_dir <- file.path(pkg_path, "vignettes")
  if (dir.exists(vignette_dir)) {
    vignette_files <- list.files(vignette_dir, pattern = "\\.(Rmd|md)$", 
                                   full.names = TRUE, recursive = FALSE)
    for (vf in vignette_files) {
      content <- paste(readLines(vf, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        source = paste0("vignettes/", basename(vf)),
        content = content,
        type = "vignette"
      )
    }
  }
  
  # Collect R source files
  r_dir <- file.path(pkg_path, "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
    for (rf in r_files) {
      content <- paste(readLines(rf, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        source = paste0("R/", basename(rf)),
        content = content,
        type = "source"
      )
    }
  }
  
  # Collect man pages (convert Rd to text)
  man_dir <- file.path(pkg_path, "man")
  if (dir.exists(man_dir)) {
    man_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE, recursive = FALSE)
    for (mf in man_files) {
      tryCatch({
        # Convert Rd to text
        content <- paste(tools::Rd2txt(mf, out = textConnection("txt", "w", local = TRUE)),
                        collapse = "\n")
        if (exists("txt") && length(txt) > 0) {
          content <- paste(txt, collapse = "\n")
        }
        
        docs[[length(docs) + 1]] <- list(
          source = paste0("man/", basename(mf)),
          content = content,
          type = "man"
        )
      }, error = function(e) {
        # Skip files that can't be converted
        message(sprintf("Warning: Could not convert %s: %s", basename(mf), e$message))
      })
    }
  }
  
  # Collect examples
  examples_dir <- file.path(pkg_path, "examples")
  if (dir.exists(examples_dir)) {
    example_files <- list.files(examples_dir, pattern = "\\.(R|r)$", 
                                 full.names = TRUE, recursive = TRUE)
    for (ef in example_files) {
      content <- paste(readLines(ef, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        source = paste0("examples/", basename(ef)),
        content = content,
        type = "example"
      )
    }
  }
  
  message(sprintf("Collected %d documents from package", length(docs)))
  return(docs)
}

#' Chunk Text Documents
#'
#' Split documents into smaller chunks for embedding and retrieval.
#'
#' @param docs List of documents with source, content, type fields
#' @param chunk_size Integer, approximate character size per chunk (default: 800)
#'
#' @return List of chunks with id, text, source, type fields
#' @export
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' chunks <- chunk_texts(docs, chunk_size = 800)
#' }
chunk_texts <- function(docs, chunk_size = 800) {
  chunks <- list()
  chunk_id <- 1
  
  for (doc in docs) {
    content <- doc$content
    
    # Split by paragraphs (double newlines)
    paragraphs <- strsplit(content, "\n\n+")[[1]]
    
    current_chunk <- ""
    for (para in paragraphs) {
      para <- trimws(para)
      if (nchar(para) == 0) next
      
      # If paragraph itself is larger than chunk_size, split it
      if (nchar(para) > chunk_size * 1.5) {
        # Save current chunk if any
        if (nchar(current_chunk) > 0) {
          chunks[[length(chunks) + 1]] <- list(
            id = chunk_id,
            text = trimws(current_chunk),
            source = doc$source,
            type = doc$type
          )
          chunk_id <- chunk_id + 1
          current_chunk <- ""
        }
        
        # Split large paragraph by sentences
        sentences <- strsplit(para, "(?<=[.!?])\\s+", perl = TRUE)[[1]]
        for (sent in sentences) {
          if (nchar(current_chunk) + nchar(sent) > chunk_size && nchar(current_chunk) > 0) {
            chunks[[length(chunks) + 1]] <- list(
              id = chunk_id,
              text = trimws(current_chunk),
              source = doc$source,
              type = doc$type
            )
            chunk_id <- chunk_id + 1
            current_chunk <- sent
          } else {
            current_chunk <- paste(current_chunk, sent)
          }
        }
      } else {
        # Check if adding this paragraph exceeds chunk_size
        if (nchar(current_chunk) + nchar(para) > chunk_size && nchar(current_chunk) > 0) {
          chunks[[length(chunks) + 1]] <- list(
            id = chunk_id,
            text = trimws(current_chunk),
            source = doc$source,
            type = doc$type
          )
          chunk_id <- chunk_id + 1
          current_chunk <- para
        } else {
          current_chunk <- paste(current_chunk, para, sep = "\n\n")
        }
      }
    }
    
    # Save remaining chunk
    if (nchar(trimws(current_chunk)) > 0) {
      chunks[[length(chunks) + 1]] <- list(
        id = chunk_id,
        text = trimws(current_chunk),
        source = doc$source,
        type = doc$type
      )
      chunk_id <- chunk_id + 1
    }
  }
  
  message(sprintf("Created %d chunks from documents", length(chunks)))
  return(chunks)
}

#' Create Embedding Index
#'
#' Generate embeddings for text chunks and create a searchable index.
#'
#' @param chunks List of text chunks with id, text, source, type fields
#' @param provider List with provider configuration
#' @param index_path Character, path to save index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param batch_size Integer, batch size for embedding requests (default: 16)
#'
#' @return Invisible path to saved index
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' chunks <- chunk_texts(collect_package_docs("."))
#' create_index(chunks, provider)
#' }
create_index <- function(chunks, provider, index_path = NULL, batch_size = 16) {
  if (is.null(index_path)) {
    index_path <- "inst/ai/iobr_embeddings.rds"
  }
  
  # Ensure directory exists
  index_dir <- dirname(index_path)
  if (!dir.exists(index_dir)) {
    dir.create(index_dir, recursive = TRUE)
  }
  
  message("Generating embeddings for chunks...")
  
  # Extract texts
  texts <- sapply(chunks, function(x) x$text)
  
  # Get embeddings
  embeddings <- send_embeddings(texts, provider, batch_size)
  
  # Build index structure
  entries <- list()
  for (i in seq_along(chunks)) {
    entries[[i]] <- list(
      id = chunks[[i]]$id,
      text = chunks[[i]]$text,
      source = chunks[[i]]$source,
      type = chunks[[i]]$type,
      embedding = embeddings[[i]]
    )
  }
  
  # For dummy provider, store vocab for consistency
  vocab <- NULL
  if (provider$name == "dummy") {
    all_text <- paste(texts, collapse = " ")
    vocab <- unique(unlist(strsplit(tolower(all_text), "[^a-z0-9]+")))
    vocab <- vocab[nzchar(vocab)]
  }
  
  index <- list(
    entries = entries,
    vocab = vocab,
    created_at = Sys.time(),
    provider_meta = list(
      name = provider$name,
      model_embeddings = provider$model_embeddings
    )
  )
  
  saveRDS(index, index_path)
  message(sprintf("Index saved to %s (%d entries)", index_path, length(entries)))
  
  invisible(index_path)
}

#' Search Index
#'
#' Search the embedding index for relevant chunks based on query.
#'
#' @param query_text Character, search query
#' @param index List or character path to index RDS file
#' @param provider List with provider configuration
#' @param top_k Integer, number of top results to return (default: 5)
#'
#' @return List of search results with id, source, type, text, score fields
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' results <- search_index("How to analyze TME?", "inst/ai/iobr_embeddings.rds", provider)
#' }
search_index <- function(query_text, index, provider, top_k = 5) {
  # Load index if path provided
  if (is.character(index)) {
    if (!file.exists(index)) {
      stop(sprintf("Index file not found: %s", index))
    }
    index <- readRDS(index)
  }
  
  if (!is.list(index) || is.null(index$entries)) {
    stop("Invalid index structure")
  }
  
  # Get query embedding
  if (provider$name == "dummy" && !is.null(index$vocab)) {
    # Use stored vocab for consistency with index
    words <- strsplit(tolower(query_text), "[^a-z0-9]+")[[1]]
    words <- words[nzchar(words)]
    query_embedding <- sapply(index$vocab, function(v) sum(words == v))
    # Normalize
    norm <- sqrt(sum(query_embedding^2))
    if (norm > 0) query_embedding <- query_embedding / norm
    query_embedding <- list(as.numeric(query_embedding))
  } else {
    query_embedding <- send_embeddings(query_text, provider, batch_size = 1)
  }
  
  query_vec <- query_embedding[[1]]
  
  # Compute cosine similarities
  scores <- sapply(index$entries, function(entry) {
    cosine_sim(query_vec, entry$embedding)
  })
  
  # Get top-k results
  top_indices <- order(scores, decreasing = TRUE)[1:min(top_k, length(scores))]
  
  results <- lapply(top_indices, function(i) {
    entry <- index$entries[[i]]
    list(
      id = entry$id,
      source = entry$source,
      type = entry$type,
      text = entry$text,
      score = scores[i]
    )
  })
  
  return(results)
}

#' Cosine Similarity
#'
#' Calculate cosine similarity between two vectors.
#'
#' @param vec1 Numeric vector
#' @param vec2 Numeric vector
#'
#' @return Numeric similarity score (0 to 1)
#' @export
#'
#' @examples
#' \dontrun{
#' similarity <- cosine_sim(c(1, 2, 3), c(4, 5, 6))
#' }
cosine_sim <- function(vec1, vec2) {
  if (length(vec1) != length(vec2)) {
    stop("Vectors must have the same length")
  }
  
  dot_product <- sum(vec1 * vec2)
  norm1 <- sqrt(sum(vec1^2))
  norm2 <- sqrt(sum(vec2^2))
  
  if (norm1 == 0 || norm2 == 0) {
    return(0)
  }
  
  dot_product / (norm1 * norm2)
}
