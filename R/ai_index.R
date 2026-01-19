#' Collect package documentation and code files
#'
#' Gathers documentation from README, vignettes, R source files, man pages, and examples.
#'
#' @param pkg_path Character, path to package root (default: current directory)
#'
#' @return List of documents, each with 'content' and 'source' fields
#' @export
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' }
collect_package_docs <- function(pkg_path = ".") {
  docs <- list()
  
  # README.md
  readme_path <- file.path(pkg_path, "README.md")
  if (file.exists(readme_path)) {
    content <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
    docs[[length(docs) + 1]] <- list(content = content, source = "README.md")
  }
  
  # Vignettes
  vignettes_dir <- file.path(pkg_path, "vignettes")
  if (dir.exists(vignettes_dir)) {
    vignette_files <- list.files(vignettes_dir, pattern = "\\.(Rmd|Rnw|md)$", 
                                  full.names = TRUE, recursive = FALSE)
    for (vfile in vignette_files) {
      content <- paste(readLines(vfile, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        content = content, 
        source = file.path("vignettes", basename(vfile))
      )
    }
  }
  
  # R source files
  r_dir <- file.path(pkg_path, "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
    for (rfile in r_files) {
      content <- paste(readLines(rfile, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        content = content,
        source = file.path("R", basename(rfile))
      )
    }
  }
  
  # Man pages (convert .Rd to text)
  man_dir <- file.path(pkg_path, "man")
  if (dir.exists(man_dir)) {
    rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE, recursive = FALSE)
    for (rdfile in rd_files) {
      tryCatch({
        # Use tools::Rd2txt to convert Rd to plain text
        temp_txt <- tempfile(fileext = ".txt")
        tools::Rd2txt(rdfile, out = temp_txt, fragment = FALSE)
        content <- paste(readLines(temp_txt, warn = FALSE), collapse = "\n")
        unlink(temp_txt)
        
        docs[[length(docs) + 1]] <- list(
          content = content,
          source = file.path("man", basename(rdfile))
        )
      }, error = function(e) {
        message(sprintf("Warning: Could not convert %s: %s", basename(rdfile), e$message))
      })
    }
  }
  
  # Examples directory (if exists)
  examples_dir <- file.path(pkg_path, "examples")
  if (dir.exists(examples_dir)) {
    example_files <- list.files(examples_dir, pattern = "\\.(R|r|Rmd|md)$", 
                                 full.names = TRUE, recursive = TRUE)
    for (efile in example_files) {
      content <- paste(readLines(efile, warn = FALSE), collapse = "\n")
      docs[[length(docs) + 1]] <- list(
        content = content,
        source = file.path("examples", basename(efile))
      )
    }
  }
  
  message(sprintf("Collected %d document(s)", length(docs)))
  return(docs)
}

#' Chunk documents into smaller pieces
#'
#' Split documents into chunks by paragraph boundaries with size limits.
#'
#' @param docs List of documents from collect_package_docs
#' @param chunk_size Integer, approximate maximum characters per chunk (default: 800)
#'
#' @return List of chunks, each with 'text' and 'source' fields
#' @export
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' chunks <- chunk_texts(docs, chunk_size = 800)
#' }
chunk_texts <- function(docs, chunk_size = 800) {
  chunks <- list()
  
  for (doc in docs) {
    content <- doc$content
    source <- doc$source
    
    # Split by double newline (paragraphs)
    paragraphs <- strsplit(content, "\n\n+")[[1]]
    paragraphs <- paragraphs[nchar(trimws(paragraphs)) > 0]
    
    current_chunk <- ""
    for (para in paragraphs) {
      para <- trimws(para)
      
      # If paragraph alone exceeds chunk_size, split it further
      if (nchar(para) > chunk_size) {
        # Save current chunk if any
        if (nchar(current_chunk) > 0) {
          chunks[[length(chunks) + 1]] <- list(
            text = trimws(current_chunk),
            source = source
          )
          current_chunk <- ""
        }
        
        # Split large paragraph by sentences or lines
        sentences <- strsplit(para, "(?<=[.!?])\\s+", perl = TRUE)[[1]]
        for (sent in sentences) {
          if (nchar(current_chunk) + nchar(sent) + 1 <= chunk_size) {
            current_chunk <- paste(current_chunk, sent)
          } else {
            if (nchar(current_chunk) > 0) {
              chunks[[length(chunks) + 1]] <- list(
                text = trimws(current_chunk),
                source = source
              )
            }
            current_chunk <- sent
          }
        }
      } else {
        # Add paragraph to current chunk
        if (nchar(current_chunk) + nchar(para) + 2 <= chunk_size) {
          current_chunk <- if (nchar(current_chunk) > 0) {
            paste(current_chunk, para, sep = "\n\n")
          } else {
            para
          }
        } else {
          # Save current chunk and start new one
          if (nchar(current_chunk) > 0) {
            chunks[[length(chunks) + 1]] <- list(
              text = trimws(current_chunk),
              source = source
            )
          }
          current_chunk <- para
        }
      }
    }
    
    # Save remaining chunk
    if (nchar(current_chunk) > 0) {
      chunks[[length(chunks) + 1]] <- list(
        text = trimws(current_chunk),
        source = source
      )
    }
  }
  
  message(sprintf("Created %d chunk(s)", length(chunks)))
  return(chunks)
}

#' Create embeddings index from document chunks
#'
#' Generate embeddings for chunks and save index to disk.
#'
#' @param chunks List of chunks from chunk_texts
#' @param provider Provider configuration list
#' @param index_path Character, path to save index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param batch_size Integer, batch size for embedding generation (default: 16)
#'
#' @return Invisible path to saved index file
#' @export
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' docs <- collect_package_docs(".")
#' chunks <- chunk_texts(docs)
#' create_index(chunks, provider)
#' }
create_index <- function(chunks, provider, index_path = NULL, batch_size = 16) {
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  # Create directory if needed
  index_dir <- dirname(index_path)
  if (!dir.exists(index_dir)) {
    dir.create(index_dir, recursive = TRUE)
  }
  
  message(sprintf("Generating embeddings for %d chunks...", length(chunks)))
  
  # Extract texts
  texts <- sapply(chunks, function(x) x$text)
  
  # Generate embeddings
  embeddings <- send_embeddings(texts, provider, batch_size)
  
  # Build index structure
  entries <- lapply(seq_along(chunks), function(i) {
    list(
      id = i,
      text = chunks[[i]]$text,
      source = chunks[[i]]$source,
      embedding = embeddings[[i]]
    )
  })
  
  # Store vocabulary for dummy provider (optional, for completeness)
  vocab <- if (tolower(provider$name) == "dummy") {
    all_words <- unique(unlist(strsplit(tolower(paste(texts, collapse = " ")), "\\W+")))
    all_words[nchar(all_words) > 0]
  } else {
    NULL
  }
  
  index <- list(
    entries = entries,
    vocab = vocab,
    created_at = Sys.time(),
    provider_meta = list(
      name = provider$name,
      model_embeddings = provider$model_embeddings %||% "default"
    )
  )
  
  # Save index
  saveRDS(index, index_path)
  message(sprintf("Index saved to: %s", index_path))
  
  invisible(index_path)
}

#' Search index for relevant chunks
#'
#' Find most similar chunks to a query using cosine similarity.
#'
#' @param query_text Character, search query
#' @param index List, index structure loaded from RDS file
#' @param provider Provider configuration list
#' @param top_k Integer, number of top results to return (default: 5)
#'
#' @return Data frame with columns: id, source, text, score
#' @export
#'
#' @examples
#' \dontrun{
#' index <- readRDS("inst/ai/iobr_embeddings.rds")
#' provider <- list(name = "dummy")
#' results <- search_index("tumor microenvironment", index, provider, top_k = 5)
#' }
search_index <- function(query_text, index, provider, top_k = 5) {
  # Generate query embedding
  query_embedding <- send_embeddings(query_text, provider, batch_size = 1)[[1]]
  
  # Compute similarities
  similarities <- sapply(index$entries, function(entry) {
    cosine_sim(query_embedding, entry$embedding)
  })
  
  # Get top-k indices
  top_indices <- order(similarities, decreasing = TRUE)[1:min(top_k, length(similarities))]
  
  # Build results
  results <- data.frame(
    id = sapply(top_indices, function(i) index$entries[[i]]$id),
    source = sapply(top_indices, function(i) index$entries[[i]]$source),
    text = sapply(top_indices, function(i) index$entries[[i]]$text),
    score = similarities[top_indices],
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' Compute cosine similarity between two vectors
#'
#' @param vec1 Numeric vector
#' @param vec2 Numeric vector
#'
#' @return Numeric, cosine similarity score
#' @export
#'
#' @examples
#' vec1 <- c(1, 2, 3)
#' vec2 <- c(4, 5, 6)
#' cosine_sim(vec1, vec2)
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

# Utility: NULL-coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
