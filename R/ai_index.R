#' AI Indexing and Retrieval Functions
#'
#' Internal functions for collecting package content, chunking, and vector search
#'
#' @name ai_index
#' @keywords internal
NULL

#' Collect package documentation and code
#'
#' @param pkg_path Path to package root
#'
#' @return List of content items with text and source
#' @keywords internal
collect_package_content <- function(pkg_path = ".") {
  content <- list()
  
  # 1. README
  readme_path <- file.path(pkg_path, "README.md")
  if (file.exists(readme_path)) {
    readme_text <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
    content[[length(content) + 1]] <- list(
      text = readme_text,
      source = "README.md",
      type = "readme"
    )
  }
  
  # 2. Vignettes
  vignette_dir <- file.path(pkg_path, "vignettes")
  if (dir.exists(vignette_dir)) {
    vignette_files <- list.files(vignette_dir, pattern = "\\.Rmd$", full.names = TRUE)
    for (vf in vignette_files) {
      vignette_text <- paste(readLines(vf, warn = FALSE), collapse = "\n")
      content[[length(content) + 1]] <- list(
        text = vignette_text,
        source = file.path("vignettes", basename(vf)),
        type = "vignette"
      )
    }
  }
  
  # 3. Man pages (Rd files)
  man_dir <- file.path(pkg_path, "man")
  if (dir.exists(man_dir)) {
    rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
    for (rd in rd_files) {
      # Simple Rd to text conversion
      rd_text <- rd_to_text(rd)
      content[[length(content) + 1]] <- list(
        text = rd_text,
        source = file.path("man", basename(rd)),
        type = "documentation"
      )
    }
  }
  
  # 4. R code files
  r_dir <- file.path(pkg_path, "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    for (rf in r_files) {
      r_text <- paste(readLines(rf, warn = FALSE), collapse = "\n")
      content[[length(content) + 1]] <- list(
        text = r_text,
        source = file.path("R", basename(rf)),
        type = "code"
      )
    }
  }
  
  return(content)
}

#' Convert Rd file to plain text
#'
#' @param rd_path Path to Rd file
#'
#' @return Plain text representation
#' @keywords internal
rd_to_text <- function(rd_path) {
  rd_lines <- readLines(rd_path, warn = FALSE)
  
  # Remove LaTeX-style commands and simplify
  text <- rd_lines
  text <- gsub("\\\\name\\{([^}]+)\\}", "Function: \\1", text)
  text <- gsub("\\\\alias\\{([^}]+)\\}", "", text)
  text <- gsub("\\\\title\\{([^}]+)\\}", "Title: \\1", text)
  text <- gsub("\\\\description\\{", "Description: ", text)
  text <- gsub("\\\\usage\\{", "Usage: ", text)
  text <- gsub("\\\\arguments\\{", "Arguments: ", text)
  text <- gsub("\\\\value\\{", "Returns: ", text)
  text <- gsub("\\\\examples\\{", "Examples: ", text)
  text <- gsub("\\\\item\\{([^}]+)\\}\\{([^}]+)\\}", "- \\1: \\2", text)
  text <- gsub("\\\\[a-zA-Z]+\\{([^}]+)\\}", "\\1", text)
  text <- gsub("\\\\[a-zA-Z]+", "", text)
  text <- gsub("^%.*$", "", text)  # Remove comments
  
  # Clean up
  text <- text[text != ""]
  text <- paste(text, collapse = "\n")
  
  return(text)
}

#' Chunk text into smaller pieces
#'
#' @param text Text to chunk
#' @param chunk_size Maximum characters per chunk
#' @param overlap Number of characters to overlap between chunks
#'
#' @return Character vector of chunks
#' @keywords internal
chunk_text <- function(text, chunk_size = 800, overlap = 100) {
  if (nchar(text) <= chunk_size) {
    return(text)
  }
  
  # Split by paragraphs first
  paragraphs <- strsplit(text, "\n\n+", perl = TRUE)[[1]]
  
  chunks <- character()
  current_chunk <- ""
  
  for (para in paragraphs) {
    # If adding this paragraph would exceed chunk size
    if (nchar(current_chunk) + nchar(para) > chunk_size) {
      if (current_chunk != "") {
        chunks <- c(chunks, current_chunk)
        # Start new chunk with overlap from previous
        if (overlap > 0 && nchar(current_chunk) > overlap) {
          overlap_text <- substr(current_chunk, 
                                 nchar(current_chunk) - overlap + 1, 
                                 nchar(current_chunk))
          current_chunk <- paste0(overlap_text, "\n\n", para)
        } else {
          current_chunk <- para
        }
      } else {
        # Paragraph itself is too long, split it
        if (nchar(para) > chunk_size) {
          words <- strsplit(para, " +")[[1]]
          temp_chunk <- ""
          for (word in words) {
            if (nchar(temp_chunk) + nchar(word) + 1 > chunk_size) {
              if (temp_chunk != "") {
                chunks <- c(chunks, temp_chunk)
              }
              temp_chunk <- word
            } else {
              temp_chunk <- if (temp_chunk == "") word else paste(temp_chunk, word)
            }
          }
          current_chunk <- temp_chunk
        } else {
          current_chunk <- para
        }
      }
    } else {
      current_chunk <- if (current_chunk == "") para else paste0(current_chunk, "\n\n", para)
    }
  }
  
  # Add remaining chunk
  if (current_chunk != "") {
    chunks <- c(chunks, current_chunk)
  }
  
  return(chunks)
}

#' Create embeddings for content chunks
#'
#' @param content_items List of content items from collect_package_content
#' @param provider_config Provider configuration
#' @param chunk_size Maximum characters per chunk
#' @param progress_callback Optional callback function for progress updates
#'
#' @return List of index entries with id, text, embedding, source
#' @keywords internal
create_index <- function(content_items, provider_config, chunk_size = 800, 
                        progress_callback = NULL) {
  index <- list()
  total_items <- length(content_items)
  
  if (!is.null(progress_callback)) {
    progress_callback(sprintf("Processing %d content items...", total_items))
  }
  
  for (i in seq_along(content_items)) {
    item <- content_items[[i]]
    
    if (!is.null(progress_callback)) {
      progress_callback(sprintf("Processing %d/%d: %s", i, total_items, item$source))
    }
    
    # Chunk the content
    chunks <- chunk_text(item$text, chunk_size = chunk_size)
    
    # Create embeddings in batch (up to 20 at a time to avoid rate limits)
    batch_size <- 20
    for (batch_start in seq(1, length(chunks), by = batch_size)) {
      batch_end <- min(batch_start + batch_size - 1, length(chunks))
      batch_chunks <- chunks[batch_start:batch_end]
      
      # Retry logic for embeddings
      max_retries <- 3
      retry_count <- 0
      embeddings <- NULL
      
      while (retry_count < max_retries && is.null(embeddings)) {
        tryCatch({
          embeddings <- send_embedding(batch_chunks, provider_config)
          
          # Add to index
          for (j in seq_along(batch_chunks)) {
            chunk_idx <- batch_start + j - 1
            index[[length(index) + 1]] <- list(
              id = sprintf("%s_chunk_%d", item$source, chunk_idx),
              text = batch_chunks[[j]],
              embedding = embeddings[[j]],
              source = item$source,
              type = item$type
            )
          }
        }, error = function(e) {
          retry_count <<- retry_count + 1
          if (retry_count >= max_retries) {
            warning(sprintf("Failed to embed chunks from %s after %d retries: %s", 
                          item$source, max_retries, e$message))
          } else {
            # Exponential backoff
            Sys.sleep(2^retry_count)
          }
        })
      }
    }
  }
  
  if (!is.null(progress_callback)) {
    progress_callback(sprintf("Index created with %d entries", length(index)))
  }
  
  return(index)
}

#' Calculate cosine similarity between two vectors
#'
#' @param v1 First vector
#' @param v2 Second vector
#'
#' @return Cosine similarity score
#' @keywords internal
cosine_similarity <- function(v1, v2) {
  if (length(v1) != length(v2)) {
    stop("Vectors must have the same length")
  }
  
  dot_product <- sum(v1 * v2)
  norm_v1 <- sqrt(sum(v1^2))
  norm_v2 <- sqrt(sum(v2^2))
  
  if (norm_v1 == 0 || norm_v2 == 0) {
    return(0)
  }
  
  similarity <- dot_product / (norm_v1 * norm_v2)
  return(similarity)
}

#' Search index with query
#'
#' @param query_text Query text
#' @param index Index list
#' @param provider_config Provider configuration
#' @param top_k Number of top results to return
#'
#' @return List of top_k results with text, source, and similarity score
#' @keywords internal
search_index <- function(query_text, index, provider_config, top_k = 5) {
  # Get query embedding
  query_embedding <- send_embedding(query_text, provider_config)[[1]]
  
  # Calculate similarities
  similarities <- sapply(index, function(entry) {
    cosine_similarity(query_embedding, entry$embedding)
  })
  
  # Get top k indices
  top_indices <- order(similarities, decreasing = TRUE)[1:min(top_k, length(similarities))]
  
  # Build results
  results <- lapply(top_indices, function(idx) {
    list(
      text = index[[idx]]$text,
      source = index[[idx]]$source,
      type = index[[idx]]$type,
      similarity = similarities[idx],
      id = index[[idx]]$id
    )
  })
  
  return(results)
}

#' Save index to disk
#'
#' @param index Index list
#' @param path Path to save RDS file
#'
#' @keywords internal
save_index <- function(index, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(index, path)
}

#' Load index from disk
#'
#' @param path Path to RDS file
#'
#' @return Index list
#' @keywords internal
load_index <- function(path) {
  if (!file.exists(path)) {
    stop("Index file not found: ", path)
  }
  readRDS(path)
}
