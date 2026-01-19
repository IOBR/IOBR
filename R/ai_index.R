#' Collect package documentation
#'
#' @description
#' Collect documentation from various sources in an R package including README,
#' vignettes, R source files, manual pages, and examples.
#'
#' @param pkg_path Character, path to package root directory (default: current directory)
#'
#' @return List of documents, each with 'content' and 'source' fields
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' length(docs)
#' }
#'
#' @export
collect_package_docs <- function(pkg_path = ".") {
  if (!dir.exists(pkg_path)) {
    stop(sprintf("Package path does not exist: %s", pkg_path))
  }
  
  docs <- list()
  
  # Collect README
  readme_files <- c("README.md", "README.Rmd", "README.txt")
  for (readme in readme_files) {
    readme_path <- file.path(pkg_path, readme)
    if (file.exists(readme_path)) {
      content <- tryCatch(
        readLines(readme_path, warn = FALSE),
        error = function(e) NULL
      )
      if (!is.null(content)) {
        docs <- c(docs, list(list(
          content = paste(content, collapse = "\n"),
          source = readme
        )))
        break  # Only include one README
      }
    }
  }
  
  # Collect vignettes
  vignettes_dir <- file.path(pkg_path, "vignettes")
  if (dir.exists(vignettes_dir)) {
    vignette_files <- list.files(vignettes_dir, pattern = "\\.(Rmd|rmd|md|Rnw)$",
                                  full.names = TRUE, recursive = FALSE)
    for (vfile in vignette_files) {
      content <- tryCatch(
        readLines(vfile, warn = FALSE),
        error = function(e) NULL
      )
      if (!is.null(content)) {
        docs <- c(docs, list(list(
          content = paste(content, collapse = "\n"),
          source = sprintf("vignettes/%s", basename(vfile))
        )))
      }
    }
  }
  
  # Collect R source files
  r_dir <- file.path(pkg_path, "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
    for (rfile in r_files) {
      content <- tryCatch(
        readLines(rfile, warn = FALSE),
        error = function(e) NULL
      )
      if (!is.null(content)) {
        docs <- c(docs, list(list(
          content = paste(content, collapse = "\n"),
          source = sprintf("R/%s", basename(rfile))
        )))
      }
    }
  }
  
  # Collect manual pages (convert Rd to text)
  man_dir <- file.path(pkg_path, "man")
  if (dir.exists(man_dir)) {
    rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE, recursive = FALSE)
    for (rdfile in rd_files) {
      content <- tryCatch({
        # Use tools::Rd2txt to convert Rd to plain text
        if (requireNamespace("tools", quietly = TRUE)) {
          tmp_file <- tempfile(fileext = ".txt")
          tools::Rd2txt(rdfile, out = tmp_file)
          txt <- readLines(tmp_file, warn = FALSE)
          unlink(tmp_file)
          paste(txt, collapse = "\n")
        } else {
          # Fallback: read raw Rd
          paste(readLines(rdfile, warn = FALSE), collapse = "\n")
        }
      }, error = function(e) NULL)
      
      if (!is.null(content) && nchar(content) > 0) {
        docs <- c(docs, list(list(
          content = content,
          source = sprintf("man/%s", basename(rdfile))
        )))
      }
    }
  }
  
  # Collect examples directory if it exists
  examples_dir <- file.path(pkg_path, "examples")
  if (dir.exists(examples_dir)) {
    example_files <- list.files(examples_dir, pattern = "\\.(R|r|Rmd|rmd)$",
                                 full.names = TRUE, recursive = TRUE)
    for (efile in example_files) {
      content <- tryCatch(
        readLines(efile, warn = FALSE),
        error = function(e) NULL
      )
      if (!is.null(content)) {
        docs <- c(docs, list(list(
          content = paste(content, collapse = "\n"),
          source = sprintf("examples/%s", basename(efile))
        )))
      }
    }
  }
  
  if (length(docs) == 0) {
    warning("No documentation files found in package")
  }
  
  message(sprintf("Collected %d documentation files", length(docs)))
  return(docs)
}


#' Chunk texts into smaller pieces
#'
#' @description
#' Split documents into smaller chunks by paragraphs with size constraints.
#'
#' @param docs List of documents with 'content' and 'source' fields
#' @param chunk_size Integer, maximum characters per chunk (default: 800)
#'
#' @return List of chunks, each with 'text', 'source', and 'chunk_id' fields
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' chunks <- chunk_texts(docs, chunk_size = 500)
#' }
#'
#' @export
chunk_texts <- function(docs, chunk_size = 800) {
  if (length(docs) == 0) {
    stop("No documents to chunk")
  }
  
  all_chunks <- list()
  chunk_counter <- 1
  
  for (doc in docs) {
    content <- doc$content
    source <- doc$source
    
    # Split by double newlines (paragraphs) first
    paragraphs <- strsplit(content, "\n\n+")[[1]]
    
    # Further split long paragraphs
    for (para in paragraphs) {
      para <- trimws(para)
      if (nchar(para) == 0) next
      
      if (nchar(para) <= chunk_size) {
        # Paragraph fits in one chunk
        all_chunks <- c(all_chunks, list(list(
          text = para,
          source = source,
          chunk_id = chunk_counter
        )))
        chunk_counter <- chunk_counter + 1
      } else {
        # Split long paragraph by sentences
        sentences <- strsplit(para, "(?<=[.!?])\\s+", perl = TRUE)[[1]]
        
        current_chunk <- ""
        for (sent in sentences) {
          if (nchar(current_chunk) + nchar(sent) + 1 <= chunk_size) {
            current_chunk <- if (nchar(current_chunk) > 0) {
              paste(current_chunk, sent)
            } else {
              sent
            }
          } else {
            # Save current chunk if not empty
            if (nchar(current_chunk) > 0) {
              all_chunks <- c(all_chunks, list(list(
                text = current_chunk,
                source = source,
                chunk_id = chunk_counter
              )))
              chunk_counter <- chunk_counter + 1
            }
            
            # Start new chunk
            if (nchar(sent) <= chunk_size) {
              current_chunk <- sent
            } else {
              # Sentence itself is too long, truncate
              current_chunk <- substr(sent, 1, chunk_size)
              all_chunks <- c(all_chunks, list(list(
                text = current_chunk,
                source = source,
                chunk_id = chunk_counter
              )))
              chunk_counter <- chunk_counter + 1
              current_chunk <- ""
            }
          }
        }
        
        # Add remaining chunk
        if (nchar(current_chunk) > 0) {
          all_chunks <- c(all_chunks, list(list(
            text = current_chunk,
            source = source,
            chunk_id = chunk_counter
          )))
          chunk_counter <- chunk_counter + 1
        }
      }
    }
  }
  
  message(sprintf("Created %d chunks from %d documents", length(all_chunks), length(docs)))
  return(all_chunks)
}


#' Create embeddings index
#'
#' @description
#' Create an embeddings index from text chunks by calling the provider's
#' embeddings API and storing the results.
#'
#' @param chunks List of chunks with 'text', 'source', and 'chunk_id' fields
#' @param provider List with provider configuration
#' @param index_path Character, path to save index RDS file (default: inst/ai/iobr_embeddings.rds)
#' @param batch_size Integer, batch size for embedding calls (default: 16)
#'
#' @return List with index metadata including path and number of entries
#'
#' @examples
#' \dontrun{
#' docs <- collect_package_docs(".")
#' chunks <- chunk_texts(docs)
#' provider <- list(name = "dummy")
#' index_meta <- create_index(chunks, provider)
#' }
#'
#' @export
create_index <- function(chunks, provider, index_path = NULL, batch_size = 16) {
  if (length(chunks) == 0) {
    stop("No chunks to index")
  }
  
  if (is.null(provider)) {
    provider <- list(name = "dummy")
  }
  
  # Determine index path
  if (is.null(index_path)) {
    index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
  }
  
  # Create directory if needed
  index_dir <- dirname(index_path)
  if (!dir.exists(index_dir)) {
    dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created directory: %s", index_dir))
  }
  
  # Extract texts
  texts <- sapply(chunks, function(ch) ch$text)
  
  message(sprintf("Generating embeddings for %d chunks...", length(texts)))
  
  # Get embeddings
  embeddings <- send_embeddings(texts, provider, batch_size = batch_size)
  
  if (length(embeddings) != length(chunks)) {
    stop(sprintf("Embeddings count (%d) does not match chunks count (%d)",
                 length(embeddings), length(chunks)))
  }
  
  # Build index structure
  entries <- lapply(seq_along(chunks), function(i) {
    list(
      id = chunks[[i]]$chunk_id,
      text = chunks[[i]]$text,
      source = chunks[[i]]$source,
      embedding = embeddings[[i]]
    )
  })
  
  # Build vocabulary for dummy provider (optional, for reference)
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
      model_embeddings = provider$model_embeddings %||% NA
    )
  )
  
  # Save index
  saveRDS(index, index_path)
  message(sprintf("Index saved to: %s", index_path))
  
  return(list(
    path = index_path,
    num_entries = length(entries),
    created_at = index$created_at,
    provider = provider$name
  ))
}


#' Search embeddings index
#'
#' @description
#' Search the embeddings index for chunks most similar to a query text.
#'
#' @param query_text Character, query string
#' @param index Either a loaded index object or path to index RDS file
#' @param provider List with provider configuration (for computing query embedding)
#' @param top_k Integer, number of top results to return (default: 5)
#'
#' @return List of top-k results with 'id', 'source', 'text', and 'score' fields
#'
#' @examples
#' \dontrun{
#' provider <- list(name = "dummy")
#' results <- search_index("How to analyze TME?", "inst/ai/iobr_embeddings.rds", provider)
#' }
#'
#' @export
search_index <- function(query_text, index, provider, top_k = 5) {
  # Load index if path provided
  if (is.character(index)) {
    if (!file.exists(index)) {
      stop(sprintf("Index file not found: %s", index))
    }
    index <- readRDS(index)
  }
  
  if (is.null(index$entries) || length(index$entries) == 0) {
    stop("Index is empty")
  }
  
  if (is.null(provider)) {
    provider <- list(name = "dummy")
  }
  
  # Get query embedding
  query_embedding <- send_embeddings(query_text, provider)[[1]]
  
  # Compute similarity scores
  scores <- sapply(index$entries, function(entry) {
    cosine_sim(query_embedding, entry$embedding)
  })
  
  # Get top-k indices
  top_indices <- order(scores, decreasing = TRUE)[1:min(top_k, length(scores))]
  
  # Build results
  results <- lapply(top_indices, function(idx) {
    entry <- index$entries[[idx]]
    list(
      id = entry$id,
      source = entry$source,
      text = entry$text,
      score = scores[idx]
    )
  })
  
  return(results)
}


#' Compute cosine similarity between two vectors
#'
#' @param vec1 Numeric vector
#' @param vec2 Numeric vector
#'
#' @return Numeric, cosine similarity score (0 to 1)
#'
#' @keywords internal
cosine_sim <- function(vec1, vec2) {
  if (length(vec1) != length(vec2)) {
    # Handle dimension mismatch by padding with zeros
    max_len <- max(length(vec1), length(vec2))
    if (length(vec1) < max_len) {
      vec1 <- c(vec1, rep(0, max_len - length(vec1)))
    }
    if (length(vec2) < max_len) {
      vec2 <- c(vec2, rep(0, max_len - length(vec2)))
    }
  }
  
  dot_product <- sum(vec1 * vec2)
  norm1 <- sqrt(sum(vec1^2))
  norm2 <- sqrt(sum(vec2^2))
  
  if (norm1 == 0 || norm2 == 0) {
    return(0)
  }
  
  return(dot_product / (norm1 * norm2))
}


#' NULL-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
