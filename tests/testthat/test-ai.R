# Tests for IOBR AI Assistant
# These tests use the dummy provider and mock HTTP requests

library(testthat)
library(IOBR)

# Helper function to create a temporary test directory
create_test_pkg <- function() {
  tmp_dir <- tempfile("test_pkg_")
  dir.create(tmp_dir)
  
  # Create minimal package structure
  dir.create(file.path(tmp_dir, "R"))
  dir.create(file.path(tmp_dir, "man"))
  
  # Create a simple README
  writeLines(
    c("# Test Package", "", "This is a test package for IOBR AI assistant.",
      "", "## Features", "", "- Feature 1: Data analysis", "- Feature 2: Visualization"),
    file.path(tmp_dir, "README.md")
  )
  
  # Create a simple R file
  writeLines(
    c("#' Test Function", "#'", "#' @description A test function for analysis",
      "#' @param x numeric vector", "#' @export",
      "test_fun <- function(x) {", "  mean(x)", "}"),
    file.path(tmp_dir, "R", "test.R")
  )
  
  return(tmp_dir)
}


# Test provider configuration ----------------------------------------------

test_that("dummy provider configuration works", {
  provider <- iobr_ai_configure_provider(list(name = "dummy"))
  
  expect_equal(provider$name, "dummy")
  expect_equal(provider$model_embeddings, "dummy-embeddings")
  expect_equal(provider$model_chat, "dummy-chat")
})

test_that("OpenAI provider gets defaults", {
  provider <- iobr_ai_configure_provider(list(
    name = "openai",
    api_key = "test-key"
  ))
  
  expect_equal(provider$name, "openai")
  expect_equal(provider$base_url, "https://api.openai.com/v1")
  expect_equal(provider$model_embeddings, "text-embedding-ada-002")
  expect_equal(provider$model_chat, "gpt-3.5-turbo")
})

test_that("Hugging Face provider gets defaults", {
  provider <- iobr_ai_configure_provider(list(
    name = "huggingface",
    api_key = "test-key"
  ))
  
  expect_equal(provider$name, "huggingface")
  expect_equal(provider$base_url, "https://api-inference.huggingface.co")
  expect_true(!is.null(provider$model_embeddings))
  expect_true(!is.null(provider$model_chat))
})

test_that("provider configuration requires name", {
  expect_error(
    iobr_ai_configure_provider(list()),
    "must have a 'name' field"
  )
})


# Test embeddings ---------------------------------------------------------

test_that("dummy embeddings work", {
  provider <- list(name = "dummy")
  texts <- c("This is a test", "Another test text", "Final test")
  
  embeddings <- send_embeddings(texts, provider)
  
  expect_equal(length(embeddings), 3)
  expect_true(is.numeric(embeddings[[1]]))
  expect_true(length(embeddings[[1]]) > 0)
  
  # Check normalization
  norm1 <- sqrt(sum(embeddings[[1]]^2))
  expect_true(abs(norm1 - 1.0) < 0.01 || norm1 == 0)
})

test_that("embeddings handle empty strings", {
  provider <- list(name = "dummy")
  texts <- c("test", "", "another")
  
  expect_error(
    embeddings <- send_embeddings(texts, provider),
    NA  # Should not error
  )
})


# Test chat ---------------------------------------------------------------

test_that("dummy chat works", {
  provider <- list(name = "dummy")
  
  response <- send_chat(
    "You are helpful",
    "What is R?",
    provider
  )
  
  expect_true("content" %in% names(response))
  expect_true("raw" %in% names(response))
  expect_true(is.character(response$content))
  expect_true(nchar(response$content) > 0)
})


# Test document collection ------------------------------------------------

test_that("collect_package_docs works", {
  tmp_dir <- create_test_pkg()
  
  docs <- collect_package_docs(tmp_dir)
  
  expect_true(length(docs) > 0)
  expect_true("content" %in% names(docs[[1]]))
  expect_true("source" %in% names(docs[[1]]))
  
  # Check that README was collected
  readme_found <- any(sapply(docs, function(d) grepl("README", d$source)))
  expect_true(readme_found)
  
  unlink(tmp_dir, recursive = TRUE)
})

test_that("collect_package_docs handles missing directory", {
  expect_error(
    collect_package_docs("/nonexistent/path"),
    "does not exist"
  )
})


# Test chunking -----------------------------------------------------------

test_that("chunk_texts works", {
  docs <- list(
    list(content = "This is a test document.\n\nWith multiple paragraphs.\n\nAnd more content.",
         source = "test.md")
  )
  
  chunks <- chunk_texts(docs, chunk_size = 50)
  
  expect_true(length(chunks) > 0)
  expect_true("text" %in% names(chunks[[1]]))
  expect_true("source" %in% names(chunks[[1]]))
  expect_true("chunk_id" %in% names(chunks[[1]]))
  
  # Check chunk sizes
  for (chunk in chunks) {
    expect_true(nchar(chunk$text) <= 50 + 10)  # Small tolerance
  }
})

test_that("chunk_texts handles long paragraphs", {
  long_text <- paste(rep("word", 200), collapse = " ")
  docs <- list(list(content = long_text, source = "long.txt"))
  
  chunks <- chunk_texts(docs, chunk_size = 100)
  
  expect_true(length(chunks) > 1)
  expect_true(all(sapply(chunks, function(c) nchar(c$text) <= 100 + 10)))
})


# Test index creation and search ------------------------------------------

test_that("create_index works with dummy provider", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  docs <- collect_package_docs(tmp_dir)
  chunks <- chunk_texts(docs, chunk_size = 200)
  provider <- list(name = "dummy")
  
  index_meta <- create_index(chunks, provider, index_path = tmp_index)
  
  expect_true(file.exists(tmp_index))
  expect_equal(index_meta$num_entries, length(chunks))
  expect_equal(index_meta$provider, "dummy")
  
  # Load and verify index structure
  index <- readRDS(tmp_index)
  expect_true("entries" %in% names(index))
  expect_true("created_at" %in% names(index))
  expect_true("provider_meta" %in% names(index))
  expect_equal(length(index$entries), length(chunks))
  
  # Verify entry structure
  entry <- index$entries[[1]]
  expect_true("id" %in% names(entry))
  expect_true("text" %in% names(entry))
  expect_true("source" %in% names(entry))
  expect_true("embedding" %in% names(entry))
  expect_true(is.numeric(entry$embedding))
  
  unlink(tmp_dir, recursive = TRUE)
  unlink(tmp_index)
})

test_that("search_index works", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  docs <- collect_package_docs(tmp_dir)
  chunks <- chunk_texts(docs, chunk_size = 200)
  provider <- list(name = "dummy")
  
  create_index(chunks, provider, index_path = tmp_index)
  
  # Search
  results <- search_index("data analysis", tmp_index, provider, top_k = 3)
  
  expect_equal(length(results), min(3, length(chunks)))
  expect_true("text" %in% names(results[[1]]))
  expect_true("source" %in% names(results[[1]]))
  expect_true("score" %in% names(results[[1]]))
  expect_true("id" %in% names(results[[1]]))
  
  # Scores should be between 0 and 1
  for (result in results) {
    expect_true(result$score >= 0 && result$score <= 1)
  }
  
  # Results should be sorted by score (descending)
  scores <- sapply(results, function(r) r$score)
  expect_equal(scores, sort(scores, decreasing = TRUE))
  
  unlink(tmp_dir, recursive = TRUE)
  unlink(tmp_index)
})


# Test cosine similarity --------------------------------------------------

test_that("cosine_sim works", {
  vec1 <- c(1, 0, 0)
  vec2 <- c(1, 0, 0)
  expect_equal(cosine_sim(vec1, vec2), 1.0)
  
  vec1 <- c(1, 0, 0)
  vec2 <- c(0, 1, 0)
  expect_equal(cosine_sim(vec1, vec2), 0.0)
  
  vec1 <- c(1, 1, 0)
  vec2 <- c(1, 1, 0)
  expect_equal(cosine_sim(vec1, vec2), 1.0)
})

test_that("cosine_sim handles zero vectors", {
  vec1 <- c(0, 0, 0)
  vec2 <- c(1, 1, 1)
  expect_equal(cosine_sim(vec1, vec2), 0)
})

test_that("cosine_sim handles dimension mismatch", {
  vec1 <- c(1, 2)
  vec2 <- c(1, 2, 3)
  
  # Should handle gracefully by padding
  result <- cosine_sim(vec1, vec2)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


# Test high-level API -----------------------------------------------------

test_that("iobr_ai_init creates valid index", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  old_wd <- getwd()
  setwd(tmp_dir)
  
  result <- iobr_ai_init(
    pkg_path = ".",
    provider = list(name = "dummy"),
    index_path = tmp_index,
    chunk_size = 200
  )
  
  expect_true(file.exists(tmp_index))
  expect_true("num_entries" %in% names(result))
  expect_true(result$num_entries > 0)
  expect_equal(result$provider, "dummy")
  
  setwd(old_wd)
  unlink(tmp_dir, recursive = TRUE)
  unlink(tmp_index)
})

test_that("iobr_ai_query returns expected structure", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  old_wd <- getwd()
  setwd(tmp_dir)
  
  # Create index
  iobr_ai_init(
    pkg_path = ".",
    provider = list(name = "dummy"),
    index_path = tmp_index
  )
  
  # Query
  response <- iobr_ai_query(
    query = "How do I use test_fun?",
    index_path = tmp_index,
    provider = list(name = "dummy"),
    top_k = 3
  )
  
  expect_true("answer" %in% names(response))
  expect_true("retrieved" %in% names(response))
  expect_true("raw" %in% names(response))
  expect_true(is.character(response$answer))
  expect_true(length(response$retrieved) > 0)
  
  # Check retrieved structure
  retrieved_item <- response$retrieved[[1]]
  expect_true("text" %in% names(retrieved_item))
  expect_true("source" %in% names(retrieved_item))
  expect_true("score" %in% names(retrieved_item))
  
  setwd(old_wd)
  unlink(tmp_dir, recursive = TRUE)
  unlink(tmp_index)
})

test_that("iobr_ai_query requires index", {
  tmp_index <- tempfile(fileext = ".rds")
  
  expect_error(
    iobr_ai_query("test query", index_path = tmp_index),
    "Index not found"
  )
})

test_that("iobr_ai_list_index works", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  old_wd <- getwd()
  setwd(tmp_dir)
  
  # Create index
  iobr_ai_init(
    pkg_path = ".",
    provider = list(name = "dummy"),
    index_path = tmp_index
  )
  
  # List index
  info <- iobr_ai_list_index(index_path = tmp_index)
  
  expect_true(!is.null(info))
  expect_true("num_entries" %in% names(info))
  expect_true("created_at" %in% names(info))
  expect_true("provider" %in% names(info))
  
  setwd(old_wd)
  unlink(tmp_dir, recursive = TRUE)
  unlink(tmp_index)
})

test_that("iobr_ai_reset_index works", {
  tmp_index <- tempfile(fileext = ".rds")
  
  # Create a dummy index file
  saveRDS(list(entries = list()), tmp_index)
  expect_true(file.exists(tmp_index))
  
  # Reset
  result <- iobr_ai_reset_index(index_path = tmp_index)
  
  expect_true(result)
  expect_false(file.exists(tmp_index))
})


# Test code extraction ----------------------------------------------------

test_that("extract_code_blocks works", {
  text <- "Here is some code:\n\n```r\nx <- 1:10\nmean(x)\n```\n\nAnd more text."
  
  # Access internal function (normally not exported)
  extract_fn <- IOBR:::extract_code_blocks
  
  blocks <- extract_fn(text)
  
  expect_equal(length(blocks), 1)
  expect_true(grepl("x <- 1:10", blocks[[1]]))
  expect_true(grepl("mean\\(x\\)", blocks[[1]]))
})

test_that("extract_code_blocks handles multiple blocks", {
  text <- "Code 1:\n```\ncode1\n```\n\nCode 2:\n```r\ncode2\n```"
  
  extract_fn <- IOBR:::extract_code_blocks
  blocks <- extract_fn(text)
  
  expect_equal(length(blocks), 2)
})

test_that("extract_code_blocks handles no code", {
  text <- "Just plain text with no code blocks"
  
  extract_fn <- IOBR:::extract_code_blocks
  blocks <- extract_fn(text)
  
  expect_equal(length(blocks), 0)
})


# Integration test --------------------------------------------------------

test_that("full workflow works end-to-end", {
  tmp_dir <- create_test_pkg()
  tmp_index <- tempfile(fileext = ".rds")
  
  old_wd <- getwd()
  setwd(tmp_dir)
  
  # Step 1: Initialize
  provider <- iobr_ai_configure_provider(list(name = "dummy"))
  index_info <- iobr_ai_init(
    pkg_path = ".",
    provider = provider,
    index_path = tmp_index
  )
  
  expect_true(file.exists(tmp_index))
  
  # Step 2: Query
  response <- iobr_ai_query(
    query = "What features does this package have?",
    index_path = tmp_index,
    provider = provider
  )
  
  expect_true(nchar(response$answer) > 0)
  expect_true(length(response$retrieved) > 0)
  
  # Step 3: List index
  info <- iobr_ai_list_index(index_path = tmp_index)
  expect_equal(info$num_entries, index_info$num_entries)
  
  # Step 4: Reset
  reset_result <- iobr_ai_reset_index(index_path = tmp_index)
  expect_true(reset_result)
  expect_false(file.exists(tmp_index))
  
  setwd(old_wd)
  unlink(tmp_dir, recursive = TRUE)
})
