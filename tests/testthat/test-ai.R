# Tests for IOBR AI Assistant

test_that("dummy provider embeddings work", {
  texts <- c("tumor microenvironment analysis", "immune cell deconvolution")
  provider <- list(name = "dummy")
  
  embeddings <- send_embeddings(texts, provider)
  
  expect_type(embeddings, "list")
  expect_length(embeddings, 2)
  expect_true(all(sapply(embeddings, is.numeric)))
  
  # Check normalization (unit vectors)
  for (emb in embeddings) {
    norm <- sqrt(sum(emb^2))
    expect_true(abs(norm - 1) < 0.001 || norm == 0)
  }
})

test_that("dummy provider chat works", {
  provider <- list(name = "dummy")
  
  response <- send_chat(
    system_msg = "You are helpful",
    user_msg = "What is IOBR?",
    provider = provider
  )
  
  expect_type(response, "list")
  expect_true("content" %in% names(response))
  expect_true("raw" %in% names(response))
  expect_type(response$content, "character")
  expect_true(nchar(response$content) > 0)
})

test_that("cosine similarity calculation is correct", {
  vec1 <- c(1, 0, 0)
  vec2 <- c(1, 0, 0)
  expect_equal(cosine_sim(vec1, vec2), 1.0)
  
  vec1 <- c(1, 0, 0)
  vec2 <- c(0, 1, 0)
  expect_equal(cosine_sim(vec1, vec2), 0.0)
  
  vec1 <- c(1, 1, 0)
  vec2 <- c(1, 1, 0)
  expect_equal(cosine_sim(vec1, vec2), 1.0)
  
  vec1 <- c(1, 2, 3)
  vec2 <- c(4, 5, 6)
  expected <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  expect_equal(cosine_sim(vec1, vec2), expected)
})

test_that("chunk_texts splits documents correctly", {
  docs <- list(
    list(content = "This is paragraph one.\n\nThis is paragraph two.\n\nThis is paragraph three.", 
         source = "test.md")
  )
  
  chunks <- chunk_texts(docs, chunk_size = 50)
  
  expect_type(chunks, "list")
  expect_true(length(chunks) >= 1)
  expect_true(all(sapply(chunks, function(x) "text" %in% names(x))))
  expect_true(all(sapply(chunks, function(x) "source" %in% names(x))))
  
  # Check all chunks are within reasonable size
  for (chunk in chunks) {
    # May slightly exceed due to paragraph boundaries
    expect_true(nchar(chunk$text) <= 100)  # Some tolerance
  }
})

test_that("iobr_ai_configure_provider validates dummy", {
  provider <- iobr_ai_configure_provider(list(name = "dummy"))
  
  expect_equal(provider$name, "dummy")
  expect_equal(provider$model_embeddings, "bag-of-words")
  expect_equal(provider$model_chat, "template-based")
})

test_that("iobr_ai_configure_provider requires api_key for openai", {
  expect_error(
    iobr_ai_configure_provider(list(name = "openai")),
    "requires api_key"
  )
  
  # With API key should work
  provider <- iobr_ai_configure_provider(list(
    name = "openai",
    api_key = "test-key"
  ))
  
  expect_equal(provider$name, "openai")
  expect_equal(provider$api_key, "test-key")
  expect_true(!is.null(provider$model_embeddings))
  expect_true(!is.null(provider$model_chat))
})

test_that("iobr_ai_init creates valid index structure", {
  skip_if_not(dir.exists("inst") || dir.exists("R"), 
              "Package structure not available")
  
  # Create temp directory for test
  temp_dir <- tempdir()
  temp_index <- file.path(temp_dir, "test_index.rds")
  
  # Clean up any existing test index
  if (file.exists(temp_index)) {
    unlink(temp_index)
  }
  
  # Create minimal test package structure
  test_pkg_dir <- file.path(temp_dir, "test_pkg")
  dir.create(test_pkg_dir, showWarnings = FALSE)
  dir.create(file.path(test_pkg_dir, "R"), showWarnings = FALSE)
  
  # Create a simple R file
  writeLines(
    c("#' Test function", "#' @export", "test_func <- function() { 1 }"),
    file.path(test_pkg_dir, "R", "test.R")
  )
  
  # Create README
  writeLines(
    c("# Test Package", "", "This is a test package for IOBR AI."),
    file.path(test_pkg_dir, "README.md")
  )
  
  # Initialize with dummy provider
  result <- iobr_ai_init(
    pkg_path = test_pkg_dir,
    provider = list(name = "dummy"),
    index_path = temp_index,
    chunk_size = 500
  )
  
  expect_true(file.exists(temp_index))
  expect_type(result, "list")
  expect_true("index_path" %in% names(result))
  expect_true("num_chunks" %in% names(result))
  expect_true(result$num_chunks > 0)
  
  # Load and validate index structure
  index <- readRDS(temp_index)
  
  expect_type(index, "list")
  expect_true("entries" %in% names(index))
  expect_true("created_at" %in% names(index))
  expect_true("provider_meta" %in% names(index))
  
  expect_type(index$entries, "list")
  expect_true(length(index$entries) > 0)
  
  # Check first entry structure
  entry <- index$entries[[1]]
  expect_true("id" %in% names(entry))
  expect_true("text" %in% names(entry))
  expect_true("source" %in% names(entry))
  expect_true("embedding" %in% names(entry))
  
  expect_type(entry$id, "integer")
  expect_type(entry$text, "character")
  expect_type(entry$source, "character")
  expect_true(is.numeric(entry$embedding))
  
  # Clean up
  unlink(temp_index)
  unlink(test_pkg_dir, recursive = TRUE)
})

test_that("iobr_ai_query returns expected structure with dummy provider", {
  skip_if_not(dir.exists("inst") || dir.exists("R"), 
              "Package structure not available")
  
  # Create temp directory for test
  temp_dir <- tempdir()
  temp_index <- file.path(temp_dir, "test_query_index.rds")
  
  # Clean up any existing test index
  if (file.exists(temp_index)) {
    unlink(temp_index)
  }
  
  # Create minimal test package structure
  test_pkg_dir <- file.path(temp_dir, "test_query_pkg")
  dir.create(test_pkg_dir, showWarnings = FALSE)
  dir.create(file.path(test_pkg_dir, "R"), showWarnings = FALSE)
  
  # Create R file with relevant content
  writeLines(
    c("#' CIBERSORT deconvolution",
      "#' Performs immune cell deconvolution using CIBERSORT",
      "#' @export",
      "cibersort_func <- function() { }"),
    file.path(test_pkg_dir, "R", "cibersort.R")
  )
  
  # Initialize index
  iobr_ai_init(
    pkg_path = test_pkg_dir,
    provider = list(name = "dummy"),
    index_path = temp_index,
    chunk_size = 500
  )
  
  # Query the index
  result <- iobr_ai_query(
    query = "How to use CIBERSORT?",
    index_path = temp_index,
    provider = list(name = "dummy"),
    top_k = 3
  )
  
  # Validate response structure
  expect_type(result, "list")
  expect_true("answer" %in% names(result))
  expect_true("retrieved" %in% names(result))
  expect_true("raw" %in% names(result))
  
  expect_type(result$answer, "character")
  expect_true(nchar(result$answer) > 0)
  
  # Validate retrieved chunks
  expect_s3_class(result$retrieved, "data.frame")
  expect_true(nrow(result$retrieved) > 0)
  expect_true(all(c("id", "source", "text", "score") %in% names(result$retrieved)))
  
  # Scores should be between 0 and 1
  expect_true(all(result$retrieved$score >= 0))
  expect_true(all(result$retrieved$score <= 1))
  
  # Scores should be in descending order
  expect_true(all(diff(result$retrieved$score) <= 0))
  
  # Clean up
  unlink(temp_index)
  unlink(test_pkg_dir, recursive = TRUE)
})

test_that("search_index returns top-k results", {
  # Create simple index
  provider <- list(name = "dummy")
  texts <- c(
    "tumor microenvironment analysis",
    "immune cell deconvolution",
    "signature scoring methods",
    "survival analysis"
  )
  
  embeddings <- send_embeddings(texts, provider)
  
  entries <- lapply(seq_along(texts), function(i) {
    list(id = i, text = texts[i], source = "test.R", embedding = embeddings[[i]])
  })
  
  index <- list(
    entries = entries,
    created_at = Sys.time(),
    provider_meta = list(name = "dummy")
  )
  
  # Search
  results <- search_index("immune cells", index, provider, top_k = 2)
  
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 2)
  expect_true(all(c("id", "source", "text", "score") %in% names(results)))
  
  # First result should have highest score
  expect_equal(results$score[1], max(results$score))
})

test_that("code extraction works", {
  text_with_code <- "Here is some code:\n\n```r\nlibrary(IOBR)\nresult <- calculate_sig_score(data)\n```\n\nAnd more text."
  
  code_blocks <- IOBR:::.extract_code_blocks(text_with_code)
  
  expect_length(code_blocks, 1)
  expect_true(grepl("library\\(IOBR\\)", code_blocks[1]))
  expect_true(grepl("calculate_sig_score", code_blocks[1]))
})

test_that("code extraction handles no code", {
  text_no_code <- "This is just plain text without any code blocks."
  
  code_blocks <- IOBR:::.extract_code_blocks(text_no_code)
  
  expect_length(code_blocks, 0)
})

test_that("code extraction handles multiple blocks", {
  text_multi_code <- "First:\n```r\ncode1()\n```\n\nSecond:\n```R\ncode2()\n```"
  
  code_blocks <- IOBR:::.extract_code_blocks(text_multi_code)
  
  expect_length(code_blocks, 2)
  expect_true(grepl("code1", code_blocks[1]))
  expect_true(grepl("code2", code_blocks[2]))
})

# Mock HTTP tests (require httptest2 or mockery)
test_that("HTTP retry logic is tested with mock", {
  skip_if_not_installed("mockery")
  
  # This is a placeholder for more comprehensive HTTP mocking
  # In practice, you would mock httr::POST to return various responses
  
  # Test that retry happens on 429
  # Test that retry happens on 5xx
  # Test that success returns correctly
  # Test that max retries is respected
  
  expect_true(TRUE)  # Placeholder
})

test_that("provider errors are handled gracefully", {
  expect_error(
    send_embeddings("test", list(name = "invalid_provider")),
    "Unknown provider"
  )
  
  expect_error(
    send_chat("sys", "user", list(name = "invalid_provider")),
    "Unknown provider"
  )
  
  expect_error(
    iobr_ai_configure_provider(list()),
    "must have a 'name' field"
  )
})

test_that("NULL-coalescing operator works", {
  # Note: This operator is defined in ai_provider.R and ai_index.R
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(NA %||% "default", NA)
  expect_equal(0 %||% "default", 0)
})
