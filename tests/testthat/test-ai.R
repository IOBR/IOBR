# Tests for IOBR AI Assistant
# These tests use the dummy provider to avoid requiring API keys

library(testthat)
library(IOBR)

test_that("send_embeddings works with dummy provider", {
  provider <- list(name = "dummy")
  texts <- c("hello world", "foo bar baz")
  
  embeddings <- send_embeddings(texts, provider)
  
  expect_type(embeddings, "list")
  expect_length(embeddings, 2)
  expect_type(embeddings[[1]], "double")
  expect_type(embeddings[[2]], "double")
})

test_that("send_chat works with dummy provider", {
  provider <- list(name = "dummy")
  
  response <- send_chat(
    system_msg = "You are helpful",
    user_msg = "Hello",
    provider = provider
  )
  
  expect_type(response, "list")
  expect_true("content" %in% names(response))
  expect_true("raw" %in% names(response))
  expect_type(response$content, "character")
  expect_gt(nchar(response$content), 0)
})

test_that("collect_package_docs gathers files", {
  # Create a temp package structure
  temp_dir <- tempdir()
  pkg_dir <- file.path(temp_dir, "testpkg")
  dir.create(pkg_dir, showWarnings = FALSE)
  dir.create(file.path(pkg_dir, "R"), showWarnings = FALSE)
  
  # Create a simple README
  writeLines("# Test Package\n\nThis is a test.", file.path(pkg_dir, "README.md"))
  
  # Create a simple R file
  writeLines("# A function\ntest_func <- function() { 1 }", 
             file.path(pkg_dir, "R", "test.R"))
  
  docs <- collect_package_docs(pkg_dir)
  
  expect_type(docs, "list")
  expect_gt(length(docs), 0)
  expect_true(all(sapply(docs, function(d) "source" %in% names(d))))
  expect_true(all(sapply(docs, function(d) "content" %in% names(d))))
  expect_true(all(sapply(docs, function(d) "type" %in% names(d))))
  
  # Cleanup
  unlink(pkg_dir, recursive = TRUE)
})

test_that("chunk_texts splits documents appropriately", {
  docs <- list(
    list(
      source = "test.txt",
      content = paste(rep("This is a test paragraph.", 50), collapse = "\n\n"),
      type = "test"
    )
  )
  
  chunks <- chunk_texts(docs, chunk_size = 100)
  
  expect_type(chunks, "list")
  expect_gt(length(chunks), 1)
  
  for (chunk in chunks) {
    expect_true("id" %in% names(chunk))
    expect_true("text" %in% names(chunk))
    expect_true("source" %in% names(chunk))
    expect_true("type" %in% names(chunk))
  }
})

test_that("create_index generates valid index structure", {
  provider <- list(name = "dummy")
  
  chunks <- list(
    list(id = 1, text = "chunk one", source = "test.R", type = "source"),
    list(id = 2, text = "chunk two", source = "test.R", type = "source")
  )
  
  temp_index <- tempfile(fileext = ".rds")
  
  result <- create_index(chunks, provider, index_path = temp_index)
  
  expect_true(file.exists(temp_index))
  
  # Load and validate index
  index <- readRDS(temp_index)
  
  expect_type(index, "list")
  expect_true("entries" %in% names(index))
  expect_true("created_at" %in% names(index))
  expect_true("provider_meta" %in% names(index))
  
  expect_length(index$entries, 2)
  
  for (entry in index$entries) {
    expect_true("id" %in% names(entry))
    expect_true("text" %in% names(entry))
    expect_true("source" %in% names(entry))
    expect_true("type" %in% names(entry))
    expect_true("embedding" %in% names(entry))
    expect_type(entry$embedding, "double")
  }
  
  # Cleanup
  unlink(temp_index)
})

test_that("search_index returns relevant results", {
  provider <- list(name = "dummy")
  
  chunks <- list(
    list(id = 1, text = "machine learning algorithms", source = "test.R", type = "source"),
    list(id = 2, text = "data visualization plots", source = "test.R", type = "source"),
    list(id = 3, text = "neural network training", source = "test.R", type = "source")
  )
  
  temp_index <- tempfile(fileext = ".rds")
  create_index(chunks, provider, index_path = temp_index)
  
  results <- search_index("machine learning", temp_index, provider, top_k = 2)
  
  expect_type(results, "list")
  expect_length(results, 2)
  
  for (result in results) {
    expect_true("id" %in% names(result))
    expect_true("source" %in% names(result))
    expect_true("type" %in% names(result))
    expect_true("text" %in% names(result))
    expect_true("score" %in% names(result))
    expect_type(result$score, "double")
    expect_gte(result$score, 0)
    expect_lte(result$score, 1)
  }
  
  # Cleanup
  unlink(temp_index)
})

test_that("cosine_sim computes correct similarity", {
  vec1 <- c(1, 0, 0)
  vec2 <- c(1, 0, 0)
  vec3 <- c(0, 1, 0)
  
  # Identical vectors
  expect_equal(cosine_sim(vec1, vec2), 1)
  
  # Orthogonal vectors
  expect_equal(cosine_sim(vec1, vec3), 0)
  
  # General case
  vec4 <- c(1, 1, 0)
  sim <- cosine_sim(vec1, vec4)
  expect_gt(sim, 0)
  expect_lt(sim, 1)
})

test_that("iobr_ai_configure_provider validates and sets defaults", {
  # Dummy provider
  provider <- iobr_ai_configure_provider(list(name = "dummy"))
  expect_equal(provider$name, "dummy")
  
  # OpenAI provider
  provider <- iobr_ai_configure_provider(list(
    name = "openai",
    api_key = "test-key"
  ))
  expect_equal(provider$name, "openai")
  expect_equal(provider$base_url, "https://api.openai.com/v1")
  expect_equal(provider$model_embeddings, "text-embedding-ada-002")
  expect_equal(provider$model_chat, "gpt-3.5-turbo")
  
  # Hugging Face provider
  provider <- iobr_ai_configure_provider(list(
    name = "huggingface",
    api_key = "test-key"
  ))
  expect_equal(provider$name, "huggingface")
  expect_true(!is.null(provider$base_url))
  expect_true(!is.null(provider$model_embeddings))
  expect_true(!is.null(provider$model_chat))
  
  # Invalid provider
  expect_error(iobr_ai_configure_provider("not a list"))
  expect_error(iobr_ai_configure_provider(list()))
})

test_that("iobr_ai_init creates index with metadata", {
  skip_if_not(file.exists("DESCRIPTION"), "Not in package root")
  
  provider <- list(name = "dummy")
  
  # Create in temp location
  temp_index <- tempfile(fileext = ".rds")
  
  metadata <- iobr_ai_init(
    pkg_path = ".",
    provider = provider,
    index_path = temp_index,
    chunk_size = 500
  )
  
  expect_type(metadata, "list")
  expect_true("index_path" %in% names(metadata))
  expect_true("provider" %in% names(metadata))
  expect_true("num_chunks" %in% names(metadata))
  expect_true("num_docs" %in% names(metadata))
  expect_true("created_at" %in% names(metadata))
  
  expect_true(file.exists(temp_index))
  
  # Cleanup
  unlink(temp_index)
})

test_that("iobr_ai_query returns structured response", {
  provider <- list(name = "dummy")
  
  # Create a minimal index
  chunks <- list(
    list(id = 1, text = "IOBR is a package for TME analysis", source = "README.md", type = "readme"),
    list(id = 2, text = "Use CIBERSORT for deconvolution", source = "vignette.Rmd", type = "vignette")
  )
  
  temp_index <- tempfile(fileext = ".rds")
  create_index(chunks, provider, index_path = temp_index)
  
  response <- iobr_ai_query(
    query = "What is IOBR?",
    index_path = temp_index,
    provider = provider,
    top_k = 2
  )
  
  expect_type(response, "list")
  expect_true("answer" %in% names(response))
  expect_true("code" %in% names(response))
  expect_true("retrieved" %in% names(response))
  expect_true("raw" %in% names(response))
  
  expect_type(response$answer, "character")
  expect_gt(nchar(response$answer), 0)
  
  expect_type(response$retrieved, "list")
  expect_length(response$retrieved, 2)
  
  # Cleanup
  unlink(temp_index)
})

test_that("iobr_ai_list_index shows index info", {
  provider <- list(name = "dummy")
  
  chunks <- list(
    list(id = 1, text = "test chunk", source = "test.R", type = "source")
  )
  
  temp_index <- tempfile(fileext = ".rds")
  create_index(chunks, provider, index_path = temp_index)
  
  info <- iobr_ai_list_index(temp_index)
  
  expect_type(info, "list")
  expect_true("index_path" %in% names(info))
  expect_true("num_entries" %in% names(info))
  expect_true("created_at" %in% names(info))
  expect_true("provider" %in% names(info))
  
  # Cleanup
  unlink(temp_index)
})

test_that("iobr_ai_reset_index deletes index", {
  provider <- list(name = "dummy")
  
  chunks <- list(
    list(id = 1, text = "test chunk", source = "test.R", type = "source")
  )
  
  temp_index <- tempfile(fileext = ".rds")
  create_index(chunks, provider, index_path = temp_index)
  
  expect_true(file.exists(temp_index))
  
  success <- iobr_ai_reset_index(temp_index)
  
  expect_true(success)
  expect_false(file.exists(temp_index))
})

# Mock HTTP tests (if httptest2 is available)
if (requireNamespace("httptest2", quietly = TRUE)) {
  
  test_that("OpenAI embeddings can be mocked", {
    skip("Mock HTTP tests require manual setup")
    
    # This is a placeholder for future HTTP mocking
    # Real implementation would use httptest2::with_mock_api()
  })
  
  test_that("OpenAI chat can be mocked", {
    skip("Mock HTTP tests require manual setup")
    
    # This is a placeholder for future HTTP mocking
    # Real implementation would use httptest2::with_mock_api()
  })
}
