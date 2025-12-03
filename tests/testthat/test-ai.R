# Tests for IOBR AI Assistant

test_that("provider configuration validation works", {
  # Valid OpenAI config
  config <- list(
    name = "openai",
    api_key = "test-key"
  )
  validated <- iobr_ai_configure_provider(config)
  expect_equal(validated$name, "openai")
  expect_equal(validated$api_key, "test-key")
  expect_equal(validated$base_url, "https://api.openai.com/v1")
  expect_equal(validated$model_embeddings, "text-embedding-ada-002")
  
  # Valid Hugging Face config
  config <- list(
    name = "huggingface",
    api_key = "test-key"
  )
  validated <- iobr_ai_configure_provider(config)
  expect_equal(validated$name, "huggingface")
  
  # Invalid: missing name
  expect_error(
    iobr_ai_configure_provider(list(api_key = "test")),
    "Provider name is required"
  )
  
  # Invalid: missing api_key
  expect_error(
    iobr_ai_configure_provider(list(name = "openai")),
    "API key is required"
  )
  
  # Invalid: unsupported provider
  expect_error(
    iobr_ai_configure_provider(list(name = "invalid", api_key = "test")),
    "Unsupported provider"
  )
  
  # Custom provider needs base_url
  expect_error(
    iobr_ai_configure_provider(list(name = "custom", api_key = "test")),
    "base_url is required"
  )
})

test_that("chunking works correctly", {
  # Short text - no chunking needed
  short_text <- "This is a short text."
  chunks <- chunk_text(short_text, chunk_size = 100)
  expect_length(chunks, 1)
  expect_equal(chunks[1], short_text)
  
  # Long text - needs chunking
  long_text <- paste(rep("This is a paragraph. ", 100), collapse = "")
  chunks <- chunk_text(long_text, chunk_size = 500)
  expect_gt(length(chunks), 1)
  
  # Each chunk should be <= chunk_size (with some tolerance for overlap)
  for (chunk in chunks) {
    expect_lte(nchar(chunk), 700)  # chunk_size + overlap
  }
})

test_that("cosine similarity works", {
  v1 <- c(1, 0, 0)
  v2 <- c(1, 0, 0)
  expect_equal(cosine_similarity(v1, v2), 1.0)
  
  v1 <- c(1, 0, 0)
  v2 <- c(0, 1, 0)
  expect_equal(cosine_similarity(v1, v2), 0.0)
  
  v1 <- c(1, 1, 0)
  v2 <- c(1, 1, 0)
  expect_equal(cosine_similarity(v1, v2), 1.0)
  
  # Different lengths should error
  expect_error(
    cosine_similarity(c(1, 2), c(1, 2, 3)),
    "same length"
  )
})

test_that("Rd to text conversion works", {
  # Create a temporary Rd file
  rd_content <- c(
    "\\name{test_func}",
    "\\alias{test_func}",
    "\\title{Test Function}",
    "\\description{This is a test function.}",
    "\\usage{test_func(x)}",
    "\\arguments{",
    "  \\item{x}{Input parameter}",
    "}",
    "\\value{Returns a value.}",
    "\\examples{",
    "test_func(1)",
    "}"
  )
  
  temp_rd <- tempfile(fileext = ".Rd")
  writeLines(rd_content, temp_rd)
  
  text <- rd_to_text(temp_rd)
  
  expect_true(grepl("Function: test_func", text))
  expect_true(grepl("Test Function", text))
  expect_true(grepl("Description:", text))
  
  unlink(temp_rd)
})

test_that("mock embedding request works", {
  skip_if_not_installed("mockery")
  
  # Skip this test as it requires complex mocking that may not work in all environments
  skip("Complex mocking test - manual verification required")
  
  # Mock the httr2 functions would go here
  # In a real test environment with proper httr2 mocking setup
})

test_that("mock chat request works", {
  skip_if_not_installed("mockery")
  
  # Skip this test as it requires complex mocking that may not work in all environments
  skip("Complex mocking test - manual verification required")
  
  # Mock the httr2 functions would go here
  # In a real test environment with proper httr2 mocking setup
})

test_that("index save and load works", {
  # Create a temporary index
  index <- list(
    list(
      id = "test_1",
      text = "Test text 1",
      embedding = c(0.1, 0.2, 0.3),
      source = "test.R",
      type = "code"
    ),
    list(
      id = "test_2",
      text = "Test text 2",
      embedding = c(0.4, 0.5, 0.6),
      source = "test.R",
      type = "code"
    )
  )
  
  temp_file <- tempfile(fileext = ".rds")
  
  # Save
  save_index(index, temp_file)
  expect_true(file.exists(temp_file))
  
  # Load
  loaded <- load_index(temp_file)
  expect_equal(length(loaded), 2)
  expect_equal(loaded[[1]]$text, "Test text 1")
  expect_equal(loaded[[2]]$embedding, c(0.4, 0.5, 0.6))
  
  unlink(temp_file)
})

test_that("collect_package_content finds files", {
  skip_if_not(dir.exists("../../R"), "Not in package root")
  
  # This test assumes we're in the tests directory
  pkg_path <- "../.."
  
  content <- collect_package_content(pkg_path)
  
  expect_gt(length(content), 0)
  
  # Check that we have different types
  types <- sapply(content, function(x) x$type)
  expect_true("readme" %in% types || "code" %in% types)
  
  # Each item should have required fields
  for (item in content) {
    expect_true(!is.null(item$text))
    expect_true(!is.null(item$source))
    expect_true(!is.null(item$type))
  }
})

test_that("iobr_ai_list_index handles missing index", {
  # Try with a non-existent path
  expect_error(
    iobr_ai_list_index(index_path = "/nonexistent/path/index.rds"),
    "Index file not found"
  )
})

test_that("iobr_ai_reset_index works", {
  temp_file <- tempfile(fileext = ".rds")
  
  # Create a dummy index
  saveRDS(list(), temp_file)
  expect_true(file.exists(temp_file))
  
  # Reset
  result <- iobr_ai_reset_index(index_path = temp_file)
  expect_true(result)
  expect_false(file.exists(temp_file))
  
  # Reset non-existent file should not error
  expect_true(iobr_ai_reset_index(index_path = temp_file))
})
