# IOBR AI Assistant - Testing Instructions

This document provides instructions for testing the IOBR AI Assistant MVP implementation.

## Overview

The AI Assistant adds semantic search and question-answering capabilities to the IOBR package. It supports:
- **Dummy provider**: Local mode for testing (no API key required)
- **OpenAI provider**: GPT models and embeddings
- **Hugging Face provider**: Open-source models

## Prerequisites

### For Local Testing (Dummy Provider)
No additional setup required! The dummy provider works offline with no API keys.

### For Real Provider Testing (Optional)
- **OpenAI**: API key from https://platform.openai.com/api-keys
- **Hugging Face**: API key from https://huggingface.co/settings/tokens

## Installation

Since this is a feature branch, install from the branch:

```r
# Install from GitHub (if package is public)
# remotes::install_github("IOBR/IOBR", ref = "ai")

# Or install locally from source
# In the package directory:
devtools::install()

# Load the package
library(IOBR)
```

## Test 1: Basic Functionality with Dummy Provider

This test verifies the AI assistant works without any API keys.

```r
library(IOBR)

# Initialize the AI assistant (creates embeddings index)
result <- iobr_ai_init()
# Expected: Messages showing document collection, chunking, and index creation
# Should complete without errors

# View index information
info <- iobr_ai_list_index()
# Expected: Shows number of entries, creation time, provider info

# Query the assistant
response <- iobr_ai_query("How do I use CIBERSORT for deconvolution?")

# View the answer
cat(response$answer)
# Expected: Text response about CIBERSORT

# View retrieved sources
print(response$retrieved)
# Expected: Data frame with id, source, text, score columns

# Check if code was extracted
if (!is.null(response$code)) {
  cat("\nExtracted code:\n")
  cat(response$code)
}
```

**Expected Results:**
- Index created successfully in `inst/ai/iobr_embeddings.rds`
- Query returns structured response with answer, retrieved chunks, and scores
- All functions execute without errors

## Test 2: Programmatic API

Test the high-level API functions.

```r
library(IOBR)

# Configure provider explicitly
provider <- iobr_ai_configure_provider(list(name = "dummy"))
print(provider)
# Expected: List with name, model_embeddings, model_chat

# Initialize with custom settings
result <- iobr_ai_init(
  pkg_path = ".",
  provider = provider,
  chunk_size = 1000
)

# Query with different parameters
response1 <- iobr_ai_query(
  "What signature calculation methods are available?",
  provider = provider,
  top_k = 3,
  temperature = 0.1
)

response2 <- iobr_ai_query(
  "Explain TME deconvolution",
  provider = provider,
  top_k = 10,
  temperature = 0.5
)

# Compare results
cat("Response 1 retrieved", nrow(response1$retrieved), "chunks\n")
cat("Response 2 retrieved", nrow(response2$retrieved), "chunks\n")

# List index
iobr_ai_list_index()

# Reset index
iobr_ai_reset_index()
# Expected: Deletes inst/ai/iobr_embeddings.rds

# Verify deletion
if (!file.exists("inst/ai/iobr_embeddings.rds")) {
  cat("Index successfully deleted\n")
}
```

**Expected Results:**
- Provider configuration works correctly
- Different parameters affect results (top_k changes number of retrieved chunks)
- Index deletion removes the file

## Test 3: Shiny App with Dummy Provider

Test the web interface.

```r
library(IOBR)
library(shiny)

# Launch the app
runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))

# Or from source directory:
# runApp("inst/shiny/iobr-ai-app.R")
```

**Manual Testing Steps:**

1. **Provider Configuration Tab:**
   - Ensure "dummy" is selected (default)
   - No API key required

2. **Create Index:**
   - Leave default settings (pkg_path = ".", chunk_size = 800)
   - Click "Create/Update Index"
   - Wait for completion message
   - Click "Show Index Info" to verify

3. **Query Tab:**
   - Enter a question: "How do I calculate signature scores?"
   - Set top_k = 5, temperature = 0.2
   - Click "Submit Query"
   - Verify response appears
   - Check retrieved sources are shown with scores
   - Expand source details to view text snippets

4. **Test Multiple Queries:**
   - "What deconvolution methods are available?"
   - "How to use TIMER?"
   - "Explain tumor microenvironment analysis"

5. **Index Management:**
   - Click "Download Index" to save the RDS file
   - Click "Reset Index" to delete
   - Upload the downloaded index to restore

6. **Verify Security Features:**
   - If code is shown, verify "Run in Sandbox" button is disabled
   - Check warning message about not executing code

**Expected Results:**
- App loads without errors
- Index creation completes successfully
- Queries return relevant responses
- Retrieved sources display with scores
- Download/upload functionality works

## Test 4: OpenAI Provider (Optional - Requires API Key)

**⚠️ Warning: This will use API credits**

```r
library(IOBR)

# Option 1: Use environment variable (recommended)
# Add to ~/.Renviron:
# OPENAI_API_KEY=sk-your-key-here

# Option 2: Set in R session
Sys.setenv(OPENAI_API_KEY = "sk-your-key-here")

# Configure OpenAI provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Create index with OpenAI embeddings
result <- iobr_ai_init(
  provider = provider,
  index_path = "inst/ai/iobr_embeddings_openai.rds"
)

# Query with OpenAI
response <- iobr_ai_query(
  "How do I perform TME deconvolution with multiple methods?",
  index_path = "inst/ai/iobr_embeddings_openai.rds",
  provider = provider
)

cat(response$answer)
print(response$retrieved)
```

**Expected Results:**
- Index creation makes HTTP requests to OpenAI API
- Embeddings are higher dimensional (1536 for ada-002)
- Chat responses are more sophisticated than dummy mode
- No errors related to authentication or rate limits

**Troubleshooting:**
- Invalid API key: Check the key is correct and active
- Rate limit errors: Wait and retry, or check API tier limits
- Network errors: Verify internet connectivity

## Test 5: Shiny App with OpenAI (Optional)

```r
runApp("inst/shiny/iobr-ai-app.R")
```

1. Select "openai" from provider dropdown
2. Paste API key in "API Key" or "API Key (secure)" field
3. Optionally customize models:
   - Embeddings: text-embedding-ada-002 (default)
   - Chat: gpt-3.5-turbo or gpt-4
4. Create index and query as in Test 3

**Expected Results:**
- Real AI responses with better quality than dummy mode
- Index may be larger due to higher-dimensional embeddings

## Test 6: Unit Tests

Run the automated test suite:

```r
# Install test dependencies if needed
install.packages(c("testthat", "mockery"))

# Run all tests
testthat::test_file("tests/testthat/test-ai.R")

# Or run specific tests
testthat::test_that("dummy provider embeddings work", {
  # Test code...
})
```

**Expected Results:**
- All tests pass
- Specific test coverage:
  - ✓ Dummy provider embeddings
  - ✓ Dummy provider chat
  - ✓ Cosine similarity calculation
  - ✓ Text chunking
  - ✓ Provider configuration
  - ✓ Index structure validation
  - ✓ Query response structure
  - ✓ Search with top-k
  - ✓ Code extraction

## Test 7: Edge Cases and Error Handling

Test error conditions:

```r
library(IOBR)

# Test 1: Query without index
tryCatch({
  iobr_ai_query("test query")
}, error = function(e) {
  cat("Expected error:", e$message, "\n")
})
# Expected: "Index not found" error

# Test 2: Invalid provider
tryCatch({
  provider <- list(name = "invalid_provider")
  send_embeddings("test", provider)
}, error = function(e) {
  cat("Expected error:", e$message, "\n")
})
# Expected: "Unknown provider" error

# Test 3: OpenAI without API key
tryCatch({
  provider <- iobr_ai_configure_provider(list(name = "openai"))
}, error = function(e) {
  cat("Expected error:", e$message, "\n")
})
# Expected: "requires api_key" error

# Test 4: Empty query
response <- iobr_ai_query("")
# Should handle gracefully

# Test 5: Very large chunk size
result <- iobr_ai_init(chunk_size = 10000)
# Should work but create fewer chunks
```

**Expected Results:**
- Appropriate error messages for invalid inputs
- No crashes or undefined behavior
- Graceful handling of edge cases

## Test 8: Performance and Scalability

Test with different package sizes:

```r
library(IOBR)
library(tictoc)

# Time index creation
tic()
result <- iobr_ai_init()
time_create <- toc()

# Time queries
tic()
response1 <- iobr_ai_query("test query 1")
toc()

tic()
response2 <- iobr_ai_query("test query 2")
toc()

# Check index size
info <- iobr_ai_list_index()
cat("Index size:", info$file_size_mb, "MB\n")
cat("Number of chunks:", info$num_entries, "\n")
```

**Expected Results:**
- Index creation takes < 5 minutes for typical package
- Queries complete in < 5 seconds with dummy provider
- Index size is reasonable (< 10 MB for most packages)

## Test 9: Documentation and Help

Verify documentation is accessible:

```r
# Function help
?iobr_ai_init
?iobr_ai_query
?send_embeddings
?send_chat
?collect_package_docs
?chunk_texts
?create_index
?search_index

# View AI README
file.show(system.file("ai", "README.md", package = "IOBR"))
```

**Expected Results:**
- All functions have Roxygen documentation
- Help pages display correctly
- README provides comprehensive guidance

## Test 10: Security Verification

Ensure security features are working:

```r
library(IOBR)

# Initialize and query
iobr_ai_init()
response <- iobr_ai_query("Show me R code examples")

# Verify code is NOT executed
if (!is.null(response$code)) {
  cat("Code extracted but NOT executed:\n")
  cat(response$code, "\n")
  cat("\nThis is correct behavior for security.\n")
}

# In Shiny app: verify "Run in Sandbox" button is disabled
```

**Expected Results:**
- Code is extracted and displayed but never executed
- No automatic evaluation of user input
- Clear warnings about reviewing code before running

## Reporting Issues

If you encounter issues, please report with:

1. **Issue description**: What you expected vs. what happened
2. **Reproducible example**: Minimal code to reproduce the issue
3. **Environment**: R version, OS, package versions
4. **Provider**: Which provider (dummy, openai, huggingface)
5. **Error messages**: Full error output if applicable

Example issue report:

```
Title: Index creation fails with large packages

Description:
iobr_ai_init() crashes when package has > 100 documentation files.

Reproducible example:
library(IOBR)
result <- iobr_ai_init(pkg_path = "/path/to/large/package")

Environment:
- R version 4.3.0
- Ubuntu 22.04
- IOBR version: 0.99.99 (ai branch)

Error:
Error in chunk_texts(...): memory allocation error

Provider: dummy
```

## Success Criteria

The implementation is considered successful if:

1. ✅ Dummy provider works without any configuration
2. ✅ Index can be created from package documentation
3. ✅ Queries return relevant results with source attribution
4. ✅ Shiny app loads and all features work
5. ✅ OpenAI provider works with valid API key
6. ✅ Unit tests pass
7. ✅ No code execution vulnerabilities
8. ✅ Clear documentation and error messages
9. ✅ Reasonable performance (index < 5 min, queries < 5 sec)
10. ✅ No API keys committed to repository

## Additional Notes

### Development Mode
When developing or testing, you can work from source:
```r
devtools::load_all()
```

### Debugging
Enable verbose output:
```r
options(verbose = TRUE)
result <- iobr_ai_init()
```

### Clean Up
Remove all AI-related files:
```r
iobr_ai_reset_index()
unlink("inst/ai", recursive = TRUE)
```

## Contact

For questions or issues:
- GitHub: https://github.com/IOBR/IOBR/issues
- Documentation: https://iobr.github.io/book/
