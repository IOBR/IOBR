# IOBR AI Assistant - Testing and Usage Guide

## Overview

This PR adds a complete MVP AI assistant to the IOBR package using Retrieval-Augmented Generation (RAG). The assistant helps users explore package documentation and get AI-powered answers about TME analysis methods.

## What's Included

### Core Implementation Files

1. **R/ai_provider.R** - Provider-agnostic client
   - Dummy provider (offline, bag-of-words embeddings)
   - OpenAI provider adapter
   - Hugging Face provider adapter
   - HTTP retry with exponential backoff
   - Configurable via runtime provider settings

2. **R/ai_index.R** - Indexing and retrieval
   - `collect_package_docs()` - Gathers README, vignettes, R files, man pages
   - `chunk_texts()` - Smart paragraph-based chunking (~800 chars)
   - `create_index()` - Generates embeddings and stores index
   - `search_index()` - Semantic search using cosine similarity
   - `cosine_sim()` - Vector similarity helper

3. **R/ai.R** - High-level API
   - `iobr_ai_init()` - Build index from package docs
   - `iobr_ai_query()` - Query assistant with RAG
   - `iobr_ai_configure_provider()` - Validate provider config
   - `iobr_ai_list_index()` - View index metadata
   - `iobr_ai_reset_index()` - Delete index

4. **inst/shiny/iobr-ai-app.R** - Interactive Shiny app
   - Provider configuration UI (OpenAI, Hugging Face, Anthropic, Custom, Dummy)
   - Index management (create, view, download, upload, reset)
   - Query interface with adjustable parameters
   - Response display with code blocks and source citations
   - Disabled "Run in Sandbox" button (safety feature)

5. **inst/ai/README.md** - Comprehensive documentation
   - Architecture overview
   - Quick start guides for all providers
   - Security best practices
   - Advanced configuration
   - Troubleshooting

6. **tests/testthat/test-ai.R** - Unit tests
   - Tests for all core functions
   - Uses dummy provider (no API keys required)
   - Validates index structure
   - Tests search and query functionality

### Updated Files

- **DESCRIPTION** - Added dependencies: shiny, httr, jsonlite, httptest2
- **NAMESPACE** - Added exports for all public functions

## Testing Instructions

### 1. Local Testing with Dummy Provider (No API Keys Required)

This is the recommended way to test the implementation without any external dependencies:

```r
# Install/load the package
library(IOBR)

# Initialize with dummy provider (default)
metadata <- iobr_ai_init()
# This will collect docs, create chunks, and build an index using local embeddings

# Query the assistant
response <- iobr_ai_query("How do I calculate TME scores using IOBR?")

# View the answer
cat(response$answer)

# View retrieved sources
print(response$retrieved)

# View extracted code (if any)
if (!is.null(response$code)) {
  cat("\nCode examples:\n")
  for (code in response$code) {
    cat(code, "\n\n")
  }
}

# View index information
iobr_ai_list_index()
```

### 2. Testing with the Shiny App (Dummy Mode)

```r
library(shiny)

# Launch the app
runApp(system.file("shiny/iobr-ai-app.R", package = "IOBR"))

# Or if testing from source:
runApp("inst/shiny/iobr-ai-app.R")
```

**Steps in the app:**
1. Keep "Dummy (Offline)" selected as the provider
2. Click "Create Index" button
3. Wait for index creation to complete
4. Go to "Query Assistant" tab
5. Enter a question like: "What deconvolution methods are available?"
6. Click "Ask Assistant"
7. Review the answer, code blocks, and sources

### 3. Testing with OpenAI (Optional - Requires API Key)

If you want to test with a real AI provider:

```r
# Set up environment variable (recommended)
Sys.setenv(OPENAI_API_KEY = "your-api-key-here")

# Or add to ~/.Renviron:
# OPENAI_API_KEY=your-api-key-here

# Configure provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Initialize (will use OpenAI embeddings)
iobr_ai_init(provider = provider)

# Query (will use OpenAI chat)
response <- iobr_ai_query(
  "Explain the CIBERSORT deconvolution method",
  provider = provider
)

cat(response$answer)
```

**In Shiny app:**
1. Select "OpenAI" from provider dropdown
2. Paste your API key
3. Leave other fields empty for defaults
4. Create index and query as before

### 4. Testing with Hugging Face (Optional - Requires Token)

```r
# Set up token
Sys.setenv(HUGGINGFACE_API_KEY = "your-hf-token-here")

provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HUGGINGFACE_API_KEY")
))

iobr_ai_init(provider = provider)
response <- iobr_ai_query("How to use EPIC?", provider = provider)
```

### 5. Running Unit Tests

```r
# Install test dependencies if needed
install.packages(c("testthat", "httptest2"))

# Run all tests
devtools::test()

# Or run just the AI tests
testthat::test_file("tests/testthat/test-ai.R")
```

**Expected test results:**
- All tests should pass using the dummy provider
- No network calls are made during tests
- Tests validate index structure, search, and query functionality

## Security Notes

⚠️ **IMPORTANT**: Never commit API keys to version control!

**Best practices:**
- Use environment variables: `Sys.getenv("OPENAI_API_KEY")`
- Store in `~/.Renviron` file (not tracked by git)
- Use secure secret management in production

**Code execution safety:**
- The assistant returns code but does NOT execute it automatically
- Always review code before running
- The Shiny app's "Run in Sandbox" button is intentionally disabled

## Expected Behavior

### Dummy Provider
- **Embeddings**: Simple bag-of-words vectors normalized by frequency
- **Chat**: Returns templated responses for offline testing
- **No network calls**: Everything runs locally
- **Purpose**: Enable testing and development without API keys

### Real Providers (OpenAI/Hugging Face)
- **Embeddings**: High-quality semantic vectors from trained models
- **Chat**: Actual AI-generated responses
- **Network calls**: Uses provider APIs with retry/backoff
- **Better quality**: More accurate retrieval and natural responses

## Index Management

### Default Location
`inst/ai/iobr_embeddings.rds`

### Index Operations
```r
# View info
iobr_ai_list_index()

# Reset (delete)
iobr_ai_reset_index()

# Custom path
iobr_ai_init(index_path = "custom/path/index.rds")
iobr_ai_query("question", index_path = "custom/path/index.rds")
```

### Index Portability
- Download from Shiny app
- Share between environments
- Upload to new installations
- Version control (with .gitignore for generated indexes)

## Troubleshooting

### "Index not found"
**Solution**: Run `iobr_ai_init()` first to create the index.

### "Package 'httr' is required"
**Solution**: Install dependencies:
```r
install.packages(c("httr", "jsonlite", "shiny"))
```

### Shiny app can't find package
**Solution**: Ensure you're running from package root, or adjust path finding:
```r
# If testing from different directory
setwd("/path/to/IOBR")
runApp("inst/shiny/iobr-ai-app.R")
```

### Poor quality responses (dummy mode)
**Expected**: Dummy provider gives generic responses for testing.
**Solution**: Use a real provider (OpenAI/Hugging Face) for actual assistance.

### OpenAI/HF API errors
**Check**:
- API key is valid
- Internet connection works
- Provider base_url is correct
- No rate limiting issues

## Example Queries to Test

1. "What is IOBR used for?"
2. "How do I perform TME deconvolution?"
3. "What's the difference between CIBERSORT and EPIC?"
4. "Show me how to calculate immune signatures"
5. "How do I visualize TME features?"
6. "What are the available deconvolution methods?"
7. "Explain the IPS calculation method"
8. "How to handle batch effects in TME data?"

## Response Structure

Every query returns:
```r
list(
  answer = "AI-generated text response",
  code = list("R code block 1", "R code block 2", ...) or NULL,
  retrieved = list(
    list(id, source, type, text, score),
    ...
  ),
  raw = list(...)  # Raw provider response
)
```

## Performance Notes

- **Index creation**: Takes ~30-60 seconds for IOBR package (depends on size)
- **Query time**: 
  - Dummy: < 1 second (local)
  - OpenAI: 2-5 seconds (API latency)
  - Hugging Face: 3-10 seconds (API + model loading)
- **Index size**: ~2-10 MB depending on package documentation

## Future Enhancements

Potential improvements for future versions:
- Conversation history/multi-turn dialogue
- Incremental index updates
- Support for images and plots
- Integration with R help system
- Sandbox code execution with safety measures
- Fine-tuning on IOBR terminology
- Caching of responses

## Questions or Issues?

- Check `inst/ai/README.md` for detailed documentation
- Open GitHub issues for bugs or feature requests
- Review test file for implementation examples

---

**Testing Checklist:**

- [ ] Install package with dependencies
- [ ] Run unit tests (all pass)
- [ ] Test dummy provider initialization
- [ ] Test dummy provider query
- [ ] View index info
- [ ] Launch Shiny app (dummy mode)
- [ ] Create index in Shiny
- [ ] Query in Shiny
- [ ] Download/upload index
- [ ] (Optional) Test with OpenAI
- [ ] (Optional) Test with Hugging Face

**All core functionality should work with dummy provider (no API keys needed).**
