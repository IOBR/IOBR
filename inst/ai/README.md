# IOBR AI Assistant

## Overview

The IOBR AI Assistant is an integrated RAG (Retrieval-Augmented Generation) system that helps users explore and use the IOBR package for immune oncology biological research. It combines semantic search over package documentation with AI-powered responses to provide accurate, context-aware answers.

## Architecture

The AI assistant consists of three main components:

### 1. Provider-Agnostic Client (`R/ai_provider.R`)

Handles communication with various AI service providers:

- **Dummy Provider**: Local bag-of-words embeddings and simple response generation for offline testing
- **OpenAI Provider**: Integration with OpenAI's embedding and chat APIs
- **Hugging Face Provider**: Integration with Hugging Face Inference API
- **Custom Providers**: Extensible architecture for additional providers

Features:
- Exponential backoff retry logic for HTTP requests
- Batch processing for embeddings
- Configurable via runtime provider configuration (no hardcoded keys)

### 2. Indexing and Retrieval (`R/ai_index.R`)

Manages the document embedding index:

- **Document Collection**: Gathers content from README, vignettes, R source files, man pages, and examples
- **Text Chunking**: Intelligently splits documents into ~800 character chunks by paragraph
- **Index Creation**: Generates embeddings for all chunks and stores them in RDS format
- **Semantic Search**: Retrieves top-K most relevant chunks using cosine similarity

### 3. High-Level API (`R/ai.R`)

User-facing functions:

- `iobr_ai_init()`: Build the index from package documentation
- `iobr_ai_query()`: Ask questions and get AI-generated answers with sources
- `iobr_ai_configure_provider()`: Validate and normalize provider configuration
- `iobr_ai_list_index()`: View index metadata
- `iobr_ai_reset_index()`: Delete the index

### 4. Shiny Application (`inst/shiny/iobr-ai-app.R`)

Interactive web interface:

- Provider configuration (OpenAI, Hugging Face, Anthropic, Custom, Dummy)
- Index management (create, view, download, upload, reset)
- Query interface with adjustable parameters
- Response display with code blocks and source citations
- Sandbox execution button (disabled for safety)

## Quick Start

### Local Dummy Mode (No API Key Required)

Perfect for testing and development without external dependencies:

```r
# Load the package
library(IOBR)

# Initialize with dummy provider (default)
iobr_ai_init()

# Query the assistant
response <- iobr_ai_query("How do I calculate TME scores?")
cat(response$answer)

# View retrieved sources
print(response$retrieved)
```

### Using OpenAI

```r
# Configure provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Initialize index
iobr_ai_init(provider = provider)

# Query
response <- iobr_ai_query(
  "What deconvolution methods are available in IOBR?",
  provider = provider
)

cat(response$answer)
```

### Using Hugging Face

```r
# Configure provider
provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HUGGINGFACE_API_KEY")
))

# Initialize and query
iobr_ai_init(provider = provider)
response <- iobr_ai_query("Explain CIBERSORT", provider = provider)
```

### Custom Provider

```r
provider <- iobr_ai_configure_provider(list(
  name = "custom",
  api_key = "your-api-key",
  base_url = "https://your-api.com",
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model",
  headers = list(
    "X-Custom-Header" = "value"
  )
))

iobr_ai_init(provider = provider)
```

## Running the Shiny App

### Launch the App

```r
library(shiny)

# From package installation
runApp(system.file("shiny/iobr-ai-app.R", package = "IOBR"))

# Or directly from source
runApp("inst/shiny/iobr-ai-app.R")
```

### Using the App

1. **Select Provider**: Choose "Dummy (Offline)" for testing or configure a real provider
2. **Create Index**: Click "Create Index" to build the embedding database
3. **Query Tab**: Enter your question and click "Ask Assistant"
4. **View Results**: See the AI-generated answer, example code, and source citations

### Provider Configuration in Shiny

For **OpenAI**:
- API Key: Your OpenAI API key
- Leave other fields empty for defaults

For **Hugging Face**:
- API Key: Your Hugging Face token
- Leave other fields empty for defaults

For **Custom**:
- Fill in all fields according to your API requirements

## Security Considerations

### API Key Safety

⚠️ **NEVER commit API keys to version control**

Best practices:
```r
# Use environment variables
Sys.setenv(OPENAI_API_KEY = "your-key")
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)

# Or use .Renviron file
# Add to ~/.Renviron:
# OPENAI_API_KEY=your-key-here
```

### Code Execution

The assistant returns code examples but **does not execute them automatically**. This is a safety feature. Always:

1. Review code before running
2. Understand what the code does
3. Test in a safe environment
4. Check for package dependencies

The Shiny app's "Run in Sandbox" button is intentionally disabled. To enable it in the future, implement:
- Isolated R process execution
- Resource limits (memory, CPU, time)
- Restricted file system access
- Output sanitization

## Index Management

### Default Index Location

`inst/ai/iobr_embeddings.rds`

### Index Structure

```r
list(
  entries = list(
    list(
      id = 1,
      text = "chunk content",
      source = "README.md",
      type = "readme",
      embedding = c(0.1, 0.2, ...)
    ),
    ...
  ),
  vocab = c("word1", "word2", ...),  # For dummy provider
  created_at = Sys.time(),
  provider_meta = list(
    name = "dummy",
    model_embeddings = NULL
  )
)
```

### Index Operations

```r
# View index info
iobr_ai_list_index()

# Reset (delete) index
iobr_ai_reset_index()

# Create with custom path
iobr_ai_init(index_path = "path/to/custom/index.rds")

# Query with custom index
iobr_ai_query("question", index_path = "path/to/custom/index.rds")
```

## Response Structure

```r
response <- iobr_ai_query("question")

# Fields:
# - answer: AI-generated text response
# - code: List of extracted R code blocks (or NULL)
# - retrieved: List of top-K source chunks with:
#   - id: Chunk identifier
#   - source: Source file path
#   - type: Document type
#   - text: Chunk content
#   - score: Similarity score (0-1)
# - raw: Raw provider response
```

## Advanced Configuration

### Adjusting Search Parameters

```r
response <- iobr_ai_query(
  query = "How to use EPIC?",
  top_k = 10,           # Retrieve more sources
  max_tokens = 1000,    # Longer responses
  temperature = 0.3     # More creative (0-1)
)
```

### Custom Chunk Size

```r
# Smaller chunks (more precise, more entries)
iobr_ai_init(chunk_size = 500)

# Larger chunks (more context, fewer entries)
iobr_ai_init(chunk_size = 1200)
```

### Batch Size for Embeddings

Controlled in `send_embeddings()`:

```r
# Larger batches (faster but more memory)
send_embeddings(texts, provider, batch_size = 32)

# Smaller batches (slower but safer)
send_embeddings(texts, provider, batch_size = 8)
```

## Extending the Assistant

### Adding a New Provider

1. Edit `R/ai_provider.R`
2. Add internal functions:
   - `.send_embeddings_yourprovider()`
   - `.send_chat_yourprovider()`
3. Update `send_embeddings()` and `send_chat()` conditionals
4. Add default configuration in `iobr_ai_configure_provider()`

### Customizing Document Collection

Edit `collect_package_docs()` in `R/ai_index.R` to include additional sources:

```r
# Add custom directory
custom_dir <- file.path(pkg_path, "custom_docs")
if (dir.exists(custom_dir)) {
  # ... collect and add to docs list
}
```

### Improving Chunking Strategy

Modify `chunk_texts()` in `R/ai_index.R`:

```r
# Add overlap between chunks
# Use different splitting strategies
# Handle code blocks specially
```

## Troubleshooting

### "Index not found"

Run `iobr_ai_init()` first to create the index.

### "Provider requires api_key"

Set your API key:
```r
provider$api_key <- Sys.getenv("YOUR_API_KEY")
```

### "HTTP request failed"

- Check internet connection
- Verify API key is valid
- Check provider base_url is correct
- Review rate limits

### Shiny app can't find package

Ensure you're running from the package root or adjust the path finding logic in the app.

### Poor quality responses

- Increase `top_k` to provide more context
- Adjust `temperature` (lower = more deterministic)
- Use a more capable model
- Improve chunking strategy

## Development and Testing

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-ai.R")
```

### Mocking HTTP Requests

Tests use mocked responses to avoid requiring API keys:

```r
library(httptest2)

with_mock_api({
  # Test code here
})
```

## Future Enhancements

Potential improvements for future versions:

- [ ] Conversation history/context
- [ ] Multi-turn dialogue
- [ ] Fine-tuning on IOBR-specific terminology
- [ ] Caching of embeddings and responses
- [ ] Incremental index updates
- [ ] Support for images and plots
- [ ] Integration with R help system
- [ ] Sandbox code execution with safety measures
- [ ] Export conversations
- [ ] Analytics and usage tracking

## References

- IOBR Package: https://iobr.github.io/book/
- OpenAI API: https://platform.openai.com/docs/
- Hugging Face Inference API: https://huggingface.co/docs/api-inference/
- RAG Overview: https://arxiv.org/abs/2005.11401

## Support

For issues, questions, or contributions:
- GitHub Issues: https://github.com/IOBR/IOBR/issues
- Package Documentation: https://iobr.github.io/book/

---

**Version**: 1.0.0 (MVP)  
**Last Updated**: 2025-12-03
