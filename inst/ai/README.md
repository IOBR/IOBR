# IOBR AI Assistant

## Overview

The IOBR AI Assistant is a Retrieval-Augmented Generation (RAG) system that helps users understand and use the IOBR R package. It combines semantic search over package documentation with Large Language Model (LLM) capabilities to provide accurate, context-aware answers.

## Architecture

### Components

1. **Provider Client (`R/ai_provider.R`)**
   - Provider-agnostic interface for embeddings and chat completions
   - Support for multiple providers: dummy (local), OpenAI, Hugging Face, and extensible to custom providers
   - Built-in retry logic with exponential backoff
   - HTTP-based communication using `httr`

2. **Indexing and Retrieval (`R/ai_index.R`)**
   - Document collection from README, vignettes, R sources, man pages, and examples
   - Text chunking with configurable size
   - Embeddings generation and storage
   - Cosine similarity-based retrieval

3. **High-Level API (`R/ai.R`)**
   - `iobr_ai_init()`: Initialize and build the documentation index
   - `iobr_ai_query()`: Query the assistant with RAG
   - `iobr_ai_configure_provider()`: Configure provider settings
   - `iobr_ai_list_index()`: View index information
   - `iobr_ai_reset_index()`: Delete the index

4. **Shiny App (`inst/shiny/iobr-ai-app.R`)**
   - Web interface for provider configuration
   - Index management (create, view, download, upload)
   - Interactive query interface
   - Source attribution and provenance display

### Data Flow

```
User Query
    ↓
1. Query Embedding (via provider)
    ↓
2. Similarity Search (cosine similarity)
    ↓
3. Retrieve Top-K Chunks
    ↓
4. Build Context Prompt
    ↓
5. LLM Completion (via provider)
    ↓
6. Parse Response (extract answer & code)
    ↓
Return to User
```

### Index Structure

The embeddings index is stored as an RDS file (default: `inst/ai/iobr_embeddings.rds`) with the following structure:

```r
list(
  entries = list(
    list(
      id = <chunk_id>,
      text = <chunk_text>,
      source = <file_path>,
      embedding = <numeric_vector>
    ),
    ...
  ),
  vocab = <vocabulary for dummy provider>,
  created_at = <timestamp>,
  provider_meta = list(
    name = <provider_name>,
    model_embeddings = <model_name>
  )
)
```

## Quick Start

### Local Testing (No API Keys Required)

The dummy provider uses a simple bag-of-words approach for local testing without any external API calls.

```r
# Load the package
library(IOBR)

# Initialize with dummy provider (default)
index_info <- iobr_ai_init()

# Query the assistant
response <- iobr_ai_query("How do I analyze the tumor microenvironment?")

# View the answer
cat(response$answer)

# View code examples (if any)
if (!is.null(response$code)) {
  cat("\nCode:\n")
  cat(response$code[[1]])
}

# View retrieved sources
for (i in seq_along(response$retrieved)) {
  source <- response$retrieved[[i]]
  cat(sprintf("\nSource %d: %s (score: %.3f)\n", 
              i, source$source, source$score))
}
```

### Using OpenAI

```r
library(IOBR)

# Configure OpenAI provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  model_embeddings = "text-embedding-ada-002",
  model_chat = "gpt-3.5-turbo"
))

# Create index with OpenAI embeddings
index_info <- iobr_ai_init(provider = provider)

# Query with OpenAI
response <- iobr_ai_query(
  "What is CIBERSORT and how do I use it in IOBR?",
  provider = provider
)

cat(response$answer)
```

### Using Hugging Face

```r
library(IOBR)

# Configure Hugging Face provider
provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HF_TOKEN"),
  model_embeddings = "sentence-transformers/all-MiniLM-L6-v2",
  model_chat = "microsoft/DialoGPT-medium"
))

# Create index
index_info <- iobr_ai_init(provider = provider)

# Query
response <- iobr_ai_query(
  "How do I perform survival analysis?",
  provider = provider
)
```

### Custom Provider

For custom API endpoints (e.g., local LLM servers, Azure OpenAI):

```r
provider <- iobr_ai_configure_provider(list(
  name = "custom",
  api_key = "your-api-key",
  base_url = "https://your-api-endpoint.com/v1",
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model",
  headers = list(
    "X-Custom-Header" = "value"
  )
))
```

## Running the Shiny App

### Local Mode (No API Keys)

```r
library(shiny)

# Run the app
runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))
```

Or if developing:

```r
shiny::runApp("inst/shiny/iobr-ai-app.R")
```

The app will start with the dummy provider by default. You can:
1. Click "Create Index" to build the documentation index
2. Enter a query and click "Submit Query"
3. View the answer and retrieved sources

### With Real Providers

1. Select provider (openai, huggingface, etc.) from dropdown
2. Enter your API key
3. (Optional) Customize models and endpoints
4. Click "Configure Provider"
5. Click "Create Index"
6. Submit queries

## Configuration Examples

### Environment Variables (Recommended)

Store API keys securely:

```bash
# In ~/.Renviron or .env
OPENAI_API_KEY=sk-...
HF_TOKEN=hf_...
```

Then in R:

```r
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)
```

### Provider-Specific Settings

**OpenAI**:
- Default base_url: `https://api.openai.com/v1`
- Default embeddings model: `text-embedding-ada-002`
- Default chat model: `gpt-3.5-turbo`

**Hugging Face**:
- Default base_url: `https://api-inference.huggingface.co`
- Default embeddings model: `sentence-transformers/all-MiniLM-L6-v2`
- Default chat model: `microsoft/DialoGPT-medium`

**Azure OpenAI**:
```r
provider <- list(
  name = "custom",
  api_key = Sys.getenv("AZURE_OPENAI_KEY"),
  base_url = "https://your-resource.openai.azure.com/openai/deployments/your-deployment",
  headers = list(
    "api-version" = "2023-05-15"
  )
)
```

## Security Best Practices

⚠️ **IMPORTANT SECURITY NOTES**:

1. **Never commit API keys to version control**
   - Use environment variables or secure configuration management
   - Add API key files to `.gitignore`

2. **Dummy provider for testing**
   - Use the dummy provider for local development and testing
   - No API keys needed, works completely offline

3. **Rate limiting and costs**
   - Be aware of API rate limits and costs
   - The dummy provider has no limits or costs

4. **Data privacy**
   - Documentation is sent to external APIs when using real providers
   - Use dummy provider if documentation contains sensitive information

5. **Sandbox execution**
   - Code examples are NOT automatically executed
   - Future versions may include sandboxed execution (currently disabled)

## Advanced Usage

### Custom Chunking

```r
# Smaller chunks for more granular retrieval
index_info <- iobr_ai_init(chunk_size = 400)

# Larger chunks for more context
index_info <- iobr_ai_init(chunk_size = 1200)
```

### Adjusting Retrieval

```r
# Retrieve more context chunks
response <- iobr_ai_query(
  "Complex question...",
  top_k = 10  # Default is 5
)

# Adjust temperature for more creative responses
response <- iobr_ai_query(
  "Question...",
  temperature = 0.7  # Default is 0.2 (more deterministic)
)
```

### Managing Multiple Indices

```r
# Create index at custom location
index_info <- iobr_ai_init(index_path = "/path/to/my_index.rds")

# Query from custom location
response <- iobr_ai_query(
  "Question...",
  index_path = "/path/to/my_index.rds"
)

# List index info
iobr_ai_list_index(index_path = "/path/to/my_index.rds")
```

### Programmatic Testing

```r
# Unit testing with dummy provider
test_that("AI query returns expected structure", {
  provider <- list(name = "dummy")
  index_info <- iobr_ai_init(provider = provider)
  
  response <- iobr_ai_query("test query", provider = provider)
  
  expect_true("answer" %in% names(response))
  expect_true("retrieved" %in% names(response))
  expect_true(length(response$retrieved) > 0)
})
```

## Troubleshooting

### "Index not found"
- Run `iobr_ai_init()` to create the index first
- Check the index path with `iobr_ai_list_index()`

### "Provider requires API key"
- Ensure API key is set in provider configuration
- Use dummy provider for testing without API keys

### "API error"
- Check API key validity
- Verify network connectivity
- Check provider base_url is correct
- Review rate limits

### Empty or poor responses
- Increase `top_k` to retrieve more context
- Adjust `chunk_size` when creating index
- Try different temperature settings
- Ensure index was created successfully

## Future Enhancements

Potential improvements for future versions:

1. **Sandboxed code execution**: Run generated code safely in isolated environment
2. **More providers**: Anthropic Claude, Cohere, local models via Ollama
3. **Hybrid search**: Combine semantic search with keyword search
4. **Conversation history**: Multi-turn conversations with context
5. **Fine-tuned models**: Custom models trained on IOBR documentation
6. **Vector databases**: Integration with Chroma, Pinecone, etc.
7. **Incremental updates**: Update index without full rebuild

## Contributing

To add support for a new provider:

1. Add embedding function in `R/ai_provider.R`:
   ```r
   send_embeddings_myprovider <- function(texts, provider, batch_size, max_retries) {
     # Implementation
   }
   ```

2. Add chat function:
   ```r
   send_chat_myprovider <- function(system_msg, user_msg, provider, ...) {
     # Implementation
   }
   ```

3. Update `send_embeddings()` and `send_chat()` dispatch logic

4. Add default configuration in `iobr_ai_configure_provider()`

5. Update documentation and tests

## Support

For issues, questions, or contributions:
- GitHub Issues: https://github.com/IOBR/IOBR/issues
- Documentation: https://iobr.github.io/book/

## License

This feature is part of IOBR and follows the same GPL-3 license.
