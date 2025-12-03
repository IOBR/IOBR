# IOBR AI Assistant

An MVP AI assistant for the IOBR R package that enables semantic search and question-answering over package documentation and code.

## Architecture

The IOBR AI assistant consists of three main components:

### 1. Provider Layer (`R/ai_provider.R`)
- **Provider-agnostic client** supporting multiple backends:
  - **Dummy**: Local bag-of-words embeddings and template-based chat (no API key required)
  - **OpenAI**: Integration with OpenAI's embedding and chat APIs
  - **Hugging Face**: Integration with Hugging Face inference API
- **HTTP with retry/backoff**: Automatic retries with exponential backoff for resilient API calls
- **Runtime configuration**: All API requests driven by provider config (api_key, base_url, models, headers)

### 2. Indexing & Retrieval Layer (`R/ai_index.R`)
- **Document collection**: Gathers README, vignettes, R source files, man pages, examples
- **Chunking**: Splits documents into ~800 character chunks by paragraph boundaries
- **Embeddings index**: Generates and stores embeddings for all chunks
- **Semantic search**: Finds relevant chunks using cosine similarity
- **Index storage**: Saves to `inst/ai/iobr_embeddings.rds` by default

### 3. High-Level API (`R/ai.R`)
- **`iobr_ai_init()`**: One-command index creation from package docs
- **`iobr_ai_query()`**: Retrieve context and generate AI responses
- **`iobr_ai_configure_provider()`**: Validate and normalize provider settings
- **`iobr_ai_list_index()`**: View index metadata
- **`iobr_ai_reset_index()`**: Delete index

### 4. Shiny App (`inst/shiny/iobr-ai-app.R`)
- Web interface for all AI assistant features
- Provider configuration UI
- Index management (create, view, download, upload)
- Interactive query interface with source display
- Disabled sandbox execution button (security-first approach)

## Quick Start

### Local Dummy Mode (No API Key Required)

This is the recommended way to test the AI assistant without any external dependencies:

```r
# Load IOBR
library(IOBR)

# Initialize the AI assistant with dummy provider
# This creates a local bag-of-words index
result <- iobr_ai_init()
# Output: Using dummy provider (local mode, no API key required)
#         Collecting package documentation...
#         Collected 85 document(s)
#         Chunking documents...
#         Created 523 chunk(s)
#         Generating embeddings for 523 chunks...
#         Index saved to: inst/ai/iobr_embeddings.rds

# Ask a question
response <- iobr_ai_query("How do I use CIBERSORT for deconvolution?")

# View the answer
cat(response$answer)

# View retrieved sources with scores
print(response$retrieved)

# If code examples are provided, they are extracted
if (!is.null(response$code)) {
  cat("\nCode examples:\n")
  cat(response$code)
}
```

### Using Real AI Providers

#### OpenAI

```r
library(IOBR)

# Configure OpenAI provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY"),  # Store in .Renviron
  model_embeddings = "text-embedding-ada-002",  # optional
  model_chat = "gpt-3.5-turbo"  # optional
))

# Create index with OpenAI embeddings
result <- iobr_ai_init(provider = provider)

# Query with OpenAI
response <- iobr_ai_query(
  "What signature calculation methods are available?",
  provider = provider
)

cat(response$answer)
```

#### Hugging Face

```r
library(IOBR)

# Configure Hugging Face provider
provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HF_API_KEY"),  # Store in .Renviron
  model_embeddings = "sentence-transformers/all-MiniLM-L6-v2",  # optional
  model_chat = "mistralai/Mistral-7B-Instruct-v0.1"  # optional
))

# Create index and query
result <- iobr_ai_init(provider = provider)
response <- iobr_ai_query("Explain TME deconvolution methods", provider = provider)
```

#### Custom Provider

```r
# For custom API endpoints
provider <- iobr_ai_configure_provider(list(
  name = "custom",
  api_key = "your-key",
  base_url = "https://your-api.com/v1",
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model",
  headers = list(
    "X-Custom-Header" = "value"
  )
))

# Note: Custom providers may require adapter implementation
# in R/ai_provider.R to handle API-specific request/response formats
```

### Running the Shiny App

```r
library(IOBR)

# Launch the Shiny app
shiny::runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))

# Or if running from package source
shiny::runApp("inst/shiny/iobr-ai-app.R")
```

The Shiny app provides:
- Provider configuration UI (select from dummy, OpenAI, Hugging Face, Anthropic, custom)
- Index creation and management
- Interactive query interface
- Display of retrieved sources with similarity scores
- Code extraction (not executed for security)
- Index download/upload functionality

## Security Best Practices

### API Key Management

**Never commit API keys to version control.** Use one of these approaches:

1. **Environment variables** (recommended):
   ```r
   # Add to ~/.Renviron
   OPENAI_API_KEY=sk-...
   HF_API_KEY=hf_...
   
   # Use in R
   provider <- list(
     name = "openai",
     api_key = Sys.getenv("OPENAI_API_KEY")
   )
   ```

2. **Secure prompts** in interactive sessions:
   ```r
   provider <- list(
     name = "openai",
     api_key = readline("Enter OpenAI API key: ")
   )
   ```

3. **Keyring package** for persistent storage:
   ```r
   keyring::key_set("openai_api_key")
   provider <- list(
     name = "openai",
     api_key = keyring::key_get("openai_api_key")
   )
   ```

### Code Execution

The AI assistant **does not automatically execute generated code**. This is a security feature.

- Generated code is displayed in the UI and returned in the `code` field
- Review all code before running
- The "Run in Sandbox" button is intentionally disabled
- Future versions may add sandboxed execution with proper isolation

## Index Management

### View Index Information
```r
info <- iobr_ai_list_index()
# Output:
# Index: inst/ai/iobr_embeddings.rds
#   Entries: 523
#   Created: 2024-12-03 10:30:15
#   Provider: dummy (bag-of-words)
#   Size: 2.45 MB
```

### Reset Index
```r
iobr_ai_reset_index()
# Deletes inst/ai/iobr_embeddings.rds
```

### Custom Index Path
```r
# Save index to custom location
result <- iobr_ai_init(index_path = "/path/to/my_index.rds")

# Query from custom index
response <- iobr_ai_query(
  "my query",
  index_path = "/path/to/my_index.rds"
)
```

### Sharing Indexes

Indexes are portable RDS files that can be:
- Downloaded from the Shiny app
- Copied to other machines
- Committed to the package (if size permits and appropriate)

**Note**: Indexes created with different providers may not be compatible for querying. Use the same provider for index creation and querying.

## Provider Configuration Details

### Dummy Provider
- **Embeddings**: Local bag-of-words using term frequency
- **Chat**: Template-based keyword matching
- **Pros**: No API key, works offline, fast
- **Cons**: Limited semantic understanding, basic responses
- **Best for**: Testing, development, offline demos

### OpenAI Provider
- **Embeddings API**: `text-embedding-ada-002` (1536 dimensions)
- **Chat API**: `gpt-3.5-turbo` or `gpt-4`
- **Cost**: Pay per token (see [OpenAI pricing](https://openai.com/pricing))
- **Rate limits**: Varies by tier
- **Best for**: Production use, high-quality responses

### Hugging Face Provider
- **Inference API**: Free tier available
- **Embeddings**: Sentence transformers (384-768 dimensions typical)
- **Chat**: Open-source models (Mistral, Llama, etc.)
- **Rate limits**: More restrictive on free tier
- **Best for**: Cost-conscious deployments, open-source preference

## Customization

### Chunk Size
Adjust the granularity of document chunks:
```r
result <- iobr_ai_init(chunk_size = 1200)  # Larger chunks
result <- iobr_ai_init(chunk_size = 400)   # Smaller chunks
```

### Number of Retrieved Contexts
Control how much context is sent to the LLM:
```r
response <- iobr_ai_query("my query", top_k = 10)  # More context
response <- iobr_ai_query("my query", top_k = 3)   # Less context
```

### Temperature
Control response randomness:
```r
response <- iobr_ai_query("my query", temperature = 0.0)  # Deterministic
response <- iobr_ai_query("my query", temperature = 0.7)  # More creative
```

## Troubleshooting

### Index not found
```
Error: Index not found at inst/ai/iobr_embeddings.rds
```
**Solution**: Run `iobr_ai_init()` first to create the index.

### API key error
```
Error: OpenAI provider requires api_key
```
**Solution**: Configure provider with valid API key:
```r
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = "sk-..."
))
```

### HTTP timeout or rate limit
The system includes automatic retry with exponential backoff, but you may encounter:
```
Request error on attempt 1: ... Retrying in 2 seconds...
```
**Solution**: Wait for retry or check API status. For persistent issues, verify:
- API key is valid
- Network connectivity
- API rate limits not exceeded
- Base URL is correct

### Large index size
Indexes can be several MB depending on package size.
**Solutions**:
- Use `.Rbuildignore` to exclude from package builds if needed
- Store indexes separately for large packages
- Reduce chunk_size to create fewer chunks

## Future Enhancements

Potential improvements for future versions:
- Sandboxed code execution environment
- Support for additional providers (Anthropic Claude, Cohere, etc.)
- Incremental index updates
- Multi-modal support (images, plots)
- Conversation history
- Fine-tuned models on IOBR-specific content
- Integration with package help system (`?iobr_ai_query`)

## Technical Details

### Index Structure
The index RDS file contains:
```r
list(
  entries = list(
    list(id = 1, text = "...", source = "README.md", embedding = c(...)),
    list(id = 2, text = "...", source = "R/function.R", embedding = c(...)),
    ...
  ),
  vocab = c("word1", "word2", ...),  # For dummy provider
  created_at = <POSIXct>,
  provider_meta = list(
    name = "dummy",
    model_embeddings = "bag-of-words"
  )
)
```

### HTTP Request Flow
1. Build request body (JSON)
2. Add headers (Authorization, Content-Type, custom)
3. Send POST request via httr
4. Check status code
5. If 429/5xx: exponential backoff retry
6. Parse JSON response
7. Extract embeddings/chat content

### Cosine Similarity
Relevance scoring uses cosine similarity:
```
similarity = (A · B) / (||A|| * ||B||)
```
Where A is the query embedding and B is a chunk embedding.

## License

This AI assistant feature is part of the IOBR package and follows the same GPL-3 license.

## Support

For issues, questions, or feature requests:
- GitHub Issues: https://github.com/IOBR/IOBR/issues
- Documentation: https://iobr.github.io/book/

## Acknowledgments

The AI assistant uses:
- httr for HTTP requests
- jsonlite for JSON parsing
- shiny for web interface
- Provider APIs: OpenAI, Hugging Face, and others
