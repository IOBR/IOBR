# IOBR AI Assistant

## Overview

The IOBR AI Assistant is an integrated Shiny-based application that provides Retrieval-Augmented Generation (RAG) capabilities for the IOBR package. It allows users to:

- Ask natural language questions about IOBR functionality
- Get AI-generated answers with code examples
- View source documentation that supports each answer
- Work with multiple LLM providers (OpenAI, Anthropic, Hugging Face, or custom endpoints)

## Architecture

### Components

1. **Provider Client** (`R/ai_provider.R`)
   - Provider-agnostic HTTP client for embeddings and chat APIs
   - Supports OpenAI, Anthropic, Hugging Face, and custom endpoints
   - Handles authentication, retry logic, and error handling

2. **Indexing and Retrieval** (`R/ai_index.R`)
   - Collects package content (README, vignettes, man pages, R code)
   - Chunks text into manageable pieces
   - Creates embeddings using provider API
   - Implements cosine similarity search
   - Saves/loads index from RDS files

3. **Core API** (`R/ai.R`)
   - Exported functions for programmatic access
   - `iobr_ai_init()`: Initialize index
   - `iobr_ai_query()`: Query the assistant
   - `iobr_ai_configure_provider()`: Configure providers
   - `iobr_ai_list_index()`: View index status
   - `iobr_ai_reset_index()`: Delete index

4. **Shiny UI** (`inst/shiny/iobr-ai-app.R`)
   - Interactive web interface
   - Provider configuration panel
   - Index management
   - Query interface with provenance display

### Data Flow

```
User Query → Embedding → Vector Search → Top-K Retrieval → 
Prompt Construction → LLM Chat → Answer + Code + Sources
```

## Getting Started

### Prerequisites

Install required packages:

```r
install.packages(c("shiny", "httr2", "jsonlite"))
```

### Provider Configuration

#### OpenAI

```r
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  # Optional overrides:
  base_url = "https://api.openai.com/v1",
  model_embeddings = "text-embedding-ada-002",
  model_chat = "gpt-3.5-turbo"
)
```

Set your API key:
```bash
export OPENAI_API_KEY="sk-..."
```

#### Anthropic (Claude)

```r
provider <- list(
  name = "anthropic",
  api_key = Sys.getenv("ANTHROPIC_API_KEY"),
  # Optional overrides:
  base_url = "https://api.anthropic.com/v1",
  model_chat = "claude-3-sonnet-20240229"
)
```

**Note:** Anthropic doesn't provide embeddings, so you'll need to use OpenAI or another provider for the embedding model.

Set your API key:
```bash
export ANTHROPIC_API_KEY="sk-ant-..."
```

#### Hugging Face

```r
provider <- list(
  name = "huggingface",
  api_key = Sys.getenv("HUGGINGFACE_API_KEY"),
  # Optional overrides:
  base_url = "https://api-inference.huggingface.co",
  model_embeddings = "sentence-transformers/all-MiniLM-L6-v2",
  model_chat = "mistralai/Mistral-7B-Instruct-v0.1"
)
```

Set your API key:
```bash
export HUGGINGFACE_API_KEY="hf_..."
```

#### Custom Provider

For any OpenAI-compatible API:

```r
provider <- list(
  name = "custom",
  api_key = "your-api-key",
  base_url = "https://your-api.com/v1/embeddings",  # For embeddings
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model",
  # Optional custom headers:
  headers = list(
    "X-Custom-Header" = "value"
  )
)
```

## Usage

### Programmatic Usage

#### Initialize Index

```r
library(IOBR)

# Configure provider
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)

# Create index
iobr_ai_init(
  pkg_path = ".",  # Path to IOBR package
  provider = provider,
  chunk_size = 800
)
```

This will:
1. Collect all package documentation and code
2. Split into chunks
3. Create embeddings for each chunk
4. Save index to `inst/ai/iobr_embeddings.rds`

#### Query the Assistant

```r
# Ask a question
result <- iobr_ai_query(
  query = "How do I perform TME deconvolution with CIBERSORT?",
  provider = provider,
  top_k = 5
)

# View answer
cat(result$answer)

# View code examples
if (length(result$code) > 0) {
  cat("Generated code:\n")
  cat(result$code[[1]])
}

# View sources
for (src in result$sources) {
  cat(sprintf("\nSource: %s (similarity: %.2f%%)\n", 
              src$source, src$similarity * 100))
}
```

#### Check Index Status

```r
# List index information
info <- iobr_ai_list_index()
print(info)

# Reset index
iobr_ai_reset_index()
```

### Shiny App Usage

#### Launch the App

```r
# From within R
shiny::runApp(system.file("shiny", package = "IOBR"))

# Or if in development
shiny::runApp("inst/shiny/iobr-ai-app.R")
```

#### Using the Interface

1. **Configure Provider**
   - Select provider (OpenAI, Anthropic, Hugging Face, or Custom)
   - Enter API key (or leave blank to use environment variable)
   - Optionally override base URL and model names
   - Click "Test Configuration" to validate

2. **Create Index**
   - Set package path (default: current directory)
   - Adjust chunk size if needed (default: 800 characters)
   - Click "Create Index"
   - Wait for completion (may take several minutes)

3. **Ask Questions**
   - Type your question in the text area
   - Adjust parameters:
     - Number of sources to retrieve (top_k)
     - Temperature for generation (0.0-1.0)
   - Click "Ask"
   - View answer, code examples, and sources

4. **Review Results**
   - Read the AI-generated answer
   - Copy any code examples for use
   - Expand source items to see relevant documentation
   - Check similarity scores to gauge relevance

## Security Considerations

### Code Execution

**IMPORTANT:** The AI assistant does NOT automatically execute generated code. All code is displayed for manual review and must be explicitly copied and run by the user.

This design ensures:
- No arbitrary code execution
- User review of all generated code
- Explicit consent before running code
- Auditability of what code is executed

### API Keys

- Never commit API keys to source control
- Use environment variables for keys
- The Shiny app supports runtime key entry
- Keys are not logged or persisted by the app

### Data Privacy

- All package content is processed locally
- Only embeddings and queries are sent to the provider API
- No user data is collected or transmitted beyond API calls
- Review your provider's data usage policy

## Troubleshooting

### "Index file not found"

You need to create an index first:
```r
iobr_ai_init(provider = your_provider_config)
```

### "API key is required"

Ensure you've set your API key:
```r
provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")  # Make sure this is set
)
```

### Rate Limit Errors

The provider client includes automatic retry with exponential backoff, but you may still hit rate limits with large indexes. Consider:
- Using a smaller chunk size
- Processing in batches
- Using a provider tier with higher rate limits

### Embedding Dimension Mismatch

Different embedding models produce different vector dimensions. If you change models, you must recreate the index.

### Poor Answer Quality

Try:
- Increasing `top_k` to retrieve more sources
- Adjusting `temperature` (lower = more focused, higher = more creative)
- Rephrasing your question
- Using a more capable chat model

## Extending the Assistant

### Adding New Providers

To add support for a new provider:

1. Edit `R/ai_provider.R`
2. Add a new case in `send_embedding()` and `send_chat()`
3. Handle the provider's specific API format
4. Update documentation

Example:
```r
} else if (provider_name == "my_new_provider") {
  endpoint <- paste0(provider_config$base_url, "/embeddings")
  # ... implement API calls
}
```

### Custom Content Sources

To index additional content:

1. Edit `collect_package_content()` in `R/ai_index.R`
2. Add logic to collect your content
3. Return list with `text`, `source`, and `type` fields

### UI Customization

The Shiny app can be customized by editing `inst/shiny/iobr-ai-app.R`. Consider:
- Adding visualization of embeddings (e.g., with UMAP)
- Implementing chat history
- Adding batch query processing
- Exporting results to files

## API Reference

See function documentation:
```r
?iobr_ai_init
?iobr_ai_query
?iobr_ai_configure_provider
?iobr_ai_list_index
?iobr_ai_reset_index
```

## License

The IOBR AI Assistant is part of the IOBR package and is distributed under the same GPL-3 license.

## Support

For issues or questions:
- GitHub Issues: https://github.com/IOBR/IOBR/issues
- Documentation: https://iobr.github.io/book/
