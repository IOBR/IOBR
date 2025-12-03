# IOBR AI Assistant - API Reference

Quick reference for all exported AI assistant functions.

## High-Level API (Start Here)

### iobr_ai_init()
Initialize AI assistant and create embeddings index.

```r
result <- iobr_ai_init(
  pkg_path = ".",                    # Package directory
  provider = NULL,                   # Provider config (default: dummy)
  index_path = NULL,                 # Index location (default: inst/ai/iobr_embeddings.rds)
  chunk_size = 800                   # Chunk size in characters
)

# Returns: list(index_path, num_chunks, num_docs, provider, created_at)
```

### iobr_ai_query()
Query the AI assistant with automatic context retrieval.

```r
response <- iobr_ai_query(
  query = "How do I use CIBERSORT?",  # User question
  index_path = NULL,                  # Index location (default: inst/ai/iobr_embeddings.rds)
  provider = NULL,                    # Provider config (default: dummy)
  top_k = 5,                          # Number of context chunks
  max_tokens = 800,                   # Max response tokens
  temperature = 0.2                   # Sampling temperature (0-2)
)

# Returns: list(answer, code, retrieved, raw)
# - answer: Text response from AI
# - code: Extracted R code blocks (if any)
# - retrieved: Data frame with id, source, text, score
# - raw: Raw API response
```

### iobr_ai_configure_provider()
Configure and validate provider settings.

```r
# Dummy provider (no API key)
provider <- iobr_ai_configure_provider(list(
  name = "dummy"
))

# OpenAI provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  model_embeddings = "text-embedding-ada-002",  # optional
  model_chat = "gpt-3.5-turbo"                  # optional
))

# Hugging Face provider
provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HF_API_KEY"),
  model_embeddings = "sentence-transformers/all-MiniLM-L6-v2",  # optional
  model_chat = "mistralai/Mistral-7B-Instruct-v0.1"             # optional
))

# Custom provider
provider <- iobr_ai_configure_provider(list(
  name = "custom",
  api_key = "your-key",
  base_url = "https://your-api.com/v1",
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model",
  headers = list("X-Custom" = "value")
))
```

### iobr_ai_list_index()
Display index metadata.

```r
info <- iobr_ai_list_index(
  index_path = NULL  # default: inst/ai/iobr_embeddings.rds
)

# Returns: list(path, num_entries, created_at, provider_name, model_embeddings, file_size_mb)
```

### iobr_ai_reset_index()
Delete the index file.

```r
iobr_ai_reset_index(
  index_path = NULL  # default: inst/ai/iobr_embeddings.rds
)

# Returns: TRUE if deleted, FALSE if not found
```

## Provider Functions (Advanced)

### send_embeddings()
Generate embeddings for texts.

```r
embeddings <- send_embeddings(
  texts = c("text1", "text2"),  # Character vector
  provider = provider,          # Provider configuration
  batch_size = 16              # Batch size for API calls
)

# Returns: List of numeric vectors (one per text)
```

### send_chat()
Get chat completion from LLM.

```r
response <- send_chat(
  system_msg = "You are a helpful assistant",
  user_msg = "What is R?",
  provider = provider,
  max_tokens = 800,
  temperature = 0.2
)

# Returns: list(content, raw)
# - content: Response text
# - raw: Full API response
```

## Indexing Functions (Advanced)

### collect_package_docs()
Gather documentation files from package.

```r
docs <- collect_package_docs(
  pkg_path = "."  # Package directory
)

# Returns: List of list(content, source)
# Collects: README.md, vignettes/*, R/*.R, man/*.Rd, examples/*
```

### chunk_texts()
Split documents into smaller chunks.

```r
chunks <- chunk_texts(
  docs = docs,        # From collect_package_docs()
  chunk_size = 800    # Target chunk size
)

# Returns: List of list(text, source)
```

### create_index()
Generate embeddings and save index.

```r
index_path <- create_index(
  chunks = chunks,                # From chunk_texts()
  provider = provider,            # Provider configuration
  index_path = NULL,              # default: inst/ai/iobr_embeddings.rds
  batch_size = 16                # Batch size for embeddings
)

# Returns: Path to saved index (invisibly)
# Creates: RDS file with list(entries, vocab, created_at, provider_meta)
```

### search_index()
Find relevant chunks for a query.

```r
results <- search_index(
  query_text = "tumor microenvironment",
  index = index,        # Loaded from readRDS(index_path)
  provider = provider,
  top_k = 5            # Number of results
)

# Returns: Data frame with id, source, text, score
```

## Utility Functions

### cosine_sim()
Calculate cosine similarity between vectors.

```r
similarity <- cosine_sim(
  vec1 = c(1, 2, 3),
  vec2 = c(4, 5, 6)
)

# Returns: Numeric value between -1 and 1
```

## Usage Patterns

### Pattern 1: Quick Start (Dummy Mode)
```r
library(IOBR)

# One-time setup
iobr_ai_init()

# Use repeatedly
response <- iobr_ai_query("your question here")
cat(response$answer)
```

### Pattern 2: OpenAI Mode
```r
library(IOBR)

# Configure provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Create index with OpenAI embeddings
iobr_ai_init(provider = provider)

# Query with OpenAI
response <- iobr_ai_query("your question", provider = provider)
```

### Pattern 3: Custom Workflow
```r
library(IOBR)

# Step-by-step control
provider <- list(name = "dummy")
docs <- collect_package_docs(".")
chunks <- chunk_texts(docs, chunk_size = 1000)
create_index(chunks, provider, index_path = "my_index.rds")

# Search directly
index <- readRDS("my_index.rds")
results <- search_index("query", index, provider, top_k = 10)

# Use results for custom processing
print(results)
```

### Pattern 4: Shiny App
```r
library(IOBR)
library(shiny)

# Launch interactive interface
runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))

# Or from source
runApp("inst/shiny/iobr-ai-app.R")
```

## Provider Comparison

| Feature | Dummy | OpenAI | Hugging Face |
|---------|-------|--------|--------------|
| API Key | No | Yes | Yes |
| Offline | Yes | No | No |
| Cost | Free | Paid | Free tier |
| Embedding Quality | Basic | High | Medium-High |
| Response Quality | Template | Excellent | Good |
| Speed | Instant | 2-5s | 3-10s |
| Best For | Testing | Production | Cost-conscious |

## Common Tasks

### Create Index
```r
iobr_ai_init()
```

### Query Assistant
```r
response <- iobr_ai_query("How do I calculate signature scores?")
cat(response$answer)
```

### View Sources
```r
print(response$retrieved)
```

### Switch Provider
```r
# From dummy to OpenAI
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Recreate index
iobr_ai_reset_index()
iobr_ai_init(provider = provider)
```

### Custom Chunk Size
```r
iobr_ai_init(chunk_size = 1200)  # Larger chunks
```

### More Context
```r
response <- iobr_ai_query("query", top_k = 10)  # More sources
```

### Lower Temperature (More Deterministic)
```r
response <- iobr_ai_query("query", temperature = 0.0)
```

### Higher Temperature (More Creative)
```r
response <- iobr_ai_query("query", temperature = 0.7)
```

## Error Handling

All functions include error handling with informative messages:

```r
# Missing index
tryCatch(
  iobr_ai_query("query"),
  error = function(e) message("Need to run iobr_ai_init() first")
)

# Invalid provider
tryCatch(
  iobr_ai_configure_provider(list(name = "invalid")),
  error = function(e) message("Provider not supported")
)

# Missing API key
tryCatch(
  iobr_ai_configure_provider(list(name = "openai")),
  error = function(e) message("OpenAI requires api_key")
)
```

## Help & Documentation

```r
# Function help
?iobr_ai_init
?iobr_ai_query
?send_embeddings
?send_chat

# Package help
help(package = "IOBR")

# AI documentation
file.show(system.file("ai", "README.md", package = "IOBR"))
```

## Tips & Best Practices

1. **Start with dummy mode** for testing
2. **Use environment variables** for API keys (never hardcode)
3. **Review code before execution** (it's not auto-executed for security)
4. **Adjust chunk_size** based on your documentation structure
5. **Use top_k=3-5** for focused answers, 10+ for comprehensive
6. **Lower temperature** for factual questions, higher for creative
7. **Check retrieved sources** to understand answer provenance
8. **Reset index** when package docs change significantly
9. **Download/backup index** before major changes
10. **Monitor API costs** when using paid providers

## Support

- Issues: https://github.com/IOBR/IOBR/issues
- Docs: https://iobr.github.io/book/
- Tests: `testthat::test_file("tests/testthat/test-ai.R")`
