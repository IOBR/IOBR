# Testing IOBR AI Assistant

This document provides instructions for testing the IOBR AI Assistant MVP.

## Prerequisites

Install the required dependencies:

```r
install.packages(c("shiny", "httr", "jsonlite", "testthat"))
```

## Local Testing (No API Keys Required)

The AI assistant includes a "dummy" provider that works completely offline for testing.

### 1. Test Programmatically

```r
# Load the package (assuming you're in the package directory)
devtools::load_all()

# Initialize with dummy provider
index_info <- iobr_ai_init()

# This will:
# - Collect documentation from README, vignettes, R files, man pages
# - Chunk the text into ~800 character pieces
# - Generate dummy embeddings using bag-of-words
# - Save index to inst/ai/iobr_embeddings.rds

# Query the assistant
response <- iobr_ai_query("How do I analyze the tumor microenvironment?")

# View the response
cat(response$answer)

# View retrieved sources
for (i in seq_along(response$retrieved)) {
  src <- response$retrieved[[i]]
  cat(sprintf("\nSource %d: %s (score: %.3f)\n", i, src$source, src$score))
  cat(substr(src$text, 1, 150), "...\n")
}
```

### 2. Test via Shiny App

```r
# Load the package
devtools::load_all()

# Run the Shiny app
shiny::runApp("inst/shiny/iobr-ai-app.R")
```

In the Shiny app:

1. **Provider Configuration**:
   - The app starts with "dummy" provider selected by default
   - Click "Configure Provider" (no API key needed for dummy)

2. **Create Index**:
   - Click "Create Index" button
   - Wait for the process to complete (may take 1-2 minutes)
   - You'll see progress messages in the app

3. **Query the Assistant**:
   - Enter a question like: "What is CIBERSORT?"
   - Click "Submit Query"
   - View the answer and retrieved sources

4. **Index Management**:
   - Click "View Index Info" to see index statistics
   - Use "Download Index" to save the index file
   - Use "Upload Index" to restore a saved index
   - Use "Reset Index" to delete and start fresh

## Testing with Real Providers

### OpenAI

```r
# Set your API key (never commit this!)
Sys.setenv(OPENAI_API_KEY = "sk-your-key-here")

# Configure OpenAI provider
provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Create index with OpenAI embeddings
index_info <- iobr_ai_init(provider = provider)

# Query with OpenAI
response <- iobr_ai_query(
  "What is CIBERSORT and how do I use it?",
  provider = provider
)

cat(response$answer)
```

In the Shiny app:
1. Select "openai" from the provider dropdown
2. Paste your API key (starts with "sk-...")
3. Click "Configure Provider"
4. Click "Create Index"
5. Submit queries

**Note**: OpenAI charges per API call. Creating an index for IOBR documentation may cost a few cents. Queries cost a few cents each.

### Hugging Face

```r
# Set your HF token
Sys.setenv(HF_TOKEN = "hf_your_token_here")

# Configure Hugging Face provider
provider <- iobr_ai_configure_provider(list(
  name = "huggingface",
  api_key = Sys.getenv("HF_TOKEN")
))

# Create index
index_info <- iobr_ai_init(provider = provider)

# Query
response <- iobr_ai_query(
  "How do I perform survival analysis?",
  provider = provider
)
```

**Note**: Hugging Face Inference API is free but has rate limits. Large requests may timeout. Use smaller `batch_size` if needed:

```r
index_info <- iobr_ai_init(
  provider = provider,
  chunk_size = 400  # Smaller chunks
)
```

## Running Unit Tests

```r
# Run all tests
devtools::test()

# Run only AI tests
testthat::test_file("tests/testthat/test-ai.R")
```

The tests use the dummy provider and don't require API keys. They cover:

- Provider configuration
- Dummy embeddings generation
- Dummy chat responses
- Document collection
- Text chunking
- Index creation and saving
- Similarity search
- High-level API functions
- Code extraction from responses
- Full end-to-end workflow

## Troubleshooting

### "Index not found"

Make sure you've run `iobr_ai_init()` first to create the index.

### "Package 'httr' is required"

Install the required packages:
```r
install.packages(c("httr", "jsonlite"))
```

### "Provider requires API key"

For testing without API keys, use the dummy provider:
```r
provider <- list(name = "dummy")
```

### Shiny app doesn't start

Make sure you have shiny installed:
```r
install.packages("shiny")
```

And run from the correct location:
```r
# If in package root:
shiny::runApp("inst/shiny/iobr-ai-app.R")

# Or if package is installed:
shiny::runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))
```

### Tests fail

Tests require the package to be loadable. Make sure you're in the package directory and run:
```r
devtools::load_all()
devtools::test()
```

## Test Scenarios

Here are some suggested test scenarios:

### Scenario 1: Basic Usage (Dummy Provider)
1. Initialize index with dummy provider
2. Query about TME analysis
3. Verify response contains relevant information from docs
4. Check that retrieved sources are listed with scores

### Scenario 2: Provider Configuration
1. Test OpenAI configuration with valid defaults
2. Test Hugging Face configuration with valid defaults
3. Test custom provider configuration
4. Verify validation errors for missing required fields

### Scenario 3: Index Management
1. Create index
2. View index info
3. Download index file
4. Reset index
5. Upload previously downloaded index
6. Verify index is restored correctly

### Scenario 4: Different Query Types
1. Query about specific functions (e.g., "How to use CIBERSORT?")
2. Query about general concepts (e.g., "What is TME?")
3. Query for code examples (e.g., "Show me code for survival analysis")
4. Verify code blocks are extracted and displayed

### Scenario 5: Edge Cases
1. Empty query (should error)
2. Very long query (should work)
3. Query before index created (should error)
4. Multiple consecutive queries (should all work)

## Security Testing

1. **API Key Handling**:
   - Verify API keys are not logged or displayed
   - Test that dummy provider works without any credentials

2. **Code Execution**:
   - Verify that code examples are NOT automatically executed
   - Confirm sandbox button is disabled
   - Test that code is only displayed, not run

3. **File Paths**:
   - Test with custom index paths
   - Verify no path traversal issues
   - Confirm files are created in expected locations

## Performance Testing

1. **Index Creation Time**:
   - Measure time to create index (dummy provider)
   - Expected: < 1 minute for IOBR docs

2. **Query Response Time**:
   - Measure time for query (dummy provider)
   - Expected: < 5 seconds

3. **Memory Usage**:
   - Monitor memory during index creation
   - Check index file size (should be < 50MB for IOBR)

## Expected Outputs

### Index Creation
```
=== IOBR AI Assistant Initialization ===
Using dummy provider (offline mode)

1. Collecting package documentation...
Collected 85 documentation files

2. Chunking texts...
Created 437 chunks from 85 documents

3. Creating embeddings index...
Generating embeddings for 437 chunks...
Index saved to: inst/ai/iobr_embeddings.rds

=== Initialization Complete ===
Index path: inst/ai/iobr_embeddings.rds
Total entries: 437
```

### Query Response
```
Retrieving relevant context...
Retrieved 5 relevant chunks
Generating response...

Answer: [AI-generated response about IOBR functionality]

Retrieved Sources:
- README.md (score: 0.892)
- R/deconvo_tme.R (score: 0.845)
- man/CIBERSORT.Rd (score: 0.823)
...
```

## Reporting Issues

If you encounter any issues:

1. Check this TESTING.md for troubleshooting steps
2. Review inst/ai/README.md for architecture details
3. Run unit tests to identify specific failures
4. Open an issue on GitHub with:
   - Steps to reproduce
   - Expected vs actual behavior
   - Provider used (dummy/openai/etc)
   - Error messages or logs

## Next Steps

After testing the MVP, consider:

1. Testing with your own domain-specific documentation
2. Experimenting with different chunk sizes and top_k values
3. Trying different LLM providers
4. Providing feedback on the interface and functionality
5. Suggesting additional features or improvements
