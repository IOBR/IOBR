# IOBR AI Assistant - Quick Start Guide

Get started with the IOBR AI Assistant in 5 minutes!

## 1. Install Package

```r
# Install from this branch (if needed)
# remotes::install_github("IOBR/IOBR", ref = "ai")

library(IOBR)
```

## 2. Create Index (One-Time Setup)

```r
# Initialize with dummy provider (no API key needed!)
result <- iobr_ai_init()

# Wait for completion...
# ✓ Collected 85 documents
# ✓ Created 523 chunks
# ✓ Index saved to: inst/ai/iobr_embeddings.rds
```

## 3. Ask Questions

```r
# Query the assistant
response <- iobr_ai_query("How do I use CIBERSORT for deconvolution?")

# View the answer
cat(response$answer)

# See which documentation was used
print(response$retrieved)
```

## 4. Use the Web Interface

```r
# Launch Shiny app
shiny::runApp(system.file("shiny", "iobr-ai-app.R", package = "IOBR"))
```

Then:
1. Make sure "dummy" provider is selected
2. Click "Create/Update Index"
3. Type your question in the Chat tab
4. Click "Submit Query"
5. View the answer and sources!

## Next Steps

### Upgrade to OpenAI (Optional)

```r
# Get API key from https://platform.openai.com/api-keys
# Add to ~/.Renviron:
# OPENAI_API_KEY=sk-your-key-here

provider <- iobr_ai_configure_provider(list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
))

# Recreate index with better embeddings
iobr_ai_init(provider = provider)

# Query with GPT
response <- iobr_ai_query(
  "Explain TME deconvolution methods",
  provider = provider
)
```

### Example Questions to Try

- "How do I calculate signature scores?"
- "What deconvolution methods are available?"
- "How to use TIMER?"
- "Explain tumor microenvironment analysis"
- "Show me survival analysis examples"
- "How to visualize TME results?"

## Troubleshooting

### Index not found
```r
# Create it first
iobr_ai_init()
```

### API key error
```r
# Use environment variable
Sys.setenv(OPENAI_API_KEY = "sk-...")

# Or check it's set
Sys.getenv("OPENAI_API_KEY")
```

### Out of memory
```r
# Use larger chunks (fewer of them)
iobr_ai_init(chunk_size = 1500)
```

## Learn More

- Full guide: `inst/ai/README.md`
- Testing: `TESTING_INSTRUCTIONS.md`
- PR details: `PR_DESCRIPTION.md`
- Help: `?iobr_ai_init`

## Tips

✅ Start with dummy mode (free, offline, instant)  
✅ Query indices are fast (< 1 second locally)  
✅ No code is executed automatically (security)  
✅ Sources show where answers came from (provenance)  
✅ Can download/share indexes as RDS files

## Support

Issues: https://github.com/IOBR/IOBR/issues

Enjoy exploring IOBR with AI! 🚀
