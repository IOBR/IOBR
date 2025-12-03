# IOBR AI Assistant - Example Usage Script
# This script demonstrates how to use the AI assistant with the dummy provider (no API keys needed)

# Load the package
library(IOBR)

# Example 1: Initialize the AI assistant with dummy provider
cat("=== Example 1: Initialize AI Assistant ===\n\n")

metadata <- iobr_ai_init()
cat("Index created with:\n")
cat("- Number of chunks:", metadata$num_chunks, "\n")
cat("- Number of documents:", metadata$num_docs, "\n")
cat("- Provider:", metadata$provider, "\n\n")

# Example 2: Query the assistant
cat("=== Example 2: Query the Assistant ===\n\n")

response <- iobr_ai_query("What is IOBR used for?")
cat("Answer:\n", response$answer, "\n\n")

cat("Retrieved sources:\n")
for (i in seq_along(response$retrieved)) {
  src <- response$retrieved[[i]]
  cat(sprintf("%d. %s (score: %.3f)\n", i, src$source, src$score))
}
cat("\n")

# Example 3: View index information
cat("=== Example 3: View Index Info ===\n\n")
iobr_ai_list_index()
cat("\n")

# Example 4: Query with different parameters
cat("=== Example 4: Query with More Sources ===\n\n")

response2 <- iobr_ai_query(
  "How do I perform deconvolution?",
  top_k = 10
)
cat("Answer:\n", response2$answer, "\n\n")

# Example 5: Using OpenAI (if API key is available)
cat("=== Example 5: Using OpenAI (Optional) ===\n\n")

if (Sys.getenv("OPENAI_API_KEY") != "") {
  cat("OpenAI API key found. Configuring provider...\n")
  
  provider_openai <- iobr_ai_configure_provider(list(
    name = "openai",
    api_key = Sys.getenv("OPENAI_API_KEY")
  ))
  
  # Re-initialize with OpenAI
  iobr_ai_init(provider = provider_openai)
  
  # Query with OpenAI
  response_openai <- iobr_ai_query(
    "Explain the CIBERSORT method",
    provider = provider_openai
  )
  
  cat("OpenAI Answer:\n", response_openai$answer, "\n")
} else {
  cat("OpenAI API key not found. Set OPENAI_API_KEY environment variable to test.\n")
  cat("Example: Sys.setenv(OPENAI_API_KEY = 'your-key-here')\n")
}

cat("\n=== Examples Complete ===\n")
cat("For more information, see inst/ai/README.md\n")
