# Example: Using the IOBR AI Assistant
# 
# This script demonstrates how to use the IOBR AI assistant programmatically.
# Make sure you have set up your API key before running this script.

library(IOBR)

# Example 1: Configure Provider (OpenAI)
# ----------------------------------------
# Set your API key as an environment variable:
# Sys.setenv(OPENAI_API_KEY = "your-key-here")

provider <- list(
  name = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY")
)

# Validate the provider configuration
validated_provider <- iobr_ai_configure_provider(provider)
print("Provider configured successfully!")

# Example 2: Initialize the Index
# --------------------------------
# This creates an index of all package documentation and code.
# Only needs to be done once (or when you want to update the index).

# Uncomment to create index (requires valid API key):
# iobr_ai_init(
#   pkg_path = ".",
#   provider = provider,
#   chunk_size = 800,
#   progress_callback = function(msg) cat(msg, "\n")
# )

# Example 3: Query the Assistant
# -------------------------------
# Ask questions about IOBR functionality

# Example query 1: General information
result1 <- iobr_ai_query(
  query = "What is IOBR and what can it do?",
  provider = provider,
  top_k = 3
)

cat("\n=== Answer ===\n")
cat(result1$answer, "\n")

cat("\n=== Sources ===\n")
for (i in seq_along(result1$sources)) {
  src <- result1$sources[[i]]
  cat(sprintf("%d. %s (similarity: %.1f%%)\n", 
              i, src$source, src$similarity * 100))
}

# Example query 2: Specific function usage
result2 <- iobr_ai_query(
  query = "How do I use CIBERSORT for TME deconvolution? Show me code.",
  provider = provider,
  top_k = 5
)

cat("\n=== Answer ===\n")
cat(result2$answer, "\n")

if (length(result2$code) > 0) {
  cat("\n=== Generated Code ===\n")
  for (i in seq_along(result2$code)) {
    cat(sprintf("\n--- Code Block %d ---\n", i))
    cat(result2$code[[i]], "\n")
  }
}

# Example 4: Using Different Providers
# -------------------------------------

# Anthropic (Claude)
provider_anthropic <- list(
  name = "anthropic",
  api_key = Sys.getenv("ANTHROPIC_API_KEY")
)

# Hugging Face
provider_hf <- list(
  name = "huggingface",
  api_key = Sys.getenv("HUGGINGFACE_API_KEY"),
  model_embeddings = "sentence-transformers/all-MiniLM-L6-v2",
  model_chat = "mistralai/Mistral-7B-Instruct-v0.1"
)

# Custom provider
provider_custom <- list(
  name = "custom",
  api_key = "your-key",
  base_url = "https://your-api.com/v1/embeddings",
  model_embeddings = "your-embedding-model",
  model_chat = "your-chat-model"
)

# Example 5: Check Index Status
# ------------------------------
index_info <- iobr_ai_list_index()
print(index_info)

# Example 6: Reset Index
# -----------------------
# Uncomment to reset the index:
# iobr_ai_reset_index()

# Example 7: Launch the Shiny App
# --------------------------------
# Uncomment to launch the interactive interface:
# shiny::runApp(system.file("shiny", package = "IOBR"))
