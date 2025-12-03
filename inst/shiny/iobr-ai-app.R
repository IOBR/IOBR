# IOBR AI Assistant - Shiny App
# This app provides a web interface for the IOBR AI assistant

library(shiny)

# UI
ui <- fluidPage(
  titlePanel("IOBR AI Assistant"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Provider Configuration"),
      selectInput("provider_name", "Provider:",
                  choices = c("dummy", "openai", "huggingface", "anthropic", "custom"),
                  selected = "dummy"),
      
      conditionalPanel(
        condition = "input.provider_name != 'dummy'",
        textInput("api_key", "API Key:", 
                  placeholder = "Enter your API key",
                  value = ""),
        passwordInput("api_key_secure", "API Key (secure):", 
                     placeholder = "Or use secure input"),
        textInput("base_url", "Base URL:", 
                  placeholder = "Optional, uses provider default"),
        textInput("model_embeddings", "Embeddings Model:", 
                  placeholder = "Optional, uses provider default"),
        textInput("model_chat", "Chat Model:", 
                  placeholder = "Optional, uses provider default"),
        textInput("extra_headers", "Extra Headers (JSON):", 
                  placeholder = '{"key": "value"}')
      ),
      
      hr(),
      
      h3("Index Management"),
      textInput("pkg_path", "Package Path:", value = "."),
      numericInput("chunk_size", "Chunk Size:", value = 800, min = 100, max = 5000),
      actionButton("create_index", "Create/Update Index", 
                   class = "btn-primary"),
      actionButton("list_index", "Show Index Info", 
                   class = "btn-info"),
      actionButton("reset_index", "Reset Index", 
                   class = "btn-danger"),
      
      hr(),
      
      downloadButton("download_index", "Download Index"),
      fileInput("upload_index", "Upload Index", 
                accept = ".rds"),
      
      hr(),
      
      p(style = "font-size: 12px; color: #666;",
        "Local dummy mode: no API key needed. For real providers, configure above.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Chat",
                 br(),
                 textAreaInput("query", "Ask a Question:", 
                              placeholder = "e.g., How do I use CIBERSORT for deconvolution?",
                              width = "100%", rows = 3),
                 
                 fluidRow(
                   column(6,
                          numericInput("top_k", "Context Chunks:", value = 5, min = 1, max = 20)
                   ),
                   column(6,
                          numericInput("temperature", "Temperature:", value = 0.2, 
                                      min = 0, max = 2, step = 0.1)
                   )
                 ),
                 
                 actionButton("submit_query", "Submit Query", 
                             class = "btn-success btn-lg"),
                 
                 hr(),
                 
                 h4("Response:"),
                 verbatimTextOutput("response_text"),
                 
                 conditionalPanel(
                   condition = "output.has_code",
                   h4("R Code:"),
                   p(style = "color: #d9534f;", 
                     "⚠️ Code is NOT automatically executed. Review before running."),
                   verbatimTextOutput("code_blocks"),
                   actionButton("run_sandbox", "Run in Sandbox", 
                               class = "btn-warning", disabled = TRUE),
                   p(style = "font-size: 12px; color: #666;",
                     "Sandbox execution is disabled. In future versions, this will allow safe code execution in an isolated environment.")
                 ),
                 
                 hr(),
                 
                 h4("Retrieved Sources:"),
                 uiOutput("retrieved_sources")
        ),
        
        tabPanel("Index Info",
                 br(),
                 verbatimTextOutput("index_info_display")
        ),
        
        tabPanel("Help",
                 br(),
                 h3("Quick Start"),
                 p("1. Choose a provider (start with 'dummy' for local testing)"),
                 p("2. Create an index by clicking 'Create/Update Index'"),
                 p("3. Ask questions in the Chat tab"),
                 
                 h3("Provider Setup"),
                 h4("Dummy Provider"),
                 p("No configuration needed. Uses local bag-of-words embeddings and template-based responses."),
                 
                 h4("OpenAI"),
                 p("Requires API key from https://platform.openai.com/api-keys"),
                 p("Default models: text-embedding-ada-002, gpt-3.5-turbo"),
                 
                 h4("Hugging Face"),
                 p("Requires API key from https://huggingface.co/settings/tokens"),
                 p("Default models: sentence-transformers/all-MiniLM-L6-v2, mistralai/Mistral-7B-Instruct-v0.1"),
                 
                 h4("Security Notes"),
                 tags$ul(
                   tags$li("API keys are stored in memory only during the session"),
                   tags$li("Never commit API keys to version control"),
                   tags$li("Use environment variables for production deployments"),
                   tags$li("Code execution is disabled by default for security")
                 ),
                 
                 h3("Programmatic Usage"),
                 p("You can also use the AI assistant programmatically:"),
                 pre('
# Initialize with dummy provider
result <- iobr_ai_init()

# Query the assistant
response <- iobr_ai_query("How do I calculate signature scores?")
cat(response$answer)

# View retrieved sources
print(response$retrieved)
                 ')
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    index_info = NULL,
    last_response = NULL,
    code_blocks = NULL
  )
  
  # Build provider config from inputs
  get_provider <- reactive({
    provider <- list(name = input$provider_name)
    
    if (input$provider_name != "dummy") {
      # Use secure password input if available, otherwise regular
      api_key <- if (nchar(input$api_key_secure) > 0) {
        input$api_key_secure
      } else {
        input$api_key
      }
      
      if (nchar(api_key) > 0) {
        provider$api_key <- api_key
      }
      
      if (nchar(input$base_url) > 0) {
        provider$base_url <- input$base_url
      }
      
      if (nchar(input$model_embeddings) > 0) {
        provider$model_embeddings <- input$model_embeddings
      }
      
      if (nchar(input$model_chat) > 0) {
        provider$model_chat <- input$model_chat
      }
      
      if (nchar(input$extra_headers) > 0) {
        tryCatch({
          provider$headers <- jsonlite::fromJSON(input$extra_headers)
        }, error = function(e) {
          showNotification(paste("Invalid headers JSON:", e$message), type = "error")
        })
      }
    }
    
    provider
  })
  
  # Create index
  observeEvent(input$create_index, {
    tryCatch({
      showNotification("Creating index... This may take a moment.", 
                      duration = NULL, id = "create_index_msg")
      
      provider <- get_provider()
      
      result <- iobr_ai_init(
        pkg_path = input$pkg_path,
        provider = provider,
        chunk_size = input$chunk_size
      )
      
      removeNotification(id = "create_index_msg")
      
      values$index_info <- result
      
      showNotification(
        sprintf("Index created: %d chunks from %d documents", 
                result$num_chunks, result$num_docs),
        type = "message"
      )
    }, error = function(e) {
      removeNotification(id = "create_index_msg")
      showNotification(paste("Error creating index:", e$message), type = "error")
    })
  })
  
  # List index
  observeEvent(input$list_index, {
    tryCatch({
      info <- iobr_ai_list_index()
      values$index_info <- info
    }, error = function(e) {
      showNotification(paste("Error listing index:", e$message), type = "error")
    })
  })
  
  # Reset index
  observeEvent(input$reset_index, {
    showModal(modalDialog(
      title = "Confirm Reset",
      "Are you sure you want to delete the index? This cannot be undone.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_reset", "Delete", class = "btn-danger")
      )
    ))
  })
  
  observeEvent(input$confirm_reset, {
    tryCatch({
      iobr_ai_reset_index()
      values$index_info <- NULL
      removeModal()
      showNotification("Index deleted", type = "message")
    }, error = function(e) {
      showNotification(paste("Error deleting index:", e$message), type = "error")
    })
  })
  
  # Submit query
  observeEvent(input$submit_query, {
    if (nchar(trimws(input$query)) == 0) {
      showNotification("Please enter a question", type = "warning")
      return()
    }
    
    tryCatch({
      showNotification("Processing query...", duration = NULL, id = "query_msg")
      
      provider <- get_provider()
      
      result <- iobr_ai_query(
        query = input$query,
        provider = provider,
        top_k = input$top_k,
        temperature = input$temperature
      )
      
      removeNotification(id = "query_msg")
      
      values$last_response <- result
      values$code_blocks <- result$code
      
      showNotification("Response received", type = "message")
      
    }, error = function(e) {
      removeNotification(id = "query_msg")
      showNotification(paste("Error processing query:", e$message), type = "error")
    })
  })
  
  # Download index
  output$download_index <- downloadHandler(
    filename = function() {
      paste0("iobr_embeddings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
      if (file.exists(index_path)) {
        file.copy(index_path, file)
      } else {
        stop("Index file not found")
      }
    }
  )
  
  # Upload index
  observeEvent(input$upload_index, {
    tryCatch({
      if (is.null(input$upload_index)) return()
      
      index_path <- file.path("inst", "ai", "iobr_embeddings.rds")
      index_dir <- dirname(index_path)
      if (!dir.exists(index_dir)) {
        dir.create(index_dir, recursive = TRUE)
      }
      
      file.copy(input$upload_index$datapath, index_path, overwrite = TRUE)
      
      showNotification("Index uploaded successfully", type = "message")
      
      # Update info
      info <- iobr_ai_list_index()
      values$index_info <- info
      
    }, error = function(e) {
      showNotification(paste("Error uploading index:", e$message), type = "error")
    })
  })
  
  # Output: Response text
  output$response_text <- renderText({
    if (!is.null(values$last_response)) {
      values$last_response$answer
    } else {
      "No response yet. Submit a query to get started."
    }
  })
  
  # Output: Code blocks
  output$code_blocks <- renderText({
    if (!is.null(values$code_blocks) && length(values$code_blocks) > 0) {
      paste(values$code_blocks, collapse = "\n\n---\n\n")
    } else {
      ""
    }
  })
  
  # Output: Has code (for conditional panel)
  output$has_code <- reactive({
    !is.null(values$code_blocks) && length(values$code_blocks) > 0
  })
  outputOptions(output, "has_code", suspendWhenHidden = FALSE)
  
  # Output: Retrieved sources
  output$retrieved_sources <- renderUI({
    if (!is.null(values$last_response) && !is.null(values$last_response$retrieved)) {
      retrieved <- values$last_response$retrieved
      
      lapply(seq_len(nrow(retrieved)), function(i) {
        div(
          style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px; border-radius: 5px;",
          h5(sprintf("Source %d: %s (Score: %.3f)", i, retrieved$source[i], retrieved$score[i])),
          tags$details(
            tags$summary("Show text"),
            pre(retrieved$text[i])
          )
        )
      })
    } else {
      p("No sources retrieved yet.")
    }
  })
  
  # Output: Index info
  output$index_info_display <- renderText({
    if (!is.null(values$index_info)) {
      info <- values$index_info
      
      sprintf(
        "Index Information:\n\n  Path: %s\n  Entries: %s\n  Documents: %s\n  Created: %s\n  Provider: %s\n  Size: %.2f MB",
        info$index_path %||% info$path,
        info$num_entries %||% info$num_chunks,
        info$num_docs %||% "N/A",
        info$created_at,
        info$provider %||% info$provider_name,
        info$file_size_mb %||% 0
      )
    } else {
      "No index information available. Click 'Show Index Info' or create an index."
    }
  })
}

# Run app
shinyApp(ui = ui, server = server)
