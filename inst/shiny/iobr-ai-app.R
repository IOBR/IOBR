# IOBR AI Assistant Shiny App
#
# This Shiny app provides a web interface for the IOBR AI assistant,
# allowing users to configure providers, create embeddings indices,
# and query the assistant with RAG.

library(shiny)

# UI Definition -----------------------------------------------------------

ui <- fluidPage(
  # Custom CSS
  tags$head(
    tags$style(HTML("
      .main-header {
        background-color: #2C3E50;
        color: white;
        padding: 20px;
        margin-bottom: 20px;
        border-radius: 5px;
      }
      .provider-config {
        background-color: #ECF0F1;
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 20px;
      }
      .code-block {
        background-color: #2C3E50;
        color: #ECF0F1;
        padding: 15px;
        border-radius: 5px;
        font-family: 'Courier New', monospace;
        white-space: pre-wrap;
        margin: 10px 0;
      }
      .source-item {
        background-color: #F8F9FA;
        border-left: 4px solid #3498DB;
        padding: 10px;
        margin: 10px 0;
        border-radius: 3px;
      }
      .score-badge {
        background-color: #3498DB;
        color: white;
        padding: 2px 8px;
        border-radius: 3px;
        font-size: 0.85em;
      }
      .warning-box {
        background-color: #FFF3CD;
        border: 1px solid #FFC107;
        border-radius: 5px;
        padding: 15px;
        margin: 10px 0;
      }
      .btn-sandbox {
        background-color: #95A5A6;
        cursor: not-allowed;
      }
    "))
  ),
  
  # Header
  div(class = "main-header",
      h1("IOBR AI Assistant"),
      p("RAG-powered documentation assistant for the IOBR package")
  ),
  
  # Main content
  fluidRow(
    # Left panel: Configuration
    column(4,
      wellPanel(
        h3("Provider Configuration"),
        div(class = "provider-config",
          selectInput("provider_name", "Provider:",
                      choices = c("dummy", "openai", "huggingface", "anthropic", "custom"),
                      selected = "dummy"),
          
          conditionalPanel(
            condition = "input.provider_name != 'dummy'",
            passwordInput("api_key", "API Key:", placeholder = "Enter your API key"),
            textInput("base_url", "Base URL:", placeholder = "Leave empty for default"),
            textInput("model_embeddings", "Embeddings Model:", placeholder = "Default model"),
            textInput("model_chat", "Chat Model:", placeholder = "Default model"),
            textAreaInput("headers", "Extra Headers (JSON):", 
                         placeholder = '{"X-Custom": "value"}',
                         rows = 2)
          ),
          
          actionButton("config_provider", "Configure Provider", 
                      class = "btn-primary btn-block")
        ),
        
        hr(),
        
        h4("Index Management"),
        actionButton("create_index", "Create Index", 
                    class = "btn-success btn-block"),
        br(),
        actionButton("view_index", "View Index Info", 
                    class = "btn-info btn-block"),
        br(),
        actionButton("reset_index", "Reset Index", 
                    class = "btn-warning btn-block"),
        br(),
        downloadButton("download_index", "Download Index", 
                      class = "btn-secondary btn-block"),
        br(),
        fileInput("upload_index", "Upload Index", 
                 accept = ".rds", 
                 buttonLabel = "Browse...",
                 placeholder = "No file chosen")
      )
    ),
    
    # Right panel: Query interface
    column(8,
      tabsetPanel(
        # Query tab
        tabPanel("Query",
          br(),
          textAreaInput("query_text", "Ask a question:", 
                       placeholder = "e.g., How do I analyze the tumor microenvironment with IOBR?",
                       rows = 3,
                       width = "100%"),
          
          fluidRow(
            column(6,
              numericInput("top_k", "Top K results:", value = 5, min = 1, max = 20)
            ),
            column(6,
              numericInput("temperature", "Temperature:", value = 0.2, min = 0, max = 1, step = 0.1)
            )
          ),
          
          actionButton("submit_query", "Submit Query", 
                      class = "btn-primary btn-lg btn-block"),
          
          hr(),
          
          h3("Response"),
          uiOutput("response_ui"),
          
          hr(),
          
          h4("Retrieved Context"),
          uiOutput("sources_ui")
        ),
        
        # About tab
        tabPanel("About",
          br(),
          h3("About IOBR AI Assistant"),
          p("This AI assistant uses Retrieval-Augmented Generation (RAG) to help you understand and use the IOBR package."),
          
          h4("How it works:"),
          tags$ol(
            tags$li("Documentation from the IOBR package is collected and chunked into smaller pieces"),
            tags$li("Each chunk is converted to an embedding vector using a language model"),
            tags$li("When you ask a question, the most relevant chunks are retrieved"),
            tags$li("The retrieved context is sent to a chat model to generate a helpful answer")
          ),
          
          h4("Providers:"),
          tags$ul(
            tags$li(tags$b("dummy:"), " Local offline mode using bag-of-words. No API key needed. Great for testing!"),
            tags$li(tags$b("openai:"), " Uses OpenAI's embedding and chat models (requires API key)"),
            tags$li(tags$b("huggingface:"), " Uses Hugging Face inference API (requires API key)"),
            tags$li(tags$b("anthropic:"), " Placeholder for future support"),
            tags$li(tags$b("custom:"), " Configure your own endpoints")
          ),
          
          div(class = "warning-box",
            h4("Security Note:"),
            p("Never commit API keys to source control. Use environment variables or secure configuration files."),
            p("The dummy provider is perfect for local testing without exposing any credentials.")
          ),
          
          h4("Sandbox Execution (Coming Soon)"),
          p("Future versions will include a sandboxed R environment for safely testing generated code."),
          actionButton("sandbox_btn", "Run in Sandbox (Disabled)", 
                      class = "btn-sandbox", disabled = TRUE),
          p(tags$small("Enable by setting up a secure sandbox environment with Docker or similar isolation."))
        ),
        
        # Logs tab
        tabPanel("Logs",
          br(),
          h3("Activity Log"),
          verbatimTextOutput("log_output")
        )
      )
    )
  )
)


# Server Logic ------------------------------------------------------------

server <- function(input, output, session) {
  
  # Reactive values
  rv <- reactiveValues(
    provider = list(name = "dummy"),
    log_messages = c("App started. Using dummy provider by default."),
    index_info = NULL,
    last_response = NULL
  )
  
  # Helper: Add log message
  add_log <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    rv$log_messages <- c(rv$log_messages, paste0("[", timestamp, "] ", msg))
  }
  
  # Configure provider
  observeEvent(input$config_provider, {
    tryCatch({
      provider <- list(name = input$provider_name)
      
      if (input$provider_name != "dummy") {
        if (nchar(input$api_key) > 0) {
          provider$api_key <- input$api_key
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
        if (nchar(input$headers) > 0) {
          provider$headers <- jsonlite::fromJSON(input$headers)
        }
      }
      
      # Validate
      provider <- iobr_ai_configure_provider(provider)
      rv$provider <- provider
      
      add_log(sprintf("Provider configured: %s", provider$name))
      showNotification(sprintf("Provider configured: %s", provider$name), 
                      type = "message")
    }, error = function(e) {
      add_log(sprintf("Error configuring provider: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Create index
  observeEvent(input$create_index, {
    tryCatch({
      add_log("Starting index creation...")
      showNotification("Creating index... This may take a few minutes.", 
                      duration = NULL, id = "index_creation")
      
      # Get package path (assume current working directory)
      pkg_path <- "."
      
      result <- iobr_ai_init(
        pkg_path = pkg_path,
        provider = rv$provider,
        index_path = NULL,  # Use default
        chunk_size = 800
      )
      
      rv$index_info <- result
      add_log(sprintf("Index created: %d entries", result$num_entries))
      
      removeNotification("index_creation")
      showNotification(sprintf("Index created successfully! (%d entries)", 
                              result$num_entries), 
                      type = "message")
    }, error = function(e) {
      removeNotification("index_creation")
      add_log(sprintf("Error creating index: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # View index info
  observeEvent(input$view_index, {
    tryCatch({
      info <- iobr_ai_list_index()
      if (!is.null(info)) {
        rv$index_info <- info
        add_log("Index info retrieved")
        showNotification(sprintf("Index has %d entries, created %s", 
                                info$num_entries, info$created_at), 
                        type = "message")
      } else {
        showNotification("No index found. Please create one first.", type = "warning")
      }
    }, error = function(e) {
      add_log(sprintf("Error viewing index: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Reset index
  observeEvent(input$reset_index, {
    tryCatch({
      result <- iobr_ai_reset_index()
      if (result) {
        rv$index_info <- NULL
        add_log("Index reset")
        showNotification("Index deleted successfully", type = "message")
      } else {
        showNotification("No index to delete", type = "warning")
      }
    }, error = function(e) {
      add_log(sprintf("Error resetting index: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
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
        add_log("Index downloaded")
      } else {
        stop("No index file found")
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
      add_log("Index uploaded")
      showNotification("Index uploaded successfully", type = "message")
    }, error = function(e) {
      add_log(sprintf("Error uploading index: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Submit query
  observeEvent(input$submit_query, {
    req(input$query_text)
    
    if (nchar(trimws(input$query_text)) == 0) {
      showNotification("Please enter a question", type = "warning")
      return()
    }
    
    tryCatch({
      add_log(sprintf("Query submitted: %s", substr(input$query_text, 1, 50)))
      showNotification("Processing query...", duration = NULL, id = "query_processing")
      
      result <- iobr_ai_query(
        query = input$query_text,
        index_path = NULL,  # Use default
        provider = rv$provider,
        top_k = input$top_k,
        max_tokens = 800,
        temperature = input$temperature
      )
      
      rv$last_response <- result
      add_log("Query completed")
      
      removeNotification("query_processing")
      showNotification("Query completed!", type = "message")
    }, error = function(e) {
      removeNotification("query_processing")
      add_log(sprintf("Error processing query: %s", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render response
  output$response_ui <- renderUI({
    if (is.null(rv$last_response)) {
      return(p("No response yet. Submit a query above."))
    }
    
    response <- rv$last_response
    
    # Split answer and code
    answer_parts <- list()
    
    # Main answer
    answer_parts[[length(answer_parts) + 1]] <- p(response$answer)
    
    # Code blocks if present
    if (!is.null(response$code) && length(response$code) > 0) {
      answer_parts[[length(answer_parts) + 1]] <- h4("Code Examples:")
      
      for (i in seq_along(response$code)) {
        answer_parts[[length(answer_parts) + 1]] <- div(
          class = "code-block",
          code(response$code[[i]])
        )
      }
      
      answer_parts[[length(answer_parts) + 1]] <- div(
        class = "warning-box",
        p(strong("Note:"), "Code is not automatically executed for safety. 
          Copy and run in your R console to test.")
      )
    }
    
    return(tagList(answer_parts))
  })
  
  # Render sources
  output$sources_ui <- renderUI({
    if (is.null(rv$last_response)) {
      return(NULL)
    }
    
    retrieved <- rv$last_response$retrieved
    
    if (is.null(retrieved) || length(retrieved) == 0) {
      return(p("No sources retrieved"))
    }
    
    source_items <- lapply(seq_along(retrieved), function(i) {
      item <- retrieved[[i]]
      div(
        class = "source-item",
        h5(
          sprintf("Source %d: %s", i, item$source),
          span(class = "score-badge", 
               sprintf("Score: %.3f", item$score))
        ),
        p(substr(item$text, 1, 300), 
          if (nchar(item$text) > 300) "..." else "")
      )
    })
    
    return(tagList(source_items))
  })
  
  # Render logs
  output$log_output <- renderText({
    paste(rev(rv$log_messages), collapse = "\n")
  })
}


# Run App -----------------------------------------------------------------

shinyApp(ui = ui, server = server)
