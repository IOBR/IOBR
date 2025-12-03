#' IOBR AI Assistant Shiny App
#'
#' Interactive interface for the IOBR AI assistant with RAG capabilities

library(shiny)

# UI Definition
ui <- fluidPage(
  titlePanel("IOBR AI Assistant"),
  
  tags$head(
    tags$style(HTML("
      .config-panel {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 20px;
      }
      .source-item {
        background-color: #e9ecef;
        padding: 10px;
        margin: 10px 0;
        border-radius: 5px;
        border-left: 3px solid #007bff;
      }
      .code-block {
        background-color: #f4f4f4;
        padding: 15px;
        border-radius: 5px;
        font-family: monospace;
        white-space: pre-wrap;
        margin: 10px 0;
      }
      .answer-box {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        border: 1px solid #dee2e6;
        margin: 20px 0;
      }
      .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        padding: 10px;
        border-radius: 5px;
        margin: 10px 0;
      }
      .similarity-score {
        display: inline-block;
        background-color: #28a745;
        color: white;
        padding: 2px 8px;
        border-radius: 3px;
        font-size: 0.85em;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Provider Configuration
      h3("Provider Configuration"),
      div(class = "config-panel",
        selectInput("provider_name", "Provider",
                   choices = c("OpenAI" = "openai",
                             "Anthropic" = "anthropic", 
                             "Hugging Face" = "huggingface",
                             "Custom" = "custom"),
                   selected = "openai"),
        
        passwordInput("api_key", "API Key",
                     placeholder = "Enter API key or use env variable"),
        
        textInput("base_url", "Base URL (optional)",
                 placeholder = "Leave empty for defaults"),
        
        textInput("model_embeddings", "Embeddings Model (optional)",
                 placeholder = "Leave empty for defaults"),
        
        textInput("model_chat", "Chat Model (optional)",
                 placeholder = "Leave empty for defaults"),
        
        actionButton("test_provider", "Test Configuration", 
                    class = "btn-primary btn-sm")
      ),
      
      # Index Management
      h3("Index Management"),
      div(class = "config-panel",
        textInput("pkg_path", "Package Path", value = "."),
        numericInput("chunk_size", "Chunk Size", value = 800, min = 200, max = 2000),
        actionButton("create_index", "Create Index", class = "btn-success"),
        hr(),
        verbatimTextOutput("index_status"),
        actionButton("refresh_status", "Refresh Status", class = "btn-sm"),
        actionButton("reset_index", "Reset Index", class = "btn-warning btn-sm")
      )
    ),
    
    mainPanel(
      width = 9,
      
      # Query Interface
      h3("Ask the AI Assistant"),
      div(class = "warning-box",
        strong("Security Notice:"),
        " Generated code is shown for review only. Do not execute untrusted code automatically."
      ),
      
      textAreaInput("query", "Your Question",
                   placeholder = "Ask anything about IOBR package...",
                   width = "100%", height = "100px"),
      
      fluidRow(
        column(6,
          numericInput("top_k", "Number of Sources to Retrieve", 
                      value = 5, min = 1, max = 20)
        ),
        column(6,
          numericInput("temperature", "Temperature", 
                      value = 0.2, min = 0, max = 1, step = 0.1)
        )
      ),
      
      actionButton("ask_btn", "Ask", class = "btn-primary btn-lg"),
      
      hr(),
      
      # Progress and status
      uiOutput("query_status"),
      
      # Answer Display
      conditionalPanel(
        condition = "output.has_answer",
        h3("Answer"),
        div(class = "answer-box",
          uiOutput("answer_display")
        ),
        
        # Code Display
        conditionalPanel(
          condition = "output.has_code",
          h3("Generated Code"),
          uiOutput("code_display"),
          actionButton("copy_code", "Copy Code", class = "btn-sm"),
          tags$small(class = "text-muted", 
                    " (Review code carefully before running in your R session)")
        ),
        
        # Sources Display
        h3("Sources Used"),
        uiOutput("sources_display")
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values
  rv <- reactiveValues(
    answer_result = NULL,
    index_info = NULL,
    provider_config = NULL,
    progress_msg = ""
  )
  
  # Get provider configuration
  get_provider_config <- reactive({
    api_key <- input$api_key
    if (api_key == "") {
      # Try to get from environment
      api_key <- switch(input$provider_name,
        "openai" = Sys.getenv("OPENAI_API_KEY"),
        "anthropic" = Sys.getenv("ANTHROPIC_API_KEY"),
        "huggingface" = Sys.getenv("HUGGINGFACE_API_KEY"),
        ""
      )
    }
    
    config <- list(
      name = input$provider_name,
      api_key = api_key
    )
    
    if (input$base_url != "") {
      config$base_url <- input$base_url
    }
    if (input$model_embeddings != "") {
      config$model_embeddings <- input$model_embeddings
    }
    if (input$model_chat != "") {
      config$model_chat <- input$model_chat
    }
    
    config
  })
  
  # Test provider configuration
  observeEvent(input$test_provider, {
    tryCatch({
      config <- get_provider_config()
      validated <- IOBR::iobr_ai_configure_provider(config)
      showNotification("Provider configuration is valid!", type = "message")
      rv$provider_config <- validated
    }, error = function(e) {
      showNotification(paste("Configuration error:", e$message), type = "error")
    })
  })
  
  # Create index
  observeEvent(input$create_index, {
    tryCatch({
      config <- get_provider_config()
      
      # Validate
      validated <- IOBR::iobr_ai_configure_provider(config)
      
      showNotification("Creating index... This may take several minutes.", 
                      type = "message", duration = NULL, id = "index_progress")
      
      # Progress callback
      progress_fn <- function(msg) {
        rv$progress_msg <- msg
        showNotification(msg, type = "message", duration = 3)
      }
      
      # Create index
      IOBR::iobr_ai_init(
        pkg_path = input$pkg_path,
        provider = validated,
        chunk_size = input$chunk_size,
        progress_callback = progress_fn
      )
      
      removeNotification("index_progress")
      showNotification("Index created successfully!", type = "message")
      
      # Refresh status
      refresh_index_status()
      
    }, error = function(e) {
      removeNotification("index_progress")
      showNotification(paste("Error creating index:", e$message), type = "error", duration = 10)
    })
  })
  
  # Refresh index status
  refresh_index_status <- function() {
    tryCatch({
      info <- IOBR::iobr_ai_list_index()
      rv$index_info <- info
    }, error = function(e) {
      rv$index_info <- NULL
    })
  }
  
  observeEvent(input$refresh_status, {
    refresh_index_status()
  })
  
  # Reset index
  observeEvent(input$reset_index, {
    showModal(modalDialog(
      title = "Confirm Reset",
      "Are you sure you want to delete the index? You will need to recreate it.",
      footer = tagList(
        actionButton("confirm_reset", "Yes, Reset"),
        modalButton("Cancel")
      )
    ))
  })
  
  observeEvent(input$confirm_reset, {
    tryCatch({
      IOBR::iobr_ai_reset_index()
      rv$index_info <- NULL
      showNotification("Index reset successfully", type = "message")
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error resetting index:", e$message), type = "error")
    })
  })
  
  # Ask question
  observeEvent(input$ask_btn, {
    req(input$query)
    
    tryCatch({
      config <- get_provider_config()
      validated <- IOBR::iobr_ai_configure_provider(config)
      
      showNotification("Processing query...", type = "message", 
                      duration = NULL, id = "query_progress")
      
      result <- IOBR::iobr_ai_query(
        query = input$query,
        provider = validated,
        top_k = input$top_k,
        temperature = input$temperature
      )
      
      rv$answer_result <- result
      removeNotification("query_progress")
      showNotification("Answer received!", type = "message")
      
    }, error = function(e) {
      removeNotification("query_progress")
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
      rv$answer_result <- NULL
    })
  })
  
  # Index status display
  output$index_status <- renderText({
    # Trigger reactive
    input$refresh_status
    input$create_index
    
    if (is.null(rv$index_info)) {
      tryCatch({
        info <- IOBR::iobr_ai_list_index()
        rv$index_info <- info
      }, error = function(e) {
        return("No index found.\nCreate an index to get started.")
      })
    }
    
    if (!is.null(rv$index_info)) {
      paste0(
        "Index Status: READY\n",
        "Total entries: ", rv$index_info$total_entries, "\n",
        "Unique sources: ", rv$index_info$unique_sources
      )
    } else {
      "No index found.\nCreate an index to get started."
    }
  })
  
  # Check if we have an answer
  output$has_answer <- reactive({
    !is.null(rv$answer_result)
  })
  outputOptions(output, "has_answer", suspendWhenHidden = FALSE)
  
  # Check if we have code
  output$has_code <- reactive({
    !is.null(rv$answer_result) && length(rv$answer_result$code) > 0
  })
  outputOptions(output, "has_code", suspendWhenHidden = FALSE)
  
  # Query status
  output$query_status <- renderUI({
    if (!is.null(rv$progress_msg) && rv$progress_msg != "") {
      div(class = "alert alert-info", rv$progress_msg)
    }
  })
  
  # Answer display
  output$answer_display <- renderUI({
    req(rv$answer_result)
    
    # Convert markdown-style formatting to HTML
    answer_text <- rv$answer_result$answer
    
    # Simple markdown conversion
    answer_html <- gsub("\n\n", "<br><br>", answer_text)
    answer_html <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", answer_html)
    answer_html <- gsub("```r\\n(.+?)\\n```", "<pre>\\1</pre>", answer_html)
    
    HTML(answer_html)
  })
  
  # Code display
  output$code_display <- renderUI({
    req(rv$answer_result)
    req(length(rv$answer_result$code) > 0)
    
    code_blocks <- lapply(seq_along(rv$answer_result$code), function(i) {
      div(
        h4(paste("Code Block", i)),
        div(class = "code-block", rv$answer_result$code[[i]]),
        tags$hr()
      )
    })
    
    do.call(tagList, code_blocks)
  })
  
  # Sources display
  output$sources_display <- renderUI({
    req(rv$answer_result)
    req(length(rv$answer_result$sources) > 0)
    
    source_items <- lapply(seq_along(rv$answer_result$sources), function(i) {
      src <- rv$answer_result$sources[[i]]
      
      div(class = "source-item",
        h5(
          paste0(i, ". ", src$source),
          tags$span(class = "similarity-score", 
                   sprintf("%.1f%%", src$similarity * 100))
        ),
        tags$small(class = "text-muted", paste("Type:", src$type)),
        tags$details(
          tags$summary("View source text"),
          div(style = "margin-top: 10px; padding: 10px; background: white; border-radius: 3px;",
            tags$pre(substr(src$text, 1, 500), 
                    if (nchar(src$text) > 500) "..." else "")
          )
        )
      )
    })
    
    do.call(tagList, source_items)
  })
  
  # Copy code functionality (using clipr or manual copy)
  observeEvent(input$copy_code, {
    req(rv$answer_result)
    req(length(rv$answer_result$code) > 0)
    
    all_code <- paste(rv$answer_result$code, collapse = "\n\n")
    
    # Show modal with code to copy
    showModal(modalDialog(
      title = "Copy Code",
      tags$p("Copy the code below:"),
      tags$textarea(
        all_code,
        style = "width: 100%; height: 300px; font-family: monospace;"
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  # Initialize
  observe({
    refresh_index_status()
  })
}

# Run the app
shinyApp(ui = ui, server = server)
