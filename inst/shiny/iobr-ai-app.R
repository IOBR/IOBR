# IOBR AI Assistant Shiny App
#
# This is a Shiny web application for the IOBR AI assistant.
# It provides a user interface to:
# - Configure AI providers (OpenAI, Hugging Face, Anthropic, Custom, or Dummy)
# - Create embedding index from package documentation
# - Query the assistant and view responses with sources
# - Download/upload index files

library(shiny)

# UI Definition
ui <- fluidPage(
  # Custom CSS
  tags$head(
    tags$style(HTML("
      .main-title {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        border-radius: 8px;
        margin-bottom: 20px;
      }
      .code-block {
        background-color: #f5f5f5;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 10px;
        margin: 10px 0;
        font-family: monospace;
        overflow-x: auto;
      }
      .source-item {
        background-color: #f9f9f9;
        border-left: 4px solid #667eea;
        padding: 10px;
        margin: 5px 0;
        border-radius: 4px;
      }
      .score-badge {
        display: inline-block;
        background-color: #667eea;
        color: white;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.85em;
        margin-left: 10px;
      }
      .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        border-radius: 4px;
        padding: 10px;
        margin: 10px 0;
      }
      .success-box {
        background-color: #d4edda;
        border: 1px solid #28a745;
        border-radius: 4px;
        padding: 10px;
        margin: 10px 0;
      }
    "))
  ),
  
  # Title
  div(class = "main-title",
    h1("IOBR AI Assistant", style = "margin: 0;"),
    p("Intelligent documentation search and Q&A powered by RAG", 
      style = "margin: 5px 0 0 0;")
  ),
  
  # Main layout
  sidebarLayout(
    # Sidebar
    sidebarPanel(
      width = 3,
      
      h4("Provider Configuration"),
      
      selectInput("provider_name", "Provider:",
                  choices = c("Dummy (Offline)" = "dummy",
                             "OpenAI" = "openai",
                             "Hugging Face" = "huggingface",
                             "Anthropic" = "anthropic",
                             "Custom" = "custom"),
                  selected = "dummy"),
      
      conditionalPanel(
        condition = "input.provider_name != 'dummy'",
        textInput("api_key", "API Key:", value = "", 
                  placeholder = "Enter your API key"),
        textInput("base_url", "Base URL (optional):", value = "",
                  placeholder = "Leave empty for defaults"),
        textInput("model_embeddings", "Embeddings Model (optional):", value = "",
                  placeholder = "Leave empty for defaults"),
        textInput("model_chat", "Chat Model (optional):", value = "",
                  placeholder = "Leave empty for defaults")
      ),
      
      hr(),
      
      h4("Index Management"),
      
      actionButton("create_index_btn", "Create Index", 
                   class = "btn-primary", width = "100%"),
      
      br(), br(),
      
      actionButton("list_index_btn", "View Index Info", 
                   class = "btn-info", width = "100%"),
      
      br(), br(),
      
      downloadButton("download_index_btn", "Download Index", 
                     class = "btn-success", style = "width: 100%;"),
      
      br(), br(),
      
      fileInput("upload_index", "Upload Index", 
                accept = ".rds", width = "100%"),
      
      br(),
      
      actionButton("reset_index_btn", "Reset Index", 
                   class = "btn-danger", width = "100%")
    ),
    
    # Main panel
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        # Query Tab
        tabPanel(
          "Query Assistant",
          br(),
          
          textAreaInput("query_input", "Ask a question:", 
                       value = "", 
                       placeholder = "E.g., How do I calculate TME scores using IOBR?",
                       height = "100px",
                       width = "100%"),
          
          fluidRow(
            column(6, 
              sliderInput("top_k", "Number of sources:", 
                         min = 1, max = 10, value = 5, step = 1)
            ),
            column(6,
              sliderInput("temperature", "Temperature:", 
                         min = 0, max = 1, value = 0.2, step = 0.1)
            )
          ),
          
          actionButton("query_btn", "Ask Assistant", 
                       class = "btn-primary btn-lg", width = "200px"),
          
          hr(),
          
          # Response area
          uiOutput("response_output"),
          
          # Code block area
          conditionalPanel(
            condition = "output.has_code",
            h4("R Code Examples"),
            div(class = "warning-box",
              strong("⚠️ Code Execution:"),
              " Review the code before running. The 'Run in Sandbox' button is disabled for safety."
            ),
            uiOutput("code_output"),
            actionButton("run_sandbox_btn", "Run in Sandbox", 
                        class = "btn-warning", disabled = TRUE),
            helpText("Sandbox execution is disabled. Copy and review code manually before running.")
          ),
          
          # Sources area
          conditionalPanel(
            condition = "output.has_sources",
            h4("Retrieved Sources"),
            uiOutput("sources_output")
          )
        ),
        
        # Info Tab
        tabPanel(
          "Information",
          br(),
          
          h3("About IOBR AI Assistant"),
          
          p("This AI assistant helps you explore and use the IOBR package for immune oncology biological research. It uses Retrieval-Augmented Generation (RAG) to provide accurate, context-aware answers based on the package documentation."),
          
          h4("Features:"),
          tags$ul(
            tags$li("Search package documentation using semantic similarity"),
            tags$li("Get AI-generated answers with source citations"),
            tags$li("View and copy example R code"),
            tags$li("Support for multiple AI providers (OpenAI, Hugging Face, etc.)"),
            tags$li("Offline dummy mode for testing without API keys")
          ),
          
          h4("Getting Started:"),
          tags$ol(
            tags$li("Choose a provider (start with 'Dummy' for testing)"),
            tags$li("Click 'Create Index' to build the search index"),
            tags$li("Enter your question in the Query tab"),
            tags$li("Click 'Ask Assistant' to get answers")
          ),
          
          h4("Using Real AI Providers:"),
          
          div(class = "warning-box",
            strong("Security Notice:"),
            " Never commit API keys to version control. Use environment variables or secure configuration."
          ),
          
          tags$b("OpenAI:"),
          tags$pre("provider <- list(\n  name = 'openai',\n  api_key = Sys.getenv('OPENAI_API_KEY')\n)"),
          
          tags$b("Hugging Face:"),
          tags$pre("provider <- list(\n  name = 'huggingface',\n  api_key = Sys.getenv('HUGGINGFACE_API_KEY')\n)"),
          
          h4("Programmatic Usage:"),
          
          p("You can also use the assistant programmatically in R:"),
          
          tags$pre(
"# Initialize
iobr_ai_init()

# Query
response <- iobr_ai_query('How do I use CIBERSORT?')
cat(response$answer)

# View sources
print(response$retrieved)
"
          ),
          
          hr(),
          
          p(em("For more information, see inst/ai/README.md in the package."))
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values
  rv <- reactiveValues(
    index_created = FALSE,
    last_response = NULL,
    status_message = NULL
  )
  
  # Get current provider configuration
  get_provider <- reactive({
    if (input$provider_name == "dummy") {
      return(list(name = "dummy"))
    }
    
    provider <- list(name = input$provider_name)
    
    if (input$api_key != "") {
      provider$api_key <- input$api_key
    }
    
    if (input$base_url != "") {
      provider$base_url <- input$base_url
    }
    
    if (input$model_embeddings != "") {
      provider$model_embeddings <- input$model_embeddings
    }
    
    if (input$model_chat != "") {
      provider$model_chat <- input$model_chat
    }
    
    return(provider)
  })
  
  # Create index
  observeEvent(input$create_index_btn, {
    showModal(modalDialog(
      title = "Creating Index",
      "Please wait while the index is being created...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      provider <- get_provider()
      
      # Try to find package root
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        # Try parent directories
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      
      metadata <- iobr_ai_init(pkg_path = pkg_path, provider = provider)
      
      rv$index_created <- TRUE
      rv$status_message <- sprintf(
        "Index created successfully! %d chunks from %d documents.",
        metadata$num_chunks, metadata$num_docs
      )
      
      removeModal()
      showModal(modalDialog(
        title = "Success",
        div(class = "success-box", rv$status_message),
        easyClose = TRUE
      ))
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Failed to create index:", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
  
  # Query assistant
  observeEvent(input$query_btn, {
    req(input$query_input)
    
    if (input$query_input == "") {
      showModal(modalDialog(
        title = "Empty Query",
        "Please enter a question.",
        easyClose = TRUE
      ))
      return()
    }
    
    showModal(modalDialog(
      title = "Processing Query",
      "Searching documentation and generating response...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      provider <- get_provider()
      
      # Find index path
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      index_path <- file.path(pkg_path, "inst/ai/iobr_embeddings.rds")
      
      response <- iobr_ai_query(
        query = input$query_input,
        index_path = index_path,
        provider = provider,
        top_k = input$top_k,
        temperature = input$temperature
      )
      
      rv$last_response <- response
      
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Failed to process query:", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
  
  # List index info
  observeEvent(input$list_index_btn, {
    tryCatch({
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      index_path <- file.path(pkg_path, "inst/ai/iobr_embeddings.rds")
      
      info <- iobr_ai_list_index(index_path)
      
      if (!is.null(info)) {
        showModal(modalDialog(
          title = "Index Information",
          div(
            p(strong("Path:"), info$index_path),
            p(strong("Entries:"), info$num_entries),
            p(strong("Created:"), format(info$created_at)),
            p(strong("Provider:"), info$provider),
            p(strong("File size:"), sprintf("%.2f MB", info$file_size / 1024^2))
          ),
          easyClose = TRUE
        ))
      } else {
        showModal(modalDialog(
          title = "No Index",
          "Index not found. Please create one first.",
          easyClose = TRUE
        ))
      }
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Failed to read index info:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Reset index
  observeEvent(input$reset_index_btn, {
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
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      index_path <- file.path(pkg_path, "inst/ai/iobr_embeddings.rds")
      
      success <- iobr_ai_reset_index(index_path)
      
      removeModal()
      
      if (success) {
        rv$index_created <- FALSE
        rv$last_response <- NULL
        
        showModal(modalDialog(
          title = "Success",
          "Index has been deleted.",
          easyClose = TRUE
        ))
      }
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Failed to reset index:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Handle index upload
  observeEvent(input$upload_index, {
    req(input$upload_index)
    
    tryCatch({
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      index_dir <- file.path(pkg_path, "inst/ai")
      
      if (!dir.exists(index_dir)) {
        dir.create(index_dir, recursive = TRUE)
      }
      
      index_path <- file.path(index_dir, "iobr_embeddings.rds")
      
      file.copy(input$upload_index$datapath, index_path, overwrite = TRUE)
      
      rv$index_created <- TRUE
      
      showModal(modalDialog(
        title = "Success",
        "Index uploaded successfully!",
        easyClose = TRUE
      ))
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Failed to upload index:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Download index handler
  output$download_index_btn <- downloadHandler(
    filename = function() {
      paste0("iobr_embeddings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      pkg_path <- getwd()
      if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
        for (i in 1:3) {
          pkg_path <- dirname(pkg_path)
          if (file.exists(file.path(pkg_path, "DESCRIPTION"))) break
        }
      }
      index_path <- file.path(pkg_path, "inst/ai/iobr_embeddings.rds")
      
      if (file.exists(index_path)) {
        file.copy(index_path, file)
      } else {
        stop("Index file not found")
      }
    }
  )
  
  # Render response
  output$response_output <- renderUI({
    req(rv$last_response)
    
    div(
      h4("Answer:"),
      div(
        style = "background-color: white; border: 1px solid #ddd; border-radius: 4px; padding: 15px; margin: 10px 0;",
        HTML(gsub("\n", "<br>", rv$last_response$answer))
      )
    )
  })
  
  # Check if response has code
  output$has_code <- reactive({
    !is.null(rv$last_response) && !is.null(rv$last_response$code)
  })
  outputOptions(output, "has_code", suspendWhenHidden = FALSE)
  
  # Render code blocks
  output$code_output <- renderUI({
    req(rv$last_response$code)
    
    code_blocks <- lapply(seq_along(rv$last_response$code), function(i) {
      div(
        class = "code-block",
        tags$pre(rv$last_response$code[[i]])
      )
    })
    
    do.call(tagList, code_blocks)
  })
  
  # Check if response has sources
  output$has_sources <- reactive({
    !is.null(rv$last_response) && !is.null(rv$last_response$retrieved)
  })
  outputOptions(output, "has_sources", suspendWhenHidden = FALSE)
  
  # Render sources
  output$sources_output <- renderUI({
    req(rv$last_response$retrieved)
    
    source_items <- lapply(rv$last_response$retrieved, function(src) {
      div(
        class = "source-item",
        p(
          strong(src$source),
          span(class = "score-badge", sprintf("%.3f", src$score))
        ),
        p(style = "margin-top: 5px; color: #666;", 
          substr(src$text, 1, 300),
          if (nchar(src$text) > 300) "..." else ""
        )
      )
    })
    
    do.call(tagList, source_items)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
