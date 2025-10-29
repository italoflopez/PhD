library(shiny)
library(ggplot2)
library(ggforce)
library(plotly)
library(dplyr)
library(tidyr)

data_long_all <- rbind(data_long, data_long_stock_and_watson, data_long_new_fast_variables)

# Parameters
ncol_page <- 6
nrow_page <- 6
plots_per_page <- ncol_page * nrow_page
total_pages <- ceiling(length(unique(data_long_all$variable)) / plots_per_page)

ui <- fluidPage(
  titlePanel("Time Series Carousel"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("prev_page", "Previous Page"),
      textOutput("page_info"),
      actionButton("next_page", "Next Page"),
      br(), br(),
      selectInput(
        "chosen_var",
        "Select a series to view:",
        choices = unique(data_long_all$variable),
        selected = unique(data_long_all$variable)[1]
      ),
      br(),
      textOutput("selected_var")
    ),
    mainPanel(
      plotlyOutput("facet_plot", height = "800px"),
      br(),
      plotlyOutput("single_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  current_page <- reactiveVal(1)
  selected_var <- reactiveVal(NULL)
  
  # Update page info
  output$page_info <- renderText({
    paste("Page", current_page(), "of", total_pages)
  })
  
  # Navigation buttons
  observeEvent(input$next_page, {
    if (current_page() < total_pages)
      current_page(current_page() + 1)
  })
  observeEvent(input$prev_page, {
    if (current_page() > 1)
      current_page(current_page() - 1)
  })
  
  # Faceted plot (paginated)
  output$facet_plot <- renderPlotly({
    facet_page <- current_page()
    
    # Get list of variable names for the current page
    vars <- unique(data_long_all$variable)
    start <- (facet_page - 1) * plots_per_page + 1
    end <- min(facet_page * plots_per_page, length(vars))
    vars_page <- vars[start:end]
    
    # Subset data to only those variables
    data_page <- data_long_all %>% filter(variable %in% vars_page)
    
    facet_plot <- ggplot(data_page, aes(x = date, y = value)) +
      geom_line(color = "steelblue") +
      facet_wrap(~ variable, ncol = ncol_page, scales = "free_y") +
      theme_minimal(base_size = 8)
    
    ggplotly(facet_plot, tooltip = "variable", source = "facet_plot") %>%
      layout(dragmode = FALSE)
  })
  
  # Handle clicks on facet plots
  observeEvent(event_data("plotly_click", source = "facet_plot"), {
    d <- event_data("plotly_click", source = "facet_plot")
    if (!is.null(d)) {
      clicked_x <- as.Date(d$x)
      clicked_y <- d$y
      var_clicked <- data_long_all %>%
        mutate(dist = abs(as.numeric(date - clicked_x)) + abs(value - clicked_y)) %>%
        arrange(dist) %>%
        slice(1) %>%
        pull(variable)
      
      selected_var(var_clicked)
      updateSelectInput(session, "chosen_var", selected = var_clicked)  # Sync dropdown
    }
  })
  
  # Handle dropdown selection
  observeEvent(input$chosen_var, {
    selected_var(input$chosen_var)
  })
  
  # Display selected variable name
  output$selected_var <- renderText({
    if (is.null(selected_var()))
      "Click on a plot or select a series from the dropdown to see it in detail"
    else
      paste("Selected variable:", selected_var())
  })
  
  # Detailed plot for selected variable
  output$single_plot <- renderPlotly({
    req(selected_var())
    ggplotly(
      ggplot(filter(data_long_all, variable == selected_var()),
             aes(x = date, y = value)) +
        geom_line(color = "tomato") +
        ggtitle(selected_var()) +
        theme_minimal(base_size = 12)
    )
  })
}

shinyApp(ui, server)
