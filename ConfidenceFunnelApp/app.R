library(shiny)


# --- UI Definition ---
ui <- fluidPage(
  titlePanel("Pharmacology Confidence Funnel Tool"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      helpText("Ensure CSV columns are group names and first row is data."),
      
      checkboxInput("log10", "Apply Log10 Transformation", value = TRUE),
      hr(),
      
      textInput("main_title", "Main Title", "Confidence Funnel"),
      textInput("x_title", "X-Axis Label", "Value"),
      
      hr(),
      h4("Download Options"),
      # Low Quality
      downloadButton("downloadLow", "Download Notebook Plot (PNG)"),
      br(), br(),
      # High Quality
      downloadButton("downloadHigh", "Download Publication Plot (PDF)"),
      br(), br(),
      # CSV Export
      downloadButton("downloadData", "Export Funnel Coordinates (CSV)")
    ),
    
    mainPanel(
      plotOutput("funnelPlot", height = "600px"),
      # The Credits Section
      hr(),
      wellPanel(
        h4("Project Information"),
        p(em("An interactive utility for generating Confidence Funnels, providing high/low-resolution visualization downloads and confidence coordinate exports for publication-ready figures.")),
        p(strong("Manuscript:"), "A Robust Workflow for Exploratory Data Analysis and Outlier Management in Small-Sample Pharmacological Datasets"),
        p(strong("Author:"), "Michael J. Marino, PhD"),
        p(strong("Status:"), "Submitted to ", em("Biochemical Pharmacology"), "(2026)"),
        p(strong("License:"), "MIT License (c) 2026. This tool is open-source. Please cite the original manuscript in any resulting publications."),
        p(strong("Repository:"), a("GitHub: pharmacology-small-sample-EDA", 
                                   href="https://github.com/mjmarino-phd/pharmacology-small-sample-EDA", 
                                   target="_blank"))
      )
    )
  )
)

# --- Server Logic ---
server <- function(input, output) {

    
    # 1. Reactive data processing
    raw_data <- reactive({
      req(input$file)
      df <- read.csv(input$file$datapath)
      
      # Check for non-positive values if Log10 is selected
      if(input$log10) {
        # We check the whole dataframe for any values <= 0
        if(any(df <= 0, na.rm = TRUE)) {
          validate(
            need(FALSE, "Error: Data contains zero or negative values. Log10 transformation cannot be performed. Please verify your input data or uncheck the Log10 box.")
          )
        }
        df <- log10(df)
      }
      
      return(df)
    })
  
  # 2. Reusable Plotting Function
  render_funnel <- function() {
    data <- raw_data()
    groups <- names(data)
    n_groups <- length(groups)
    
    line_colors <- hcl.colors(n_groups, palette = "Set 2")
    fill_colors <- adjustcolor(line_colors, alpha.f = 0.15)
    point_colors <- adjustcolor(line_colors, alpha.f = 0.6)
    
    x_min <- min(unlist(data), na.rm = TRUE)
    x_max <- max(unlist(data), na.rm = TRUE)
    pad <- (x_max - x_min) * 0.1
    
    plot(NULL, xlim = c(x_min - pad, x_max + pad), ylim = c(50, 100),
         xlab = input$x_title, ylab = "Confidence Level (%)",
         main = input$main_title, frame.plot = FALSE)
    
    grid(nx = NULL, ny = NULL, col = "gray90", lty = "solid")
    
    for (i in 1:n_groups) {
      vals <- na.omit(data[[groups[i]]])
      n <- length(vals)
      if(n < 3) next
      
      m <- mean(vals)
      se <- sd(vals) / sqrt(n)
      conf_levels <- seq(0.50, 0.999, length.out = 100)
      
      t_crit <- qt(conf_levels + (1 - conf_levels) / 2, df = n - 1)
      lower_bound <- m - t_crit * se
      upper_bound <- m + t_crit * se
      
      polygon(c(lower_bound, rev(upper_bound)), 
              c(conf_levels * 100, rev(conf_levels * 100)), 
              col = fill_colors[i], border = NA)
      
      points(x = vals, y = rep(50, n), pch = 19, col = point_colors[i], cex = 1.3)
      abline(v = m, col = line_colors[i], lty = 2, lwd = 2)
      
      error_95 <- qt(0.975, df = n - 1) * se
      segments(m - error_95, 95, m + error_95, 95, col = line_colors[i], lwd = 2)
    }
    
    # Plot Legend - Positioned at Y=60 to clear the points at Y=50
    legend(x = x_max, y = 60, 
           xjust = 1, 
           legend = groups,
           fill = fill_colors, 
           border = line_colors, 
           bty = "n", 
           cex = 0.9, 
           y.intersp = 1.2)
    
  }
  
  # 3. Output: Screen Plot
  output$funnelPlot <- renderPlot({
    req(raw_data())
    render_funnel()
  })
  
  # 4. Download: Low Quality PNG
  output$downloadLow <- downloadHandler(
    filename = function() { paste("funnel_notebook_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600, res = 72) # Standard screen res
      render_funnel()
      dev.off()
    }
  )
  
  # 5. Download: High Quality PDF (Vector format is best for publication)
  output$downloadHigh <- downloadHandler(
    filename = function() { paste("funnel_pub_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 8, height = 6) # Vector format, infinite resolution
      render_funnel()
      dev.off()
    }
  )
  
  # 6. Export: CSV Coordinates
  output$downloadData <- downloadHandler(
    filename = function() { paste("funnel_coords_", Sys.Date(), ".csv", sep="") },
    content = function(file) {
      data <- raw_data()
      export_list <- list()
      
      for(name in names(data)) {
        vals <- na.omit(data[[name]])
        n <- length(vals)
        m <- mean(vals)
        se <- sd(vals) / sqrt(n)
        conf_levels <- seq(0.50, 0.99, by = 0.05) # Export key intervals
        t_crit <- qt(conf_levels + (1 - conf_levels) / 2, df = n - 1)
        
        temp_df <- data.frame(
          Group = name,
          Confidence = conf_levels * 100,
          Lower = m - t_crit * se,
          Upper = m + t_crit * se
        )
        export_list[[name]] <- temp_df
      }
      write.csv(do.call(rbind, export_list), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
