#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggtree)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scatterpie)
library(scales)
library(cowplot)
library(ape)

source("tree_functions.R")

# Load data once at startup
processed_data <- process_tree_data(
  tree_file = "tree.nw",
  busco_file = if(file.exists("Busco_combined.tsv")) "Busco_combined.tsv" else NULL,
  quast_file = if(file.exists("Quast_to_plot.tsv")) "Quast_to_plot.tsv" else NULL,
  genes_file = if(file.exists("gene_stats.tsv")) "gene_stats.tsv" else NULL,
  nseqs_file = if(file.exists("n_seqs_above_x_buscos_output.tsv")) "n_seqs_above_x_buscos_output.tsv" else NULL,
  ortho_file = if(file.exists("species_orthologous_chromosomes_output.tsv")) "species_orthologous_chromosomes_output.tsv" else NULL
)

# UI
ui <- fluidPage(
  titlePanel("Dynamic Tree Plot Adjuster"),

  tabsetPanel(
    tabPanel("Plot Controls",
             br(),

             # File Status
#             wellPanel(
#               h4("📁 Loaded Data Files"),
#               fluidRow(
#                 column(2, div(style = "text-align: center;", icon("tree", "fa-2x", style = "color: green;"), br(), "Tree")),
#                 column(2, div(style = "text-align: center;",
#                               if(file.exists("busco_file.tsv")) list(icon("check-circle", "fa-2x", style = "color: green;"), br(), "BUSCO")
#                               else list(icon("times-circle", "fa-2x", style = "color: red;"), br(), "No BUSCO"))),
#                 column(2, div(style = "text-align: center;",
#                               if(file.exists("quast_file.tsv")) list(icon("check-circle", "fa-2x", style = "color: green;"), br(), "Quast")
#                               else list(icon("times-circle", "fa-2x", style = "color: red;"), br(), "No Quast"))),
#                 column(2, div(style = "text-align: center;",
#                               if(file.exists("genes_file.tsv")) list(icon("check-circle", "fa-2x", style = "color: green;"), br(), "Genes")
#                               else list(icon("times-circle", "fa-2x", style = "color: red;"), br(), "No Genes"))),
#                 column(2, div(style = "text-align: center;",
#                               if(file.exists("nseqs_file.tsv")) list(icon("check-circle", "fa-2x", style = "color: green;"), br(), "NSeqs")
#                               else list(icon("times-circle", "fa-2x", style = "color: red;"), br(), "No NSeqs"))),
#                 column(2, div(style = "text-align: center;",
#                               if(file.exists("ortho_file.tsv")) list(icon("check-circle", "fa-2x", style = "color: green;"), br(), "Ortho")
#                               else list(icon("times-circle", "fa-2x", style = "color: red;"), br(), "No Ortho")))
#               )
#             ),

             # Controls
             fluidRow(
               column(4,
                      wellPanel(
                        h4("Margin Controls"),
                        sliderInput("tree_space_ratio", "Tree space ratio:", value = 1.3, min = -100, max = 100, step = 0.1),
                        sliderInput("top_margin", "Top margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                        sliderInput("right_margin", "Right margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                        sliderInput("bottom_margin", "Bottom Margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                        sliderInput("left_margin", "Left Margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                        sliderInput("tree_margin", "Tree Margin:", min = 0, max = 100, value = 15, step = 1)
                      )
               ),

               column(4,
                      wellPanel(
                        h4("Plot Parameters"),
                        numericInput("text_size", "Tip Text Size:", value = 3, min = 1, max = 10, step = 0.1),
                        numericInput("tree_scale", "Tree Scale:", value = 0.0005, min = 0.0001, max = 0.01, step = 0.0001),
                        #numericInput("tree_space_ratio", "Tree space ratio:", value = 1.3, min = -100, max = 100, step = 0.1),
                        numericInput("bar_width", "Bar Width:", value = 0.7, min = 0.1, max = 2, step = 0.1),
                        numericInput("rad_width", "Pie Radius:", value = 0.4, min = 0.1, max = 1, step = 0.05)
                      )
               ),

               column(4,
                      wellPanel(
                        h4("Display Options"),
                        checkboxGroupInput("skip_stats", "Skip Statistics:",
                                           choices = list(
                                             "Sequence Count" = "ch_plot",
                                             "NSeqs Plot" = "nseqs_plot",
                                             "Ortho Plot" = "ortho_plot",
                                             "Genome Length" = "len_plot",
                                             "Gene Statistics" = "gene_plot",
                                             "N50 Statistics" = "n50_plot",
                                             "BUSCO Pies" = "pies_plot"
                                           )),
                        #selectInput("plot_type", "Plot Type:",
                        #            choices = list("Genome + Annotation" = "genome_anno", "Genome Only" = "genome_only"),
                        #            selected = "genome_anno"),
                        sliderInput("plot_height", "Plot Height (px):", min = 400, max = 1200, value = 600, step = 50)
                      )
               )
             ),

             # Plot Display
             wellPanel(
               h4("Tree Plot Preview"),
               plotOutput("tree_plot", height = "auto"),
               br(),
               fluidRow(
                 column(4, actionButton("refresh_plot", "Refresh Plot", class = "btn-primary"), style = "float: center"),
               )
             )
    ),

    tabPanel("Export Settings",
             h4("Export preview"),
             uiOutput("export_preview"),
             br(),
             fluidRow(
               column(6,
                      wellPanel(
                        h4("Export Settings"),
                        numericInput("export_width", "Width (inches):", value = 12, min = 4, max = 20),
                        numericInput("export_height", "Height (inches):", value = 8, min = 4, max = 16),
                        numericInput("export_dpi", "DPI:", value = 300, min = 72, max = 600),
                        textInput("export_filename", "Filename:", value = "tree_plot")
                      )
               )
             ),
             fluidRow(
                column(6,
                      wellPanel(
                        fluidRow(
                          column(4, downloadButton("download_pdf", "Download PNG", class = "btn-success")),
                          column(4, downloadButton("download_svg", "Download SVG", class = "btn-info"), style = "float: right")
                        )
                      )
               )
             )
    )
  )
)



# Server
server <- function(input, output, session) {

  # Reactive values for the current plot
  current_plot <- reactiveVal(NULL)

  # Generate plot when parameters change or refresh button is clicked
  observeEvent(c(input$refresh_plot, input$text_size, input$tree_scale,
                 input$bar_width, input$rad_width, input$skip_stats,
                 input$tree_space_ratio, input$top_margin, input$right_margin,
                 input$bottom_margin, input$left_margin, input$tree_margin, input$skip_stats,
                 input$export_width, input$export_height, input$export_dpi), {

                   # Show loading notification
                   showNotification("Generating plot...", type = "message", duration = 2)

                   # Generate the plot
                   plot_result <- generate_complete_plot(
                     processed_data,
                     text_size = input$text_size,
                     tree_scale = input$tree_scale,
                     bar_width = input$bar_width,
                     rad_width = input$rad_width,
                     top_margin = input$top_margin,
                     right_margin = input$right_margin,
                     bottom_margin = input$bottom_margin,
                     left_margin = input$left_margin,
                     tree_margin = input$tree_margin,
                     skip_stats = input$skip_stats,
                     tree_space_ratio = input$tree_space_ratio

                   )

                   # Store the plot
                   current_plot(plot_result)

                   # Render the current plot in the UI
                   output$tree_plot <- renderPlot({
                     req(current_plot())  # Wait until current_plot() is not NULL
                     current_plot()
                   }, height = function() input$plot_height)

                  # Render plot to be exported
                  output$export_preview <- renderUI({
                    req(current_plot())  # Wait until current_plot() is not NULL
                    tmp <- tempfile(fileext = ".png")
                    # Save the plot as PNG using the export settings
                    ggsave(
                      filename = tmp,
                      plot = current_plot(),
                      device = "png",
                      width = input$export_width,
                      height = input$export_height,
                      units = "in",
                      dpi = input$export_dpi
                    )
                    # Convert to base64 to embed directly in Shiny
                    encoded <- base64enc::dataURI(file = tmp, mime = "image/png")

                    tags$img(src = encoded, style = "max-width:100%; height:auto; border:1px solid #ccc;")
                }
                )


                   # PNG download
                   output$download_pdf <- downloadHandler(
                     filename = function() {
                       paste0("tree_plot_", Sys.Date(), ".png")
                     },
                     content = function(file) {
                       req(current_plot())  # make sure plot exists
                       ggsave(
                         filename = file,
                         plot = current_plot(),
                         device = "png",
                         width = input$export_width, height = input$export_height,    # adjust or make reactive to input
                         units = "in",
                         dpi = input$export_dpi
                       )
                     }
                   )

                   # SVG download
                   output$download_svg <- downloadHandler(
                     filename = function() {
                       paste0("tree_plot_", Sys.Date(), ".svg")
                     },
                     content = function(file) {
                       req(current_plot())
                       ggsave(
                         filename = file,
                         plot = current_plot(),
                         device = "svg",
                         width = input$export_width, height = input$export_height,
                         units = "in",
                         dpi = input$export_dpi
                       )
                     }
                   )


                   #showNotification("Plot updated!", type = "success", duration = 1)
                 }, ignoreNULL = FALSE)


}

# Run the application
shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 8000, launch.browser = FALSE))
