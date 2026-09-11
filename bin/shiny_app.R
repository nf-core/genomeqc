# Written by Fernando Duarte and released under the MIT license.
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
library(ggtreeExtra)
library(ggnewscale)
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
  busco_file_geno = if(file.exists("Busco_combined_geno.tsv")) "Busco_combined_geno.tsv" else NULL,
  busco_file_prot = if(file.exists("Busco_combined_prot.tsv")) "Busco_combined_prot.tsv" else NULL,
  quast_file = if(file.exists("Quast_to_plot.tsv")) "Quast_to_plot.tsv" else NULL,
  genes_file = if(file.exists("gene_stats.tsv")) "gene_stats.tsv" else NULL,
  nseqs_file = if(file.exists("n_seqs_above_x_buscos_output.tsv")) "n_seqs_above_x_buscos_output.tsv" else NULL,
  ortho_file = if(file.exists("species_ortho_seq_count_output.tsv")) "species_ortho_seq_count_output.tsv" else NULL,
  te_file    = if(file.exists("te_table_output.tsv")) "te_table_output.tsv" else NULL,
  fcs_file   = if(file.exists("fcs_table_output.tsv")) "fcs_table_output.tsv" else NULL
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
                      # Margins / branch spacing apply to the non-circular layouts only
                      conditionalPanel(
                        condition = "input.tree_style != 'circular'",
                        wellPanel(
                          h4("Margin Controls"),
                          helpText("Applies to roundrect / ellipse / rectangular layouts."),
                          sliderInput("tree_space_ratio", "Branch length:", value = 1.3, min = -100, max = 100, step = 0.1),
                          sliderInput("top_margin", "Top margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                          sliderInput("right_margin", "Right margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                          sliderInput("bottom_margin", "Bottom Margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                          sliderInput("left_margin", "Left Margin:", min = 0, max = 100, value = 5.5, step = 0.5),
                          sliderInput("tree_margin", "Tree Margin:", min = 0, max = 100, value = 15, step = 1)
                        )
                      ),
                      # Ring options apply to the circular layout only
                      conditionalPanel(
                        condition = "input.tree_style == 'circular'",
                        wellPanel(
                          h4("Circular Ring Options"),
                          helpText("Rings show every statistic that is not skipped (right-hand panel)."),
                          sliderInput("ring_width", "Ring thickness:", value = 0.13, min = 0.05, max = 0.4, step = 0.01),
                          selectInput("quality_preset", "Quality thresholds (phylo group):",
                                      choices = c("Generic (lenient)" = "generic",
                                                  "Vertebrate" = "vertebrate",
                                                  "Insect" = "insect",
                                                  "Plant" = "plant",
                                                  "Fungi" = "fungi",
                                                  "Bacteria" = "bacteria",
                                                  "Custom..." = "custom"),
                                      selected = "generic"),
                          helpText("Quality rings are scored Good/Warn/Poor against these thresholds. N50 and sequence-count cut-offs are clade-dependent - pick the group matching your taxa, or set them yourself with 'Custom'."),
                          conditionalPanel(
                            condition = "input.quality_preset == 'custom'",
                            helpText("Direction is fixed per metric: BUSCO complete and N50 are higher-is-better; BUSCO duplicated and sequence count are lower-is-better."),
                            numericInput("thr_busco_good", "BUSCO complete % - Good at or above:", value = 95, min = 0, max = 100, step = 1),
                            numericInput("thr_busco_warn", "BUSCO complete % - Warn at or above:", value = 90, min = 0, max = 100, step = 1),
                            numericInput("thr_dup_good",   "BUSCO duplicated % - Good at or below:", value = 5, min = 0, max = 100, step = 1),
                            numericInput("thr_dup_warn",   "BUSCO duplicated % - Warn at or below:", value = 10, min = 0, max = 100, step = 1),
                            numericInput("thr_n50_good",   "N50 (bp) - Good at or above:", value = 1e6, min = 0, step = 1e5),
                            numericInput("thr_n50_warn",   "N50 (bp) - Warn at or above:", value = 1e5, min = 0, step = 1e4),
                            numericInput("thr_seq_good",   "Sequence count - Good at or below:", value = 1000, min = 0, step = 10),
                            numericInput("thr_seq_warn",   "Sequence count - Warn at or below:", value = 10000, min = 0, step = 100),
                            numericInput("thr_fcs_good",   "FCS non-contaminant % - Good at or above:", value = 99.5, min = 0, max = 100, step = 0.1),
                            numericInput("thr_fcs_warn",   "FCS non-contaminant % - Warn at or above:", value = 98, min = 0, max = 100, step = 0.1)
                          ),
                          checkboxInput("show_values", "Print values on rings", value = FALSE)
                        )
                      )
               ),

               column(4,
                      wellPanel(
                        h4("Plot Parameters"),
                        selectInput("tree_style", "Tree Style:",
                                    choices = c("Rounded (roundrect)" = "roundrect",
                                                "Curved (ellipse)" = "ellipse",
                                                "Rectangular (legacy)" = "rectangular",
                                                "Circular (rings)" = "circular"),
                                    selected = "roundrect"),
                        # Text size applies to all layouts (tip labels / circular tip numbers)
                        numericInput("text_size", "Tip Text Size:", value = 3, min = 1, max = 10, step = 0.1),
                        # These only affect the non-circular layouts
                        conditionalPanel(
                          condition = "input.tree_style != 'circular'",
                          numericInput("tree_scale", "Tree Scale:", value = 0.0005, min = 0.0001, max = 0.01, step = 0.0001),
                          numericInput("bar_width", "Bar Width:", value = 0.7, min = 0.1, max = 2, step = 0.1),
                          numericInput("rad_width", "Pie Radius:", value = 0.4, min = 0.1, max = 1, step = 0.05)
                        )
                      )
               ),

               column(4,
                      wellPanel(
                        h4("Display Options"),
                        helpText("Skip a statistic to drop its panel (non-circular) or ring (circular). Applies to all layouts."),
                        checkboxGroupInput("skip_stats", "Skip Statistics:",
                                           choices = list(
                                             "Sequence Count" = "ch_plot",
                                             "NSeqs" = "nseqs_plot",
                                             "Ortho Seqs" = "ortho_plot",
                                             "Genome Size" = "len_plot",
                                             "Gene Number" = "gene_plot",
                                             "N50" = "n50_plot",
                                             "BUSCO Genome" = "busco_gen_plot",
                                             "BUSCO Protein" = "busco_prot_plot",
                                             "TE Composition" = "te_plot",
                                             "FCS Contamination" = "fcs_plot"
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
                 input$bar_width, input$rad_width, input$skip_stats, input$tree_style,
                 input$ring_width, input$quality_preset, input$show_values,
                 input$thr_busco_good, input$thr_busco_warn, input$thr_dup_good, input$thr_dup_warn,
                 input$thr_n50_good, input$thr_n50_warn, input$thr_seq_good, input$thr_seq_warn,
                 input$thr_fcs_good, input$thr_fcs_warn,
                 input$tree_space_ratio, input$top_margin, input$right_margin,
                 input$bottom_margin, input$left_margin, input$tree_margin, input$skip_stats,
                 input$export_width, input$export_height, input$export_dpi), {

                   # Show loading notification
                   showNotification("Generating plot...", type = "message", duration = 2)

                   # "Custom" thresholds override the phylo-group preset. Directions are
                   # fixed per metric and must match those in QUALITY_PRESETS.
                   is_custom <- identical(input$quality_preset, "custom")
                   custom_thresholds <- if (is_custom) list(
                     busco_complete   = list(direction = "higher", good = input$thr_busco_good, warn = input$thr_busco_warn),
                     busco_duplicated = list(direction = "lower",  good = input$thr_dup_good,   warn = input$thr_dup_warn),
                     n50              = list(direction = "higher", good = input$thr_n50_good,   warn = input$thr_n50_warn),
                     seq_number       = list(direction = "lower",  good = input$thr_seq_good,   warn = input$thr_seq_warn),
                     fcs_noncontam    = list(direction = "higher", good = input$thr_fcs_good,   warn = input$thr_fcs_warn)
                   ) else NULL

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
                     tree_space_ratio = input$tree_space_ratio,
                     tree_style = input$tree_style,
                     ring_width = if (is.null(input$ring_width)) 0.13 else input$ring_width,
                     quality_preset = if (is_custom || is.null(input$quality_preset)) "generic" else input$quality_preset,
                     thresholds = custom_thresholds,
                     show_values = isTRUE(input$show_values)

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
                       # ggsave(..., device = "svg") requires the svglite package, which
                       # isn't installed in this app's container - use base R's svg()
                       # device directly instead, matching plot_tree_summary.R's own
                       # SVG export (no extra dependency, and DPI doesn't apply to a
                       # vector format so it's dropped here).
                       svg(filename = file, width = input$export_width, height = input$export_height)
                       on.exit(dev.off())
                       print(current_plot())
                     }
                   )


                   #showNotification("Plot updated!", type = "success", duration = 1)
                 }, ignoreNULL = FALSE)


}

# Run the application
shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 8000, launch.browser = FALSE))
