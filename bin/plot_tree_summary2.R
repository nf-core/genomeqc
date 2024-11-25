#!/usr/bin/Rscript

# Written by Chris Wyatt and released under the MIT license. 
# Plots the phylogenetic tree with BUSCO result in pie charts

# Load necessary libraries
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

library(ggtree)
library(ggplot2)
library(cowplot)
library(argparse)
library(dplyr)
library(tidyr)

# Function to extract legend from a ggplot
extract_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  return(legend)
}

# Parse command-line arguments
parser <- ArgumentParser(description = 'Plot phylogenetic tree with statistics and true/false data')
parser$add_argument('tree_file', type = 'character', help = 'Path to the Newick formatted tree file')
parser$add_argument('data_file', type = 'character', help = 'Path to the CSV file with statistic and true/false data')
parser$add_argument('--text_size', type = 'double', default = 3, help = 'Text size for the plots')
parser$add_argument('--tree_size', type = 'double', default = 0.3, help = 'Proportion of the plot width for the tree')
args <- parser$parse_args()

# Read the Newick tree from the file
tree <- read.tree(args$tree_file)

# Clean tree tip labels
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- tolower(tree$tip.label)

# Read the data table from the file, ensuring species column is read as character
data <- read.csv(args$data_file, sep = "\t", colClasses = c("species" = "character"))

# Clean data species names
data$species <- trimws(data$species)
data$species <- tolower(data$species)

# Extract the column headers and the second row to determine plot types
plot_types <- as.character(data[1, ])
data <- data[-1, ]  # Remove the second row used for plot types

# Debugging: Print species names from the tree and the data
cat("Species names in the tree based on nw:\n")
print(tree$tip.label)
cat("\nSpecies names in the data table before processing:\n")
print(data$species)

# Debugging: Print the data frame to check its structure
cat("\nData frame:\n")
print(data)

# Ensure species names match the tree tips
missing_in_tree <- setdiff(data$species, tree$tip.label)
missing_in_data <- setdiff(tree$tip.label, data$species)

if (length(missing_in_tree) > 0) {
  cat("Species in data but not in tree:\n")
  print(missing_in_tree)
}

if (length(missing_in_data) > 0) {
  cat("Species in tree but not in data:\n")
  print(missing_in_data)
}

if (length(missing_in_tree) > 0 || length(missing_in_data) > 0) {
  stop("Species names in the data do not match the tree tips")
}

# Plot the phylogenetic tree
tree_plot <- ggtree(tree) + 
  geom_tiplab() + 
  coord_cartesian(clip="off") +
  ggtitle("Phylogenetic Tree") +
  theme(plot.margin = margin(10, 70, 10, 10))  # Increase margins

pdf ("Tree_only.pdf")
tree_plot
dev.off()


# Extract the exact tree tip names:axis.text
tree$tip.label <- get_taxa_name(tree_plot)

# Reorder the data based on the order of species in the tree
data <- data %>% arrange(match(species, tree$tip.label))

# Debugging: Print plot types
cat("\nPlot types:\n")
print(plot_types)

# Debugging: Print the reordered data frame to check its structure
cat("\nReordered data frame:\n")
print(data)

write.table(data, "Reordered_output_tree.tsv", sep="\t", quote=F)

# Check if the number of unique species matches the number of tree tips
if (length(unique(data$species)) != length(tree$tip.label)) {
  warning("The number of unique species in the data does not match the number of tree tips.")
}

# Find column headers
column_headers <- colnames(data)

# Create plots based on the plot types
plots <- list()
legend_plot <- NULL
for (i in 2:length(column_headers)) {
  column_name <- column_headers[i]
  plot_type <- plot_types[i]
  
  if (plot_type == "bar") {
    bar_plot <- ggplot(data, aes(x = as.numeric(!!sym(column_name)), y = factor(species, levels = rev(tree$tip.label)))) + 
      geom_bar(stat = "identity", fill = "darkblue") + 
      geom_text(aes(label = !!sym(column_name)), hjust = -0.1, size = args$text_size) +  # Add text labels
      theme_minimal() + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
            panel.grid.minor.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0)) + 
      labs(x = column_name) +
      ggtitle(column_name)
    plots[[length(plots) + 1]] <- bar_plot
  } else if (plot_type == "text") {
    text_plot <- ggplot(data, aes(y = factor(species, levels = rev(tree$tip.label)), x = 1, label = !!sym(column_name))) + 
      geom_text(size = args$text_size, hjust = 0) + 
      theme_void() + 
      labs(x = NULL, y = NULL) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)) +
      ggtitle(column_name)
    plots[[length(plots) + 1]] <- text_plot
  } else if (plot_type == "stacked") {
    # Split the stacked values into separate columns
    stacked_data <- data %>%
      separate(column_name, into = paste0(column_name, "_", 1:4), sep = ",", convert = TRUE, extra = "drop") %>%
      pivot_longer(cols = starts_with(column_name), names_to = "stack", values_to = "value") %>%
      mutate(stack = factor(stack, levels = paste0(column_name, "_", 1:4)))
    
    stacked_plot <- ggplot(stacked_data, aes(x = value, y = factor(species, levels = rev(tree$tip.label)), fill = stack)) + 
      geom_bar(stat = "identity", position = "fill") + 
      theme_minimal() + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
            panel.grid.minor.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0)) + 
      labs(x = column_name) +
      ggtitle(column_name)
    
    # Extract the legend from the stacked plot
    legend_plot <- extract_legend(stacked_plot)
    
    # Remove the legend from the stacked plot
    stacked_plot <- stacked_plot + theme(legend.position = "none")
    
    plots[[length(plots) + 1]] <- stacked_plot
  } else if (plot_type == "pie") {
    # Create a list to hold the pie charts for each row
    pie_plots <- list()
    
    for (j in 1:nrow(data)) {
      # Extract the values for the current row (species)
      pie_data <- data[j, ] %>%
        separate(column_name, into = paste0(column_name, "_", 1:4), sep = ",", convert = TRUE, extra = "drop") %>%
        pivot_longer(cols = starts_with(column_name), names_to = "slice", values_to = "value") %>%
        mutate(slice = factor(slice, levels = paste0(column_name, "_", 1:4)))

      # Create the pie plot for this row (species)
      pie_plot <- ggplot(pie_data, aes(x = "", y = value, fill = slice)) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_polar(theta = "y") +  # Convert the bar plot to a pie chart
        theme_void() +  # Remove axes and labels
        theme(legend.position = "none")  # Remove the legend

      # Save the pie chart for this species
      pie_plots[[length(pie_plots) + 1]] <- pie_plot
    }
    
    # Combine pie charts in a vertical grid (aligned with tree tips)
    combined_pie_plot <- plot_grid(plotlist = pie_plots, ncol = 1, align = "v")
    plots[[length(plots) + 1]] <- combined_pie_plot
  } else {
    cat("Unknown plot type:", plot_type, "for column:", column_name, "\n")
  }
}


# Combine the tree and the plots
if (length(plots) > 0) {
  combined_plot <- plot_grid(tree_plot, plot_grid(plotlist = plots, ncol = 1), ncol = 2, rel_widths = c(args$tree_size, 1 - args$tree_size))
  
  # Save the combined plot to a PDF file
  ggsave("Phyloplot_busco.pdf", plot = combined_plot, width = 10, height = 8)
  
  # Save the legend to a separate PDF file if present
  if (!is.null(legend_plot)) {
    legend_plot <- cowplot::plot_grid(legend_plot)
    ggsave("Legend.pdf", plot = legend_plot, width = 10, height = 2)
  }
}




# Debugging: Print the list of plots
cat("List of plots:\n")
print(plots)

# Combine the plots if there are any
if (length(plots) > 0) {
  combined_plot <- plot_grid(tree_plot, plot_grid(plotlist = plots, ncol = length(plots)), ncol = 2, rel_widths = c(args$tree_size, 1 - args$tree_size))  # Adjust widths based on tree_size

  # Save the combined plot to a PDF file
  ggsave("Phyloplot.pdf", plot = combined_plot, width = 10, height = 8)
  
  # Save the legend to a separate PDF file
  if (!is.null(legend_plot)) {
    legend_plot <- cowplot::plot_grid(legend_plot)
    ggsave("Legend.pdf", plot = legend_plot, width = 10, height = 2)
  }
} else {
  cat("No plots to combine.\n")
}

# Print warnings
warnings()