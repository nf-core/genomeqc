#!/usr/bin/Rscript

# Written by Chris Wyatt and Fernando Duarte and released under the MIT license.
# Plots the phylogenetic tree with BUSCO, Quast and gene stats results

# Function to plot tree and plots
# Improved function to plot tree and plots
build_tree_plot <- function(tree, plots, legends, xlimit, right_margin = 15, bottom_margin = 60, tree_space_ratio = 1.3) {

  # Calculate tree width dynamically based on actual rendered plot
  tree_built <- ggplot_build(tree)
  tree_data <- tree_built$data[[1]]  # Get the tree data

  # Find the maximum x position and estimate label width
  max_x <- max(tree_data$x, na.rm = TRUE)

  # Get tip labels and find the longest one
  tip_labels <- tree$data$label[!is.na(tree$data$label)]
  max_label_chars <- max(nchar(tip_labels), na.rm = TRUE)

  # More balanced width calculation
  # Base it on the text size and character count, but be more conservative
  text_size_pts <- tree$theme$text$size %||% 11  # Default ggplot text size
  char_width_estimate <- text_size_pts * 0.015  # Reduced from 0.02 to give more space to tree
  label_padding <- max_label_chars * char_width_estimate

  # Set tree x-limit with balanced padding
  # Keep the tree structure prominent while ensuring labels fit
  tree_xlim <- max_x * tree_space_ratio + label_padding  # Adjustable tree space + label padding

  # Update tree with calculated xlim
  tree <- tree + xlim(0, tree_xlim)

  # Initialize combined plot with the tree
  combined_plots <- tree

  # Calculate widths based on number of plots
  n_plots <- length(plots)
  tree_width <- max(5, 10 - n_plots)  # Ensure minimum tree width
  plot_widths <- rep(1, n_plots)
  widths <- c(tree_width, plot_widths)

  # Add each additional plot
  for (plot in plots) {
    combined_plots <- combined_plots | plot
  }

  # Initialize combined legends with empty plot aligned with tree
  combined_legends <- plot_spacer() + xlimit

  # Add each additional legend
  for (legend in legends) {
    if (!is.null(legend) && length(legend) != 0) {
      legend_plot <- legend
    } else {
      legend_plot <- plot_spacer() + xlimit
    }
    combined_legends <- combined_legends | legend_plot
  }

  # Apply the layout widths
  combined_plots <- combined_plots + plot_layout(widths = widths)

  combined_legends <- combined_legends +
    plot_layout(widths = widths) +
    theme(plot.margin = margin(0, right_margin, bottom_margin, 0))

  # Combine plots and legends
  final_plot <- combined_plots / combined_legends +
    plot_layout(heights = c(0.99, 0.01))

  return(final_plot)
}

# Load libraries
library(ggtree)
library(ggplot2)
library(patchwork)
library(argparse)
library(dplyr)
library(tidyr)
library(scatterpie)
library(scales)

# Parse command-line arguments
parser <- ArgumentParser(description = 'Plot phylogenetic tree with statistics and true/false data')
parser$add_argument('tree_file', type = 'character', help = 'Path to the Newick formatted tree file')
parser$add_argument('busco_file', type = 'character', help = 'Path to processed BUSCO output file')
parser$add_argument('quast_file', type = 'character', help = 'Path to processed Quast output file')
parser$add_argument('genes_file', type = 'character', help = 'Path to gene stats output file')
parser$add_argument('nseqs_file', type = 'character', help = 'Path to number sequences with at least x number of complete BUSCOs file')
parser$add_argument('--ortho_file', type = 'character', default = NULL, help = 'Path to number of orthologous sequences file')
parser$add_argument('--text_size', type = 'double', default = 3, help = 'Text size for the tree plot')
parser$add_argument('--tree_scale', type = 'double', default = 0.0005, help = 'x axis limits scaling for tree plot (useful when tree labels appear truncated)')
parser$add_argument('--tree_margin', type = 'double', default = 15, help = "Tree's right margin size")
parser$add_argument('--bar_width', type = 'double', default = 0.7, help = 'Width of bar plots')
parser$add_argument('--rad_width', type = 'double', default = 0.4, help = 'Radius of pie charts')
parser$add_argument('--skip_stats', type = 'character', default = NULL, help = "Don't plot these stats (comma separated list)")
parser$add_argument('--type', type = 'character', choices = c('genome_only', 'genome_anno'), default = 'genome_anno', help = 'Select stats for genome only or for both genome and annotation')

args <- parser$parse_args()

# Avoid scientific notation in all plots
options(scipen = 999)

# Skipt these plots (parse skip arguments, thanks chat gpt)
skip <- if (!is.null(args$skip_stats)) strsplit(args$skip_stats, ",")[[1]] else character(0)

print(skip)

# Read the Newick tree from the file
tree <- read.tree(args$tree_file)

# Clean tree tip labels
tree$tip.label <- trimws(tree$tip.label)
#tree$tip.label <- tolower(tree$tip.label)

# If radious of pie charts is to big, it can
# mess the position of the pies, make them
# smaller
if (length(tree$tip.label) < 7) {
  args$bar_width <- args$bar_width/2
  args$rad_width <- args$rad_width/3
}

# Get order of tips (useful for data transformation of stats)
tree_plot <- ggtree(tree) # Temporary plot for get_taxa_name()
tips_order <- rev(get_taxa_name(tree_plot))

# --- Helper function to load BUSCO data ---
load_busco <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Read the data table from the file, ensuring species column is read as character
    # Load BUSCO
    data_busco <- read.csv(file, sep = "\t", colClasses = c("Input_file" = "character"))
    # Prepare BUSCO data (tidy)
    data_busco <- data_busco %>%
      # Remove extension from Input_file
      mutate(Input_file = tools::file_path_sans_ext(Input_file)) %>%
      # Rename 'Input_file' to 'species'
      rename(species = Input_file) %>%
    # Arrange data according to tree labels
      arrange(match(species, tree_tips)) %>%
    # Add node column
      mutate(node = 1:length(species)) # Node number needed for nodpie
  }, error = function(e) {
    warning("Failed to load BUSCO file: ", conditionMessage(e))
    NULL
  })
}

# --- Helper function to load Quast data ---
load_quast <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Load Quast
    data_quast <- read.csv(file, sep = "\t")
    # Change header of GC% and contigs column
    colnames(data_quast)[5] <- "GC"
    colnames(data_quast)[6] <- "Sequences"
    #Prepare Quast data (tidy)
    data_quast <- data_quast %>%
      # Remove the row where species is NA
      filter(!is.na(species)) %>%
      # Remove any remaining "bar" rows if necessary (check Chris script)
      filter(N50 != "bar") %>%
      # Total length values to Mb
      mutate(Total.length = (as.numeric(Total.length)/1000000)) %>%
      #Change Sequence values to integers
      mutate(Sequences = as.integer(Sequences)) %>%
      # Create new col with numbers of GC bp
      mutate(GC = as.numeric(Total.length)*as.numeric(GC)/100) %>%
      # Rename column to make it shorter
      rename(Length = Total.length)
      # Arrange data according to tree labels
      data_quast <- data_quast %>%
        arrange(match(species, tree_tips)) %>%
        mutate(node = 1:length(species))
    # Tidy Quast data
    # For N50/N90
    n5090 <- data_quast %>%
    # Convert wide to long format
      pivot_longer(cols = c(N50, N90),
                   names_to = "metric",
                   values_to = "value") %>%
      # Convert value column to numeric if needed
      mutate(value = as.numeric(value)) %>%
      mutate(value = (as.numeric(value)/1000000)) # Values in Mb
    # For GC content and length
    len <- data_quast %>%
      pivot_longer(cols = c(GC, Length),
                   names_to = "metric",
                   values_to = "value")
    list(full = data_quast, n5090 = n5090, len = len)
  }, error = function(e) {
    warning("Failed to load Quast file: ", conditionMessage(e))
    NULL
  })
}

# --- Helper function to load gene count data ---
load_genes <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Load gene stats
    data_genes <- read.csv(file, sep = "\t")
    # Prepare gene stats
    data_genes <- data_genes %>%
    # Rename columns
      rename(species = File) %>%
      rename(Total = Total_genes) %>%
      rename(Overlapping = Overlapping_genes) %>%
      # Remove ".counts.tsv"
      mutate(species = gsub("\\.counts\\.tsv", "", species))
    # Arrange data according to tree labels
    data_genes <- data_genes %>%
      arrange(match(species, tree_tips)) %>%
      mutate(node = 1:length(species))
    # Tidy gene stats data
    data_genes <- data_genes %>%
      pivot_longer(cols = c(Total, Overlapping),
                   names_to = "stat",
                   values_to = "value")
  }, error = function(e) {
    warning("Failed to load gene stats file: ", conditionMessage(e))
    NULL
  })
}

# --- Helper function to load gene count data ---
load_nseqs <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Load gene stats
    data_nseqs <- read.csv(file, sep = "\t")
    data_nseqs
    # Arrange data according to tree labels
    data_nseqs <- data_nseqs %>%
      arrange(match(.data[[names(.)[1]]], tree_tips)) %>%
      mutate(node = 1:n())
  }, error = function(e) {
    warning("Failed to load nseqs file: ", conditionMessage(e))
    NULL
  })
}

# --- Load optional input files ---
data_busco <- load_busco(args$busco_file, tips_order)
data_quast <- load_quast(args$quast_file, tips_order)
data_genes <- load_genes(args$genes_file, tips_order)
data_nseqs <- load_nseqs(args$nseqs_file, tips_order)
data_ortho <- load_nseqs(args$ortho_file, tips_order)

# Extract names for debugging
# tree_sp <- sort(tree$tip.label)
# quast_sp <- sort(unique(data_quast$species))
# busco_sp <- sort(data_busco$species)
# gene_sp <- sort(unique(data_genes$species))

# Debugging: Print species names from the tree and the data
# cat("Species names in the tree based on nw:\n")
# print(tree_sp )
# cat("\nSpecies names in data tables:\n")
# cat("BUSCO:", busco_sp, "\nQUAST:", quast_sp, "\ngene_stats:", gene_sp)

# Debugging: Check if there are any mismatches in species names
# datasets <- list(BUSCO = busco_sp, Quast = quast_sp, GeneStats = gene_sp)
# for (name in names(datasets)) {
#   if (any(tree_sp != datasets[[name]])) {
#    stop(paste("Species names in", name, "and tree labels do not match"))
#  }
#}


# Match names with new tree tips (only necessary for Quast)
# data_quast$species <- gsub("_", " ", data_quast$species)



# This is for the synteny paper, remove "_" and changes the first letter to uppercase
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Set standard theme for all barplots
barplots_theme <- theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 6),
    axis.ticks.y=element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_blank()
  )

# Ssequences with single copy orthologues plot
if (!is.null(data_nseqs)) {
# Plot number of chromosomes/sequences
  nseqs_plot <- ggplot(data_nseqs, aes(x=1, y=node)) +
    geom_text(aes(label = data_nseqs[,2])) +
    theme_void() +
    ggtitle("Seqs ≥5\nBUSCOs") +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))
}  else {
  nseqs_plot <- NULL
}

# Orthologous sequences
if (!is.null(data_ortho)) {
# Plot number of chromosomes/sequences
  ortho_plot <- ggplot(data_ortho, aes(x=1, y=node)) +
    geom_text(aes(label = data_ortho[,2])) +
    theme_void() +
    ggtitle("Ortho \nSeqs") +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))
}  else {
  ortho_plot <- NULL
}

# Quast plots
if (!is.null(data_quast)) {
# Plot number of chromosomes/sequences
  ch_plot <- ggplot(data_quast$full, aes(x=1, y=node)) +
    geom_text(aes(label = Sequences)) +
    theme_void() +
    ggtitle("Seq\nNumber") +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))
  # Plot Quast data genome size
  len_plot <- ggplot(
    data_quast$len,
    aes(y=value, x=node)
  ) +
    geom_col(
      aes(fill=metric),
      position = position_stack(reverse = TRUE),
      width = args$bar_width
    ) +
    scale_fill_manual(labels = c("GC %", "Length"), values = c("brown1", "cornflowerblue")) +
    ggtitle("Genome\nsize (Mb)") +
    barplots_theme +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -5)) +
    coord_flip() + #Flip plot
    xlab(NULL) +
    ylab(NULL)
  # Extract legend
  legend_len <- cowplot::get_legend(
    len_plot +
      theme(legend.position = "right",
            legend.justification = c(0, 1.2), # This is what actually move the legend, play with it, default position is c(1,0.5)
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.background = element_rect(fill = NA),
            legend.text = element_text(size = 8))
  )
  # Remove legend
  len_plot <- len_plot + guides(fill="none")
  # Prepare Quast data for plotting
  data_quast_n50 <- data_quast$n5090[data_quast$n5090$metric %in% "N50",]
  #data_quast_n90 <- data_quast_n5090[data_quast_n5090$metric %in% "N90",]
  # Plot Quast data N50
  n50_plot <- ggplot(
    data_quast_n50,
    aes(y=value, x=node)
  ) +
   geom_col(
      position = position_stack(reverse = TRUE), # For GC%
     width = args$bar_width,
     fill = "steelblue"
    ) +
    ggtitle("N50 (Mb)") +
    barplots_theme +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -0.4)) +
    coord_flip() +
    xlab(NULL) +
    ylab(NULL)
# Remove legend
n50_plot <- n50_plot + guides(fill="none")
}  else {
  n50_plot <- NULL
  len_plot <- NULL
  legend_len <- NULL
}

if (!is.null(data_busco)) {
  # Create the scatterpie plot
  pies_plot <- ggplot() +
    geom_scatterpie(
     aes(x = 0, y = node, group = species, r = args$rad_width),  # r determines the radius of the pies
      data = data_busco,
      cols = c("Single", "Duplicated", "Fragmented", "Missing"),
     color = NA
    ) +
    scale_fill_manual(values = c("deepskyblue", "orange", "darkorchid4", "firebrick1")) +
    coord_fixed() +
    theme_void() +
   ggtitle("BUSCO") +
   theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05))

  # Extract legend
  legend_busco <- cowplot::get_legend(
    pies_plot +
      #guides(fill=guide_legend(ncol=2)) +
     theme(legend.position = "right",
           legend.justification = c(0, 1.08),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, "cm"),
           #legend.background = element_rect(fill = NA), # I don't know why this doesn't work, if I set this to NA an outline appears around the legend
           legend.text = element_text(size = 8))
  )

  # Display the legend alone
  #cowplot::ggdraw() + cowplot::draw_grob(legend_busco)

  # Remove lenged for pieplot
  pies_plot <- pies_plot + guides(fill="none")
} else {
  pies_plot <- NULL
  legend_busco <- NULL
}

# Display the legend alone
#cowplot::ggdraw() + cowplot::draw_grob(legend_len)

if (!is.null(data_genes)) {
# Plot gene stats
  gene_plot <- ggplot(
    data_genes,
   aes(y=value, x=node)
  ) +
    geom_col(
      aes(fill=stat),
      position = position_stack(reverse = TRUE),
      width = args$bar_width
    ) +
    scale_fill_manual(values = c("indianred1", "lightsteelblue")) +
    ggtitle("Gene\nnumber") +
    barplots_theme +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -5)) +
    coord_flip() + #Flip plot
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) #+
  # Extract legend
  legend_gene <- cowplot::get_legend(
    gene_plot +
      theme(legend.position = "right",
            legend.justification = c(0, 1.2),
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.background = element_rect(fill = NA),
            legend.text = element_text(size = 8))
  )
  # Remove legend from pie plot
  gene_plot <- gene_plot + guides(fill="none")
} else {
  gene_plot <- NULL
  legend_gene <- NULL
}

# Display the legend alone
#cowplot::ggdraw() + cowplot::draw_grob(legend_gene)

# Helper function to safely extract axis ranges
get_plot_range <- function(plot, axis = "y") {
  # For error message in case plot is null
  plot_name <- deparse(substitute(plot))
  # To match tips with plots it's necessary to set the same ylim for all plots
  # Select the biggest range to avoid cropping
  tryCatch({
    built <- ggplot_build(plot)
    if (axis == "y") {
      return(built$layout$panel_scales_y[[1]]$range$range)
    } else if (axis == "x") {
      return(built$layout$panel_scales_x[[1]]$range$range)
    } else {
      stop("Invalid axis specified.")
   }
  }, error = function(e) {
    warning("Failed to load plot ", plot_name, ": ", conditionMessage(e))
    numeric(0)
  })
}

# Collect ranges safely
all_ranges <- c(
  get_plot_range(ch_plot, "y"),
  get_plot_range(nseqs_plot, "y"),
  get_plot_range(ortho_plot, "y"),
  get_plot_range(pies_plot, "y"),
  get_plot_range(len_plot, "x"),
  get_plot_range(n50_plot, "x"),
  get_plot_range(gene_plot, "x")
)

# Set new ylim based on the highest value taking into account both plots
# Add -0,1 and 0.5 to avoid the cropping of first and last pies
new_ylim <- ylim(c(min(all_ranges), max(all_ranges)))
# A new xlim is needed for barplots (equivalent to ylim), as these are flipped
# using coord_flip()
new_xlim <- xlim(c(min(all_ranges), max(all_ranges)))

# Set new ylim for sequnces plots (numbers)
if (!is.null(nseqs_plot)) nseqs_plot <- nseqs_plot + new_ylim
if (!is.null(ortho_plot)) ortho_plot <- ortho_plot + new_ylim
if (!is.null(ch_plot))    ch_plot   <- ch_plot + new_ylim

# Set new xlim for Quast genome size (equivalent to ylim)
if (!is.null(len_plot))  len_plot  <- len_plot + new_xlim

# Set new xlim for Quast N50 (equivalent to ylim)
if (!is.null(n50_plot))  n50_plot  <- n50_plot + new_xlim

# Set new ylim for Quast pies
if (!is.null(pies_plot)) pies_plot <- pies_plot + new_ylim

# Set new xlim for gene stats (equivalent to ylim)
if (!is.null(gene_plot)) gene_plot <- gene_plot + new_xlim

# Build tree
# Usage - replace your current tree building section with:
tree_plot <- ggtree(tree) +
  geom_tiplab(size = args$text_size, fontface = "italic", align = TRUE, hjust = -0.05) +
  theme(plot.margin = margin(10, 30, 10, 10)) +  # Increased right margin
  coord_cartesian(clip = "off") +
  new_ylim

# Set new ylim and xlim for tree
tree_plot <- tree_plot + new_ylim

# Set value for tree xlim to avoid the truncation of labels:
# Why "^2*0.001"? ^2 is because the relatin between number of characters and the number
# of pixels is close to beexponential, not proportional. 0.001 would be the length
# per character in the x axis scale. Script should allow to change this value
m = max(tree_plot$data$x) + max(nchar(tree_plot$data$label))^2*args$tree_scale

# Define named plot and legend lists (thanks to chat gpt)
all_plots <- list(
  ch_plot    = ch_plot,
  nseqs_plot = nseqs_plot,
  ortho_plot = ortho_plot,
  len_plot   = len_plot,
  gene_plot  = gene_plot,
  n50_plot   = n50_plot,
  pies_plot  = pies_plot
)

all_legends <- list(
  ch_plot    = NULL,
  nseqs_plot = NULL,
  ortho_plot = NULL,
  len_plot   = legend_len,
  gene_plot  = legend_gene,
  n50_plot   = NULL,
  pies_plot  = legend_busco
)

# Keep only plots and legends not in the skip list (thanks to chat gpt)
plots <- all_plots[!names(all_plots) %in% skip & !sapply(all_plots, is.null)]
legends <- all_legends[names(plots)]  # Re-align legends to plots

print("plots")
plots
print("legends")
legends

# Call the function
if (args$type == 'genome_anno') {
  final_plot <- build_tree_plot(
    tree = tree_plot,
    #n = m, # Only affects tree_plot
    plots = plots,
    legends = legends,
    new_xlim,
    15,
    60,
    args$tree_margin
  )
} else if (args$type == 'genome_only') {
  final_plot <- build_tree_plot(
    tree = tree_plot,
    #n = m, # Only affects tree_plot
    plots = plots,
    legends = legends,
    new_xlim,
    15,
    60,
    args$tree_margin
  )
}

pdf("tree_plot.pdf")
final_plot
dev.off()

svg("tree_plot.svg")
final_plot
dev.off()
