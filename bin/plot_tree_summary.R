#!/usr/bin/env Rscript

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
library(ggtreeExtra)
library(ggnewscale)
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
parser$add_argument('quast_file', type = 'character', help = 'Path to processed Quast output file')
parser$add_argument('genes_file', type = 'character', help = 'Path to gene stats output file')
parser$add_argument('nseqs_file', type = 'character', help = 'Path to number sequences with at least x number of complete BUSCOs file')
parser$add_argument('--busco_geno', type = 'character', help = 'Path to processed BUSCO genome output file')
parser$add_argument('--busco_prot', type = 'character', help = 'Path to processed BUSCO protein output file')
parser$add_argument('--ortho_file', type = 'character', default = NULL, help = 'Path to number of orthologous sequences file')
parser$add_argument('--te_file', type = 'character', help = 'Path to Transposable Elementes stats table')
parser$add_argument('--fcs_file', type = 'character', help = 'Path to FCS-GX contamination stats table')
parser$add_argument('--text_size', type = 'double', default = 3, help = 'Text size for the tree plot')
parser$add_argument('--tree_scale', type = 'double', default = 0.0005, help = 'x axis limits scaling for tree plot (useful when tree labels appear truncated)')
parser$add_argument('--tree_margin', type = 'double', default = 15, help = "Tree's right margin size")
parser$add_argument('--bar_width', type = 'double', default = 0.7, help = 'Width of bar plots')
parser$add_argument('--rad_width', type = 'double', default = 0.4, help = 'Radius of pie charts')
parser$add_argument('--skip_stats', type = 'character', default = NULL, help = "Don't plot these stats (comma separated list)")
parser$add_argument('--type', type = 'character', choices = c('genome_only', 'genome_anno'), default = 'genome_anno', help = 'Select stats for genome only or for both genome and annotation')
parser$add_argument('--tree_style', type = 'character', choices = c('roundrect', 'ellipse', 'rectangular', 'circular'), default = 'roundrect', help = 'Tree layout style: roundrect (rounded branches, default), ellipse (curved branches with node points), rectangular (legacy look with dotted leader lines), or circular (fan tree with the stats drawn as concentric coloured rings)')
parser$add_argument('--circular_rings', type = 'character', default = 'ch_plot,n50_plot,busco_gen_plot,busco_prot_plot,fcs_plot,len_plot,gene_plot', help = "Circular layout only: comma-separated, ordered (inner->outer) list of stats to draw as rings, or 'all' to show every available stat. Descriptive rings are always grouped nearest the tree and quality (traffic-light) rings on the outer rim, so the two stay visually separate. Keys: ch_plot, len_plot, n50_plot, gene_plot, busco_gen_plot, busco_prot_plot, busco_dup_plot, fcs_plot, nseqs_plot, ortho_plot")
parser$add_argument('--quality_preset', type = 'character', choices = c('generic', 'vertebrate', 'insect', 'plant', 'fungi', 'bacteria'), default = 'generic', help = "Circular layout only: phylogenetic-group thresholds used to score the quality rings (traffic light). 'generic' is deliberately lenient; pick the group matching your taxa. The N50/sequence-count cut-offs in particular are starting points and should be tuned per project.")
parser$add_argument('--quality_thresholds', type = 'character', default = NULL, help = "Circular layout only: override individual --quality_preset cut-offs, as comma-separated 'metric=good:warn' pairs. Metrics: busco_complete, busco_duplicated, n50, seq_number, fcs_noncontam (direction is fixed per metric, so only the cut-offs are given). Example: 'n50=2e6:5e5,seq_number=500:5000'. Unmentioned metrics keep the preset's cut-offs.")
parser$add_argument('--show_ring_values', action = 'store_true', help = 'Circular layout only: print each value on its ring (redundant encoding, so the figure is not colour-only). Best for small trees.')
parser$add_argument('--len_pos_x', type = 'double', default = 0.5, help = 'Position (legend.justification x-anchor) of the BUSCO pie and TE bar legends')

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

# --- Helper function to load BUSCO n seqs data ---
load_nseqs <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Load n seqs
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

# --- Helper function to load TE data ---
load_te <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    # Load gene stats
    data_te <- read.csv(file, sep = "\t")
    # node = this species' real index in tree_tips, not its row position -
    # TE annotation can be missing for some species (e.g. a species-level
    # failure), and 1:length(species) would silently misassign the row that
    # follows a gap onto the wrong tip. Drop any species tree_tips doesn't
    # recognise (shouldn't happen, but never plot at a fabricated position).
    data_te <- data_te %>%
      mutate(node = match(species, tree_tips)) %>%
      filter(!is.na(node))
  }, error = function(e) {
    warning("Failed to load TE file: ", conditionMessage(e))
    NULL
  })
}

# --- Helper function to load FCS-GX contamination data ---
load_fcs <- function(file, tree_tips) {
  if (is.null(file)) return(NULL)
  tryCatch({
    data_fcs <- read.csv(file, sep = "\t")
    # node = this species' real index in tree_tips, not its row position -
    # FCS-GX only runs for samples with a taxid set (a documented, per-sample
    # opt-in), so this table is routinely a subset of all species. Using
    # 1:length(species) here would silently misassign a present species'
    # pie onto a DIFFERENT tip's row rather than just omitting the absent
    # ones. Drop any species tree_tips doesn't recognise (shouldn't happen,
    # but never plot at a fabricated position).
    data_fcs <- data_fcs %>%
      mutate(node = match(species, tree_tips)) %>%
      filter(!is.na(node))
  }, error = function(e) {
    warning("Failed to load FCS file: ", conditionMessage(e))
    NULL
  })
}

# --- Load optional input files ---
data_busco_geno <- load_busco(args$busco_geno, tips_order)
data_busco_prot <- load_busco(args$busco_prot, tips_order)
data_quast <- load_quast(args$quast_file, tips_order)
data_genes <- load_genes(args$genes_file, tips_order)
data_nseqs <- load_nseqs(args$nseqs_file, tips_order)
data_ortho <- load_nseqs(args$ortho_file, tips_order)
data_te <- load_te(args$te_file, tips_order)
data_fcs <- load_fcs(args$fcs_file, tips_order)

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

# Helper function to plot BUSCO pies
make_busco_scatterpie <- function(data_busco,
                                  rad_width,
                                  len_pos_x = 0,
                                  type = c("genome", "protein")) {

  if (is.null(data_busco)) {
    return(list(
      pies_plot = NULL,
      legend_busco = NULL
    ))
  }

  type <- match.arg(type)

  # Create the scatterpie plot
  pies_plot <- ggplot() +
    geom_scatterpie(
      aes(x = 0, y = node, group = species, r = rad_width),
      data = data_busco,
      cols = c("Single", "Duplicated", "Fragmented", "Missing"),
      color = NA
    ) +
    scale_fill_manual(
      values = c(
        "Single"     = "deepskyblue",
        "Duplicated" = "orange",
        "Fragmented" = "darkorchid4",
        "Missing"    = "firebrick1"
      )
    ) +
    coord_fixed() +
    theme_void() +
    ggtitle(paste0("BUSCO\n", type)) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05)
    )

  # Extract legend
  legend_busco <- cowplot::get_legend(
    pies_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8)
      )
  )

  # Remove legend from pie plot
  pies_plot <- pies_plot + guides(fill = "none")

  list(
    plot   = pies_plot,
    legend = legend_busco
  )
}

# Helper function to plot TE composition as a 100%-stacked horizontal bar
make_te_barplot <- function(data_te,
                             bar_width,
                             len_pos_x = 0) {

  if (is.null(data_te)) {
    return(list(
      plot   = NULL,
      legend = NULL
    ))
  }

  te_cols <- c("SINE", "LINE", "LTR", "Penelope", "DNA", "Rolling_Circle", "Unclassified", "Other", "Non_Repeat")

  # Reshape one column per TE category into long format for a stacked bar
  data_te_long <- data_te %>%
    select(node, all_of(te_cols)) %>%
    pivot_longer(cols = all_of(te_cols), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = te_cols))

  bars_plot <- ggplot(data_te_long, aes(x = node, y = value, fill = metric)) +
    geom_col(
      position = position_fill(reverse = TRUE), # Rescales each bar to 100% (like a pie), stacked in legend order
      width = bar_width,
      color = NA
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(
      values = c(
        "SINE"     = "deepskyblue",
        "LINE"     = "orange",
        "LTR"      = "darkorchid4",
        "Penelope" = "firebrick1",
        "DNA"      = "forestgreen",
        "Rolling_Circle" = "goldenrod",
        "Unclassified" = "purple",
        "Other" = "indianred",
        "Non_Repeat" = "darkgray"
      )
    ) +
    ggtitle("TE") +
    barplots_theme +
    # vjust=-5 (used by the two-line titles like "Genome\nsize (Mb)") pushes a
    # single-line title like this one below the panel instead of above it -
    # match N50 (Mb)'s single-line calibration instead.
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -0.4)) +
    coord_flip() +
    xlab(NULL) +
    ylab(NULL)

  # Extract legend
  legend_te <- cowplot::get_legend(
    bars_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8)
      )
  )

  # Remove legend from bar plot
  bars_plot <- bars_plot + guides(fill = "none")

  list(
    plot   = bars_plot,
    legend = legend_te
  )
}

# Helper function to plot FCS-GX contamination as a 2-category pie chart
make_fcs_piechart <- function(data_fcs,
                               rad_width,
                               len_pos_x = 0) {

  if (is.null(data_fcs)) {
    return(list(
      plot   = NULL,
      legend = NULL
    ))
  }

  pies_plot <- ggplot() +
    geom_scatterpie(
      aes(x = 0, y = node, group = species, r = rad_width),
      data = data_fcs,
      cols = c("non_contaminant_pct", "contaminant_pct"),
      color = NA
    ) +
    scale_fill_manual(
      labels = c("Non-contaminant", "Contaminant"),
      values = c(
        "non_contaminant_pct" = "deepskyblue",
        "contaminant_pct"     = "firebrick1"
      )
    ) +
    coord_fixed() +
    theme_void() +
    ggtitle("FCS") +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05)
    )

  # Extract legend
  legend_fcs <- cowplot::get_legend(
    pies_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8)
      )
  )

  # Remove legend from pie plot
  pies_plot <- pies_plot + guides(fill = "none")

  list(
    plot   = pies_plot,
    legend = legend_fcs
  )
}

# BUSCO plots
# -- if both genome and proteome busco datasets are present,
# change legend x position so that it's not skewed --
# len_pos_x <- args$len_pos_x * (!is.null(data_busco_geno) && !is.null(data_busco_prot)) # very smart chatgpt
len_pos_x <- args$len_pos_x

# Plot both genome and proteome BUSCO pies
busco_gen_plot <- make_busco_scatterpie(
  data_busco = data_busco_geno,
  rad_width  = args$rad_width,
  len_pos_x = len_pos_x,
  type       = "genome"
)

busco_prot_plot <- make_busco_scatterpie(
  data_busco = data_busco_prot,
  rad_width  = args$rad_width,
  len_pos_x = len_pos_x,
  type       = "protein"
)

# TE plots
te_plot <- make_te_barplot(
  data_te = data_te,
  bar_width = args$bar_width,
  len_pos_x = len_pos_x
)

# FCS-GX contamination plot
fcs_plot <- make_fcs_piechart(
  data_fcs = data_fcs,
  rad_width = args$rad_width,
  len_pos_x = len_pos_x
)

#if (!is.null(data_busco)) {
  # Create the scatterpie plot
#  pies_plot <- ggplot() +
#    geom_scatterpie(
#     aes(x = 0, y = node, group = species, r = args$rad_width),  # r determines the radius of the pies
#      data = data_busco,
#      cols = c("Single", "Duplicated", "Fragmented", "Missing"),
#     color = NA
#    ) +
#    scale_fill_manual(values = c("deepskyblue", "orange", "darkorchid4", "firebrick1")) +
#    coord_fixed() +
#    theme_void() +
#   ggtitle("BUSCO\ngenome") +
#   theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05))

  # Extract legend
#  legend_busco <- cowplot::get_legend(
#    pies_plot +
      #guides(fill=guide_legend(ncol=2)) +
#     theme(legend.position = "right",
#           legend.justification = c(0, 1.08),
#           legend.title = element_blank(),
#           legend.key.size = unit(0.2, "cm"),
           #legend.background = element_rect(fill = NA), # I don't know why this doesn't work, if I set this to NA an outline appears around the legend
#           legend.text = element_text(size = 8))
#  )

  # Display the legend alone
  #cowplot::ggdraw() + cowplot::draw_grob(legend_busco)

  # Remove lenged for pieplot
#  pies_plot <- pies_plot + guides(fill="none")
#} else {
#  pies_plot <- NULL
#  legend_busco <- NULL
#}

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
  get_plot_range(busco_gen_plot$plot, "y"),
  get_plot_range(busco_prot_plot$plot, "y"),
  get_plot_range(fcs_plot$plot, "y"),
  get_plot_range(len_plot, "x"),
  get_plot_range(n50_plot, "x"),
  get_plot_range(gene_plot, "x"),
  get_plot_range(te_plot$plot, "x") # te_plot is now a coord_flip()'d bar (like len_plot/n50_plot), not a pie - node range lives on x, not y
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

# Set new ylim for BUSCO pies
if (!is.null(busco_gen_plot$plot)) busco_gen_plot$plot <- busco_gen_plot$plot + new_ylim
if (!is.null(busco_prot_plot$plot)) busco_prot_plot$plot <- busco_prot_plot$plot + new_ylim

# Set new ylim for FCS pie
if (!is.null(fcs_plot$plot)) fcs_plot$plot <- fcs_plot$plot + new_ylim

# Set new xlim for TE bars (equivalent to ylim, since te_plot is coord_flip()'d)
if (!is.null(te_plot$plot)) te_plot$plot <- te_plot$plot + new_xlim


# Set new xlim for gene stats (equivalent to ylim)
if (!is.null(gene_plot)) gene_plot <- gene_plot + new_xlim

# --- Quality scoring ----------------------------------------------------------
# Quality metrics are drawn as a discrete "traffic light" instead of a sequential
# ramp, because a light->dark ramp implies "dark = good", which is wrong for
# lower-is-better stats (e.g. sequence count) and meaningless for descriptive
# ones (genome size, gene number). Colours are the colour-vision-safe Okabe-Ito
# triple, so they stay distinguishable under red-green colour blindness.
QUALITY_COLOURS <- c(Good = "#009E73", Warn = "#E69F00", Poor = "#D55E00")

# Threshold presets by phylogenetic group.
# NOTE: the BUSCO cut-offs follow community practice; the N50 and sequence-count
# cut-offs are clade-dependent STARTING POINTS and should be tuned per project.
# 'generic' is deliberately lenient - eukaryote/vertebrate standards are far too
# strict for, say, a bacterial assembly.
# fcs_noncontam (FCS-GX non-contaminant %) is the same across every preset: unlike
# BUSCO/N50/sequence-count, "how much of the genome is foreign contamination" has
# no taxon-dependent expectation, so there's no reason to vary it by clade.
QUALITY_PRESETS <- list(
  generic = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  vertebrate = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e7,    warn = 1e6),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  insect = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  plant = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 10,     warn = 20),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 5000,   warn = 50000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  fungi = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 100,    warn = 1000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  bacteria = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 5e5,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 10,     warn = 100),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  )
)

# Bin values into Good / Warn / Poor given a threshold spec
classify_quality <- function(values, thr) {
  if (is.null(thr) || is.null(values)) return(NULL)
  v <- as.numeric(values)
  out <- rep(NA_character_, length(v))
  if (identical(thr$direction, "higher")) {
    out[v >= thr$good]                  <- "Good"
    out[v <  thr$good & v >= thr$warn]  <- "Warn"
    out[v <  thr$warn]                  <- "Poor"
  } else {
    out[v <= thr$good]                  <- "Good"
    out[v >  thr$good & v <= thr$warn]  <- "Warn"
    out[v >  thr$warn]                  <- "Poor"
  }
  factor(out, levels = names(QUALITY_COLOURS))
}

# Parse --quality_thresholds ('metric=good:warn' pairs) into the same shape as
# a QUALITY_PRESETS entry. Direction is looked up per metric (never taken from
# the CLI) so a custom threshold cannot silently invert a metric's meaning.
QUALITY_METRIC_DIRECTIONS <- vapply(QUALITY_PRESETS$generic, `[[`, character(1), "direction")

parse_quality_thresholds <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(NULL)
  out <- list()
  for (part in strsplit(spec, ",")[[1]]) {
    part <- trimws(part)
    if (!nzchar(part)) next
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      stop("Malformed --quality_thresholds entry '", part, "': expected 'metric=good:warn'")
    }
    metric <- trimws(kv[1])
    if (!metric %in% names(QUALITY_METRIC_DIRECTIONS)) {
      stop("Unknown metric '", metric, "' in --quality_thresholds. Expected one of: ",
           paste(names(QUALITY_METRIC_DIRECTIONS), collapse = ", "))
    }
    gw <- strsplit(trimws(kv[2]), ":", fixed = TRUE)[[1]]
    if (length(gw) != 2 || anyNA(suppressWarnings(as.numeric(gw)))) {
      stop("Malformed cut-offs for '", metric, "' in --quality_thresholds: expected 'good:warn' numbers")
    }
    out[[metric]] <- list(direction = QUALITY_METRIC_DIRECTIONS[[metric]],
                           good = as.numeric(gw[1]), warn = as.numeric(gw[2]))
  }
  out
}

# --- Circular layout helper ---------------------------------------------------
# The circular ("fan") layout is a separate plotting path: instead of the
# concatenated side panels, each stat is drawn as a concentric coloured ring
# around a fan tree, tips are numbered, and a species key is shown alongside.
build_circular_plot <- function(tree, tips_order, data_quast = NULL, data_genes = NULL,
                                data_busco_geno = NULL, data_busco_prot = NULL,
                                data_nseqs = NULL, data_ortho = NULL, data_fcs = NULL,
                                text_size = 3, skip = NULL, rings = NULL,
                                quality_preset = "generic", thresholds = NULL, show_values = FALSE,
                                open_angle = 14, ring_width = 0.13) {

  if (is.null(skip)) skip <- character(0)

  thr_set <- QUALITY_PRESETS[[quality_preset]]
  if (is.null(thr_set)) {
    warning("Unknown quality_preset '", quality_preset, "', falling back to 'generic'")
    thr_set <- QUALITY_PRESETS[["generic"]]
  }
  # Explicit per-metric thresholds (e.g. from the Shiny app) override the preset
  if (!is.null(thresholds)) {
    for (m in names(thresholds)) {
      if (!is.null(thresholds[[m]])) thr_set[[m]] <- thresholds[[m]]
    }
  }

  # Single source of truth for key typography: headers one step larger than body
  # text. Used by the hand-drawn keys (geom text, mm) and - via .pt - the ggplot
  # legends (theme text, points), so the two cannot drift apart.
  key_title_size <- text_size * 1.15
  key_text_size  <- text_size * 0.95

  n_tips <- length(tips_order)
  labs <- gsub("_", " ", tips_order)              # tree tip labels in node order
  num_df  <- data.frame(label = labs, num = seq_len(n_tips), stringsAsFactors = FALSE)
  ring_df <- data.frame(label = labs, stringsAsFactors = FALSE)

  ring_specs <- list()
  add_ring <- function(spec, values) {
    if (is.null(values) || all(is.na(values))) return(invisible())
    name <- spec$name
    if (identical(spec$scale, "quality")) {
      thr <- thr_set[[spec$metric]]
      grade <- classify_quality(values, thr)
      if (is.null(grade)) return(invisible())
      ring_df[[name]] <<- grade
      spec$values <- values          # keep raw values for the printed labels
    } else {
      ring_df[[name]] <<- values
      spec$values <- values
    }
    ring_specs[[length(ring_specs) + 1]] <<- spec
  }
  # pull a stat's per-tip values in node (i.e. labs) order
  # Aligns df[[col]] to the ring's full n_tips positions by node index, NA
  # everywhere df has no row. df is routinely a subset (e.g. FCS-GX only runs
  # for taxid-bearing samples) - returning just the present rows in df's own
  # order would give a short vector that R silently RECYCLES across every
  # tip when it's assigned into ring_df, fabricating data for species that
  # have none at all instead of leaving them blank.
  col_by_node <- function(df, col) {
    out <- rep(NA_real_, n_tips)
    out[df$node] <- as.numeric(df[[col]])
    out
  }

  # Registry of every stat that can be a ring. Each key matches its
  # rectangular-layout panel (and the Shiny "Skip Statistics" checkboxes).
  # `scale` is "quality" (discrete traffic light, needs `metric` for thresholds)
  # or "descriptive" (sequential ramp, no good/bad implied). `get` returns the
  # per-tip values or NULL when that data was not supplied.
  ring_registry <- list(
    ch_plot         = list(name = "Seq number",     scale = "quality", metric = "seq_number",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "Sequences")),
    len_plot        = list(name = "Genome size",    scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "Length")),
    n50_plot        = list(name = "N50",            scale = "quality", metric = "n50",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "N50")),
    gene_plot       = list(name = "Gene number",    scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() {
                             if (is.null(data_genes)) return(NULL)
                             tot <- data_genes[data_genes$stat == "Total", ]
                             as.numeric(tot[order(tot$node), ]$value)
                           }),
    busco_gen_plot  = list(name = "BUSCO genome",   scale = "quality", metric = "busco_complete",
                           get = function() if (!is.null(data_busco_geno)) data_busco_geno$Single + data_busco_geno$Duplicated),
    busco_prot_plot = list(name = "BUSCO protein",  scale = "quality", metric = "busco_complete",
                           get = function() if (!is.null(data_busco_prot)) data_busco_prot$Single + data_busco_prot$Duplicated),
    busco_dup_plot  = list(name = "BUSCO duplicated", scale = "quality", metric = "busco_duplicated",
                           get = function() if (!is.null(data_busco_geno)) data_busco_geno$Duplicated),
    fcs_plot        = list(name = "FCS non-contaminant %", scale = "quality", metric = "fcs_noncontam",
                           get = function() if (!is.null(data_fcs)) col_by_node(data_fcs, "non_contaminant_pct")),
    nseqs_plot      = list(name = "Seqs ≥5 BUSCOs", scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_nseqs)) col_by_node(data_nseqs, names(data_nseqs)[2])),
    ortho_plot      = list(name = "Ortho seqs",     scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_ortho)) col_by_node(data_ortho, names(data_ortho)[2]))
  )

  # `rings` fixes which stats to show and their inner->outer order (the curated
  # static default). When NULL (e.g. the Shiny app) every available stat is
  # shown. Either way `skip` is still honoured, so a ring is drawn only if it is
  # selected, not skipped, and its data is present.
  order_keys <- if (!is.null(rings)) rings else names(ring_registry)
  for (key in order_keys) {
    entry <- ring_registry[[key]]
    if (is.null(entry) || key %in% skip) next
    add_ring(entry, entry$get())
  }

  # Keep the traffic-light (quality) rings visually separate from the neutral
  # descriptive rings: descriptive inner (nearest the tree), quality outer (on the
  # rim, where they read most clearly), preserving the requested order within each
  # group.
  is_quality <- vapply(ring_specs, function(s) identical(s$scale, "quality"), logical(1))
  ring_specs <- c(ring_specs[!is_quality], ring_specs[is_quality])

  # Fan tree
  p <- ggtree(tree, layout = "fan", open.angle = open_angle, size = 0.5, colour = "grey30")

  n_rings <- length(ring_specs)

  ring_offset <- 0.055

  if (n_rings > 0) {
    # PRIMER RING: geom_fruit()'s very first call on a fresh ggtree allocates a
    # one-off, disproportionately large slice of radius regardless of width/
    # offset (a ggtreeExtra quirk, not something width/offset can compensate
    # for) - every ring after the first is sized consistently. A fully
    # transparent zero-content ring absorbs that anomaly so every REAL ring
    # below gets uniform, comparable thickness. Deliberately uses only width/
    # offset (never pwidth: explicitly passing it - at any value - corrupts
    # the whole plot's scale in the ggtreeExtra version this pipeline pins).
    p <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1), fill = NA, colour = NA,
        width = ring_width, offset = 0.10
      ) +
      ggnewscale::new_scale_fill()

    # CALIBRATION: geom_fruit's `width` is not the ring's actual rendered
    # thickness - rings are drawn as overlapping tiles staggered by `offset`,
    # each one visually clipped by the next ring drawn on top of it down to
    # the gap between their start positions. The last ring has nothing after
    # it to clip its trailing edge, so it alone renders at its full declared
    # `width` - several times wider than every other (clipped) ring. Measure
    # the true stagger `ring_offset` produces on this tree with a throwaway
    # invisible ring, then use that measured value (not ring_width) as every
    # real ring's width below, so every ring - including the last - renders
    # the same actual thickness.
    range_before_calib <- suppressWarnings(max(vapply(
      ggplot_build(p)$data,
      function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
      numeric(1)), na.rm = TRUE))
    p_calib <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1), fill = NA, colour = NA,
        width = ring_width, offset = ring_offset
      )
    range_after_calib <- suppressWarnings(max(vapply(
      ggplot_build(p_calib)$data,
      function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
      numeric(1)), na.rm = TRUE))
    measured_width <- range_after_calib - range_before_calib
  }

  # Add each stat as a concentric ring. Quality rings share one discrete
  # Good/Warn/Poor scale (so the legend is shown only once - the ring key on the
  # left says which ring is which); descriptive rings keep a sequential ramp,
  # each with its own legend ordered outer-ring-first.
  shown_quality_legend <- FALSE
  for (i in seq_along(ring_specs)) {
    spec <- ring_specs[[i]]
    p <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1, fill = .data[[spec$name]]),
        width = measured_width, offset = ring_offset,
        color = "white", linewidth = 0.2
      )
    if (identical(spec$scale, "quality")) {
      # The Good/Warn/Poor key is drawn manually in the side panel: ggplot only
      # renders legend keys for levels present in the layer that owns the legend,
      # so Warn/Poor vanish whenever the first quality ring happens to be all-Good.
      p <- p + scale_fill_manual(
        values = QUALITY_COLOURS, limits = names(QUALITY_COLOURS),
        drop = FALSE, na.value = "grey90", name = "Quality", guide = "none"
      )
      shown_quality_legend <- TRUE
    } else {
      p <- p + scale_fill_gradient(
        low = spec$low, high = spec$high, name = spec$name,
        guide = guide_colourbar(order = n_rings - i + 1)
      )
    }
    p <- p + ggnewscale::new_scale_fill()
  }

  # Outermost labels: tip numbers, kept upright. geom_fruit/geom_tiplab rotate
  # text tangentially in a fan layout, so instead we place a plain geom_text
  # (angle = 0 keeps the glyphs horizontal) just beyond the outer ring radius.
  ring_max_x <- suppressWarnings(max(vapply(
    ggplot_build(p)$data,
    function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
    numeric(1)), na.rm = TRUE))
  tip_pos <- p$data[p$data$isTip, ]
  num_pos <- data.frame(y = tip_pos$y, num = match(tip_pos$label, labs))

  # All circular text (tip numbers, legends, species key) scales with text_size.
  p <- p +
    geom_text(data = num_pos, aes(x = ring_max_x * 1.06, y = y, label = num),
              angle = 0, size = text_size * 0.9, inherit.aes = FALSE) +
    theme(legend.position = "right",
          legend.title = element_text(size = key_title_size * .pt, face = "bold"),
          legend.text = element_text(size = key_text_size * .pt),
          legend.key.width = unit(0.3, "cm"),
          legend.key.height = unit(0.35, "cm"))

  # Optional redundant encoding: print each value on its ring, so the figure is
  # not colour-only (important for the traffic-light rings). Kept upright.
  if (isTRUE(show_values) && n_rings > 0) {
    fmt_val <- function(v) {
      v <- as.numeric(v)
      m <- suppressWarnings(max(abs(v), na.rm = TRUE))
      if (!is.finite(m))                                    rep("", length(v))
      else if (m >= 1e6)                                    sprintf("%.1fM", v / 1e6)
      else if (m >= 1000)                                   format(round(v), big.mark = ",", trim = TRUE)
      else if (all(abs(v - round(v)) < 1e-8, na.rm = TRUE)) as.character(round(v))
      else                                                  sprintf("%.1f", v)
    }
    # Radial centre of each ring, read back from the rendered tile layers.
    # tile_layers[[1]] is the invisible primer ring (always added above when
    # n_rings > 0), not a real ring - skip it so values land on the ring they
    # actually describe instead of the one drawn just inside it.
    tile_layers <- Filter(function(dd) all(c("xmin", "xmax") %in% names(dd)),
                          ggplot_build(p)$data)
    if (length(tile_layers) >= n_rings + 1) {
      for (i in seq_len(n_rings)) {
        dd <- tile_layers[[i + 1]]
        # NB: the radius must live in the data, not the aes expression - aes() is
        # evaluated lazily, so aes(x = rad) would resolve every layer to the last
        # value of the loop variable.
        val_pos <- data.frame(
          x   = mean(c(dd$xmin, dd$xmax), na.rm = TRUE),
          y   = tip_pos$y,
          lab = fmt_val(ring_specs[[i]]$values)[match(tip_pos$label, labs)]
        )
        p <- p + geom_text(data = val_pos, aes(x = x, y = y, label = lab),
                           angle = 0, size = text_size * 0.6, colour = "grey15",
                           inherit.aes = FALSE)
      }
    }
  }

  # Species key (number -> italic name): a compact, top-aligned list on the left
  sp_txt <- paste0(num_df$num, "  ", num_df$label)
  sp_block <- paste(sp_txt, collapse = "\n")

  # Ring key (inner -> outer). Needed because all quality rings share the same
  # Good/Warn/Poor colours, so colour alone cannot identify a ring.
  ring_names <- vapply(ring_specs, function(s) s$name, character(1))
  ring_block <- paste0(seq_along(ring_names), ". ", ring_names, collapse = "\n")
  n_ring_lines <- length(ring_names)

  has_quality <- any(vapply(ring_specs, function(s) identical(s$scale, "quality"), logical(1)))
  lh <- 0.035                                   # one text line, in panel units
  y  <- 1.00
  sp_leg <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
    annotate("text", x = 0, y = y, label = "Species",
             hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
  y <- y - lh * 1.6
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = sp_block,
             hjust = 0, vjust = 1, size = key_text_size, fontface = "italic", lineheight = 1.2)
  y <- y - lh * n_tips - lh * 0.8
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = "Rings (inner -> outer)",
             hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
  y <- y - lh * 1.6
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = ring_block,
             hjust = 0, vjust = 1, size = key_text_size, lineheight = 1.2)
  y <- y - lh * length(ring_names) - lh * 0.8

  # Manual Good / Warn / Poor key - always shows all three swatches
  if (has_quality) {
    sp_leg <- sp_leg +
      annotate("text", x = 0, y = y, label = "Quality",
               hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
    y <- y - lh * 1.4
    for (k in seq_along(QUALITY_COLOURS)) {
      yy <- y - (k - 1) * lh
      sp_leg <- sp_leg +
        annotate("rect", xmin = 0, xmax = 0.05,
                 ymin = yy - lh * 0.62, ymax = yy - lh * 0.08,
                 fill = QUALITY_COLOURS[[k]], colour = NA) +
        annotate("text", x = 0.075, y = yy - lh * 0.35,
                 label = names(QUALITY_COLOURS)[k],
                 hjust = 0, vjust = 0.5, size = key_text_size)
    }
    y <- y - lh * length(QUALITY_COLOURS) - lh * 0.8

    # Spell out what earns each grade, per quality ring, so the figure is
    # self-documenting rather than relying on the reader knowing the preset.
    fmt_thr <- function(v) {
      if (!is.finite(v)) return("?")
      if (abs(v) >= 1e6)       sprintf("%.1fM", v / 1e6)
      else if (abs(v) >= 1000) format(round(v), big.mark = ",", trim = TRUE)
      else                     as.character(round(v, 1))
    }
    rules <- character(0)
    for (i in seq_along(ring_specs)) {
      sp  <- ring_specs[[i]]
      if (!identical(sp$scale, "quality")) next
      thr <- thr_set[[sp$metric]]
      if (is.null(thr)) next
      rule <- if (identical(thr$direction, "higher")) {
        sprintf("%s+ / %s+ / <%s", fmt_thr(thr$good), fmt_thr(thr$warn), fmt_thr(thr$warn))
      } else {
        sprintf("<=%s / <=%s / >%s", fmt_thr(thr$good), fmt_thr(thr$warn), fmt_thr(thr$warn))
      }
      rules <- c(rules, paste0(i, ". ", sp$name, ": ", rule))
    }
    if (length(rules) > 0) {
      sp_leg <- sp_leg +
        annotate("text", x = 0, y = y, label = "Thresholds (Good / Warn / Poor)",
                 hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
      y <- y - lh * 1.4
      sp_leg <- sp_leg +
        annotate("text", x = 0, y = y, label = paste(rules, collapse = "\n"),
                 hjust = 0, vjust = 1, size = key_text_size, lineheight = 1.25)
    }
  }

  sp_leg + p + plot_layout(widths = c(0.35, 1))
}

if (args$tree_style == "circular") {
  # 'all' shows every available stat (respecting skip); otherwise use the
  # curated, ordered ring list from --circular_rings.
  circular_rings <- if (identical(tolower(args$circular_rings), "all")) {
    NULL
  } else {
    trimws(strsplit(args$circular_rings, ",")[[1]])
  }
  final_plot <- build_circular_plot(
    tree            = tree,
    tips_order      = tips_order,
    data_quast      = data_quast,
    data_genes      = data_genes,
    data_busco_geno = data_busco_geno,
    data_busco_prot = data_busco_prot,
    data_nseqs      = data_nseqs,
    data_ortho      = data_ortho,
    data_fcs        = data_fcs,
    text_size       = args$text_size,
    skip            = skip,
    rings           = circular_rings,
    quality_preset  = args$quality_preset,
    thresholds      = parse_quality_thresholds(args$quality_thresholds),
    show_values     = isTRUE(args$show_ring_values)
  )
} else {

# Build tree according to the selected style
if (args$tree_style == "rectangular") {
  # Legacy look: thin black branches with dotted alignment leader lines
  tree_plot <- ggtree(tree) +
    geom_tiplab(size = args$text_size, fontface = "italic", align = TRUE, hjust = -0.05)
} else {
  # Modern look: thicker grey branches, aligned labels without dotted leaders
  tree_plot <- ggtree(tree, layout = args$tree_style, size = 0.7, colour = "grey30") +
    geom_tiplab(size = args$text_size, fontface = "italic", align = TRUE,
                linetype = NA, hjust = -0.05)
  if (args$tree_style == "ellipse") {
    # Subtle node markers to accentuate the curved layout
    tree_plot <- tree_plot + geom_nodepoint(colour = "steelblue", size = 1.2, alpha = 0.75)
  }
}
tree_plot <- tree_plot +
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
  busco_gen_plot  = busco_gen_plot$plot,
  busco_prot_plot = busco_prot_plot$plot,
  te_plot = te_plot$plot,
  fcs_plot = fcs_plot$plot
)

all_legends <- list(
  ch_plot    = NULL,
  nseqs_plot = NULL,
  ortho_plot = NULL,
  len_plot   = legend_len,
  gene_plot  = legend_gene,
  n50_plot   = NULL,
  busco_gen_plot  = busco_gen_plot$legend,
  busco_prot_plot = if (!is.null(busco_gen_plot$legend)) NULL else busco_prot_plot$legend, # Only plot legend once
  te_plot = te_plot$legend,
  fcs_plot = fcs_plot$legend
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

} # end of non-circular (side-panel) layout branch

# Circular plots read better on a larger, squarer canvas
plot_w <- if (args$tree_style == "circular") 12 else 10
plot_h <- if (args$tree_style == "circular") 8.5 else 7

pdf("tree_plot.pdf", width = plot_w, height = plot_h)
final_plot
dev.off()

svg("tree_plot.svg", width = plot_w, height = plot_h)
final_plot
dev.off()
