#!/usr/bin/Rscript

# Written by Chris Wyatt and Fernando Duarte and released under the MIT license.
# Plots the phylogenetic tree with BUSCO, Quast and gene stats results

# Function to plot tree and plots
build_tree_plot <- function(tree, n, plots, legends, xlimit, rigth_margin, bottom_margin) { #xlim for legends, use same xlim as barplots (new_xlim)
  # Update xlim for the tree plot
  tree <- tree + ggplot2::xlim(0, n)

  # Initialize combined plot with the tree
  combined_plots <- tree
  widths <- c(10 - length(plots)) # Dynamic width for the tree

  # Add each additional plot
  for (plot in plots) {
    combined_plots <- combined_plots | plot
    widths <- c(widths, 1) # Same widths for plots and legends
  }

  # Initialize combined legends with empty plot aligned w/ tree
  combined_legends <- plot_spacer() + xlimit

  # Add each additional legend
  for (legend in legends) {
    if (length(legend) != 0) {
      legend <- legend#wrap_elements(rotate_grob(legend, -45)) + xlimit
    } else {
      legend <- plot_spacer() + xlimit
    }
    combined_legends <- combined_legends | legend
  }

  # Apply the layout widths
  combined_plots <- combined_plots + plot_layout(widths = widths)

  combined_legends <- combined_legends +
    plot_layout(widths = widths) +
    theme(plot.margin = margin(0, rigth_margin, bottom_margin, 0))

  combined_plots <- combined_plots / combined_legends +
    plot_layout(heights = c(0.99, 0.01))

  return(combined_plots)
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
parser$add_argument('--text_size', type = 'double', default = 3, help = 'Text size for the tree plot')
parser$add_argument('--tree_scale', type = 'double', default = 0.0005, help = 'x axis limits scaling for tree plot (useful when tree labels appear truncated)')
parser$add_argument('--bar_width', type = 'double', default = 0.7, help = 'Width of bar plots')
parser$add_argument('--rad_width', type = 'double', default = 0.4, help = 'Radius of pie charts')
args <- parser$parse_args()

# Avoid scientific notation in all plots
options(scipen = 999)

# Read the Newick tree from the file
tree <- read.tree(args$tree_file)

# Clean tree tip labels
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- tolower(tree$tip.label)

# If radious of pie charts is to big, it can
# mess the position of the pies, make them
# smaller
if (length(tree$tip.label) < 7) {
  args$bar_width <- args$bar_width/1.5
  args$rad_width <- args$rad_width/2
}

# Capitalize first letter of tip labels
tree$tip.label <- sub("^(\\w)(.*)", "\\U\\1\\L\\2", tree$tip.label, perl = TRUE)

# Read the data table from the file, ensuring species column is read as character
# Load BUSCO
data_busco <- read.csv(args$busco_file, sep = "\t", colClasses = c("Input_file" = "character"))

# Prepare BUSCO data
data_busco <- data_busco %>%
  # Remove extension from Input_file
  mutate(Input_file = tools::file_path_sans_ext(Input_file)) %>%
  # Rename 'Input_file' to 'species'
  rename(species = Input_file)

# Load Quast
data_quast <- read.csv(args$quast_file, sep = "\t")

# Change header of GC% and contigs column
colnames(data_quast)[5] <- "GC"
colnames(data_quast)[6] <- "Sequences"

#Prepare Quast data
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

# Load gene stats
data_genes <- read.csv(args$genes_file, sep = "\t")

# Prepare gene stats
data_genes <- data_genes %>%
  # Rename columns
  rename(species = File) %>%
  rename(Total = Total_genes) %>%
  rename(Overlapping = Overlapping_genes) %>%
  # Remove prefix and suffix
  mutate(species = gsub("Count\\.|\\.tsv", "", species))

# Extract names for debugging
tree_sp <- sort(tree$tip.label)
quast_sp <- sort(unique(data_quast$species))
busco_sp <- sort(data_busco$species)
gene_sp <- sort(unique(data_genes$species))

# Debugging: Print species names from the tree and the data
cat("Species names in the tree based on nw:\n")
print(tree_sp )
cat("\nSpecies names in data tables:\n")
cat("BUSCO:", busco_sp, "\nQUAST:", quast_sp, "\ngene_stats:", gene_sp)

# Debugging: Check if there are any mismatches in species names
datasets <- list(BUSCO = busco_sp, Quast = quast_sp, GeneStats = gene_sp)
for (name in names(datasets)) {
  if (any(tree_sp != datasets[[name]])) {
    stop(paste("Species names in", name, "and tree labels do not match"))
  }
}

# Get order of tips
tree_plot <- ggtree(tree) # Temporary plot for get_taxa_name()
tips_order <- rev(get_taxa_name(tree_plot))

# Arrange data according to tree labels
data_busco <- data_busco %>%
  arrange(match(species, tips_order))
data_quast <- data_quast %>%
  arrange(match(species, tips_order)) %>%
  mutate(node = 1:length(species))
data_genes <- data_genes %>%
  arrange(match(species, tips_order)) %>%
  mutate(node = 1:length(species))

# Match names with new tree tips (only necessary for Quast)
# data_quast$species <- gsub("_", " ", data_quast$species)

# Tidy Quast data
# For N50/N90
data_quast_n5090 <- data_quast %>%
  # Convert wide to long format
  pivot_longer(cols = c(N50, N90),
               names_to = "metric",
               values_to = "value") %>%
  # Convert value column to numeric if needed
  mutate(value = as.numeric(value)) %>%
  mutate(value = (as.numeric(value)/1000000)) # Values in Mb

# For GC content and length
data_quast_len <- data_quast %>%
  pivot_longer(cols = c(GC, Length),
               names_to = "metric",
               values_to = "value")

# Tidy gene stats data
data_genes <- data_genes %>%
  pivot_longer(cols = c(Total, Overlapping),
               names_to = "stat",
               values_to = "value")

# Add node column
data_busco <- data_busco %>%
  mutate(node = 1:length(species)) # Node number needed for nodpie

# This is for the synteny paper, remove "_" and changes the first letter to uppercase
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Plot number of chromosomes/sequences
ch_plot <- ggplot(data_quast, aes(x=1, y=node)) +
  geom_text(aes(label = Sequences)) +
  theme_void() +
  ggtitle("Sequence\nnumber") +
  theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))

# Now with the bar plots. Set standard theme for all barplots
barplots_theme <- theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 6),
    axis.ticks.y=element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_blank()
  )

# Prepare Quast data for plotting
data_quast_n50 <- data_quast_n5090[data_quast_n5090$metric %in% "N50",]
#data_quast_n90 <- data_quast_n5090[data_quast_n5090$metric %in% "N90",]

# Plot Quast data genome size
len_plot <- ggplot(
  data_quast_len,
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

# Display the legend alone
#cowplot::ggdraw() + cowplot::draw_grob(legend_len)

# Remove legend
len_plot <- len_plot + guides(fill="none")

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

# Display the legend alone
#cowplot::ggdraw() + cowplot::draw_grob(legend_gene)

# Remove legend from pie plot
pies_plot <- pies_plot + guides(fill="none")

# Remove legend from gene plot
gene_plot <- gene_plot + guides(fill="none")

# To match tips with plots it's necessary to set the same ylim for all plots
# Select the biggest range to avoid cropping
p_ranges_y <- c(
  ggplot_build(ch_plot)$layout$panel_scales_y[[1]]$range$range,
  ggplot_build(len_plot)$layout$panel_scales_x[[1]]$range$range,
  ggplot_build(n50_plot)$layout$panel_scales_x[[1]]$range$range,
  ggplot_build(pies_plot)$layout$panel_scales_y[[1]]$range$range
)

# Set new ylim based on the highest value taking into account both plots
# Add -0,1 and 0.5 to avoid the cropping of first and last pies
new_ylim <- ylim(c(min(p_ranges_y), max(p_ranges_y)))
# A new xlim is needed for barplots (equivalent to ylim), as these are flipped
# using coord_flip()
new_xlim <- xlim(c(min(p_ranges_y), max(p_ranges_y)))

# Set new ylim for sequnces
ch_plot <- ch_plot + new_ylim

# Set new xlim for Quast genome size (equivalent to ylim)
len_plot <- len_plot + new_xlim

# Set new xlim for Quast N50 (equivalent to ylim)
n50_plot <- n50_plot + new_xlim

# Set new ylim for Quast pies
pies_plot <- pies_plot + new_ylim

# Set new xlim for gene stats (equivalent to ylim)
gene_plot <- gene_plot + new_xlim

# Build tree
tree_plot <- ggtree(tree) +
  # Tip font size, should be an arg
  geom_tiplab(size=3, fontface = "italic", align = TRUE) +
  theme(plot.margin = margin(10, 10, 10, 10)) + # Increase margins
  coord_cartesian(clip="off")

# Set new ylim and xlim for tree
tree_plot <- tree_plot + new_ylim

# Set value for tree xlim to avoid the truncation of labels:
# Why "^2*0.001"? ^2 is because the relatin between number of characters and the number
# of pixels is close to beexponential, not proportional. 0.001 would be the length
# per character in the x axis scale. Script should allow to change this value
m = max(tree_plot$data$x) + max(nchar(tree_plot$data$label))^2*args$tree_scale

# Call the function
final_plot <- build_tree_plot(
  tree = tree_plot,
  n = m, # Only affects tree_plot
  plots = list(ch_plot, len_plot, gene_plot, n50_plot, pies_plot),
  legends = list(NULL, legend_len, legend_gene, NULL, legend_busco),
  new_xlim,
  15,
  60
)

pdf("tree_plot.pdf")
final_plot
dev.off()
