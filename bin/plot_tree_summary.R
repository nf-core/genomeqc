#!/usr/bin/Rscript

# Written by Chris Wyatt and Fernando Duarte and released under the MIT license.
# Plots a phylogenetic tree alongside BUSCO and Quast stats

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
if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
  BiocManager::install("ggtreeExtra")
}
if (!requireNamespace("ggimage", quietly = TRUE)) {
  install.packages("ggimage")
}

library(ggtree)
library(ggplot2)
library(cowplot)
library(argparse)
library(dplyr)
library(tidyr)
library(ggimage)
library(ggtreeExtra)

# Parse command-line arguments
parser <- ArgumentParser(description = 'Plot phylogenetic tree with statistics and true/false data')
parser$add_argument('tree_file', type = 'character', help = 'Path to the Newick formatted tree file')
parser$add_argument('busco_file', type = 'character', help = 'Path to processed BUSCO output file')
parser$add_argument('quast_file', type = 'character', help = 'Path to processed Quast output file')
parser$add_argument('--text_size', type = 'double', default = 2.5, help = 'Tree tip text size for the plots')
parser$add_argument('--pie_size', type = 'double', default = 0, help = 'Increase pie size of each tip by this fold')
parser$add_argument('--pie_hjust', type = 'double', default = 0.4, help = 'Horizontally adjust pie positions')
#parser$add_argument('--tree_size', type = 'double', default = 0.3, help = 'Proportion of the plot width for the tree')
args <- parser$parse_args()

# Create a function so that, if the Quast data is off limits and warnings related
# to this pop up , margins will be increased by the value "step" on each iteration
increase_scale_until_no_warnings <- function(plot, initial_limit = 5, step = 0.1, max_iterations = 50) {
  # Initialize limit and warning flag
  current_limit <- initial_limit
  warning_flag <- TRUE
  iteration <- 0

  while (warning_flag && iteration < max_iterations) {
    iteration <- iteration + 1
    # Update the plot with the current limit
    updated_plot <- plot + scale_x_continuous(limits = c(0, current_limit))

    # Capture warnings
    warnings <- tryCatch({
      print(updated_plot)  # Render the plot
      NULL  # No warnings
    }, warning = function(w) {
      w$message  # Capture warning message
    })

    # Check if there are warnings
    if (is.null(warnings)) {
      warning_flag <- FALSE  # No warnings, exit loop
    } else {
      # Increase the limit for the next iteration
      current_limit <- current_limit + step
    }
  }

  # Return the final plot
  if (iteration >= max_iterations) {
    warning("Maximum iterations reached. Plot may still have issues.")
  }
  return(updated_plot)
}

# Read the Newick tree from the file
tree <- read.tree(args$tree_file)

# Clean tree tip labels
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- tolower(tree$tip.label)

# Read the data table from the file, ensuring species column is read as character
# Load BUSCO
data_busco <- read.csv(args$busco_file, sep = "\t", colClasses = c("Input_file" = "character"))

# Prepare BUSCO data
data_busco <- data_busco %>%
  # Remove extension from Input_file
  mutate(Input_file = tools::file_path_sans_ext(Input_file)) %>%
  # Rename 'Input_file' to 'species'
  rename(species = Input_file) %>%
  # Make everything lowercase
  mutate(species = tolower(species))

# Load Quast data
data_quast <- read.csv(args$quast_file, sep = "\t")

#Tidy Quast data
data_quast <- data_quast %>%
  # Remove the row where species is NA
  filter(!is.na(species)) %>%
  # Make species lowercase
  mutate(species = tolower(species)) %>%
  # Convert wide to long format
  pivot_longer(cols = c(N50, N90),
               names_to = "metric",
               values_to = "value") %>%
  # Remove any remaining "bar" rows if necessary (check Chris script)
  filter(value != "bar") %>%
  # Convert value column to numeric if needed
  mutate(value = as.numeric(value))

# Debugging: Print species names from the tree and the data
cat("Species names in the tree based on nw:\n")
cat(tree$tip.label)
cat("\nSpecies names in data tables:\n")
quast_sp <- unique(data_quast$species)
cat("BUSCO:", data_busco$species, "\nQUAST:", quast_sp, "\n")

# Debugging: Check if there are any mismatches in species names (might need to be
# changed in the future to support optional stats)
if (all(sort(tree$tip.label) != sort(data_busco$species))) {
  stop("Species names in BUSCO and tree labels do not match")
} else if (all(sort(quast_sp) != sort(data_busco$species))) {
  stop("Species names in Quast and tree labels do not match")
}

# Arrange data according to tree labels, so that they are plotted correctly
# For BUSCO
data_busco <- data_busco %>%
  arrange(match(species, tree$tip.label))

# For Quast
data_quast <- data_quast %>%
  arrange(match(species, tree$tip.label))

# Node number needed in BUSCO dataframe for nodpie. This is a generic 1-n
# numbering (n=number of species)
data_busco <- data_busco %>%
  mutate(node = seq_len(n()))

# Plot the phylogenetic tree
tree_plot <- ggtree(tree, branch.length = "none") + # No branch length, only topology
  geom_tiplab(size=args$text_size) + # Tip font size
  ggplot2::xlim(0, 5) +
  ggtitle("Phylogenetic Tree") +
  theme(plot.margin = margin(10, 10, 10, 10))  # Increase margins

# Save tree only
ggsave("tree_only.pdf", tree_plot)

# Create list of pies from BUSCO data
pies <- nodepie(data_busco,
                cols = 4:7,
                color = c(
                  "Single" = "steelblue",
                  "Duplicated" = "orange",
                  "Fragmented" = "purple",
                  "Missing" = "red"
                )
)

print("Check 1")

# Create a dummy pie as for the legend it seems impossible to extract it
# from nodpie
dummy_pie_data <- data.frame(
  Type = factor(c("Single", "Duplicated", "Fragmented", "Missing")),
  value = c(1, 1, 1, 1)  # Dummy values for equal slices
)

dummy_pie <- ggplot(dummy_pie_data) +
  geom_bar(
    aes(x = "", y = value, fill = Type),
    stat = "identity",
    width = 1
  ) +
  coord_polar("y") +
  scale_fill_manual(
    values = c(
      "Single" = "steelblue",
      "Duplicated" = "orange",
      "Fragmented" = "purple",
      "Missing" = "red"
    ),
    name = "BUSCO"
  ) +
  theme_void() +
  theme(legend.position = "right")  # Show the legend

print("Check 2")

# Create variable for standarization of legends format
legends_theme <- theme(legend.title = element_text(size = 10, face = "bold"),
                       legend.text = element_text(size = 8),
                       legend.key.size = unit(0.5, "cm"))

# Extract legend of dummy BUSCO pie chart
legend_busco <- cowplot::get_legend(dummy_pie + legends_theme)

print("Check 2.5")

# Plot BUSCO pies next to tree tips
p <- ggtree::inset(tree_plot,
                   pies,
                   width=(args$pie_size * .16 + .16),
                   height=(args$pie_size * .6 + .6),
                   hjust=-1.2)

print("Check 3")

# Plot Quast data
p2 <- p + geom_fruit(data=data_quast,
                     geom = geom_bar,
                     mapping = aes(y=species, x=value, fill=metric),
                     position=position_dodgex(),
                     pwidth=0.38,
                     orientation="y",
                     stat="identity",
                     axis.params = list(
                                        axis = "x",
                                        text.size  = 1.8,
                                        hjust = 1,
                                        vjust = 0.5,
                                        nbreak = 3
                                        ),
                     grid.params = list(),
                     offset = 1,
                     width = 0.4) +
                     labs(fill = "Quast") +
                     scale_x_continuous(limits=c(0,5)) # limits=c(0,5) by default for p2

print("Check 4")

#
p2 <- increase_scale_until_no_warnings(p2)

# Save plot without legend
ggsave("Phyloplot_no_legend.pdf", p2 + theme(legend.position="none"))

# Extract Quast legend
legend_quast <- cowplot::get_legend(p2 + legends_theme)

# Wrap the legends in fixed-sized containers and adjust margins so that legends
# are perfectly aligned
legend_busco <- ggdraw(legend_busco) + theme(plot.margin = margin(0, 0, 0, 0))
legend_quast <- ggdraw(legend_quast) + theme(plot.margin = margin(0, 29, 60, 0))

# Combine both BUSCO and Quast legends
combined_legends <- plot_grid(legend_busco,
                              legend_quast,
                              NULL, # Dummy rows for better positioning
                              NULL,
                              nrow = 4,
                              rel_widths = c(1, 1))

# Save legend
ggsave("legend.pdf", plot = combined_legends)

# Add combined legends to tree plot
p3 <- plot_grid(p2 + theme(legend.position = "none"),
                combined_legends,
                ncol = 2,
                rel_widths = c(4, 1))

# Save full plot
ggsave("Phyloplot_complete.pdf", plot = p3)
