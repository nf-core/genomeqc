#!/usr/bin/env Rscript

# Script adapted from: https://gitlab.com/ezlab/busco_protocol/-/blob/main/support_protocol2/plot_markers/scripts/plot_markers2.R?ref_type=heads

# Load required libraries
library(optparse)
library(RIdeogram)
library(ggplot2)

# Parse command line arguments
option_list <- list(
    make_option(c("--busco_output"), type="character", default=NULL, help="Path to BUSCO output file", metavar="file"),
    make_option(c("--karyotype"), type="character", default=NULL, help="Path to karyotype file", metavar="file"),
    make_option(c("--prefix"), type="character", default="output", help="Prefix for output files", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Load BUSCO output data
busco_mappings <- read.table(opt$busco_output, header = FALSE, stringsAsFactors = FALSE)
colnames(busco_mappings) <- c("Status", "Chr", "Start", "End")

# Load karyotype file
karyotype <- read.table(opt$karyotype, header = TRUE, stringsAsFactors = FALSE)

# Process BUSCO mappings
busco_mappings$Type <- busco_mappings$Status
busco_mappings$Shape <- "circle" 

# Change status to color
busco_mappings$Status <- gsub("Complete", '2dacd6', busco_mappings$Status)
busco_mappings$Status <- gsub("Duplicated", '0a0e1a', busco_mappings$Status)
busco_mappings$Status <- gsub("Fragmented", 'ff0000', busco_mappings$Status)  # Added color for Fragmented
colnames(busco_mappings)[1] <- "color"

busco_mappings <- busco_mappings[, c("Type", "Shape", "Chr", "Start", "End", "color")]

head(busco_mappings)

# Ensure karyotype has required columns
required_columns <- c("Chr", "Start", "End")
if (!all(required_columns %in% colnames(karyotype))) {
    stop("Karyotype file must contain columns: Chr, Start, End")
}

# Use only required columns from karyotype
karyotype <- karyotype[, required_columns]

# Create a vector with all the chromosomes that contain markers
chr_w_markers = unique(sort(busco_mappings$Chr))

# Plot only those chromosomes in which markers where found
filtered_karyotype = karyotype[karyotype$Chr %in% chr_w_markers, ]

# Generate ideogram
ideogram(karyotype = filtered_karyotype, label = busco_mappings, label_type = "marker", output = paste0(opt$prefix, ".svg"))

# Convert to png
convertSVG(paste0(opt$prefix, ".svg"), file = opt$prefix, device = "png")

cat("Ideogram has been generated and saved as", paste0(opt$prefix, ".svg"), "and", paste0(opt$prefix, ".png"), "\n")
