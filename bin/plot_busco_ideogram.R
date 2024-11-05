#!/usr/bin/env Rscript

# Script adapted from: https://gitlab.com/ezlab/busco_protocol/-/blob/main/support_protocol2/plot_markers/scripts/plot_markers2.R?ref_type=heads

library(optparse)
library(RIdeogram)

# Parse command line arguments
option_list <- list(
    make_option(c("--busco_full_table"), type="character", default=NULL, help="Path to BUSCO full table", metavar="file"),
    make_option(c("--prefix"), type="character", default="output", help="Prefix for output files", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Load BUSCO full table and process it
busco_full <- read.csv(opt$busco_full_table, sep="\t", header = FALSE, stringsAsFactors = F)
busco_mappings <- busco_full[busco_full$V2 %in% c("Complete", "Duplicated"), c(2,3,4,5)]
colnames(busco_mappings) <- c("Status", "Chr", "Start", "End")

# Generate karyotype from BUSCO data
karyotype <- data.frame(
  Chr = unique(busco_mappings$Chr),
  Start = 0,
  End = sapply(unique(busco_mappings$Chr), function(chr) max(busco_mappings$End[busco_mappings$Chr == chr]))
)

busco_mappings$Type <- "BUSCO_marker"
busco_mappings$Shape <- "circle" 

# Change status to color
busco_mappings$Status <- gsub("Complete", '2dacd6', busco_mappings$Status)
busco_mappings$Status <- gsub("Duplicated", '0a0e1a', busco_mappings$Status)
colnames(busco_mappings)[1] <- "color"

busco_mappings <- busco_mappings[, c(5,6,2,3,4,1)]

# Generate ideogram
ideogram(karyotype = karyotype, label = busco_mappings, label_type = "marker", output = paste0(opt$prefix, ".svg"))

# Convert to png
convertSVG(paste0(opt$prefix, ".svg"), device = "png")