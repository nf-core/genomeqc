#!/usr/bin/env Rscript

# load libraries
library(RIdeogram)

# Script adapted from: https://gitlab.com/ezlab/busco_protocol/-/blob/main/support_protocol2/plot_markers/scripts/plot_markers2.R?ref_type=heads

args <- commandArgs(trailingOnly = TRUE)
file_name_for_karyotype <- args[1]
out_name <- args[2] # how to name the karyotype
 
##### karyotype-------
karyotype1 <- read.csv(file=file_name_for_karyotype, sep="\t", header = FALSE, stringsAsFactors = F)
karyotype1$start <- 1
karyotype1$genome <- out_name
karyotype1$size <- 12
karyotype1$color <- "25252"

# fix columns names
colnames(karyotype1) <- c("Chr", "End", "Start", "species", "size", "color")
# reorder columns
karyotype <- karyotype1[, c(1,3,2,4,5,6)]

write.table(karyotype, file = paste0(out_name, "_karyotype.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


