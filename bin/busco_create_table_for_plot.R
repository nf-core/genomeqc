#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript match_busco_gff.R <busco_file> <gff_file> <output_file>")
}

busco_file <- args[1]
gff_file <- args[2]
output_file <- args[3]

# Step 1: Read the BUSCO file line-by-line, filter out comment and "Missing" lines
busco_raw <- readLines(busco_file)
busco_filtered <- busco_raw[!grepl("^#|Missing", busco_raw)]

# Step 2: Parse the remaining lines as a TSV without column names, then rename columns
busco_data <- read_delim(
  I(busco_filtered),
  delim = "\t",
  col_names = FALSE,
  show_col_types = FALSE
)

# Check if the expected 7 columns are present
if (ncol(busco_data) != 7) {
  stop("Expected 7 columns in BUSCO data after filtering, but found ", ncol(busco_data), ". Please check the input file format.")
}

# Rename columns
colnames(busco_data) <- c("Busco_id", "Status", "Sequence", "Score", "Length", "OrthoDB_url", "Description")

# Read the GFF file
gff_data <- read_tsv(
  gff_file, 
  comment = "#", 
  col_names = FALSE, 
  col_types = cols(
    X1 = col_character(), X2 = col_character(), X3 = col_character(),
    X4 = col_integer(), X5 = col_integer(), X6 = col_character(),
    X7 = col_character(), X8 = col_character(), X9 = col_character()
  ),
  show_col_types = FALSE,
  skip_empty_rows = TRUE
)

# Extract the gene name from the 9th column in GFF, looking for ID=<value> up to the first ;
gff_data <- gff_data %>%
  mutate(gene_name = str_extract(X9, "ID=([^;]+)")) %>%
  mutate(gene_name = str_replace(gene_name, "ID=", "")) %>%  # Remove the "ID=" prefix
  filter(!is.na(gene_name))

# Perform the join on gene name from both data frames
result <- inner_join(
  busco_data,
  gff_data,
  by = c("Sequence" = "gene_name")
)

# Select and rename the columns we need
output_data <- result %>%
  select(Status, Scaffold = X1, Start = X4, End = X5) %>%
  distinct()  # Remove any potential duplicates

# Write the output in the requested format
write.table(output_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Print a message to confirm the output has been written
cat("Output has been written to", output_file, "\n")