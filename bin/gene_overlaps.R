#!/usr/bin/env Rscript

# Code written by Chris Wyatt with some editing by ChatGPT

library(dplyr)
library(GenomicRanges)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: overlap_analysis.R <input_gff_file> <output_table>")
}

input_gff_file <- args[1]
output_table <- args[2]
output_count <- args[3]

# Helper function to read a GFF file and convert it into a GRanges object
read_gff_to_granges <- function(gff_file) {
    gff_data <- read.delim(gff_file, header = FALSE, comment.char = "#")
    colnames(gff_data) <- c("seqname", "source", "feature", "start", "end",
                            "score", "strand", "frame", "attribute")
    gff_data <- gff_data %>%
      filter(strand %in% c("-","+")) 

    gr <- GRanges(
        seqnames = gff_data$seqname,
        ranges = IRanges(start = gff_data$start, end = gff_data$end),
        strand = gff_data$strand,
        feature = gff_data$feature,
        attribute = gff_data$attribute
    )
    return(gr)
}

# Load GFF file into GRanges
gr <- read_gff_to_granges(input_gff_file)

# Filter only "gene" features
genes <- gr[gr$feature == "gene"]
total_genes <- length(genes)

# Find overlaps
overlap_results <- findOverlaps(genes, genes, ignore.strand = FALSE)

# Initialize counts
sense_count_within <- 0
antisense_count_within <- 0

# Initialize a results data frame
results <- data.frame(
    query_gene = character(),
    subject_gene = character(),
    query_start = integer(),
    query_end = integer(),
    subject_start = integer(),
    subject_end = integer(),
    query_strand = character(),
    subject_strand = character(),
    overlap_type = character(),
    query_overlap_pct = numeric(),
    subject_overlap_pct = numeric(),
    stringsAsFactors = FALSE
)

# Iterate through overlaps
for (i in seq_len(length(overlap_results))) {
    query_idx <- queryHits(overlap_results)[i]
    subject_idx <- subjectHits(overlap_results)[i]

    if (query_idx != subject_idx) {
        query_gene <- genes[query_idx]
        subject_gene <- genes[subject_idx]
        overlap_range <- intersect(ranges(query_gene), ranges(subject_gene))
        overlap_length <- width(overlap_range)

        query_length <- width(ranges(query_gene))
        subject_length <- width(ranges(subject_gene))

        query_overlap_pct <- (overlap_length / query_length) * 100
        subject_overlap_pct <- (overlap_length / subject_length) * 100

        # Determine overlap type
        query_strand <- as.character(strand(query_gene))
        subject_strand <- as.character(strand(subject_gene))

        if (!is.na(query_strand) && !is.na(subject_strand)) {
            if (query_strand == subject_strand) {
                overlap_type <- "sense"
            } else {
                overlap_type <- "antisense"
            }
        } else {
            overlap_type <- "unknown"
        }


        # Increment counters for fully contained overlaps
        if (query_overlap_pct == 100 && overlap_type == "sense") {
            sense_count_within <- sense_count_within + 1
        } else if (query_overlap_pct == 100 && overlap_type == "antisense") {
            antisense_count_within <- antisense_count_within + 1
        }

        # Append to results
        results <- rbind(
            results,
            data.frame(
                query_gene = query_gene$attribute,
                subject_gene = subject_gene$attribute,
                query_start = start(query_gene),
                query_end = end(query_gene),
                subject_start = start(subject_gene),
                subject_end = end(subject_gene),
                query_strand = as.character(strand(query_gene)),
                subject_strand = as.character(strand(subject_gene)),
                overlap_type = overlap_type,
                query_overlap_pct = query_overlap_pct,
                subject_overlap_pct = subject_overlap_pct,
                stringsAsFactors = FALSE
            )
        )
    }
}

# Write results to output file
write.table(results, file = output_table, sep = "\t", row.names = FALSE, quote = FALSE)

# Print the summary
cat("Number of genes fully contained in sense direction:", sense_count_within, "\n")
cat("Number of genes fully contained in antisense direction:", antisense_count_within, "\n")
cat("Overlap analysis complete. Results saved to:", output_table, "\n")

# Calculate total number of overlapping genes
total_overlapping_genes <- length(unique(c(results$query_gene, results$subject_gene)))

# Create summary statistics table
summary_stats <- data.frame(
    Statistic = c("Total number of genes",
                  "Number of genes fully contained in sense direction",
                  "Number of genes fully contained in antisense direction",
                  "Total number of overlapping genes"),
    Count = c(total_genes,
              sense_count_within,
              antisense_count_within,
              total_overlapping_genes)
)

# Write summary statistics table
write.table(summary_stats, file=output_count, sep="\t", quote=FALSE, row.names=FALSE)

# Print the summary
cat("Total number of genes:", total_genes, "\n")
cat("Number of genes fully contained in sense direction:", sense_count_within, "\n")
cat("Number of genes fully contained in antisense direction:", antisense_count_within, "\n")
cat("Total number of overlapping genes:", total_overlapping_genes, "\n")
cat("Overlap analysis complete. Results saved to:", output_table, "\n")
cat("Summary statistics saved to:", output_count, "\n")
