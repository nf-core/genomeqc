#!/usr/bin/env Rscript

# Code written by Chris Wyatt with some editing by ChatGPT
#!/usr/bin/env Rscript

library(GenomicRanges)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: overlap_analysis.R <input_gff_file> <output_table>")
}

input_gff_file <- args[1]
output_table <- args[2]

# Helper function to read a GFF file and convert it into a GRanges object
read_gff_to_granges <- function(gff_file) {
    gff_data <- read.delim(gff_file, header = FALSE, comment.char = "#")
    colnames(gff_data) <- c("seqname", "source", "feature", "start", "end",
                            "score", "strand", "frame", "attribute")
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

# Find overlaps
overlap_results <- findOverlaps(genes, genes, ignore.strand = FALSE)

# Initialize a results data frame
results <- data.frame(
    gene1_id = character(),
    gene1_start = integer(),
    gene1_end = integer(),
    gene2_id = character(),
    gene2_start = integer(),
    gene2_end = integer(),
    overlap_percent = numeric(),
    stringsAsFactors = FALSE
)

# Iterate through overlaps
for (i in seq_len(length(overlap_results))) {
    query_idx <- queryHits(overlap_results)[i]
    subject_idx <- subjectHits(overlap_results)[i]

    if (query_idx != subject_idx) {
        gene1 <- genes[query_idx]
        gene2 <- genes[subject_idx]

        overlap_range <- intersect(ranges(gene1), ranges(gene2))
        overlap_length <- width(overlap_range)

        gene1_length <- width(ranges(gene1))
        gene2_length <- width(ranges(gene2))

        overlap_percent <- (overlap_length / gene1_length) * 100

        results <- rbind(
            results,
            data.frame(
                gene1_id = gene1$attribute,
                gene1_start = start(gene1),
                gene1_end = end(gene1),
                gene2_id = gene2$attribute,
                gene2_start = start(gene2),
                gene2_end = end(gene2),
                overlap_percent = overlap_percent,
                stringsAsFactors = FALSE
            )
        )
    }
}

# Write results to output file
write.table(results, file = output_table, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Overlap analysis complete. Results saved to:", output_table, "\n")
