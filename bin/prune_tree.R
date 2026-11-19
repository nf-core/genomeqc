#!/usr/bin/env Rscript

# Written by Chris Wyatt and released under the MIT license.
# Prepares a user-supplied species tree (--species_tree) for the tree summary
# plot: keeps only the samples in this run, so the tree's tips match the QC
# tables exactly, and stops with a clear message if a sample is not in the tree.
#
# Usage: prune_tree.R <tree.nwk> <samples.txt> <out.nw>
#   samples.txt: one sample id per line (the samplesheet's `assembly` column)

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("usage: prune_tree.R <tree.nwk> <samples.txt> <out.nw>")
}
tree_file <- args[1]

tree <- tryCatch(suppressWarnings(read.tree(tree_file)), error = function(e) NULL)
if (is.null(tree)) {
    stop("--species_tree: could not read a Newick tree from '", basename(tree_file), "' ",
         "(a Newick tree looks like ((A,B),C); and must end with a semicolon)")
}
if (inherits(tree, "multiPhylo")) {
    stop("--species_tree: '", basename(tree_file), "' contains ", length(tree), " trees; give a file with one tree")
}

tree$tip.label <- trimws(tree$tip.label)
samples <- unique(trimws(readLines(args[2])))
samples <- samples[samples != ""]

dup <- unique(tree$tip.label[duplicated(tree$tip.label)])
if (length(dup) > 0) {
    stop("--species_tree: these tip labels appear more than once: ", paste(dup, collapse = ", "))
}

missing <- setdiff(samples, tree$tip.label)
if (length(missing) > 0) {
    stop(
        "--species_tree: these samples are not tips in the tree: ", paste(missing, collapse = ", "), ".\n",
        "Tip labels must match the samplesheet `assembly` column exactly. ",
        "The tree's tips are: ", paste(sort(tree$tip.label), collapse = ", ")
    )
}

extra <- setdiff(tree$tip.label, samples)
if (length(extra) > 0) {
    message("Dropping ", length(extra), " tip(s) not in this run: ", paste(sort(extra), collapse = ", "))
    tree <- drop.tip(tree, extra)
}

write.tree(tree, args[3])
message("Species tree with ", length(tree$tip.label), " tips written to ", args[3])
