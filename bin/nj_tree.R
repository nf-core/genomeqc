#!/usr/bin/env Rscript

# Written by Chris Wyatt and released under the MIT license.
# Builds a neighbour-joining tree from a `mash triangle` distance matrix
# (--tree_method mash) and roots it at its midpoint, for the tree summary plot.
#
# Usage: nj_tree.R <mash_distances.tsv> <out.nwk>
#   The input is mash triangle's lower-triangle (PHYLIP-like) output: the number
#   of genomes on the first line, then one line per genome with its name and
#   its distances to the genomes listed before it.

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("usage: nj_tree.R <mash_distances.tsv> <out.nwk>")
}

# Read the lower-triangle matrix
lines <- readLines(args[1])
lines <- lines[nzchar(trimws(lines))]
n <- as.integer(trimws(lines[1]))
rows <- strsplit(lines[-1], "\t")
if (length(rows) != n) {
    stop("expected ", n, " genomes in ", args[1], " but found ", length(rows))
}
labels <- vapply(rows, function(r) basename(r[1]), character(1))
d <- matrix(0, n, n, dimnames = list(labels, labels))
for (i in seq_len(n)) {
    vals <- as.numeric(rows[[i]][-1])
    if (length(vals) != i - 1) stop("row ", i, " of ", args[1], " has ", length(vals), " distances, expected ", i - 1)
    if (i > 1) {
        d[i, seq_len(i - 1)] <- vals
        d[seq_len(i - 1), i] <- vals
    }
}
if (n < 3) stop("a tree needs at least 3 genomes; got ", n)

# Neighbour joining; NJ can give slightly negative branch lengths, set those to 0
tree <- nj(as.dist(d))
tree$edge.length[tree$edge.length < 0] <- 0

# Root at the midpoint of the longest tip-to-tip path
midpoint_root <- function(tr) {
    dn <- dist.nodes(tr)
    ntip <- length(tr$tip.label)
    tips <- dn[seq_len(ntip), seq_len(ntip)]
    if (max(tips) == 0) return(tr)                       # all genomes identical: leave as is
    ab <- which(tips == max(tips), arr.ind = TRUE)[1, ]
    a <- ab[[1]]; b <- ab[[2]]
    half <- tips[a, b] / 2
    path <- nodepath(tr, a, b)                           # nodes from tip a to tip b
    from_a <- dn[a, path]
    k <- max(which(from_a <= half))                      # midpoint lies between path[k] and path[k + 1]
    u <- path[k]; v <- path[k + 1]
    # The edge joining u and v, and which end is the child
    e <- which((tr$edge[, 1] == u & tr$edge[, 2] == v) | (tr$edge[, 1] == v & tr$edge[, 2] == u))
    child <- tr$edge[e, 2]
    len <- tr$edge.length[e]
    to_child <- abs(dn[a, child] - half)                 # distance from the midpoint down to the child node
    out <- if (child <= ntip) tr$tip.label[child] else extract.clade(tr, child)$tip.label
    rooted <- root(tr, outgroup = out, resolve.root = TRUE, edgelabel = FALSE)
    # Split the original edge at the midpoint between the two root branches
    root_node <- ntip + 1
    root_edges <- which(rooted$edge[, 1] == root_node)
    clade_of <- function(node) if (node <= ntip) rooted$tip.label[node] else extract.clade(rooted, node)$tip.label
    for (re in root_edges) {
        if (setequal(clade_of(rooted$edge[re, 2]), out)) {
            rooted$edge.length[re] <- to_child
        } else {
            rooted$edge.length[re] <- len - to_child
        }
    }
    rooted
}
tree <- ladderize(midpoint_root(tree))

write.tree(tree, args[2])
message("Neighbour-joining tree of ", n, " genomes written to ", args[2])
