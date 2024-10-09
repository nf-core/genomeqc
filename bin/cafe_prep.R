#!/usr/bin/Rscript
library(ape)
library(data.table)

tre <- read.tree('pruned_tree')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre)
}
write.tree(utre, 'SpeciesTree_rooted_ultra.txt')

hog <- fread('N0.tsv')
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 1]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')
