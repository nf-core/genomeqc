
# load libraries
require(RIdeogram)

# From the Rstudio console set working directory "plot_markers/" folder:
# setwd("/path/to/support_protocol2/plot_markers/")


karyotype_file <- "./data/karyotype/Dmel_karyotype.txt"  # provide path to karyotype file 
coordinates <- "./plot_results/Dmel_busco_coordinates.txt" # provide path to markers' coordinates
out_name <- "Dmel" # sample name

# load karyotype
karyotype <- read.csv(karyotype_file, sep="\t", header = TRUE, stringsAsFactors = F)

# load mapping busco_coordinates.txt
busco_mappings <- read.csv(coordinates, sep="\t", header = FALSE, stringsAsFactors = F)

colnames(busco_mappings)[1] <- "Status"
colnames(busco_mappings)[2] <- "Chr"
colnames(busco_mappings)[3] <- "Start"
colnames(busco_mappings)[4] <- "End"

busco_mappings$Type <- "BUSCO_marker"
busco_mappings$Shape <- "circle" 

# change status in color
busco_mappings$Status <- gsub("Complete", '2dacd6', busco_mappings$Status)
busco_mappings$Status <- gsub("Duplicated", '0a0e1a', busco_mappings$Status)
busco_mappings$Status <- gsub("Fragmented", 'eded13', busco_mappings$Status)
colnames(busco_mappings)[1] <- "color"

busco_mappings <- busco_mappings[, c(5,6,2,3,4,1)]

ideogram(karyotype = karyotype, label = busco_mappings, label_type = "marker", output = paste0(out_name, ".svg"))
# convert to png
convertSVG(paste0(out_name, ".svg"), device = "png")

