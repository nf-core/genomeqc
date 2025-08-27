#!/usr/bin/env Rscript

# Written by Chris Wyatt and Fernando Duarte and released under the MIT license.
# Plots the phylogenetic tree with BUSCO, Quast and gene stats results

# Functions mirroring nf-core/genomeqc (gff_valid) plot_tree_summary.R behavior
# Same function names as your prior modular version

# -- Packages are expected to be loaded by the caller (keep as comments if you prefer) --
# library(ggtree)
# library(ggplot2)
# library(patchwork)
# library(dplyr)
# library(tidyr)
# library(scatterpie)
# library(scales)
# library(cowplot)
# library(ape)

# Avoid scientific notation in all plots
options(scipen = 999)

# ---------------------------
# Helper: avoid mistmatches between species name in data frame and tree tips
# ---------------------------
check_match <- function(data, col, reference) {
  vals <- data[[col]]
  if (!setequal(vals, reference)) {
    stop("Values in column `", col, "` do not match reference")
  }
  data
}

# ---------------------------
# Helper: safely extract axis ranges (unchanged name)
# ---------------------------
get_plot_range <- function(plot, axis = "y") {
  plot_name <- deparse(substitute(plot))
  tryCatch({
    built <- ggplot_build(plot)
    if (axis == "y") {
      return(built$layout$panel_scales_y[[1]]$range$range)
    } else if (axis == "x") {
      return(built$layout$panel_scales_x[[1]]$range$range)
    } else {
      stop("Invalid axis specified.")
    }
  }, error = function(e) {
    warning("Failed to load plot ", plot_name, ": ", conditionMessage(e))
    numeric(0)
  })
}

# ---------------------------
# Data loaders (same names)
# ---------------------------
load_busco <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    data_busco <- read.csv(file, sep = "\t", colClasses = c("Input_file" = "character"))
    data_busco <- data_busco %>%
      mutate(Input_file = tools::file_path_sans_ext(Input_file)) %>%
      rename(species = Input_file) %>%
      check_match("species", tree_tips) %>%
      arrange(match(species, tree_tips)) %>%
      mutate(node = 1:length(species))
    data_busco
  }, error = function(e) {
    warning("Failed to load BUSCO file: ", conditionMessage(e))
    NULL
  })
}

load_quast <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    data_quast <- read.csv(file, sep = "\t")
    colnames(data_quast)[5] <- "GC"
    colnames(data_quast)[6] <- "Sequences"
    data_quast <- data_quast %>%
      filter(!is.na(species)) %>%
      filter(N50 != "bar") %>%
      mutate(Total.length = (as.numeric(Total.length)/1000000)) %>%   # Mb
      mutate(Sequences = as.integer(Sequences)) %>%
      mutate(GC = as.numeric(Total.length) * as.numeric(GC) / 100) %>%
      rename(Length = Total.length) %>%
      check_match("species", tree_tips) %>%
      arrange(match(species, tree_tips)) %>%
      mutate(node = 1:length(species))

    n5090 <- data_quast %>%
      pivot_longer(cols = c(N50, N90), names_to = "metric", values_to = "value") %>%
      mutate(value = as.numeric(value)/1000000)                      # Mb

    len <- data_quast %>%
      pivot_longer(cols = c(GC, Length), names_to = "metric", values_to = "value")

    list(full = data_quast, n5090 = n5090, len = len)
  }, error = function(e) {
    warning("Failed to load Quast file: ", conditionMessage(e))
    NULL
  })
}

load_genes <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    data_genes <- read.csv(file, sep = "\t")
    data_genes <- data_genes %>%
      rename(species = File, Total = Total_genes, Overlapping = Overlapping_genes) %>%
      mutate(species = gsub("\\.counts\\.tsv", "", species)) %>%
      check_match("species", tree_tips) %>%
      arrange(match(species, tree_tips)) %>%
      mutate(node = 1:length(species)) %>%
      pivot_longer(cols = c(Total, Overlapping), names_to = "stat", values_to = "value")
    data_genes
  }, error = function(e) {
    warning("Failed to load gene stats file: ", conditionMessage(e), ". Ignore this message if there's no gene stats file (e.g. genome only).")
    NULL
  })
}

load_nseqs <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    data_nseqs <- read.csv(file, sep = "\t")
    #print(data_nseqs[,1])
    print("nseqs in")
    print(tree_tips)
    print(data_nseqs[,1])
    print("nseqs out")
    data_nseqs <- data_nseqs %>%
      rename(species = 1) %>%
      check_match("species", tree_tips) %>% # Ensure species match tree tips
      arrange(match(species, tree_tips)) %>% # What is this match here for then?
      mutate(node = 1:n())
    data_nseqs
  }, error = function(e) {
    warning("Failed to load nseqs file: ", conditionMessage(e), ". This will always fail once (loads both ortho and nseq, but one of them will always be empty). Worry if it fails more than once")
    NULL
  })
}

# ---------------------------
# Main data ingest (same name)
# ---------------------------
process_tree_data <- function(tree_file, busco_file = NULL, quast_file = NULL,
                              genes_file = NULL, nseqs_file = NULL, ortho_file = NULL) {

  tree <- read.tree(tree_file)
  tree$tip.label <- trimws(tree$tip.label)

  tree_plot <- ggtree(tree)               # temp plot for get_taxa_name()
  tips_order <- rev(get_taxa_name(tree_plot))

  data_busco <- load_busco(busco_file, tips_order)
  data_quast <- load_quast(quast_file, tips_order)
  data_genes <- load_genes(genes_file, tips_order)
  data_nseqs <- load_nseqs(nseqs_file, tips_order)
  data_ortho <- load_nseqs(ortho_file, tips_order)

  tree$tip.label <- gsub("_", " ", tree$tip.label)  # matches gff_valid behavior

  list(
    tree = tree,
    tips_order = tips_order,
    data_busco = data_busco,
    data_quast = data_quast,
    data_genes = data_genes,
    data_nseqs = data_nseqs,
    data_ortho = data_ortho
  )
}

# ---------------------------
# Plot generator (same name)
# ---------------------------
generate_plots <- function(processed_data, text_size = 3, tree_scale = 0.0005,
                           bar_width = 0.7, rad_width = 0.4, skip_stats = NULL) {

  tree <- processed_data$tree
  data_busco <- processed_data$data_busco
  data_quast <- processed_data$data_quast
  data_genes <- processed_data$data_genes
  data_nseqs <- processed_data$data_nseqs
  data_ortho <- processed_data$data_ortho

  # Small trees: reduce bar/pie sizes (as in gff_valid)
  if (length(tree$tip.label) < 7) {
    bar_width <- bar_width/2
    rad_width <- rad_width/3
  }

  # Theme used in gff_valid
  barplots_theme <- theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 6),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(),
      axis.line.y = element_blank()
    )

  # --- Plots ---
  nseqs_plot <- if (!is.null(data_nseqs)) {
    ggplot(data_nseqs, aes(x = 1, y = node)) +
      geom_text(aes(label = data_nseqs[, 2])) +
      theme_void() +
      ggtitle("Seqs ≥5\nBUSCOs") +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))
  } else NULL

  ortho_plot <- if (!is.null(data_ortho)) {
    ggplot(data_ortho, aes(x = 1, y = node)) +
      geom_text(aes(label = data_ortho[, 2])) +
      theme_void() +
      ggtitle("Ortho \nSeqs") +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))
  } else NULL

  ch_plot <- NULL
  len_plot <- NULL
  n50_plot <- NULL
  legend_len <- NULL

  if (!is.null(data_quast)) {
    ch_plot <- ggplot(data_quast$full, aes(x = 1, y = node)) +
      geom_text(aes(label = Sequences)) +
      theme_void() +
      ggtitle("Seq\nNumber") +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -2.2))

    len_plot <- ggplot(data_quast$len, aes(y = value, x = node)) +
      geom_col(aes(fill = metric), position = position_stack(reverse = TRUE), width = bar_width) +
      scale_fill_manual(labels = c("GC %", "Length"), values = c("brown1", "cornflowerblue")) +
      ggtitle("Genome\nsize (Mb)") +
      barplots_theme +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -5)) +
      coord_flip() +
      xlab(NULL) + ylab(NULL)

    legend_len <- cowplot::get_legend(
      len_plot +
        theme(legend.position = "right",
              legend.justification = c(0, 1.2),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.background = element_rect(fill = NA),
              legend.text = element_text(size = 8))
    )
    len_plot <- len_plot + guides(fill = "none")

    data_quast_n50 <- data_quast$n5090[data_quast$n5090$metric %in% "N50", ]
    n50_plot <- ggplot(data_quast_n50, aes(y = value, x = node)) +
      geom_col(position = position_stack(reverse = TRUE), width = bar_width, fill = "steelblue") +
      ggtitle("N50 (Mb)") +
      barplots_theme +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -0.4)) +
      coord_flip() +
      xlab(NULL) + ylab(NULL) +
      guides(fill = "none")
  }

  pies_plot <- NULL
  legend_busco <- NULL
  if (!is.null(data_busco)) {
    pies_plot <- ggplot() +
      geom_scatterpie(
        aes(x = 0, y = node, group = species, r = rad_width),
        data = data_busco,
        cols = c("Single", "Duplicated", "Fragmented", "Missing"),
        color = NA
      ) +
      scale_fill_manual(values = c("deepskyblue", "orange", "darkorchid4", "firebrick1")) +
      coord_fixed() +
      theme_void() +
      ggtitle("BUSCO") +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05))

    legend_busco <- cowplot::get_legend(
      pies_plot +
        theme(legend.position = "right",
              legend.justification = c(0, 1.08),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.text = element_text(size = 8))
    )
    pies_plot <- pies_plot + guides(fill = "none")
  }

  gene_plot <- NULL
  legend_gene <- NULL
  if (!is.null(data_genes)) {
    gene_plot <- ggplot(data_genes, aes(y = value, x = node)) +
      geom_col(aes(fill = stat), position = position_stack(reverse = TRUE), width = bar_width) +
      scale_fill_manual(values = c("indianred1", "lightsteelblue")) +
      ggtitle("Gene\nnumber") +
      barplots_theme +
      theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -5)) +
      coord_flip() +
      scale_y_continuous(breaks = pretty_breaks(n = 3)) +
      xlab(NULL) + ylab(NULL)

    legend_gene <- cowplot::get_legend(
      gene_plot +
        theme(legend.position = "right",
              legend.justification = c(0, 1.2),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              legend.background = element_rect(fill = NA),
              legend.text = element_text(size = 8))
    )
    gene_plot <- gene_plot + guides(fill = "none")
  }

  # Build the tree base, as in gff_valid
  tree_plot <- ggtree(tree) +
    geom_tiplab(size = text_size, fontface = "italic", align = TRUE, hjust = -0.05) +
    theme(plot.margin = margin(10, 30, 10, 10)) +
    coord_cartesian(clip = "off")

  # Collect ranges (y for non-flipped, x for flipped)
  all_ranges <- c(
    get_plot_range(ch_plot, "y"),
    get_plot_range(nseqs_plot, "y"),
    get_plot_range(ortho_plot, "y"),
    get_plot_range(pies_plot, "y"),
    get_plot_range(len_plot, "x"),
    get_plot_range(n50_plot, "x"),
    get_plot_range(gene_plot, "x")
  )

  if (length(all_ranges) > 0) {
    new_ylim <- ylim(c(min(all_ranges, na.rm = TRUE), max(all_ranges, na.rm = TRUE)))
    new_xlim <- xlim(c(min(all_ranges, na.rm = TRUE), max(all_ranges, na.rm = TRUE)))

    if (!is.null(nseqs_plot)) nseqs_plot <- nseqs_plot + new_ylim
    if (!is.null(ortho_plot)) ortho_plot <- ortho_plot + new_ylim
    if (!is.null(ch_plot))    ch_plot    <- ch_plot    + new_ylim
    if (!is.null(pies_plot))  pies_plot  <- pies_plot  + new_ylim

    if (!is.null(len_plot))   len_plot   <- len_plot   + new_xlim
    if (!is.null(n50_plot))   n50_plot   <- n50_plot   + new_xlim
    if (!is.null(gene_plot))  gene_plot  <- gene_plot  + new_xlim

    tree_plot <- tree_plot + new_ylim
  } else {
    new_xlim <- xlim(0, 10)
  }

  # Build lists, then apply skipping exactly like the script’s intent
  all_plots <- list(
    ch_plot   = ch_plot,
    nseqs_plot = nseqs_plot,
    ortho_plot = ortho_plot,
    len_plot   = len_plot,
    gene_plot  = gene_plot,
    n50_plot   = n50_plot,
    pies_plot  = pies_plot
  )
  all_legends <- list(
    ch_plot   = NULL,
    nseqs_plot = NULL,
    ortho_plot = NULL,
    len_plot   = legend_len,
    gene_plot  = legend_gene,
    n50_plot   = NULL,
    pies_plot  = legend_busco
  )

  if (!is.null(skip_stats) && length(skip_stats) > 0) {
    plots   <- all_plots[!names(all_plots) %in% skip_stats & !sapply(all_plots, is.null)]
    legends <- all_legends[names(plots)]
  } else {
    plots   <- all_plots[!sapply(all_plots, is.null)]
    legends <- all_legends[names(plots)]
  }

  list(
    tree_plot = tree_plot,
    plots = plots,
    legends = legends,
    xlim = new_xlim
  )
}

# ---------------------------
# Layout combiner (same name)
# ---------------------------
build_tree_plot <- function(tree, plots, legends, xlimit, top_margin = 5.5,
                            right_margin = 5.5, bottom_margin = 5.5, # 5.5 is the default value for margins
                            left_margin = 5.5, tree_space_ratio = 1.3) {

  # Dynamic tree width exactly as in gff_valid's "improved" version
  tree_built <- ggplot_build(tree)
  tree_data  <- tree_built$data[[1]]

  max_x <- max(tree_data$x, na.rm = TRUE)
  tip_labels <- tree$data$label[!is.na(tree$data$label)]
  max_label_chars <- max(nchar(tip_labels), na.rm = TRUE)

  text_size_pts <- tree$theme$text$size %||% 11
  char_width_estimate <- text_size_pts * 0.015
  label_padding <- max_label_chars * char_width_estimate

  tree_xlim <- max_x * tree_space_ratio + label_padding
  tree <- tree + xlim(0, tree_xlim)

  # Build the patchwork layout
  combined_plots <- tree
  n_plots <- length(plots)
  tree_width <- max(5, 10 - n_plots)
  widths <- c(tree_width, rep(1, n_plots))

  for (p in plots) combined_plots <- combined_plots | p

  combined_legends <- plot_spacer() + xlimit
  for (lg in legends) {
    legend_plot <- if (!is.null(lg) && length(lg) != 0) lg else plot_spacer() + xlimit
    combined_legends <- combined_legends | legend_plot
  }

  combined_plots <- combined_plots + plot_layout(widths = widths)
  combined_legends <- combined_legends +
    plot_layout(widths = widths) +
    theme(plot.margin = margin(0, 15, 60, 0))

  final_plot <- combined_plots / combined_legends + plot_layout(heights = c(0.99, 0.01)) +
                  plot_annotation(theme = theme( # Had to use plot_annotation to set margins as + theme doesn't work with patchwork
                                  plot.margin = margin(top_margin, right_margin, bottom_margin, left_margin)
                                  ))

  final_plot
}

# ---------------------------
# Final wrapper (same name)
# ---------------------------
generate_complete_plot <- function(processed_data, text_size = 3, tree_scale = 0.0005,
                                   bar_width = 0.7, rad_width = 0.4, skip_stats = NULL,
                                   top_margin = 5.5, right_margin = 5.5, bottom_margin = 5.5,
                                   left_margin = 5.5, tree_margin = 15, tree_space_ratio = 1.3) {

  plot_results <- generate_plots(
    processed_data = processed_data,
    text_size = text_size,
    tree_scale = tree_scale,
    bar_width = bar_width,
    rad_width = rad_width,
    skip_stats = skip_stats
  )

  if (length(plot_results$plots) > 0) {
    build_tree_plot(
      tree     = plot_results$tree_plot,
      plots    = plot_results$plots,
      legends  = plot_results$legends,
      xlimit   = plot_results$xlim,
      top_margin = top_margin,
      right_margin = right_margin,
      bottom_margin = bottom_margin,
      left_margin = left_margin,
      tree_space_ratio = tree_space_ratio
    )
  } else {
    plot_results$tree_plot
  }
}
