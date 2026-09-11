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

load_te <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    # node = this species' real index in tree_tips, not its row position -
    # TE annotation can be missing for some species (e.g. a species-level
    # failure), and 1:length(species) would silently misassign the row that
    # follows a gap onto the wrong tip. Drop any species tree_tips doesn't
    # recognise (shouldn't happen, but never plot at a fabricated position).
    data_te <- read.csv(file, sep = "\t") %>%
      mutate(node = match(species, tree_tips)) %>%
      filter(!is.na(node))
    data_te
  }, error = function(e) {
    warning("Failed to load TE file: ", conditionMessage(e))
    NULL
  })
}

load_fcs <- function(file, tree_tips) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch({
    # node = this species' real index in tree_tips, not its row position -
    # FCS-GX only runs for samples with a taxid set (a documented, per-sample
    # opt-in), so this table is routinely a subset of all species. Using
    # 1:length(species) here would silently misassign a present species'
    # pie onto a DIFFERENT tip's row rather than just omitting the absent
    # ones. Drop any species tree_tips doesn't recognise (shouldn't happen,
    # but never plot at a fabricated position).
    data_fcs <- read.csv(file, sep = "\t") %>%
      mutate(node = match(species, tree_tips)) %>%
      filter(!is.na(node))
    data_fcs
  }, error = function(e) {
    warning("Failed to load FCS file: ", conditionMessage(e))
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
process_tree_data <- function(tree_file, busco_file_geno = NULL, busco_file_prot = NULL, quast_file = NULL,
                              genes_file = NULL, nseqs_file = NULL, ortho_file = NULL, te_file = NULL, fcs_file = NULL) {

  tree <- read.tree(tree_file)
  tree$tip.label <- trimws(tree$tip.label)

  tree_plot <- ggtree(tree)               # temp plot for get_taxa_name()
  tips_order <- rev(get_taxa_name(tree_plot))

  data_busco_geno <- load_busco(busco_file_geno, tips_order)
  data_busco_prot <- load_busco(busco_file_prot, tips_order)
  data_quast <- load_quast(quast_file, tips_order)
  data_genes <- load_genes(genes_file, tips_order)
  data_nseqs <- load_nseqs(nseqs_file, tips_order)
  data_ortho <- load_nseqs(ortho_file, tips_order)
  data_te    <- load_te(te_file, tips_order)
  data_fcs   <- load_fcs(fcs_file, tips_order)

  tree$tip.label <- gsub("_", " ", tree$tip.label)  # matches gff_valid behavior

  list(
    tree = tree,
    tips_order = tips_order,
    data_busco_geno = data_busco_geno,
    data_busco_prot = data_busco_prot,
    data_quast = data_quast,
    data_genes = data_genes,
    data_nseqs = data_nseqs,
    data_ortho = data_ortho,
    data_te    = data_te,
    data_fcs   = data_fcs
  )
}

# ---------------------------
# Plot generator (same name)
# ---------------------------

# But first let's create a function that generates BUSCO pie plots separately for genome and proteome
plot_busco_pies <- function(data_busco,
                            title = "BUSCO",
                            rad_width = NULL,
                            len_pos_x = 0,
                            fill_colors = c(
                              Single     = "deepskyblue",
                              Duplicated = "orange",
                              Fragmented = "darkorchid4",
                              Missing    = "firebrick1"
                            )) {

  # Default return
  out <- list(
    plot   = NULL,
    legend = NULL
  )

  if (is.null(data_busco)) {
    return(out)
  }

  pies_plot <- ggplot() +
    geom_scatterpie(
      aes(x = 0, y = node, group = species, r = rad_width),
      data = data_busco,
      cols = names(fill_colors),
      color = NA
    ) +
    scale_fill_manual(values = fill_colors) +
    coord_fixed() +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05)
    )

  legend_busco <- cowplot::get_legend(
    pies_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        #legend.box.margin = margin(l = 6)
      )
  )

  pies_plot <- pies_plot + guides(fill = "none")

  out$plot   <- pies_plot
  out$legend <- legend_busco

  return(out)
}

plot_te_bars <- function(data_te, bar_width = NULL, len_pos_x = 0) {
  out <- list(plot = NULL, legend = NULL)
  if (is.null(data_te)) return(out)

  te_colors <- c(
    "SINE"           = "deepskyblue",
    "LINE"           = "orange",
    "LTR"            = "darkorchid4",
    "Penelope"       = "firebrick1",
    "DNA"            = "forestgreen",
    "Rolling_Circle" = "goldenrod",
    "Unclassified"   = "purple",
    "Other"          = "indianred",
    "Non_Repeat"     = "darkgray"
  )

  data_te_long <- data_te %>%
    select(node, all_of(names(te_colors))) %>%
    pivot_longer(cols = all_of(names(te_colors)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(te_colors)))

  te_plot <- ggplot(data_te_long, aes(x = node, y = value, fill = metric)) +
    geom_col(
      position = position_fill(reverse = TRUE), # Rescales each bar to 100% (like a pie), stacked in legend order
      width = bar_width,
      color = NA
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = te_colors) +
    ggtitle("TE") +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 6),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(),
      axis.line.y = element_blank(),
      # vjust=-5 (used by the two-line titles like "Genome\nsize (Mb)") pushes a
      # single-line title like this one below the panel instead of above it -
      # match N50 (Mb)'s single-line calibration instead.
      plot.title = element_text(size = 9, hjust = 0.5, vjust = -0.4)
    ) +
    coord_flip() +
    xlab(NULL) + ylab(NULL)

  legend_te <- cowplot::get_legend(
    te_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8)
      )
  )

  out$plot   <- te_plot + guides(fill = "none")
  out$legend <- legend_te
  out
}

plot_fcs_pie <- function(data_fcs, rad_width = NULL, len_pos_x = 0) {
  out <- list(plot = NULL, legend = NULL)
  if (is.null(data_fcs)) return(out)

  pies_plot <- ggplot() +
    geom_scatterpie(
      aes(x = 0, y = node, group = species, r = rad_width),
      data = data_fcs,
      cols = c("non_contaminant_pct", "contaminant_pct"),
      color = NA
    ) +
    scale_fill_manual(
      labels = c("Non-contaminant", "Contaminant"),
      values = c(
        "non_contaminant_pct" = "deepskyblue",
        "contaminant_pct"     = "firebrick1"
      )
    ) +
    coord_fixed() +
    theme_void() +
    ggtitle("FCS") +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, vjust = 0.05)
    )

  legend_fcs <- cowplot::get_legend(
    pies_plot +
      theme(
        legend.position = "right",
        legend.justification = c(len_pos_x, 1.08),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8)
      )
  )

  out$plot   <- pies_plot + guides(fill = "none")
  out$legend <- legend_fcs
  out
}

# Now for the plot generator
generate_plots <- function(processed_data, text_size = 3, tree_scale = 0.0005,
                           bar_width = 0.7, rad_width = 0.4, skip_stats = NULL, busco_len_pos_x = 0.5,
                           tree_style = "roundrect") {

  tree <- processed_data$tree
  data_busco_geno <- processed_data$data_busco_geno
  data_busco_prot <- processed_data$data_busco_prot
  data_quast <- processed_data$data_quast
  data_genes <- processed_data$data_genes
  data_nseqs <- processed_data$data_nseqs
  data_ortho <- processed_data$data_ortho
  data_te    <- processed_data$data_te
  data_fcs   <- processed_data$data_fcs

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
      ggtitle("Seq\nnumber") +
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

  # BUSCO plots
  # -- if both genome and proteome busco datasets are present,
  # change legend x position so that it's not skewed --
  # len_pos_x <- busco_len_pos_x * (!is.null(data_busco_geno) && !is.null(data_busco_prot)) # very smart chatgpt
  len_pos_x <- busco_len_pos_x
  # Plot both genome and proteome BUSCO pies
  busco_gen_plot  <- plot_busco_pies(data_busco_geno,
                                      rad_width = rad_width,
                                      title = "BUSCO\ngenome",
                                      len_pos_x = len_pos_x)
  busco_prot_plot <- plot_busco_pies(data_busco_prot,
                                      rad_width = rad_width,
                                      title = "BUSCO\nprotein",
                                      len_pos_x = len_pos_x)
  te_result <- plot_te_bars(data_te, bar_width = bar_width, len_pos_x = len_pos_x)
  fcs_result <- plot_fcs_pie(data_fcs, rad_width = rad_width, len_pos_x = len_pos_x)


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

  # Build the tree base according to the selected style
  if (tree_style == "rectangular") {
    # Legacy look: thin black branches with dotted alignment leader lines
    tree_plot <- ggtree(tree) +
      geom_tiplab(size = text_size, fontface = "italic", align = TRUE, hjust = -0.05)
  } else {
    # Modern look: thicker grey branches, aligned labels without dotted leaders
    tree_plot <- ggtree(tree, layout = tree_style, size = 0.7, colour = "grey30") +
      geom_tiplab(size = text_size, fontface = "italic", align = TRUE,
                  linetype = NA, hjust = -0.05)
    if (tree_style == "ellipse") {
      # Subtle node markers to accentuate the curved layout
      tree_plot <- tree_plot + geom_nodepoint(colour = "steelblue", size = 1.2, alpha = 0.75)
    }
  }
  tree_plot <- tree_plot +
    theme(plot.margin = margin(10, 30, 10, 10)) +
    coord_cartesian(clip = "off")

  # Collect ranges (y for non-flipped, x for flipped)
  all_ranges <- c(
    get_plot_range(ch_plot, "y"),
    get_plot_range(nseqs_plot, "y"),
    get_plot_range(ortho_plot, "y"),
    get_plot_range(busco_gen_plot$plot, "y"),
    get_plot_range(busco_prot_plot$plot, "y"),
    get_plot_range(fcs_result$plot, "y"),
    get_plot_range(len_plot, "x"),
    get_plot_range(n50_plot, "x"),
    get_plot_range(gene_plot, "x"),
    get_plot_range(te_result$plot, "x") # te_result is now a coord_flip()'d bar (like len_plot/n50_plot), not a pie - node range lives on x, not y
  )

  if (length(all_ranges) > 0) {
    new_ylim <- ylim(c(min(all_ranges, na.rm = TRUE), max(all_ranges, na.rm = TRUE)))
    new_xlim <- xlim(c(min(all_ranges, na.rm = TRUE), max(all_ranges, na.rm = TRUE)))

    if (!is.null(nseqs_plot)) nseqs_plot <- nseqs_plot + new_ylim
    if (!is.null(ortho_plot)) ortho_plot <- ortho_plot + new_ylim
    if (!is.null(ch_plot))    ch_plot    <- ch_plot    + new_ylim
    if (!is.null(busco_gen_plot$plot))  busco_gen_plot$plot  <- busco_gen_plot$plot  + new_ylim
    if (!is.null(busco_prot_plot$plot))  busco_prot_plot$plot  <- busco_prot_plot$plot  + new_ylim
    if (!is.null(fcs_result$plot)) fcs_result$plot <- fcs_result$plot + new_ylim

    if (!is.null(len_plot))   len_plot   <- len_plot   + new_xlim
    if (!is.null(n50_plot))   n50_plot   <- n50_plot   + new_xlim
    if (!is.null(gene_plot))  gene_plot  <- gene_plot  + new_xlim
    if (!is.null(te_result$plot)) te_result$plot <- te_result$plot + new_xlim # te_result is coord_flip()'d

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
    busco_gen_plot  = busco_gen_plot$plot,
    busco_prot_plot = busco_prot_plot$plot,
    te_plot         = te_result$plot,
    fcs_plot        = fcs_result$plot
  )
  all_legends <- list(
    ch_plot   = NULL,
    nseqs_plot = NULL,
    ortho_plot = NULL,
    len_plot   = legend_len,
    gene_plot  = legend_gene,
    n50_plot   = NULL,
    busco_gen_plot  = busco_gen_plot$legend,
    busco_prot_plot = if (!is.null(busco_gen_plot$legend)) NULL else busco_prot_plot$legend, # Only plot legend once
    te_plot         = te_result$legend,
    fcs_plot        = fcs_result$legend
  )

  # What's this for? To filter out skipped plots and legends?
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
# Quality scoring for the circular layout
# ---------------------------
# Quality metrics are drawn as a discrete "traffic light" instead of a sequential
# ramp, because a light->dark ramp implies "dark = good", which is wrong for
# lower-is-better stats (e.g. sequence count) and meaningless for descriptive
# ones (genome size, gene number). Colour-vision-safe Okabe-Ito triple.
QUALITY_COLOURS <- c(Good = "#009E73", Warn = "#E69F00", Poor = "#D55E00")

# Threshold presets by phylogenetic group.
# NOTE: the BUSCO cut-offs follow community practice; the N50 and sequence-count
# cut-offs are clade-dependent STARTING POINTS and should be tuned per project.
# 'generic' is deliberately lenient.
# fcs_noncontam (FCS-GX non-contaminant %) is the same across every preset: unlike
# BUSCO/N50/sequence-count, "how much of the genome is foreign contamination" has
# no taxon-dependent expectation, so there's no reason to vary it by clade.
QUALITY_PRESETS <- list(
  generic = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  vertebrate = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e7,    warn = 1e6),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  insect = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 1000,   warn = 10000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  plant = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 10,     warn = 20),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 5000,   warn = 50000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  fungi = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 1e6,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 100,    warn = 1000),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  ),
  bacteria = list(
    busco_complete   = list(direction = "higher", good = 95,     warn = 90),
    busco_duplicated = list(direction = "lower",  good = 5,      warn = 10),
    n50              = list(direction = "higher", good = 5e5,    warn = 1e5),
    seq_number       = list(direction = "lower",  good = 10,     warn = 100),
    fcs_noncontam    = list(direction = "higher", good = 99.5,   warn = 98)
  )
)

# Bin values into Good / Warn / Poor given a threshold spec
classify_quality <- function(values, thr) {
  if (is.null(thr) || is.null(values)) return(NULL)
  v <- as.numeric(values)
  out <- rep(NA_character_, length(v))
  if (identical(thr$direction, "higher")) {
    out[v >= thr$good]                  <- "Good"
    out[v <  thr$good & v >= thr$warn]  <- "Warn"
    out[v <  thr$warn]                  <- "Poor"
  } else {
    out[v <= thr$good]                  <- "Good"
    out[v >  thr$good & v <= thr$warn]  <- "Warn"
    out[v >  thr$warn]                  <- "Poor"
  }
  factor(out, levels = names(QUALITY_COLOURS))
}

# ---------------------------
# Circular ("fan") layout: stats as concentric rings instead of side panels
# ---------------------------
build_circular_plot <- function(tree, tips_order, data_quast = NULL, data_genes = NULL,
                                data_busco_geno = NULL, data_busco_prot = NULL,
                                data_nseqs = NULL, data_ortho = NULL, data_fcs = NULL,
                                text_size = 3, skip = NULL, rings = NULL,
                                quality_preset = "generic", thresholds = NULL, show_values = FALSE,
                                open_angle = 14, ring_width = 0.13) {

  if (is.null(skip)) skip <- character(0)

  thr_set <- QUALITY_PRESETS[[quality_preset]]
  if (is.null(thr_set)) {
    warning("Unknown quality_preset '", quality_preset, "', falling back to 'generic'")
    thr_set <- QUALITY_PRESETS[["generic"]]
  }
  # Explicit per-metric thresholds (e.g. from the Shiny app) override the preset
  if (!is.null(thresholds)) {
    for (m in names(thresholds)) {
      if (!is.null(thresholds[[m]])) thr_set[[m]] <- thresholds[[m]]
    }
  }

  # Single source of truth for key typography: headers one step larger than body
  # text. Used by the hand-drawn keys (geom text, mm) and - via .pt - the ggplot
  # legends (theme text, points), so the two cannot drift apart.
  key_title_size <- text_size * 1.15
  key_text_size  <- text_size * 0.95

  n_tips <- length(tips_order)
  labs <- gsub("_", " ", tips_order)              # tree tip labels in node order
  num_df  <- data.frame(label = labs, num = seq_len(n_tips), stringsAsFactors = FALSE)
  ring_df <- data.frame(label = labs, stringsAsFactors = FALSE)

  ring_specs <- list()
  add_ring <- function(spec, values) {
    if (is.null(values) || all(is.na(values))) return(invisible())
    name <- spec$name
    if (identical(spec$scale, "quality")) {
      grade <- classify_quality(values, thr_set[[spec$metric]])
      if (is.null(grade)) return(invisible())
      ring_df[[name]] <<- grade
    } else {
      ring_df[[name]] <<- values
    }
    spec$values <- values             # keep raw values for the printed labels
    ring_specs[[length(ring_specs) + 1]] <<- spec
  }
  # Aligns df[[col]] to the ring's full n_tips positions by node index, NA
  # everywhere df has no row. df is routinely a subset (e.g. FCS-GX only runs
  # for taxid-bearing samples) - returning just the present rows in df's own
  # order would give a short vector that R silently RECYCLES across every
  # tip when it's assigned into ring_df, fabricating data for species that
  # have none at all instead of leaving them blank.
  col_by_node <- function(df, col) {
    out <- rep(NA_real_, n_tips)
    out[df$node] <- as.numeric(df[[col]])
    out
  }

  # Registry of every stat that can be a ring. Each key matches its
  # rectangular-layout panel (and the Shiny "Skip Statistics" checkboxes).
  # `scale` is "quality" (discrete traffic light, needs `metric` for thresholds)
  # or "descriptive" (sequential ramp, no good/bad implied).
  ring_registry <- list(
    ch_plot         = list(name = "Seq number",     scale = "quality", metric = "seq_number",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "Sequences")),
    len_plot        = list(name = "Genome size",    scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "Length")),
    n50_plot        = list(name = "N50",            scale = "quality", metric = "n50",
                           get = function() if (!is.null(data_quast)) col_by_node(data_quast$full, "N50")),
    gene_plot       = list(name = "Gene number",    scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() {
                             if (is.null(data_genes)) return(NULL)
                             tot <- data_genes[data_genes$stat == "Total", ]
                             as.numeric(tot[order(tot$node), ]$value)
                           }),
    busco_gen_plot  = list(name = "BUSCO genome",   scale = "quality", metric = "busco_complete",
                           get = function() if (!is.null(data_busco_geno)) data_busco_geno$Single + data_busco_geno$Duplicated),
    busco_prot_plot = list(name = "BUSCO protein",  scale = "quality", metric = "busco_complete",
                           get = function() if (!is.null(data_busco_prot)) data_busco_prot$Single + data_busco_prot$Duplicated),
    busco_dup_plot  = list(name = "BUSCO duplicated", scale = "quality", metric = "busco_duplicated",
                           get = function() if (!is.null(data_busco_geno)) data_busco_geno$Duplicated),
    fcs_plot        = list(name = "FCS non-contaminant %", scale = "quality", metric = "fcs_noncontam",
                           get = function() if (!is.null(data_fcs)) col_by_node(data_fcs, "non_contaminant_pct")),
    nseqs_plot      = list(name = "Seqs ≥5 BUSCOs", scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_nseqs)) col_by_node(data_nseqs, names(data_nseqs)[2])),
    ortho_plot      = list(name = "Ortho seqs",     scale = "descriptive", low = "#f0f0f0", high = "#4d4d4d",
                           get = function() if (!is.null(data_ortho)) col_by_node(data_ortho, names(data_ortho)[2]))
  )

  # `rings` fixes which stats to show and their inner->outer order (the curated
  # static default). When NULL (e.g. the Shiny app) every available stat is
  # shown. Either way `skip` is still honoured.
  order_keys <- if (!is.null(rings)) rings else names(ring_registry)
  for (key in order_keys) {
    entry <- ring_registry[[key]]
    if (is.null(entry) || key %in% skip) next
    add_ring(entry, entry$get())
  }

  # Keep the traffic-light (quality) rings visually separate from the neutral
  # descriptive rings: descriptive inner (nearest the tree), quality outer (on the
  # rim, where they read most clearly), preserving the requested order within each
  # group.
  is_quality <- vapply(ring_specs, function(s) identical(s$scale, "quality"), logical(1))
  ring_specs <- c(ring_specs[!is_quality], ring_specs[is_quality])

  p <- ggtree(tree, layout = "fan", open.angle = open_angle, size = 0.5, colour = "grey30")

  n_rings <- length(ring_specs)

  ring_offset <- 0.055

  if (n_rings > 0) {
    # PRIMER RING: geom_fruit()'s very first call on a fresh ggtree allocates a
    # one-off, disproportionately large slice of radius regardless of width/
    # offset (a ggtreeExtra quirk, not something width/offset can compensate
    # for) - every ring after the first is sized consistently. A fully
    # transparent zero-content ring absorbs that anomaly so every REAL ring
    # below gets uniform, comparable thickness. Deliberately uses only width/
    # offset (never pwidth: explicitly passing it - at any value - corrupts
    # the whole plot's scale in the ggtreeExtra version this pipeline pins).
    p <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1), fill = NA, colour = NA,
        width = ring_width, offset = 0.10
      ) +
      ggnewscale::new_scale_fill()

    # CALIBRATION: geom_fruit's `width` is not the ring's actual rendered
    # thickness - rings are drawn as overlapping tiles staggered by `offset`,
    # each one visually clipped by the next ring drawn on top of it down to
    # the gap between their start positions. The last ring has nothing after
    # it to clip its trailing edge, so it alone renders at its full declared
    # `width` - several times wider than every other (clipped) ring. Measure
    # the true stagger `ring_offset` produces on this tree with a throwaway
    # invisible ring, then use that measured value (not ring_width) as every
    # real ring's width below, so every ring - including the last - renders
    # the same actual thickness.
    range_before_calib <- suppressWarnings(max(vapply(
      ggplot_build(p)$data,
      function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
      numeric(1)), na.rm = TRUE))
    p_calib <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1), fill = NA, colour = NA,
        width = ring_width, offset = ring_offset
      )
    range_after_calib <- suppressWarnings(max(vapply(
      ggplot_build(p_calib)$data,
      function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
      numeric(1)), na.rm = TRUE))
    measured_width <- range_after_calib - range_before_calib
  }

  # Quality rings share one discrete Good/Warn/Poor scale (legend shown once -
  # the ring key on the left says which ring is which); descriptive rings keep a
  # sequential ramp, each with its own legend ordered outer-ring-first.
  shown_quality_legend <- FALSE
  for (i in seq_along(ring_specs)) {
    spec <- ring_specs[[i]]
    p <- p +
      ggtreeExtra::geom_fruit(
        data = ring_df, geom = geom_tile,
        mapping = aes(y = label, x = 1, fill = .data[[spec$name]]),
        width = measured_width, offset = ring_offset,
        color = "white", linewidth = 0.2
      )
    if (identical(spec$scale, "quality")) {
      # The Good/Warn/Poor key is drawn manually in the side panel: ggplot only
      # renders legend keys for levels present in the layer that owns the legend,
      # so Warn/Poor vanish whenever the first quality ring happens to be all-Good.
      p <- p + scale_fill_manual(
        values = QUALITY_COLOURS, limits = names(QUALITY_COLOURS),
        drop = FALSE, na.value = "grey90", name = "Quality", guide = "none"
      )
      shown_quality_legend <- TRUE
    } else {
      p <- p + scale_fill_gradient(
        low = spec$low, high = spec$high, name = spec$name,
        guide = guide_colourbar(order = n_rings - i + 1)
      )
    }
    p <- p + ggnewscale::new_scale_fill()
  }

  # Outermost labels: tip numbers, kept upright. geom_fruit/geom_tiplab rotate
  # text tangentially in a fan layout, so instead we place a plain geom_text
  # (angle = 0 keeps the glyphs horizontal) just beyond the outer ring radius.
  ring_max_x <- suppressWarnings(max(vapply(
    ggplot_build(p)$data,
    function(dd) if ("x" %in% names(dd)) max(dd$x, na.rm = TRUE) else NA_real_,
    numeric(1)), na.rm = TRUE))
  tip_pos <- p$data[p$data$isTip, ]
  num_pos <- data.frame(y = tip_pos$y, num = match(tip_pos$label, labs))

  # All circular text (tip numbers, legends, species key) scales with text_size.
  p <- p +
    geom_text(data = num_pos, aes(x = ring_max_x * 1.06, y = y, label = num),
              angle = 0, size = text_size * 0.9, inherit.aes = FALSE) +
    theme(legend.position = "right",
          legend.title = element_text(size = key_title_size * .pt, face = "bold"),
          legend.text = element_text(size = key_text_size * .pt),
          legend.key.width = unit(0.3, "cm"),
          legend.key.height = unit(0.35, "cm"))

  # Optional redundant encoding: print each value on its ring, so the figure is
  # not colour-only (important for the traffic-light rings). Kept upright.
  if (isTRUE(show_values) && n_rings > 0) {
    fmt_val <- function(v) {
      v <- as.numeric(v)
      m <- suppressWarnings(max(abs(v), na.rm = TRUE))
      if (!is.finite(m))                                    rep("", length(v))
      else if (m >= 1e6)                                    sprintf("%.1fM", v / 1e6)
      else if (m >= 1000)                                   format(round(v), big.mark = ",", trim = TRUE)
      else if (all(abs(v - round(v)) < 1e-8, na.rm = TRUE)) as.character(round(v))
      else                                                  sprintf("%.1f", v)
    }
    tile_layers <- Filter(function(dd) all(c("xmin", "xmax") %in% names(dd)),
                          ggplot_build(p)$data)
    # tile_layers[[1]] is the invisible primer ring (always added above when
    # n_rings > 0), not a real ring - skip it so values land on the ring they
    # actually describe instead of the one drawn just inside it.
    if (length(tile_layers) >= n_rings + 1) {
      for (i in seq_len(n_rings)) {
        dd <- tile_layers[[i + 1]]
        # NB: the radius must live in the data, not the aes expression - aes() is
        # evaluated lazily, so aes(x = rad) would resolve every layer to the last
        # value of the loop variable.
        val_pos <- data.frame(
          x   = mean(c(dd$xmin, dd$xmax), na.rm = TRUE),
          y   = tip_pos$y,
          lab = fmt_val(ring_specs[[i]]$values)[match(tip_pos$label, labs)]
        )
        p <- p + geom_text(data = val_pos, aes(x = x, y = y, label = lab),
                           angle = 0, size = text_size * 0.6, colour = "grey15",
                           inherit.aes = FALSE)
      }
    }
  }

  # Species key (number -> italic name): a compact, top-aligned list on the left
  sp_txt <- paste0(num_df$num, "  ", num_df$label)
  sp_block <- paste(sp_txt, collapse = "\n")

  # Ring key (inner -> outer). Needed because all quality rings share the same
  # Good/Warn/Poor colours, so colour alone cannot identify a ring.
  ring_names <- vapply(ring_specs, function(s) s$name, character(1))
  ring_block <- paste0(seq_along(ring_names), ". ", ring_names, collapse = "\n")

  has_quality <- any(vapply(ring_specs, function(s) identical(s$scale, "quality"), logical(1)))
  lh <- 0.035                                   # one text line, in panel units
  y  <- 1.00
  sp_leg <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
    annotate("text", x = 0, y = y, label = "Species",
             hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
  y <- y - lh * 1.6
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = sp_block,
             hjust = 0, vjust = 1, size = key_text_size, fontface = "italic", lineheight = 1.2)
  y <- y - lh * n_tips - lh * 0.8
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = "Rings (inner -> outer)",
             hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
  y <- y - lh * 1.6
  sp_leg <- sp_leg +
    annotate("text", x = 0, y = y, label = ring_block,
             hjust = 0, vjust = 1, size = key_text_size, lineheight = 1.2)
  y <- y - lh * length(ring_names) - lh * 0.8

  # Manual Good / Warn / Poor key - always shows all three swatches
  if (has_quality) {
    sp_leg <- sp_leg +
      annotate("text", x = 0, y = y, label = "Quality",
               hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
    y <- y - lh * 1.4
    for (k in seq_along(QUALITY_COLOURS)) {
      yy <- y - (k - 1) * lh
      sp_leg <- sp_leg +
        annotate("rect", xmin = 0, xmax = 0.05,
                 ymin = yy - lh * 0.62, ymax = yy - lh * 0.08,
                 fill = QUALITY_COLOURS[[k]], colour = NA) +
        annotate("text", x = 0.075, y = yy - lh * 0.35,
                 label = names(QUALITY_COLOURS)[k],
                 hjust = 0, vjust = 0.5, size = key_text_size)
    }
    y <- y - lh * length(QUALITY_COLOURS) - lh * 0.8

    # Spell out what earns each grade, per quality ring, so the figure is
    # self-documenting rather than relying on the reader knowing the preset.
    fmt_thr <- function(v) {
      if (!is.finite(v)) return("?")
      if (abs(v) >= 1e6)       sprintf("%.1fM", v / 1e6)
      else if (abs(v) >= 1000) format(round(v), big.mark = ",", trim = TRUE)
      else                     as.character(round(v, 1))
    }
    rules <- character(0)
    for (i in seq_along(ring_specs)) {
      sp  <- ring_specs[[i]]
      if (!identical(sp$scale, "quality")) next
      thr <- thr_set[[sp$metric]]
      if (is.null(thr)) next
      rule <- if (identical(thr$direction, "higher")) {
        sprintf("%s+ / %s+ / <%s", fmt_thr(thr$good), fmt_thr(thr$warn), fmt_thr(thr$warn))
      } else {
        sprintf("<=%s / <=%s / >%s", fmt_thr(thr$good), fmt_thr(thr$warn), fmt_thr(thr$warn))
      }
      rules <- c(rules, paste0(i, ". ", sp$name, ": ", rule))
    }
    if (length(rules) > 0) {
      sp_leg <- sp_leg +
        annotate("text", x = 0, y = y, label = "Thresholds (Good / Warn / Poor)",
                 hjust = 0, vjust = 1, fontface = "bold", size = key_title_size)
      y <- y - lh * 1.4
      sp_leg <- sp_leg +
        annotate("text", x = 0, y = y, label = paste(rules, collapse = "\n"),
                 hjust = 0, vjust = 1, size = key_text_size, lineheight = 1.25)
    }
  }

  sp_leg + p + plot_layout(widths = c(0.35, 1))
}

# ---------------------------
# Final wrapper (same name)
# ---------------------------
generate_complete_plot <- function(processed_data, text_size = 3, tree_scale = 0.0005,
                                   bar_width = 0.7, rad_width = 0.4, skip_stats = NULL,
                                   top_margin = 5.5, right_margin = 5.5, bottom_margin = 5.5,
                                   left_margin = 5.5, tree_margin = 15, tree_space_ratio = 1.3,
                                   tree_style = "roundrect", circular_rings = NULL,
                                   ring_width = 0.13, quality_preset = "generic",
                                   thresholds = NULL, show_values = FALSE) {

  # Circular layout is a separate plotting path (rings instead of side panels).
  # circular_rings = NULL shows every available stat (the interactive Shiny
  # default); pass an ordered vector of keys to curate/fix the rings.
  # Note: only text_size, skip_stats, circular_rings and ring_width affect this
  # layout; the bar/pie/margin/tree-scale controls apply to the other styles.
  if (tree_style == "circular") {
    return(build_circular_plot(
      tree            = processed_data$tree,
      tips_order      = processed_data$tips_order,
      data_quast      = processed_data$data_quast,
      data_genes      = processed_data$data_genes,
      data_busco_geno = processed_data$data_busco_geno,
      data_busco_prot = processed_data$data_busco_prot,
      data_nseqs      = processed_data$data_nseqs,
      data_ortho      = processed_data$data_ortho,
      data_fcs        = processed_data$data_fcs,
      text_size       = text_size,
      skip            = skip_stats,
      rings           = circular_rings,
      ring_width      = ring_width,
      quality_preset  = quality_preset,
      thresholds      = thresholds,
      show_values     = show_values
    ))
  }

  plot_results <- generate_plots(
    processed_data = processed_data,
    text_size = text_size,
    tree_scale = tree_scale,
    bar_width = bar_width,
    rad_width = rad_width,
    skip_stats = skip_stats,
    tree_style = tree_style
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
