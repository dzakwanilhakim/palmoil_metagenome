# =============================================================================
# R/functions_qc.R  —  Stage 3, Steps 1-2 (taxonomic filter + pre-QC assessment)
# Per-marker. Steps 3-5 (threshold application, sync, rarefy+CLR) are wired
# AFTER the user inspects the pre-QC visuals and supplies thresholds.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ----------------------------------------------------------------------------
# STEP 1. TAXONOMIC FILTERING  (per marker)
#   Input: long_qc (genus rows carry the `tax` lineage string:
#          "Superkingdom;Kingdom;Phylum;Class;Order;Family;Genus")
#   16S : keep superkingdom in {Bacteria, Archaea}
#   ITS : keep kingdom == Fungi
#   both: drop order == Chloroplast OR family == Mitochondria
#   Off-target / Unassigned / blank are dropped by the keep-rules above.
# ----------------------------------------------------------------------------

# split the `tax` lineage into ranks. ONT/EPI2ME lineage is 7 levels:
# superkingdom; kingdom; phylum; class; order; family; genus
.split_lineage <- function(tax) {
  parts <- stringr::str_split(tax, ";", simplify = TRUE)
  # pad to 7 columns if some rows are shorter
  if (ncol(parts) < 7) {
    parts <- cbind(parts,
                   matrix("", nrow = nrow(parts), ncol = 7 - ncol(parts)))
  }
  colnames(parts) <- c("superkingdom","kingdom","phylum_l","class",
                       "order","family","genus_l")[seq_len(ncol(parts))]
  tibble::as_tibble(parts[, 1:7, drop = FALSE])
}

taxonomic_filter <- function(long_qc) {
  # operate on genus-rank rows (analysis rank); phylum rows handled separately
  dat <- dplyr::filter(long_qc, rank == "genus")

  lin <- .split_lineage(dat$tax)
  dat <- dplyr::bind_cols(dat, lin)

  is_blank <- function(x) is.na(x) | stringr::str_trim(x) == "" |
                          stringr::str_detect(x, "(?i)^unassigned|^unclassified|^unknown$")

  keep_16s <- with(dat, marker == "16S" &
                        superkingdom %in% c("Bacteria", "Archaea"))
  keep_its <- with(dat, marker == "ITS" & kingdom == "Fungi")

  # organelle removal (both markers)
  organelle <- with(dat, stringr::str_detect(order,  "(?i)chloroplast") |
                          stringr::str_detect(family, "(?i)mitochondria"))
  organelle[is.na(organelle)] <- FALSE

  keep <- (keep_16s | keep_its) & !organelle

  # diagnostics
  dropped <- dat[!keep, ]
  msg <- dat |>
    dplyr::mutate(.kept = keep) |>
    dplyr::group_by(marker) |>
    dplyr::summarise(otu_rows = dplyr::n(),
                     kept = sum(.kept),
                     dropped = sum(!.kept), .groups = "drop")
  message("Taxonomic filter (genus rows):")
  print(msg)

  dplyr::filter(dat, keep) |>
    dplyr::select(-dplyr::any_of(c("phylum_l","genus_l","class")))
}

# ----------------------------------------------------------------------------
# STEP 2. PRE-QC DEPTH ASSESSMENT  (per marker)  — then STOP for thresholds
# ----------------------------------------------------------------------------

# per-sample sequencing depth (sum of counts over kept taxa), per marker
sample_depth <- function(filtered_long) {
  filtered_long |>
    dplyr::group_by(marker, id_sampel, stage) |>
    dplyr::summarise(depth = sum(count, na.rm = TRUE), .groups = "drop")
}

# 5-number summary per marker
depth_summary <- function(depth_tbl) {
  depth_tbl |>
    dplyr::group_by(marker) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      Min = min(depth), Q1 = quantile(depth, .25),
      Median = median(depth), Q3 = quantile(depth, .75),
      Max = max(depth), .groups = "drop")
}

# boxplot of sample depths, faceted by marker, colored by stage
plot_depth_boxplot <- function(depth_tbl,
                               out_path = "results/preqc_depth_boxplot.png",
                               thresholds = NULL) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  # drop zero-depth samples from the log plot but report them
  n_zero <- sum(depth_tbl$depth <= 0)
  pdat <- dplyr::filter(depth_tbl, depth > 0)

  # n per marker for subtitle
  n_lab <- pdat |> dplyr::count(marker) |>
    dplyr::mutate(lab = paste0(marker, " (n=", n, ")"))

  pal <- c(Nursery = "#E07A5F", TM = "#3D8C8C")

  p <- ggplot2::ggplot(pdat,
                       ggplot2::aes(x = marker, y = depth)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = stage),
                          width = .55, outlier.shape = NA, alpha = .85,
                          position = ggplot2::position_dodge(.7),
                          linewidth = .4, colour = "grey25") +
    ggplot2::geom_point(ggplot2::aes(fill = stage),
                        position = ggplot2::position_jitterdodge(
                          jitter.width = .12, dodge.width = .7),
                        shape = 21, size = 1.6, stroke = .2,
                        colour = "grey20", alpha = .55) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::annotation_logticks(sides = "l", colour = "grey60") +
    ggplot2::scale_fill_manual(values = pal, name = "Stage") +
    ggplot2::labs(
      title = "Pre-QC sample sequencing depth",
      subtitle = "Read depth per sample after taxonomic filtering, before threshold cuts",
      x = NULL, y = "Total reads (log scale)") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 10),
      legend.position = "right",
      axis.text = ggplot2::element_text(colour = "grey20"))

  # optional threshold reference lines (per marker) if provided
  if (!is.null(thresholds)) {
    th <- tibble::tibble(
      marker = names(thresholds),
      min_depth = unlist(lapply(thresholds, function(x) x$min_sample_depth)))
    th <- dplyr::filter(th, !is.na(min_depth))
    if (nrow(th))
      p <- p + ggplot2::geom_hline(data = th,
                ggplot2::aes(yintercept = min_depth),
                linetype = "dashed", colour = "firebrick", linewidth = .5)
  }

  if (n_zero > 0)
    p <- p + ggplot2::labs(caption = paste0(n_zero,
              " zero-depth sample(s) omitted from log plot (will be dropped)."))

  ggplot2::ggsave(out_path, p, width = 7.5, height = 5.5, dpi = 300)
  message("Wrote ", out_path, if (n_zero) paste0(" (", n_zero, " zero-depth omitted)"))
  out_path
}

# ----------------------------------------------------------------------------
# INTERACTIVE THRESHOLD EXPLORER  (run in console, NOT inside tar_make)
#   How many samples survive a candidate min-depth, per marker & stage?
#   Usage:
#     depth <- tar_read(preqc)$depth_table
#     explore_threshold(depth, marker = "16S", min_depth = 30000)
#     survival_table(depth, marker = "ITS", candidates = c(5e3,1e4,2e4,5e4))
# ----------------------------------------------------------------------------
explore_threshold <- function(depth_tbl, marker, min_depth) {
  d <- dplyr::filter(depth_tbl, marker == !!marker)
  res <- d |>
    dplyr::mutate(survives = depth >= min_depth) |>
    dplyr::group_by(stage) |>
    dplyr::summarise(total = dplyr::n(),
                     kept  = sum(survives),
                     dropped = sum(!survives),
                     min_kept = ifelse(any(survives), min(depth[survives]), NA),
                     .groups = "drop")
  cat(sprintf("\n%s  @ min_depth = %s\n", marker, format(min_depth, big.mark=",")))
  print(res)
  cat(sprintf("TOTAL kept: %d of %d  | rarefy depth would be: %s\n",
              sum(res$kept), sum(res$total),
              format(min(d$depth[d$depth >= min_depth]), big.mark=",")))
  invisible(res)
}

# sweep several candidate thresholds at once -> survival table
survival_table <- function(depth_tbl, marker, candidates) {
  d <- dplyr::filter(depth_tbl, marker == !!marker)
  purrr::map_dfr(candidates, function(th) {
    tibble::tibble(
      min_depth = th,
      kept = sum(d$depth >= th),
      dropped = sum(d$depth < th),
      pct_kept = round(100 * mean(d$depth >= th), 1),
      rarefy_depth = ifelse(any(d$depth >= th),
                            min(d$depth[d$depth >= th]), NA))
  }) |> dplyr::mutate(marker = marker, .before = 1)
}

# rarefaction curves (vegan::rarecurve) per marker, BEFORE any dropping.
# Builds a sample x taxon count matrix per marker, colored by stage.
plot_rarefaction <- function(filtered_long, marker,
                             step = 200, max_depth = NULL,
                             out_dir = "results") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir, paste0("preqc_rarefaction_", tolower(marker), ".png"))

  dat <- dplyr::filter(filtered_long, marker == !!marker)
  if (nrow(dat) == 0) { message("No data for ", marker); return(NA_character_) }

  wide <- dat |>
    dplyr::group_by(id_sampel, taxon) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    tidyr::pivot_wider(names_from = taxon, values_from = count, values_fill = 0)

  samp_ids <- wide$id_sampel
  m <- as.matrix(wide[, -1]); rownames(m) <- samp_ids
  storage.mode(m) <- "integer"

  stage_of <- dat |> dplyr::distinct(id_sampel, stage) |>
    dplyr::arrange(match(id_sampel, samp_ids))
  cols <- as.integer(factor(stage_of$stage))

  if (is.null(max_depth)) max_depth <- max(rowSums(m))

  grDevices::png(out_path, width = 1400, height = 900, res = 150)
  vegan::rarecurve(m, step = step, label = FALSE, col = cols,
                   xlab = "Sequencing depth", ylab = "Observed genera",
                   main = paste0("Pre-QC rarefaction — ", marker),
                   xlim = c(0, max_depth))
  graphics::legend("bottomright", legend = levels(factor(stage_of$stage)),
                   col = seq_along(levels(factor(stage_of$stage))), lty = 1, bty = "n")
  grDevices::dev.off()
  message("Wrote ", out_path)
  out_path
}

# convenience: run the whole pre-QC bundle, return paths + summary
preqc_assess <- function(filtered_long) {
  depth <- sample_depth(filtered_long)
  summ  <- depth_summary(depth)
  message("\n=== PRE-QC DEPTH SUMMARY (per marker) ===")
  print(summ)
  list(
    depth_table = depth,
    summary     = summ,
    boxplot     = plot_depth_boxplot(depth),
    rarefaction = purrr::map_chr(unique(filtered_long$marker),
                                 ~ plot_rarefaction(filtered_long, .x))
  )
}
