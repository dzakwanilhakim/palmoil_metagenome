# =============================================================================
# R/functions_legends.R — Results/legends/: shared keys for every figure
#
# Built from docs/supplementary.md (treatment codes, biofertilizer contents)
# plus the pipeline's own conventions (timepoint colours, plot types):
#   treatment_TM.png, treatment_Nursery.png  plot label -> treatment
#   biofertilizer_contents.png               product -> microbes
#   timepoint_colours.png                    T0..T4 colour key
#   figure_legends.md                        caption per image type
# Plots label treatments PH1..PH6 (TM) and PHN01..PHN13 (Nursery); the
# supplementary tables use PH01 / PHN1, so both codes are shown.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# Every markdown pipe table in the file, as tibbles (all columns character).
.read_md_tables <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  is_row <- grepl("^\\s*\\|", lines)
  runs <- rle(is_row)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  blocks <- Map(function(s, e) lines[s:e], starts[runs$values], ends[runs$values])
  lapply(blocks, function(b) {
    cells <- lapply(b, function(x)
      trimws(strsplit(sub("^\\s*\\|", "", sub("\\|\\s*$", "", x)), "|", fixed = TRUE)[[1]]))
    header <- cells[[1]]
    body <- cells[-(1:2)]   # drop header and |---| separator
    tibble::as_tibble(stats::setNames(
      as.data.frame(do.call(rbind, body), stringsAsFactors = FALSE), header))
  })
}

# PH01 -> PH1 (TM plots); PHN1 -> PHN01 (Nursery plots)
.plot_label <- function(code) {
  n <- as.integer(stringr::str_extract(code, "[0-9]+$"))
  ifelse(startsWith(code, "PHN"), sprintf("PHN%02d", n), paste0("PH", n))
}

# Draw a table as a PNG; cells may use markdown *italic* (ggtext).
.legend_table_png <- function(df, title, out_path, wrap, style) {
  wrap <- rep_len(wrap, ncol(df))
  wrapped <- purrr::map2(df, wrap, function(col, w)
    gsub("\n", "<br>", stringr::str_wrap(col, w)))
  n_lines <- do.call(pmax, lapply(wrapped, function(x) stringr::str_count(x, "<br>") + 1))
  row_h <- n_lines + 0.6
  row_top <- c(0, cumsum(row_h))[seq_len(nrow(df))]
  x_left <- c(0, cumsum(wrap + 2))[seq_len(ncol(df))]
  total_w <- sum(wrap + 2)
  header_h <- 1.6

  cells <- purrr::imap_dfr(wrapped, function(col, nm) tibble::tibble(
    column = nm, row = seq_along(col), label = col))
  cells$x <- x_left[match(cells$column, names(df))] + 0.5
  cells$y <- -(header_h + row_top[cells$row] + row_h[cells$row] / 2)
  stripes <- tibble::tibble(row = seq_len(nrow(df)),
                            ymax = -(header_h + row_top), ymin = -(header_h + row_top + row_h)) |>
    dplyr::filter(row %% 2 == 0)
  header <- tibble::tibble(label = paste0("**", names(df), "**"),
                           x = x_left + 0.5, y = -header_h / 2)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = stripes, ggplot2::aes(xmin = 0, xmax = total_w,
                                                    ymin = ymin, ymax = ymax),
                       fill = "grey94") +
    ggplot2::annotate("rect", xmin = 0, xmax = total_w, ymin = -header_h, ymax = 0,
                      fill = "grey85") +
    ggtext::geom_richtext(data = header, ggplot2::aes(x, y, label = label),
                          hjust = 0, size = 4.2, fill = NA, label.colour = NA) +
    ggtext::geom_richtext(data = cells, ggplot2::aes(x, y, label = label),
                          hjust = 0, size = 3.8, fill = NA, label.colour = NA,
                          lineheight = 1.1) +
    ggplot2::scale_x_continuous(limits = c(0, total_w), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-(header_h + sum(row_h)), 0),
                                expand = ggplot2::expansion(add = c(0.4, 0))) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = style$base_size %||% 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.margin = ggplot2::margin(10, 14, 10, 14))
  ggplot2::ggsave(out_path, p, width = 0.088 * total_w + 0.5,
                  height = 0.19 * (header_h + sum(row_h)) + 0.9, dpi = 200)
  out_path
}

.timepoint_key_png <- function(out_path, style) {
  pal <- waktu_palette()
  d <- tibble::tibble(waktu = factor(names(pal), levels = names(pal)), x = seq_along(pal))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, 1, fill = waktu)) +
    ggplot2::geom_point(shape = 21, size = 9, colour = "grey20") +
    ggplot2::geom_text(ggplot2::aes(y = 0.55, label = waktu), size = 5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0.3, 1.3)) +
    ggplot2::scale_x_continuous(limits = c(0.5, length(pal) + 0.5)) +
    ggplot2::labs(title = "Timepoint colours (all figures)",
                  subtitle = "T0 = baseline before treatment; T1..T4 = later sampling rounds") +
    ggplot2::theme_void(base_size = style$base_size %||% 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.margin = ggplot2::margin(10, 14, 10, 14))
  ggplot2::ggsave(out_path, p, width = 7, height = 2.4, dpi = 200)
  out_path
}

# One caption per image type. {cohort} = 16S_TM / ITS_TM / 16S_Nursery,
# {kebun} = Kode Kebun, {rank} = phylum / genus / species, {Tx}/{Ty} = timepoints.
FIGURE_LEGENDS <- list(
  "Shared keys" = c(
    "Treatment codes on x axes and in titles: see `treatment_TM.png` (PH1-PH6) and `treatment_Nursery.png` (PHN01-PHN13); biofertilizer products PH01-PH05: `biofertilizer_contents.png`. `Untreated` = T0 sample without a treatment code (pre-treatment baseline).",
    "Timepoint colours are the same in every figure: see `timepoint_colours.png`.",
    "Cohorts are analysed separately and never compared: 16S_TM, ITS_TM, 16S_Nursery. 16S = bacteria/archaea, ITS = fungi; TM = mature plantation, Nursery = nursery.",
    "Fertilizers are pooled within a kebun for the temporal analyses (ANCOM-BC2, FAPROTAX): those figures compare timepoints, not treatments.",
    "Sequencing/extraction batch overlaps with timepoint, and 16S batches were classified with two different reference databases; read temporal differences with that caveat."),
  "Results/dashboard" = c(
    "`qc_status_pie.png` - Share of sequenced samples per QC class (HIGH PASS = mapped reads at or above the cutoff, 16S 50,000 / ITS 10,000; LOW PASS = below cutoff; NO DATA = no run).",
    "`prepost_qc_depth.png` - Read depth per sample before and after the taxonomic and depth filters.",
    "`rarefaction_{cohort}.png` - Rarefaction curves: observed taxa as reads are subsampled; a flattening curve means sequencing depth captured most of the diversity.",
    "`pcoa_{cohort}.png` - Principal coordinates of between-sample distances; points coloured by timepoint.",
    "`shannon_barplot.png` - Shannon diversity per sample."),
  "Results/data_counts/{cohort}" = c(
    "`{kebun}_counts.png` - Number of analysed biological replicates per treatment (rows) and timepoint (columns) in one kebun."),
  "Goal A - one kebun, one timepoint (`Goal_A_Intra_Snapshot/{kebun}`)" = c(
    "`alpha_box_*.png` - Alpha diversity (Observed richness, Pielou evenness, Shannon) per treatment at one timepoint; boxes = median and interquartile range, points = samples. Kruskal-Wallis / Wilcoxon where a group has n >= 2; markers use raw p, CSVs carry raw and BH-adjusted p."),
  "Goal B - one kebun over time (`Goal_B_Intra_Longitudinal/{kebun}`)" = c(
    "`alpha_traj_*.png` - Alpha diversity per treatment across timepoints.",
    "`beta_ord_*.png` - Ordination of beta diversity (Aitchison, Bray-Curtis or Jaccard as named in the file); colours = timepoint; PERMANOVA in the CSVs.",
    "`beta_dendro_*.png` - Ward.D2 hierarchical clustering of samples on the same distance.",
    "`stackbar_{rank}_topN_{kebun}.png` - Mean relative abundance per treatment (x) and timepoint (panels) from the rarefied species table. The top 50 taxa are selected independently at each rank (phylum, genus, species), ranked by summed group mean relative abundance; the remaining taxa are summed as Other (light grey); Unknown (dark grey) = unresolved name."),
  "Goal C - all kebun, one timepoint (`Goal_C_Cross_Snapshot`)" = c(
    "`alpha_box_crossfield_*.png` - Alpha diversity per treatment at one timepoint with all kebun pooled."),
  "Goal D - all kebun over time (`Goal_D_Cross_Longitudinal`)" = c(
    "`alpha_traj_pooled_*.png`, `alpha_traj_faceted_*.png` - Alpha diversity over time, pooled across kebun or one panel per kebun.",
    "`beta_ord_*_faceted.png`, `beta_dendro_*_pooled.png` - Beta diversity ordination and clustering across kebun.",
    "`stackbar_{rank}_topN_pooled.png` - As Goal B stacked bars, all kebun pooled."),
  "ANCOM-BC2 differential abundance (`Differential_Abundance_Temporal/{kebun}`)" = c(
    "`{Tx}_vs_{Ty}/{rank}/dumbbell.png` - Taxa whose abundance differs between two timepoints (BH q < 0.05, passing the pseudo-count sensitivity check); grey dot = 0, coloured dot = bias-corrected log fold change (positive = higher at the later timepoint). Pairs are consecutive (T0 vs T1, T1 vs T2, ...) plus first vs last timepoint (e.g. T0 vs T3).",
    "`Trajectory_{T0}_to_{Tn}/{rank}/trajectory.png` - Log fold change of each taxon versus the first timepoint (0 = no change) for the top 20 taxa, ranked by the ANCOM-BC2 global test across all timepoints. Lines = taxa (labelled at the right); point fill = timepoint; large point = significant versus the first timepoint (BH q < 0.05)."),
  "FAPROTAX functional profiles (`Functional_FAPROTAX`, 16S only)" = c(
    "Functions are predicted from taxon names (FAPROTAX 1.2.12) and are putative. Clinical/host-associated groups (human pathogens, human/mammal gut, animal parasites) are excluded from every Top-N plot.",
    "`{kebun}/{Tx}_vs_{Ty}/dumbbell.png` - Mean functional score at two timepoints for functions that differ (Wilcoxon, BH q < 0.05), top 20 by absolute change; dots coloured by timepoint. All dumbbells of one kebun share one x axis (square-root scale).",
    "`faprotax_effect_heatmap.png` - Change in mean score (later minus earlier, percentage points) per function and comparison; red = increase, blue = decrease; * = BH q < 0.05; x = insufficient samples.",
    "`faprotax_trajectories.png` - Selected functions over time: grey points = samples, coloured points = means with standard error.",
    "`{kebun}/faprotax_trajectory_panel.png` - One panel, first to last timepoint: log2 fold change of the mean score versus the first timepoint for the top 20 functions (Kruskal-Wallis across time, BH). Lines = functions; point fill = timepoint; large point = significant versus the first timepoint (Wilcoxon, BH q < 0.05).",
    "`faprotax_annotation_coverage.png` - Fraction of each sample's reads assigned to at least one FAPROTAX function.")
)

.write_figure_legends_md <- function(out_path) {
  body <- purrr::imap_chr(FIGURE_LEGENDS, function(items, heading)
    paste0("## ", heading, "\n\n", paste0("- ", items, collapse = "\n"), "\n"))
  writeLines(c("# Figure legends",
               "",
               "Captions for every image type in `Results/`. Placeholders: `{cohort}`, `{kebun}`, `{rank}`, `{Tx}`/`{Ty}` (timepoints), `topN` (number of named taxa).",
               "",
               body), out_path, useBytes = TRUE)
  out_path
}

build_legend_folder <- function(supplementary_path = "docs/supplementary.md",
                                style = load_plot_style(),
                                out_dir = "Results/legends") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- .read_md_tables(supplementary_path)
  is_product <- vapply(tables, function(t) "Produk" %in% names(t), logical(1))
  codes <- lapply(tables, function(t) t[[1]])
  is_nursery <- !is_product & vapply(codes, function(x) all(startsWith(x, "PHN")), logical(1))
  is_tm <- !is_product & !is_nursery

  treatment <- function(t, stage) {
    t |>
      dplyr::transmute(`Plot label` = .plot_label(.data[[names(t)[1]]]),
                       Kode = .data[[names(t)[1]]],
                       Perlakuan = .data[[names(t)[2]]])
  }
  written <- character()
  if (any(is_tm))
    written <- c(written, .legend_table_png(
      treatment(tables[[which(is_tm)[1]]]), "Treatments - TM (mature plantation)",
      file.path(out_dir, "treatment_TM.png"), c(11, 7, 60), style))
  if (any(is_nursery))
    written <- c(written, .legend_table_png(
      treatment(tables[[which(is_nursery)[1]]]), "Treatments - Nursery",
      file.path(out_dir, "treatment_Nursery.png"), c(11, 7, 40), style))
  if (any(is_product)) {
    prod <- tables[[which(is_product)[1]]]
    written <- c(written, .legend_table_png(
      dplyr::select(prod, -dplyr::any_of("No.")), "Commercial biofertilizer products",
      file.path(out_dir, "biofertilizer_contents.png"), c(7, 22, 70), style))
  }
  written <- c(written,
               .timepoint_key_png(file.path(out_dir, "timepoint_colours.png"), style),
               .write_figure_legends_md(file.path(out_dir, "figure_legends.md")))
  message("Legends written to ", out_dir, ": ", paste(basename(written), collapse = ", "))
  written
}
