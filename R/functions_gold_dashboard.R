# =============================================================================
# R/functions_gold_dashboard.R  —  Dashboard panels (gold-tier pipeline)
#   Universe = marker x stage (Jenis Kebun): 16S_TM, 16S_Nursery, ITS_TM,
#   ITS_Nursery. Built incrementally:
#     - qc_status pie chart (HIGH PASS / LOW PASS / NO DATA)
#     - Pre-QC vs Post-QC sample depth boxplot
#   Rarefaction curve + PCoA panels follow once these are verified.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

DASHBOARD_DIR <- "Results/dashboard"

.universe_qc_counts <- function(gold_qc, marker_label) {
  # NOTE: gold_qc already has its own "marker" column inherited from
  # silver_wf_runs (NA for metadata-only/unmatched rows) -- assigning to
  # `marker` inside dplyr::mutate() would resolve the bare symbol on the RHS
  # against THAT existing column, not this function's argument. Aggregate
  # first, then attach marker_label afterward to sidestep the collision.
  gold_qc |>
    dplyr::mutate(qc_status = dplyr::coalesce(qc_status, "NO DATA")) |>
    dplyr::filter(!is.na(`Jenis Kebun`)) |>
    dplyr::count(stage = `Jenis Kebun`, qc_status) |>
    dplyr::mutate(marker = marker_label, .before = 1)
}

# ----------------------------------------------------------------------------
# Pie chart: HIGH PASS / LOW PASS / NO DATA counts, faceted by universe
# ----------------------------------------------------------------------------
build_qc_status_pie <- function(gold_qc_16s, gold_qc_its, style = load_plot_style(),
                                out_path = file.path(DASHBOARD_DIR, "qc_status_pie.png")) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  # NO DATA excluded entirely (not just hidden) -- percentages are
  # recomputed over HIGH PASS + LOW PASS only, so the remaining slices
  # still fill the circle to 100%.
  d <- dplyr::bind_rows(.universe_qc_counts(gold_qc_16s, "16S"),
                        .universe_qc_counts(gold_qc_its, "ITS")) |>
    dplyr::filter(qc_status != "NO DATA") |>
    dplyr::mutate(universe = paste0(marker, "_", stage)) |>
    dplyr::group_by(universe) |>
    dplyr::mutate(pct = 100 * n / sum(n),
                  label = paste0(n, "\n(", sprintf("%.1f", pct), "%)")) |>
    dplyr::ungroup()

  pal <- c("HIGH PASS" = "#2E8B57", "LOW PASS" = "#E8A33D")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = "", y = n, fill = qc_status)) +
    ggplot2::geom_col(width = 1, colour = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3.4,
                       fontface = "bold", lineheight = .9) +
    ggplot2::scale_fill_manual(values = pal, name = "QC status") +
    # scales="free": each pie fills its OWN 360 degrees. Default "fixed"
    # trains one shared y-scale off the largest per-universe total, so a
    # universe with a smaller total (fewer samples) only fills a fraction
    # of the circle under coord_polar, leaving a blank gap -- looks exactly
    # like a missing slice.
    ggplot2::facet_wrap(~ universe, scales = "free") +
    ggplot2::labs(title = "Sample QC status by universe",
                  subtitle = "HIGH PASS / LOW PASS (NO DATA excluded), gold-tier read-mapping gate") +
    gold_plot_theme(style) +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  ggplot2::ggsave(out_path, p, width = 9, height = 7, dpi = 200)
  message("Wrote ", out_path)
  out_path
}

# ----------------------------------------------------------------------------
# Pre-QC vs Post-QC sample depth boxplot, per universe
#   Pre-QC  = Mapped reads from gold_qc (all samples, before taxonomic/OTU
#             filtering, i.e. right after the read-mapping gate)
#   Post-QC = surviving per-sample depth after taxonomic + threshold filter
#             (gold_processed_matrix, species rank)
# ----------------------------------------------------------------------------
.processed_matrix_depth <- function(mat, marker) {
  taxon_col <- names(mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(mat))
  sample_cols <- setdiff(names(mat), c(taxon_col, meta_cols))
  tibble::tibble(`Sample alias` = sample_cols,
                 depth = colSums(mat[, sample_cols, drop = FALSE]),
                 marker = marker)
}

build_prepost_depth_plot <- function(gold_qc_16s, gold_qc_its,
                                     gold_processed_matrix_16s_species,
                                     gold_processed_matrix_its_species,
                                     style = load_plot_style(),
                                     out_path = file.path(DASHBOARD_DIR, "prepost_qc_depth.png")) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  pre <- dplyr::bind_rows(
    dplyr::transmute(gold_qc_16s, marker = "16S", stage = `Jenis Kebun`,
                     depth = Mapped),
    dplyr::transmute(gold_qc_its, marker = "ITS", stage = `Jenis Kebun`,
                     depth = Mapped)) |>
    dplyr::filter(!is.na(stage), !is.na(depth), depth > 0) |>
    dplyr::mutate(phase = "Pre-QC")

  # NOTE: "Sample alias" (e.g. barcode01_169) is only unique WITHIN one
  # marker's own barcode numbering -- a 16S and an ITS run can independently
  # reuse the same barcode index for a different sample, so this lookup must
  # join on (marker, Sample alias), never Sample alias alone.
  stage_lookup <- dplyr::bind_rows(
    dplyr::transmute(gold_qc_16s, marker = "16S", `Sample alias`,
                     stage = `Jenis Kebun`),
    dplyr::transmute(gold_qc_its, marker = "ITS", `Sample alias`,
                     stage = `Jenis Kebun`))

  post <- dplyr::bind_rows(
    .processed_matrix_depth(gold_processed_matrix_16s_species, "16S"),
    .processed_matrix_depth(gold_processed_matrix_its_species, "ITS")) |>
    dplyr::left_join(stage_lookup, by = c("marker", "Sample alias")) |>
    dplyr::filter(!is.na(stage), depth > 0) |>
    dplyr::mutate(phase = "Post-QC")

  d <- dplyr::bind_rows(dplyr::select(pre, marker, stage, depth, phase),
                        dplyr::select(post, marker, stage, depth, phase))
  d$phase <- factor(d$phase, levels = c("Pre-QC", "Post-QC"))
  d$universe <- paste0(d$marker, "_", d$stage)

  p <- ggplot2::ggplot(d, ggplot2::aes(phase, depth, fill = phase)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = .8, width = .55) +
    ggplot2::geom_jitter(width = .12, size = 1, alpha = .5) +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::scale_fill_manual(values = c("Pre-QC" = "#8DA0CB", "Post-QC" = "#66C2A5"),
                               guide = "none") +
    ggplot2::facet_wrap(~ universe, scales = "free_y") +
    ggplot2::labs(title = "Sample sequencing depth: Pre-QC vs Post-QC",
                  subtitle = "Pre-QC = Mapped reads (gold_qc, before filtering) | Post-QC = species-rank depth after taxonomic + threshold filter",
                  x = NULL, y = "Depth (log scale)") +
    gold_plot_theme(style)

  ggplot2::ggsave(out_path, p, width = 9, height = 7, dpi = 200)
  message("Wrote ", out_path)
  out_path
}

# ----------------------------------------------------------------------------
# Rarefaction curve, per universe — on the TAXONOMICALLY-filtered species
# matrix (gold_matrix, HIGH PASS samples), deliberately BEFORE the OTU
# min-reads/prevalence threshold cut. A rarefaction curve exists to show
# whether depth is sufficient to have discovered rare taxa; running it on
# data that's already had rare OTUs pruned would artificially flatten the
# curve. Mirrors the established pipeline's functions_qc.R::plot_rarefaction(),
# which likewise uses tax_filtered (pre-threshold), not the thresholded table.
# `step` is scaled to each sample's own depth (~target_points evaluations per
# curve) rather than fixed, since 16S depths (~1M) vs ITS (~10-50k) differ by
# two orders of magnitude and a fixed small step is very slow at high depth.
# ----------------------------------------------------------------------------
plot_rarefaction_curve <- function(gold_matrix_species, universe_lookup,
                                   rarefied_universe, marker, stage,
                                   out_dir = DASHBOARD_DIR, target_points = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir,
    paste0("rarefaction_", tolower(marker), "_", tolower(stage), ".png"))

  tf <- taxonomic_filter_gold(gold_matrix_species, marker)
  sm <- .gold_matrix_to_sample_matrix(tf)

  ids <- intersect(universe_lookup$`Sample alias`[universe_lookup$stage == stage],
                   rownames(sm))
  sub <- sm[ids, , drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  if (nrow(sub) == 0) {
    message("  rarefaction curve skipped for ", marker, "/", stage, ": no samples")
    return(NA_character_)
  }
  storage.mode(sub) <- "integer"

  rdepth <- if (!is.null(rarefied_universe) && nrow(rarefied_universe) > 0)
    unique(rowSums(rarefied_universe))[1] else min(rowSums(sub))
  step <- max(200, floor(max(rowSums(sub)) / target_points))

  grDevices::png(out_path, width = 1400, height = 900, res = 150)
  vegan::rarecurve(sub, step = step, label = FALSE, col = "#3D8C8C",
                   xlab = "Sequencing depth", ylab = "Observed species",
                   main = paste0(marker, "_", stage, " — rarefaction curve"))
  graphics::abline(v = rdepth, lty = 2, col = "firebrick", lwd = 2)
  grDevices::dev.off()
  message("  Wrote ", out_path, " (rarefy depth line at ",
          format(rdepth, big.mark = ","), ", step=", step, ")")
  out_path
}

build_rarefaction_curves <- function(gold_matrix_16s_species,
                                     gold_universe_lookup_16s,
                                     gold_rarefied_species_16s,
                                     gold_matrix_its_species,
                                     gold_universe_lookup_its,
                                     gold_rarefied_species_its) {
  written <- c(
    purrr::map_chr(names(gold_rarefied_species_16s), ~ plot_rarefaction_curve(
      gold_matrix_16s_species, gold_universe_lookup_16s,
      gold_rarefied_species_16s[[.x]], "16S", .x)),
    purrr::map_chr(names(gold_rarefied_species_its), ~ plot_rarefaction_curve(
      gold_matrix_its_species, gold_universe_lookup_its,
      gold_rarefied_species_its[[.x]], "ITS", .x)))
  written[!is.na(written)]
}

# ----------------------------------------------------------------------------
# PCoA (Bray-Curtis, rarefied species matrix), per universe.
#   colour = Waktu (time), shape = Kode Kebun (field)
# ----------------------------------------------------------------------------
plot_pcoa_universe <- function(rarefied_mat, universe_lookup, marker, stage,
                               style = load_plot_style(), out_dir = DASHBOARD_DIR) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir,
    paste0("pcoa_", tolower(marker), "_", tolower(stage), ".png"))
  if (is.null(rarefied_mat) || nrow(rarefied_mat) < 3) {
    message("  PCoA skipped for ", marker, "/", stage, ": <3 samples")
    return(NA_character_)
  }

  d <- vegan::vegdist(rarefied_mat, method = "bray")
  pcoa <- stats::cmdscale(d, k = 2, eig = TRUE)
  scores <- as.data.frame(pcoa$points); names(scores) <- c("Axis1", "Axis2")
  scores$`Sample alias` <- rownames(rarefied_mat)
  eig <- pcoa$eig[pcoa$eig > 0]
  ve <- round(100 * (pcoa$eig / sum(eig))[1:2], 1)

  meta <- universe_lookup[universe_lookup$stage == stage, ]
  df <- dplyr::left_join(scores, meta, by = "Sample alias")

  fields <- sort(unique(df$field))
  shapes <- rep(c(16, 17, 15, 18, 3, 7, 8, 9, 10, 11, 12, 13, 14),
               length.out = length(fields))

  p <- ggplot2::ggplot(df, ggplot2::aes(Axis1, Axis2, colour = waktu, shape = field)) +
    ggplot2::geom_point(size = 3, alpha = .85) +
    ggplot2::scale_shape_manual(values = stats::setNames(shapes, fields),
                                name = "Kode Kebun") +
    ggplot2::scale_colour_manual(values = waktu_palette(), name = "Waktu") +
    ggplot2::labs(
      title = paste0(marker, "_", stage, " — Bray-Curtis PCoA (species, rarefied)"),
      x = paste0("PCoA Axis 1 (", ve[1], "%)"),
      y = paste0("PCoA Axis 2 (", ve[2], "%)")) +
    gold_plot_theme(style)

  ggplot2::ggsave(out_path, p, width = 8, height = 6.5, dpi = 200)
  message("  Wrote ", out_path)
  out_path
}

build_pcoa_plots <- function(gold_rarefied_species_16s, gold_universe_lookup_16s,
                             gold_rarefied_species_its, gold_universe_lookup_its,
                             style = load_plot_style()) {
  written <- c(
    purrr::map_chr(names(gold_rarefied_species_16s), ~ plot_pcoa_universe(
      gold_rarefied_species_16s[[.x]], gold_universe_lookup_16s, "16S", .x, style)),
    purrr::map_chr(names(gold_rarefied_species_its), ~ plot_pcoa_universe(
      gold_rarefied_species_its[[.x]], gold_universe_lookup_its, "ITS", .x, style)))
  written[!is.na(written)]
}

# ----------------------------------------------------------------------------
# Shannon diversity barplot, per universe (16S_Nursery, 16S_TM, ITS_TM only --
# ITS_Nursery has no data, matching the rest of the pipeline's convention).
#   x = Waktu, dodged bars = Kode Pupuk (fertilizer), bar height = mean
#   Shannon, error bar = +/- 1 SD (drawn only when n>1). y-axis fixed 0-5.
# ----------------------------------------------------------------------------
build_shannon_dashboard_barplot <- function(alpha_shannon, style = load_plot_style(),
                                            out_path = file.path(DASHBOARD_DIR,
                                              "shannon_barplot.png")) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  d <- alpha_shannon |>
    dplyr::mutate(universe = paste0(marker, "_", stage)) |>
    dplyr::filter(universe %in% c("16S_Nursery", "16S_TM", "ITS_TM"),
                  !is.na(fertilizer), !is.na(waktu))
  if (nrow(d) == 0) {
    message("Shannon dashboard barplot skipped: no data"); return(NA_character_)
  }

  summ <- d |>
    dplyr::group_by(universe, waktu, fertilizer) |>
    dplyr::summarise(mean_shannon = mean(Shannon, na.rm = TRUE),
                     sd_shannon = stats::sd(Shannon, na.rm = TRUE),
                     n = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(sd_shannon = ifelse(n > 1, sd_shannon, NA_real_))

  p <- ggplot2::ggplot(summ, ggplot2::aes(waktu, mean_shannon, fill = fertilizer)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(.8), width = .7,
                      colour = "grey30", linewidth = .2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_shannon - sd_shannon, ymax = mean_shannon + sd_shannon),
      position = ggplot2::position_dodge(.8), width = .2, na.rm = TRUE) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +   # flexible top, fixed 0 floor
    ggplot2::facet_wrap(~ universe, nrow = 1) +
    ggplot2::labs(title = "Shannon diversity by Waktu and Kode Pupuk",
                  subtitle = "Bars = mean Shannon; error bars = ± 1 SD (shown when n>1)",
                  x = "Waktu", y = "Shannon index", fill = "Kode Pupuk") +
    gold_plot_theme(style)

  ggplot2::ggsave(out_path, p, width = 13, height = 6, dpi = 200)
  message("Wrote ", out_path)
  out_path
}
