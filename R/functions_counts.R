# =============================================================================
# R/functions_counts.R — replicate-count tables per field (post-QC), as PNGs
#   Rows = fertilizers (TM: PH1-PH6; Nursery: PHN01-PHN13)
#   Cols = timepoints T0..T6 (T2..T6 forced to exist, pre-filled 0)
#   Cells = number of clean biological replicates (samples) after QC.
#   Zero cells -> red font. Routed to Results/data_counts/<universe>/<field>.png
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# fixed scaffolds (future-proofed)
TIMEPOINTS_ALL <- paste0("T", 0:4)
FERT_TM        <- paste0("PH", 1:6)
FERT_NURSERY   <- sprintf("PHN%02d", 1:13)

fert_rows_for_stage <- function(stage) {
  if (toupper(stage) == "NURSERY") FERT_NURSERY else FERT_TM
}

# build a complete fertilizer x timepoint count grid for one field
field_count_grid <- function(master_samples, marker, stage, field) {
  ferts <- fert_rows_for_stage(stage)

  obs <- master_samples |>
    dplyr::filter(marker == !!marker, stage == !!stage, field == !!field) |>
    dplyr::count(fertilizer, timepoint, name = "n")

  grid <- tidyr::expand_grid(fertilizer = ferts, timepoint = TIMEPOINTS_ALL) |>
    dplyr::left_join(obs, by = c("fertilizer", "timepoint")) |>
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) |>
    tidyr::pivot_wider(names_from = timepoint, values_from = n)

  # force all T0..T6 columns to exist (missing -> 0), in order
  for (tp in TIMEPOINTS_ALL) if (!tp %in% names(grid)) grid[[tp]] <- 0L
  grid <- grid[, c("fertilizer", TIMEPOINTS_ALL)]
  grid <- grid[match(ferts, grid$fertilizer), , drop = FALSE]
  grid$fertilizer <- ferts
  grid[is.na(grid)] <- 0L
  grid
}

# render one styled table as a PNG using pure ggplot2 (no V8/browser needed).
# red bold font where value == 0; dark blue otherwise. Header + title.
render_count_png <- function(grid, title, subtitle, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  ferts <- grid$fertilizer
  long <- grid |>
    tidyr::pivot_longer(-fertilizer, names_to = "timepoint", values_to = "n") |>
    dplyr::mutate(
      fertilizer = factor(fertilizer, levels = rev(ferts)),   # top row first
      timepoint  = factor(timepoint, levels = TIMEPOINTS_ALL),
      txt_col    = ifelse(n == 0, "#C0392B", "#1B4F72"),
      txt_face   = ifelse(n == 0, "bold", "plain"))

  ncol <- length(TIMEPOINTS_ALL); nrow <- length(ferts)

  p <- ggplot2::ggplot(long, ggplot2::aes(timepoint, fertilizer)) +
    ggplot2::geom_tile(fill = "white", colour = "grey80", linewidth = .4) +
    ggplot2::geom_text(ggplot2::aes(label = n, colour = txt_col,
                                    fontface = txt_face), size = 4) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 0.9) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9))

  ggplot2::ggsave(out_path, p,
                  width = 1.0 + 0.7 * ncol,
                  height = 1.4 + 0.45 * nrow, dpi = 200, limitsize = FALSE)
  out_path
}

# ============================================================================
# MAIN: one table per field, per universe, into Results/data_counts/<universe>/
# ============================================================================
build_count_tables <- function(master_samples, root = "Results") {
  universes <- master_samples |> dplyr::distinct(marker, stage) |>
    dplyr::filter(!is.na(stage)) |> dplyr::arrange(marker, stage)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    outdir <- file.path(root, "data_counts", ucode)

    fields <- master_samples |>
      dplyr::filter(marker == mk, stage == st) |>
      dplyr::distinct(field) |> dplyr::pull(field) |> sort()
    if (length(fields) == 0) next
    message("=== COUNTS ", ucode, " (", length(fields), " fields) ===")

    for (fld in fields) {
      grid <- field_count_grid(master_samples, mk, st, fld)
      out  <- file.path(outdir, paste0(fld, "_counts.png"))
      w <- tryCatch(
        render_count_png(grid,
          title = paste0(ucode, " — ", fld, " — replicate counts (post-QC)"),
          subtitle = "Clean biological replicates per fertilizer x timepoint",
          out_path = out),
        error = function(e) { message("count table failed for ", fld, ": ",
                                       conditionMessage(e)); NA_character_ })
      if (!is.na(w)) written <- c(written, w)
    }
  }
  message("Count tables complete: ", length(written), " PNGs.")
  written
}
