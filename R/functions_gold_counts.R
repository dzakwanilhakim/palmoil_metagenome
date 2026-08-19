# =============================================================================
# R/functions_gold_counts.R  —  Synchronized: data_counts (gold-tier)
#   Replicate-count tables per field (post gold-tier QC), from
#   gold_universe_lookup_{16s,its} (Kode Pupuk x Waktu grid) instead of the
#   established pipeline's master_samples. Reuses render_count_png() and the
#   TIMEPOINTS_ALL/fert_rows_for_stage() scaffolding from functions_counts.R
#   unchanged -- only the grid-building step needs a gold-native source.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

gold_field_count_grid <- function(universe_lookup, stage, field) {
  ferts <- fert_rows_for_stage(stage)

  obs <- universe_lookup |>
    dplyr::filter(stage == !!stage, field == !!field) |>
    dplyr::count(fertilizer, waktu, name = "n")

  grid <- tidyr::expand_grid(fertilizer = ferts, waktu = TIMEPOINTS_ALL) |>
    dplyr::left_join(obs, by = c("fertilizer", "waktu")) |>
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) |>
    tidyr::pivot_wider(names_from = waktu, values_from = n)

  for (tp in TIMEPOINTS_ALL) if (!tp %in% names(grid)) grid[[tp]] <- 0L
  grid <- grid[, c("fertilizer", TIMEPOINTS_ALL)]
  grid <- grid[match(ferts, grid$fertilizer), , drop = FALSE]
  grid$fertilizer <- ferts
  grid[is.na(grid)] <- 0L
  grid
}

build_gold_count_tables <- function(gold_universe_lookup_16s, gold_universe_lookup_its,
                                    root = "Results") {
  lookups <- list("16S" = gold_universe_lookup_16s, "ITS" = gold_universe_lookup_its)
  written <- character(0)

  for (marker in names(lookups)) {
    lookup <- lookups[[marker]]
    stages <- sort(unique(lookup$stage))
    for (st in stages) {
      ucode <- paste0(marker, "_", st)
      outdir <- file.path(root, "data_counts", ucode)
      fields <- lookup |> dplyr::filter(stage == st) |>
        dplyr::distinct(field) |> dplyr::pull(field) |> sort()
      if (length(fields) == 0) next
      message("=== GOLD COUNTS ", ucode, " (", length(fields), " fields) ===")

      for (fld in fields) {
        grid <- gold_field_count_grid(lookup, st, fld)
        out  <- file.path(outdir, paste0(fld, "_counts.png"))
        w <- tryCatch(
          render_count_png(grid,
            title = paste0(ucode, " — ", fld, " — replicate counts (gold, HIGH PASS)"),
            subtitle = "Clean biological replicates per Kode Pupuk x Waktu",
            out_path = out),
          error = function(e) { message("count table failed for ", fld, ": ",
                                        conditionMessage(e)); NA_character_ })
        if (!is.na(w)) written <- c(written, w)
      }
    }
  }
  message("Gold count tables complete: ", length(written), " PNGs.")
  written
}
