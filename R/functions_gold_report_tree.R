# =============================================================================
# R/functions_gold_report_tree.R  —  Synchronized: Goal A/B/C/D alpha report
#   Shannon-only, species rank, sourced from alpha_shannon (functions_gold_
#   alpha.R) instead of the established pipeline's alpha_div/master_samples.
#   Reuses GOAL_DIR, fert_palette(), list_universes(), universe_alpha_stats()
#   /run_alpha_stats()/expand_goal() (functions_alpha.R + functions_report_
#   tree.R) completely unchanged, and the plot_goalA_field()/
#   plot_trajectory_metric()/plot_goalC_pooled()/plot_goalD_pooled() plotters
#   via their new metrics/style parameters (generalized, not duplicated).
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# reshape gold alpha_shannon (Sample alias/stage/field/fertilizer/waktu/
# Shannon) into the column names the established report-tree engine expects
# (timepoint instead of waktu; id_sampel for parity/future-proofing).
gold_alpha_for_report <- function(alpha_shannon) {
  alpha_shannon |>
    dplyr::rename(id_sampel = `Sample alias`, timepoint = waktu) |>
    dplyr::filter(!is.na(stage), !is.na(field), !is.na(fertilizer), !is.na(timepoint))
}

build_gold_report_tree <- function(alpha_shannon, comparisons,
                                   style = load_plot_style(), root = "Results") {
  d_all <- gold_alpha_for_report(alpha_shannon)
  universes <- list_universes(d_all)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)
    d_u   <- dplyr::filter(d_all, marker == mk, stage == st)
    if (nrow(d_u) == 0) next
    pal   <- fert_palette(d_u$fertilizer)
    message("=== GOLD Universe ", ucode, " (", nrow(d_u), " samples, Shannon only) ===")

    stats_u <- universe_alpha_stats(d_all, comparisons, mk, st)
    fields  <- sort(unique(d_u$field))

    ## Goal A — one combined-timepoint figure per field
    for (fld in fields) {
      leaf <- file.path(udir, GOAL_DIR["A"], fld)
      dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
      dsub <- dplyr::filter(d_u, field == fld)
      w <- plot_goalA_field(dsub, pal,
            title = paste0(ucode, " — Goal A — ", fld),
            out_path = file.path(leaf, paste0("alpha_box_", fld, ".png")),
            metrics = "Shannon", style = style,
            chart_type = "bar", y_limits = c(0, NA))
      if (!is.na(w)) written <- c(written, w)
      sl <- dplyr::filter(stats_u$global, goal == "A",
                          grepl(paste0("field=", fld, "(;|$)"), context))
      readr::write_csv(sl, file.path(leaf, paste0("stats_A_", fld, ".csv")))
    }

    ## Goal B — per field, Shannon trajectory + p-table
    for (fld in fields) {
      leaf <- file.path(udir, GOAL_DIR["B"], fld)
      dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
      dsub <- dplyr::filter(d_u, field == fld)
      w <- plot_trajectory_metric(dsub, "Shannon", pal,
            title = paste0(ucode, " — Goal B — ", fld, " — Shannon"),
            out_path = file.path(leaf, paste0("alpha_traj_", fld, "_shannon.png")),
            facet_field = FALSE, style = style, y_limits = c(0, NA))
      if (!is.na(w)) written <- c(written, w)
      sl <- dplyr::filter(stats_u$global, goal %in% c("B", "B_per_fertilizer"),
                          grepl(paste0("field=", fld, "(;| |$)"), context))
      readr::write_csv(sl, file.path(leaf, paste0("stats_B_", fld, ".csv")))
    }

    ## Goal C — pooled cross-field snapshot (no field subfolders)
    leafC <- file.path(udir, GOAL_DIR["C"])
    dir.create(leafC, recursive = TRUE, showWarnings = FALSE)
    w <- plot_goalC_pooled(d_u, pal,
          title = paste0(ucode, " — Goal C — fertilizers across fields"),
          out_path = file.path(leafC, "alpha_box_crossfield.png"),
          metrics = "Shannon", style = style,
          chart_type = "bar", y_limits = c(0, NA))
    if (!is.na(w)) written <- c(written, w)
    readr::write_csv(dplyr::filter(stats_u$global, goal %in% c("C", "C_per_field")),
                     file.path(leafC, "stats_C.csv"))

    ## Goal D — faceted-by-field trajectory + master pooled plot
    leafD <- file.path(udir, GOAL_DIR["D"])
    dir.create(leafD, recursive = TRUE, showWarnings = FALSE)
    w1 <- plot_trajectory_metric(d_u, "Shannon", pal,
          title = paste0(ucode, " — Goal D (faceted) — Shannon"),
          out_path = file.path(leafD, "alpha_traj_faceted_shannon.png"),
          facet_field = TRUE, style = style, y_limits = c(0, NA))
    if (!is.na(w1)) written <- c(written, w1)
    w2 <- plot_goalD_pooled(d_u, "Shannon", pal,
          title = paste0(ucode, " — Goal D (pooled) — Shannon"),
          out_path = file.path(leafD, "alpha_traj_pooled_shannon.png"),
          style = style, y_limits = c(0, NA))
    if (!is.na(w2)) written <- c(written, w2)
    readr::write_csv(dplyr::filter(stats_u$global, goal %in% c("D", "D_per_fertilizer")),
                     file.path(leafD, "stats_D.csv"))
  }
  message("Gold report tree complete: ", length(written), " plots under ", root, "/")
  written
}
