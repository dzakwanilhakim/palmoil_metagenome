# =============================================================================
# R/functions_alpha_plots.R — Stage 4: alpha diversity visualisations
#   A & C (snapshots) -> boxplots of metric by group_by
#   B & D (longitudinal) -> trajectory plots: individual points (incl. unpaired
#         standalone T0) + line connecting the MEDIAN of each timepoint group
#   Three metrics each: Observed, Pielou, Shannon. Per marker.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

.alpha_long <- function(alpha_df, marker) {
  dplyr::filter(alpha_df, marker == !!marker) |>
    tidyr::pivot_longer(c(Observed, Pielou, Shannon),
                        names_to = "metric", values_to = "value")
}

PAL_FERT <- c(PH1="#1b9e77", PH2="#d95f02", PH3="#7570b3",
              PH4="#e7298a", PH5="#66a61e", PH6="#e6ab02")

# ----------------------------------------------------------------------------
# A & C — SNAPSHOT BOXPLOTS
#   A: facet per (field, timepoint), x = fertilizer
#   C: facet per timepoint, x = fertilizer (pooled across fields);
#      plus a field-faceted variant (C_per_field)
# ----------------------------------------------------------------------------
plot_snapshot_box <- function(alpha_df, marker, goal,
                              out_dir = "results/alpha") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  al <- .alpha_long(alpha_df, marker)
  if (nrow(al) == 0) return(NA_character_)

  if (goal == "A") {
    al <- dplyr::mutate(al, panel = paste(field, timepoint, sep = " / "))
    facet <- ggplot2::facet_grid(metric ~ panel, scales = "free_y")
    sub <- "Goal A — fertilizers within each field x timepoint"
  } else { # C
    al <- dplyr::mutate(al, panel = timepoint)
    facet <- ggplot2::facet_grid(metric ~ panel, scales = "free_y")
    sub <- "Goal C — fertilizers across fields, per timepoint (pooled)"
  }

  p <- ggplot2::ggplot(al, ggplot2::aes(fertilizer, value, fill = fertilizer)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = .7, linewidth = .3) +
    ggplot2::geom_jitter(width = .15, size = 1.1, alpha = .6, shape = 21,
                         colour = "grey20") +
    ggplot2::scale_fill_manual(values = PAL_FERT, guide = "none") +
    facet +
    ggplot2::labs(title = paste0(marker, " — alpha diversity"),
                  subtitle = sub, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  n_panel <- dplyr::n_distinct(al$panel)
  out <- file.path(out_dir, paste0("alpha_box_goal", goal, "_", tolower(marker), ".png"))
  ggplot2::ggsave(out, p, width = max(7, 1.6 * n_panel), height = 8,
                  dpi = 200, limitsize = FALSE)
  message("Wrote ", out); out
}

# ----------------------------------------------------------------------------
# B & D — TRAJECTORY PLOTS
#   individual points (all samples, incl. unpaired standalone T0) +
#   line connecting the MEDIAN of each timepoint group.
#   B: facet per field. D: facet per field (same layout; D is "all fields").
#   shape = fertilizer; one figure per metric.
# ----------------------------------------------------------------------------
plot_trajectory <- function(alpha_df, marker, goal, metric,
                            out_dir = "results/alpha") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- dplyr::filter(alpha_df, marker == !!marker)
  if (nrow(d) == 0) return(NA_character_)
  d$value <- d[[metric]]
  d$timepoint <- factor(d$timepoint, levels = sort(unique(d$timepoint)))

  # median per (field, timepoint) for the trend line (D & B both facet by field)
  med <- d |>
    dplyr::group_by(field, timepoint) |>
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(d, ggplot2::aes(timepoint, value)) +
    ggplot2::geom_point(ggplot2::aes(shape = fertilizer, colour = fertilizer),
                        size = 2.2, alpha = .85,
                        position = ggplot2::position_jitter(width = .06, height = 0)) +
    ggplot2::geom_line(data = med, ggplot2::aes(timepoint, med, group = 1),
                       colour = "grey20", linewidth = .7) +
    ggplot2::geom_point(data = med, ggplot2::aes(timepoint, med),
                        colour = "grey20", size = 2.4) +
    ggplot2::scale_colour_manual(values = PAL_FERT, name = "Fertilizer") +
    ggplot2::scale_shape_manual(values = c(PH1=16,PH2=17,PH3=15,
                                           PH4=18,PH5=3,PH6=8), name = "Fertilizer") +
    ggplot2::facet_wrap(~ field, scales = "free_y") +
    ggplot2::labs(
      title = paste0(marker, " — ", metric, " trajectory (Goal ", goal, ")"),
      subtitle = "Points = samples (unpaired T0 shown standalone); line = median per timepoint",
      x = "Timepoint", y = metric) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right")

  n_field <- dplyr::n_distinct(d$field)
  out <- file.path(out_dir,
    paste0("alpha_traj_goal", goal, "_", tolower(metric), "_", tolower(marker), ".png"))
  ggplot2::ggsave(out, p, width = max(7, 2.2 * ceiling(sqrt(n_field))),
                  height = max(5, 2.0 * ceiling(sqrt(n_field))),
                  dpi = 200, limitsize = FALSE)
  message("Wrote ", out); out
}

# convenience driver: all alpha plots for one marker
make_alpha_plots <- function(alpha_df, marker, metrics = c("Observed","Pielou","Shannon")) {
  out <- c(
    plot_snapshot_box(alpha_df, marker, "A"),
    plot_snapshot_box(alpha_df, marker, "C"),
    purrr::map_chr(metrics, ~ plot_trajectory(alpha_df, marker, "B", .x)),
    purrr::map_chr(metrics, ~ plot_trajectory(alpha_df, marker, "D", .x))
  )
  out[!is.na(out)]
}
