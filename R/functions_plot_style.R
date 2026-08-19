# =============================================================================
# R/functions_plot_style.R  —  shared ggplot theme for Dashboard + Analysis
#   plots, driven by config/plot_style.yaml (font size, bold titles/axes).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

load_plot_style <- function(path = "config/plot_style.yaml") {
  yaml::read_yaml(path)
}

# Global time (Waktu) colour palette, keyed off the FULL possible timepoint
# range (TIMEPOINTS_ALL = T0..T4, from functions_counts.R) rather than
# whichever levels happen to be present in one specific plot's subset --
# so T0/T1/etc. always map to the SAME colour everywhere in the pipeline
# (Dashboard PCoA, dendrograms, Goal B/D beta plots, any time-coloured
# alpha plot), instead of each plot independently re-deriving a viridis
# scale over just its own observed levels.
waktu_palette <- function() {
  stats::setNames(viridisLite::viridis(length(TIMEPOINTS_ALL)), TIMEPOINTS_ALL)
}

gold_plot_theme <- function(style = load_plot_style()) {
  bs <- style$base_size %||% 14
  face_title <- if (isTRUE(style$bold_titles)) "bold" else "plain"
  face_axis  <- if (isTRUE(style$bold_axis_text)) "bold" else "plain"

  ggplot2::theme_bw(base_size = bs) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = face_title),
      plot.subtitle = ggplot2::element_text(face = face_title, colour = "grey35"),
      axis.title    = ggplot2::element_text(face = face_axis),
      axis.text     = ggplot2::element_text(face = face_axis),
      strip.text    = ggplot2::element_text(face = face_axis),
      legend.title  = ggplot2::element_text(face = face_axis),
      legend.text   = ggplot2::element_text(face = face_axis))
}
