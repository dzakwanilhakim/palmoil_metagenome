# =============================================================================
# R/functions_gold_alpha.R  —  Analysis output: alpha diversity
#   Shannon only, from the rarefied species matrix, per universe (marker x
#   stage). Per project spec: "Alpha diversity using only rarefied species
#   matrix with output only Shannon diversity."
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(yaml)
})

# shared by alpha here and by the ANCOM-BC DA/co-occurrence steps that follow
load_analysis_thresholds <- function(path = "config/analysis_thresholds.yaml") {
  yaml::read_yaml(path)
}

compute_shannon_universe <- function(rarefied_universe) {
  if (is.null(rarefied_universe) || nrow(rarefied_universe) == 0)
    return(tibble::tibble())
  tibble::tibble(`Sample alias` = rownames(rarefied_universe),
                 Shannon = vegan::diversity(rarefied_universe, index = "shannon"))
}

build_alpha_shannon <- function(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                gold_rarefied_species_its, gold_universe_lookup_its) {
  one_marker <- function(rarefied_list, lookup, marker) {
    purrr::map_dfr(names(rarefied_list), function(st) {
      compute_shannon_universe(rarefied_list[[st]]) |>
        dplyr::mutate(marker = marker, stage = st)
    }) |>
      dplyr::left_join(lookup, by = c("Sample alias", "stage"))
  }
  dplyr::bind_rows(
    one_marker(gold_rarefied_species_16s, gold_universe_lookup_16s, "16S"),
    one_marker(gold_rarefied_species_its, gold_universe_lookup_its, "ITS"))
}

plot_alpha_shannon <- function(alpha_df, style = load_plot_style(),
                               out_path = "Results/analysis/alpha_shannon.png") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  alpha_df <- dplyr::mutate(alpha_df, universe = paste0(marker, "_", stage))

  p <- ggplot2::ggplot(alpha_df, ggplot2::aes(fertilizer, Shannon, fill = fertilizer)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = .75) +
    ggplot2::geom_jitter(width = .15, size = 1.3, alpha = .6) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +   # flexible top, fixed 0 floor
    ggplot2::facet_grid(universe ~ waktu, scales = "free_x") +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::labs(title = "Shannon diversity (species, rarefied)",
                  subtitle = "x = Kode Pupuk (fertilizer), facet rows = universe, facet cols = Waktu (time)",
                  x = "Kode Pupuk (fertilizer)", y = "Shannon index") +
    gold_plot_theme(style) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(out_path, p, width = 12, height = 10, dpi = 200, limitsize = FALSE)
  message("Wrote ", out_path)
  out_path
}

write_alpha_shannon_csv <- function(alpha_df,
                                    out_path = "Results/analysis/alpha_shannon.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(alpha_df, out_path)
  message("Wrote ", out_path)
  out_path
}
