# =============================================================================
# R/functions_report_tree.R — Stage 4 output (per-universe depth-first tree)
#   4 isolated universes: {marker} x {stage}. Empty universes skipped.
#   Results/<MARKER>_<STAGE>/Goal_{A,B,C,D}.../[FIELD]/  plots + sliced CSVs
#
#   - dynamic fertilizer palette via scales::hue_pal() (handles >=13)
#   - n=1 groups -> single dot (no flat boxplot) GLOBALLY
#   - trajectory lines = MEDIAN per fertilizer per timepoint (colored)
#   - embedded dynamic p-value tables (consecutive-timepoint raw Wilcoxon,
#     rows = Overall + each fertilizer, cols = T0vsT1, T1vsT2, ...)
#   - Goal D adds a master POOLED plot: median +/- IQR per fertilizer over time
#   - significance from RAW p<0.05 (asterisks); CSVs hold raw + BH
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

GOAL_DIR <- c(
  A = "Goal_A_Intra_Snapshot",
  B = "Goal_B_Intra_Longitudinal",
  C = "Goal_C_Cross_Snapshot",
  D = "Goal_D_Cross_Longitudinal")

MIN_N <- 2  # global stats gate

sig_stars <- function(p) dplyr::case_when(
  is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
  p < 0.05 ~ "*", TRUE ~ "ns")

# dynamic categorical palette for any number of fertilizers
fert_palette <- function(levels) {
  levels <- sort(unique(levels))
  setNames(scales::hue_pal()(length(levels)), levels)
}

list_universes <- function(alpha_df) {
  alpha_df |> dplyr::distinct(marker, stage) |>
    dplyr::filter(!is.na(stage)) |> dplyr::arrange(marker, stage)
}

.alpha_long_u <- function(d) {
  tidyr::pivot_longer(d, c(Observed, Pielou, Shannon),
                      names_to = "metric", values_to = "value")
}

# add a flag: groups with n==1 plotted as dots, n>=2 as boxplots
.box_or_dot <- function(al, group_cols) {
  al |> dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::mutate(.n = dplyr::n()) |> dplyr::ungroup()
}

# overall fertilizer test (raw p) per metric on a subset
.overall_fert_p <- function(d, min_n = MIN_N) {
  purrr::map_dfr(c("Observed","Pielou","Shannon"), function(m){
    g <- as.factor(d$fertilizer); v <- d[[m]]
    keep <- !is.na(v); v <- v[keep]; g <- droplevels(g[keep])
    tab <- table(g); p <- NA_real_
    if (length(tab) >= 2 && all(tab >= min_n))
      p <- if (nlevels(g)==2) suppressWarnings(stats::wilcox.test(v~g)$p.value)
           else stats::kruskal.test(v~g)$p.value
    tibble::tibble(metric = m, p = p)
  })
}

# ---- DYNAMIC P-VALUE TABLE -------------------------------------------------
# consecutive-timepoint raw unpaired Wilcoxon per fertilizer + Overall row.
# returns a tibble: row (Overall/PH..) x each transition column.
pvalue_transition_table <- function(d, metric, min_n = MIN_N) {
  d <- d[!is.na(d[[metric]]), ]
  tps <- sort(unique(d$timepoint))
  if (length(tps) < 2) return(NULL)
  transitions <- paste0(tps[-length(tps)], "vs", tps[-1])

  test_pair <- function(sub, t1, t2) {
    s <- sub[sub$timepoint %in% c(t1,t2), ]
    tab <- table(s$timepoint)
    if (length(tab) < 2 || any(tab < min_n)) return("ns")
    p <- suppressWarnings(stats::wilcox.test(s[[metric]] ~ s$timepoint)$p.value)
    if (is.na(p)) "ns" else paste0(signif(p,2), sig_stars(p))
  }

  rows <- c("Overall", sort(unique(d$fertilizer)))
  mat <- lapply(rows, function(r){
    sub <- if (r == "Overall") d else d[d$fertilizer == r, ]
    vapply(seq_along(transitions), function(i)
      test_pair(sub, tps[i], tps[i+1]), character(1))
  })
  out <- as.data.frame(do.call(rbind, mat), stringsAsFactors = FALSE)
  names(out) <- transitions
  tibble::tibble(group = rows) |> dplyr::bind_cols(out)
}

# render a p-value table tibble as a ggplot grob (for patchwork)
pvalue_table_grob <- function(tbl, title = "Raw Wilcoxon p (consecutive timepoints)") {
  if (is.null(tbl) || nrow(tbl) == 0)
    return(patchwork::wrap_elements(grid::textGrob("No temporal pairs available")))
  long <- tidyr::pivot_longer(tbl, -group, names_to = "transition", values_to = "p")
  long$group <- factor(long$group, levels = rev(tbl$group))
  long$transition <- factor(long$transition, levels = names(tbl)[-1])
  ggplot2::ggplot(long, ggplot2::aes(transition, group)) +
    ggplot2::geom_tile(fill = "grey97", colour = "grey80") +
    ggplot2::geom_text(ggplot2::aes(label = p), size = 3) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 9, face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# ---- GOAL A: combined-timepoint snapshot, one figure per field -------------
# 3 rows (metric) x dynamic cols (timepoint); n=1 -> dot only.
plot_goalA_field <- function(d, pal, title, out_path, min_n = MIN_N) {
  al <- .alpha_long_u(d)
  if (nrow(al) == 0) return(NA_character_)
  al <- .box_or_dot(al, c("metric","timepoint","fertilizer"))
  box <- dplyr::filter(al, .n >= 2); dot <- dplyr::filter(al, .n == 1)

  # overall fert p per (metric, timepoint)
  pann <- al |> dplyr::distinct(metric, timepoint) |>
    purrr::pmap_dfr(function(metric, timepoint){
      sub <- dplyr::filter(d, timepoint == !!timepoint)
      pp <- .overall_fert_p(sub, min_n) |> dplyr::filter(metric == !!metric)
      tibble::tibble(metric=metric, timepoint=timepoint,
        lab = ifelse(is.na(pp$p), "ns", paste0("p=",signif(pp$p,2)," ",sig_stars(pp$p))))
    })

  p <- ggplot2::ggplot(mapping = ggplot2::aes(fertilizer, value, fill = fertilizer)) +
    { if (nrow(box)) ggplot2::geom_boxplot(data = box, outlier.shape = NA,
        alpha=.7, linewidth=.3) } +
    { if (nrow(box)) ggplot2::geom_jitter(data = box, width=.12, size=1.1,
        alpha=.6, shape=21, colour="grey20") } +
    { if (nrow(dot)) ggplot2::geom_point(data = dot, size=2.4, shape=21,
        colour="grey20") } +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::facet_grid(metric ~ timepoint, scales = "free_y") +
    ggplot2::geom_text(data = pann, ggplot2::aes(x=Inf, y=Inf, label=lab),
        inherit.aes=FALSE, hjust=1.05, vjust=1.4, size=2.8) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face="bold"),
        axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  n_tp <- dplyr::n_distinct(al$timepoint)
  ggplot2::ggsave(out_path, p, width = max(7, 3.2*n_tp), height = 8, dpi = 200,
                  limitsize = FALSE)
  out_path
}

# ---- GOAL B/D-faceted: colored median trajectory + embedded p-table --------
# one figure per metric. `facet_field` TRUE => facet_wrap(~field) (Goal D faceted)
plot_trajectory_metric <- function(d, metric, pal, title, out_path,
                                    facet_field = FALSE, min_n = MIN_N) {
  if (nrow(d) == 0) return(NA_character_)
  d$value <- d[[metric]]
  d$timepoint <- factor(d$timepoint, levels = sort(unique(d$timepoint)))

  grp <- if (facet_field) c("field","fertilizer","timepoint") else c("fertilizer","timepoint")
  med <- d |> dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(med = median(value, na.rm=TRUE), .groups="drop")

  p <- ggplot2::ggplot(d, ggplot2::aes(timepoint, value, colour = fertilizer)) +
    ggplot2::geom_point(ggplot2::aes(shape = fertilizer), size=2, alpha=.8,
        position = ggplot2::position_jitter(width=.05, height=0)) +
    ggplot2::geom_line(data = med, ggplot2::aes(timepoint, med,
        group = fertilizer, colour = fertilizer), linewidth=.8) +
    ggplot2::geom_point(data = med, ggplot2::aes(timepoint, med, colour = fertilizer),
        size=2.2) +
    ggplot2::scale_colour_manual(values = pal, name="Fertilizer") +
    ggplot2::scale_shape_manual(values = rep(c(16,17,15,18,3,8,7,9,10,11,12,13,14),3),
        name="Fertilizer") +
    ggplot2::labs(title = title, x="Timepoint", y=metric) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face="bold"))
  if (facet_field) p <- p + ggplot2::facet_wrap(~ field, scales="free_y")

  # p-value table (pool fields if faceted-D uses overall; here per-plot scope)
  tbl <- pvalue_transition_table(d, metric, min_n)
  grob <- pvalue_table_grob(tbl, paste0("Raw Wilcoxon p — ", metric))

  combined <- p / grob + patchwork::plot_layout(heights = c(3, 1.4))
  nf <- if (facet_field) dplyr::n_distinct(d$field) else 1
  ggplot2::ggsave(out_path, combined,
                  width = max(7, 2.4*ceiling(sqrt(nf))),
                  height = max(7, 2.2*ceiling(sqrt(nf)) + 2), dpi = 200,
                  limitsize = FALSE)
  out_path
}

# ---- GOAL C: pooled cross-field snapshot (boxplots, n=1 -> dot) ------------
plot_goalC_pooled <- function(d, pal, title, out_path, min_n = MIN_N) {
  al <- .alpha_long_u(d)
  if (nrow(al) == 0) return(NA_character_)
  al <- .box_or_dot(al, c("metric","timepoint","fertilizer"))
  box <- dplyr::filter(al, .n >= 2); dot <- dplyr::filter(al, .n == 1)
  pann <- al |> dplyr::distinct(metric, timepoint) |>
    purrr::pmap_dfr(function(metric, timepoint){
      sub <- dplyr::filter(d, timepoint == !!timepoint)
      pp <- .overall_fert_p(sub, min_n) |> dplyr::filter(metric == !!metric)
      tibble::tibble(metric=metric, timepoint=timepoint,
        lab = ifelse(is.na(pp$p),"ns",paste0("p=",signif(pp$p,2)," ",sig_stars(pp$p))))
    })
  p <- ggplot2::ggplot(mapping = ggplot2::aes(fertilizer, value, fill=fertilizer)) +
    { if (nrow(box)) ggplot2::geom_boxplot(data=box, outlier.shape=NA, alpha=.7, linewidth=.3) } +
    { if (nrow(box)) ggplot2::geom_jitter(data=box, width=.12, size=1, alpha=.5, shape=21, colour="grey20") } +
    { if (nrow(dot)) ggplot2::geom_point(data=dot, size=2.2, shape=21, colour="grey20") } +
    ggplot2::scale_fill_manual(values = pal, guide="none") +
    ggplot2::facet_grid(metric ~ timepoint, scales="free_y") +
    ggplot2::geom_text(data=pann, ggplot2::aes(x=Inf, y=Inf, label=lab),
        inherit.aes=FALSE, hjust=1.05, vjust=1.4, size=2.8) +
    ggplot2::labs(title=title, x=NULL, y=NULL) +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
        plot.title=ggplot2::element_text(face="bold"),
        axis.text.x=ggplot2::element_text(angle=45, hjust=1))
  ggplot2::ggsave(out_path, p, width=9, height=8, dpi=200)
  out_path
}

# ---- GOAL D master POOLED plot: median +/- IQR per fertilizer over time ----
plot_goalD_pooled <- function(d, metric, pal, title, out_path, min_n = MIN_N) {
  if (nrow(d) == 0) return(NA_character_)
  d$value <- d[[metric]]
  d$timepoint <- factor(d$timepoint, levels = sort(unique(d$timepoint)))
  summ <- d |> dplyr::group_by(fertilizer, timepoint) |>
    dplyr::summarise(med = median(value, na.rm=TRUE),
                     lo = quantile(value, .25, na.rm=TRUE),
                     hi = quantile(value, .75, na.rm=TRUE), .groups="drop")
  p <- ggplot2::ggplot(summ, ggplot2::aes(timepoint, med, colour=fertilizer, group=fertilizer)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lo, ymax=hi), width=.12, linewidth=.5,
        position = ggplot2::position_dodge(.2)) +
    ggplot2::geom_line(linewidth=.8, position=ggplot2::position_dodge(.2)) +
    ggplot2::geom_point(size=2.4, position=ggplot2::position_dodge(.2)) +
    ggplot2::scale_colour_manual(values=pal, name="Fertilizer") +
    ggplot2::labs(title=title, subtitle="Median ± IQR, pooled across all fields",
        x="Timepoint", y=metric) +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
        plot.title=ggplot2::element_text(face="bold"))
  tbl <- pvalue_transition_table(d, metric, min_n)
  grob <- pvalue_table_grob(tbl, paste0("Raw Wilcoxon p (pooled) — ", metric))
  combined <- p / grob + patchwork::plot_layout(heights = c(3, 1.4))
  ggplot2::ggsave(out_path, combined, width=8, height=8, dpi=200)
  out_path
}

# ============================================================================
# MAIN: depth-first tree across universes
# ============================================================================
build_report_tree <- function(alpha_df, comparisons, root = "Results") {
  universes <- list_universes(alpha_df)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)
    d_u   <- dplyr::filter(alpha_df, marker==mk, stage==st)
    if (nrow(d_u) == 0) next
    pal   <- fert_palette(d_u$fertilizer)
    message("=== Universe ", ucode, " (", nrow(d_u), " samples) ===")

    stats_u <- universe_alpha_stats(alpha_df, comparisons, mk, st)
    fields  <- sort(unique(d_u$field))

    ## Goal A — one combined-timepoint figure per field
    for (fld in fields) {
      leaf <- file.path(udir, GOAL_DIR["A"], fld)
      dir.create(leaf, recursive=TRUE, showWarnings=FALSE)
      dsub <- dplyr::filter(d_u, field==fld)
      w <- plot_goalA_field(dsub, pal,
            title=paste0(ucode," — Goal A — ",fld),
            out_path=file.path(leaf, paste0("alpha_box_",fld,".png")))
      if (!is.na(w)) written <- c(written, w)
      sl <- dplyr::filter(stats_u$global, goal=="A",
                          grepl(paste0("field=",fld,"(;|$)"), context))
      readr::write_csv(sl, file.path(leaf, paste0("stats_A_",fld,".csv")))
    }

    ## Goal B — per field, 3 metric trajectories + p-table
    for (fld in fields) {
      leaf <- file.path(udir, GOAL_DIR["B"], fld)
      dir.create(leaf, recursive=TRUE, showWarnings=FALSE)
      dsub <- dplyr::filter(d_u, field==fld)
      for (met in c("Observed","Pielou","Shannon")) {
        w <- plot_trajectory_metric(dsub, met, pal,
              title=paste0(ucode," — Goal B — ",fld," — ",met),
              out_path=file.path(leaf, paste0("alpha_traj_",fld,"_",tolower(met),".png")),
              facet_field=FALSE)
        if (!is.na(w)) written <- c(written, w)
      }
      sl <- dplyr::filter(stats_u$global, goal %in% c("B","B_per_fertilizer"),
                          grepl(paste0("field=",fld,"(;| |$)"), context))
      readr::write_csv(sl, file.path(leaf, paste0("stats_B_",fld,".csv")))
    }

    ## Goal C — pooled cross-field snapshot (no field subfolders)
    leafC <- file.path(udir, GOAL_DIR["C"]); dir.create(leafC, recursive=TRUE, showWarnings=FALSE)
    w <- plot_goalC_pooled(d_u, pal,
          title=paste0(ucode," — Goal C — fertilizers across fields"),
          out_path=file.path(leafC, "alpha_box_crossfield.png"))
    if (!is.na(w)) written <- c(written, w)
    readr::write_csv(dplyr::filter(stats_u$global, goal %in% c("C","C_per_field")),
                     file.path(leafC, "stats_C.csv"))

    ## Goal D — faceted-by-field trajectories + master pooled plot
    leafD <- file.path(udir, GOAL_DIR["D"]); dir.create(leafD, recursive=TRUE, showWarnings=FALSE)
    for (met in c("Observed","Pielou","Shannon")) {
      w1 <- plot_trajectory_metric(d_u, met, pal,
            title=paste0(ucode," — Goal D (faceted) — ",met),
            out_path=file.path(leafD, paste0("alpha_traj_faceted_",tolower(met),".png")),
            facet_field=TRUE)
      if (!is.na(w1)) written <- c(written, w1)
      w2 <- plot_goalD_pooled(d_u, met, pal,
            title=paste0(ucode," — Goal D (pooled) — ",met),
            out_path=file.path(leafD, paste0("alpha_traj_pooled_",tolower(met),".png")))
      if (!is.na(w2)) written <- c(written, w2)
    }
    readr::write_csv(dplyr::filter(stats_u$global, goal %in% c("D","D_per_fertilizer")),
                     file.path(leafD, "stats_D.csv"))
  }
  message("Report tree complete: ", length(written), " plots under ", root, "/")
  written
}

# per-universe stats wrapper (uses the alpha engine on the stage subset)
universe_alpha_stats <- function(alpha_df, comparisons, marker, stage) {
  mk <- marker; st <- stage
  d <- dplyr::filter(alpha_df, marker == mk, stage == st)
  if (nrow(d)==0) return(list(global=tibble::tibble(), posthoc=tibble::tibble()))
  run_alpha_stats(d, comparisons, mk)
}
