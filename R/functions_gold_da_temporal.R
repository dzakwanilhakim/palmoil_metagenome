# =============================================================================
# R/functions_gold_da_temporal.R — consecutive temporal ANCOM-BC2 per kebun
#
# Unit of analysis:
#   cohort (marker x stage) x field (Kode Kebun) x rank (genus/species)
#
# All fertilizers are pooled within a field. One model is fitted across the
# usable timepoints in that field, then consecutive contrasts are exported
# (T0 vs T1, T1 vs T2, ...) plus one initial-vs-final contrast (e.g. T0 vs T3)
# between the first and last usable timepoints. Each contrast gets an
# all-taxa CSV, a significant-taxa CSV, and a dumbbell plot.
#
# Trajectory: when a field has >= 3 usable timepoints, the same model's
# LFCs vs the first timepoint are drawn as one panel (T0 -> final T) for the
# top-N taxa, ranked by ANCOM-BC2's global test (any change across time),
# then by significance vs the first timepoint, then by max |LFC|.
#
# The post-gold-QC cohort gate in _targets.R determines which cohorts reach
# this driver; `da.fields` in analysis_thresholds.yaml optionally limits kebun.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
})

.da_cfg <- function(x, default) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) default else x[[1]]
}

.ordered_timepoints <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  number <- suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+$")))
  x[order(number, x, na.last = TRUE)]
}

# Only nominally consecutive timepoints are compared. If T1 is absent, this
# deliberately does not reinterpret T0 vs T2 as a consecutive comparison.
consecutive_time_pairs <- function(timepoints) {
  tp <- .ordered_timepoints(timepoints)
  if (length(tp) < 2)
    return(tibble::tibble(earlier = character(), later = character()))

  number <- suppressWarnings(as.integer(stringr::str_extract(tp, "[0-9]+$")))
  keep <- !is.na(number[-length(number)]) & !is.na(number[-1]) &
    number[-1] == number[-length(number)] + 1L
  tibble::tibble(earlier = tp[-length(tp)][keep], later = tp[-1][keep])
}

.temporal_status_result <- function(cohort, marker, stage, field, rank,
                                    earlier, later, n_earlier, n_later,
                                    status, note) {
  tibble::tibble(
    cohort = cohort, marker = marker, stage = stage, field = field, rank = rank,
    earlier = earlier, later = later,
    contrast = paste0(later, " vs ", earlier),
    n_earlier = as.integer(n_earlier), n_later = as.integer(n_later),
    taxon = NA_character_, lfc = NA_real_, se = NA_real_, W = NA_real_,
    p_val = NA_real_, q_val = NA_real_, diff_abn = NA,
    passed_ss = NA, significant = FALSE, status = status, note = note)
}

.result_column <- function(res, prefix, suffix, default = NA) {
  nm <- paste0(prefix, suffix)
  if (nm %in% names(res)) res[[nm]] else rep(default, nrow(res))
}

# Extract one directional contrast. Positive LFC means higher abundance at the
# later timepoint. With T0 as the model reference, T0->T1 is `waktuT1`; later
# pairs in ANCOM-BC2's pairwise table are `waktuT2_waktuT1`, etc.
.extract_temporal_contrast <- function(out, cohort, marker, stage, field, rank,
                                       earlier, later, n_earlier, n_later,
                                       cutoff, pseudo_sensitivity) {
  earlier_term <- paste0("waktu", make.names(earlier))
  later_term <- paste0("waktu", make.names(later))

  pair_res <- out$res_pair
  if (!is.null(pair_res)) {
    suffix <- paste0(later_term, "_", earlier_term)
    # Comparisons against the reference timepoint retain only the later term.
    if (!paste0("lfc_", suffix) %in% names(pair_res) &&
        paste0("lfc_", later_term) %in% names(pair_res))
      suffix <- later_term
    res <- pair_res
  } else {
    suffix <- later_term
    res <- out$res
  }

  if (is.null(res) || !paste0("lfc_", suffix) %in% names(res)) {
    available <- if (is.null(res)) "none" else
      paste(grep("^lfc_", names(res), value = TRUE), collapse = ", ")
    return(.temporal_status_result(
      cohort, marker, stage, field, rank, earlier, later, n_earlier, n_later,
      "parse_failed", paste0("Expected ANCOM-BC2 contrast '", suffix,
                              "'; available LFC columns: ", available)))
  }

  diff_abn <- as.logical(.result_column(res, "diff_", suffix, FALSE))
  passed_ss <- as.logical(.result_column(res, "passed_ss_", suffix, NA))
  robust_col <- paste0("diff_robust_", suffix)
  significant <- if (isTRUE(pseudo_sensitivity) && robust_col %in% names(res)) {
    as.logical(res[[robust_col]])
  } else if (isTRUE(pseudo_sensitivity) &&
             paste0("passed_ss_", suffix) %in% names(res)) {
    diff_abn & passed_ss
  } else {
    diff_abn
  }
  significant[is.na(significant)] <- FALSE

  tibble::tibble(
    cohort = cohort, marker = marker, stage = stage, field = field, rank = rank,
    earlier = earlier, later = later,
    contrast = paste0(later, " vs ", earlier),
    n_earlier = as.integer(n_earlier), n_later = as.integer(n_later),
    taxon = as.character(res$taxon),
    lfc = as.numeric(.result_column(res, "lfc_", suffix, NA_real_)),
    se = as.numeric(.result_column(res, "se_", suffix, NA_real_)),
    W = as.numeric(.result_column(res, "W_", suffix, NA_real_)),
    p_val = as.numeric(.result_column(res, "p_", suffix, NA_real_)),
    q_val = as.numeric(.result_column(res, "q_", suffix, NA_real_)),
    diff_abn = diff_abn, passed_ss = passed_ss,
    significant = significant, status = "ok", note = NA_character_)
}

.processed_matrix_sample_columns <- function(processed_mat) {
  taxon_col <- names(processed_mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(processed_mat))
  setdiff(names(processed_mat), c(taxon_col, meta_cols))
}

.temporal_phyloseq <- function(processed_mat, meta, time_levels) {
  sample_cols <- .processed_matrix_sample_columns(processed_mat)
  keep <- intersect(sample_cols, meta$`Sample alias`)
  meta <- meta[match(keep, meta$`Sample alias`), , drop = FALSE]

  taxon_col <- names(processed_mat)[1]
  otu <- as.matrix(processed_mat[, keep, drop = FALSE])
  storage.mode(otu) <- "double"

  taxa <- as.character(processed_mat[[taxon_col]])
  missing_taxa <- is.na(taxa) | taxa == ""
  taxa[missing_taxa] <- paste0("Unknown_taxon_", which(missing_taxa))
  otu <- rowsum(otu, group = taxa, reorder = FALSE)
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]

  meta$waktu <- factor(as.character(meta$waktu), levels = time_levels)
  samp <- as.data.frame(meta)
  rownames(samp) <- keep

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(samp))
}

# Initial-vs-final pair: first and last timepoints in this field that meet
# min_samples (e.g. T0 vs T3). Not added when it equals a consecutive pair.
.initial_final_pair <- function(counts, min_samples, pairs) {
  times <- .ordered_timepoints(names(counts)[counts >= min_samples])
  if (length(times) < 2)
    return(tibble::tibble(earlier = character(), later = character(),
                          comparison = character()))
  tibble::tibble(earlier = times[[1]], later = times[[length(times)]],
                 comparison = "initial_vs_final") |>
    dplyr::anti_join(pairs, by = c("earlier", "later"))
}

# Trajectory table: LFC of every model timepoint vs the model reference (the
# first usable timepoint, LFC = 0), with significance vs the reference and the
# global test across all timepoints. `trajectory_rank` marks the top-N taxa.
.extract_temporal_trajectory <- function(out, cohort, marker, stage, field, rank,
                                         model_times, counts, pseudo_sensitivity,
                                         max_taxa) {
  res <- out$res
  ref <- model_times[[1]]
  long <- purrr::map_dfr(model_times[-1], function(tp) {
    term <- paste0("waktu", make.names(tp))
    if (is.null(res) || !paste0("lfc_", term) %in% names(res)) return(NULL)
    diff_abn <- as.logical(.result_column(res, "diff_", term, FALSE))
    robust_col <- paste0("diff_robust_", term)
    sig <- if (isTRUE(pseudo_sensitivity) && robust_col %in% names(res))
      as.logical(res[[robust_col]]) else diff_abn
    tibble::tibble(
      taxon = as.character(res$taxon), waktu = tp,
      lfc = as.numeric(.result_column(res, "lfc_", term, NA_real_)),
      se = as.numeric(.result_column(res, "se_", term, NA_real_)),
      q_val = as.numeric(.result_column(res, "q_", term, NA_real_)),
      significant_vs_ref = tidyr::replace_na(sig, FALSE))
  })
  if (nrow(long) == 0) return(NULL)

  glob <- out$res_global
  glob <- if (is.null(glob)) {
    tibble::tibble(taxon = unique(long$taxon), global_q = NA_real_,
                   global_significant = NA)
  } else {
    gsig <- if (isTRUE(pseudo_sensitivity) && "diff_robust_abn" %in% names(glob))
      glob$diff_robust_abn else glob$diff_abn
    tibble::tibble(taxon = as.character(glob$taxon),
                   global_q = as.numeric(glob$q_val),
                   global_significant = tidyr::replace_na(as.logical(gsig), FALSE))
  }

  ranking <- long |>
    dplyr::group_by(taxon) |>
    dplyr::summarise(any_sig_vs_ref = any(significant_vs_ref),
                     max_abs_lfc = suppressWarnings(max(abs(lfc), na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::left_join(glob, by = "taxon") |>
    dplyr::mutate(selected_by = dplyr::case_when(
      global_significant %in% TRUE ~ "global test",
      any_sig_vs_ref ~ paste0("vs ", ref),
      TRUE ~ NA_character_)) |>
    dplyr::filter(!is.na(selected_by)) |>
    dplyr::arrange(selected_by != "global test", dplyr::desc(max_abs_lfc)) |>
    dplyr::slice_head(n = max_taxa) |>
    dplyr::mutate(trajectory_rank = dplyr::row_number()) |>
    dplyr::select(taxon, selected_by, trajectory_rank)

  base <- tibble::tibble(taxon = unique(long$taxon), waktu = ref, lfc = 0,
                         se = 0, q_val = NA_real_, significant_vs_ref = FALSE)
  dplyr::bind_rows(base, long) |>
    dplyr::left_join(glob, by = "taxon") |>
    dplyr::left_join(ranking, by = "taxon") |>
    dplyr::mutate(cohort = cohort, marker = marker, stage = stage,
                  field = field, rank = rank, reference = ref,
                  final = model_times[[length(model_times)]],
                  time_order = match(waktu, model_times),
                  n = as.integer(counts[waktu]), .before = 1) |>
    dplyr::arrange(trajectory_rank, taxon, time_order)
}

.run_temporal_field_rank <- function(processed_mat, lookup, marker, stage,
                                     field, field_members, rank, pairs, cutoff,
                                     min_samples, pseudo_sensitivity, n_cl,
                                     initial_vs_final = TRUE, global_test = TRUE,
                                     max_taxa_trajectory = 20L) {
  cohort <- paste0(marker, "_", stage)
  sample_cols <- .processed_matrix_sample_columns(processed_mat)
  meta <- lookup |>
    dplyr::filter(stage == !!stage, .data$field %in% field_members,
                  `Sample alias` %in% sample_cols) |>
    dplyr::mutate(field = !!field)
  counts <- table(as.character(meta$waktu))
  pairs <- dplyr::mutate(pairs, comparison = "consecutive")
  if (isTRUE(initial_vs_final))
    pairs <- dplyr::bind_rows(pairs, .initial_final_pair(counts, min_samples, pairs))
  pair_counts <- pairs |>
    dplyr::mutate(
      n_earlier = as.integer(counts[earlier]),
      n_later = as.integer(counts[later]),
      n_earlier = tidyr::replace_na(n_earlier, 0L),
      n_later = tidyr::replace_na(n_later, 0L),
      usable = n_earlier >= min_samples & n_later >= min_samples)
  tag <- function(x) dplyr::left_join(
    x, dplyr::select(pairs, earlier, later, comparison), by = c("earlier", "later"))
  done <- function(x, trajectory = NULL) list(result = tag(x), trajectory = trajectory)

  skipped <- purrr::pmap_dfr(
    dplyr::filter(pair_counts, !usable),
    function(earlier, later, n_earlier, n_later, ...) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "skipped_insufficient_replication",
        paste0("Requires >=", min_samples, " samples at each timepoint"))
    })

  usable_pairs <- dplyr::filter(pair_counts, usable)
  if (nrow(usable_pairs) == 0)
    return(done(skipped))

  model_times <- .ordered_timepoints(c(usable_pairs$earlier, usable_pairs$later))
  model_meta <- dplyr::filter(meta, waktu %in% model_times)
  ps <- .temporal_phyloseq(processed_mat, model_meta, model_times)

  if (phyloseq::ntaxa(ps) < 2) {
    failed <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                     n_later, ...) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "skipped_insufficient_taxa", "<2 non-zero taxa")
    })
    return(done(dplyr::bind_rows(failed, skipped)))
  }

  multi <- length(model_times) > 2
  context <- paste(cohort, field, rank, paste(model_times, collapse = "-"), sep = "|")
  message("ANCOM-BC2 temporal [", context, "]: samples=", phyloseq::nsamples(ps),
          ", taxa=", phyloseq::ntaxa(ps), ", fertilizers pooled",
          if (multi && isTRUE(global_test)) ", global test on" else "")

  out <- tryCatch(
    ANCOMBC::ancombc2(
      data = ps, fix_formula = "waktu", group = "waktu",
      p_adj_method = "BH", pseudo_sens = pseudo_sensitivity,
      prv_cut = 0, lib_cut = 0, struc_zero = FALSE, neg_lb = FALSE,
      alpha = cutoff, n_cl = n_cl, verbose = FALSE,
      global = multi && isTRUE(global_test), pairwise = multi,
      dunnet = FALSE, trend = FALSE),
    error = function(e) e)

  if (inherits(out, "error")) {
    failed <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                     n_later, ...) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "ancombc2_failed", conditionMessage(out))
    })
    return(done(dplyr::bind_rows(failed, skipped)))
  }

  tested <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                   n_later, ...) {
    .extract_temporal_contrast(
      out, cohort, marker, stage, field, rank, earlier, later,
      n_earlier, n_later, cutoff, pseudo_sensitivity)
  })
  trajectory <- if (multi)
    .extract_temporal_trajectory(out, cohort, marker, stage, field, rank,
                                 model_times, counts, pseudo_sensitivity,
                                 max_taxa_trajectory)
  done(dplyr::bind_rows(tested, skipped), trajectory)
}

.safe_da_path <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

.da_placeholder_plot <- function(reason, title, subtitle, style) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = stringr::str_wrap(reason, width = 90),
                      size = 5, fontface = "bold", colour = "grey35") +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    gold_plot_theme(style) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}

plot_temporal_da_dumbbell <- function(result, cutoff, max_taxa,
                                      style = load_plot_style(), out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  first <- result[1, , drop = FALSE]
  sig <- result |>
    dplyr::filter(status == "ok", significant, !is.na(taxon), is.finite(lfc)) |>
    dplyr::arrange(dplyr::desc(abs(lfc))) |>
    dplyr::slice_head(n = max_taxa)

  title <- paste0(first$cohort, " — ", first$field, " — ",
                  stringr::str_to_title(first$rank), " — ", first$contrast)
  subtitle <- paste0("ANCOM-BC2 temporal DA; all fertilizers pooled within kebun; ",
                     "BH q < ", cutoff)

  if (nrow(sig) == 0) {
    reason <- if (any(result$status == "ok"))
      "No significant taxa" else paste(unique(stats::na.omit(result$note)), collapse = "; ")
    if (!nzchar(reason)) reason <- "Analysis did not produce significant taxa"
    p <- .da_placeholder_plot(reason, title, subtitle, style)
  } else {
    sig <- dplyr::arrange(sig, lfc)
    sig$taxon <- factor(sig$taxon, levels = unique(sig$taxon))
    sig$direction <- ifelse(sig$lfc >= 0,
                            paste0("Higher in ", first$later),
                            paste0("Higher in ", first$earlier))
    pal <- stats::setNames(c("#E74C3C", "#2980B9"),
                           c(paste0("Higher in ", first$later),
                             paste0("Higher in ", first$earlier)))

    p <- ggplot2::ggplot(sig, ggplot2::aes(y = taxon)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = lfc, yend = taxon,
                                         colour = direction), linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(x = 0), colour = "grey45", size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = lfc, colour = direction), size = 3) +
      ggplot2::scale_colour_manual(values = pal, name = NULL) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = paste0("Log fold change (", first$later,
                               " relative to ", first$earlier, ")"), y = NULL) +
      gold_plot_theme(style) + ggplot2::theme(legend.position = "bottom")
  }

  ggplot2::ggsave(out_path, p, width = 10,
                  height = max(5, min(16, 0.34 * max(1, nrow(sig)) + 3)),
                  dpi = 200, limitsize = FALSE)
  out_path
}

.write_temporal_pair <- function(result, cutoff, max_taxa, style, root) {
  first <- result[1, , drop = FALSE]
  leaf <- file.path(root, first$cohort, "Differential_Abundance_Temporal",
                    .safe_da_path(first$field),
                    paste0(.safe_da_path(first$earlier), "_vs_",
                           .safe_da_path(first$later)), first$rank)
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

  all_csv <- file.path(leaf, "ancombc_all_taxa.csv")
  sig_csv <- file.path(leaf, "ancombc_significant.csv")
  plot_path <- file.path(leaf, "dumbbell.png")
  readr::write_csv(result, all_csv)
  readr::write_csv(dplyr::filter(result, status == "ok", significant), sig_csv)
  plot_temporal_da_dumbbell(result, cutoff, max_taxa, style, plot_path)

  summary <- first |>
    dplyr::select(cohort, marker, stage, field, rank, comparison, earlier, later,
                  contrast, n_earlier, n_later) |>
    dplyr::mutate(
      status = if (any(result$status == "ok")) "ok" else result$status[[1]],
      note = paste(unique(stats::na.omit(result$note)), collapse = "; "),
      taxa_tested = sum(result$status == "ok" & !is.na(result$taxon)),
      significant_taxa = sum(result$status == "ok" & result$significant, na.rm = TRUE),
      all_taxa_csv = all_csv, significant_csv = sig_csv, dumbbell_plot = plot_path)
  list(files = c(all_csv, sig_csv, plot_path), summary = summary)
}

# One panel, T0 -> final T: LFC vs the first timepoint for the top-N taxa.
# Filled points = significant vs the first timepoint at that timepoint.
plot_temporal_da_trajectory <- function(trajectory, cutoff, max_taxa,
                                        style = load_plot_style(), out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  first <- trajectory[1, , drop = FALSE]
  times <- dplyr::distinct(trajectory, time_order, waktu, n) |>
    dplyr::arrange(time_order)
  title <- paste0(first$cohort, " — ", first$field, " — ",
                  stringr::str_to_title(first$rank), " — trajectory ",
                  first$reference, " to ", first$final)
  subtitle <- paste0("ANCOM-BC2, all fertilizers pooled within kebun; top ", max_taxa,
                     " taxa by global test, then vs ", first$reference,
                     ", then max |LFC|; BH q < ", cutoff) |>
    stringr::str_wrap(width = 95)
  d <- dplyr::filter(trajectory, !is.na(trajectory_rank), is.finite(lfc))

  if (nrow(d) == 0) {
    p <- .da_placeholder_plot("No taxa significant across time", title, subtitle, style)
    height <- 5
  } else {
    sig_lab <- paste0("q < ", cutoff, " vs ", first$reference)
    d <- d |>
      dplyr::mutate(
        waktu = factor(waktu, levels = times$waktu),
        taxon = factor(taxon, levels = unique(taxon[order(trajectory_rank)])),
        point = dplyr::case_when(time_order == 1 ~ "Reference",
                                 significant_vs_ref ~ sig_lab,
                                 TRUE ~ "Not significant"))
    ends <- dplyr::filter(d, time_order == max(time_order))
    n_tp <- nrow(times)
    pal <- stats::setNames(grDevices::hcl.colors(nlevels(d$taxon), "Dark 3"),
                           levels(d$taxon))

    p <- ggplot2::ggplot(d, ggplot2::aes(waktu, lfc, group = taxon, colour = taxon)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
      ggplot2::geom_line(linewidth = 0.8, alpha = 0.85) +
      ggplot2::geom_point(ggplot2::aes(shape = point), size = 2.4,
                          fill = "white", stroke = 0.9) +
      ggrepel::geom_text_repel(
        data = ends, ggplot2::aes(label = taxon), size = 3, hjust = 0,
        direction = "y", nudge_x = 0.2, xlim = c(n_tp + 0.15, NA),
        min.segment.length = 0, segment.colour = "grey70", show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = pal, guide = "none") +
      ggplot2::scale_shape_manual(
        values = stats::setNames(c(18, 16, 21), c("Reference", sig_lab, "Not significant")),
        name = NULL) +
      ggplot2::scale_x_discrete(
        labels = paste0(times$waktu, "\n(n=", times$n, ")"),
        expand = ggplot2::expansion(add = c(0.3, 1.2))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title = title, subtitle = subtitle, x = NULL,
                    y = paste0("Log fold change vs ", first$reference)) +
      gold_plot_theme(style) + ggplot2::theme(legend.position = "bottom")
    height <- max(6, min(12, 0.25 * nlevels(d$taxon) + 4))
  }

  ggplot2::ggsave(out_path, p, width = 12, height = height, dpi = 200,
                  limitsize = FALSE)
  out_path
}

.write_temporal_trajectory <- function(trajectory, cutoff, max_taxa, style, root) {
  first <- trajectory[1, , drop = FALSE]
  leaf <- file.path(root, first$cohort, "Differential_Abundance_Temporal",
                    .safe_da_path(first$field),
                    paste0("Trajectory_", .safe_da_path(first$reference), "_to_",
                           .safe_da_path(first$final)), first$rank)
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

  all_csv <- file.path(leaf, "trajectory_all_taxa.csv")
  top_csv <- file.path(leaf, "trajectory_top_taxa.csv")
  plot_path <- file.path(leaf, "trajectory.png")
  readr::write_csv(trajectory, all_csv)
  readr::write_csv(dplyr::filter(trajectory, !is.na(trajectory_rank)), top_csv)
  plot_temporal_da_trajectory(trajectory, cutoff, max_taxa, style, plot_path)

  times <- dplyr::distinct(trajectory, time_order, waktu, n) |> dplyr::arrange(time_order)
  per_taxon <- dplyr::distinct(trajectory, taxon, global_significant, trajectory_rank)
  summary <- first |>
    dplyr::select(cohort, marker, stage, field, rank, reference, final) |>
    dplyr::mutate(
      timepoints = paste0(times$waktu, " (n=", times$n, ")", collapse = ", "),
      taxa_tested = nrow(per_taxon),
      global_significant_taxa = sum(per_taxon$global_significant %in% TRUE),
      taxa_plotted = sum(!is.na(per_taxon$trajectory_rank)),
      all_taxa_csv = all_csv, top_taxa_csv = top_csv, trajectory_plot = plot_path)
  list(files = c(all_csv, top_csv, plot_path), summary = summary)
}

build_gold_temporal_da <- function(gold_processed_matrix_16s_genus,
                                   gold_processed_matrix_16s_species,
                                   gold_universe_lookup_16s,
                                   gold_processed_matrix_its_genus,
                                   gold_processed_matrix_its_species,
                                   gold_universe_lookup_its,
                                   analysis_thresholds,
                                   style = load_plot_style(), root = "Results",
                                   combined_kebun = NULL) {
  da_cfg <- analysis_thresholds$da
  cutoff <- as.numeric(.da_cfg(da_cfg$adj_pval_cutoff, 0.05))
  min_samples <- as.integer(.da_cfg(da_cfg$min_samples_per_timepoint, 2L))
  pseudo_sensitivity <- as.logical(.da_cfg(da_cfg$pseudo_sensitivity, TRUE))
  n_cl <- as.integer(.da_cfg(da_cfg$n_cl, 1L))
  max_taxa <- as.integer(.da_cfg(da_cfg$max_taxa_dumbbell, 30L))
  initial_vs_final <- as.logical(.da_cfg(da_cfg$initial_vs_final, TRUE))
  global_test <- as.logical(.da_cfg(da_cfg$trajectory_global_test, TRUE))
  max_taxa_traj <- as.integer(.da_cfg(da_cfg$max_taxa_trajectory, 20L))
  only_fields <- as.character(unlist(da_cfg$fields %||% character()))

  if (!is.finite(cutoff) || cutoff <= 0 || cutoff >= 1)
    stop("Temporal DA: adj_pval_cutoff must be between 0 and 1.", call. = FALSE)
  if (is.na(min_samples) || min_samples < 2)
    stop("Temporal DA: min_samples_per_timepoint must be >= 2.", call. = FALSE)
  if (is.na(pseudo_sensitivity))
    stop("Temporal DA: pseudo_sensitivity must be true or false.", call. = FALSE)
  if (is.na(n_cl) || n_cl < 1)
    stop("Temporal DA: n_cl must be >= 1.", call. = FALSE)
  if (is.na(max_taxa) || max_taxa < 1)
    stop("Temporal DA: max_taxa_dumbbell must be >= 1.", call. = FALSE)
  if (is.na(initial_vs_final) || is.na(global_test))
    stop("Temporal DA: initial_vs_final and trajectory_global_test must be true or false.",
         call. = FALSE)
  if (is.na(max_taxa_traj) || max_taxa_traj < 1)
    stop("Temporal DA: max_taxa_trajectory must be >= 1.", call. = FALSE)

  message("Temporal DA settings: BH cutoff=", cutoff,
          ", min samples/timepoint=", min_samples,
          ", pseudo sensitivity=", pseudo_sensitivity,
          ", ANCOM-BC2 cores=", n_cl,
          ", initial vs final=", initial_vs_final,
          ", trajectory top N=", max_taxa_traj,
          " (global test=", global_test, ")",
          if (length(only_fields)) paste0(", fields=", paste(only_fields, collapse = ","))
          else "")

  specs <- list(
    list(marker = "16S", rank = "genus", matrix = gold_processed_matrix_16s_genus,
         lookup = gold_universe_lookup_16s),
    list(marker = "16S", rank = "species", matrix = gold_processed_matrix_16s_species,
         lookup = gold_universe_lookup_16s),
    list(marker = "ITS", rank = "genus", matrix = gold_processed_matrix_its_genus,
         lookup = gold_universe_lookup_its),
    list(marker = "ITS", rank = "species", matrix = gold_processed_matrix_its_species,
         lookup = gold_universe_lookup_its))

  written <- character(); summaries <- list(); idx <- 0L; traj_summaries <- list()
  for (spec in specs) {
    stages <- sort(unique(spec$lookup$stage))
    for (st in stages) {
      meta <- dplyr::filter(spec$lookup, stage == st)
      pairs <- consecutive_time_pairs(meta$waktu)
      if (nrow(pairs) == 0) next
      definitions <- analysis_field_definitions(
        spec$lookup, spec$marker, st, combined_kebun)
      for (field_i in seq_len(nrow(definitions))) {
        fld <- definitions$field[[field_i]]
        if (length(only_fields) && !fld %in% only_fields) next
        members <- as.character(unlist(
          definitions$members[[field_i]], use.names = FALSE))
        run <- .run_temporal_field_rank(
          spec$matrix, spec$lookup, spec$marker, st, fld, members,
          spec$rank, pairs,
          cutoff, min_samples, pseudo_sensitivity, n_cl,
          initial_vs_final, global_test, max_taxa_traj)
        field_pairs <- dplyr::distinct(run$result, comparison, earlier, later) |>
          dplyr::arrange(comparison != "consecutive", earlier)
        for (pair_i in seq_len(nrow(field_pairs))) {
          pair_result <- dplyr::filter(
            run$result, earlier == field_pairs$earlier[[pair_i]],
            later == field_pairs$later[[pair_i]])
          output <- .write_temporal_pair(pair_result, cutoff, max_taxa, style, root)
          idx <- idx + 1L
          written <- c(written, output$files)
          summaries[[idx]] <- output$summary
        }
        if (!is.null(run$trajectory) && nrow(run$trajectory) > 0) {
          traj <- .write_temporal_trajectory(run$trajectory, cutoff, max_taxa_traj,
                                             style, root)
          written <- c(written, traj$files)
          traj_summaries[[length(traj_summaries) + 1L]] <- traj$summary
        }
      }
    }
  }

  summary_tbl <- dplyr::bind_rows(summaries)
  summary_path <- file.path(root, "analysis", "da_temporal_summary.csv")
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary_tbl, summary_path)
  traj_tbl <- dplyr::bind_rows(traj_summaries)
  traj_path <- file.path(root, "analysis", "da_trajectory_summary.csv")
  readr::write_csv(traj_tbl, traj_path)
  written <- c(written, summary_path, traj_path)
  message("Temporal DA complete: ", nrow(summary_tbl), " field/rank/time comparisons, ",
          nrow(traj_tbl), " trajectories.")
  list(files = unique(written), summary = summary_tbl, summary_csv = summary_path,
       trajectory_summary = traj_tbl, trajectory_csv = traj_path)
}
