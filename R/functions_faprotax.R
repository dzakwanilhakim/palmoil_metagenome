# =============================================================================
# R/functions_faprotax.R — FAPROTAX annotation + temporal functional analysis
#
# FAPROTAX is applicable to prokaryotic 16S profiles, not ITS fungi. Species
# counts are converted to per-sample functional scores by the bundled
# collapse_table.py script. Functional groups can overlap, so their scores do
# not form a mutually exclusive composition. Consecutive timepoints, plus the
# first vs last usable timepoint (e.g. T0 vs T3), are therefore tested with an
# unpaired Wilcoxon rank-sum test (BH correction), separately for each cohort
# and kebun. Kebun with >= 3 usable timepoints also get a one-panel trajectory
# (top-N functions, Kruskal-Wallis across time). Fertilizers are pooled
# analytically: individual samples are retained, but fertilizer is not a
# grouping variable.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

.fap_cfg <- function(x, default) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) default else x[[1]]
}

.fap_safe_path <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

.fap_sample_columns <- function(processed_mat) {
  taxon_col <- names(processed_mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(processed_mat))
  setdiff(names(processed_mat), c(taxon_col, meta_cols))
}

.fap_empty_scores <- function() {
  tibble::tibble(function_name = character())
}

.fap_write_status <- function(path, status, note, n_taxa = 0L,
                              n_samples = 0L, n_functions = 0L) {
  out <- tibble::tibble(
    status = status,
    note = note,
    n_input_taxa = as.integer(n_taxa),
    n_samples = as.integer(n_samples),
    n_functions = as.integer(n_functions),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  readr::write_csv(out, path)
  path
}

# Run the locally supplied FAPROTAX program on the cohort-filtered 16S species
# matrix. Normalization occurs before collapsing, yielding sample-level
# relative functional scores while preserving overlapping function groups.
run_faprotax_annotation <- function(
    processed_16s_species,
    collapse_script = "FAPROTAX_1.2.12/collapse_table.py",
    database = "FAPROTAX_1.2.12/FAPROTAX.txt",
    python_bin = "python3",
    out_dir = "data/gold/faprotax") {
  if (!file.exists(collapse_script))
    stop("FAPROTAX collapse script not found: ", collapse_script, call. = FALSE)
  if (!file.exists(database))
    stop("FAPROTAX database not found: ", database, call. = FALSE)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  input_path <- file.path(out_dir, "faprotax_input_16s.tsv")
  scores_path <- file.path(out_dir, "faprotax_function_scores_16s.tsv")
  report_path <- file.path(out_dir, "faprotax_report.txt")
  mapping_path <- file.path(out_dir, "faprotax_taxon_function_mapping.tsv")
  log_path <- file.path(out_dir, "faprotax_run.log")
  status_path <- file.path(out_dir, "faprotax_status.csv")
  output_files <- c(input_path, scores_path, report_path, mapping_path,
                    log_path, status_path)

  sample_cols <- .fap_sample_columns(processed_16s_species)
  if (length(sample_cols) == 0 || nrow(processed_16s_species) == 0) {
    readr::write_tsv(tibble::tibble(taxonomy = character()), input_path)
    readr::write_tsv(.fap_empty_scores(), scores_path)
    writeLines("FAPROTAX skipped: no selected 16S samples.", report_path)
    writeLines("record\tgroup", mapping_path)
    writeLines("FAPROTAX skipped: no selected 16S samples.", log_path)
    .fap_write_status(status_path, "skipped_no_16s_samples",
                      "No selected 16S samples reached functional analysis.")
    return(list(scores = .fap_empty_scores(), files = output_files,
                status_csv = status_path, report = report_path,
                mapping = mapping_path))
  }

  taxonomy <- if ("tax" %in% names(processed_16s_species)) {
    as.character(processed_16s_species[["tax"]])
  } else {
    lineage <- setdiff(intersect(LINEAGE_COLS, names(processed_16s_species)),
                       c("total", "tax"))
    apply(processed_16s_species[, lineage, drop = FALSE], 1, function(x) {
      paste(stats::na.omit(as.character(x[x != ""])), collapse = ";")
    })
  }
  missing_taxonomy <- is.na(taxonomy) | !nzchar(taxonomy)
  taxonomy[missing_taxonomy] <- paste0(
    "Unclassified_taxon_", which(missing_taxonomy))

  abundance <- as.matrix(processed_16s_species[, sample_cols, drop = FALSE])
  storage.mode(abundance) <- "double"
  if (any(!is.finite(abundance)) || any(abundance < 0))
    stop("FAPROTAX input contains non-finite or negative counts.", call. = FALSE)
  abundance <- rowsum(abundance, group = taxonomy, reorder = FALSE)
  abundance <- abundance[rowSums(abundance) > 0, , drop = FALSE]

  input <- tibble::as_tibble(abundance, rownames = "taxonomy")
  readr::write_tsv(input, input_path)

  args <- c(
    collapse_script,
    "-i", input_path,
    "-g", database,
    "-o", scores_path,
    "-r", report_path,
    "-l", log_path,
    "--out_groups2records_table_dense", mapping_path,
    "-d", "taxonomy",
    "--column_names_are_in", "first_data_line",
    "--non_numeric", "ignore",
    "--normalize_collapsed", "columns_before_collapsing",
    "--group_leftovers_as", "unassigned",
    "--omit_unrepresented_groups",
    "-f", "-v")

  command_output <- tryCatch(
    system2(python_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 1L))
  exit_status <- attr(command_output, "status") %||% 0L
  # collapse_table.py also writes its own log; append captured Python output so
  # environment warnings and the exit context are retained in one audit file.
  cat(paste0("\n# Captured stdout/stderr\n",
             paste(command_output, collapse = "\n"), "\n"),
      file = log_path, append = TRUE)

  if (!identical(as.integer(exit_status), 0L) || !file.exists(scores_path)) {
    .fap_write_status(status_path, "failed", paste0(
      "collapse_table.py exited with status ", exit_status),
      n_taxa = nrow(abundance), n_samples = ncol(abundance))
    stop("FAPROTAX annotation failed (exit status ", exit_status,
         "). See ", log_path, call. = FALSE)
  }

  scores <- readr::read_tsv(scores_path, comment = "#", show_col_types = FALSE,
                            name_repair = "minimal")
  if (ncol(scores) == 0)
    scores <- .fap_empty_scores()
  if (ncol(scores) > 0)
    names(scores)[1] <- "function_name"
  score_cols <- intersect(sample_cols, names(scores))
  scores <- dplyr::select(scores, function_name, dplyr::all_of(score_cols))
  scores <- dplyr::mutate(scores,
                          dplyr::across(dplyr::all_of(score_cols), as.numeric))

  .fap_write_status(
    status_path, "ok",
    paste0("FAPROTAX annotation completed using ", basename(database), "."),
    n_taxa = nrow(abundance), n_samples = length(score_cols),
    n_functions = sum(scores$function_name != "unassigned"))
  message("FAPROTAX annotation: taxa=", nrow(abundance),
          ", samples=", length(score_cols),
          ", represented functions=",
          sum(scores$function_name != "unassigned"))

  list(scores = scores, files = output_files, status_csv = status_path,
       report = report_path, mapping = mapping_path)
}

.fap_status_result <- function(cohort, stage, field, earlier, later,
                               n_earlier, n_later, status, note) {
  tibble::tibble(
    cohort = cohort, marker = "16S", stage = stage, field = field,
    earlier = earlier, later = later, contrast = paste0(later, " vs ", earlier),
    n_earlier = as.integer(n_earlier), n_later = as.integer(n_later),
    function_name = NA_character_, prevalence = NA_real_,
    mean_earlier = NA_real_, mean_later = NA_real_, delta = NA_real_,
    log2_fc = NA_real_, p_val = NA_real_, q_val = NA_real_,
    significant = FALSE, status = status, note = note)
}

.fap_test_pair <- function(scores, meta, stage, field, earlier, later,
                           cutoff, min_samples, min_prevalence, pseudocount) {
  cohort <- paste0("16S_", stage)
  earlier_ids <- intersect(meta$`Sample alias`[meta$waktu == earlier], names(scores))
  later_ids <- intersect(meta$`Sample alias`[meta$waktu == later], names(scores))
  n_earlier <- length(earlier_ids)
  n_later <- length(later_ids)

  if (n_earlier < min_samples || n_later < min_samples) {
    return(.fap_status_result(
      cohort, stage, field, earlier, later, n_earlier, n_later,
      "skipped_insufficient_replication",
      paste0("Requires >=", min_samples, " samples at each timepoint.")))
  }

  functions <- which(!is.na(scores$function_name) &
                       scores$function_name != "unassigned")
  if (length(functions) == 0) {
    return(.fap_status_result(
      cohort, stage, field, earlier, later, n_earlier, n_later,
      "skipped_no_annotated_functions", "No represented FAPROTAX functions."))
  }

  tested <- purrr::map_dfr(functions, function(i) {
    x <- as.numeric(unlist(scores[i, earlier_ids, drop = FALSE],
                           use.names = FALSE))
    y <- as.numeric(unlist(scores[i, later_ids, drop = FALSE],
                           use.names = FALSE))
    prevalence <- mean(c(x, y) > 0, na.rm = TRUE)
    if (!is.finite(prevalence) || prevalence < min_prevalence)
      return(NULL)
    p <- if (length(unique(c(x, y))) < 2) {
      1
    } else {
      tryCatch(suppressWarnings(stats::wilcox.test(
        x, y, paired = FALSE, exact = FALSE)$p.value), error = function(e) NA_real_)
    }
    mean_x <- mean(x, na.rm = TRUE)
    mean_y <- mean(y, na.rm = TRUE)
    tibble::tibble(
      cohort = cohort, marker = "16S", stage = stage, field = field,
      earlier = earlier, later = later, contrast = paste0(later, " vs ", earlier),
      n_earlier = n_earlier, n_later = n_later,
      function_name = as.character(scores$function_name[[i]]),
      prevalence = prevalence,
      mean_earlier = mean_x, mean_later = mean_y,
      delta = mean_y - mean_x,
      log2_fc = log2((mean_y + pseudocount) / (mean_x + pseudocount)),
      p_val = p)
  })

  if (nrow(tested) == 0) {
    return(.fap_status_result(
      cohort, stage, field, earlier, later, n_earlier, n_later,
      "skipped_no_prevalent_functions",
      paste0("No functions met prevalence >=", min_prevalence, ".")))
  }

  tested |>
    dplyr::mutate(
      q_val = stats::p.adjust(p_val, method = "BH"),
      significant = !is.na(q_val) & q_val < cutoff,
      status = "ok", note = NA_character_) |>
    dplyr::arrange(q_val, dplyr::desc(abs(delta)), function_name)
}

# Rows a dumbbell shows: significant, not excluded, top N by |mean change|.
.fap_dumbbell_rows <- function(result, max_functions) {
  result |>
    dplyr::filter(status == "ok", significant, !is.na(function_name),
                  is.finite(mean_earlier), is.finite(mean_later),
                  !(excluded_from_top %in% TRUE)) |>
    dplyr::arrange(dplyr::desc(abs(delta)), q_val) |>
    dplyr::slice_head(n = max_functions)
}

# Upper x limit shared by every dumbbell of one kebun folder (0 to a rounded
# max of all plotted means), so the T0_vs_T1, T1_vs_T2, ... plots compare.
.fap_shared_x_max <- function(results, max_functions) {
  shown <- dplyr::bind_rows(lapply(results, .fap_dumbbell_rows, max_functions))
  if (nrow(shown) == 0) return(NULL)
  max(pretty(c(0, shown$mean_earlier, shown$mean_later)))
}

# x axis of the dumbbells. "sqrt" keeps small (0.1%) and large (40%) scores
# readable on one shared axis; "linear" is the plain percentage axis.
.fap_dumbbell_x_scale <- function(x_max, x_transform) {
  pct <- function(x) paste0(signif(100 * x, 2), "%")
  if (x_transform == "sqrt") {
    steps <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1)
    top <- x_max %||% 1
    return(ggplot2::scale_x_continuous(
      transform = "sqrt", labels = pct, limits = c(0, x_max),
      breaks = steps[steps <= top]))
  }
  ggplot2::scale_x_continuous(
    labels = pct, limits = if (!is.null(x_max)) c(0, x_max),
    breaks = if (!is.null(x_max)) pretty(c(0, x_max)) else ggplot2::waiver())
}

plot_faprotax_dumbbell <- function(result, cutoff, max_functions,
                                   style = load_plot_style(), out_path,
                                   x_max = NULL, x_transform = "sqrt") {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  first <- result[1, , drop = FALSE]
  sig <- .fap_dumbbell_rows(result, max_functions)

  title <- paste0(first$cohort, " — ", first$field, " — ", first$contrast)
  subtitle <- stringr::str_wrap(paste0(
    "FAPROTAX functional scores; fertilizers pooled within kebun; BH q < ",
    cutoff, "; top ", max_functions, " by absolute mean change",
    "; clinical/host-associated groups excluded",
    if (!is.null(x_max)) "; x axis shared within kebun" else "",
    if (x_transform == "sqrt") " (square-root scale)" else ""), width = 80)

  if (nrow(sig) == 0) {
    reason <- if (any(result$status == "ok")) {
      "No significant functions"
    } else {
      paste(unique(stats::na.omit(result$note)), collapse = "; ")
    }
    if (!nzchar(reason)) reason <- "Analysis did not produce significant functions"
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = stringr::str_wrap(reason, 90),
                        size = 5, fontface = "bold", colour = "grey35") +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
      gold_plot_theme(style) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
  } else {
    sig <- dplyr::arrange(sig, delta)
    sig$function_name <- factor(sig$function_name,
                                levels = unique(sig$function_name))
    long <- sig |>
      dplyr::select(function_name, earlier, later, mean_earlier, mean_later) |>
      tidyr::pivot_longer(c(mean_earlier, mean_later),
                          names_to = "endpoint", values_to = "mean_score") |>
      dplyr::mutate(timepoint = ifelse(endpoint == "mean_earlier", earlier, later))
    time_pal <- waktu_palette()
    missing_times <- setdiff(unique(long$timepoint), names(time_pal))
    if (length(missing_times) > 0) {
      time_pal <- c(time_pal, stats::setNames(
        grDevices::hcl.colors(length(missing_times), "Dark 3"), missing_times))
    }

    p <- ggplot2::ggplot(sig, ggplot2::aes(y = function_name)) +
      ggplot2::geom_segment(ggplot2::aes(x = mean_earlier, xend = mean_later,
                                         yend = function_name),
                            colour = "grey65", linewidth = 1.1) +
      ggplot2::geom_point(data = long,
                          ggplot2::aes(x = mean_score, y = function_name,
                                       colour = timepoint), size = 3.2) +
      ggplot2::scale_colour_manual(values = time_pal, name = "Timepoint") +
      .fap_dumbbell_x_scale(x_max, x_transform) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = "Mean normalized functional score", y = NULL) +
      gold_plot_theme(style) +
      ggplot2::theme(legend.position = "bottom", plot.title.position = "plot")
  }

  ggplot2::ggsave(out_path, p, width = 10,
                  height = max(5, min(13, 0.38 * max(1, nrow(sig)) + 3)),
                  dpi = 200, limitsize = FALSE)
  out_path
}

.write_fap_pair <- function(result, cutoff, max_functions, style, root,
                            x_max = NULL, x_transform = "sqrt") {
  first <- result[1, , drop = FALSE]
  leaf <- file.path(root, first$cohort, "Functional_FAPROTAX",
                    .fap_safe_path(first$field),
                    paste0(.fap_safe_path(first$earlier), "_vs_",
                           .fap_safe_path(first$later)))
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  all_csv <- file.path(leaf, "faprotax_all_functions.csv")
  sig_csv <- file.path(leaf, "faprotax_significant.csv")
  plot_path <- file.path(leaf, "dumbbell.png")
  readr::write_csv(result, all_csv)
  readr::write_csv(dplyr::filter(result, status == "ok", significant), sig_csv)
  plot_faprotax_dumbbell(result, cutoff, max_functions, style, plot_path,
                         x_max, x_transform)

  summary <- first |>
    dplyr::select(cohort, marker, stage, field, comparison, earlier, later,
                  contrast, n_earlier, n_later) |>
    dplyr::mutate(
      status = if (any(result$status == "ok")) "ok" else result$status[[1]],
      note = paste(unique(stats::na.omit(result$note)), collapse = "; "),
      functions_tested = sum(result$status == "ok" &
                               !is.na(result$function_name)),
      significant_functions = sum(result$status == "ok" & result$significant,
                                  na.rm = TRUE),
      all_functions_csv = all_csv, significant_csv = sig_csv,
      dumbbell_plot = plot_path)
  list(files = c(all_csv, sig_csv, plot_path), summary = summary)
}

.fap_coverage <- function(scores, lookup) {
  sample_cols <- intersect(names(scores), lookup$`Sample alias`)
  if (length(sample_cols) == 0)
    return(tibble::tibble())
  unassigned_i <- which(scores$function_name == "unassigned")
  unassigned <- if (length(unassigned_i) == 0) {
    stats::setNames(rep(0, length(sample_cols)), sample_cols)
  } else {
    stats::setNames(
      as.numeric(unlist(scores[unassigned_i[[1]], sample_cols, drop = FALSE],
                        use.names = FALSE)),
      sample_cols)
  }
  tibble::tibble(`Sample alias` = sample_cols,
                 assigned_fraction = pmax(0, pmin(1, 1 - unassigned))) |>
    dplyr::left_join(lookup, by = "Sample alias") |>
    dplyr::mutate(cohort = paste0("16S_", stage), .before = 1)
}

.fap_select_functions <- function(results, n) {
  eligible <- results |>
    dplyr::filter(status == "ok", !is.na(function_name), is.finite(delta),
                  !(excluded_from_top %in% TRUE))
  # dplyr may evaluate summary expressions once for type discovery even when
  # there are no groups. Return early so skipped kebun do not call max() on an
  # empty vector while their placeholder figures are being prepared.
  if (nrow(eligible) == 0)
    return(character())

  ranked <- eligible |>
    dplyr::group_by(function_name) |>
    dplyr::summarise(
      significant_any = any(significant, na.rm = TRUE),
      max_abs_change = {
        values <- abs(delta[is.finite(delta)])
        if (length(values) == 0) NA_real_ else max(values)
      },
      min_q = {
        values <- q_val[is.finite(q_val)]
        if (length(values) == 0) NA_real_ else min(values)
      },
      .groups = "drop") |>
    dplyr::arrange(dplyr::desc(significant_any),
                   dplyr::desc(max_abs_change), min_q, function_name)
  if (any(ranked$significant_any))
    ranked <- dplyr::filter(ranked, significant_any)
  utils::head(ranked$function_name, n)
}

.fap_placeholder_plot <- function(title, subtitle, note, style, out_path,
                                  width = 10, height = 6) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0,
                      label = stringr::str_wrap(note, 90), size = 5,
                      fontface = "bold", colour = "grey35") +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    gold_plot_theme(style) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(out_path, p, width = width, height = height, dpi = 200)
  out_path
}

# Primary overview figure: functions are ranked by their largest absolute
# change in mean score, never by a sum of q-values. The diverging scale shows
# later-minus-earlier percentage-point change; stars encode BH significance.
plot_faprotax_effect_heatmap <- function(results, max_functions, cutoff,
                                         style = load_plot_style(), out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  cohort <- unique(stats::na.omit(results$cohort))[[1]] %||% "16S"
  fields_in_plot <- unique(stats::na.omit(results$field))
  scope <- if (length(fields_in_plot) == 1)
    paste0(" — ", fields_in_plot[[1]]) else ""
  selected <- .fap_select_functions(results, max_functions)
  if (length(selected) == 0) {
    return(.fap_placeholder_plot(
      paste0(cohort, scope, " — FAPROTAX temporal effect overview"),
      "Later minus earlier mean functional score",
      "No testable FAPROTAX functions.", style, out_path))
  }

  pair_meta <- results |>
    dplyr::group_by(cohort, field, earlier, later) |>
    dplyr::summarise(
      pair_status = ifelse(any(status == "ok"), "ok", dplyr::first(status)),
      pair_type = dplyr::first(comparison),
      .groups = "drop") |>
    dplyr::arrange(field, pair_type != "consecutive", earlier) |>
    dplyr::mutate(
      pair_id = paste(field, earlier, later, sep = "__"),
      comparison = paste0(field, "\n", earlier, "→", later))

  values <- results |>
    dplyr::filter(function_name %in% selected) |>
    dplyr::mutate(pair_id = paste(field, earlier, later, sep = "__")) |>
    dplyr::select(pair_id, function_name, delta, q_val, significant, status)
  plot_data <- tidyr::expand_grid(
    function_name = selected, pair_id = pair_meta$pair_id) |>
    dplyr::left_join(values, by = c("function_name", "pair_id")) |>
    dplyr::left_join(dplyr::select(pair_meta, pair_id, comparison, pair_status),
                     by = "pair_id") |>
    dplyr::mutate(
      change_pp = 100 * delta,
      label = dplyr::case_when(
        pair_status != "ok" ~ "×",
        !is.na(significant) & significant ~ "*",
        TRUE ~ ""))

  pair_order <- pair_meta$comparison
  plot_data$comparison <- factor(plot_data$comparison, levels = pair_order)
  plot_data$function_name <- factor(
    plot_data$function_name, levels = rev(selected))
  max_change <- max(abs(plot_data$change_pp), na.rm = TRUE)
  if (!is.finite(max_change) || max_change == 0) max_change <- 1

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = comparison, y = function_name,
                                    fill = change_pp)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 5,
                       fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-max_change, max_change),
      name = "Mean change\n(percentage points)") +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_replace_all(x, "_", " ")) +
    ggplot2::labs(
      title = paste0(cohort, scope, " — FAPROTAX temporal effect overview"),
      subtitle = paste0("Top ", length(selected),
                        " functions by maximum absolute mean change; ",
                        "* BH q < ", cutoff, "; × insufficient data; ",
                        "clinical/host-associated groups excluded"),
      x = "Kebun and comparison (consecutive, then first to last)", y = NULL) +
    gold_plot_theme(style) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "right")

  ggplot2::ggsave(out_path, p,
                  width = max(11, 1.1 * nrow(pair_meta) + 5),
                  height = max(7, 0.34 * length(selected) + 3.5),
                  dpi = 200, limitsize = FALSE)
  out_path
}

plot_faprotax_trajectories <- function(scores, lookup, functions, valid_fields,
                                       style = load_plot_style(), out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  cohort <- paste0("16S_", unique(stats::na.omit(lookup$stage))[[1]] %||% "")
  fields_in_plot <- unique(stats::na.omit(lookup$field))
  scope <- if (length(fields_in_plot) == 1)
    paste0(" — ", fields_in_plot[[1]]) else ""
  sample_cols <- intersect(names(scores), lookup$`Sample alias`)
  if (length(functions) == 0 || length(sample_cols) == 0 ||
      length(valid_fields) == 0) {
    return(.fap_placeholder_plot(
      paste0(cohort, scope, " — selected FAPROTAX trajectories"),
      "Mean ± standard error; fertilizers pooled",
      "No functions or kebun had a testable temporal comparison.",
      style, out_path))
  }

  long <- scores |>
    dplyr::filter(function_name %in% functions) |>
    tidyr::pivot_longer(dplyr::all_of(sample_cols),
                        names_to = "Sample alias", values_to = "score") |>
    dplyr::left_join(lookup, by = "Sample alias") |>
    dplyr::filter(field %in% valid_fields, !is.na(waktu)) |>
    dplyr::mutate(score_pct = 100 * as.numeric(score),
                  waktu = factor(waktu, levels = .ordered_timepoints(waktu)),
                  function_name = factor(function_name, levels = functions))
  means <- long |>
    dplyr::group_by(function_name, field, waktu) |>
    dplyr::summarise(
      mean_score = mean(score_pct, na.rm = TRUE),
      se = ifelse(dplyr::n() > 1, stats::sd(score_pct, na.rm = TRUE) /
                    sqrt(dplyr::n()), NA_real_),
      .groups = "drop")

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = long, ggplot2::aes(x = waktu, y = score_pct),
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      colour = "grey65", alpha = 0.65, size = 1.3) +
    ggplot2::geom_line(
      data = means, ggplot2::aes(x = waktu, y = mean_score, group = 1),
      colour = "grey35", linewidth = 0.7) +
    ggplot2::geom_errorbar(
      data = means,
      ggplot2::aes(x = waktu, ymin = mean_score - se, ymax = mean_score + se),
      width = 0.12, colour = "grey35", na.rm = TRUE) +
    ggplot2::geom_point(
      data = means,
      ggplot2::aes(x = waktu, y = mean_score, colour = waktu), size = 2.6) +
    ggplot2::scale_colour_manual(values = waktu_palette(), name = "Timepoint") +
    ggplot2::facet_grid(function_name ~ field, scales = "free_y",
                        labeller = ggplot2::labeller(
                          function_name = function(x)
                            stringr::str_replace_all(x, "_", " "))) +
    ggplot2::labs(
      title = paste0(cohort, scope,
                     " — selected FAPROTAX temporal trajectories"),
      subtitle = "Points are individual samples; connected points are means ± SE; fertilizers pooled",
      x = "Timepoint", y = "Normalized functional score (%)") +
    gold_plot_theme(style) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(out_path, p,
                  width = max(11, 2.2 * length(valid_fields) + 4),
                  height = max(7, 1.65 * length(functions) + 3),
                  dpi = 200, limitsize = FALSE)
  out_path
}

# Trajectory table for one kebun: every prevalent function x usable timepoint
# (>= min_samples). value = log2 fold change of the mean score vs the first
# timepoint; Kruskal-Wallis across all timepoints = global test (BH across
# functions); Wilcoxon vs the first timepoint per timepoint (BH per timepoint).
# trajectory_rank marks the top N: global significant first, then significant
# vs the first timepoint, then max |log2 FC|; excluded groups never ranked.
.fap_trajectory <- function(scores, meta, stage, field, cutoff, min_samples,
                            min_prevalence, pseudocount, exclude, max_n) {
  meta <- dplyr::filter(meta, `Sample alias` %in% names(scores), !is.na(waktu))
  counts <- table(as.character(meta$waktu))
  times <- .ordered_timepoints(names(counts)[counts >= min_samples])
  if (length(times) < 3) return(NULL)
  meta <- dplyr::filter(meta, waktu %in% times)
  ids <- meta$`Sample alias`
  grp <- factor(as.character(meta$waktu), levels = times)
  ref_i <- grp == times[[1]]
  functions <- which(!is.na(scores$function_name) & scores$function_name != "unassigned")

  long <- purrr::map_dfr(functions, function(i) {
    v <- as.numeric(unlist(scores[i, ids, drop = FALSE], use.names = FALSE))
    if (mean(v > 0, na.rm = TRUE) < min_prevalence) return(NULL)
    global_p <- if (length(unique(v)) < 2) 1 else
      tryCatch(stats::kruskal.test(v, grp)$p.value, error = function(e) NA_real_)
    p_ref <- vapply(times[-1], function(t) {
      x <- v[ref_i]; y <- v[grp == t]
      if (length(unique(c(x, y))) < 2) return(1)
      tryCatch(suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
               error = function(e) NA_real_)
    }, numeric(1))
    means <- as.numeric(tapply(v, grp, mean))
    tibble::tibble(
      function_name = as.character(scores$function_name[[i]]),
      waktu = times, time_order = seq_along(times),
      n = as.integer(table(grp)), mean_score = means,
      value = log2((means + pseudocount) / (means[[1]] + pseudocount)),
      p_vs_ref = c(NA_real_, p_ref), global_p = global_p)
  })
  if (nrow(long) == 0) return(NULL)

  global <- dplyr::distinct(long, function_name, global_p) |>
    dplyr::mutate(global_q = stats::p.adjust(global_p, method = "BH"))
  long <- long |>
    dplyr::group_by(waktu) |>
    dplyr::mutate(q_vs_ref = stats::p.adjust(p_vs_ref, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::left_join(dplyr::select(global, function_name, global_q), by = "function_name") |>
    dplyr::mutate(significant_vs_ref = !is.na(q_vs_ref) & q_vs_ref < cutoff,
                  global_significant = !is.na(global_q) & global_q < cutoff,
                  excluded_from_top = function_name %in% exclude)

  ranking <- long |>
    dplyr::filter(!excluded_from_top) |>
    dplyr::group_by(function_name) |>
    dplyr::summarise(global_sig = any(global_significant),
                     any_sig = any(significant_vs_ref),
                     max_abs = max(abs(value), na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(global_sig | any_sig) |>
    dplyr::arrange(dplyr::desc(global_sig), dplyr::desc(max_abs)) |>
    dplyr::slice_head(n = max_n) |>
    dplyr::mutate(trajectory_rank = dplyr::row_number(),
                  selected_by = ifelse(global_sig, "global test", paste0("vs ", times[[1]])))

  long |>
    dplyr::left_join(dplyr::select(ranking, function_name, selected_by, trajectory_rank),
                     by = "function_name") |>
    dplyr::mutate(cohort = paste0("16S_", stage), stage = stage, field = field,
                  reference = times[[1]], final = times[[length(times)]], .before = 1) |>
    dplyr::rename(log2_fc_vs_ref = value) |>
    dplyr::arrange(trajectory_rank, function_name, time_order)
}

.write_fap_trajectory <- function(trajectory, cutoff, max_n, style, field_dir) {
  first <- trajectory[1, , drop = FALSE]
  csv <- file.path(field_dir, "faprotax_trajectory_panel.csv")
  png <- file.path(field_dir, "faprotax_trajectory_panel.png")
  readr::write_csv(trajectory, csv)
  plot_trajectory_panel(
    dplyr::rename(trajectory, label = function_name, value = log2_fc_vs_ref) |>
      dplyr::mutate(label = stringr::str_replace_all(label, "_", " ")),
    title = paste0(first$cohort, " — ", first$field, " — FAPROTAX trajectory ",
                   first$reference, " to ", first$final),
    subtitle = paste0("Mean functional score, fertilizers pooled within kebun; top ", max_n,
                      " functions by Kruskal-Wallis across time, then Wilcoxon vs ",
                      first$reference, ", then max |log2 FC|; BH q < ", cutoff,
                      "; clinical/host-associated groups excluded"),
    y_label = paste0("log2 fold change of mean score vs ", first$reference),
    reference = first$reference, cutoff = cutoff,
    empty_note = "No functions significant across time", style = style, out_path = png)

  per_fn <- dplyr::distinct(trajectory, function_name, global_significant, trajectory_rank)
  times <- dplyr::distinct(trajectory, time_order, waktu, n) |> dplyr::arrange(time_order)
  summary <- first |>
    dplyr::select(cohort, stage, field, reference, final) |>
    dplyr::mutate(
      timepoints = paste0(times$waktu, " (n=", times$n, ")", collapse = ", "),
      functions_tested = nrow(per_fn),
      global_significant = sum(per_fn$global_significant),
      functions_plotted = sum(!is.na(per_fn$trajectory_rank)),
      trajectory_csv = csv, trajectory_plot = png)
  list(files = c(csv, png), summary = summary)
}

plot_faprotax_coverage <- function(coverage, style = load_plot_style(), out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  if (nrow(coverage) == 0) {
    return(.fap_placeholder_plot(
      "FAPROTAX annotation coverage", NULL, "No selected 16S samples.",
      style, out_path))
  }
  coverage <- coverage |>
    dplyr::mutate(waktu = factor(waktu, levels = .ordered_timepoints(waktu)))
  cohorts_in_plot <- unique(stats::na.omit(coverage$cohort))
  fields_in_plot <- unique(stats::na.omit(coverage$field))
  scope <- paste0(
    if (length(cohorts_in_plot) == 1) paste0(cohorts_in_plot[[1]], " — ") else "",
    if (length(fields_in_plot) == 1) paste0(fields_in_plot[[1]], " — ") else "")
  p <- ggplot2::ggplot(
    coverage, ggplot2::aes(x = waktu, y = assigned_fraction, colour = waktu)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.15, width = 0.55) +
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.8, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = waktu_palette(), name = "Timepoint") +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
    ggplot2::facet_grid(stage ~ field, scales = "free_x", space = "free_x") +
    ggplot2::labs(
      title = paste0(scope, "FAPROTAX annotation coverage"),
      subtitle = "Fraction of each 16S sample assigned to at least one FAPROTAX function",
      x = "Timepoint", y = "Assigned fraction") +
    gold_plot_theme(style) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(out_path, p, width = 14,
                  height = max(6, 3.2 * length(unique(coverage$stage)) + 2),
                  dpi = 200, limitsize = FALSE)
  out_path
}

build_faprotax_temporal <- function(faprotax_annotation, lookup_16s,
                                    analysis_thresholds,
                                    style = load_plot_style(),
                                    root = "Results",
                                    combined_kebun = NULL) {
  cfg <- analysis_thresholds$faprotax
  cutoff <- as.numeric(.fap_cfg(cfg$adj_pval_cutoff, 0.05))
  min_samples <- as.integer(.fap_cfg(cfg$min_samples_per_timepoint, 2L))
  min_prevalence <- as.numeric(.fap_cfg(cfg$min_prevalence_frac, 0.10))
  pseudocount <- as.numeric(.fap_cfg(cfg$pseudocount, 1e-6))
  max_functions <- as.integer(.fap_cfg(cfg$max_functions_dumbbell, 20L))
  max_heatmap <- as.integer(.fap_cfg(cfg$max_functions_heatmap, 25L))
  max_trajectory <- as.integer(.fap_cfg(cfg$max_functions_trajectory, 8L))
  exclude_top <- as.character(unlist(cfg$exclude_from_top %||% character()))
  initial_vs_final <- as.logical(.fap_cfg(cfg$initial_vs_final, TRUE))
  max_panel <- as.integer(.fap_cfg(cfg$max_functions_panel, 20L))
  x_transform <- as.character(.fap_cfg(cfg$dumbbell_x_scale, "sqrt"))

  if (!is.finite(cutoff) || cutoff <= 0 || cutoff >= 1)
    stop("FAPROTAX: adj_pval_cutoff must be between 0 and 1.", call. = FALSE)
  if (is.na(min_samples) || min_samples < 2)
    stop("FAPROTAX: min_samples_per_timepoint must be >= 2.", call. = FALSE)
  if (!is.finite(min_prevalence) || min_prevalence < 0 || min_prevalence > 1)
    stop("FAPROTAX: min_prevalence_frac must be between 0 and 1.", call. = FALSE)
  if (!is.finite(pseudocount) || pseudocount <= 0)
    stop("FAPROTAX: pseudocount must be > 0.", call. = FALSE)
  if (is.na(max_functions) || max_functions < 1)
    stop("FAPROTAX: max_functions_dumbbell must be >= 1.", call. = FALSE)
  if (is.na(max_heatmap) || max_heatmap < 1)
    stop("FAPROTAX: max_functions_heatmap must be >= 1.", call. = FALSE)
  if (is.na(max_trajectory) || max_trajectory < 1)
    stop("FAPROTAX: max_functions_trajectory must be >= 1.", call. = FALSE)
  if (is.na(initial_vs_final))
    stop("FAPROTAX: initial_vs_final must be true or false.", call. = FALSE)
  if (!x_transform %in% c("sqrt", "linear"))
    stop("FAPROTAX: dumbbell_x_scale must be sqrt or linear.", call. = FALSE)
  if (is.na(max_panel) || max_panel < 1)
    stop("FAPROTAX: max_functions_panel must be >= 1.", call. = FALSE)

  message("FAPROTAX temporal settings: BH cutoff=", cutoff,
          ", min samples/timepoint=", min_samples,
          ", prevalence=", min_prevalence,
          ", top functions/plot=", max_functions,
          ", heatmap=", max_heatmap, ", trajectories=", max_trajectory,
          ", trajectory panel=", max_panel, ", initial vs final=", initial_vs_final,
          "; fertilizers pooled; excluded from Top-N: ",
          if (length(exclude_top)) paste(exclude_top, collapse = ", ") else "none")

  scores <- faprotax_annotation$scores
  analysis_dir <- file.path(root, "analysis")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  coverage <- .fap_coverage(scores, lookup_16s)
  coverage_path <- file.path(analysis_dir, "faprotax_annotation_coverage.csv")
  coverage_plot_path <- file.path(analysis_dir, "faprotax_annotation_coverage.png")
  readr::write_csv(coverage, coverage_path)
  plot_faprotax_coverage(coverage, style, coverage_plot_path)

  written <- c(coverage_path, coverage_plot_path)
  summaries <- list()
  all_results <- list()
  idx <- 0L
  stages <- sort(unique(stats::na.omit(lookup_16s$stage)))
  for (st in stages) {
    stage_meta <- dplyr::filter(lookup_16s, stage == st)
    pairs <- consecutive_time_pairs(stage_meta$waktu)
    if (nrow(pairs) == 0) next
    definitions <- analysis_field_definitions(
      lookup_16s, "16S", st, combined_kebun)
    for (field_i in seq_len(nrow(definitions))) {
      definition <- definitions[field_i, , drop = FALSE]
      fld <- definition$field[[1]]
      field_meta <- analysis_field_subset(stage_meta, definition)
      field_pairs <- dplyr::mutate(pairs, comparison = "consecutive")
      if (initial_vs_final) {
        counts <- table(as.character(
          field_meta$waktu[field_meta$`Sample alias` %in% names(scores)]))
        field_pairs <- dplyr::bind_rows(
          field_pairs, .initial_final_pair(counts, min_samples, field_pairs))
      }
      field_results <- lapply(seq_len(nrow(field_pairs)), function(pair_i) {
        .fap_test_pair(
          scores, field_meta, st, fld, field_pairs$earlier[[pair_i]],
          field_pairs$later[[pair_i]], cutoff, min_samples,
          min_prevalence, pseudocount) |>
          # clinical/host groups stay in the CSVs but never take a Top-N slot
          dplyr::mutate(comparison = field_pairs$comparison[[pair_i]],
                        excluded_from_top = function_name %in% exclude_top)
      })
      x_max <- .fap_shared_x_max(field_results, max_functions)
      for (result in field_results) {
        output <- .write_fap_pair(result, cutoff, max_functions, style, root,
                                  x_max, x_transform)
        idx <- idx + 1L
        written <- c(written, output$files)
        summaries[[idx]] <- output$summary
        all_results[[idx]] <- result
      }
    }
  }

  summary_tbl <- dplyr::bind_rows(summaries)
  results_tbl <- dplyr::bind_rows(all_results)
  if (nrow(summary_tbl) == 0) {
    summary_tbl <- tibble::tibble(
      cohort = character(), marker = character(), stage = character(),
      field = character(), comparison = character(), earlier = character(),
      later = character(), contrast = character(), n_earlier = integer(), n_later = integer(),
      status = character(), note = character(), functions_tested = integer(),
      significant_functions = integer(), all_functions_csv = character(),
      significant_csv = character(), dumbbell_plot = character())
  }
  summary_path <- file.path(analysis_dir, "faprotax_temporal_summary.csv")
  readr::write_csv(summary_tbl, summary_path)
  written <- c(written, summary_path)

  # Cohort-level overview figures consolidate all kebun and transitions.
  panel_summaries <- list()
  for (st in stages) {
    cohort <- paste0("16S_", st)
    cohort_results <- dplyr::filter(results_tbl, stage == st)
    if (nrow(cohort_results) == 0) next
    definitions <- analysis_field_definitions(
      lookup_16s, "16S", st, combined_kebun)
    original_fields <- definitions$field[!definitions$is_combined]
    overview_results <- dplyr::filter(
      cohort_results, field %in% original_fields)
    overview_dir <- file.path(root, cohort, "Functional_FAPROTAX")
    dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)
    heatmap_path <- file.path(overview_dir, "faprotax_effect_heatmap.png")
    trajectory_path <- file.path(overview_dir, "faprotax_trajectories.png")
    selected <- .fap_select_functions(overview_results, max_trajectory)
    valid_fields <- summary_tbl |>
      dplyr::filter(stage == st, field %in% original_fields, status == "ok") |>
      dplyr::pull(field) |>
      unique()
    plot_faprotax_effect_heatmap(
      overview_results, max_heatmap, cutoff, style, heatmap_path)
    plot_faprotax_trajectories(
      scores, dplyr::filter(lookup_16s, stage == st), selected, valid_fields,
      style, trajectory_path)
    written <- c(written, heatmap_path, trajectory_path)

    # Repeat all three overview visualizations within each kebun. Function
    # selection is recalculated locally so a strong signal in one kebun does
    # not determine what is displayed for another kebun.
    for (field_i in seq_len(nrow(definitions))) {
      definition <- definitions[field_i, , drop = FALSE]
      fld <- definition$field[[1]]
      field_results <- dplyr::filter(cohort_results, field == fld)
      if (nrow(field_results) == 0) next
      field_dir <- file.path(overview_dir, .fap_safe_path(fld))
      dir.create(field_dir, recursive = TRUE, showWarnings = FALSE)
      field_heatmap <- file.path(field_dir, "faprotax_effect_heatmap.png")
      field_trajectory <- file.path(field_dir, "faprotax_trajectories.png")
      field_coverage <- file.path(field_dir, "faprotax_annotation_coverage.png")
      field_functions <- .fap_select_functions(field_results, max_trajectory)
      field_has_test <- any(
        summary_tbl$stage == st & summary_tbl$field == fld &
          summary_tbl$status == "ok")

      plot_faprotax_effect_heatmap(
        field_results, max_heatmap, cutoff, style, field_heatmap)
      plot_faprotax_trajectories(
        scores,
        analysis_field_subset(
          dplyr::filter(lookup_16s, stage == st), definition),
        field_functions,
        if (field_has_test) fld else character(),
        style, field_trajectory)
      plot_faprotax_coverage(
        analysis_field_subset(
          dplyr::filter(coverage, stage == st), definition),
        style, field_coverage)
      written <- c(written, field_heatmap, field_trajectory, field_coverage)

      panel <- .fap_trajectory(
        scores, analysis_field_subset(dplyr::filter(lookup_16s, stage == st), definition),
        st, fld, cutoff, min_samples, min_prevalence, pseudocount, exclude_top, max_panel)
      if (!is.null(panel)) {
        out <- .write_fap_trajectory(panel, cutoff, max_panel, style, field_dir)
        written <- c(written, out$files)
        panel_summaries[[length(panel_summaries) + 1L]] <- out$summary
      }
    }
  }

  panel_tbl <- dplyr::bind_rows(panel_summaries)
  panel_path <- file.path(analysis_dir, "faprotax_trajectory_summary.csv")
  readr::write_csv(panel_tbl, panel_path)
  written <- c(written, panel_path)
  message("FAPROTAX temporal analysis complete: ", nrow(summary_tbl),
          " kebun/time comparisons, ", nrow(panel_tbl), " trajectory panels.")
  list(files = unique(written), summary = summary_tbl,
       summary_csv = summary_path, coverage_csv = coverage_path)
}
