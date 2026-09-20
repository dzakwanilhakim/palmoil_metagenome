# =============================================================================
# R/functions_gold_da_temporal.R — consecutive temporal ANCOM-BC2 per kebun
#
# Unit of analysis:
#   cohort (marker x stage) x field (Kode Kebun) x rank (genus/species)
#
# All fertilizers are pooled within a field. One model is fitted across the
# usable timepoints in that field, then only consecutive contrasts are
# exported (T0 vs T1, T1 vs T2, ...). Each contrast gets an all-taxa CSV, a
# significant-taxa CSV, and a dumbbell plot. The post-gold-QC cohort gate in
# _targets.R determines which cohorts reach this driver.
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

.run_temporal_field_rank <- function(processed_mat, lookup, marker, stage,
                                     field, field_members, rank, pairs, cutoff,
                                     min_samples, pseudo_sensitivity, n_cl) {
  cohort <- paste0(marker, "_", stage)
  sample_cols <- .processed_matrix_sample_columns(processed_mat)
  meta <- lookup |>
    dplyr::filter(stage == !!stage, .data$field %in% field_members,
                  `Sample alias` %in% sample_cols) |>
    dplyr::mutate(field = !!field)
  counts <- table(as.character(meta$waktu))
  pair_counts <- pairs |>
    dplyr::mutate(
      n_earlier = as.integer(counts[earlier]),
      n_later = as.integer(counts[later]),
      n_earlier = tidyr::replace_na(n_earlier, 0L),
      n_later = tidyr::replace_na(n_later, 0L),
      usable = n_earlier >= min_samples & n_later >= min_samples)

  skipped <- purrr::pmap_dfr(
    dplyr::filter(pair_counts, !usable),
    function(earlier, later, n_earlier, n_later, usable) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "skipped_insufficient_replication",
        paste0("Requires >=", min_samples, " samples at each timepoint"))
    })

  usable_pairs <- dplyr::filter(pair_counts, usable)
  if (nrow(usable_pairs) == 0)
    return(skipped)

  model_times <- .ordered_timepoints(c(usable_pairs$earlier, usable_pairs$later))
  model_meta <- dplyr::filter(meta, waktu %in% model_times)
  ps <- .temporal_phyloseq(processed_mat, model_meta, model_times)

  if (phyloseq::ntaxa(ps) < 2) {
    failed <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                     n_later, usable) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "skipped_insufficient_taxa", "<2 non-zero taxa")
    })
    return(dplyr::bind_rows(failed, skipped))
  }

  context <- paste(cohort, field, rank, paste(model_times, collapse = "-"), sep = "|")
  message("ANCOM-BC2 temporal [", context, "]: samples=", phyloseq::nsamples(ps),
          ", taxa=", phyloseq::ntaxa(ps), ", fertilizers pooled")

  out <- tryCatch(
    ANCOMBC::ancombc2(
      data = ps, fix_formula = "waktu", group = "waktu",
      p_adj_method = "BH", pseudo_sens = pseudo_sensitivity,
      prv_cut = 0, lib_cut = 0, struc_zero = FALSE, neg_lb = FALSE,
      alpha = cutoff, n_cl = n_cl, verbose = FALSE,
      global = FALSE, pairwise = length(model_times) > 2,
      dunnet = FALSE, trend = FALSE),
    error = function(e) e)

  if (inherits(out, "error")) {
    failed <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                     n_later, usable) {
      .temporal_status_result(
        cohort, marker, stage, field, rank, earlier, later,
        n_earlier, n_later, "ancombc2_failed", conditionMessage(out))
    })
    return(dplyr::bind_rows(failed, skipped))
  }

  tested <- purrr::pmap_dfr(usable_pairs, function(earlier, later, n_earlier,
                                                   n_later, usable) {
    .extract_temporal_contrast(
      out, cohort, marker, stage, field, rank, earlier, later,
      n_earlier, n_later, cutoff, pseudo_sensitivity)
  })
  dplyr::bind_rows(tested, skipped)
}

.safe_da_path <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

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
    reason <- stringr::str_wrap(reason, width = 90)
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = reason,
                        size = 5, fontface = "bold", colour = "grey35") +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
      gold_plot_theme(style) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
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
    dplyr::select(cohort, marker, stage, field, rank, earlier, later, contrast,
                  n_earlier, n_later) |>
    dplyr::mutate(
      status = if (any(result$status == "ok")) "ok" else result$status[[1]],
      note = paste(unique(stats::na.omit(result$note)), collapse = "; "),
      taxa_tested = sum(result$status == "ok" & !is.na(result$taxon)),
      significant_taxa = sum(result$status == "ok" & result$significant, na.rm = TRUE),
      all_taxa_csv = all_csv, significant_csv = sig_csv, dumbbell_plot = plot_path)
  list(files = c(all_csv, sig_csv, plot_path), summary = summary)
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

  message("Temporal DA settings: BH cutoff=", cutoff,
          ", min samples/timepoint=", min_samples,
          ", pseudo sensitivity=", pseudo_sensitivity,
          ", ANCOM-BC2 cores=", n_cl)

  specs <- list(
    list(marker = "16S", rank = "genus", matrix = gold_processed_matrix_16s_genus,
         lookup = gold_universe_lookup_16s),
    list(marker = "16S", rank = "species", matrix = gold_processed_matrix_16s_species,
         lookup = gold_universe_lookup_16s),
    list(marker = "ITS", rank = "genus", matrix = gold_processed_matrix_its_genus,
         lookup = gold_universe_lookup_its),
    list(marker = "ITS", rank = "species", matrix = gold_processed_matrix_its_species,
         lookup = gold_universe_lookup_its))

  written <- character(); summaries <- list(); idx <- 0L
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
        members <- as.character(unlist(
          definitions$members[[field_i]], use.names = FALSE))
        result <- .run_temporal_field_rank(
          spec$matrix, spec$lookup, spec$marker, st, fld, members,
          spec$rank, pairs,
          cutoff, min_samples, pseudo_sensitivity, n_cl)
        for (pair_i in seq_len(nrow(pairs))) {
          pair_result <- dplyr::filter(
            result, earlier == pairs$earlier[[pair_i]], later == pairs$later[[pair_i]])
          if (nrow(pair_result) == 0) next
          output <- .write_temporal_pair(pair_result, cutoff, max_taxa, style, root)
          idx <- idx + 1L
          written <- c(written, output$files)
          summaries[[idx]] <- output$summary
        }
      }
    }
  }

  summary_tbl <- dplyr::bind_rows(summaries)
  summary_path <- file.path(root, "analysis", "da_temporal_summary.csv")
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary_tbl, summary_path)
  written <- c(written, summary_path)
  message("Temporal DA complete: ", nrow(summary_tbl), " field/rank/time comparisons.")
  list(files = unique(written), summary = summary_tbl, summary_csv = summary_path)
}
