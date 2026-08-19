# =============================================================================
# R/functions_gold_da.R  —  Analysis output: differential abundance (ANCOM-BC)
#   Runs on the UNRAREFIED, thresholded counts (gold_processed_matrix) --
#   ANCOM-BC does its own bias-correction normalization and should not be fed
#   rarefied data. Independently at genus and species rank.
#
#   Both cases test the WAKTU (time) effect -- the difference is pooling
#   scope, not the variable being tested:
#   Case 1: within one field (Kode Kebun), POOL all fertilizers (Kode Pupuk)
#           together, test WAKTU -- the broad field-level temporal trend.
#   Case 2: within one field AND ONE fertilizer, test WAKTU -- the
#           fertilizer-specific temporal trend within that field.
#
#   Significant = adj p (BH) < config's da$adj_pval_cutoff (0.05 default).
#
#   NOTE: this has not been run against a live ANCOMBC install yet (the
#   package wasn't present in this environment as of writing) -- the exact
#   column-naming convention in ancombc2()'s $res output can shift slightly
#   between versions, so run_ancombc_gold() is written defensively (regex
#   column discovery + a diagnostic dump if nothing matches) and should be
#   treated as a first pass to verify against real output.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
})

# gold_processed_matrix (wide: taxon, <barcodes>, total, lineage cols) ->
# phyloseq object restricted to a set of sample ids, with sample_data from
# universe_lookup (field/fertilizer/waktu).
.processed_matrix_to_phyloseq <- function(processed_mat, sample_ids, meta) {
  taxon_col <- names(processed_mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(processed_mat))
  sample_cols <- setdiff(names(processed_mat), c(taxon_col, meta_cols))
  keep <- intersect(sample_cols, sample_ids)
  if (length(keep) < 4) return(NULL)

  otu <- as.matrix(processed_mat[, keep, drop = FALSE])
  rownames(otu) <- processed_mat[[taxon_col]]
  storage.mode(otu) <- "double"
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  if (nrow(otu) < 2) return(NULL)

  samp <- as.data.frame(meta[match(keep, meta$`Sample alias`), , drop = FALSE])
  rownames(samp) <- keep

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(samp))
}

# run ancombc2, return tidy significant taxon x contrast rows (q < cutoff).
# Column names in ancombc2()$res follow lfc_<term>/q_<term> per contrast --
# discovered by regex rather than hardcoded, since exact suffixes depend on
# factor levels and whether `pairwise` mode was used.
run_ancombc_gold <- function(ps, fix_formula, group_var, adj_pval_cutoff, context) {
  if (is.null(ps) || phyloseq::nsamples(ps) < 4)
    return(tibble::tibble(context = context, taxon = NA, contrast = NA,
                          lfc = NA, q_val = NA, note = "skipped: <4 samples"))
  n_levels <- dplyr::n_distinct(phyloseq::sample_data(ps)[[group_var]])
  if (n_levels < 2)
    return(tibble::tibble(context = context, taxon = NA, contrast = NA,
                          lfc = NA, q_val = NA,
                          note = paste0("skipped: <2 levels in ", group_var)))

  out <- tryCatch(
    ANCOMBC::ancombc2(
      data = ps, fix_formula = fix_formula, group = group_var,
      p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
      global = n_levels > 2, pairwise = n_levels > 2,
      alpha = adj_pval_cutoff, verbose = FALSE),
    error = function(e) {
      message("  ancombc2 failed [", context, "]: ", conditionMessage(e))
      NULL
    })
  if (is.null(out))
    return(tibble::tibble(context = context, taxon = NA, contrast = NA,
                          lfc = NA, q_val = NA, note = "ancombc2 failed"))

  res <- out$res
  lfc_cols <- grep(paste0("^lfc_", group_var), names(res), value = TRUE)
  q_cols   <- grep(paste0("^q_",   group_var), names(res), value = TRUE)
  if (length(lfc_cols) == 0 || length(q_cols) == 0) {
    message("  [", context, "] no lfc_/q_ columns matched '", group_var,
            "' -- res columns: ", paste(names(res), collapse = ", "))
    return(tibble::tibble(context = context, taxon = NA, contrast = NA,
                          lfc = NA, q_val = NA,
                          note = "no matching lfc_/q_ columns (see message log)"))
  }

  purrr::map2_dfr(lfc_cols, q_cols, function(lc, qc) {
    tibble::tibble(context = context, taxon = res$taxon,
                   contrast = sub("^lfc_", "", lc),
                   lfc = res[[lc]], q_val = res[[qc]], note = NA_character_) |>
      dplyr::filter(!is.na(q_val), q_val < adj_pval_cutoff)
  })
}

# ----------------------------------------------------------------------------
# Case 1: per field, ALL fertilizers pooled, waktu (time) effect
# ----------------------------------------------------------------------------
build_da_case1 <- function(processed_mat, universe_lookup, marker, stage, rank,
                           adj_pval_cutoff) {
  meta <- dplyr::filter(universe_lookup, stage == !!stage)
  fields <- sort(unique(meta$field))
  purrr::map_dfr(fields, function(fld) {
    ids <- meta$`Sample alias`[meta$field == fld]   # every fertilizer pooled
    ps <- .processed_matrix_to_phyloseq(processed_mat, ids, meta)
    run_ancombc_gold(ps, fix_formula = "waktu", group_var = "waktu",
                     adj_pval_cutoff = adj_pval_cutoff,
                     context = paste(marker, stage, rank, "Case1", fld, sep = "|"))
  })
}

# ----------------------------------------------------------------------------
# Case 2: per field x fertilizer, waktu (time) effect
# ----------------------------------------------------------------------------
build_da_case2 <- function(processed_mat, universe_lookup, marker, stage, rank,
                           adj_pval_cutoff) {
  meta <- dplyr::filter(universe_lookup, stage == !!stage)
  combos <- dplyr::distinct(meta, field, fertilizer)
  purrr::pmap_dfr(combos, function(field, fertilizer) {
    ids <- meta$`Sample alias`[meta$field == field & meta$fertilizer == fertilizer]
    ps <- .processed_matrix_to_phyloseq(processed_mat, ids, meta)
    run_ancombc_gold(ps, fix_formula = "waktu", group_var = "waktu",
                     adj_pval_cutoff = adj_pval_cutoff,
                     context = paste(marker, stage, rank, "Case2", field, fertilizer,
                                     sep = "|"))
  })
}

# ----------------------------------------------------------------------------
# Log-fold-change dumbbell plot for either case -- one per (marker, stage,
# rank, case), pooling all its significant-taxon contexts. Dumbbell runs
# from 0 (no change) to the LFC point estimate; colour = direction of change
# over time (matches the established DA volcano plot's up/down palette in
# functions_network.R).
#   case_label: "Case 1" (context = marker|stage|rank|Case1|field, no
#               fertilizer segment -- all fertilizers were pooled) or
#               "Case 2" (context = ...|Case2|field|fertilizer).
# ----------------------------------------------------------------------------
plot_da_dumbbell <- function(da_rows, marker, stage, rank, case_label,
                             style = load_plot_style(), out_dir = "Results/analysis") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  case_tag <- gsub("[^A-Za-z0-9]", "", case_label)   # "Case1" / "Case2"
  out_path <- file.path(out_dir,
    paste0("da_dumbbell_", tolower(case_tag), "_", tolower(marker), "_",
          tolower(stage), "_", rank, ".png"))

  d <- dplyr::filter(da_rows, !is.na(taxon))
  if (nrow(d) == 0) {
    message("  DA dumbbell skipped [", marker, "/", stage, "/", rank, "/", case_label,
            "]: no significant taxa")
    return(NA_character_)
  }

  # contrast distinguishes multiple waktu pairs when >2 timepoints exist,
  # so it's part of the label too, to avoid collisions.
  ctx_parts <- stringr::str_split(d$context, "\\|")
  d$field      <- purrr::map_chr(ctx_parts, ~ .x[5])
  d$fertilizer <- purrr::map_chr(ctx_parts,
                                 ~ if (length(.x) >= 6) .x[6] else NA_character_)
  d$context_label <- ifelse(is.na(d$fertilizer), d$field,
                            paste0(d$field, "/", d$fertilizer))
  d$direction  <- ifelse(d$lfc >= 0, "Higher at later timepoint",
                                     "Lower at later timepoint")
  d$label <- paste0(d$taxon, "  (", d$context_label, ", ", d$contrast, ")")
  d <- dplyr::arrange(d, lfc)
  d$label <- factor(d$label, levels = unique(d$label))

  pal <- c("Higher at later timepoint" = "#E74C3C",
           "Lower at later timepoint"  = "#2980B9")

  subt <- if (case_tag == "Case1")
    "All Kode Pupuk pooled per field -- broad field-level temporal trend" else
    "Within one Kode Pupuk per field -- fertilizer-specific temporal trend"

  p <- ggplot2::ggplot(d, ggplot2::aes(y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = lfc, yend = label,
                                       colour = direction), linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(x = 0), size = 2, colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(x = lfc, colour = direction), size = 3.2) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(
      title = paste0(marker, "_", stage, " — ", stringr::str_to_title(rank),
                     " — ", case_label, " (temporal) ANCOM-BC log fold change"),
      subtitle = paste0(subt, " | grey dot = no change (0); coloured dot = ",
                        "significant LFC (q < cutoff)"),
      x = "Log fold change (ANCOM-BC)", y = NULL) +
    gold_plot_theme(style) +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(out_path, p, width = 9,
                  height = max(4, 0.35 * dplyr::n_distinct(d$label) + 2),
                  dpi = 200, limitsize = FALSE)
  message("  Wrote ", out_path)
  out_path
}

# ----------------------------------------------------------------------------
# MAIN driver: Case 1 + Case 2, genus + species, both markers, plus a
# dumbbell plot for each case
# ----------------------------------------------------------------------------
build_gold_da_analysis <- function(gold_processed_matrix_16s_genus,
                                   gold_processed_matrix_16s_species,
                                   gold_universe_lookup_16s,
                                   gold_processed_matrix_its_genus,
                                   gold_processed_matrix_its_species,
                                   gold_universe_lookup_its,
                                   analysis_thresholds, style = load_plot_style(),
                                   root = "Results/analysis") {
  cutoff <- analysis_thresholds$da$adj_pval_cutoff
  message("DA (ANCOM-BC) thresholds: adj_pval_cutoff = ", cutoff)

  run_one <- function(processed_mat, lookup, marker, rank) {
    stages <- sort(unique(lookup$stage))
    purrr::map(stages, function(st) {
      list(marker = marker, stage = st, rank = rank,
           case1 = build_da_case1(processed_mat, lookup, marker, st, rank, cutoff),
           case2 = build_da_case2(processed_mat, lookup, marker, st, rank, cutoff))
    })
  }

  all_runs <- c(
    run_one(gold_processed_matrix_16s_genus,   gold_universe_lookup_16s, "16S", "genus"),
    run_one(gold_processed_matrix_16s_species, gold_universe_lookup_16s, "16S", "species"),
    run_one(gold_processed_matrix_its_genus,   gold_universe_lookup_its, "ITS", "genus"),
    run_one(gold_processed_matrix_its_species, gold_universe_lookup_its, "ITS", "species"))

  res <- dplyr::bind_rows(purrr::map(all_runs, "case1"), purrr::map(all_runs, "case2"))
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  csv_path <- file.path(root, "da_ancombc_significant.csv")
  readr::write_csv(res, csv_path)
  message("Wrote ", csv_path, " (", nrow(res), " significant taxon-contrast rows)")

  plots <- c(
    purrr::map_chr(all_runs, function(r)
      plot_da_dumbbell(r$case1, r$marker, r$stage, r$rank, "Case 1", style, out_dir = root)),
    purrr::map_chr(all_runs, function(r)
      plot_da_dumbbell(r$case2, r$marker, r$stage, r$rank, "Case 2", style, out_dir = root)))
  plots <- plots[!is.na(plots)]

  list(csv = csv_path, plots = plots)
}
