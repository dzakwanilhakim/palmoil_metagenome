# =============================================================================
# R/functions_qc_apply.R — Stage 3, Steps 3-5
#   Step 3: apply thresholds (sample depth, OTU min reads, prevalence)
#   Step 4: post-QC assessment (summary + boxplot)
#   Step 5: synchronized Rarefied + CLR tables (identical samples AND OTUs)
# Per marker. Reads config/thresholds.yaml.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

load_thresholds <- function(path = "config/thresholds.yaml") {
  yaml::read_yaml(path)
}

# ----------------------------------------------------------------------------
# STEP 3. APPLY THRESHOLDS  (per marker, strict order)
#   3a drop samples below min_sample_depth
#   3b on surviving samples: drop OTUs with global total reads < otu_min_reads
#   3c drop OTUs present in < prevalence_frac of surviving samples (count > 0)
# Returns a long table (filtered) + records the master sample list per marker.
# ----------------------------------------------------------------------------
apply_thresholds <- function(tax_filtered, thresholds) {
  markers <- unique(tax_filtered$marker)

  out <- purrr::map_dfr(markers, function(mk) {
    th <- thresholds$markers[[mk]]
    if (is.null(th)) stop("No thresholds defined for marker ", mk)

    dat <- dplyr::filter(tax_filtered, marker == mk)

    # 3a — sample depth
    depth <- dat |>
      dplyr::group_by(id_sampel) |>
      dplyr::summarise(depth = sum(count, na.rm = TRUE), .groups = "drop")
    keep_samples <- depth$id_sampel[depth$depth >= th$min_sample_depth]
    dat <- dplyr::filter(dat, id_sampel %in% keep_samples)
    message(sprintf("[%s] 3a sample depth >= %s : kept %d samples",
                    mk, format(th$min_sample_depth, big.mark=","),
                    length(keep_samples)))

    # 3b — global OTU minimum reads (on surviving samples)
    otu_tot <- dat |>
      dplyr::group_by(taxon) |>
      dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop")
    keep_otu1 <- otu_tot$taxon[otu_tot$total >= th$otu_min_reads]
    dat <- dplyr::filter(dat, taxon %in% keep_otu1)

    # 3c — prevalence (>0 in >= prevalence_frac of surviving samples)
    n_samp <- length(keep_samples)
    min_prev <- ceiling(th$prevalence_frac * n_samp)
    prev <- dat |>
      dplyr::filter(count > 0) |>
      dplyr::group_by(taxon) |>
      dplyr::summarise(prev = dplyr::n_distinct(id_sampel), .groups = "drop")
    keep_otu2 <- prev$taxon[prev$prev >= min_prev]
    dat <- dplyr::filter(dat, taxon %in% keep_otu2)

    message(sprintf("[%s] 3b OTU reads >= %d & 3c prevalence >= %d/%d (%.0f%%) : kept %d OTUs",
                    mk, th$otu_min_reads, min_prev, n_samp,
                    100*th$prevalence_frac, dplyr::n_distinct(dat$taxon)))
    dat
  })
  out
}

# Master sample list (the synchronized survivors), per marker.
# Drops samples with missing/blank fertilizer (cannot enter any A/B/C/D test).
master_sample_list <- function(thresholded_long) {
  ms <- thresholded_long |>
    dplyr::distinct(marker, id_sampel, stage, field, fertilizer, timepoint)
  bad <- is.na(ms$fertilizer) |
         stringr::str_trim(as.character(ms$fertilizer)) %in% c("", "NA")
  if (any(bad)) {
    message("master_sample_list: dropped ", sum(bad),
            " sample(s) with NA/blank fertilizer.")
    print(ms[bad, c("marker","id_sampel","stage","field","timepoint")])
  }
  ms[!bad, , drop = FALSE]
}

# ----------------------------------------------------------------------------
# QC DROPPED-SAMPLE REPORT
#   Compares samples present BEFORE thresholding (tax_filtered) against the
#   survivors (master list), per marker, and records why each was dropped:
#     - "below_min_depth" : sample depth < min_sample_depth
#     - "lost_after_otu_filter" : passed depth but had no reads left once
#        low-count / low-prevalence OTUs were removed (rare, but possible)
#   Writes a tidy CSV and returns the path (file target).
# ----------------------------------------------------------------------------
write_qc_drop_report <- function(tax_filtered, thresholded_long, thresholds,
                                 out_path = "results/dropped_samples_qc.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  # depth per sample BEFORE thresholding, with metadata
  pre <- tax_filtered |>
    dplyr::group_by(marker, id_sampel, stage, field, fertilizer, timepoint) |>
    dplyr::summarise(depth = sum(count, na.rm = TRUE), .groups = "drop")

  survivors <- thresholded_long |>
    dplyr::distinct(marker, id_sampel)

  dropped <- dplyr::anti_join(pre, survivors, by = c("marker", "id_sampel"))

  # annotate reason using each marker's min_sample_depth
  min_depth_of <- function(mk) thresholds$markers[[mk]]$min_sample_depth
  dropped <- dropped |>
    dplyr::rowwise() |>
    dplyr::mutate(
      min_sample_depth = min_depth_of(marker),
      reason = dplyr::if_else(depth < min_sample_depth,
                              "below_min_depth",
                              "lost_after_otu_filter")) |>
    dplyr::ungroup() |>
    dplyr::arrange(marker, reason, depth) |>
    dplyr::select(marker, id_sampel, stage, field, fertilizer, timepoint,
                  depth, min_sample_depth, reason)

  readr::write_csv(dropped, out_path)
  message("QC drop report: ", nrow(dropped), " samples dropped -> ", out_path)
  # brief per-marker/reason tally to console
  if (nrow(dropped))
    print(dplyr::count(dropped, marker, reason))
  out_path
}

# ----------------------------------------------------------------------------
# UNIFIED FILTERED-BARCODE REPORT
#   One CSV per the whole pipeline listing every barcode that was filtered out,
#   at any stage, with metadata attached where available and the reason.
#   Reasons:
#     no_metadata_match     - barcode in matrix but no metadata row
#     no_matrix_match       - barcode in metadata but no matrix data
#     flag_qc               - metadata Flag == 1
#     below_min_depth       - dropped at Stage 3 sample-depth threshold
#     lost_after_otu_filter - passed depth but emptied by OTU filtering
#   For barcodes with no metadata, only barcode + marker + reason are written.
# ----------------------------------------------------------------------------
write_filtered_barcode_report <- function(long_counts, metadata_all, joined_counts,
                                          tax_filtered, thresholded, schema, thresholds,
                                          out_path = "results/filtered_barcodes_all.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  registry <- build_sample_registry(metadata_all, schema)

  # helper: metadata lookup per marker keyed by barcode string
  meta_by_barcode <- purrr::map_dfr(names(schema$marker_columns), function(mk) {
    bc <- schema$marker_columns[[mk]]$barcode
    fl <- schema$marker_columns[[mk]]$flag
    if (!bc %in% names(registry)) return(tibble::tibble())
    registry |>
      dplyr::filter(!is.na(.data[[bc]]), .data[[bc]] != "NA") |>
      dplyr::transmute(marker = mk, barcode = .data[[bc]],
                       id_sampel, stage, field, fertilizer, timepoint,
                       extraction_batch,
                       flag = if (fl %in% names(registry)) .data[[fl]] else NA_character_)
  })

  matched <- if (nrow(joined_counts))
    dplyr::distinct(joined_counts, marker, barcode) else
    tibble::tibble(marker=character(), barcode=character())

  # --- 1. matrix barcodes with NO metadata (Stage 2) ---
  matrix_bc <- long_counts |> dplyr::distinct(marker, barcode)
  no_meta <- dplyr::anti_join(matrix_bc, matched, by = c("marker","barcode")) |>
    dplyr::mutate(reason = "no_metadata_match")

  # --- 2. metadata barcodes with NO matrix (Stage 2) ---
  no_mat <- dplyr::anti_join(meta_by_barcode, matrix_bc, by = c("marker","barcode")) |>
    dplyr::mutate(reason = "no_matrix_match")

  # --- 3. flagged at QC (Flag == 1), among joined ---
  drop_val <- as.character(schema$qc$drop_when_flag_equals)
  flagged <- joined_counts |>
    dplyr::distinct(marker, barcode, id_sampel, stage, field, fertilizer,
                    timepoint, extraction_batch, flag) |>
    dplyr::filter(!is.na(flag) & as.character(flag) == drop_val) |>
    dplyr::mutate(reason = "flag_qc")

  # --- 4. dropped at Stage 3 thresholding ---
  # samples present pre-threshold (tax_filtered) but not in survivors
  pre <- tax_filtered |>
    dplyr::group_by(marker, id_sampel, stage, field, fertilizer, timepoint) |>
    dplyr::summarise(depth = sum(count, na.rm = TRUE), .groups = "drop")
  survivors <- dplyr::distinct(thresholded, marker, id_sampel)
  stage3 <- dplyr::anti_join(pre, survivors, by = c("marker","id_sampel")) |>
    dplyr::rowwise() |>
    dplyr::mutate(reason = dplyr::if_else(
      depth < thresholds$markers[[marker]]$min_sample_depth,
      "below_min_depth", "lost_after_otu_filter")) |>
    dplyr::ungroup()
  # attach the barcode for these (from joined_counts)
  bc_lookup <- joined_counts |> dplyr::distinct(marker, id_sampel, barcode,
                                                extraction_batch, flag)
  stage3 <- stage3 |>
    dplyr::left_join(bc_lookup, by = c("marker","id_sampel"))

  # --- combine, unify columns ---
  cols <- c("marker","barcode","id_sampel","stage","field","fertilizer",
            "timepoint","extraction_batch","depth","flag","reason")
  pad <- function(df) {
    for (c in setdiff(cols, names(df))) df[[c]] <- NA
    df[, cols]
  }

  report <- dplyr::bind_rows(
    pad(no_meta), pad(no_mat), pad(flagged), pad(stage3)
  ) |>
    dplyr::arrange(marker, reason, barcode)

  readr::write_csv(report, out_path)
  message("Filtered-barcode report: ", nrow(report), " entries -> ", out_path)
  print(dplyr::count(report, marker, reason))
  out_path
}

# ----------------------------------------------------------------------------
# STEP 4. POST-QC ASSESSMENT
# ----------------------------------------------------------------------------
postqc_summary <- function(thresholded_long) {
  d <- thresholded_long |>
    dplyr::group_by(marker, id_sampel, stage) |>
    dplyr::summarise(depth = sum(count, na.rm = TRUE), .groups = "drop")
  summ <- d |>
    dplyr::group_by(marker) |>
    dplyr::summarise(n_samples = dplyr::n(),
                     Min = min(depth), Q1 = quantile(depth, .25),
                     Median = median(depth), Q3 = quantile(depth, .75),
                     Max = max(depth), .groups = "drop")
  message("\n=== POST-QC DEPTH SUMMARY (per marker) ==="); print(summ)
  list(depth = d, summary = summ)
}

plot_postqc_boxplot <- function(depth_tbl,
                                out_path = "results/postqc_depth_boxplot.png") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  pal <- c(Nursery = "#E07A5F", TM = "#3D8C8C")
  p <- ggplot2::ggplot(depth_tbl, ggplot2::aes(marker, depth)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = stage), width = .55,
                          outlier.shape = NA, alpha = .85,
                          position = ggplot2::position_dodge(.7),
                          linewidth = .4, colour = "grey25") +
    ggplot2::geom_point(ggplot2::aes(fill = stage),
                        position = ggplot2::position_jitterdodge(.12, dodge.width=.7),
                        shape = 21, size = 1.6, stroke = .2, colour="grey20", alpha=.55) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::annotation_logticks(sides = "l", colour = "grey60") +
    ggplot2::scale_fill_manual(values = pal, name = "Stage") +
    ggplot2::labs(title = "Post-QC sample sequencing depth",
                  subtitle = "After threshold filtering (samples + OTUs)",
                  x = NULL, y = "Total reads (log scale)") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour="grey40", size=10))
  ggplot2::ggsave(out_path, p, width = 7.5, height = 5.5, dpi = 300)
  message("Wrote ", out_path)
  out_path
}

# ----------------------------------------------------------------------------
# STEP 5. SYNCHRONIZED RAREFIED + CLR TABLES
#   Build a sample x OTU count matrix per marker from the SAME thresholded data.
#   - Rarefied: vegan::rrarefy to rarefy_depth (samples below it are dropped;
#     those dropped samples are ALSO removed from the CLR table -> strict sync).
#   - CLR: pseudocount + log-ratio, on the EXACT same surviving samples & OTUs.
# Returns, per marker: $rarefied (matrix), $clr (matrix), $samples (kept ids).
# ----------------------------------------------------------------------------
.long_to_matrix <- function(dat) {
  wide <- dat |>
    dplyr::group_by(id_sampel, taxon) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = taxon, values_from = count, values_fill = 0)
  m <- as.matrix(wide[, -1]); rownames(m) <- wide$id_sampel
  storage.mode(m) <- "double"
  m
}

clr_transform <- function(mat, pseudocount = 1) {
  x <- mat + pseudocount
  logx <- log(x)
  gm <- rowMeans(logx)              # geometric mean per sample (in log space)
  sweep(logx, 1, gm, "-")           # CLR = log(x) - mean(log(x)) per sample
}

# post-QC rarefaction curve per marker, with a vertical line at the rarefy
# depth actually used (so you can see whether curves saturate at/below it).
plot_postqc_rarefaction <- function(thresholded_long, thresholds, marker,
                                     step = 200, out_dir = "results") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir,
                        paste0("postqc_rarefaction_", tolower(marker), ".png"))
  dat <- dplyr::filter(thresholded_long, marker == !!marker)
  if (nrow(dat) == 0) { message("No data for ", marker); return(NA_character_) }

  wide <- dat |>
    dplyr::group_by(id_sampel, taxon) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    tidyr::pivot_wider(names_from = taxon, values_from = count, values_fill = 0)
  m <- as.matrix(wide[, -1]); rownames(m) <- wide$id_sampel
  storage.mode(m) <- "integer"

  stage_of <- dat |> dplyr::distinct(id_sampel, stage) |>
    dplyr::arrange(match(id_sampel, wide$id_sampel))
  cols <- as.integer(factor(stage_of$stage))

  th <- thresholds$markers[[marker]]
  rdepth <- th$rarefy_depth
  if (is.null(rdepth) || is.na(rdepth)) rdepth <- min(rowSums(m))

  grDevices::png(out_path, width = 1400, height = 900, res = 150)
  vegan::rarecurve(m, step = step, label = FALSE, col = cols,
                   xlab = "Sequencing depth", ylab = "Observed genera",
                   main = paste0("Post-QC rarefaction - ", marker))
  graphics::abline(v = rdepth, lty = 2, col = "firebrick", lwd = 2)
  graphics::legend("bottomright",
                   legend = levels(factor(stage_of$stage)),
                   col = seq_along(levels(factor(stage_of$stage))),
                   lty = 1, bty = "n")
  grDevices::dev.off()
  message("Wrote ", out_path, "  (rarefy depth line at ", format(rdepth, big.mark=","), ")")
  out_path
}

build_rarefied_and_clr <- function(thresholded_long, thresholds, seed = 42) {
  markers <- unique(thresholded_long$marker)
  set.seed(seed)

  purrr::map(markers, function(mk) {
    th  <- thresholds$markers[[mk]]
    dat <- dplyr::filter(thresholded_long, marker == mk)
    m   <- .long_to_matrix(dat)
    m_int <- m; storage.mode(m_int) <- "integer"

    depths <- rowSums(m)
    rdepth <- th$rarefy_depth
    if (is.null(rdepth) || is.na(rdepth)) rdepth <- min(depths)  # = min surviving

    # samples that meet the rarefaction depth (should be all, since min-depth
    # threshold >= rarefy depth, but guard strictly)
    keep <- names(depths)[depths >= rdepth]

    rar <- vegan::rrarefy(m_int[keep, , drop = FALSE], sample = rdepth)
    # drop OTUs that became all-zero after rarefying (keep tables consistent)
    nz  <- colSums(rar) > 0
    rar <- rar[, nz, drop = FALSE]

    # CLR on the SAME samples and SAME OTU set as the rarefied table
    clr <- clr_transform(m[keep, colnames(rar), drop = FALSE],
                         pseudocount = thresholds$clr$pseudocount)

    message(sprintf("[%s] rarefy depth=%s | samples=%d | OTUs=%d (rarefied & CLR synchronized)",
                    mk, format(rdepth, big.mark=","), length(keep), ncol(rar)))

    list(marker = mk, rarefy_depth = rdepth,
         samples = keep, rarefied = rar, clr = clr)
  }) |> setNames(markers)
}
