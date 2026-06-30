# =============================================================================
# R/functions_ingest.R  —  Stage 2: Ingest & Validation
# Pure functions. No targets here; the DAG in _targets.R calls these.
# Everything is driven by config/schema.yaml — no hardcoded column names.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

# ----------------------------------------------------------------------------
# 0. Schema
# ----------------------------------------------------------------------------
load_schema <- function(path = "config/schema.yaml") {
  yaml::read_yaml(path)
}

# ----------------------------------------------------------------------------
# 1. METADATA: read all per-batch tsv, rename to canonical, bind
#    `files` = character vector of metadata tsv paths (16s and/or its).
# ----------------------------------------------------------------------------
compile_metadata <- function(files, schema) {
  cmap <- schema$column_map

  read_one <- function(f) {
    df <- readr::read_tsv(f, show_col_types = FALSE,
                          name_repair = "minimal", progress = FALSE)
    # rename only the columns we know about; leave the rest as-is
    present <- intersect(names(cmap), names(df))
    rename_vec <- setNames(present, unlist(cmap[present]))  # new = old
    df <- dplyr::rename(df, dplyr::all_of(rename_vec))
    df$.source_file <- basename(f)
    df
  }

  meta <- purrr::map(files, read_one) |> dplyr::bind_rows()

  # normalise obvious whitespace issues that silently break joins/validation
  chr_cols <- c("stage", "field", "fertilizer", "timepoint")
  for (c in intersect(chr_cols, names(meta))) {
    meta[[c]] <- stringr::str_trim(as.character(meta[[c]]))
  }
  meta
}

# ----------------------------------------------------------------------------
# 2. MATRIX: parse filename -> (marker, seq_batch, rank); melt wide -> long
#    Returns long table with columns: barcode, taxon, rank, marker,
#    seq_batch, count, tax (lineage), phylum (for the phylum files).
# ----------------------------------------------------------------------------
parse_matrix_filename <- function(path, schema) {
  rx <- schema$join$matrix_filename_regex
  m  <- stringr::str_match(basename(path), rx)
  if (any(is.na(m)))
    stop("Matrix filename does not match schema regex: ", basename(path))
  list(marker = toupper(m[2]), seq_batch = as.integer(m[3]), rank = m[4])
}

melt_matrix <- function(path, schema) {
  info <- parse_matrix_filename(path, schema)
  non_sample <- schema$join$non_sample_columns

  mat <- readr::read_tsv(path, show_col_types = FALSE,
                         name_repair = "minimal", progress = FALSE)

  # the taxon label column is the rank name itself (genus / phylum)
  taxon_col <- info$rank
  if (!taxon_col %in% names(mat))
    stop("Expected taxon column '", taxon_col, "' not found in ", basename(path))

  # sample (barcode) columns = everything not in the known non-sample set
  sample_cols <- setdiff(names(mat), non_sample)
  if (length(sample_cols) == 0)
    stop("No barcode/sample columns detected in ", basename(path))

  keep_tax <- intersect(c("tax", "phylum"), names(mat))

  long <- mat |>
    dplyr::select(dplyr::all_of(c(taxon_col, keep_tax, sample_cols))) |>
    dplyr::rename(taxon = !!taxon_col) |>
    tidyr::pivot_longer(dplyr::all_of(sample_cols),
                        names_to = "barcode", values_to = "count") |>
    dplyr::mutate(marker    = info$marker,
                  seq_batch = info$seq_batch,
                  rank      = info$rank,
                  count     = as.numeric(count))
  long
}

compile_matrices <- function(files, schema) {
  purrr::map(files, melt_matrix, schema = schema) |> dplyr::bind_rows()
}

# ----------------------------------------------------------------------------
# 3. JOIN: attach metadata to melted counts on BARCODE STRING (per marker).
#    The full barcode string (e.g. "barcode01_97") is GLOBALLY UNIQUE across
#    the whole study (verified), so it — not the batch number — is the join
#    key. Batch numbers in metadata filenames do NOT correspond to batch
#    numbers in matrix filenames, so seq_batch is kept only as provenance,
#    never used to join.
#
#    First we collapse metadata to ONE row per id_sampel, carrying whichever
#    barcode_16s / barcode_its values are recorded for that sample (they may
#    have been entered in different files). Then per marker we join matrix
#    barcodes to that pool on the barcode string alone.
# ----------------------------------------------------------------------------

# Build a combined, de-duplicated sample registry keyed by id_sampel.
# For each sample we take the first non-NA value of each field across all
# metadata files (handles barcodes split across separate per-batch files).
build_sample_registry <- function(meta, schema) {
  first_non_na <- function(x) {
    x <- x[!is.na(x) & x != "NA"]
    if (length(x)) x[1] else NA_character_
  }

  bc_cols  <- vapply(schema$marker_columns, function(m) m$barcode, character(1))
  flag_cols <- vapply(schema$marker_columns, function(m) m$flag,    character(1))
  keep <- c("id_sampel", "stage", "field", "fertilizer", "timepoint",
            "extraction_batch", unname(bc_cols), unname(flag_cols))
  keep <- intersect(keep, names(meta))

  meta |>
    dplyr::select(dplyr::all_of(keep)) |>
    dplyr::group_by(id_sampel) |>
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~ first_non_na(as.character(.x))),
                     .groups = "drop")
}

join_counts_metadata <- function(long_counts, meta, schema) {
  registry <- build_sample_registry(meta, schema)
  markers  <- names(schema$marker_columns)

  joined <- purrr::map_dfr(markers, function(mk) {
    bc_col   <- schema$marker_columns[[mk]]$barcode    # barcode_16s / barcode_its
    flag_col <- schema$marker_columns[[mk]]$flag
    if (!bc_col %in% names(registry)) return(tibble::tibble())

    meta_mk <- registry |>
      dplyr::filter(!is.na(.data[[bc_col]]), .data[[bc_col]] != "NA") |>
      dplyr::transmute(
        id_sampel, stage, field, fertilizer, timepoint, extraction_batch,
        barcode = .data[[bc_col]],
        flag    = if (flag_col %in% names(registry)) .data[[flag_col]] else NA_character_,
        marker  = mk
      )

    cnt_mk <- dplyr::filter(long_counts, marker == mk)

    # join on barcode string only; seq_batch stays as matrix provenance
    dplyr::inner_join(cnt_mk, meta_mk, by = c("marker", "barcode"))
  })

  joined
}

# Kept for _targets.R compatibility. With a globally-unique barcode key there
# is no per-batch prefix fallback needed; this is now a thin wrapper.
join_with_fallback <- function(long_counts, meta, schema) {
  join_counts_metadata(long_counts, meta, schema)
}

# ----------------------------------------------------------------------------
# 4. VALIDATION (STRICT — stop() on any violation)
#    Runs on metadata + the join result. Collects ALL errors, then halts.
# ----------------------------------------------------------------------------
validate_data <- function(meta, joined, long_counts, schema) {
  errs <- character(0)
  add <- function(...) errs <<- c(errs, paste0(...))

  # -- 4a. required columns present after rename
  miss <- setdiff(schema$required_columns, names(meta))
  if (length(miss)) add("Missing required metadata columns: ",
                        paste(miss, collapse = ", "))

  # -- 4b. stage allowed values
  if ("stage" %in% names(meta)) {
    bad <- setdiff(unique(meta$stage), schema$allowed$stage)
    if (length(bad)) add("Invalid stage value(s): ", paste(bad, collapse = ", "),
                         " (allowed: ", paste(schema$allowed$stage, collapse = "/"), ")")
  }

  # -- 4c. timepoint pattern
  if ("timepoint" %in% names(meta)) {
    tp_rx <- schema$allowed$timepoint_pattern
    bad <- unique(meta$timepoint[!stringr::str_detect(meta$timepoint, tp_rx)])
    bad <- bad[!is.na(bad)]
    if (length(bad)) add("Invalid timepoint value(s): ", paste(bad, collapse = ", "))
  }

  # -- 4d. fertilizer pattern, PER STAGE
  if (all(c("stage", "fertilizer") %in% names(meta))) {
    for (st in intersect(unique(meta$stage), names(schema$allowed$fertilizer_pattern))) {
      rx <- schema$allowed$fertilizer_pattern[[st]]
      vals <- meta$fertilizer[meta$stage == st]
      bad <- unique(vals[!stringr::str_detect(vals, rx)])
      bad <- bad[!is.na(bad)]
      if (length(bad)) add("Stage '", st, "' has fertilizer not matching '", rx,
                           "': ", paste(bad, collapse = ", "))
    }
  }

  # -- 4e. blanks in extraction_batch
  if ("extraction_batch" %in% names(meta)) {
    if (any(is.na(meta$extraction_batch)))
      add(sum(is.na(meta$extraction_batch)), " rows have blank extraction_batch")
  }

  # -- 4f. duplicate id_sampel within a marker (same sample joined twice)
  if (nrow(joined) > 0) {
    dup <- joined |>
      dplyr::distinct(marker, id_sampel, barcode) |>
      dplyr::count(marker, id_sampel) |>
      dplyr::filter(n > 1)
    if (nrow(dup))
      add("Duplicate id_sampel within marker (one sample on multiple barcodes): ",
          paste(sprintf("%s/%s", dup$marker, dup$id_sampel), collapse = "; "))
  }

  # -- 4g. Report unmatched matrix barcodes as INFO (not an error).
  #    Inner join already drops barcodes present on only one side. This is
  #    expected and intended behaviour, so we only message, never halt.
  attempted <- long_counts |> dplyr::distinct(marker, barcode)
  matched   <- if (nrow(joined)) dplyr::distinct(joined, marker, barcode) else
                 attempted[0, ]
  unmatched <- dplyr::anti_join(attempted, matched, by = c("marker", "barcode"))
  if (nrow(unmatched)) {
    by_marker <- unmatched |>
      dplyr::group_by(marker) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    message("INFO: dropped unmatched matrix barcodes (no metadata): ",
            paste(sprintf("%s=%d", by_marker$marker, by_marker$n), collapse = ", "))
  }

  if (length(errs)) {
    stop("VALIDATION FAILED (strict mode):\n  - ",
         paste(errs, collapse = "\n  - "), call. = FALSE)
  }
  message("Validation passed: ", nrow(joined), " count rows across ",
          dplyr::n_distinct(joined$marker), " marker(s).")
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# 5. QC: drop samples whose marker flag == 1
# ----------------------------------------------------------------------------
apply_qc <- function(joined, schema) {
  drop_val <- schema$qc$drop_when_flag_equals
  before <- dplyr::n_distinct(joined$id_sampel, joined$marker)

  kept <- joined |>
    dplyr::filter(is.na(flag) | as.character(flag) != as.character(drop_val))

  after <- dplyr::n_distinct(kept$id_sampel, kept$marker)
  message("QC: dropped ", before - after, " flagged sample-units (kept ", after, ").")
  kept
}

# ----------------------------------------------------------------------------
# 6. DROP REPORT: record every barcode that did NOT intersect, both directions
#    - matrix_only : barcode in a matrix but no metadata row  (count dropped)
#    - metadata_only: barcode registered in metadata but absent from matrices
#    Writes a tidy CSV to `out_path` and returns it (as a file target).
# ----------------------------------------------------------------------------
write_drop_report <- function(long_counts, meta, joined, schema,
                              out_path = "results/dropped_barcodes.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  # matched (marker, barcode) pairs
  matched <- if (nrow(joined)) dplyr::distinct(joined, marker, barcode) else
               tibble::tibble(marker = character(), barcode = character())

  # --- side 1: matrix barcodes with no metadata ---
  matrix_bc <- long_counts |> dplyr::distinct(marker, barcode)
  matrix_only <- dplyr::anti_join(matrix_bc, matched, by = c("marker", "barcode")) |>
    dplyr::mutate(side = "matrix_only",
                  reason = "barcode in matrix but no matching metadata row")

  # --- side 2: metadata barcodes with no matrix column ---
  registry <- build_sample_registry(meta, schema)
  meta_long <- purrr::map_dfr(names(schema$marker_columns), function(mk) {
    bc_col <- schema$marker_columns[[mk]]$barcode
    if (!bc_col %in% names(registry)) return(tibble::tibble())
    registry |>
      dplyr::filter(!is.na(.data[[bc_col]]), .data[[bc_col]] != "NA") |>
      dplyr::transmute(marker = mk, barcode = .data[[bc_col]], id_sampel)
  })
  meta_only <- dplyr::anti_join(meta_long, matched, by = c("marker", "barcode")) |>
    dplyr::mutate(side = "metadata_only",
                  reason = "barcode in metadata but absent from matrix")

  report <- dplyr::bind_rows(
    matrix_only |> dplyr::mutate(id_sampel = NA_character_),
    meta_only
  ) |>
    dplyr::select(marker, side, barcode, id_sampel, reason) |>
    dplyr::arrange(marker, side, barcode)

  readr::write_csv(report, out_path)

  message("Drop report: ", nrow(matrix_only), " matrix-only + ",
          nrow(meta_only), " metadata-only barcodes -> ", out_path)
  out_path
}
