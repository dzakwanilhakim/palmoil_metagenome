# =============================================================================
# R/functions_raw_export.R
#   Aggregate ALL raw batch TSVs -> dated raw CSVs (no QC / filtering):
#     <yymmdd>_phylum_long.csv / _genus_long.csv / _species_long.csv
#     <yymmdd>_<rank>_<marker>.csv  (wide: taxa x samples)
#     <yymmdd>_metadata.csv
#   Self-contained: parses filenames directly (does NOT depend on schema regex).
#   Two ways to run:
#     (a) targets:  tar_target(raw_exports, export_raw_data(schema))
#     (b) standalone: source(); export_raw_data(schema)
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(yaml) })

export_raw_data <- function(schema = NULL,
                            mat_dirs  = c("data/raw/raw_mat_16s",
                                          "data/raw/raw_mat_its"),
                            meta_dirs = c("data/raw/metadata_16s",
                                          "data/raw/metadata_its"),
                            out_dir   = "Results/raw_exports",
                            date_tag  = format(Sys.Date(), "%y%m%d")) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pfx <- function(name) file.path(out_dir, paste0(date_tag, "_", name))
  written <- character(0)

  # taxonomy/lineage columns that are NOT sample barcodes
  non_sample <- c("genus","phylum","species","class","order","family",
                  "tax","total","superkingdom","kingdom")

  # parse marker + rank straight from filename (no schema dependency)
  parse_fn <- function(path) {
    b <- basename(path)
    m <- regmatches(b, regexec(
      "wf_(16s|its)_batch_([0-9]+)_(genus|phylum|species)\\.tsv", b))[[1]]
    if (length(m) < 4) return(NULL)
    list(marker = toupper(m[2]), seq_batch = m[3], rank = m[4])
  }

  # ---- 1. read + melt every matrix file --------------------------------------
  mat_files <- unlist(lapply(mat_dirs, function(d)
    list.files(d, pattern = "\\.tsv$", full.names = TRUE)))
  message("Found ", length(mat_files), " matrix files.")

  read_one <- function(f) {
    info <- parse_fn(f)
    if (is.null(info)) { message("  skip (name): ", basename(f)); return(tibble::tibble()) }
    mat  <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE,
                            name_repair = "minimal")
    label_col <- if (info$rank %in% names(mat)) info$rank else names(mat)[1]
    sample_cols <- setdiff(names(mat), c(non_sample, label_col))
    tax_col <- if ("tax" %in% names(mat)) "tax" else NULL
    mat |>
      dplyr::select(dplyr::all_of(c(label_col, tax_col, sample_cols))) |>
      dplyr::rename(taxon = !!label_col) |>
      tidyr::pivot_longer(dplyr::all_of(sample_cols),
                          names_to = "barcode", values_to = "count") |>
      dplyr::mutate(marker = info$marker, rank = info$rank,
                    count = as.numeric(count))
  }
  long_all <- purrr::map_dfr(mat_files, read_one)
  message("Total long rows: ", nrow(long_all))

  # ---- 2. metadata: read all batch files, canonical rename -------------------
  for (d in meta_dirs)
    message("  meta_dir '", d, "' exists=", dir.exists(d),
            " files=", length(list.files(d, pattern = "\\.tsv$")))
  meta_files <- unlist(lapply(meta_dirs, function(d)
    list.files(d, pattern = "\\.tsv$", full.names = TRUE, ignore.case = TRUE)))
  message("Metadata files found: ", length(meta_files))
  if (length(meta_files) == 0)
    stop("No metadata .tsv files found in: ", paste(meta_dirs, collapse=", "),
         " | working dir = ", getwd())
  meta_raw <- purrr::map_dfr(meta_files, function(f)
    readr::read_tsv(f, show_col_types = FALSE, progress = FALSE,
                    name_repair = "minimal") |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)))

  # canonical rename via schema$column_map if provided, else built-in map
  default_map <- c(
    "No/Kode_Sampel"="row_code","ID_Sampel"="id_sampel","Nama_Sampel"="sample_name",
    "Jenis_Kebun"="stage","Kelompok_Kebun"="field","Kelompok_Pupuk"="fertilizer",
    "Fase_Treatment"="timepoint","Barcode_16S"="barcode_16s","Barcode_ITS"="barcode_its",
    "Batch_Ekstraksi"="extraction_batch","Hasil_Sekuensing_16S"="seq_result_16s",
    "Flag_16S"="flag_16s","Hasil_Sekuensing_ITS"="seq_result_its","Flag_ITS"="flag_its")
  cmap <- if (!is.null(schema) && !is.null(schema$column_map))
    unlist(schema$column_map) else default_map
  rename_vec <- setNames(names(cmap), unname(cmap))
  meta <- meta_raw |> dplyr::rename(dplyr::any_of(rename_vec))

  # ---- 2b. aggregated metadata csv (dated) -----------------------------------
  canon_cols <- c("row_code","id_sampel","sample_name","stage","field",
                  "fertilizer","timepoint","barcode_16s","barcode_its",
                  "extraction_batch","seq_result_16s","flag_16s",
                  "seq_result_its","flag_its")
  meta_agg <- meta |>
    dplyr::select(dplyr::any_of(canon_cols)) |>
    dplyr::distinct()
  mp <- pfx("metadata.csv")
  readr::write_csv(meta_agg, mp); written <- c(written, mp)
  message("wrote ", mp, " (", nrow(meta_agg), " rows, ",
          dplyr::n_distinct(meta_agg$id_sampel), " unique samples)")

  # ---- 3. barcode -> sample map (per marker) ---------------------------------
  bc_map <- purrr::map_dfr(c("16S","ITS"), function(mk) {
    bc <- if (mk == "16S") "barcode_16s" else "barcode_its"
    if (!bc %in% names(meta)) return(tibble::tibble())
    meta |>
      dplyr::filter(!is.na(.data[[bc]]), .data[[bc]] != "NA", .data[[bc]] != "") |>
      dplyr::transmute(marker = mk, barcode = .data[[bc]],
                       id_sampel,
                       stage = dplyr::coalesce(stage, NA_character_),
                       field, fertilizer, timepoint)
  })
  message("Barcode map rows: ", nrow(bc_map))

  joined <- dplyr::inner_join(long_all, bc_map, by = c("marker", "barcode"))
  message("After join: ", nrow(joined), " rows (",
          dplyr::n_distinct(joined$id_sampel), " samples)")

  agg <- joined |>
    dplyr::group_by(marker, rank, taxon, tax, id_sampel,
                    stage, field, fertilizer, timepoint) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  # ---- 4. dated per-rank raw csvs (long + wide) ------------------------------
  for (rk in c("phylum", "genus", "species")) {
    d <- dplyr::filter(agg, rank == rk)
    if (nrow(d) == 0) { message("No data for rank ", rk); next }
    lp <- pfx(paste0(rk, "_long.csv"))
    readr::write_csv(d, lp); written <- c(written, lp)
    for (mk in unique(d$marker)) {
      dm <- dplyr::filter(d, marker == mk)
      wide <- dm |>
        dplyr::select(taxon, tax, id_sampel, count) |>
        tidyr::pivot_wider(names_from = id_sampel, values_from = count,
                           values_fill = 0)
      wp <- pfx(paste0(rk, "_", mk, ".csv"))
      readr::write_csv(wide, wp); written <- c(written, wp)
      message("wrote ", wp, " (", nrow(wide), " taxa x ", ncol(wide) - 2, " samples)")
    }
    message("wrote ", lp, " (", nrow(d), " rows)")
  }

  message("\nDone. Dated raw exports in: ", out_dir)
  written
}
