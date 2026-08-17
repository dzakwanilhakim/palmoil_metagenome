# =============================================================================
# R/functions_matrix_prep.R  —  TASK 2: abundance_matrix_prep
#   Reads data/source/wf_{16s,its}_batch_<N>.html (EPI2ME wf-metagenomics /
#   wf-16s sequencing reports) and produces:
#     BRONZE: data/bronze/bronze_wf_runs_{16s,its}.csv        (per-sample read stats)
#             data/bronze/bronze_matrix_{16s,its}_{rank}.csv  (merged across all
#               batches; rank in superkingdom/phylum/genus/species)
#     SILVER: data/silver/silver_wf_runs_{16s,its}.csv        (deduped by latest run)
#             data/silver/silver_matrix_{16s,its}_{rank}.csv  (synced to survivors)
#
#   Reports are large (10-30MB) self-contained HTML with abundance tables
#   embedded as DataTables. Tables are isolated by string search on their
#   `id="DataTable_..._inner"` anchor (not full-document DOM parsing) and each
#   isolated fragment is parsed individually via rvest — this keeps parsing
#   fast even though the surrounding document is huge.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rvest)
  library(xml2)
})

RANKS_SCRAPED <- c("phylum", "genus", "species")   # native DataTables
RANKS_BRONZE  <- c("superkingdom", "phylum", "genus", "species")
LINEAGE_COLS  <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "total", "tax")

# ----------------------------------------------------------------------------
# LOW-LEVEL HTML extraction
# ----------------------------------------------------------------------------
parse_wf_html_filename <- function(path) {
  b <- basename(path)
  m <- regmatches(b, regexec("wf_(16s|its)_batch_([0-9]+)\\.html", b))[[1]]
  if (length(m) < 3) stop("Unexpected wf html filename: ", b)
  list(marker = toupper(m[2]), seq_batch = as.integer(m[3]))
}

.datatable_ids <- function(html) {
  hits <- stringr::str_match_all(html, 'id="(DataTable_[a-f0-9]+_inner)"')[[1]]
  unique(hits[, 2])
}

# substring spanning [start of id="..."] .. [end of the following </table>]
.datatable_html <- function(html, tid, window = 8e6) {
  idx <- stringr::str_locate(html, stringr::fixed(paste0('id="', tid, '"')))[1, 1]
  if (is.na(idx)) return(NA_character_)
  seg <- substr(html, idx, min(nchar(html), idx + window))
  end_rel <- regexpr("</table>", seg, fixed = TRUE)
  if (end_rel == -1) return(NA_character_)
  substr(seg, 1, end_rel + attr(end_rel, "match.length") - 1)
}

# nearest <p>...</p> caption BEFORE the table — distinguishes
# "Abundance table for the X rank." from "Rarefied abundance table ..."
.datatable_caption <- function(html, tid, lookback = 2000) {
  idx <- stringr::str_locate(html, stringr::fixed(paste0('id="', tid, '"')))[1, 1]
  if (is.na(idx)) return(NA_character_)
  ctx <- substr(html, max(1, idx - lookback), idx)
  caps <- stringr::str_match_all(ctx,
            stringr::regex("<p[^>]*>(.*?)</p>", dotall = TRUE))[[1]]
  if (nrow(caps) == 0) return(NA_character_)
  stringr::str_squish(caps[nrow(caps), 2])
}

.datatable_first_header <- function(tbl_html) {
  m <- stringr::str_match(tbl_html, '<th[^>]*scope="col"[^>]*>([^<]*)</th>')
  if (is.na(m[1, 2])) NA_character_ else stringr::str_trim(m[1, 2])
}

.parse_datatable <- function(tbl_html) {
  node <- xml2::read_html(paste0("<table>", tbl_html, "</table>"))
  rvest::html_table(node, header = TRUE)[[1]]
}

# raw (non-rarefied) abundance table for one rank: phylum / genus / species
.find_abundance_table <- function(html, rank) {
  for (tid in .datatable_ids(html)) {
    tbl_html <- .datatable_html(html, tid)
    if (is.na(tbl_html) || !identical(.datatable_first_header(tbl_html), rank)) next
    cap <- .datatable_caption(html, tid)
    if (!is.na(cap) && stringr::str_detect(cap, "(?i)^rarefied")) next
    return(.parse_datatable(tbl_html))
  }
  stop("No raw abundance table found for rank '", rank, "'")
}

# the "Number of reads" table: Sample alias, Reads, Unclassified|Unmapped, ...
.find_reads_table <- function(html) {
  for (tid in .datatable_ids(html)) {
    tbl_html <- .datatable_html(html, tid)
    if (!is.na(tbl_html) &&
        identical(.datatable_first_header(tbl_html), "Sample alias"))
      return(.parse_datatable(tbl_html))
  }
  stop("No 'Number of reads' table found")
}

# report generation date badge, right after the report title
.find_report_date <- function(html) {
  m <- stringr::str_match(html,
        stringr::regex("Sequencing Report.{0,1000}?(20[0-9]{2}-[0-9]{2}-[0-9]{2})",
                       dotall = TRUE))
  if (is.na(m[1, 2])) NA_character_ else m[1, 2]
}

# ----------------------------------------------------------------------------
# PER-REPORT parse: one wf_<marker>_batch_<N>.html -> reads table + rank matrices
# ----------------------------------------------------------------------------
parse_wf_report <- function(path) {
  info <- parse_wf_html_filename(path)
  message("Parsing ", basename(path), " ...")
  html <- readr::read_file(path)

  date_run <- .find_report_date(html)

  reads <- .find_reads_table(html) |>
    dplyr::mutate(
      `No/Kode Sampel` = stringr::str_extract(`Sample alias`, "(?<=_)[0-9]+$"),
      `date run` = date_run,
      marker = info$marker, seq_batch = info$seq_batch) |>
    dplyr::select(`No/Kode Sampel`, `Sample alias`, `date run`, Reads,
                  `Unclassified|Unmapped`, `Unclassified|Unmapped (%)`,
                  marker, seq_batch)

  matrices <- purrr::map(RANKS_SCRAPED, ~ .find_abundance_table(html, .x)) |>
    stats::setNames(RANKS_SCRAPED)

  list(marker = info$marker, seq_batch = info$seq_batch, date_run = date_run,
       reads = reads, matrices = matrices)
}

# ----------------------------------------------------------------------------
# BRONZE — merge across all batches of one marker
# ----------------------------------------------------------------------------
compile_wf_runs <- function(reports) {
  purrr::map_dfr(reports, "reads")
}

# collapse one batch's phylum-level wide table to a synthetic superkingdom-
# level wide table (superkingdom, <barcodes...>, total) by summing phyla.
.superkingdom_from_phylum <- function(phylum_tbl) {
  sample_cols <- setdiff(names(phylum_tbl),
                         c("phylum", "total", "superkingdom", "kingdom", "tax"))
  out <- phylum_tbl |>
    dplyr::group_by(superkingdom) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), sum), .groups = "drop")
  out$total <- rowSums(out[, sample_cols, drop = FALSE])
  out
}

.wide_to_long <- function(wide_tbl, taxon_col) {
  meta_cols <- setdiff(intersect(LINEAGE_COLS, names(wide_tbl)), taxon_col)
  sample_cols <- setdiff(names(wide_tbl), c(taxon_col, meta_cols))
  wide_tbl |>
    dplyr::rename(taxon = dplyr::all_of(taxon_col)) |>
    tidyr::pivot_longer(dplyr::all_of(sample_cols), names_to = "barcode",
                        values_to = "count") |>
    dplyr::select(taxon, dplyr::all_of(setdiff(meta_cols, "total")), barcode, count)
}

# merge one rank's per-batch wide tables into a single wide matrix
# (taxa x ALL barcodes across ALL batches; gaps filled with 0).
merge_rank_across_batches <- function(reports, rank) {
  long <- purrr::map_dfr(reports, function(r) {
    wide <- if (rank == "superkingdom")
      .superkingdom_from_phylum(r$matrices$phylum) else r$matrices[[rank]]
    .wide_to_long(wide, rank)
  })

  meta_cols <- setdiff(names(long), c("taxon", "barcode", "count"))
  lineage <- long |>
    dplyr::distinct(taxon, dplyr::across(dplyr::all_of(meta_cols))) |>
    dplyr::group_by(taxon) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::first), .groups = "drop")

  wide_out <- long |>
    dplyr::group_by(taxon, barcode) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = barcode, values_from = count, values_fill = 0)

  sample_cols <- setdiff(names(wide_out), "taxon")
  out <- dplyr::left_join(lineage, wide_out, by = "taxon")
  out$total <- rowSums(out[, sample_cols, drop = FALSE])
  names(out)[names(out) == "taxon"] <- rank
  dplyr::relocate(out, dplyr::all_of(rank), dplyr::all_of(sample_cols), total)
}

# add all-zero column(s) for any Sample alias present in wf_runs but absent
# from this rank's merged matrix (e.g. a barcode with no classified reads).
add_missing_sample_columns <- function(mat, reads_all) {
  taxon_col <- names(mat)[1]
  present <- setdiff(names(mat), c(taxon_col, intersect(LINEAGE_COLS, names(mat))))
  expected <- unique(reads_all$`Sample alias`)
  missing <- setdiff(expected, present)
  if (length(missing)) {
    message("add_missing_sample_columns [", taxon_col, "]: adding ", length(missing),
            " all-zero sample column(s) (in wf_runs, absent from matrix): ",
            paste(utils::head(missing, 10), collapse = ", "),
            if (length(missing) > 10) ", ..." else "")
    for (m in missing) mat[[m]] <- 0L
  }
  mat
}

compile_bronze_matrices <- function(reports) {
  reads_all <- compile_wf_runs(reports)
  purrr::map(RANKS_BRONZE, function(rk) {
    m <- merge_rank_across_batches(reports, rk)
    add_missing_sample_columns(m, reads_all)
  }) |> stats::setNames(RANKS_BRONZE)
}

write_csv_out <- function(df, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, out_path)
  message("Wrote ", out_path, " (", nrow(df), " rows x ", ncol(df), " cols)")
  out_path
}

# ----------------------------------------------------------------------------
# SILVER
# ----------------------------------------------------------------------------

# 1. dedup by "No/Kode Sampel": keep the row with the latest "date run";
#    if the latest date is tied between >1 rows, drop the whole group
#    (ambiguous — can't determine which re-run is authoritative).
dedup_wf_runs <- function(bronze_runs) {
  d <- dplyr::mutate(bronze_runs,
                     .date_parsed = suppressWarnings(as.Date(`date run`)))

  resolved <- d |>
    dplyr::group_by(`No/Kode Sampel`) |>
    dplyr::group_modify(function(sub, key) {
      if (nrow(sub) == 1) return(sub)
      mx <- max(sub$.date_parsed, na.rm = TRUE)
      at_max <- sub[!is.na(sub$.date_parsed) & sub$.date_parsed == mx, , drop = FALSE]
      if (nrow(at_max) == 1) return(at_max)
      tibble::tibble()
    }) |>
    dplyr::ungroup() |>
    dplyr::select(-.date_parsed)

  n_in  <- dplyr::n_distinct(bronze_runs$`No/Kode Sampel`)
  n_out <- dplyr::n_distinct(resolved$`No/Kode Sampel`)
  message("dedup_wf_runs: ", nrow(bronze_runs), " rows / ", n_in,
          " sample(s) -> ", nrow(resolved), " rows / ", n_out,
          " sample(s) (", n_in - n_out, " dropped: repeat run with no later date).")
  resolved
}

# Restrict a taxon x sample matrix to a set of surviving sample (barcode)
# ids, recomputing `total` over the surviving columns. Shared by the silver
# sync step below and the gold HIGH-PASS filter in functions_gold_qc.R.
filter_matrix_to_samples <- function(mat, keep_ids, context_label = "") {
  taxon_col <- names(mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(mat))
  sample_cols <- setdiff(names(mat), c(taxon_col, meta_cols))
  keep <- intersect(sample_cols, keep_ids)
  dropped <- setdiff(sample_cols, keep)
  if (length(dropped))
    message("filter_matrix_to_samples [", taxon_col,
            if (nzchar(context_label)) paste0(" | ", context_label) else "",
            "]: dropping ", length(dropped), " sample column(s): ",
            paste(utils::head(dropped, 10), collapse = ", "),
            if (length(dropped) > 10) ", ..." else "")

  out <- mat[, c(taxon_col, setdiff(meta_cols, "total"), keep), drop = FALSE]
  out$total <- rowSums(out[, keep, drop = FALSE])
  out
}

# 2. drop matrix sample columns whose barcode was dropped at wf_runs dedup.
sync_matrix_to_silver <- function(bronze_matrix, silver_runs) {
  filter_matrix_to_samples(bronze_matrix, silver_runs$`Sample alias`,
                           "silver wf_runs dedup")
}
