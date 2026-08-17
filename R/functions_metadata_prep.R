# =============================================================================
# R/functions_metadata_prep.R  —  TASK 1: metadata_prep
#   Reads data/source/<date>_metadata.xlsx (row 1 = merged section headers,
#   dropped; row 2 = real column headers; data starts row 3) and produces:
#     BRONZE: data/bronze/bronze_metadata.csv        (raw column selection)
#     SILVER: data/silver/silver_metadata_16s.csv,
#             data/silver/silver_metadata_its.csv    (parsed, QC'd, split)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# ----------------------------------------------------------------------------
# SOURCE — read + combine all dated source xlsx files
# ----------------------------------------------------------------------------
read_source_metadata_one <- function(path) {
  dat <- readxl::read_excel(path, skip = 1, col_types = "text",
                            .name_repair = "minimal")
  dat$.source_file <- basename(path)
  dat
}

compile_source_metadata <- function(files) {
  purrr::map(files, read_source_metadata_one) |> dplyr::bind_rows()
}

# ----------------------------------------------------------------------------
# BRONZE — select the raw columns, no transformation
# ----------------------------------------------------------------------------
BRONZE_META_COLS <- c("No/Kode Sampel", "Nama Sampel", "Jenis Kebun",
                      "Jenis Analisis", "Batch Ekstraksi")

build_bronze_metadata <- function(source_metadata) {
  missing <- setdiff(BRONZE_META_COLS, names(source_metadata))
  if (length(missing))
    stop("bronze_metadata: missing expected column(s): ",
         paste(missing, collapse = ", "))
  out <- dplyr::select(source_metadata, dplyr::all_of(BRONZE_META_COLS))
  message("bronze_metadata: ", nrow(out), " rows.")
  out
}

write_bronze_metadata <- function(bronze, out_path = "data/bronze/bronze_metadata.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(bronze, out_path)
  message("Wrote ", out_path, " (", nrow(bronze), " rows)")
  out_path
}

# ----------------------------------------------------------------------------
# SILVER — parse Batch Ekstraksi, QC, derive columns, split by marker
# ----------------------------------------------------------------------------

# "9 (31/07/2026)" -> batch = 9L, date_ekstraksi = "260731" (YYMMDD)
parse_batch_ekstraksi <- function(x) {
  m <- stringr::str_match(x,
        "^\\s*([0-9]+)\\s*\\(([0-9]{1,2})/([0-9]{1,2})/([0-9]{4})\\)\\s*$")
  batch <- suppressWarnings(as.integer(m[, 2]))
  day <- m[, 3]; mon <- m[, 4]; yr <- m[, 5]
  date_ekstraksi <- ifelse(
    is.na(day) | is.na(mon) | is.na(yr), NA_character_,
    paste0(substr(yr, 3, 4),
           sprintf("%02d", as.integer(mon)),
           sprintf("%02d", as.integer(day))))
  tibble::tibble(batch = batch, date_ekstraksi = date_ekstraksi)
}

# first "-" delimited token of Nama Sampel, e.g. "GSDI-PHN04-T0-06-SSUM" -> "GSDI"
extract_kode_kebun <- function(nama) stringr::str_extract(nama, "^[^-]+")

# "-" delimited token matching PH<n> or PHN<nn>, e.g. "...-PH1-..." -> "PH1"
extract_kode_pupuk <- function(nama) {
  toks <- stringr::str_split(nama, "-")
  purrr::map_chr(toks, function(t) {
    hit <- t[stringr::str_detect(t, "^PHN?[0-9]+$")]
    if (length(hit)) hit[1] else NA_character_
  })
}

# "-" delimited token matching T<n>, e.g. "...-T2-..." -> "T2"
extract_waktu <- function(nama) {
  toks <- stringr::str_split(nama, "-")
  purrr::map_chr(toks, function(t) {
    hit <- t[stringr::str_detect(t, "^T[0-9]+$")]
    if (length(hit)) hit[1] else NA_character_
  })
}

# Step 5 QC gate — STRICT, mirrors validate_data() in functions_ingest.R
validate_silver_metadata <- function(df) {
  errs <- character(0)
  add <- function(...) errs <<- c(errs, paste0(...))
  blank <- function(x) is.na(x) | stringr::str_trim(as.character(x)) == ""

  ks <- df$`No/Kode Sampel`
  if (any(blank(ks))) add("No/Kode Sampel has NA/empty value(s)")
  not_int <- !stringr::str_detect(ks, "^[0-9]+$")
  if (any(not_int & !blank(ks)))
    add("No/Kode Sampel has non-integer value(s): ",
        paste(unique(ks[not_int & !blank(ks)]), collapse = ", "))
  dup <- unique(ks[duplicated(ks)])
  if (length(dup))
    add("No/Kode Sampel has duplicate value(s): ", paste(dup, collapse = ", "))

  ns <- df$`Nama Sampel`
  if (any(blank(ns))) add("Nama Sampel has NA/empty value(s)")
  dupn <- unique(ns[duplicated(ns)])
  if (length(dupn))
    add("Nama Sampel has duplicate value(s): ", paste(dupn, collapse = ", "))

  jk <- df$`Jenis Kebun`
  if (any(blank(jk))) add("Jenis Kebun has NA/empty value(s)")
  bad_jk <- setdiff(unique(jk[!blank(jk)]), c("Nursery", "TM"))
  if (length(bad_jk))
    add("Jenis Kebun has value(s) outside {Nursery, TM}: ",
        paste(bad_jk, collapse = ", "))

  ja <- df$`Jenis Analisis`
  if (any(blank(ja))) add("Jenis Analisis has NA/empty value(s)")
  bad_ja <- setdiff(unique(ja[!blank(ja)]), c("16S", "ITS", "16S + ITS"))
  if (length(bad_ja))
    add("Jenis Analisis has value(s) outside {16S, ITS, 16S + ITS}: ",
        paste(bad_ja, collapse = ", "))

  if (any(is.na(df$batch))) add("Batch Ekstraksi has NA value(s)")
  if (any(blank(df$date_ekstraksi))) add("Date Ekstraksi has NA/empty value(s)")

  if (length(errs))
    stop("SILVER METADATA QC FAILED:\n  - ", paste(errs, collapse = "\n  - "),
         call. = FALSE)
  message("Silver metadata QC passed: ", nrow(df), " rows.")
  invisible(TRUE)
}

build_silver_metadata <- function(bronze) {
  df <- bronze

  # 4. drop rows with blank/NA Batch Ekstraksi; parse into batch + date
  before <- nrow(df)
  df <- dplyr::filter(df, !is.na(`Batch Ekstraksi`) &
                          stringr::str_trim(`Batch Ekstraksi`) != "")
  message("silver_metadata: dropped ", before - nrow(df),
          " row(s) with blank Batch Ekstraksi (kept ", nrow(df), ").")

  parsed <- parse_batch_ekstraksi(df$`Batch Ekstraksi`)
  bad <- is.na(parsed$batch) | is.na(parsed$date_ekstraksi)
  if (any(bad))
    stop("silver_metadata: unparseable Batch Ekstraksi value(s): ",
         paste(unique(df$`Batch Ekstraksi`[bad]), collapse = ", "), call. = FALSE)
  df$batch <- parsed$batch
  df$date_ekstraksi <- parsed$date_ekstraksi
  df$`Batch Ekstraksi` <- NULL   # superseded by batch / date_ekstraksi

  # 6-8. derived columns from Nama Sampel
  df$`Kode Kebun` <- extract_kode_kebun(df$`Nama Sampel`)

  df$`Kode Pupuk` <- extract_kode_pupuk(df$`Nama Sampel`)
  n_pupuk <- sum(is.na(df$`Kode Pupuk`))
  if (n_pupuk > 0)
    message(n_pupuk, " sample(s) don't have Kode Pupuk")

  df$Waktu <- extract_waktu(df$`Nama Sampel`)
  n_waktu <- sum(is.na(df$Waktu))
  if (n_waktu > 0)
    message(n_waktu, " sample(s) don't have Waktu")

  # 9. 16S / ITS presence indicators from Jenis Analisis
  df$`16S` <- as.integer(stringr::str_detect(df$`Jenis Analisis`, "16S"))
  df$ITS   <- as.integer(stringr::str_detect(df$`Jenis Analisis`, "ITS"))

  # 5. QC gate (after parsing, before split)
  validate_silver_metadata(df)

  # 10. split by marker
  list(
    all          = df,
    metadata_16s = dplyr::filter(df, `16S` == 1),
    metadata_its = dplyr::filter(df, ITS == 1)
  )
}

write_silver_metadata <- function(df, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, out_path)
  message("Wrote ", out_path, " (", nrow(df), " rows)")
  out_path
}
