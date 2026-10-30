# =============================================================================
# R/functions_gold_data_audit.R — data audit / LOW PASS / data report workbooks
#   Built from the complete bronze -> silver -> gold QC layers (all cohorts,
#   before the analysis cohort gate), written to Results/dashboard/:
#     data_audit.xlsx           metadata <-> wf mismatches (+ possible matches),
#                               repeat sequencing, queued rows, name issues
#     low_pass_hypothesis.xlsx  LOW PASS samples with rule-based hypotheses
#                               and recommendations, by batch and by kebun
#     data_report.xlsx          overview, QC by cohort, sequencing batches,
#                               analysable coverage per kebun x waktu
#   Hypotheses are heuristics from read counts, re-runs, the other marker and
#   batch / kebun failure rates -- leads to check in the lab, not diagnoses.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# rule thresholds for the LOW PASS hypotheses
LOWPASS_RULES <- list(
  failed_reads      = 1000,  # Mapped below this = library/amplification failure
  borderline_frac   = 0.60,  # Mapped >= this fraction of cutoff = borderline
  high_unclassified = 10,    # Unclassified|Unmapped (%) at or above this
  batch_low_rate    = 0.40,  # share of a seq batch below cutoff = weak run
  kebun_low_rate    = 0.30,  # share of a kebun's runs below cutoff = site effect
  min_group_n       = 5      # min runs in a batch / kebun before rating it
)

# one recommendation per primary hypothesis
LOWPASS_RECOMMENDATION <- c(
  "Weak sequencing run (borderline depth)" = "Top-up resequencing of the existing library (better-loaded run)",
  "Borderline depth" = "Top-up resequencing of the existing library",
  "Sample-level problem (repeated failure)" = "Re-extract DNA (inhibitor clean-up / more input soil), then re-amplify and resequence",
  "DNA extract problem" = "Re-extract DNA (inhibitor clean-up / more input soil), then re-amplify and resequence",
  "Marker-specific amplification failure" = "Re-amplify this marker (check PCR / primers), then resequence",
  "Weak sequencing run" = "Resequence in a better-loaded run (check pooling / flow cell yield)",
  "Site effect (soil matrix / low biomass)" = "Re-extract with inhibitor removal or more input soil, then resequence",
  "Library/amplification failure" = "Re-amplify and resequence",
  "Low depth, cause unclear" = "Re-amplify and resequence")

.blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
.base_name <- function(x) sub("-II$", "", x)

# ----------------------------------------------------------------------------
# Shared audit tables (one target; the three workbooks are views of it)
# ----------------------------------------------------------------------------
build_data_audit <- function(bronze_metadata, silver_metadata_all,
                             bronze_wf_runs_16s, bronze_wf_runs_its,
                             silver_wf_runs_16s, silver_wf_runs_its,
                             gold_qc_16s, gold_qc_its, reads_cutoff) {
  markers <- c("16S", "ITS")
  bronze_wf  <- list("16S" = bronze_wf_runs_16s, ITS = bronze_wf_runs_its)
  silver_wf  <- list("16S" = silver_wf_runs_16s, ITS = silver_wf_runs_its)
  gold_qc    <- list("16S" = gold_qc_16s, ITS = gold_qc_its)
  bm <- dplyr::mutate(bronze_metadata, `No/Kode Sampel` = as.character(`No/Kode Sampel`))
  meta <- dplyr::mutate(silver_metadata_all,
                        `No/Kode Sampel` = as.character(`No/Kode Sampel`),
                        base_name = .base_name(`Nama Sampel`))

  # ---- queued: no Batch Ekstraksi (excluded from QC, informational) --------
  queued <- bm |>
    dplyr::filter(.blank(`Batch Ekstraksi`)) |>
    dplyr::mutate(status = "Queued: no extraction batch yet (excluded from QC)")

  # ---- name issues ---------------------------------------------------------
  name_issues <- meta |>
    dplyr::filter(is.na(`Kode Pupuk`) | is.na(Waktu)) |>
    dplyr::transmute(`No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`,
                     `Jenis Analisis`, `Kode Kebun`, `Kode Pupuk`, Waktu,
                     handling = dplyr::case_when(
                       is.na(Waktu) ~ "Excluded from analysis: no Waktu (T<n>) in name",
                       Waktu == "T0" ~ paste0("Included as '", T0_UNTREATED_LABEL,
                                              "' (T0 pre-treatment baseline)"),
                       TRUE ~ "Excluded from analysis: no Kode Pupuk in name"))

  per_marker <- purrr::map(markers, function(m) {
    other <- setdiff(markers, m)
    qc <- gold_qc[[m]] |> dplyr::mutate(`No/Kode Sampel` = as.character(`No/Kode Sampel`))
    qc_other <- gold_qc[[other]] |>
      dplyr::mutate(`No/Kode Sampel` = as.character(`No/Kode Sampel`)) |>
      dplyr::filter(!is.na(`Sample alias`)) |>
      dplyr::transmute(`No/Kode Sampel`,
                       other_marker_run = paste0(`Sample alias`, " (batch ", seq_batch,
                                                 ", ", qc_status, ", ",
                                                 format(Mapped, big.mark = ",", trim = TRUE),
                                                 " mapped)"))
    cutoff <- reads_cutoff[[m]]

    # repeat sequencing (all runs of samples sequenced more than once)
    bw <- bronze_wf[[m]] |> dplyr::mutate(`No/Kode Sampel` = as.character(`No/Kode Sampel`))
    sw <- silver_wf[[m]] |> dplyr::mutate(`No/Kode Sampel` = as.character(`No/Kode Sampel`))
    repeats <- bw |>
      dplyr::group_by(`No/Kode Sampel`) |> dplyr::filter(dplyr::n() > 1) |>
      dplyr::mutate(n_runs = dplyr::n()) |> dplyr::ungroup() |>
      dplyr::mutate(
        Mapped = Reads - `Unclassified|Unmapped`,
        run_vs_cutoff = ifelse(Mapped >= cutoff, "HIGH", "LOW"),
        selected = paste(`Sample alias`, seq_batch, `date run`) %in%
          paste(sw$`Sample alias`, sw$seq_batch, sw$`date run`)) |>
      dplyr::left_join(dplyr::select(meta, `No/Kode Sampel`, `Nama Sampel`,
                                     `Jenis Kebun`),
                       by = "No/Kode Sampel") |>
      dplyr::transmute(marker = m, `No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`,
                       n_runs, `Sample alias`, seq_batch, `date run`, Reads,
                       Mapped, run_vs_cutoff,
                       selected = ifelse(selected, "SELECTED (latest)", "not used")) |>
      dplyr::arrange(as.integer(`No/Kode Sampel`), `date run`)

    # metadata rows with no wf run, + possible matches
    sequenced_batches <- qc$batch[!is.na(qc$`Sample alias`) & !is.na(qc$batch)]
    max_seq_batch <- if (length(sequenced_batches)) max(sequenced_batches) else -Inf
    partners <- meta |>
      dplyr::select(partner_code = `No/Kode Sampel`, partner_name = `Nama Sampel`,
                    base_name) |>
      dplyr::left_join(dplyr::select(qc, partner_code = `No/Kode Sampel`,
                                     partner_qc = qc_status),
                       by = "partner_code")
    meta_no_wf <- qc |>
      dplyr::filter(flag == 1, is.na(`Sample alias`)) |>
      dplyr::mutate(base_name = .base_name(`Nama Sampel`)) |>
      dplyr::left_join(qc_other, by = "No/Kode Sampel") |>
      dplyr::mutate(likely_status = dplyr::case_when(
        batch > max_seq_batch ~ "Likely queued: extraction batch not sequenced yet",
        batch == max_seq_batch ~ "Possibly queued: extraction batch partly sequenced",
        TRUE ~ "Missing: extraction batch already sequenced, check barcode / sample sheet"))
    pm <- meta_no_wf |>
      dplyr::select(`No/Kode Sampel`, base_name) |>
      dplyr::inner_join(partners, by = "base_name", relationship = "many-to-many") |>
      dplyr::filter(partner_code != `No/Kode Sampel`) |>
      dplyr::group_by(`No/Kode Sampel`) |>
      dplyr::summarise(replicate_partner = paste0(
        partner_code, " ", partner_name, " (", dplyr::coalesce(partner_qc, paste("no", m, "run")), ")",
        collapse = "; "), .groups = "drop")
    meta_no_wf <- meta_no_wf |>
      dplyr::left_join(pm, by = "No/Kode Sampel") |>
      dplyr::transmute(marker = m, `No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`,
                       `Jenis Analisis`, `Batch Ekstraksi` = batch, date_ekstraksi,
                       likely_status,
                       `possible match: same code, other marker` = other_marker_run,
                       `possible match: replicate (-II) partner` = replicate_partner) |>
      dplyr::arrange(as.integer(`No/Kode Sampel`))

    # wf runs with no metadata row for this marker, + why
    wf_no_meta <- qc |>
      dplyr::filter(flag == 1, is.na(`Nama Sampel`)) |>
      dplyr::select(`No/Kode Sampel`, `Sample alias`, seq_batch, `date run`, Reads,
                    Mapped) |>
      dplyr::left_join(dplyr::select(bm, `No/Kode Sampel`,
                                     `possible match: Nama Sampel` = `Nama Sampel`,
                                     `possible match: Jenis Analisis` = `Jenis Analisis`,
                                     `possible match: Batch Ekstraksi` = `Batch Ekstraksi`),
                       by = "No/Kode Sampel") |>
      dplyr::mutate(reason = dplyr::case_when(
        is.na(`possible match: Nama Sampel`) ~ "Code not in metadata at all (typo in sample sheet / barcode?)",
        .blank(`possible match: Batch Ekstraksi`) ~ "In metadata but Batch Ekstraksi blank",
        !stringr::str_detect(`possible match: Jenis Analisis`, m) ~
          paste0("In metadata but Jenis Analisis has no ", m),
        TRUE ~ "Other")) |>
      dplyr::mutate(marker = m, .before = 1)

    list(repeats = repeats, meta_no_wf = meta_no_wf, wf_no_meta = wf_no_meta,
         qc = qc, qc_other = qc_other, bronze_wf = bw, cutoff = cutoff)
  }) |> stats::setNames(markers)

  list(markers = markers, per_marker = per_marker, queued = queued,
       name_issues = name_issues, meta = meta, reads_cutoff = reads_cutoff)
}

# ----------------------------------------------------------------------------
# LOW PASS hypotheses
# ----------------------------------------------------------------------------
.low_pass_table <- function(audit) {
  R <- LOWPASS_RULES
  purrr::map_dfr(audit$markers, function(m) {
    pmk <- audit$per_marker[[m]]
    qc <- pmk$qc; cutoff <- pmk$cutoff
    seqd <- dplyr::filter(qc, !is.na(qc_status), qc_status %in% c("HIGH PASS", "LOW PASS"))

    batch_rate <- seqd |>
      dplyr::group_by(seq_batch) |>
      dplyr::summarise(batch_n = dplyr::n(),
                       batch_low_rate = mean(qc_status == "LOW PASS"), .groups = "drop")
    kebun_rate <- seqd |>
      dplyr::group_by(`Jenis Kebun`, `Kode Kebun`) |>
      dplyr::summarise(kebun_n = dplyr::n(),
                       kebun_low_rate = mean(qc_status == "LOW PASS"), .groups = "drop")
    other_status <- audit$per_marker[[setdiff(audit$markers, m)]]$qc |>
      dplyr::filter(!is.na(qc_status)) |>
      dplyr::select(`No/Kode Sampel`, other_marker_qc = qc_status)
    prev_runs <- pmk$repeats |>
      dplyr::filter(selected == "not used") |>
      dplyr::group_by(`No/Kode Sampel`) |>
      dplyr::summarise(previous_runs = paste0("batch ", seq_batch, ": ",
                                              format(Mapped, big.mark = ",", trim = TRUE),
                                              " mapped", collapse = "; "),
                       previous_all_low = all(run_vs_cutoff == "LOW"),
                       .groups = "drop")

    qc |>
      dplyr::filter(qc_status == "LOW PASS") |>
      dplyr::left_join(batch_rate, by = "seq_batch") |>
      dplyr::left_join(kebun_rate, by = c("Jenis Kebun", "Kode Kebun")) |>
      dplyr::left_join(other_status, by = "No/Kode Sampel") |>
      dplyr::left_join(prev_runs, by = "No/Kode Sampel") |>
      dplyr::mutate(
        marker = m, cutoff = cutoff,
        pct_of_cutoff = round(100 * Mapped / cutoff, 1),
        other_marker_qc = dplyr::coalesce(other_marker_qc,
                                          ifelse(`Jenis Kebun` == "Nursery",
                                                 "not requested", "no run")),
        h_failed     = Mapped < R$failed_reads,
        h_borderline = Mapped >= R$borderline_frac * cutoff,
        h_repeat     = !is.na(previous_all_low) & previous_all_low,
        h_both_low   = other_marker_qc == "LOW PASS",
        h_marker     = other_marker_qc == "HIGH PASS",
        h_batch      = batch_n >= R$min_group_n & batch_low_rate >= R$batch_low_rate,
        h_kebun      = kebun_n >= R$min_group_n & kebun_low_rate >= R$kebun_low_rate,
        h_unclass    = `Unclassified|Unmapped (%)` >= R$high_unclassified)
  }) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      hypotheses = paste(c(
        if (h_failed) "H1 library/amplification failure (almost no reads)",
        if (h_repeat) "H2 failed again on re-run -> sample-level problem (low DNA yield / PCR inhibitors)",
        if (h_both_low) "H3 both markers low -> DNA extract quantity/quality",
        if (h_marker) "H4 other marker passed -> marker-specific PCR/primer issue",
        if (h_batch) sprintf("H5 weak run: %.0f%% of seq batch %s below cutoff",
                             100 * batch_low_rate, seq_batch),
        if (h_kebun) sprintf("H6 site effect: %.0f%% of %s %s runs below cutoff (soil matrix / low biomass)",
                             100 * kebun_low_rate, `Kode Kebun`, `Jenis Kebun`),
        if (h_borderline) sprintf("H7 borderline depth (%.0f%% of cutoff)", pct_of_cutoff),
        if (h_unclass) sprintf("H8 high unclassified reads (%.1f%%)", `Unclassified|Unmapped (%)`)),
        collapse = "; "),
      primary_hypothesis = dplyr::case_when(
        h_borderline & h_batch ~ "Weak sequencing run (borderline depth)",
        h_borderline ~ "Borderline depth",
        h_repeat ~ "Sample-level problem (repeated failure)",
        h_both_low ~ "DNA extract problem",
        h_marker & h_failed ~ "Marker-specific amplification failure",
        h_batch ~ "Weak sequencing run",
        h_kebun ~ "Site effect (soil matrix / low biomass)",
        h_failed ~ "Library/amplification failure",
        TRUE ~ "Low depth, cause unclear"),
      recommendation = LOWPASS_RECOMMENDATION[[primary_hypothesis]]) |>
    dplyr::ungroup() |>
    dplyr::transmute(marker, `No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`,
                     `Kode Kebun`, `Kode Pupuk`, Waktu, seq_batch, `date run`,
                     Reads, `Unclassified|Unmapped`, `Unclassified|Unmapped (%)`,
                     Mapped, cutoff, pct_of_cutoff, other_marker_qc,
                     previous_runs, batch_low_rate = round(batch_low_rate, 2),
                     kebun_low_rate = round(kebun_low_rate, 2),
                     primary_hypothesis, hypotheses, recommendation) |>
    dplyr::arrange(marker, Mapped)
}

# ----------------------------------------------------------------------------
# workbook helpers
# ----------------------------------------------------------------------------
.xl_styles <- function() list(
  hdr   = openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                                border = "TopBottom", wrapText = TRUE),
  title = openxlsx::createStyle(textDecoration = "bold", fontSize = 13),
  note  = openxlsx::createStyle(fontColour = "#595959", textDecoration = "italic"),
  bold  = openxlsx::createStyle(textDecoration = "bold"))

# one titled table per sheet (or several, stacked, via `row`)
.xl_table <- function(wb, sheet, df, title = NULL, note = NULL, row = 1,
                      widths = NULL) {
  st <- .xl_styles()
  if (!sheet %in% names(wb)) openxlsx::addWorksheet(wb, sheet)
  if (!is.null(title)) {
    openxlsx::writeData(wb, sheet, title, startRow = row)
    openxlsx::addStyle(wb, sheet, st$title, rows = row, cols = 1)
    row <- row + 1
  }
  if (!is.null(note)) {
    openxlsx::writeData(wb, sheet, note, startRow = row)
    openxlsx::addStyle(wb, sheet, st$note, rows = row, cols = 1)
    row <- row + 1
  }
  df <- as.data.frame(df, check.names = FALSE)
  openxlsx::writeData(wb, sheet, df, startRow = row, headerStyle = st$hdr)
  if (nrow(df) == 0)
    openxlsx::writeData(wb, sheet, "(none)", startRow = row + 1)
  if (row <= 3) openxlsx::freezePane(wb, sheet, firstActiveRow = row + 1)
  if (ncol(df) > 0) {
    w <- widths %||% pmin(60, pmax(10, purrr::map_int(seq_along(df), function(i)
      max(nchar(c(names(df)[i], as.character(df[[i]]))), na.rm = TRUE)) + 2))
    openxlsx::setColWidths(wb, sheet, cols = seq_along(df), widths = w)
  }
  row + max(nrow(df), 1) + 3
}

.xl_readme <- function(wb, lines) {
  openxlsx::addWorksheet(wb, "README")
  openxlsx::writeData(wb, "README", tibble::tibble(` ` = lines), colNames = FALSE)
  openxlsx::addStyle(wb, "README", .xl_styles()$title, rows = 1, cols = 1)
  openxlsx::setColWidths(wb, "README", cols = 1, widths = 120)
}

.xl_save <- function(wb, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Wrote ", out_path)
  out_path
}

# ----------------------------------------------------------------------------
# 1. data_audit.xlsx
# ----------------------------------------------------------------------------
build_data_audit_xlsx <- function(audit,
                                  out_path = "Results/dashboard/data_audit.xlsx") {
  pm <- audit$per_marker
  meta_no_wf <- purrr::map_dfr(pm, "meta_no_wf")
  wf_no_meta <- purrr::map_dfr(pm, "wf_no_meta")
  repeats    <- purrr::map_dfr(pm, "repeats")

  summary <- purrr::map_dfr(audit$markers, function(m) tibble::tibble(
    marker = m,
    `metadata rows needing marker` = nrow(pm[[m]]$qc) - nrow(pm[[m]]$wf_no_meta),
    `in metadata, no wf run` = nrow(pm[[m]]$meta_no_wf),
    `  of which likely/possibly queued` =
      sum(stringr::str_detect(pm[[m]]$meta_no_wf$likely_status, "queued")),
    `  of which with a possible match` = sum(
      !is.na(pm[[m]]$meta_no_wf$`possible match: same code, other marker`) |
        !is.na(pm[[m]]$meta_no_wf$`possible match: replicate (-II) partner`)),
    `in wf, no metadata` = nrow(pm[[m]]$wf_no_meta),
    `samples sequenced >1 time` = dplyr::n_distinct(pm[[m]]$repeats$`No/Kode Sampel`),
    `wf runs not used (older repeat)` = sum(pm[[m]]$repeats$selected == "not used")))

  wb <- openxlsx::createWorkbook()
  .xl_readme(wb, c(
    "Data audit: metadata vs EPI2ME wf runs",
    paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    "Summary          - counts per marker",
    "Metadata_no_wf   - samples in metadata that need this marker but have no wf run.",
    "                   likely_status compares the sample's extraction batch with the latest extraction batch already sequenced.",
    "                   Possible matches: the same code sequenced for the other marker, and the '-II' replicate partner (same name without -II).",
    "Wf_no_metadata   - wf runs whose code has no metadata row for this marker, with the metadata row for that code if it exists.",
    "Repeat_sequencing- every run of samples sequenced more than once; only the latest run (date, then seq batch) is used.",
    "Queued_no_batch  - metadata rows with no Batch Ekstraksi: not yet extracted/sequenced, excluded from QC.",
    paste0("Name_issues      - names without Kode Pupuk or Waktu. T0 samples are kept as '",
           T0_UNTREATED_LABEL, "' (pre-treatment baseline).")))
  .xl_table(wb, "Summary", summary, title = "Data audit summary")
  .xl_table(wb, "Metadata_no_wf", meta_no_wf,
            title = "In metadata but no wf sequencing run")
  .xl_table(wb, "Wf_no_metadata", wf_no_meta,
            title = "wf sequencing run but no metadata row for this marker")
  .xl_table(wb, "Repeat_sequencing", repeats,
            title = "Samples sequenced more than once (latest run selected)")
  .xl_table(wb, "Queued_no_batch", audit$queued,
            title = "Metadata rows with no Batch Ekstraksi (queued)")
  .xl_table(wb, "Name_issues", audit$name_issues,
            title = "Nama Sampel without Kode Pupuk / Waktu")
  .xl_save(wb, out_path)
}

# ----------------------------------------------------------------------------
# 2. low_pass_hypothesis.xlsx
# ----------------------------------------------------------------------------
build_low_pass_xlsx <- function(audit,
                                out_path = "Results/dashboard/low_pass_hypothesis.xlsx") {
  R <- LOWPASS_RULES
  low <- .low_pass_table(audit)

  rate_tbl <- function(group_cols) purrr::map_dfr(audit$markers, function(m) {
    audit$per_marker[[m]]$qc |>
      dplyr::filter(qc_status %in% c("HIGH PASS", "LOW PASS")) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(runs = dplyr::n(), low_pass = sum(qc_status == "LOW PASS"),
                       low_rate = round(low_pass / runs, 2),
                       median_mapped = stats::median(Mapped),
                       min_mapped = min(Mapped), .groups = "drop") |>
      dplyr::mutate(marker = m, .before = 1)
  })
  by_batch <- rate_tbl(c("seq_batch", "date run")) |> dplyr::arrange(marker, seq_batch)
  by_kebun <- rate_tbl(c("Jenis Kebun", "Kode Kebun")) |>
    dplyr::arrange(marker, dplyr::desc(low_rate))
  by_hyp <- low |> dplyr::count(marker, primary_hypothesis, recommendation,
                                name = "samples") |>
    dplyr::arrange(marker, dplyr::desc(samples))

  wb <- openxlsx::createWorkbook()
  .xl_readme(wb, c(
    "LOW PASS samples and rule-based hypotheses",
    paste0("LOW PASS = Mapped (Reads - Unclassified|Unmapped) below the marker cutoff: ",
           paste(names(audit$reads_cutoff), unlist(audit$reads_cutoff), sep = " ",
                 collapse = ", "), " (config/gold_qc_thresholds.yaml)."),
    "Hypotheses are heuristics to guide lab follow-up, not diagnoses. Rules:",
    sprintf("H1 library/amplification failure : Mapped < %s", R$failed_reads),
    "H2 sample-level problem             : an earlier run of the same sample was also below cutoff",
    "H3 DNA extract problem              : the other marker (TM: 16S + ITS) is also LOW PASS",
    "H4 marker-specific PCR issue        : the other marker is HIGH PASS",
    sprintf("H5 weak sequencing run            : >= %.0f%% of the seq batch below cutoff (batch >= %d runs)",
            100 * R$batch_low_rate, R$min_group_n),
    sprintf("H6 site effect                    : >= %.0f%% of the kebun's runs below cutoff (kebun >= %d runs)",
            100 * R$kebun_low_rate, R$min_group_n),
    sprintf("H7 borderline depth               : Mapped >= %.0f%% of cutoff", 100 * R$borderline_frac),
    sprintf("H8 high unclassified              : Unclassified|Unmapped >= %s%%", R$high_unclassified),
    "primary_hypothesis / recommendation pick the most actionable rule in the order: borderline, repeated failure, both markers low, marker-specific, weak run, site, failure.",
    "Thresholds live in LOWPASS_RULES (R/functions_gold_data_audit.R)."))
  .xl_table(wb, "Low_pass", low, title = "LOW PASS samples")
  .xl_table(wb, "By_hypothesis", by_hyp, title = "LOW PASS samples by primary hypothesis")
  .xl_table(wb, "By_batch", by_batch, title = "Below-cutoff rate per sequencing batch (selected runs)")
  .xl_table(wb, "By_kebun", by_kebun, title = "Below-cutoff rate per kebun (selected runs)")
  .xl_save(wb, out_path)
}

# ----------------------------------------------------------------------------
# 3. data_report.xlsx
# ----------------------------------------------------------------------------
build_data_report_xlsx <- function(audit,
                                   out_path = "Results/dashboard/data_report.xlsx") {
  pm <- audit$per_marker
  lookups <- purrr::map(pm, function(x) build_universe_lookup(x$qc))

  overview <- purrr::map_dfr(audit$markers, function(m) {
    qc <- pm[[m]]$qc; lk <- lookups[[m]]
    tibble::tibble(
      marker = m,
      `read cutoff (mapped)` = pm[[m]]$cutoff,
      `metadata samples needing marker` = sum(!is.na(qc$`Nama Sampel`)),
      `queued (no extraction batch)` = sum(stringr::str_detect(
        audit$queued$`Jenis Analisis`, m)),
      `wf batches` = dplyr::n_distinct(pm[[m]]$bronze_wf$seq_batch),
      `wf runs` = nrow(pm[[m]]$bronze_wf),
      `samples sequenced >1 time` = dplyr::n_distinct(pm[[m]]$repeats$`No/Kode Sampel`),
      `HIGH PASS` = sum(qc$qc_status == "HIGH PASS", na.rm = TRUE),
      `LOW PASS` = sum(qc$qc_status == "LOW PASS", na.rm = TRUE),
      `NO DATA (no wf run or no metadata)` = sum(qc$qc_status == "NO DATA", na.rm = TRUE),
      `T0 included as Untreated` = sum(lk$fertilizer == T0_UNTREATED_LABEL),
      `analysable samples` = nrow(lk))
  })

  qc_by_cohort <- purrr::map_dfr(audit$markers, function(m) {
    pm[[m]]$qc |>
      dplyr::filter(!is.na(`Jenis Kebun`)) |>
      dplyr::count(cohort = paste0(m, "_", `Jenis Kebun`), qc_status) |>
      tidyr::pivot_wider(names_from = qc_status, values_from = n, values_fill = 0)
  }) |>
    dplyr::mutate(total = rowSums(dplyr::across(-cohort)),
                  `HIGH PASS %` = round(100 * `HIGH PASS` / total, 1))

  batches <- purrr::map_dfr(audit$markers, function(m) {
    sel <- pm[[m]]$qc |> dplyr::filter(!is.na(`Sample alias`), !is.na(`Nama Sampel`))
    pm[[m]]$bronze_wf |>
      dplyr::mutate(Mapped = Reads - `Unclassified|Unmapped`) |>
      dplyr::group_by(seq_batch, `date run`) |>
      dplyr::summarise(runs = dplyr::n(), median_reads = stats::median(Reads),
                       median_mapped = stats::median(Mapped),
                       min_mapped = min(Mapped), max_mapped = max(Mapped),
                       .groups = "drop") |>
      dplyr::left_join(sel |> dplyr::group_by(seq_batch) |>
                         dplyr::summarise(`selected HIGH` = sum(qc_status == "HIGH PASS"),
                                          `selected LOW` = sum(qc_status == "LOW PASS"),
                                          .groups = "drop"),
                       by = "seq_batch") |>
      dplyr::mutate(marker = m, .before = 1)
  }) |> dplyr::arrange(marker, seq_batch)

  coverage <- purrr::imap_dfr(lookups, function(lk, m) {
    lk |>
      dplyr::count(cohort = paste0(m, "_", stage), `Kode Kebun` = field, waktu) |>
      tidyr::pivot_wider(names_from = waktu, values_from = n, values_fill = 0)
  })
  tcols <- sort(grep("^T[0-9]+$", names(coverage), value = TRUE))
  coverage <- coverage |>
    dplyr::mutate(dplyr::across(dplyr::all_of(tcols), ~ tidyr::replace_na(.x, 0L))) |>
    dplyr::select(cohort, `Kode Kebun`, dplyr::all_of(tcols)) |>
    dplyr::mutate(total = rowSums(dplyr::across(dplyr::all_of(tcols)))) |>
    dplyr::arrange(cohort, `Kode Kebun`)

  coverage_detail <- purrr::imap_dfr(lookups, function(lk, m) {
    lk |>
      dplyr::mutate(replicate = ifelse(grepl("-II$", `Nama Sampel`), "II", "I")) |>
      dplyr::count(cohort = paste0(m, "_", stage), `Kode Kebun` = field,
                   `Kode Pupuk` = fertilizer, Waktu = waktu, replicate) |>
      tidyr::pivot_wider(names_from = replicate, values_from = n, values_fill = 0,
                         names_prefix = "replicate ")
  }) |>
    dplyr::mutate(n = rowSums(dplyr::across(dplyr::starts_with("replicate ")))) |>
    dplyr::arrange(cohort, `Kode Kebun`, `Kode Pupuk`, Waktu)

  notes <- tibble::tibble(Note = c(
    "Sources: data/source/*_metadata.xlsx and wf_{16s,its}_batch_<N>.html, via bronze -> silver -> gold QC.",
    "Join key: No/Kode Sampel (metadata) = numeric suffix of the wf Sample alias (barcodeNN_<code>).",
    "Repeat sequencing: only the latest run per sample is used (date run, then seq batch).",
    "HIGH PASS = mapped reads (Reads - Unclassified|Unmapped) >= marker cutoff; LOW PASS below it; NO DATA = no wf run or no metadata.",
    paste0("Analysable = HIGH PASS with Kode Pupuk; T0 samples without Kode Pupuk are kept as '",
           T0_UNTREATED_LABEL, "' (pre-treatment baseline)."),
    "'-II' samples are biological replicates added in the latest sampling; they are pooled with the original replicate of the same Kode Kebun x Kode Pupuk x Waktu.",
    "Rows with no Batch Ekstraksi are queued (not yet extracted/sequenced) and are excluded from QC.",
    "Extraction/sequencing batch is confounded with timepoint: temporal comparisons are descriptive.",
    "ITS is analysed for TM only (Nursery is 16S only).",
    "Details: data_audit.xlsx (mismatches, repeats) and low_pass_hypothesis.xlsx (LOW PASS follow-up)."))

  wb <- openxlsx::createWorkbook()
  r <- .xl_table(wb, "Overview", overview, title = "Data report: overview per marker",
                 note = paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")))
  .xl_table(wb, "Overview", notes, title = "Notes", row = r)
  .xl_table(wb, "QC_by_cohort", qc_by_cohort, title = "QC status per cohort")
  .xl_table(wb, "Sequencing_batches", batches,
            title = "Sequencing batches (all runs; selected = run used after repeat dedup)")
  .xl_table(wb, "Coverage", coverage,
            title = "Analysable samples per cohort x Kode Kebun x Waktu")
  .xl_table(wb, "Coverage_detail", coverage_detail,
            title = "Analysable samples per Kode Kebun x Kode Pupuk x Waktu (replicate I / II)")
  .xl_save(wb, out_path)
}

# ----------------------------------------------------------------------------
# 4. final_metadata.xlsx — sheet "16S" then "ITS", one row per sample that
#   requests that marker (Jenis Analisis): analysed or not, and why
#   "Analysed" = the sample's run is in the final processed matrix that every
#   analysis reads (HIGH PASS, selected cohort, usable name, enough reads left
#   after the taxonomic filter). Reasons are checked in pipeline order, so each
#   sample gets the first step that removed it.
# ----------------------------------------------------------------------------
.processed_depth_after_taxfilter <- function(gold_matrix_species, marker) {
  tf <- taxonomic_filter_gold(gold_matrix_species, marker)
  ids <- setdiff(names(tf), c(names(tf)[1], intersect(LINEAGE_COLS, names(tf))))
  colSums(tf[, ids, drop = FALSE])
}

build_final_metadata_xlsx <- function(audit, analysis_cohorts,
                                      lookup_16s, lookup_its,
                                      processed_16s_species, processed_its_species,
                                      selected_gold_matrix_16s, selected_gold_matrix_its,
                                      processing_thresholds,
                                      out_path = "Results/dashboard/final_metadata.xlsx") {
  lookups   <- list("16S" = lookup_16s, ITS = lookup_its)
  processed <- list("16S" = processed_16s_species, ITS = processed_its_species)
  selected  <- list("16S" = selected_gold_matrix_16s$species,
                    ITS = selected_gold_matrix_its$species)
  low <- .low_pass_table(audit) |>
    dplyr::select(marker, `No/Kode Sampel`, lp_hypothesis = primary_hypothesis,
                  lp_action = recommendation)

  per_marker <- purrr::map_dfr(audit$markers, function(m) {
    pmk <- audit$per_marker[[m]]
    cutoff <- pmk$cutoff
    min_depth <- processing_thresholds$sample_depth[[m]]
    proc_ids <- setdiff(names(processed[[m]]),
                        c(names(processed[[m]])[1], LINEAGE_COLS))
    depth_tf <- .processed_depth_after_taxfilter(selected[[m]], m)
    lk <- lookups[[m]]

    pmk$qc |>
      dplyr::left_join(dplyr::filter(low, marker == m) |> dplyr::select(-marker),
                       by = "No/Kode Sampel") |>
      dplyr::left_join(dplyr::select(lk, `Sample alias`, fert_analysis = fertilizer),
                       by = "Sample alias") |>
      dplyr::mutate(
        depth_tf = unname(depth_tf[`Sample alias`]),
        cohort_ok = paste0(m, "_", `Jenis Kebun`) %in% analysis_cohorts,
        reason = dplyr::case_when(
          is.na(`Sample alias`) ~ "Extracted, not yet sequenced (no wf run)",
          is.na(`Nama Sampel`) ~ "wf run has no metadata row for this marker",
          qc_status == "LOW PASS" ~ sprintf(
            "LOW PASS: %s mapped reads < %s cutoff. Likely: %s. Action: %s",
            format(Mapped, big.mark = ",", trim = TRUE), format(cutoff, big.mark = ",", trim = TRUE),
            lp_hypothesis, lp_action),
          !cohort_ok ~ "Cohort not selected in config/analysis_cohorts.yaml",
          is.na(Waktu) ~ "Nama Sampel has no Waktu (T<n>)",
          is.na(fert_analysis) ~ "Nama Sampel has no Kode Pupuk (T1+ needs one)",
          !`Sample alias` %in% proc_ids ~ sprintf(
            "Too few target reads after taxonomic filter (%s): %s < %s",
            if (m == "16S") "Bacteria/Archaea" else "Fungi",
            format(dplyr::coalesce(depth_tf, 0), big.mark = ",", trim = TRUE),
            format(min_depth, big.mark = ",", trim = TRUE)),
          TRUE ~ "Analysed"),
        # same checks, short category for filtering / pivoting
        general = dplyr::case_when(
          is.na(`Sample alias`) ~ "Queued",
          is.na(`Nama Sampel`) ~ "No metadata",
          qc_status == "LOW PASS" ~ "LOW PASS",
          !cohort_ok ~ "Cohort not selected",
          is.na(Waktu) | is.na(fert_analysis) ~ "Name issue",
          !`Sample alias` %in% proc_ids ~ "LOW PASS after taxonomic filter",
          TRUE ~ "Analysed"),
        status = dplyr::case_when(general == "Analysed" ~ "Analysed",
                                  general == "Queued" ~ "Queued",
                                  TRUE ~ "Not analysed")) |>
      dplyr::transmute(marker = m, `No/Kode Sampel`, wf_only = is.na(`Nama Sampel`),
                       fert_analysis, status, general, reason, `Sample alias`, seq_batch,
                       `date run`, Reads, Mapped,
                       `reads after taxonomic filter` = depth_tf)
  })

  # base rows: every metadata sample (incl. queued) + wf runs without metadata
  queued <- audit$queued |>
    dplyr::transmute(`No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`, `Jenis Analisis`,
                     batch = NA_integer_, date_ekstraksi = NA_character_,
                     `Kode Kebun` = extract_kode_kebun(`Nama Sampel`),
                     `Kode Pupuk` = extract_kode_pupuk(`Nama Sampel`),
                     Waktu = extract_waktu(`Nama Sampel`), queued = TRUE)
  base <- audit$meta |>
    dplyr::transmute(`No/Kode Sampel`, `Nama Sampel`, `Jenis Kebun`, `Jenis Analisis`,
                     batch, date_ekstraksi, `Kode Kebun`, `Kode Pupuk`, Waktu,
                     queued = FALSE) |>
    dplyr::bind_rows(queued)

  fills <- c("Analysed" = "#E2EFDA", "Not analysed" = "#FCE4D6", "Queued" = "#FFF2CC")
  wb <- openxlsx::createWorkbook()
  for (m in audit$markers) {
    # one sheet per marker: only samples whose Jenis Analisis requests it
    # (16S + ITS samples appear on both sheets, single-marker samples on one)
    pm_m <- dplyr::filter(per_marker, marker == m)
    requested <- base[stringr::str_detect(dplyr::coalesce(base$`Jenis Analisis`, ""), m), ]
    wf_only <- pm_m |>
      dplyr::filter(wf_only, !`No/Kode Sampel` %in% requested$`No/Kode Sampel`) |>
      dplyr::select(`No/Kode Sampel`) |>
      dplyr::mutate(queued = FALSE)
    sheet <- dplyr::bind_rows(requested, wf_only) |>
      dplyr::left_join(dplyr::select(pm_m, -marker, -wf_only), by = "No/Kode Sampel") |>
      dplyr::mutate(
        Status = dplyr::case_when(queued ~ "Queued", !is.na(status) ~ status,
                                  TRUE ~ "Not analysed"),
        `General reason` = dplyr::case_when(
          queued ~ "Queued", !is.na(general) ~ general, TRUE ~ "Unknown"),
        Reason = dplyr::case_when(
          queued ~ "No Batch Ekstraksi yet: not yet extracted",
          !is.na(reason) ~ reason, TRUE ~ "Unknown"),
        `Kode Pupuk (analysis)` = dplyr::coalesce(fert_analysis, `Kode Pupuk`),
        replicate = ifelse(grepl("-II$", `Nama Sampel`), "II", "I")) |>
      dplyr::select(`No/Kode Sampel`, `Nama Sampel`, Status, `General reason`, Reason,
                    `Jenis Kebun`, `Jenis Analisis`, `Kode Kebun`,
                    `Kode Pupuk (analysis)`, Waktu, replicate,
                    `Batch Ekstraksi` = batch, date_ekstraksi, `Sample alias`,
                    seq_batch, `date run`, Reads, Mapped,
                    `reads after taxonomic filter`) |>
      dplyr::arrange(suppressWarnings(as.integer(`No/Kode Sampel`)))

    counts <- dplyr::count(sheet, Status)
    note <- paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | ",
                   paste(counts$Status, counts$n, sep = ": ", collapse = " | "),
                   " | Samples whose Jenis Analisis includes ", m, ".",
                   " Analysed = in the final processed matrix used by all analyses;",
                   " General reason = category; Reason = detail of the first",
                   " pipeline step that removed the sample.")
    .xl_table(wb, m, sheet, title = paste0(m, ": analysed / not analysed"), note = note)
    for (k in names(fills)) {
      rows <- which(sheet$Status == k) + 3
      if (length(rows))
        openxlsx::addStyle(wb, m, openxlsx::createStyle(fgFill = fills[[k]]),
                           rows = rows, cols = 3, stack = TRUE)
    }
    openxlsx::addFilter(wb, m, rows = 3, cols = seq_along(sheet))
  }
  .xl_save(wb, out_path)
}
