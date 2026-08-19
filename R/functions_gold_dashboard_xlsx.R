# =============================================================================
# R/functions_gold_dashboard_xlsx.R  —  Dashboard summary workbook
#   3 sheets: 16S_Nursery, 16S_TM, ITS_TM (ITS_Nursery has no data). HIGH PASS
#   samples only -- gold_universe_lookup_{16s,its} is already restricted to
#   HIGH PASS by build_universe_lookup(), so no extra filtering needed here.
#   Each sheet has:
#     1. sample count per Kode Kebun (field), + a Total row
#     2. distinct Kode Pupuk (fertilizer) count per Kode Kebun, + a Total row
#        (Total = sum of each field's own distinct-fertilizer count, i.e. a
#        simple column footer -- NOT the deduplicated count of unique
#        fertilizers across the whole stage, which would be smaller since
#        fertilizers repeat across fields)
#     3. one Kode Pupuk x Waktu replicate-count grid PER Kode Kebun, with a
#        Total row (samples per Waktu, across all Kode Pupuk) and a Total
#        column (samples per Kode Pupuk, across all Waktu) -- answers "how
#        many samples in Kode Pupuk A at T0 for Kode Kebun X?" and "how many
#        samples total in T0 for Kode Kebun X?" directly. Reuses
#        gold_field_count_grid() (functions_gold_counts.R) unchanged -- same
#        fixed PH1-6/PHN01-13 x T0-T4 scaffold as the data_counts PNGs.
#
#   NOTE: openxlsx is referenced by renv.lock only as someone else's
#   dependency, not pinned/installed itself as of writing -- install with
#   install.packages("openxlsx") then renv::snapshot(). Called namespaced
#   (openxlsx::...) rather than library()'d in tar_option_set, same reason
#   as ANCOMBC in functions_gold_da.R: an uninstalled package listed there
#   would break every other target, not just this one.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

.dashboard_summary_tables <- function(lookup, stage) {
  d <- dplyr::filter(lookup, stage == !!stage)

  t1 <- d |>
    dplyr::count(field, name = "Jumlah Sampel") |>
    dplyr::rename(`Kode Kebun` = field) |>
    dplyr::arrange(`Kode Kebun`)
  t1 <- dplyr::bind_rows(t1,
    tibble::tibble(`Kode Kebun` = "Total", `Jumlah Sampel` = sum(t1$`Jumlah Sampel`)))

  t2 <- d |>
    dplyr::group_by(field) |>
    dplyr::summarise(`Jumlah Kode Pupuk` = dplyr::n_distinct(fertilizer, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(`Kode Kebun` = field) |>
    dplyr::arrange(`Kode Kebun`)
  t2 <- dplyr::bind_rows(t2,
    tibble::tibble(`Kode Kebun` = "Total", `Jumlah Kode Pupuk` = sum(t2$`Jumlah Kode Pupuk`)))

  list(samples = t1, fertilizers = t2)
}

# Kode Pupuk x Waktu replicate grid for one field, + Total row and column.
.dashboard_field_grid_with_totals <- function(lookup, stage, field) {
  grid <- gold_field_count_grid(lookup, stage, field)
  names(grid)[1] <- "Kode Pupuk"
  tp_cols <- setdiff(names(grid), "Kode Pupuk")

  grid$Total <- rowSums(grid[, tp_cols, drop = FALSE])

  total_row <- tibble::tibble(`Kode Pupuk` = "Total")
  for (col in tp_cols) total_row[[col]] <- sum(grid[[col]])
  total_row$Total <- sum(grid$Total)

  dplyr::bind_rows(grid, total_row)
}

build_dashboard_xlsx <- function(gold_universe_lookup_16s, gold_universe_lookup_its,
                                 out_path = "Results/dashboard/dashboard_summary.xlsx") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  lookups <- list("16S_Nursery" = list(lookup = gold_universe_lookup_16s, stage = "Nursery"),
                  "16S_TM"      = list(lookup = gold_universe_lookup_16s, stage = "TM"),
                  "ITS_TM"      = list(lookup = gold_universe_lookup_its, stage = "TM"))

  wb <- openxlsx::createWorkbook()
  hdr_style   <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                                       border = "TopBottom")
  total_style <- openxlsx::createStyle(textDecoration = "bold")
  title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12)
  field_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 11,
                                       fgFill = "#F2F2F2")

  for (sheet_name in names(lookups)) {
    lookup <- lookups[[sheet_name]]$lookup
    stage  <- lookups[[sheet_name]]$stage
    tabs   <- .dashboard_summary_tables(lookup, stage)
    t1 <- tabs$samples; t2 <- tabs$fertilizers

    openxlsx::addWorksheet(wb, sheet_name)
    row <- 1

    openxlsx::writeData(wb, sheet_name, "Jumlah Sampel per Kode Kebun",
                        startRow = row, startCol = 1)
    openxlsx::addStyle(wb, sheet_name, title_style, rows = row, cols = 1)
    row <- row + 1
    openxlsx::writeData(wb, sheet_name, t1, startRow = row, startCol = 1,
                        headerStyle = hdr_style)
    openxlsx::addStyle(wb, sheet_name, total_style,
                       rows = row + nrow(t1), cols = 1:2, gridExpand = TRUE)
    row <- row + nrow(t1) + 3

    openxlsx::writeData(wb, sheet_name, "Jumlah Kode Pupuk per Kode Kebun",
                        startRow = row, startCol = 1)
    openxlsx::addStyle(wb, sheet_name, title_style, rows = row, cols = 1)
    row <- row + 1
    openxlsx::writeData(wb, sheet_name, t2, startRow = row, startCol = 1,
                        headerStyle = hdr_style)
    openxlsx::addStyle(wb, sheet_name, total_style,
                       rows = row + nrow(t2), cols = 1:2, gridExpand = TRUE)
    row <- row + nrow(t2) + 3

    openxlsx::writeData(wb, sheet_name,
      "Jumlah Sampel per Kode Kebun x Kode Pupuk x Waktu",
      startRow = row, startCol = 1)
    openxlsx::addStyle(wb, sheet_name, title_style, rows = row, cols = 1)
    row <- row + 2

    fields <- sort(unique(lookup$field[lookup$stage == stage]))
    for (fld in fields) {
      g <- .dashboard_field_grid_with_totals(lookup, stage, fld)
      n_cols <- ncol(g)

      openxlsx::writeData(wb, sheet_name, paste0("Kode Kebun: ", fld),
                          startRow = row, startCol = 1)
      openxlsx::addStyle(wb, sheet_name, field_style, rows = row, cols = 1:n_cols,
                         gridExpand = TRUE)
      row <- row + 1
      openxlsx::writeData(wb, sheet_name, g, startRow = row, startCol = 1,
                          headerStyle = hdr_style)
      openxlsx::addStyle(wb, sheet_name, total_style,
                         rows = row + nrow(g), cols = 1:n_cols, gridExpand = TRUE)
      openxlsx::addStyle(wb, sheet_name, total_style,
                         rows = (row + 1):(row + nrow(g)), cols = n_cols)
      row <- row + nrow(g) + 3
    }

    openxlsx::setColWidths(wb, sheet_name, cols = 1:8, widths = c(18, rep(9, 6), 10))
  }

  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Wrote ", out_path)
  out_path
}
