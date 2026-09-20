# =============================================================================
# R/functions_cohort_selection.R — cohort gate after gold QC
#
# Gold QC and recommendation outputs always retain every sample. This gate
# limits only the downstream processed matrices, normalization, dashboards,
# and analyses. The default `selected: all` expands to every supported cohort.
# =============================================================================

SUPPORTED_ANALYSIS_COHORTS <- c("16S_TM", "ITS_TM", "16S_Nursery")

validate_analysis_cohorts <- function(selected = "all") {
  if (is.null(selected) || length(selected) == 0)
    selected <- "all"

  selected <- trimws(as.character(unlist(selected, use.names = FALSE)))
  selected <- selected[nzchar(selected)]

  if (length(selected) == 0)
    selected <- "all"

  if ("all" %in% tolower(selected)) {
    if (length(selected) > 1)
      stop("analysis_cohorts: 'all' cannot be combined with named cohorts.",
           call. = FALSE)
    selected <- SUPPORTED_ANALYSIS_COHORTS
  }

  invalid <- setdiff(selected, SUPPORTED_ANALYSIS_COHORTS)
  if (length(invalid))
    stop("analysis_cohorts: unsupported cohort(s): ",
         paste(invalid, collapse = ", "),
         ". Allowed values: all, ",
         paste(SUPPORTED_ANALYSIS_COHORTS, collapse = ", "),
         call. = FALSE)

  selected <- SUPPORTED_ANALYSIS_COHORTS[
    SUPPORTED_ANALYSIS_COHORTS %in% unique(selected)]
  message("Downstream analysis cohorts: ", paste(selected, collapse = ", "))
  selected
}

load_analysis_cohorts <- function(path = "config/analysis_cohorts.yaml") {
  cfg <- yaml::read_yaml(path)
  validate_analysis_cohorts(cfg$selected)
}

write_analysis_cohorts <- function(selected_cohorts,
                                   out_path = "Results/analysis/selected_cohorts.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(data.frame(cohort = selected_cohorts), out_path)
  message("Wrote ", out_path)
  out_path
}

# Restrict one marker's gold-QC table by marker x stage cohort. Rows with no
# resolved stage cannot belong to an analysis cohort and are excluded here;
# they remain present in the complete gold_qc CSV and recommendation table.
subset_gold_qc_cohorts <- function(gold_qc, marker, selected_cohorts) {
  if (!marker %in% c("16S", "ITS"))
    stop("subset_gold_qc_cohorts: marker must be '16S' or 'ITS'.", call. = FALSE)

  cohort <- paste(marker, gold_qc[["Jenis Kebun"]], sep = "_")
  out <- gold_qc[!is.na(gold_qc[["Jenis Kebun"]]) &
                   cohort %in% selected_cohorts, , drop = FALSE]
  message("subset_gold_qc_cohorts [", marker, "]: kept ", nrow(out), "/",
          nrow(gold_qc), " QC row(s).")
  out
}
