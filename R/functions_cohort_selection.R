# =============================================================================
# R/functions_cohort_selection.R — cohort gate after gold QC
#
# Gold QC and recommendation outputs always retain every sample. This gate
# limits only the downstream processed matrices, normalization, dashboards,
# and analyses. The default `selected: all` expands to every supported cohort.
# `combined_kebun` adds overlapping virtual fields to per-kebun analyses only;
# pooled analyses continue using original fields to prevent double-counting.
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

# Parse optional, overlapping virtual kebun definitions. One row represents
# one virtual group in one marker-stage cohort; `members` is a list-column of
# original Kode Kebun values. Groups not belonging to a selected cohort are
# retained in YAML but do not enter the current target graph.
load_combined_kebun <- function(path = "config/analysis_cohorts.yaml",
                                selected_cohorts = validate_analysis_cohorts()) {
  cfg <- yaml::read_yaml(path)
  groups <- cfg$combined_kebun
  if (is.null(groups) || length(groups) == 0) {
    return(tibble::tibble(name = character(), cohort = character(),
                          members = list()))
  }
  if (is.null(names(groups)) || any(!nzchar(names(groups))))
    stop("analysis_cohorts: every combined_kebun entry must have a name.",
         call. = FALSE)
  bad_names <- names(groups)[!grepl("^[A-Za-z0-9_.-]+$", names(groups))]
  if (length(bad_names))
    stop("analysis_cohorts: combined_kebun names may contain only letters, ",
         "numbers, underscore, dot, and hyphen: ",
         paste(bad_names, collapse = ", "), call. = FALSE)

  rows <- purrr::imap_dfr(groups, function(spec, group_name) {
    cohorts <- unique(trimws(as.character(unlist(spec$cohorts,
                                                  use.names = FALSE))))
    members <- unique(trimws(as.character(unlist(spec$members,
                                                  use.names = FALSE))))
    cohorts <- cohorts[nzchar(cohorts)]
    members <- members[nzchar(members)]
    invalid <- setdiff(cohorts, SUPPORTED_ANALYSIS_COHORTS)
    if (length(invalid))
      stop("analysis_cohorts: combined_kebun '", group_name,
           "' has unsupported cohort(s): ", paste(invalid, collapse = ", "),
           call. = FALSE)
    if (length(cohorts) == 0)
      stop("analysis_cohorts: combined_kebun '", group_name,
           "' must specify at least one cohort.", call. = FALSE)
    if (length(members) < 2)
      stop("analysis_cohorts: combined_kebun '", group_name,
           "' must contain at least two original kebun.", call. = FALSE)
    tibble::tibble(name = group_name, cohort = cohorts,
                   members = rep(list(members), length(cohorts)))
  }) |>
    dplyr::filter(cohort %in% selected_cohorts)

  duplicates <- rows |>
    dplyr::count(name, cohort) |>
    dplyr::filter(n > 1)
  if (nrow(duplicates))
    stop("analysis_cohorts: duplicate combined_kebun name/cohort definitions.",
         call. = FALSE)
  if (nrow(rows))
    message("Additional combined kebun: ", paste0(
      rows$name, " [", rows$cohort, ": ",
      purrr::map_chr(rows$members, paste, collapse = "+"), "]",
      collapse = ", "))
  rows
}

# Return original and applicable virtual field definitions for one universe.
# The returned list-column always contains original field names; downstream
# code subsets those samples and then relabels only the analysis copy.
analysis_field_definitions <- function(lookup, marker, stage,
                                       combined_kebun = NULL) {
  original <- sort(unique(stats::na.omit(
    lookup$field[lookup$stage == stage])))
  base <- tibble::tibble(field = original, members = as.list(original),
                         is_combined = FALSE)
  if (is.null(combined_kebun) || nrow(combined_kebun) == 0)
    return(base)

  cohort_code <- paste0(marker, "_", stage)
  extra <- dplyr::filter(combined_kebun, .data$cohort == .env$cohort_code)
  if (nrow(extra) == 0) return(base)
  collision <- intersect(extra$name, original)
  if (length(collision))
    stop("analysis_cohorts: combined_kebun name conflicts with original ",
         "Kode Kebun in ", cohort_code, ": ",
         paste(collision, collapse = ", "),
         call. = FALSE)

  extra <- extra |>
    dplyr::mutate(
      requested_members = members,
      members = purrr::map(members, intersect, y = original),
      missing = purrr::map2_chr(
        requested_members, members,
        function(requested, present) paste(setdiff(requested, present),
                                           collapse = ","))) |>
    dplyr::select(-requested_members)
  missing_rows <- dplyr::filter(extra, nzchar(missing))
  if (nrow(missing_rows))
    warning("analysis_cohorts: combined kebun member(s) absent from ",
            cohort_code,
            ": ", paste0(missing_rows$name, "=", missing_rows$missing,
                           collapse = "; "), call. = FALSE)
  empty <- extra$name[purrr::map_int(extra$members, length) == 0]
  if (length(empty))
    stop("analysis_cohorts: combined kebun has no available members in ",
         cohort_code, ": ", paste(empty, collapse = ", "), call. = FALSE)

  dplyr::bind_rows(
    base,
    dplyr::transmute(extra, field = name, members = members,
                     is_combined = TRUE))
}

analysis_field_subset <- function(data, definition) {
  stopifnot(nrow(definition) == 1, "field" %in% names(data))
  label <- definition$field[[1]]
  members <- as.character(unlist(definition$members[[1]], use.names = FALSE))
  data |>
    dplyr::filter(.data$field %in% members) |>
    dplyr::mutate(field = label)
}

write_analysis_cohorts <- function(selected_cohorts, combined_kebun = NULL,
                                   out_path = "Results/analysis/selected_cohorts.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  selected <- tibble::tibble(type = "selected_cohort",
                             cohort = selected_cohorts,
                             name = selected_cohorts, members = NA_character_)
  groups <- if (is.null(combined_kebun) || nrow(combined_kebun) == 0) {
    tibble::tibble(type = character(), cohort = character(),
                   name = character(), members = character())
  } else {
    combined_kebun |>
      dplyr::transmute(type = "combined_kebun", cohort, name,
                       members = purrr::map_chr(members, paste, collapse = "+"))
  }
  readr::write_csv(dplyr::bind_rows(selected, groups), out_path)
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
