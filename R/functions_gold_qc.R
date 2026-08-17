# =============================================================================
# R/functions_gold_qc.R  —  TASK 3: metadata_matrix (gold layer)
#   Joins silver_metadata_{16s,its} to silver_wf_runs_{16s,its} on
#   "No/Kode Sampel", derives QC status + recommendation, and synchronizes
#   the silver abundance matrices down to HIGH PASS samples only.
#     data/gold/gold_qc_{16s,its}.csv                          (full audit table)
#     data/gold/gold_matrix_{16s,its}_{rank}.csv                (HIGH PASS only)
#     data/gold/recommendation_{16s,its}.csv                    (not HIGH PASS)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

load_gold_thresholds <- function(path = "config/gold_qc_thresholds.yaml") {
  yaml::read_yaml(path)
}

# ----------------------------------------------------------------------------
# 1-4. join + flag + reads_status + qc_status + recommendation
# ----------------------------------------------------------------------------
build_gold_qc <- function(silver_metadata, silver_wf_runs, cutoff) {
  in_meta <- silver_metadata$`No/Kode Sampel`
  in_wf   <- silver_wf_runs$`No/Kode Sampel`

  joined <- dplyr::full_join(silver_metadata, silver_wf_runs,
                             by = "No/Kode Sampel")

  # flag = 1 unless the sample matched on BOTH sides of the join
  joined$flag <- as.integer(!(joined$`No/Kode Sampel` %in% in_meta &
                              joined$`No/Kode Sampel` %in% in_wf))

  # 3. Mapped = Reads - Unclassified|Unmapped; reads_status: HIGH if Mapped
  #    >= marker cutoff, else LOW (NA if Reads itself is unknown, i.e. no
  #    wf_runs match)
  joined$Mapped <- joined$Reads - joined$`Unclassified|Unmapped`
  joined$reads_status <- dplyr::case_when(
    is.na(joined$Mapped)     ~ NA_character_,
    joined$Mapped >= cutoff  ~ "HIGH",
    TRUE                     ~ "LOW")

  # 4. qc_status
  joined$qc_status <- dplyr::case_when(
    joined$flag == 1                    ~ "NO DATA",
    joined$reads_status == "HIGH"       ~ "HIGH PASS",
    joined$reads_status == "LOW"        ~ "LOW PASS",
    TRUE                                ~ NA_character_)

  # 5. recommendation (only for NOT "HIGH PASS")
  joined$recommendation <- dplyr::case_when(
    joined$qc_status == "HIGH PASS" ~ NA_character_,
    joined$qc_status == "LOW PASS"  ~ "RESEQUENCING",
    joined$qc_status == "NO DATA" & joined$reads_status == "LOW" ~
      "RESEQUENCING;CLARIFY",
    joined$qc_status == "NO DATA" ~ "NO DATA",
    TRUE ~ NA_character_)

  message("build_gold_qc: ", nrow(joined), " rows | qc_status: ",
          paste(names(table(joined$qc_status, useNA = "ifany")),
                table(joined$qc_status, useNA = "ifany"), sep = "=",
                collapse = ", "))
  joined
}

# ----------------------------------------------------------------------------
# gold matrix: silver matrix restricted to HIGH PASS Sample alias values
# ----------------------------------------------------------------------------
gold_matrix_high_pass <- function(silver_matrix, gold_qc) {
  keep_ids <- gold_qc$`Sample alias`[
    !is.na(gold_qc$qc_status) & gold_qc$qc_status == "HIGH PASS"]
  filter_matrix_to_samples(silver_matrix, keep_ids, "gold HIGH PASS filter")
}

# ----------------------------------------------------------------------------
# recommendation table: subset of gold_qc that is NOT HIGH PASS
# ----------------------------------------------------------------------------
build_recommendation <- function(gold_qc) {
  out <- dplyr::filter(gold_qc, is.na(qc_status) | qc_status != "HIGH PASS")
  message("build_recommendation: ", nrow(out), " sample(s) not HIGH PASS.")
  out
}
