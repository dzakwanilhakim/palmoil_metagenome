# =============================================================================
# R/functions_gold_processing.R  —  QC: gold -> processed matrix
#   Synchronizes the gold-tier abundance matrices (already restricted to
#   HIGH PASS samples) with the same taxonomic + threshold filtering the
#   established pipeline applied to data/raw, so downstream analysis
#   (alpha/beta/relabund/ANCOM-BC) has a single, QC'd entry point:
#     data/gold/gold_processed_matrix_{16s,its}_{genus,species}.csv
#   Applied independently at genus and species rank (both are DA/co-occurrence
#   targets downstream); phylum is a view-only rollup of processed genus, not
#   given its own independent threshold pass. Operates on the WIDE gold_matrix
#   shape directly (taxon rows already carry lineage columns — no lineage
#   string to parse, unlike the original long-format Stage 3 pipeline).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

load_gold_processing_thresholds <- function(path = "config/gold_processing_thresholds.yaml") {
  yaml::read_yaml(path)
}

# ----------------------------------------------------------------------------
# Taxonomic filter (genus/species rank; both carry superkingdom/kingdom and
# order/family lineage columns straight from the gold matrix).
#   16S: keep superkingdom in {Bacteria, Archaea}
#   ITS: keep kingdom == Fungi
#   both: drop order == Chloroplast OR family == Mitochondria
# ----------------------------------------------------------------------------
taxonomic_filter_gold <- function(mat, marker) {
  keep_16s <- marker == "16S" & mat$superkingdom %in% c("Bacteria", "Archaea")
  keep_its <- marker == "ITS" & mat$kingdom == "Fungi"

  organelle <- rep(FALSE, nrow(mat))
  if ("order" %in% names(mat))
    organelle <- organelle | stringr::str_detect(mat$order, "(?i)chloroplast")
  if ("family" %in% names(mat))
    organelle <- organelle | stringr::str_detect(mat$family, "(?i)mitochondria")
  organelle[is.na(organelle)] <- FALSE

  keep <- (keep_16s | keep_its) & !organelle
  out <- mat[keep, , drop = FALSE]
  message("  taxonomic_filter_gold [", marker, "]: kept ", nrow(out), "/",
          nrow(mat), " taxa")
  out
}

# ----------------------------------------------------------------------------
# Threshold filter, strict order, on the WIDE gold matrix:
#   3a drop SAMPLE COLUMNS below min_sample_depth
#   3b on surviving samples: drop taxon ROWS with total reads < otu_min_reads
#   3c drop taxon ROWS present in < prevalence_frac of surviving samples
# ----------------------------------------------------------------------------
threshold_filter_gold <- function(mat, min_sample_depth, otu_min_reads,
                                  prevalence_frac) {
  taxon_col <- names(mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(mat))
  sample_cols <- setdiff(names(mat), c(taxon_col, meta_cols))

  # 3a — sample depth
  depth <- colSums(mat[, sample_cols, drop = FALSE])
  keep_samples <- names(depth)[depth >= min_sample_depth]
  message(sprintf("  3a sample depth >= %s : kept %d/%d samples",
                  format(min_sample_depth, big.mark = ","),
                  length(keep_samples), length(sample_cols)))
  mat <- mat[, c(taxon_col, setdiff(meta_cols, "total"), keep_samples),
            drop = FALSE]

  # 3b — OTU global reads (on surviving samples)
  otu_tot <- rowSums(mat[, keep_samples, drop = FALSE])
  mat <- mat[otu_tot >= otu_min_reads, , drop = FALSE]

  # 3c — prevalence (>0 in >= prevalence_frac of surviving samples)
  n_samp <- length(keep_samples)
  min_prev <- ceiling(prevalence_frac * n_samp)
  prev <- rowSums(mat[, keep_samples, drop = FALSE] > 0)
  mat <- mat[prev >= min_prev, , drop = FALSE]

  message(sprintf(
    "  3b/3c OTU reads >= %d & prevalence >= %d/%d (%.0f%%) : kept %d taxa",
    otu_min_reads, min_prev, n_samp, 100 * prevalence_frac, nrow(mat)))

  mat$total <- rowSums(mat[, keep_samples, drop = FALSE])
  mat
}

# ----------------------------------------------------------------------------
# Driver: taxonomic filter -> print thresholds -> threshold filter, for one
# (marker, rank) combination.
# ----------------------------------------------------------------------------
build_gold_processed_matrix <- function(gold_matrix_rank, marker, rank,
                                        thresholds) {
  message("=== gold_processed_matrix [", marker, " / ", rank, "] ===")
  th_depth <- thresholds$sample_depth[[marker]]
  th_otu   <- thresholds$otu[[marker]][[rank]]
  message(sprintf(
    "  thresholds: min_sample_depth=%s, otu_min_reads=%s, prevalence_frac=%s",
    format(th_depth, big.mark = ","), th_otu$min_reads, th_otu$prevalence_frac))

  tf <- taxonomic_filter_gold(gold_matrix_rank, marker)
  threshold_filter_gold(tf, th_depth, th_otu$min_reads, th_otu$prevalence_frac)
}
