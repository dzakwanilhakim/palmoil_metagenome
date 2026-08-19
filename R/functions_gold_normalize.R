# =============================================================================
# R/functions_gold_normalize.R  —  rarefied species matrix, per universe
#   Universe = (marker x stage), matching this project's established "four
#   isolated universes, never compared" convention (16S_TM, 16S_Nursery,
#   ITS_TM, ITS_Nursery). Species is the sole rarefaction anchor: alpha/beta
#   diversity and the PCoA/rarefaction-curve dashboard panels all run on this;
#   genus/phylum relative-abundance bars will later roll UP from it. ANCOM-BC
#   (genus + species DA/co-occurrence) instead uses the unrarefied
#   gold_processed_matrix directly — see functions_gold_processing.R.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# HIGH PASS sample -> stage/field/fertilizer/time lookup, from gold_qc.
# Also drops samples with unresolved fertilizer (NA Kode Pupuk, e.g. a
# malformed name like "PSG-PHN-T0-08-SSUM") HERE, at the shared source, not
# in individual downstream consumers -- every consumer of this lookup
# (rarefaction, alpha, beta, relabund, dashboard, DA) needs a valid
# fertilizer to be meaningfully grouped/adjusted-for, and filtering
# per-consumer risks exactly the inconsistency this fixes: alpha's
# gold_alpha_for_report() used to filter NA fertilizer out while beta's
# .beta_meta() didn't, so the same "HIGH PASS" sample set silently differed
# between the two -- e.g. a lone T0 sample showing up in a Goal B PCoA but
# never in the matching Shannon trajectory for the same field.
build_universe_lookup <- function(gold_qc) {
  d <- gold_qc |>
    dplyr::filter(!is.na(qc_status), qc_status == "HIGH PASS") |>
    dplyr::transmute(
      `Sample alias`,
      `Nama Sampel`,
      stage      = `Jenis Kebun`,
      field      = `Kode Kebun`,
      fertilizer = `Kode Pupuk`,
      waktu      = Waktu)

  n_bad <- sum(is.na(d$fertilizer))
  if (n_bad > 0)
    message("build_universe_lookup: dropping ", n_bad,
            " HIGH PASS sample(s) with unresolved Kode Pupuk (can't be ",
            "grouped by fertilizer in any downstream analysis): ",
            paste(d$`Nama Sampel`[is.na(d$fertilizer)], collapse = ", "))

  dplyr::filter(d, !is.na(fertilizer))
}

# taxon x sample wide table -> sample x taxon numeric matrix
.gold_matrix_to_sample_matrix <- function(mat) {
  taxon_col <- names(mat)[1]
  meta_cols <- intersect(LINEAGE_COLS, names(mat))
  sample_cols <- setdiff(names(mat), c(taxon_col, meta_cols))
  m <- as.matrix(mat[, sample_cols, drop = FALSE])
  rownames(m) <- mat[[taxon_col]]
  storage.mode(m) <- "double"
  t(m)   # samples x taxa
}

# rarefy one (marker's) processed species matrix, independently per stage
# (universe), to that universe's own min surviving sample depth. Returns a
# named list keyed by stage ("TM","Nursery"), each a sample x species matrix.
rarefy_species_by_universe <- function(processed_species, universe_lookup,
                                       seed = 42) {
  sm <- .gold_matrix_to_sample_matrix(processed_species)
  stages <- sort(unique(universe_lookup$stage))

  set.seed(seed)
  purrr::map(stages, function(st) {
    ids <- intersect(
      universe_lookup$`Sample alias`[universe_lookup$stage == st],
      rownames(sm))
    sub <- sm[ids, , drop = FALSE]
    depths <- rowSums(sub)
    keep <- rownames(sub)[depths > 0]
    sub <- sub[keep, , drop = FALSE]
    if (nrow(sub) == 0) {
      message("  rarefy_species_by_universe [", st, "]: 0 samples with depth > 0, skipped")
      return(NULL)
    }
    rdepth <- min(rowSums(sub))
    sub_int <- sub; storage.mode(sub_int) <- "integer"
    rar <- vegan::rrarefy(sub_int, sample = rdepth)
    nz <- colSums(rar) > 0
    rar <- rar[, nz, drop = FALSE]
    message(sprintf("  rarefy_species_by_universe [%s]: depth=%s | samples=%d | species=%d",
                    st, format(rdepth, big.mark = ","), nrow(rar), ncol(rar)))
    rar
  }) |> stats::setNames(stages)
}
