# =============================================================================
# R/functions_gold_relabund.R  —  Synchronized: relative abundance stacked
#   bars, Goal B (per field) + Goal D (pooled), sourced from the RAREFIED
#   SPECIES matrix per project spec ("using rarefied species matrix, genus,
#   and phylum"). Unlike the established pipeline (genus is the base rank,
#   phylum rolls up from genus, species is a separate track with its own
#   normalization), here SPECIES is the base and genus/phylum both roll up
#   from it directly via the lineage columns already carried in
#   gold_processed_matrix_*_species.
#   Reuses plot_stacked_bar() from functions_relabund.R unchanged (generic:
#   taxon/mean_rel/fertilizer/timepoint columns only) and fert_rows_for_stage()
#   / UNCLASS_RX from the established pipeline.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# rarefied species matrix (one universe) -> long format, joined to species'
# own genus/phylum lineage (from the processed species table) and to sample
# metadata (field/fertilizer/waktu) from universe_lookup.
gold_rarefied_species_long <- function(rarefied_universe, processed_species,
                                       universe_lookup, stage) {
  if (is.null(rarefied_universe) || nrow(rarefied_universe) == 0)
    return(tibble::tibble())

  # one row per species (first-wins on any lineage ambiguity) -- mirrors
  # merge_rank_across_batches()'s dplyr::first()-per-taxon convention in
  # functions_matrix_prep.R. Precautionary: current data has zero ambiguous
  # species, but distinct(species, genus, phylum) alone would silently fan
  # out the join below if that ever changes.
  lineage <- processed_species |>
    dplyr::distinct(species, genus, phylum) |>
    dplyr::group_by(species) |>
    dplyr::summarise(genus = dplyr::first(genus),
                     phylum = dplyr::first(phylum), .groups = "drop")

  df <- as.data.frame(rarefied_universe)
  df$`Sample alias` <- rownames(rarefied_universe)
  long <- tidyr::pivot_longer(df, -`Sample alias`, names_to = "species",
                              values_to = "count")

  meta <- universe_lookup[universe_lookup$stage == stage, ]

  long |>
    dplyr::left_join(lineage, by = "species") |>
    dplyr::inner_join(meta, by = "Sample alias")
}

# mean relative abundance at a given rank (species/genus/phylum), Top-N +
# Other. Mirrors the established mean_relabund()'s per-sample-then-mean
# logic, generalized to read taxon directly off the long table's own rank
# column (species/genus/phylum are ALL already present per row via lineage)
# rather than needing a separate genus->phylum join.
mean_relabund_gold <- function(long, rank, top_n = 10, full = FALSE) {
  d <- long
  d$taxon <- d[[rank]]
  # unresolved taxa become their OWN visible "Unknown" bucket -- they used
  # to be silently dropped here, which meant Top-N + Other never accounted
  # for 100% of a sample's reads (bars fell short of 100%, and worse: with
  # a large unresolved fraction gone, "top N by abundance" among what was
  # left could be decided by ties among genuinely low-abundance taxa,
  # producing an effectively arbitrary/alphabetical "top N").
  d$taxon[is.na(d$taxon) | d$taxon == "" |
          stringr::str_detect(d$taxon, UNCLASS_RX)] <- "Unknown"

  per_sample <- d |>
    dplyr::group_by(fertilizer, timepoint = waktu, `Sample alias`, taxon) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop_last") |>
    dplyr::group_by(fertilizer, timepoint, `Sample alias`) |>
    # NOTE: sum(count) is a per-group SCALAR; base ifelse(test, yes, no)
    # returns a result the length of `test`, so ifelse(sum(count)>0, ...)
    # collapsed to ONE value (the first row's ratio) and mutate() recycled
    # it to every row in the group -- every taxon in a sample silently got
    # the same `rel`. Materialize the total as its own column first, then
    # use dplyr::if_else (strictly vectorized, errors loudly instead of
    # silently corrupting data if this pattern is ever broken again).
    dplyr::mutate(.samp_total = sum(count),
                  rel = dplyr::if_else(.samp_total > 0, count / .samp_total, 0)) |>
    dplyr::select(-.samp_total) |>
    dplyr::ungroup()

  grp_mean <- per_sample |>
    dplyr::group_by(fertilizer, timepoint, taxon) |>
    dplyr::summarise(mean_rel = mean(rel), .groups = "drop")

  if (full) return(grp_mean)   # un-truncated, keeps Unknown -> for CSV

  # Top-N among CLASSIFIED taxa only -- Unknown never competes for a slot,
  # it's always its own bucket (see below). with_ties=FALSE: many taxa can
  # tie at mean_rel=0 for a small subset, and slice_max()'s default
  # with_ties=TRUE would return ALL of them, blowing the Top-N cap open.
  classified <- dplyr::filter(grp_mean, taxon != "Unknown")
  top_taxa <- classified |> dplyr::group_by(taxon) |>
    dplyr::summarise(tot = sum(mean_rel), .groups = "drop") |>
    dplyr::slice_max(tot, n = top_n, with_ties = FALSE) |> dplyr::pull(taxon)

  bucketed <- grp_mean |>
    dplyr::mutate(taxon = dplyr::case_when(
      taxon == "Unknown"    ~ "Unknown",
      taxon %in% top_taxa   ~ taxon,
      TRUE                  ~ "Other")) |>
    dplyr::group_by(fertilizer, timepoint, taxon) |>
    dplyr::summarise(mean_rel = sum(mean_rel), .groups = "drop")

  # stretch every (fertilizer, timepoint) bar to sum to exactly 100% --
  # Top-N + Other + Unknown already accounts for the full sample (nothing
  # was dropped above), so this mainly guards against floating-point drift,
  # but guarantees every bar visually reaches the same 100% maximum.
  # same scalar-ifelse() bug as per_sample above -- fixed the same way.
  bucketed |>
    dplyr::group_by(fertilizer, timepoint) |>
    dplyr::mutate(.grp_total = sum(mean_rel),
                  mean_rel = dplyr::if_else(.grp_total > 0,
                                            mean_rel / .grp_total, 0)) |>
    dplyr::select(-.grp_total) |>
    dplyr::ungroup()
}

# ============================================================================
# MAIN: Goal B (per field) + Goal D (pooled), ranks = phylum/genus/species
#   phylum & genus: Top-10 and Top-15 (matches the established pipeline's
#   top_ns default). species: Top-15 only (matches the established species
#   track's own convention).
# ============================================================================
build_gold_relabund_tree <- function(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                     gold_processed_matrix_16s_species,
                                     gold_rarefied_species_its, gold_universe_lookup_its,
                                     gold_processed_matrix_its_species,
                                     style = load_plot_style(), root = "Results",
                                     top_ns = c(10, 15)) {
  written <- character(0)

  run_marker <- function(rarefied_list, lookup, processed_species, marker) {
    for (st in names(rarefied_list)) {
      mat <- rarefied_list[[st]]
      if (is.null(mat) || nrow(mat) == 0) next
      ucode <- paste0(marker, "_", st)
      udir  <- file.path(root, ucode)

      long <- gold_rarefied_species_long(mat, processed_species, lookup, st)
      if (nrow(long) == 0) next
      fields <- sort(unique(long$field))
      all_ferts <- fert_rows_for_stage(st)
      message("=== GOLD RELABUND ", ucode, " (", length(fields), " fields) ===")

      for (rank in c("phylum", "genus", "species")) {
        top_n_this_rank <- if (rank == "species") 15 else top_ns

        ## Goal B — per field
        for (fld in fields) {
          lf <- dplyr::filter(long, field == fld)
          if (nrow(lf) == 0) next
          leaf <- file.path(udir, "Goal_B_Intra_Longitudinal", fld)
          dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

          for (tn in top_n_this_rank) {
            plt <- mean_relabund_gold(lf, rank, top_n = tn)
            w <- plot_stacked_bar(plt, rank,
                   title = paste0(ucode, " — Goal B — ", fld, " — ", rank,
                                  " (Top ", tn, ")"),
                   out_path = file.path(leaf,
                     paste0("stackbar_", rank, "_top", tn, "_", fld, ".png")),
                   all_ferts = all_ferts, style = style)
            if (!is.na(w)) written <<- c(written, w)
          }
          full <- mean_relabund_gold(lf, rank, full = TRUE)
          readr::write_csv(full, file.path(leaf,
            paste0("relabund_", rank, "_full_", fld, ".csv")))
        }

        ## Goal D — pooled all fields
        leafD <- file.path(udir, "Goal_D_Cross_Longitudinal")
        dir.create(leafD, recursive = TRUE, showWarnings = FALSE)
        for (tn in top_n_this_rank) {
          pltD <- mean_relabund_gold(long, rank, top_n = tn)
          wD <- plot_stacked_bar(pltD, rank,
                 title = paste0(ucode, " — Goal D (pooled) — ", rank,
                                " (Top ", tn, ")"),
                 out_path = file.path(leafD,
                   paste0("stackbar_", rank, "_top", tn, "_pooled.png")),
                 all_ferts = all_ferts, style = style)
          if (!is.na(wD)) written <<- c(written, wD)
        }
        fullD <- mean_relabund_gold(long, rank, full = TRUE)
        readr::write_csv(fullD, file.path(leafD,
          paste0("relabund_", rank, "_full_pooled.csv")))
      }
    }
  }
  run_marker(gold_rarefied_species_16s, gold_universe_lookup_16s,
            gold_processed_matrix_16s_species, "16S")
  run_marker(gold_rarefied_species_its, gold_universe_lookup_its,
            gold_processed_matrix_its_species, "ITS")

  message("Gold relabund tree complete: ", length(written), " plots.")
  written
}
