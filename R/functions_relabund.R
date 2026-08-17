# =============================================================================
# R/functions_relabund.R — Stage 4: relative abundance stacked bars
#   Source: rarefied genus table. Phylum = agglomerated from genus via lineage.
#   Strict taxonomy already applied upstream; restricted to master samples.
#   Plots: x=fertilizer, facet=timepoint, MEAN relative abundance per group,
#          Top-10 taxa + "Other", "Unclassified" EXCLUDED from plots.
#   CSVs : full un-truncated mean relative abundance per rank.
#   Goal B: per field (subset). Goal D: pooled all fields. Both same layout.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

UNCLASS_RX <- "(?i)^unclassified|^unassigned|^unknown|_incertae_sedis$|^$|^NA$"

# ---- SPECIES INGEST (parallel to genus, restricted to master samples) -------
# Reads wf_<marker>_batch_N_species.tsv, melts, joins to master_samples by the
# barcode string (same join logic as genus). Returns long: marker, id_sampel,
# stage, field, fertilizer, timepoint, species, genus (from lineage), count.
load_species_long <- function(species_files, metadata_all, master_samples,
                              schema) {
  if (length(species_files) == 0) {
    message("load_species_long: no species files found."); return(tibble::tibble())
  }
  message("load_species_long: ", length(species_files), " species files.")
  registry <- build_sample_registry(metadata_all, schema)

  # parse marker/rank DIRECTLY from filename (independent of schema regex,
  # which may be stale). e.g. wf_16s_batch_1_species.tsv -> marker=16S
  parse_species_fn <- function(path) {
    b <- basename(path)
    m <- regmatches(b, regexec(
      "wf_(16s|its)_batch_([0-9]+)_species\\.tsv", b))[[1]]
    if (length(m) < 3) return(NULL)
    list(marker = toupper(m[2]), seq_batch = m[3])
  }
  # non-sample (lineage) columns; fall back to a built-in set if schema lacks it
  ns_cols <- if (!is.null(schema$join$non_sample_columns))
    schema$join$non_sample_columns else
    c("genus","phylum","species","class","order","family","tax","total",
      "superkingdom","kingdom")

  long <- purrr::map_dfr(species_files, function(f) {
    info <- parse_species_fn(f)
    if (is.null(info)) {
      message("  skip (name): ", basename(f)); return(tibble::tibble())
    }
    mat <- readr::read_tsv(f, show_col_types = FALSE,
                           name_repair = "minimal", progress = FALSE)
    label_col <- if ("species" %in% names(mat)) "species" else NA_character_
    if (is.na(label_col)) {
      message("  skip (no 'species' column): ", basename(f),
              " | cols: ", paste(utils::head(names(mat), 6), collapse=", "))
      return(tibble::tibble())
    }
    keep_tax <- intersect(c("tax"), names(mat))
    sample_cols <- setdiff(names(mat), ns_cols)
    sample_cols <- setdiff(sample_cols, c(label_col, "species"))
    mat |>
      dplyr::select(dplyr::all_of(c(label_col, keep_tax, sample_cols))) |>
      dplyr::rename(species = !!label_col) |>
      tidyr::pivot_longer(dplyr::all_of(sample_cols),
                          names_to = "barcode", values_to = "count") |>
      dplyr::mutate(marker = info$marker, count = as.numeric(count))
  })
  if (nrow(long) == 0) {
    message("load_species_long: no species rows after read."); return(tibble::tibble())
  }

  markers <- names(schema$marker_columns)
  joined <- purrr::map_dfr(markers, function(mk) {
    bc <- schema$marker_columns[[mk]]$barcode
    if (!bc %in% names(registry)) return(tibble::tibble())
    meta_mk <- registry |>
      dplyr::filter(!is.na(.data[[bc]]), .data[[bc]] != "NA") |>
      dplyr::transmute(marker = mk, barcode = .data[[bc]], id_sampel)
    dplyr::inner_join(dplyr::filter(long, marker == mk), meta_mk,
                      by = c("marker", "barcode"))
  })
  if (nrow(joined) == 0) {
    message("load_species_long: 0 barcodes matched metadata."); return(tibble::tibble())
  }

  joined <- joined |>
    dplyr::mutate(genus = stringr::str_trim(
                    stringr::str_split_fixed(tax, ";", 8)[, 7]))

  out <- dplyr::inner_join(joined,
    dplyr::select(master_samples, marker, id_sampel, stage, field,
                  fertilizer, timepoint),
    by = c("marker", "id_sampel"))
  message("load_species_long: ", dplyr::n_distinct(out$id_sampel, out$marker),
          " sample-units, ", nrow(out), " rows.")
  out
}

# ---- SPECIES stacked bar for one subset --------------------------------------
# 1) find Top-N genera in THIS subset (by mean rel abundance, genus level)
# 2) keep only species whose parent genus is in that Top-N
# 3) Top-N species + "Other"; mean rel abundance per (fertilizer, timepoint)
species_bar_df <- function(sp_long, genus_long, gp_map, marker, top_n = 15) {
  # Top-N genera for this subset (reuse the genus mean-relabund Top-N logic)
  g_top <- mean_relabund(genus_long, gp_map, marker, "genus",
                         group_cols = c("fertilizer","timepoint"),
                         top_n = top_n, exclude_unclassified = TRUE, full = FALSE)
  top_genera <- setdiff(unique(g_top$taxon), "Other")

  # keep species within the Top-N genera
  d <- dplyr::filter(sp_long, genus %in% top_genera)
  message("  species_bar_df: ", length(top_genera), " top genera, ",
          dplyr::n_distinct(d$species), " species matched (",
          nrow(d), " rows)")
  d$taxon <- d$species
  d$taxon[is.na(d$taxon) | d$taxon == ""] <- "Unclassified"

  # AGGREGATE per (sample, species) FIRST to collapse any duplicate rows
  # (same species via multiple barcodes/lineage variants) -> prevents >100%.
  per_ss <- d |>
    dplyr::group_by(fertilizer, timepoint, id_sampel, taxon) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  # per-sample proportion within the Top-genera species set (sums to 100%)
  per_sample <- per_ss |>
    dplyr::group_by(fertilizer, timepoint, id_sampel) |>
    dplyr::mutate(samp_total = sum(count),
                  rel = ifelse(samp_total > 0, count / samp_total, 0)) |>
    dplyr::ungroup()

  grp_mean <- per_sample |>
    dplyr::group_by(fertilizer, timepoint, taxon) |>
    dplyr::summarise(mean_rel = mean(rel), .groups = "drop")

  grp_mean <- dplyr::filter(grp_mean, !stringr::str_detect(taxon, UNCLASS_RX))
  top_sp <- grp_mean |> dplyr::group_by(taxon) |>
    dplyr::summarise(t = sum(mean_rel), .groups = "drop") |>
    dplyr::slice_max(t, n = top_n) |> dplyr::pull(taxon)

  grp_mean |>
    dplyr::mutate(taxon = ifelse(taxon %in% top_sp, taxon, "Other")) |>
    dplyr::group_by(fertilizer, timepoint, taxon) |>
    dplyr::summarise(mean_rel = sum(mean_rel), .groups = "drop")
}

# genus -> phylum map from the thresholded long table (carries `tax` lineage)
# lineage: superkingdom;kingdom;phylum;class;order;family;genus
build_genus_phylum_map <- function(thresholded_long) {
  thresholded_long |>
    dplyr::filter(rank == "genus") |>
    dplyr::distinct(marker, taxon, tax) |>
    dplyr::mutate(phylum = stringr::str_split_fixed(tax, ";", 7)[, 3]) |>
    dplyr::select(marker, genus = taxon, phylum)
}

# long rel-abundance per sample from the rarefied matrix, joined to metadata
# returns: marker, id_sampel, stage, field, fertilizer, timepoint, genus, count
rarefied_long <- function(norm_tables, master_samples, marker) {
  rar <- norm_tables[[marker]]$rarefied
  if (is.null(rar) || nrow(rar) == 0) return(tibble::tibble())
  df <- as.data.frame(rar)
  df$id_sampel <- rownames(rar)
  long <- tidyr::pivot_longer(df, -id_sampel, names_to = "genus",
                              values_to = "count")
  long$marker <- marker
  dplyr::inner_join(long,
    dplyr::filter(master_samples, marker == !!marker),
    by = c("marker", "id_sampel"))
}

# compute MEAN relative abundance per group, at a given rank, with Top-N + Other
# group_cols defines the bar groups (fertilizer, timepoint [, field])
# returns tidy: group_cols..., taxon, mean_rel (0-1). full = no Top-N truncation.
mean_relabund <- function(long, gp_map, marker, rank, group_cols,
                          top_n = 10, exclude_unclassified = TRUE,
                          full = FALSE) {
  d <- long
  if (rank == "phylum") {
    d <- dplyr::left_join(d, dplyr::filter(gp_map, marker == !!marker),
                          by = c("marker", "genus"))
    d$taxon <- d$phylum
  } else {
    d$taxon <- d$genus
  }
  d$taxon[is.na(d$taxon) | d$taxon == ""] <- "Unclassified"

  # per-sample relative abundance, then mean across samples within group
  per_sample <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "id_sampel"))), taxon) |>
    dplyr::summarise(count = sum(count), .groups = "drop_last") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "id_sampel")))) |>
    dplyr::mutate(rel = count / sum(count)) |>
    dplyr::ungroup()

  grp_mean <- per_sample |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), taxon) |>
    dplyr::summarise(mean_rel = mean(rel), .groups = "drop")

  if (full) return(grp_mean)   # un-truncated, keeps Unclassified -> for CSV

  # ---- plotting version: optionally drop Unclassified, then Top-N + Other ---
  if (exclude_unclassified)
    grp_mean <- dplyr::filter(grp_mean,
                  !stringr::str_detect(taxon, UNCLASS_RX))

  top_taxa <- grp_mean |>
    dplyr::group_by(taxon) |>
    dplyr::summarise(tot = sum(mean_rel), .groups = "drop") |>
    dplyr::slice_max(tot, n = top_n) |> dplyr::pull(taxon)

  grp_mean |>
    dplyr::mutate(taxon = ifelse(taxon %in% top_taxa, taxon, "Other")) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), taxon) |>
    dplyr::summarise(mean_rel = sum(mean_rel), .groups = "drop")
}

# stacked bar: x=fertilizer, facet=timepoint, fill=taxon
plot_stacked_bar <- function(plot_df, rank, title, out_path, all_ferts = NULL) {
  if (nrow(plot_df) == 0) return(NA_character_)
  # fixed fertilizer x-axis if provided
  if (!is.null(all_ferts))
    plot_df$fertilizer <- factor(plot_df$fertilizer, levels = all_ferts)
  # GLOBAL order: rank taxa by total mean abundance across all groups in plot.
  # ggplot stacks the FIRST factor level at the TOP, so to put "Other" on top
  # we make it the first level, followed by taxa in ascending abundance
  # (so the most abundant sits at the BOTTOM of the stack, a common convention).
  ord_desc <- plot_df |> dplyr::group_by(taxon) |>
    dplyr::summarise(t = sum(mean_rel), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(t)) |> dplyr::pull(taxon)
  real <- setdiff(ord_desc, "Other")
  # legend reads most-abundant first; stack: Other on top, then descending down
  legend_levels <- c(if ("Other" %in% ord_desc) "Other", real)
  plot_df$taxon <- factor(plot_df$taxon, levels = legend_levels)

  n_tax <- length(legend_levels)
  pal <- setNames(scales::hue_pal()(length(real)), real)
  if ("Other" %in% legend_levels) pal <- c(Other = "grey75", pal)

  p <- ggplot2::ggplot(plot_df,
        ggplot2::aes(fertilizer, mean_rel, fill = taxon)) +
    ggplot2::geom_col(width = .8, colour = "grey30", linewidth = .1) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal, name = stringr::str_to_title(rank),
                               breaks = legend_levels) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                expand = ggplot2::expansion(c(0, .02))) +
    ggplot2::facet_wrap(~ timepoint, nrow = 1) +
    ggplot2::labs(title = title, x = "Fertilizer",
                  y = "Mean relative abundance") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.text = ggplot2::element_text(size = 8))
  n_tp <- dplyr::n_distinct(plot_df$timepoint)
  ggplot2::ggsave(out_path, p, width = max(7, 2.6 * n_tp), height = 6, dpi = 200,
                  limitsize = FALSE)
  out_path
}

# ============================================================================
# MAIN: relabund bars into the tree (Goal B per field, Goal D pooled)
# ============================================================================
build_relabund_tree <- function(thresholded_long, norm_tables, master_samples,
                                 root = "Results", top_ns = c(10, 15)) {
  gp_map    <- build_genus_phylum_map(thresholded_long)
  universes <- master_samples |> dplyr::distinct(marker, stage) |>
    dplyr::filter(!is.na(stage)) |> dplyr::arrange(marker, stage)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)

    long <- rarefied_long(norm_tables, master_samples, mk)
    long <- dplyr::filter(long, stage == st)
    if (nrow(long) == 0) next
    message("=== RELABUND universe ", ucode, " ===")
    fields <- sort(unique(long$field))
    all_ferts <- fert_rows_for_stage(st)   # PH1-6 (TM) or PHN01-13 (Nursery)

    for (rank in c("phylum", "genus")) {

      ## Goal B — per field
      for (fld in fields) {
        lf <- dplyr::filter(long, field == fld)
        if (nrow(lf) == 0) next
        leaf <- file.path(udir, "Goal_B_Intra_Longitudinal", fld)
        dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

        for (tn in top_ns) {
          plt <- mean_relabund(lf, gp_map, mk, rank,
                   group_cols = c("fertilizer","timepoint"),
                   top_n = tn, exclude_unclassified = TRUE, full = FALSE)
          w <- plot_stacked_bar(plt, rank,
                 title = paste0(ucode," — Goal B — ",fld," — ",rank," (Top ",tn,")"),
                 out_path = file.path(leaf,
                   paste0("stackbar_", rank, "_top", tn, "_", fld, ".png")),
                 all_ferts = all_ferts)
          if (!is.na(w)) written <- c(written, w)
        }
        full <- mean_relabund(lf, gp_map, mk, rank,
                  group_cols = c("fertilizer","timepoint"), full = TRUE)
        readr::write_csv(full, file.path(leaf,
                  paste0("relabund_", rank, "_full_", fld, ".csv")))
      }

      ## Goal D — pooled all fields
      leafD <- file.path(udir, "Goal_D_Cross_Longitudinal")
      dir.create(leafD, recursive = TRUE, showWarnings = FALSE)
      for (tn in top_ns) {
        pltD <- mean_relabund(long, gp_map, mk, rank,
                  group_cols = c("fertilizer","timepoint"),
                  top_n = tn, exclude_unclassified = TRUE, full = FALSE)
        wD <- plot_stacked_bar(pltD, rank,
               title = paste0(ucode," — Goal D (pooled) — ",rank," (Top ",tn,")"),
               out_path = file.path(leafD,
                 paste0("stackbar_", rank, "_top", tn, "_pooled.png")),
               all_ferts = all_ferts)
        if (!is.na(wD)) written <- c(written, wD)
      }
      fullD <- mean_relabund(long, gp_map, mk, rank,
                 group_cols = c("fertilizer","timepoint"), full = TRUE)
      readr::write_csv(fullD, file.path(leafD,
                 paste0("relabund_", rank, "_full_pooled.csv")))
    }
  }
  message("Relabund tree complete: ", length(written), " plots.")
  written
}

# ============================================================================
# SPECIES stacked bars (Goal B per field, Goal D pooled) — Top-15 only
#   Top-15 genera computed per-plot; species filtered to those genera; then
#   Top-15 species + "Other". Same samples as genus QC (master list).
# ============================================================================
build_species_tree <- function(thresholded_long, norm_tables, master_samples,
                                species_long, root = "Results", top_n = 15) {
  if (nrow(species_long) == 0) {
    message("No species data ingested; skipping species bars."); return(character(0))
  }
  gp_map <- build_genus_phylum_map(thresholded_long)
  universes <- master_samples |> dplyr::distinct(marker, stage) |>
    dplyr::filter(!is.na(stage)) |> dplyr::arrange(marker, stage)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)

    g_long <- rarefied_long(norm_tables, master_samples, mk)
    g_long <- dplyr::filter(g_long, stage == st)
    sp_u <- dplyr::filter(species_long, marker == mk, stage == st)
    if (nrow(g_long) == 0 || nrow(sp_u) == 0) next
    message("=== SPECIES bars ", ucode, " ===")
    fields <- sort(unique(g_long$field))
    all_ferts <- fert_rows_for_stage(st)

    ## Goal B — per field
    for (fld in fields) {
      gl <- dplyr::filter(g_long, field == fld)
      sl <- dplyr::filter(sp_u,  field == fld)
      if (nrow(sl) == 0) next
      leaf <- file.path(udir, "Goal_B_Intra_Longitudinal", fld)
      dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
      df <- species_bar_df(sl, gl, gp_map, mk, top_n)
      w <- plot_stacked_bar(df, "species",
            title = paste0(ucode," — Goal B — ",fld," — species (Top ",top_n,
                           ", within Top ",top_n," genera)"),
            out_path = file.path(leaf, paste0("stackbar_species_top",top_n,"_",fld,".png")),
            all_ferts = all_ferts)
      if (!is.na(w)) written <- c(written, w)
      readr::write_csv(df, file.path(leaf,
            paste0("relabund_species_top",top_n,"_",fld,".csv")))
    }

    ## Goal D — pooled
    leafD <- file.path(udir, "Goal_D_Cross_Longitudinal")
    dir.create(leafD, recursive = TRUE, showWarnings = FALSE)
    dfD <- species_bar_df(sp_u, g_long, gp_map, mk, top_n)
    wD <- plot_stacked_bar(dfD, "species",
           title = paste0(ucode," — Goal D (pooled) — species (Top ",top_n,
                          ", within Top ",top_n," genera)"),
           out_path = file.path(leafD, paste0("stackbar_species_top",top_n,"_pooled.png")),
           all_ferts = all_ferts)
    if (!is.na(wD)) written <- c(written, wD)
    readr::write_csv(dfD, file.path(leafD,
           paste0("relabund_species_top",top_n,"_pooled.csv")))
  }
  message("Species bar tree complete: ", length(written), " plots.")
  written
}
