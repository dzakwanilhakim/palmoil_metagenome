# =============================================================================
# R/functions_gold_beta.R  —  Analysis output: beta diversity
#   Bray-Curtis only, from the rarefied species matrix, per universe.
#   Per project spec: "Beta diversity using only rarefied species matrix with
#   output only Bray Curtis ordination pcoa, dendrogram (color based on
#   time), and PERMANOVA."
#   PCoA reuses plot_pcoa_universe() from functions_gold_dashboard.R (same
#   plot the Dashboard already produces) into a separate Results/analysis/
#   copy, rather than re-deriving it.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
})

ANALYSIS_DIR <- "Results/analysis"

# ---- dendrogram (Ward.D2), leaves colored by Waktu (time) ------------------
# Mirrors the established pipeline's functions_beta.R::plot_dendrogram() bug
# fix (leaflab="none" + colored mtext labels), grouped by Waktu instead of
# fertilizer.
plot_beta_dendrogram <- function(rarefied_mat, universe_lookup, marker, stage,
                                 style = load_plot_style(),
                                 out_dir = ANALYSIS_DIR, fname_suffix = "") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  sfx <- if (nzchar(fname_suffix)) paste0("_", fname_suffix) else ""
  out_path <- file.path(out_dir,
    paste0("beta_dendro_", tolower(marker), "_", tolower(stage), sfx, ".png"))
  if (is.null(rarefied_mat) || nrow(rarefied_mat) < 3) {
    message("  dendrogram skipped for ", marker, "/", stage, ": <3 samples")
    return(NA_character_)
  }

  d <- vegan::vegdist(rarefied_mat, method = "bray")
  hc <- stats::hclust(d, method = "ward.D2")
  lab_ids <- hc$labels[hc$order]

  meta <- universe_lookup[universe_lookup$stage == stage, ]
  waktu <- meta$waktu[match(lab_ids, meta$`Sample alias`)]
  # display label = Nama Sampel (sample name), not the barcode; lab_ids
  # itself stays barcode-based since that's what the distance matrix and
  # hclust object are keyed on internally.
  nama <- meta$`Nama Sampel`[match(lab_ids, meta$`Sample alias`)]
  nama[is.na(nama)] <- lab_ids[is.na(nama)]   # fallback if no metadata match
  # global time palette (waktu_palette(), functions_plot_style.R) -- same
  # colour scale used everywhere else in the pipeline, not re-derived per plot.
  pal <- waktu_palette()
  cols <- pal[as.character(waktu)]; cols[is.na(cols)] <- "grey50"

  bs <- style$base_size %||% 14
  grDevices::png(out_path, width = max(1400, 34 * length(lab_ids)),
                 height = 1000, res = 150)
  op <- graphics::par(mar = c(16, 4, 3, 1), cex.main = bs / 11, font.main = 2)
  dend <- stats::as.dendrogram(hc)
  plot(dend,
       main = paste0(marker, "_", stage, " — Bray-Curtis dendrogram (Ward.D2)"),
       ylab = "Ward.D2 distance", xlab = "", leaflab = "none")
  graphics::mtext(side = 1, at = seq_along(lab_ids), text = nama,
                  col = cols, las = 2, line = .5, cex = bs / 26, font = 2)
  present <- intersect(names(pal), sort(unique(waktu)))
  graphics::legend("topright", legend = present, fill = pal[present], bty = "n",
                   cex = bs / 16, title = "Waktu")
  graphics::par(op); grDevices::dev.off()
  message("  Wrote ", out_path)
  out_path
}

# ---- PERMANOVA (Bray-Curtis ~ waktu + fertilizer, strata = field) ----------
run_beta_permanova <- function(rarefied_mat, universe_lookup, marker, stage,
                               formula_rhs = "waktu + fertilizer",
                               strata_var = "field", min_n = 2) {
  if (is.null(rarefied_mat) || nrow(rarefied_mat) < 3)
    return(tibble::tibble(marker = marker, stage = stage, term = NA,
                          note = "skipped: <3 samples"))

  meta <- universe_lookup[universe_lookup$stage == stage, ]
  meta <- meta[match(rownames(rarefied_mat), meta$`Sample alias`), ]

  factors <- trimws(strsplit(formula_rhs, "\\+")[[1]])
  for (f in factors) {
    tab <- table(meta[[f]])
    if (length(tab) < 2 || all(tab < min_n))
      return(tibble::tibble(marker = marker, stage = stage, term = f,
             note = sprintf("skipped: factor '%s' lacks replication (n>=%d, >=2 levels)",
                            f, min_n)))
  }

  d <- vegan::vegdist(rarefied_mat, method = "bray")
  fml <- stats::as.formula(paste("d ~", formula_rhs))
  set.seed(42)
  ad <- tryCatch({
    if (!is.null(strata_var) && strata_var %in% names(meta) &&
        dplyr::n_distinct(meta[[strata_var]]) > 1) {
      perm <- permute::how(nperm = 999)
      permute::setBlocks(perm) <- meta[[strata_var]]
      vegan::adonis2(fml, data = meta, permutations = perm, by = "terms")
    } else {
      vegan::adonis2(fml, data = meta, permutations = 999, by = "terms")
    }
  }, error = function(e) NULL)
  if (is.null(ad))
    return(tibble::tibble(marker = marker, stage = stage, term = NA,
                          note = "adonis2 failed"))

  tibble::tibble(marker = marker, stage = stage, term = rownames(ad),
                 Df = ad$Df, R2 = ad$R2, F = ad$F, p = ad$`Pr(>F)`,
                 note = NA_character_)
}

# ---- driver: PCoA + dendrogram + PERMANOVA, all universes ------------------
build_beta_analysis <- function(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                gold_rarefied_species_its, gold_universe_lookup_its,
                                style = load_plot_style()) {
  run_marker <- function(rarefied_list, lookup, marker) {
    purrr::map(names(rarefied_list), function(st) {
      mat <- rarefied_list[[st]]
      list(
        pcoa   = plot_pcoa_universe(mat, lookup, marker, st, style,
                                    out_dir = ANALYSIS_DIR),
        dendro = plot_beta_dendrogram(mat, lookup, marker, st, style),
        perm   = run_beta_permanova(mat, lookup, marker, st))
    })
  }
  all_res <- c(run_marker(gold_rarefied_species_16s, gold_universe_lookup_16s, "16S"),
              run_marker(gold_rarefied_species_its, gold_universe_lookup_its, "ITS"))

  written <- c(purrr::map_chr(all_res, "pcoa"), purrr::map_chr(all_res, "dendro"))
  written <- written[!is.na(written)]

  perm_tbl <- dplyr::bind_rows(purrr::map(all_res, "perm"))
  perm_path <- file.path(ANALYSIS_DIR, "beta_permanova.csv")
  dir.create(dirname(perm_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(perm_tbl, perm_path)
  message("Wrote ", perm_path)
  print(perm_tbl)

  list(plots = written, permanova_csv = perm_path)
}

# ============================================================================
# Goal B/D beta breakdown (synchronized from established functions_beta.R)
#   Goal B: per field, within-field ordination + dendrogram + PERMANOVA ~waktu
#   Goal D: pooled (all fields), faceted-by-field ordination + dendrogram +
#           PERMANOVA ~ waktu + fertilizer, strata = field
#   Reuses plot_ordination() from functions_beta.R unchanged (generic: beta
#   object + meta + color/shape/facet vars); compute_beta_gold() below
#   produces that same $dist/$ord/$method/$ordination_type/$var_explained
#   shape from a plain rarefied matrix, sidestepping the established
#   compute_beta()'s dependency on the old norm_tables[[marker]] structure.
# ============================================================================

compute_beta_gold <- function(rarefied_mat) {
  if (is.null(rarefied_mat) || nrow(rarefied_mat) < 3) return(NULL)
  d <- vegan::vegdist(rarefied_mat, method = "bray")
  pcoa <- stats::cmdscale(d, k = 2, eig = TRUE)
  scores <- as.data.frame(pcoa$points); names(scores) <- c("Axis1", "Axis2")
  scores$id_sampel <- rownames(rarefied_mat)
  eig <- pcoa$eig[pcoa$eig > 0]
  ve  <- (pcoa$eig / sum(eig))[1:2]
  list(dist = d, ord = scores, method = "bray_curtis",
       ordination_type = "PCoA", var_explained = ve)
}

# universe_lookup subset reshaped for plot_ordination(), which expects
# id_sampel/timepoint column names (established pipeline convention).
.beta_meta <- function(universe_lookup, stage) {
  universe_lookup |>
    dplyr::filter(stage == !!stage) |>
    dplyr::rename(id_sampel = `Sample alias`, timepoint = waktu)
}

build_gold_beta_goal_tree <- function(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                      gold_rarefied_species_its, gold_universe_lookup_its,
                                      style = load_plot_style(), root = "Results") {
  written <- character(0); perm_rows <- list()

  run_marker <- function(rarefied_list, lookup, marker) {
    for (st in names(rarefied_list)) {
      mat <- rarefied_list[[st]]
      if (is.null(mat) || nrow(mat) == 0) next
      ucode <- paste0(marker, "_", st)
      udir  <- file.path(root, ucode)
      meta  <- .beta_meta(lookup, st)

      # global time palette (waktu_palette()) -- ggplot's scale_colour_manual
      # only draws legend keys for values actually present in the data, so
      # passing the full T0..T4 palette here is safe (no unused-level clutter).
      time_pal <- waktu_palette()
      fields <- sort(unique(meta$field[meta$id_sampel %in% rownames(mat)]))
      message("=== GOLD BETA GOALS ", ucode, " (", length(fields), " fields) ===")

      ## Goal B — per field
      for (fld in fields) {
        ids <- meta$id_sampel[meta$field == fld]
        sub_mat <- mat[intersect(ids, rownames(mat)), , drop = FALSE]
        if (nrow(sub_mat) < 3) next
        leaf <- file.path(udir, GOAL_DIR["B"], fld)
        dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

        beta <- compute_beta_gold(sub_mat)
        o <- plot_ordination(beta, meta, color_var = "timepoint", shape_var = NULL,
              facet_var = NULL, pal = time_pal,
              title = paste0(ucode, " — Goal B — ", fld, " — Bray-Curtis PCoA"),
              out_path = file.path(leaf, "beta_ord_bray_curtis.png"))
        if (!is.na(o)) written <<- c(written, o)

        g <- plot_beta_dendrogram(sub_mat, lookup, marker, st, style, out_dir = leaf)
        if (!is.na(g)) written <<- c(written, g)

        pm <- run_beta_permanova(sub_mat, lookup, marker, st,
                                 formula_rhs = "waktu", strata_var = NULL)
        perm_rows[[paste0(ucode, "_B_", fld)]] <<-
          dplyr::mutate(pm, field = fld, goal = "B", .before = 1)
      }

      ## Goal D — pooled, faceted by field
      leafD <- file.path(udir, GOAL_DIR["D"])
      dir.create(leafD, recursive = TRUE, showWarnings = FALSE)
      betaD <- compute_beta_gold(mat)
      oD <- plot_ordination(betaD, meta, color_var = "timepoint", shape_var = NULL,
             facet_var = "field", pal = time_pal,
             title = paste0(ucode, " — Goal D (faceted) — Bray-Curtis PCoA"),
             out_path = file.path(leafD, "beta_ord_bray_curtis_faceted.png"))
      if (!is.na(oD)) written <<- c(written, oD)

      gD <- plot_beta_dendrogram(mat, lookup, marker, st, style, out_dir = leafD,
                                 fname_suffix = "pooled")
      if (!is.na(gD)) written <<- c(written, gD)

      pmD <- run_beta_permanova(mat, lookup, marker, st,
                                formula_rhs = "waktu + fertilizer", strata_var = "field")
      perm_rows[[paste0(ucode, "_D")]] <<-
        dplyr::mutate(pmD, field = "pooled", goal = "D", .before = 1)
    }
  }
  run_marker(gold_rarefied_species_16s, gold_universe_lookup_16s, "16S")
  run_marker(gold_rarefied_species_its, gold_universe_lookup_its, "ITS")

  perm_tbl <- dplyr::bind_rows(perm_rows)
  perm_path <- file.path(root, "beta_goal_permanova.csv")
  dir.create(dirname(perm_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(perm_tbl, perm_path)
  message("Wrote ", perm_path)
  print(perm_tbl)

  message("Gold beta goal tree complete: ", length(written), " plots.")
  list(plots = written, permanova_csv = perm_path)
}
