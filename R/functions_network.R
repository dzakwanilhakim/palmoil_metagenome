# =============================================================================
# R/functions_network.R
#   Co-occurrence network analysis + Active Keystone identification.
#   Strictly per universe (16S_TM, 16S_Nursery, ITS_TM, ITS_Nursery).
#
#   Keystone logic (3-step intersection):
#     1. Global network -> rank all genera by Node Degree
#     2. Welch t-test CLR -> genera significant (p<0.05) in >=1 comparison
#     3. Intersect -> re-rank by degree -> Top 10 = Active Keystones
#
#   Outputs per universe:
#     cooccurrence/global/
#       network_hubs.jpg, heatmap.png, keystone_active10.csv,
#       volcano_genus/species, dumbbell_genus.png, dumbbell_species.png
#     cooccurrence/per_field/<field>/
#       network_hubs.jpg, heatmap.png, keystone_top10.csv
#     cooccurrence/per_fertilizer/<fert>/
#       (same)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(igraph); library(ggraph); library(pheatmap)
})

MIN_N_NET   <- 5      # minimum n to attempt a network (n-gate)
HUB_TOP_PCT <- 0.20   # top 20% by degree = hub
HUB_MIN_DEG <- 3      # absolute degree floor for hub plots
R_THRESH    <- 0.6    # |r| threshold for edges
TOP_K       <- 10     # active keystones to report
DA_P_THRESH <- 0.05   # significance threshold for DA filter

# ============================================================================
# A.  CLR helper
# ============================================================================

# CLR from a sample x taxon count matrix. Pseudocount +1. Drops zero-var cols.
clr_matrix <- function(mat, pseudo = 1L) {
  m <- mat + pseudo
  m <- m[, colSums(m) > 0, drop = FALSE]
  log_m <- log(m)
  sweep(log_m, 1, rowMeans(log_m), "-")
}

# Extract genus CLR for one universe from norm_tables (uses the stored CLR).
genus_clr_for_universe <- function(norm_tables, master_samples, marker, stage) {
  clr <- norm_tables[[marker]]$clr
  if (is.null(clr)) return(NULL)
  ids <- master_samples |>
    dplyr::filter(marker == !!marker, stage == !!stage) |>
    dplyr::pull(id_sampel)
  ids <- intersect(ids, rownames(clr))
  if (length(ids) < 2) return(NULL)
  clr[ids, , drop = FALSE]
}

# ============================================================================
# B.  Network construction
# ============================================================================

# Spearman |r| >= R_THRESH network from a CLR matrix.
# Returns igraph object or NULL if n < MIN_N_NET or no edges pass threshold.
build_network <- function(clr_mat, group_label,
                          r_thresh = R_THRESH, min_n = MIN_N_NET) {
  n <- if (is.null(clr_mat)) 0L else nrow(clr_mat)
  if (n < min_n) {
    message("  SKIP (n=", n, " < ", min_n, "): ", group_label)
    return(NULL)
  }
  keep <- apply(clr_mat, 2, var, na.rm = TRUE) > 0
  m <- clr_mat[, keep, drop = FALSE]
  if (ncol(m) < 2) { message("  SKIP (<2 variable taxa): ", group_label); return(NULL) }

  cr <- suppressWarnings(cor(m, method = "spearman"))
  ut <- upper.tri(cr)
  idx <- which(ut & abs(cr) >= r_thresh, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    message("  SKIP (0 edges at |r|>=", r_thresh, "): ", group_label)
    return(NULL)
  }
  edges <- tibble::tibble(
    from  = rownames(cr)[idx[,1]],
    to    = colnames(cr)[idx[,2]],
    r     = cr[idx],
    color = ifelse(cr[idx] >= 0, "#27AE60", "#E74C3C"))

  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
         vertices = data.frame(name = colnames(m)))
  # drop isolated single nodes (no edges) for efficiency
  iso <- igraph::V(g)[igraph::degree(g) == 0]
  if (length(iso) > 0) g <- igraph::delete_vertices(g, iso)
  igraph::V(g)$degree      <- igraph::degree(g)
  igraph::V(g)$betweenness <- igraph::betweenness(g, normalized = TRUE)
  igraph::E(g)$color       <- edges$color[match(
    paste(igraph::as_data_frame(g,"edges")$from, igraph::as_data_frame(g,"edges")$to),
    paste(edges$from, edges$to))]
  igraph::E(g)$r           <- edges$r[match(
    paste(igraph::as_data_frame(g,"edges")$from, igraph::as_data_frame(g,"edges")$to),
    paste(edges$from, edges$to))]
  message("  OK n=", n, " | nodes=", igraph::vcount(g),
          " edges=", igraph::ecount(g), ": ", group_label)
  g
}

# ============================================================================
# C.  Network visualization: full + hub-only
# ============================================================================

.draw_network <- function(g, title, out_path,
                          hub_only = FALSE,
                          hub_pct  = HUB_TOP_PCT,
                          hub_min  = HUB_MIN_DEG) {
  if (is.null(g) || igraph::ecount(g) == 0) return(NA_character_)

  if (hub_only) {
    deg   <- igraph::V(g)$degree
    thr   <- max(hub_min, stats::quantile(deg, 1 - hub_pct, names = FALSE))
    keep  <- igraph::V(g)$name[deg >= thr]
    if (length(keep) < 2) {
      message("  hub plot skipped (<2 hub nodes)"); return(NA_character_)
    }
    g <- igraph::induced_subgraph(g, keep)
  }

  n_nodes <- igraph::vcount(g)
  set.seed(42)
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(colour = color), alpha = .55, linewidth = .5,
      show.legend = FALSE) +
    ggraph::scale_edge_colour_identity() +
    ggraph::geom_node_point(
      ggplot2::aes(size = degree), colour = "#2C3E50", alpha = .82) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name, size = pmax(degree * .30, 2)),
      repel = TRUE, max.overlaps = 25,
      colour = "black", show.legend = FALSE) +
    ggplot2::scale_size_continuous(range = c(2, 9), name = "Degree") +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(n_nodes, " nodes | ",
                        igraph::ecount(g), " edges | green=positive, red=negative")) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      legend.position = "right")

  # cap canvas: ragg max is 50000px; at 180dpi keep well under (<= ~40in)
  w <- min(40, max(9, 1.2 + 0.28 * n_nodes))
  h <- min(40, max(7, 0.9 + 0.25 * n_nodes))
  ok <- tryCatch({ ggplot2::ggsave(out_path, p, width = w, height = h,
                  dpi = 150, limitsize = FALSE); TRUE },
                  error = function(e) { message("  network plot error: ",
                  conditionMessage(e)); FALSE })
  if (ok && file.exists(out_path)) out_path else NA_character_
}

plot_network_full <- function(g, title, out_path)
  .draw_network(g, title, out_path, hub_only = FALSE)

plot_network_hubs <- function(g, title, out_path)
  .draw_network(g, title, out_path, hub_only = TRUE)

# ============================================================================
# D.  Correlation heatmap
# ============================================================================

plot_corr_heatmap <- function(g, title, out_path, max_nodes = 60) {
  if (is.null(g) || igraph::vcount(g) < 2) return(NA_character_)

  # restrict to the most-connected nodes: a 1270x1270 heatmap is unreadable
  # and slow. Keep the top `max_nodes` by degree (hub-focused heatmap).
  deg <- igraph::degree(g)
  if (length(deg) > max_nodes) {
    keep_nodes <- names(sort(deg, decreasing = TRUE))[seq_len(max_nodes)]
    g <- igraph::induced_subgraph(g, keep_nodes)
  }
  nodes <- igraph::V(g)$name
  if (length(nodes) < 2) return(NA_character_)

  # rebuild Spearman r matrix from edge list
  adj <- matrix(0, length(nodes), length(nodes),
                dimnames = list(nodes, nodes))
  diag(adj) <- 1
  edf <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edf) > 0) for (i in seq_len(nrow(edf))) {
    adj[edf$from[i], edf$to[i]] <- edf$r[i]
    adj[edf$to[i], edf$from[i]] <- edf$r[i]
  }

  cell_px <- max(10, min(40, 1000 %/% length(nodes)))
  px_w <- min(8000, max(800, cell_px * length(nodes) + 250))
  px_h <- min(8000, max(700, cell_px * length(nodes) + 200))

  ok <- tryCatch({
    grDevices::png(out_path, width = px_w, height = px_h, res = 130)
    pheatmap::pheatmap(
      adj, main = title,
      fontsize = max(5, min(9, 600 %/% length(nodes))),
      color = grDevices::colorRampPalette(c("#E74C3C","white","#27AE60"))(100),
      breaks = seq(-1, 1, length.out = 101),
      treeheight_row = 18, treeheight_col = 18,
      border_color = NA, silent = TRUE)
    grDevices::dev.off()
    TRUE
  }, error = function(e) {
    message("  heatmap error: ", conditionMessage(e))
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    FALSE
  })

  # only return the path if the file was actually created
  if (ok && file.exists(out_path)) out_path else NA_character_
}

# ============================================================================
# E.  Keystone table from network (topology only — used for per-field/fert)
# ============================================================================

keystone_table <- function(g, top_k = TOP_K) {
  if (is.null(g)) return(tibble::tibble())
  tibble::tibble(
    genus       = igraph::V(g)$name,
    degree      = igraph::V(g)$degree,
    betweenness = igraph::V(g)$betweenness) |>
    dplyr::arrange(dplyr::desc(degree), dplyr::desc(betweenness)) |>
    dplyr::slice_head(n = top_k)
}

# ============================================================================
# F.  Differential Abundance (Welch t-test, consecutive timepoints)
#     Returns tidy data frame of all genera that are significant in >=1 pair.
# ============================================================================

da_significant_genera <- function(clr_mat, meta, p_thresh = DA_P_THRESH) {
  df <- as.data.frame(clr_mat)
  df$id_sampel <- rownames(clr_mat)
  df <- dplyr::inner_join(df,
    dplyr::select(meta, id_sampel, timepoint), by = "id_sampel")
  tps  <- sort(unique(df$timepoint))
  taxa <- setdiff(names(df), c("id_sampel", "timepoint"))
  if (length(tps) < 2 || length(taxa) == 0) return(character(0))

  sig_genera <- character(0)
  for (i in seq_len(length(tps) - 1)) {
    t0 <- tps[i]; t1 <- tps[i + 1]
    d0 <- dplyr::filter(df, timepoint == t0)
    d1 <- dplyr::filter(df, timepoint == t1)
    if (nrow(d0) < 2 || nrow(d1) < 2) next
    for (tx in taxa) {
      tt <- tryCatch(stats::t.test(d1[[tx]], d0[[tx]], var.equal = FALSE),
                     error = function(e) NULL)
      if (!is.null(tt) && tt$p.value < p_thresh)
        sig_genera <- c(sig_genera, tx)
    }
  }
  unique(sig_genera)
}

# ============================================================================
# G.  Active Keystone identification (3-step intersection)
# ============================================================================

active_keystones <- function(g, clr_mat, meta,
                             top_k = TOP_K, p_thresh = DA_P_THRESH) {
  if (is.null(g)) return(tibble::tibble())

  # Step 1: rank all genera in the network by degree
  topo <- tibble::tibble(
    genus       = igraph::V(g)$name,
    degree      = igraph::V(g)$degree,
    betweenness = igraph::V(g)$betweenness) |>
    dplyr::arrange(dplyr::desc(degree))

  # Step 2: DA filter — genera significant in >=1 consecutive comparison
  sig <- da_significant_genera(clr_mat, meta, p_thresh)
  message("  DA significant genera: ", length(sig))

  # Step 3: intersect + re-rank by degree -> top_k
  ks <- topo |>
    dplyr::filter(genus %in% sig) |>
    dplyr::arrange(dplyr::desc(degree)) |>
    dplyr::slice_head(n = top_k) |>
    dplyr::mutate(keystone_rank = dplyr::row_number())

  message("  Active keystones after intersection: ", nrow(ks))
  ks
}

# ============================================================================
# H.  Per-subset network helper (shared by field, fertilizer, global)
# ============================================================================

.run_network_subset <- function(clr_sub, label, leaf, tag,
                                compute_keystones = FALSE,
                                clr_full = NULL, meta_full = NULL,
                                r_thresh = R_THRESH,
                                species_long = NULL, marker = NULL, stage = NULL) {
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  g <- build_network(clr_sub, label, r_thresh = r_thresh)

  if (is.null(g)) {
    readr::write_csv(
      tibble::tibble(note = paste0("skipped n<", MIN_N_NET, ": ", label)),
      file.path(leaf, "network_skipped.csv"))
    return(list(g = NULL, written = character(0),
                keystones = tibble::tibble()))
  }

  written <- character(0)

  # only hub network (full network omitted per spec)
  w <- plot_network_hubs(g, paste0(tag, " | hub network"),
                          file.path(leaf, "network_hubs.jpg"))
  if (!is.na(w)) written <- c(written, w)

  w <- plot_corr_heatmap(g, paste0(tag, " | correlation heatmap"),
                          file.path(leaf, "heatmap.png"))
  if (!is.na(w)) written <- c(written, w)

  # keystones: Active (topology x DA intersection)
  if (compute_keystones && !is.null(clr_full) && !is.null(meta_full)) {
    ks <- active_keystones(g, clr_full, meta_full)
    readr::write_csv(ks, file.path(leaf, "keystone_active10.csv"))

    # ---- volcano plots: genus ------------------------------------------
    written <- c(written,
      volcano_for_subset(clr_full, meta_full, "genus", leaf, tag))

    # ---- volcano plots: species (built from species_long for these ids) -
    if (!is.null(species_long) && !is.null(marker) && !is.null(stage)) {
      sp_clr <- species_clr_for_ids(species_long, meta_full$id_sampel, marker, stage)
      if (!is.null(sp_clr))
        written <- c(written,
          volcano_for_subset(sp_clr, meta_full, "species", leaf, tag))
    } else sp_clr <- NULL

    # ---- dumbbell: genus (Top-10 keystone genera) ----------------------
    if (nrow(ks) > 0) {
      w <- plot_genus_dumbbell(clr_full, meta_full, ks,
             paste0(tag, " | keystone genera dumbbell"),
             file.path(leaf, "dumbbell_genus.png"))
      if (!is.na(w)) written <- c(written, w)

      # ---- dumbbell: species (top-1 species per keystone genus) --------
      if (!is.null(sp_clr)) {
        w <- plot_species_dumbbell(sp_clr, meta_full, ks,
               paste0(tag, " | keystone species dumbbell"),
               file.path(leaf, "dumbbell_species.png"))
        if (!is.na(w)) written <- c(written, w)
      }
    }
  } else {
    ks <- keystone_table(g)
    readr::write_csv(ks, file.path(leaf, "keystone_top10.csv"))
  }
  list(g = g, written = written, keystones = ks)
}

# ============================================================================
# I.  MAIN DRIVER  build_network_tree()
# ============================================================================

build_network_tree <- function(norm_tables, master_samples, species_long = NULL,
                                root = "Results") {
  universes <- master_samples |>
    dplyr::distinct(marker, stage) |>
    dplyr::filter(!is.na(stage)) |>
    dplyr::arrange(marker, stage)
  all_written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)
    message("\n====== NETWORK universe ", ucode, " ======")

    meta <- dplyr::filter(master_samples, marker == mk, stage == st)
    clr_g <- genus_clr_for_universe(norm_tables, master_samples, mk, st)
    if (is.null(clr_g)) {
      message("  no CLR data — skipping"); next
    }
    meta <- meta[meta$id_sampel %in% rownames(clr_g), , drop = FALSE]

    # ---- Per-field networks (Active Keystone per field) ---------------------
    # Pool ALL fertilizers + ALL timepoints within each field (ignore fertilizer)
    fields <- sort(unique(meta$field))
    message(" -- Per-field networks (", length(fields), " fields, all ferts pooled) --")
    for (fld in fields) {
      ids_f  <- meta$id_sampel[meta$field == fld]   # all ferts, all timepoints
      leaf_f <- file.path(udir, "cooccurrence", "per_field", fld)
      meta_f <- meta[meta$id_sampel %in% ids_f, , drop = FALSE]
      res_f  <- .run_network_subset(
        clr_g[ids_f, , drop = FALSE],
        paste0(ucode, "|field:", fld), leaf_f,
        paste0(ucode, " | field: ", fld, " (all ferts pooled)"),
        compute_keystones = TRUE,
        clr_full = clr_g[ids_f, , drop = FALSE], meta_full = meta_f,
        species_long = species_long, marker = mk, stage = st)
      all_written <- c(all_written, res_f$written)
    }

    # ---- 3. Per-fertilizer networks -----------------------------------------
    # Pool ALL fields + ALL timepoints for each fertilizer (ignore field)
    ferts <- sort(unique(meta$fertilizer))
    message(" -- Per-fertilizer networks (", length(ferts), " ferts, all fields pooled) --")
    for (ft in ferts) {
      ids_ft  <- meta$id_sampel[meta$fertilizer == ft]   # all fields, all timepoints
      leaf_ft <- file.path(udir, "cooccurrence", "per_fertilizer", ft)
      meta_ft <- meta[meta$id_sampel %in% ids_ft, , drop = FALSE]
      res_ft  <- .run_network_subset(
        clr_g[ids_ft, , drop = FALSE],
        paste0(ucode, "|fert:", ft), leaf_ft,
        paste0(ucode, " | fertilizer: ", ft, " (all fields pooled)"),
        compute_keystones = TRUE,
        clr_full = clr_g[ids_ft, , drop = FALSE], meta_full = meta_ft,
        species_long = species_long, marker = mk, stage = st)
      all_written <- c(all_written, res_ft$written)
    }

    message(" Universe ", ucode, " done: ",
            length(all_written), " total outputs so far.")
  }

  message("\nNetwork tree complete: ", length(all_written), " outputs.")
  all_written
}

# ============================================================================
# J.  Species CLR helper (for volcano + zoom-in within a subset)
# ============================================================================

# Build species CLR matrix for a given set of sample ids, with genus map attr.
species_clr_for_ids <- function(species_long, ids, marker, stage) {
  d <- species_long |>
    dplyr::filter(marker == !!marker, stage == !!stage, id_sampel %in% ids) |>
    dplyr::group_by(id_sampel, species, genus) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  if (nrow(d) == 0) return(NULL)
  wide <- d |>
    dplyr::select(id_sampel, species, count) |>
    tidyr::pivot_wider(names_from = species, values_from = count,
                       values_fill = 0)
  mat <- as.matrix(dplyr::select(wide, -id_sampel))
  rownames(mat) <- wide$id_sampel
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)
  clr_m <- clr_matrix(mat)
  attr(clr_m, "genus_map") <- dplyr::distinct(d, species, genus)
  clr_m
}

# ============================================================================
# K.  Volcano plots (Welch t-test on CLR, consecutive timepoints)
#     rank = "genus" or "species". Returns written file paths.
# ============================================================================

.volcano_plot <- function(res, title, out_path) {
  if (nrow(res) == 0 || all(is.na(res$p))) return(NA_character_)
  res$col <- dplyr::case_when(
    res$sig & res$clr_diff >  0 ~ "up",
    res$sig & res$clr_diff <  0 ~ "down",
    TRUE                        ~ "ns")
  pal <- c(up = "#E74C3C", down = "#2980B9", ns = "grey70")

  p <- ggplot2::ggplot(res, ggplot2::aes(clr_diff, -log10(p), colour = col)) +
    ggplot2::geom_point(size = 1.8, alpha = .75) +
    ggplot2::scale_colour_manual(values = pal, name = NULL,
      labels = c(up = "Higher at T+1", down = "Lower at T+1", ns = "n.s."),
      breaks = c("up","down","ns")) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = .4) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                        colour = "grey50", linewidth = .4) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label),
      size = 2.3, max.overlaps = 12, na.rm = TRUE, colour = "black") +
    ggplot2::labs(title = title,
                  x = "CLR difference (mean T+1 \u2212 mean T0)",
                  y = expression(-log[10](p))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10))

  ok <- tryCatch({ ggplot2::ggsave(out_path, p, width = 7, height = 5.5, dpi = 170); TRUE },
                  error = function(e) { message("  volcano error: ",
                  conditionMessage(e)); FALSE })
  if (ok && file.exists(out_path)) out_path else NA_character_
}

volcano_for_subset <- function(clr_mat, meta, rank_label, leaf,
                               universe_label) {
  if (is.null(clr_mat) || nrow(clr_mat) < 4) return(character(0))
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(clr_mat)
  df$id_sampel <- rownames(clr_mat)
  df <- dplyr::inner_join(df,
    dplyr::select(meta, id_sampel, timepoint), by = "id_sampel")
  tps  <- sort(unique(df$timepoint))
  taxa <- setdiff(names(df), c("id_sampel","timepoint"))
  if (length(tps) < 2) return(character(0))
  written <- character(0)

  for (i in seq_len(length(tps) - 1)) {
    t0 <- tps[i]; t1 <- tps[i + 1]
    d0 <- dplyr::filter(df, timepoint == t0)
    d1 <- dplyr::filter(df, timepoint == t1)
    if (nrow(d0) < 2 || nrow(d1) < 2) next

    res <- purrr::map_dfr(taxa, function(tx) {
      x0 <- d0[[tx]]; x1 <- d1[[tx]]
      tt <- tryCatch(stats::t.test(x1, x0, var.equal = FALSE),
                     error = function(e) NULL)
      if (is.null(tt)) return(tibble::tibble(taxon = tx, clr_diff = NA, p = NA))
      tibble::tibble(taxon = tx, clr_diff = mean(x1) - mean(x0), p = tt$p.value)
    }) |>
      dplyr::filter(!is.na(p)) |>
      dplyr::mutate(p_adj = p.adjust(p, "BH"),
                    sig   = p < 0.05,                 # RAW p (per spec)
                    label = ifelse(sig, taxon, NA_character_))

    nm <- paste0(t0, "_vs_", t1)
    readr::write_csv(res, file.path(leaf,
      paste0("diffabund_", rank_label, "_", nm, ".csv")))
    w <- .volcano_plot(res,
      title = paste0(universe_label, " | ", rank_label, " | ", nm),
      out_path = file.path(leaf, paste0("volcano_", rank_label, "_", nm, ".png")))
    if (!is.na(w)) written <- c(written, w)
  }
  written
}

# ============================================================================
# L.  Dumbbell plots (genus + species): rows = taxa, x = mean CLR,
#     dots coloured pale->dark blue across timepoints, line connects them.
# ============================================================================

# Generic dumbbell. long_df: taxon, timepoint, mean_clr.
.dumbbell_plot <- function(long_df, title, out_path, rank_label = "Genus") {
  if (nrow(long_df) == 0) return(NA_character_)
  tps <- sort(unique(long_df$timepoint))
  long_df$timepoint <- factor(long_df$timepoint, levels = tps)

  # order taxa rows by mean CLR at the LAST timepoint (lowest at bottom)
  last_tp <- tps[length(tps)]
  ord <- long_df |>
    dplyr::group_by(taxon) |>
    dplyr::summarise(
      key = {
        v <- mean_clr[timepoint == last_tp]
        if (length(v) == 0) mean(mean_clr) else v[1]
      }, .groups = "drop") |>
    dplyr::arrange(key) |> dplyr::pull(taxon)
  long_df$taxon <- factor(long_df$taxon, levels = ord)

  p <- ggplot2::ggplot(long_df, ggplot2::aes(mean_clr, taxon)) +
    ggplot2::geom_line(ggplot2::aes(group = taxon),
                       colour = "grey55", linewidth = .8) +
    ggplot2::geom_point(ggplot2::aes(fill = timepoint),
                        shape = 21, size = 4, colour = "grey30", stroke = .3) +
    ggplot2::scale_fill_manual(
      values = grDevices::colorRampPalette(c("#cfe3f5", "#08306b"))(length(tps)),
      name = "Timepoint") +
    ggplot2::labs(title = title, x = "Mean CLR abundance", y = rank_label) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = if (rank_label == "Species") "italic" else "plain", size = 9))

  ok <- tryCatch({
    ggplot2::ggsave(out_path, p, width = 8.5,
                    height = max(4, 0.45 * dplyr::n_distinct(long_df$taxon) + 1.5),
                    dpi = 170, limitsize = FALSE); TRUE },
    error = function(e){ message("  dumbbell error: ", conditionMessage(e)); FALSE })
  if (ok && file.exists(out_path)) out_path else NA_character_
}

# Genus dumbbell: Top-10 keystone genera.
plot_genus_dumbbell <- function(clr_mat, meta, keystones, title, out_path) {
  if (nrow(keystones) == 0 || is.null(clr_mat)) return(NA_character_)
  genera <- intersect(keystones$genus, colnames(clr_mat))
  if (length(genera) == 0) return(NA_character_)
  df <- as.data.frame(clr_mat[, genera, drop = FALSE])
  df$id_sampel <- rownames(clr_mat)
  long <- dplyr::inner_join(df,
    dplyr::select(meta, id_sampel, timepoint), by = "id_sampel") |>
    tidyr::pivot_longer(-c(id_sampel, timepoint),
                        names_to = "taxon", values_to = "clr") |>
    dplyr::group_by(taxon, timepoint) |>
    dplyr::summarise(mean_clr = mean(clr), .groups = "drop")
  .dumbbell_plot(long, title, out_path, rank_label = "Genus")
}

# Species dumbbell: top-1 species per keystone genus -> up to 10 species.
plot_species_dumbbell <- function(sp_clr, meta, keystones, title, out_path) {
  if (is.null(sp_clr) || nrow(keystones) == 0) return(NA_character_)
  gmap <- attr(sp_clr, "genus_map")
  if (is.null(gmap)) return(NA_character_)
  genera <- keystones$genus

  df <- as.data.frame(sp_clr)
  df$id_sampel <- rownames(sp_clr)
  long <- dplyr::inner_join(df,
    dplyr::select(meta, id_sampel, timepoint), by = "id_sampel") |>
    tidyr::pivot_longer(-c(id_sampel, timepoint),
                        names_to = "species", values_to = "clr") |>
    dplyr::left_join(gmap, by = "species") |>
    dplyr::filter(genus %in% genera, !is.na(genus))
  if (nrow(long) == 0) return(NA_character_)

  top1 <- long |>
    dplyr::group_by(genus, species) |>
    dplyr::summarise(m = mean(clr), .groups = "drop") |>
    dplyr::group_by(genus) |>
    dplyr::slice_max(m, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |> dplyr::pull(species)

  long <- long |>
    dplyr::filter(species %in% top1) |>
    dplyr::group_by(species, timepoint) |>
    dplyr::summarise(mean_clr = mean(clr), .groups = "drop") |>
    dplyr::rename(taxon = species)
  .dumbbell_plot(long, title, out_path, rank_label = "Species")
}
