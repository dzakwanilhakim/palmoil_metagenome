# =============================================================================
# R/functions_beta.R — Stage 4: beta diversity (per universe, into the tree)
#   Three metrics: Aitchison (CLR->euclidean, PCA), Bray-Curtis (rarefied,PCoA),
#                  Jaccard (rarefied, binary, PCoA).
#   Outputs per goal: ordination (95% ellipse where n permits), Ward.D2
#                     dendrogram (leaves colored by grouping var), PERMANOVA.
#   Goal placement:
#     A: gated (>=3 groups w/ n>=2, or total n>=10) -> ~fertilizer
#     B: ordination colored by timepoint + PERMANOVA ~timepoint (per field)
#     C: full beta, color=fertilizer shape=field, PERMANOVA ~fertilizer strata=field
#     D: full beta, color=fertilizer shape=field, PERMANOVA ~timepoint+fertilizer strata=field
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(vegan) })

MIN_N_BETA <- 2

# ---- distance computation ---------------------------------------------------
# returns a list with $dist (dist obj), $ord (data.frame of axis1/2 + ids),
# $method, $ordination_type, $var_explained (named vec for PCA/PCoA axes 1-2)
compute_beta <- function(norm_tables, marker, ids, metric) {
  nt <- norm_tables[[marker]]
  if (is.null(nt)) return(NULL)

  if (metric == "aitchison") {
    m <- nt$clr[rownames(nt$clr) %in% ids, , drop = FALSE]
    if (nrow(m) < 3) return(NULL)
    d <- stats::dist(m, method = "euclidean")          # Aitchison
    pco <- stats::prcomp(m, center = TRUE, scale. = FALSE)
    scores <- as.data.frame(pco$x[, 1:2]); names(scores) <- c("Axis1","Axis2")
    ve <- (pco$sdev^2 / sum(pco$sdev^2))[1:2]
    otype <- "PCA"
  } else {
    m <- nt$rarefied[rownames(nt$rarefied) %in% ids, , drop = FALSE]
    if (nrow(m) < 3) return(NULL)
    binary <- identical(metric, "jaccard")
    meth   <- if (binary) "jaccard" else "bray"
    d <- vegan::vegdist(m, method = meth, binary = binary)
    pcoa <- stats::cmdscale(d, k = 2, eig = TRUE)
    scores <- as.data.frame(pcoa$points); names(scores) <- c("Axis1","Axis2")
    eig <- pcoa$eig[pcoa$eig > 0]
    ve  <- (pcoa$eig / sum(eig))[1:2]
    otype <- "PCoA"
  }
  scores$id_sampel <- rownames(m)
  list(dist = d, ord = scores, method = metric,
       ordination_type = otype, var_explained = ve)
}

# group sizes ok for an ellipse (need >=4 points typically)
.can_ellipse <- function(g) {
  tab <- table(g); any(tab >= 4)
}

# ---- ordination plot --------------------------------------------------------
# color_var always; shape_var optional (field for C/D). meta has grouping cols.
plot_ordination <- function(beta, meta, color_var, shape_var = NULL,
                            pal, title, out_path) {
  if (is.null(beta)) return(NA_character_)
  df <- dplyr::left_join(beta$ord, meta, by = "id_sampel")
  ve <- round(100 * beta$var_explained, 1)

  aes_base <- if (is.null(shape_var))
    ggplot2::aes(Axis1, Axis2, colour = .data[[color_var]]) else
    ggplot2::aes(Axis1, Axis2, colour = .data[[color_var]], shape = .data[[shape_var]])

  p <- ggplot2::ggplot(df, aes_base) +
    ggplot2::geom_point(size = 2.6, alpha = .85) +
    ggplot2::scale_colour_manual(values = pal, name = color_var) +
    ggplot2::labs(title = title,
      x = paste0(beta$ordination_type, " Axis 1 (", ve[1], "%)"),
      y = paste0(beta$ordination_type, " Axis 2 (", ve[2], "%)")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))

  # 95% ellipse where a color-group has >=4 points
  if (.can_ellipse(df[[color_var]]))
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = .data[[color_var]]),
                                   level = .95, linewidth = .4, na.rm = TRUE)
  if (!is.null(shape_var))
    p <- p + ggplot2::scale_shape_manual(
      values = rep(c(15,16,17,18,3,7,8,9,10,12,13,14), 3)[seq_len(dplyr::n_distinct(df[[shape_var]]))],
      name = shape_var)

  ggplot2::ggsave(out_path, p, width = 7, height = 5.5, dpi = 200)
  out_path
}

# ---- dendrogram (Ward.D2, leaves colored by grouping var) -------------------
plot_dendrogram <- function(beta, meta, group_var, pal, title, out_path) {
  if (is.null(beta) || attr(beta$dist, "Size") < 3) return(NA_character_)
  hc <- stats::hclust(beta$dist, method = "ward.D2")
  lab_ids <- hc$labels[hc$order]
  grp <- meta[[group_var]][match(lab_ids, meta$id_sampel)]
  cols <- pal[as.character(grp)]; cols[is.na(cols)] <- "grey50"

  grDevices::png(out_path, width = max(1200, 30 * length(lab_ids)),
                 height = 800, res = 150)
  op <- graphics::par(mar = c(8, 4, 3, 1))
  dend <- stats::as.dendrogram(hc)
  plot(dend, main = title, ylab = "Ward.D2 distance", xlab = "")
  # color tick labels
  graphics::mtext(side = 1, at = seq_along(lab_ids), text = lab_ids,
                  col = cols, las = 2, line = .5, cex = .6)
  graphics::legend("topright", legend = names(pal), fill = pal, bty = "n", cex = .8)
  graphics::par(op); grDevices::dev.off()
  out_path
}

# ---- PERMANOVA --------------------------------------------------------------
# formula_rhs e.g. "fertilizer" or "timepoint + fertilizer"; strata optional.
run_permanova <- function(beta, meta, formula_rhs, strata_var = NULL,
                          min_n = MIN_N_BETA) {
  if (is.null(beta)) return(tibble::tibble(term=NA, note="no ordination"))
  md <- meta[match(labels(beta$dist), meta$id_sampel), , drop = FALSE]
  factors <- trimws(strsplit(formula_rhs, "\\+")[[1]])

  # n-gate: each factor must have >=2 levels each with >= min_n
  for (f in factors) {
    tab <- table(md[[f]])
    if (length(tab) < 2 || all(tab < min_n))
      return(tibble::tibble(term = f,
             note = sprintf("skipped: factor '%s' lacks replication (n>=%d, >=2 levels)",
                            f, min_n)))
  }

  fml <- stats::as.formula(paste("beta$dist ~", formula_rhs))
  set.seed(42)
  ad <- tryCatch({
    if (!is.null(strata_var) && strata_var %in% names(md)) {
      perm <- permute::how(nperm = 999)
      permute::setBlocks(perm) <- md[[strata_var]]
      vegan::adonis2(fml, data = md, permutations = perm, by = "terms")
    } else {
      vegan::adonis2(fml, data = md, permutations = 999, by = "terms")
    }
  }, error = function(e) {
    # fallback: legacy strata argument
    tryCatch(vegan::adonis2(fml, data = md, permutations = 999,
              strata = if (!is.null(strata_var)) md[[strata_var]] else NULL,
              by = "terms"),
             error = function(e2) NULL)
  })
  if (is.null(ad)) return(tibble::tibble(term=NA, note="adonis2 failed"))

  tibble::tibble(term = rownames(ad),
                 Df = ad$Df, R2 = ad$R2, F = ad$F, p = ad$`Pr(>F)`,
                 strata = strata_var %||% NA, note = NA_character_)
}

# ---- Goal A gate ------------------------------------------------------------
beta_gate_A <- function(d) {
  tab <- table(d$fertilizer)
  (sum(tab >= 2) >= 3) || (nrow(d) >= 10)
}

# ---- driver: all beta outputs for a subset, into a leaf dir -----------------
# group_color = fertilizer/timepoint; group_shape = field or NULL
beta_for_subset <- function(norm_tables, d, marker, leaf, tag,
                            color_var, shape_var, perm_rhs, strata_var,
                            fname_suffix = "",
                            metrics = c("aitchison","bray_curtis","jaccard")) {
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  ids  <- d$id_sampel
  meta <- d
  pal  <- fert_palette(d$fertilizer)
  if (color_var == "timepoint")
    pal <- setNames(scales::hue_pal()(dplyr::n_distinct(d$timepoint)),
                    sort(unique(d$timepoint)))
  sfx <- if (nzchar(fname_suffix)) paste0("_", fname_suffix) else ""

  written <- character(0); perm_rows <- list()
  for (met in metrics) {
    beta <- compute_beta(norm_tables, marker, ids, met)
    if (is.null(beta)) next

    o <- plot_ordination(beta, meta, color_var, shape_var, pal,
          title = paste0(tag, " — ", met, " ordination"),
          out_path = file.path(leaf, paste0("beta_ord_", met, sfx, ".png")))
    if (!is.na(o)) written <- c(written, o)

    g <- plot_dendrogram(beta, meta, color_var, pal,
          title = paste0(tag, " — ", met, " dendrogram (Ward.D2)"),
          out_path = file.path(leaf, paste0("beta_dendro_", met, sfx, ".png")))
    if (!is.na(g)) written <- c(written, g)

    pm <- run_permanova(beta, meta, perm_rhs, strata_var)
    perm_rows[[met]] <- dplyr::mutate(pm, metric = met, .before = 1)
  }
  if (length(perm_rows))
    readr::write_csv(dplyr::bind_rows(perm_rows),
                     file.path(leaf, paste0("permanova", sfx, ".csv")))
  written
}

# ============================================================================
# MAIN: beta tree across universes (depth-first), into the SAME Results/ tree
# ============================================================================
build_beta_tree <- function(alpha_df, norm_tables, root = "Results") {
  universes <- list_universes(alpha_df)
  written <- character(0)

  for (u in seq_len(nrow(universes))) {
    mk <- universes$marker[u]; st <- universes$stage[u]
    ucode <- paste0(mk, "_", st)
    udir  <- file.path(root, ucode)
    d_u   <- dplyr::filter(alpha_df, marker == mk, stage == st)
    if (nrow(d_u) == 0) next
    message("=== BETA universe ", ucode, " ===")
    fields <- sort(unique(d_u$field))

    ## Goal A — gated, per field x timepoint, ~fertilizer
    for (fld in fields) {
      dsub_f <- dplyr::filter(d_u, field == fld)
      for (tp in sort(unique(dsub_f$timepoint))) {
        dsub <- dplyr::filter(dsub_f, timepoint == tp)
        if (!beta_gate_A(dsub)) next   # strict gate; skip silently
        leaf <- file.path(udir, GOAL_DIR["A"], fld)
        w <- beta_for_subset(norm_tables, dsub, mk, leaf,
              tag = paste0(ucode," A ",fld," ",tp),
              color_var = "fertilizer", shape_var = NULL,
              perm_rhs = "fertilizer", strata_var = NULL,
              fname_suffix = tp)
        written <- c(written, w)
      }
    }

    ## Goal B — per field, ordination colored by timepoint, ~timepoint
    for (fld in fields) {
      dsub <- dplyr::filter(d_u, field == fld)
      if (dplyr::n_distinct(dsub$timepoint) < 2) next
      leaf <- file.path(udir, GOAL_DIR["B"], fld)
      w <- beta_for_subset(norm_tables, dsub, mk, leaf,
            tag = paste0(ucode," B ",fld),
            color_var = "timepoint", shape_var = NULL,
            perm_rhs = "timepoint", strata_var = NULL)
      written <- c(written, w)
    }

    ## Goal C — pooled per timepoint, color=fertilizer shape=field, strata=field
    leafC <- file.path(udir, GOAL_DIR["C"])
    for (tp in sort(unique(d_u$timepoint))) {
      dsub <- dplyr::filter(d_u, timepoint == tp)
      sub_leaf <- file.path(leafC, paste0("timepoint_", tp))
      w <- beta_for_subset(norm_tables, dsub, mk, sub_leaf,
            tag = paste0(ucode," C ",tp),
            color_var = "fertilizer", shape_var = "field",
            perm_rhs = "fertilizer", strata_var = "field")
      written <- c(written, w)
    }

    ## Goal D — all data, color=fertilizer shape=field, ~timepoint+fertilizer strata=field
    leafD <- file.path(udir, GOAL_DIR["D"])
    w <- beta_for_subset(norm_tables, d_u, mk, leafD,
          tag = paste0(ucode," D"),
          color_var = "fertilizer", shape_var = "field",
          perm_rhs = "timepoint + fertilizer", strata_var = "field")
    written <- c(written, w)
  }
  message("Beta tree complete: ", length(written), " plots.")
  written
}
