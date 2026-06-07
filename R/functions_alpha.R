# =============================================================================
# R/functions_alpha.R — Stage 4: alpha diversity + A/B/C/D engine + stats
#   - alpha metrics (Observed, Pielou, Shannon) from the RAREFIED table
#   - generic goal expander: each goal -> many (subset, group_by) analyses
#   - unpaired stats: KW + Wilcoxon post-hoc (kw_wilcox);
#                     unpaired Wilcoxon rank-sum T0 vs T1 (wilcox_time)
#   - n>=3 gate per compared group; BH adjustment; skipped tests logged
# Runs per marker, independently (16S and ITS never co-analysed).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

load_comparisons <- function(path = "config/comparisons.yaml") {
  yaml::read_yaml(path)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# ALPHA METRICS  (per sample, from rarefied matrix)  per marker
#   Observed = richness; Shannon = vegan::diversity; Pielou = Shannon/log(Obs)
# ----------------------------------------------------------------------------
compute_alpha <- function(norm_tables, master_samples) {
  out <- purrr::map_dfr(names(norm_tables), function(mk) {
    rar <- norm_tables[[mk]]$rarefied
    if (is.null(rar) || nrow(rar) == 0) return(tibble::tibble())

    obs <- vegan::specnumber(rar)
    sha <- vegan::diversity(rar, index = "shannon")
    pie <- ifelse(obs > 1, sha / log(obs), NA_real_)

    tibble::tibble(marker = mk, id_sampel = rownames(rar),
                   Observed = obs, Shannon = sha, Pielou = pie)
  }) |>
    dplyr::left_join(master_samples, by = c("marker", "id_sampel"))

  # drop samples with missing/blank fertilizer (cannot be grouped)
  bad <- is.na(out$fertilizer) | stringr::str_trim(as.character(out$fertilizer)) %in% c("", "NA")
  if (any(bad)) {
    message("Dropped ", sum(bad), " sample(s) with NA/blank fertilizer:")
    print(out[bad, c("marker","id_sampel","stage","field","timepoint")])
  }
  out[!bad, , drop = FALSE]
}

# ----------------------------------------------------------------------------
# GOAL EXPANDER
#   Turn one goal spec into a list of concrete analysis units, each with:
#     $goal, $marker, $context (named list of the `per` values),
#     $data (subset of alpha rows), $group_by
#   `per` values enumerated from the data. "global" = no subsetting.
# ----------------------------------------------------------------------------
expand_goal <- function(alpha_df, goal_name, spec, marker) {
  d <- dplyr::filter(alpha_df, marker == !!marker)
  if (nrow(d) == 0) return(list())

  per <- spec$per
  per_vec <- unlist(per)
  if (length(per_vec) == 1 && identical(as.character(per_vec), "global")) {
    combos <- tibble::tibble(.dummy = 1)
  } else {
    combos <- dplyr::distinct(d, dplyr::across(dplyr::all_of(per_vec)))
  }

  purrr::map(seq_len(nrow(combos)), function(i) {
    ctx <- combos[i, , drop = FALSE]
    sub <- d
    if (!".dummy" %in% names(ctx)) {
      for (col in names(ctx)) sub <- dplyr::filter(sub, .data[[col]] == ctx[[col]][1])
    }
    list(goal = goal_name, marker = marker,
         context = if (".dummy" %in% names(ctx)) list() else as.list(ctx),
         group_by = spec$group_by, spec = spec, data = sub)
  })
}

# ----------------------------------------------------------------------------
# STATS — unpaired, per analysis unit, per metric
# ----------------------------------------------------------------------------
.gate_ok <- function(groups, min_n) {
  tab <- table(groups)
  length(tab) >= 2 && all(tab >= min_n)
}

# KW (>2 groups) + Wilcoxon rank-sum post-hoc; or just Wilcoxon if 2 groups
run_kw_wilcox <- function(values, groups, metric, min_n, padj) {
  groups <- as.factor(groups)
  keep <- !is.na(values)
  values <- values[keep]; groups <- droplevels(groups[keep])
  res <- list(metric = metric, test = NA, p_global = NA_real_,
              posthoc = NULL, skipped = NA_character_)

  if (!.gate_ok(groups, min_n)) {
    res$skipped <- sprintf("n<%d in a group or <2 groups", min_n); return(res)
  }
  k <- nlevels(groups)
  if (k == 2) {
    w <- suppressWarnings(stats::wilcox.test(values ~ groups))
    res$test <- "wilcoxon_ranksum"; res$p_global <- w$p.value
  } else {
    kw <- stats::kruskal.test(values ~ groups)
    res$test <- "kruskal_wallis"; res$p_global <- kw$p.value
    # pairwise Wilcoxon post-hoc, BH-adjusted
    pw <- suppressWarnings(
      stats::pairwise.wilcox.test(values, groups, p.adjust.method = padj))
    res$posthoc <- broom_pw(pw)
  }
  res
}

# unpaired Wilcoxon rank-sum specifically for T0 vs T1
run_wilcox_time <- function(values, timepoints, metric, min_n) {
  tp <- as.factor(timepoints)
  keep <- !is.na(values)
  values <- values[keep]; tp <- droplevels(tp[keep])
  res <- list(metric = metric, test = NA, p_global = NA_real_,
              skipped = NA_character_)
  if (!.gate_ok(tp, min_n) || nlevels(tp) < 2) {
    res$skipped <- sprintf("need 2 timepoints with n>=%d", min_n); return(res)
  }
  w <- suppressWarnings(stats::wilcox.test(values ~ tp))
  res$test <- "wilcoxon_ranksum_T0vsT1"; res$p_global <- w$p.value
  res
}

# tidy a pairwise.wilcox.test matrix into long form
broom_pw <- function(pw) {
  m <- pw$p.value
  if (is.null(m)) return(NULL)
  out <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(out) <- c("group1", "group2", "p_adj")
  out[!is.na(out$p_adj), ]
}

# ----------------------------------------------------------------------------
# DRIVER: run all goals x metrics for one marker, return tidy stats table
# ----------------------------------------------------------------------------
run_alpha_stats <- function(alpha_df, comparisons, marker) {
  ss  <- comparisons$stats_settings
  mets <- comparisons$alpha_metrics
  goals <- comparisons$goals

  rows <- list(); ph <- list()

  for (gname in names(goals)) {
    spec  <- goals[[gname]]
    units <- expand_goal(alpha_df, gname, spec, marker)

    for (u in units) {
      ctx_lab <- if (length(u$context))
        paste(names(u$context), unlist(u$context), sep = "=", collapse = "; ") else "global"

      for (met in mets) {
        vals <- u$data[[met]]

        if (identical(spec$stats, "kw_wilcox")) {
          r <- run_kw_wilcox(vals, u$data[[u$group_by]], met,
                             ss$min_n_per_group, ss$p_adjust_method)
        } else if (identical(spec$stats, "wilcox_time")) {
          r <- run_wilcox_time(vals, u$data$timepoint, met, ss$min_n_per_group)
        } else next

        rows[[length(rows)+1]] <- tibble::tibble(
          marker = marker, goal = gname, context = ctx_lab,
          group_by = u$group_by, metric = met,
          test = r$test %||% NA, p_global = r$p_global %||% NA_real_,
          skipped = r$skipped %||% NA_character_)

        if (!is.null(r$posthoc) && nrow(r$posthoc))
          ph[[length(ph)+1]] <- tibble::tibble(
            marker = marker, goal = gname, context = ctx_lab, metric = met,
            r$posthoc)
      }

      # stats_extra: per_fertilizer time effect (B, D) / per_field (C)
      if (!is.null(spec$stats_extra)) {
        extra <- run_stats_extra(u, spec, mets, ss, marker, gname, ctx_lab)
        rows <- c(rows, extra$rows); ph <- c(ph, extra$ph)
      }
    }
  }

  global_tbl <- if (length(rows)) dplyr::bind_rows(rows) else tibble::tibble()
  if (nrow(global_tbl)) global_tbl$marker <- as.character(global_tbl$marker)

  # BH-adjust p_global WITHIN each goal family (only over tests that ran)
  if (nrow(global_tbl)) {
    global_tbl <- global_tbl |>
      dplyr::group_by(goal) |>
      dplyr::mutate(p_global_BH = ifelse(is.na(p_global), NA_real_,
                      stats::p.adjust(p_global, method = ss$p_adjust_method))) |>
      dplyr::ungroup() |>
      dplyr::relocate(p_global_BH, .after = p_global)
  }

  list(
    global = global_tbl,
    posthoc = if (length(ph)) dplyr::bind_rows(ph) else tibble::tibble()
  )
}

# helper for the per_fertilizer / per_field secondary tests
run_stats_extra <- function(u, spec, mets, ss, marker, gname, ctx_lab) {
  rows <- list(); ph <- list()
  split_col <- dplyr::case_when(
    spec$stats_extra == "per_fertilizer" ~ "fertilizer",
    spec$stats_extra == "per_field"      ~ "field",
    TRUE ~ NA_character_)
  if (is.na(split_col) || !split_col %in% names(u$data)) return(list(rows=rows, ph=ph))

  for (lv in unique(u$data[[split_col]])) {
    sub <- dplyr::filter(u$data, .data[[split_col]] == lv)
    for (met in mets) {
      vals <- sub[[met]]
      if (identical(spec$stats, "wilcox_time")) {
        r <- run_wilcox_time(vals, sub$timepoint, met, ss$min_n_per_group)
      } else {
        r <- run_kw_wilcox(vals, sub[[u$group_by]], met,
                           ss$min_n_per_group, ss$p_adjust_method)
      }
      rows[[length(rows)+1]] <- tibble::tibble(
        marker = marker, goal = paste0(gname, "_", spec$stats_extra),
        context = paste0(ctx_lab, " | ", split_col, "=", lv),
        group_by = u$group_by, metric = met,
        test = r$test %||% NA, p_global = r$p_global %||% NA_real_,
        skipped = r$skipped %||% NA_character_)
      if (!is.null(r$posthoc) && nrow(r$posthoc))
        ph[[length(ph)+1]] <- tibble::tibble(
          marker = marker, goal = paste0(gname,"_",spec$stats_extra),
          context = paste0(ctx_lab," | ",split_col,"=",lv), metric = met, r$posthoc)
    }
  }
  list(rows = rows, ph = ph)
}

# write the stats tables to CSV
write_alpha_stats <- function(stats_list, marker,
                              out_dir = "results/alpha") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  g <- file.path(out_dir, paste0("alpha_stats_", tolower(marker), ".csv"))
  p <- file.path(out_dir, paste0("alpha_posthoc_", tolower(marker), ".csv"))
  readr::write_csv(stats_list$global, g)
  readr::write_csv(stats_list$posthoc, p)
  message("Wrote ", g, " and ", p)
  c(g, p)
}
