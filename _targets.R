# =============================================================================
# _targets.R  —  Palm Oil Soil Metagenomics pipeline (gold-tier)
#   data/source/*.xlsx + wf_*.html -> bronze -> silver -> gold -> QC/normalize
#   -> Dashboard + Analysis output (alpha/beta/relabund/DA), synchronized with
#   the originally-established Goal A/B/C/D + data_counts structure.
# =============================================================================

library(targets)
library(tarchetypes)

tar_option_set(
  # NOTE: ANCOMBC deliberately excluded here -- tar_option_set's packages
  # are library()'d before ANY target runs, so listing an uninstalled
  # package here would break the whole pipeline, not just the DA target.
  # functions_gold_da.R calls it namespaced (ANCOMBC::ancombc2()), which
  # only requires the package to exist when that specific target runs.
  packages = c("tidyverse", "yaml", "vegan", "patchwork", "permute",
               "igraph", "ggraph", "pheatmap", "ggrepel", "readxl",
               "rvest", "xml2", "phyloseq"),
  format   = "rds"     # default for all targets
)

# load all function files
lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

list(

  # ===== STAGE 1 — metadata_prep (TASK 1: source xlsx -> bronze -> silver) ==
  tar_target(source_metadata_files,
             list.files("data/source", pattern = "_metadata\\.xlsx$",
                        full.names = TRUE),
             format = "file"),

  tar_target(source_metadata_raw, compile_source_metadata(source_metadata_files)),

  tar_target(bronze_metadata, build_bronze_metadata(source_metadata_raw)),
  tar_target(bronze_metadata_csv,
             write_bronze_metadata(bronze_metadata,
                                   "data/bronze/bronze_metadata.csv"),
             format = "file"),

  tar_target(silver_metadata, build_silver_metadata(bronze_metadata)),
  tar_target(silver_metadata_16s_csv,
             write_silver_metadata(silver_metadata$metadata_16s,
                                   "data/silver/silver_metadata_16s.csv"),
             format = "file"),
  tar_target(silver_metadata_its_csv,
             write_silver_metadata(silver_metadata$metadata_its,
                                   "data/silver/silver_metadata_its.csv"),
             format = "file"),

  # ===== STAGE 2 — abundance_matrix_prep (TASK 2: wf html -> bronze -> silver)
  tar_target(wf_html_16s_files,
             list.files("data/source", pattern = "^wf_16s_batch_[0-9]+\\.html$",
                        full.names = TRUE),
             format = "file"),
  tar_target(wf_html_its_files,
             list.files("data/source", pattern = "^wf_its_batch_[0-9]+\\.html$",
                        full.names = TRUE),
             format = "file"),

  tar_target(wf_reports_16s, purrr::map(wf_html_16s_files, parse_wf_report)),
  tar_target(wf_reports_its, purrr::map(wf_html_its_files, parse_wf_report)),

  # ---- bronze: wf run stats ------------------------------------------------
  tar_target(bronze_wf_runs_16s, compile_wf_runs(wf_reports_16s)),
  tar_target(bronze_wf_runs_its, compile_wf_runs(wf_reports_its)),
  tar_target(bronze_wf_runs_16s_csv,
             write_csv_out(bronze_wf_runs_16s, "data/bronze/bronze_wf_runs_16s.csv"),
             format = "file"),
  tar_target(bronze_wf_runs_its_csv,
             write_csv_out(bronze_wf_runs_its, "data/bronze/bronze_wf_runs_its.csv"),
             format = "file"),

  # ---- bronze: abundance matrices (superkingdom/phylum/genus/species) -----
  tar_target(bronze_matrix_16s, compile_bronze_matrices(wf_reports_16s)),
  tar_target(bronze_matrix_its, compile_bronze_matrices(wf_reports_its)),

  tar_target(bronze_matrix_16s_superkingdom_csv,
             write_csv_out(bronze_matrix_16s$superkingdom,
                           "data/bronze/bronze_matrix_16s_superkingdom.csv"),
             format = "file"),
  tar_target(bronze_matrix_16s_phylum_csv,
             write_csv_out(bronze_matrix_16s$phylum,
                           "data/bronze/bronze_matrix_16s_phylum.csv"),
             format = "file"),
  tar_target(bronze_matrix_16s_genus_csv,
             write_csv_out(bronze_matrix_16s$genus,
                           "data/bronze/bronze_matrix_16s_genus.csv"),
             format = "file"),
  tar_target(bronze_matrix_16s_species_csv,
             write_csv_out(bronze_matrix_16s$species,
                           "data/bronze/bronze_matrix_16s_species.csv"),
             format = "file"),

  tar_target(bronze_matrix_its_superkingdom_csv,
             write_csv_out(bronze_matrix_its$superkingdom,
                           "data/bronze/bronze_matrix_its_superkingdom.csv"),
             format = "file"),
  tar_target(bronze_matrix_its_phylum_csv,
             write_csv_out(bronze_matrix_its$phylum,
                           "data/bronze/bronze_matrix_its_phylum.csv"),
             format = "file"),
  tar_target(bronze_matrix_its_genus_csv,
             write_csv_out(bronze_matrix_its$genus,
                           "data/bronze/bronze_matrix_its_genus.csv"),
             format = "file"),
  tar_target(bronze_matrix_its_species_csv,
             write_csv_out(bronze_matrix_its$species,
                           "data/bronze/bronze_matrix_its_species.csv"),
             format = "file"),

  # ---- silver: dedup wf runs, sync matrices to survivors --------------------
  tar_target(silver_wf_runs_16s, dedup_wf_runs(bronze_wf_runs_16s)),
  tar_target(silver_wf_runs_its, dedup_wf_runs(bronze_wf_runs_its)),
  tar_target(silver_wf_runs_16s_csv,
             write_csv_out(silver_wf_runs_16s, "data/silver/silver_wf_runs_16s.csv"),
             format = "file"),
  tar_target(silver_wf_runs_its_csv,
             write_csv_out(silver_wf_runs_its, "data/silver/silver_wf_runs_its.csv"),
             format = "file"),

  tar_target(silver_matrix_16s,
             purrr::map(bronze_matrix_16s, sync_matrix_to_silver, silver_wf_runs_16s)),
  tar_target(silver_matrix_its,
             purrr::map(bronze_matrix_its, sync_matrix_to_silver, silver_wf_runs_its)),

  tar_target(silver_matrix_16s_superkingdom_csv,
             write_csv_out(silver_matrix_16s$superkingdom,
                           "data/silver/silver_matrix_16s_superkingdom.csv"),
             format = "file"),
  tar_target(silver_matrix_16s_phylum_csv,
             write_csv_out(silver_matrix_16s$phylum,
                           "data/silver/silver_matrix_16s_phylum.csv"),
             format = "file"),
  tar_target(silver_matrix_16s_genus_csv,
             write_csv_out(silver_matrix_16s$genus,
                           "data/silver/silver_matrix_16s_genus.csv"),
             format = "file"),
  tar_target(silver_matrix_16s_species_csv,
             write_csv_out(silver_matrix_16s$species,
                           "data/silver/silver_matrix_16s_species.csv"),
             format = "file"),

  tar_target(silver_matrix_its_superkingdom_csv,
             write_csv_out(silver_matrix_its$superkingdom,
                           "data/silver/silver_matrix_its_superkingdom.csv"),
             format = "file"),
  tar_target(silver_matrix_its_phylum_csv,
             write_csv_out(silver_matrix_its$phylum,
                           "data/silver/silver_matrix_its_phylum.csv"),
             format = "file"),
  tar_target(silver_matrix_its_genus_csv,
             write_csv_out(silver_matrix_its$genus,
                           "data/silver/silver_matrix_its_genus.csv"),
             format = "file"),
  tar_target(silver_matrix_its_species_csv,
             write_csv_out(silver_matrix_its$species,
                           "data/silver/silver_matrix_its_species.csv"),
             format = "file"),

  # ===== STAGE 3 — metadata_matrix (TASK 3: gold layer) =====================
  tar_target(gold_thresholds_file, "config/gold_qc_thresholds.yaml", format = "file"),
  tar_target(gold_thresholds, load_gold_thresholds(gold_thresholds_file)),

  tar_target(gold_qc_16s, build_gold_qc(silver_metadata$metadata_16s,
                                        silver_wf_runs_16s,
                                        gold_thresholds$reads_cutoff[["16S"]])),
  tar_target(gold_qc_its, build_gold_qc(silver_metadata$metadata_its,
                                        silver_wf_runs_its,
                                        gold_thresholds$reads_cutoff[["ITS"]])),
  tar_target(gold_qc_16s_csv,
             write_csv_out(gold_qc_16s, "data/gold/gold_qc_16s.csv"),
             format = "file"),
  tar_target(gold_qc_its_csv,
             write_csv_out(gold_qc_its, "data/gold/gold_qc_its.csv"),
             format = "file"),

  tar_target(recommendation_16s, build_recommendation(gold_qc_16s)),
  tar_target(recommendation_its, build_recommendation(gold_qc_its)),
  tar_target(recommendation_16s_csv,
             write_csv_out(recommendation_16s, "data/gold/recommendation_16s.csv"),
             format = "file"),
  tar_target(recommendation_its_csv,
             write_csv_out(recommendation_its, "data/gold/recommendation_its.csv"),
             format = "file"),

  tar_target(gold_matrix_16s,
             purrr::map(silver_matrix_16s, gold_matrix_high_pass, gold_qc_16s)),
  tar_target(gold_matrix_its,
             purrr::map(silver_matrix_its, gold_matrix_high_pass, gold_qc_its)),

  tar_target(gold_matrix_16s_superkingdom_csv,
             write_csv_out(gold_matrix_16s$superkingdom,
                           "data/gold/gold_matrix_16s_superkingdom.csv"),
             format = "file"),
  tar_target(gold_matrix_16s_phylum_csv,
             write_csv_out(gold_matrix_16s$phylum,
                           "data/gold/gold_matrix_16s_phylum.csv"),
             format = "file"),
  tar_target(gold_matrix_16s_genus_csv,
             write_csv_out(gold_matrix_16s$genus,
                           "data/gold/gold_matrix_16s_genus.csv"),
             format = "file"),
  tar_target(gold_matrix_16s_species_csv,
             write_csv_out(gold_matrix_16s$species,
                           "data/gold/gold_matrix_16s_species.csv"),
             format = "file"),

  tar_target(gold_matrix_its_superkingdom_csv,
             write_csv_out(gold_matrix_its$superkingdom,
                           "data/gold/gold_matrix_its_superkingdom.csv"),
             format = "file"),
  tar_target(gold_matrix_its_phylum_csv,
             write_csv_out(gold_matrix_its$phylum,
                           "data/gold/gold_matrix_its_phylum.csv"),
             format = "file"),
  tar_target(gold_matrix_its_genus_csv,
             write_csv_out(gold_matrix_its$genus,
                           "data/gold/gold_matrix_its_genus.csv"),
             format = "file"),
  tar_target(gold_matrix_its_species_csv,
             write_csv_out(gold_matrix_its$species,
                           "data/gold/gold_matrix_its_species.csv"),
             format = "file"),

  # ===== QC — gold -> processed matrix (taxonomic + threshold filter) =======
  #   Independent at genus and species rank (both are DA/co-occurrence
  #   targets downstream); pre-rarefaction (ANCOM-BC wants raw counts).
  tar_target(gold_processing_thresholds_file,
             "config/gold_processing_thresholds.yaml", format = "file"),
  tar_target(gold_processing_thresholds,
             load_gold_processing_thresholds(gold_processing_thresholds_file)),

  tar_target(gold_processed_matrix_16s_genus,
             build_gold_processed_matrix(gold_matrix_16s$genus, "16S", "genus",
                                         gold_processing_thresholds)),
  tar_target(gold_processed_matrix_16s_species,
             build_gold_processed_matrix(gold_matrix_16s$species, "16S", "species",
                                         gold_processing_thresholds)),
  tar_target(gold_processed_matrix_its_genus,
             build_gold_processed_matrix(gold_matrix_its$genus, "ITS", "genus",
                                         gold_processing_thresholds)),
  tar_target(gold_processed_matrix_its_species,
             build_gold_processed_matrix(gold_matrix_its$species, "ITS", "species",
                                         gold_processing_thresholds)),

  tar_target(gold_processed_matrix_16s_genus_csv,
             write_csv_out(gold_processed_matrix_16s_genus,
                           "data/gold/gold_processed_matrix_16s_genus.csv"),
             format = "file"),
  tar_target(gold_processed_matrix_16s_species_csv,
             write_csv_out(gold_processed_matrix_16s_species,
                           "data/gold/gold_processed_matrix_16s_species.csv"),
             format = "file"),
  tar_target(gold_processed_matrix_its_genus_csv,
             write_csv_out(gold_processed_matrix_its_genus,
                           "data/gold/gold_processed_matrix_its_genus.csv"),
             format = "file"),
  tar_target(gold_processed_matrix_its_species_csv,
             write_csv_out(gold_processed_matrix_its_species,
                           "data/gold/gold_processed_matrix_its_species.csv"),
             format = "file"),

  # ===== NORMALIZE — rarefied species matrix, per universe (marker x stage) =
  tar_target(gold_universe_lookup_16s, build_universe_lookup(gold_qc_16s)),
  tar_target(gold_universe_lookup_its, build_universe_lookup(gold_qc_its)),

  tar_target(gold_rarefied_species_16s,
             rarefy_species_by_universe(gold_processed_matrix_16s_species,
                                        gold_universe_lookup_16s)),
  tar_target(gold_rarefied_species_its,
             rarefy_species_by_universe(gold_processed_matrix_its_species,
                                        gold_universe_lookup_its)),

  # ===== DASHBOARD ============================================================
  tar_target(plot_style_file, "config/plot_style.yaml", format = "file"),
  tar_target(plot_style, load_plot_style(plot_style_file)),

  tar_target(dashboard_qc_status_pie,
             build_qc_status_pie(gold_qc_16s, gold_qc_its, plot_style),
             format = "file"),
  tar_target(dashboard_prepost_depth,
             build_prepost_depth_plot(gold_qc_16s, gold_qc_its,
                                      gold_processed_matrix_16s_species,
                                      gold_processed_matrix_its_species,
                                      plot_style),
             format = "file"),
  tar_target(dashboard_rarefaction_curves,
             build_rarefaction_curves(gold_matrix_16s$species,
                                      gold_universe_lookup_16s,
                                      gold_rarefied_species_16s,
                                      gold_matrix_its$species,
                                      gold_universe_lookup_its,
                                      gold_rarefied_species_its),
             format = "file"),
  tar_target(dashboard_pcoa,
             build_pcoa_plots(gold_rarefied_species_16s, gold_universe_lookup_16s,
                              gold_rarefied_species_its, gold_universe_lookup_its,
                              plot_style),
             format = "file"),
  tar_target(dashboard_shannon_barplot,
             build_shannon_dashboard_barplot(alpha_shannon, plot_style),
             format = "file"),
  tar_target(dashboard_xlsx,
             build_dashboard_xlsx(gold_universe_lookup_16s, gold_universe_lookup_its),
             format = "file"),

  # ===== ANALYSIS OUTPUT ======================================================
  tar_target(analysis_thresholds_file, "config/analysis_thresholds.yaml",
             format = "file"),
  tar_target(analysis_thresholds, load_analysis_thresholds(analysis_thresholds_file)),

  # ---- alpha diversity (Shannon only, rarefied species matrix) -------------
  tar_target(alpha_shannon,
             build_alpha_shannon(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                 gold_rarefied_species_its, gold_universe_lookup_its)),
  tar_target(alpha_shannon_plot,
             plot_alpha_shannon(alpha_shannon, plot_style), format = "file"),
  tar_target(alpha_shannon_csv,
             write_alpha_shannon_csv(alpha_shannon), format = "file"),

  # ---- beta diversity (Bray-Curtis only, rarefied species matrix) ----------
  tar_target(beta_analysis,
             build_beta_analysis(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                 gold_rarefied_species_its, gold_universe_lookup_its,
                                 plot_style)),
  tar_target(beta_analysis_plots, beta_analysis$plots, format = "file"),
  tar_target(beta_analysis_permanova_csv, beta_analysis$permanova_csv, format = "file"),

  # ---- data_counts (synchronized from established pipeline) ----------------
  tar_target(gold_count_tables,
             build_gold_count_tables(gold_universe_lookup_16s, gold_universe_lookup_its),
             format = "file"),

  # ---- Goal A/B/C/D alpha report tree (synchronized, Shannon-only) ---------
  #   `comparisons` is declared in the established-pipeline Stage 4 block
  #   below; targets resolves the dependency by reference, not list order.
  tar_target(gold_report_tree,
             build_gold_report_tree(alpha_shannon, comparisons, plot_style),
             format = "file"),

  # ---- Goal B/D beta breakdown (synchronized, Bray-Curtis only) -----------
  tar_target(beta_goal_tree,
             build_gold_beta_goal_tree(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                       gold_rarefied_species_its, gold_universe_lookup_its,
                                       plot_style)),
  tar_target(beta_goal_tree_plots, beta_goal_tree$plots, format = "file"),
  tar_target(beta_goal_tree_permanova_csv, beta_goal_tree$permanova_csv, format = "file"),

  # ---- Goal B/D relative abundance stacked bars (synchronized) -------------
  #   species is the base rank; genus/phylum roll up from it.
  tar_target(gold_relabund_tree,
             build_gold_relabund_tree(gold_rarefied_species_16s, gold_universe_lookup_16s,
                                      gold_processed_matrix_16s_species,
                                      gold_rarefied_species_its, gold_universe_lookup_its,
                                      gold_processed_matrix_its_species,
                                      plot_style),
             format = "file"),

  # ---- DA (ANCOM-BC): Case 1 (fertilizer within field) + Case 2 (waktu ----
  #      within field x fertilizer), genus + species, + Case 2 LFC dumbbells
  tar_target(gold_da_analysis,
             build_gold_da_analysis(gold_processed_matrix_16s_genus,
                                    gold_processed_matrix_16s_species,
                                    gold_universe_lookup_16s,
                                    gold_processed_matrix_its_genus,
                                    gold_processed_matrix_its_species,
                                    gold_universe_lookup_its,
                                    analysis_thresholds, plot_style)),
  tar_target(gold_da_analysis_csv, gold_da_analysis$csv, format = "file"),
  tar_target(gold_da_analysis_plots, gold_da_analysis$plots, format = "file"),

  # ---- shared config: comparisons.yaml (Goal A/B/C/D definitions) ---------
  #   Read by gold_report_tree above.
  #
  #   The REST of the originally-established data/raw-based pipeline (Stage 2
  #   ingest/validation, Stage 3 taxonomic/threshold filter + rarefy+CLR,
  #   Stage 4 alpha/beta/relabund/count report trees, Stage 5 co-occurrence
  #   networks, raw exports) has been fully synchronized onto the gold-tier
  #   pipeline above and its targets removed from this DAG: data/raw/ no
  #   longer exists in this project (superseded by data/source -> bronze ->
  #   silver -> gold), so list.files() against it always returned character(0),
  #   which propagated to an empty metadata_all/long_counts and a hard
  #   validate_data() stop() -- that aborted the ENTIRE tar_make() run, even
  #   though every gold-tier target above it was otherwise fine on its own.
  #   The original implementations are untouched and still callable directly
  #   (functions_ingest.R, functions_qc.R, functions_qc_apply.R, and the
  #   build_report_tree()/build_beta_tree()/build_relabund_tree()/
  #   build_count_tables()/build_network_tree() drivers) -- they're just no
  #   longer wired into the DAG. Their still-shared building blocks (GOAL_DIR,
  #   fert_palette(), plot_ordination(), plot_stacked_bar(),
  #   render_count_png(), universe_alpha_stats(), etc.) remain in active use
  #   by the gold-tier targets above.
  tar_target(comparisons_file, "config/comparisons.yaml", format = "file"),
  tar_target(comparisons, load_comparisons(comparisons_file))
)
