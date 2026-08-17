# =============================================================================
# _targets.R  —  Palm Oil Soil Metagenomics pipeline
# Stage 2 implemented (ingest + validation). Stage 3+ targets appended later.
# =============================================================================

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "yaml", "arrow", "vegan", "patchwork", "permute",
               "igraph", "ggraph", "pheatmap", "ggrepel", "readxl",
               "rvest", "xml2"),
  format   = "rds"     # default; long_qc overridden to parquet below
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

  # ---- config -------------------------------------------------------------
  tar_target(schema_file, "config/schema.yaml", format = "file"),
  tar_target(schema,      load_schema(schema_file)),

  # ---- raw file tracking (re-runs when ANY file changes/added) ------------
  tar_target(meta_files,
             list.files(c("data/raw/metadata_16s", "data/raw/metadata_its"),
                        pattern = "\\.tsv$", full.names = TRUE),
             format = "file"),

  tar_target(mat_files_16s,
            list.files("data/raw/raw_mat_16s",
                        pattern = "_(genus|phylum)\\.tsv$", full.names = TRUE),
            format = "file"),

  tar_target(mat_files_its,
            list.files("data/raw/raw_mat_its",
                        pattern = "_(genus|phylum)\\.tsv$", full.names = TRUE),
            format = "file"),
  # ---- compile ------------------------------------------------------------
  tar_target(metadata_all, compile_metadata(meta_files, schema)),

  tar_target(long_counts,
             compile_matrices(c(mat_files_16s, mat_files_its), schema)),

  # ---- join (full string, prefix fallback) --------------------------------
  tar_target(joined_counts,
             join_with_fallback(long_counts, metadata_all, schema)),

  # ---- STRICT validation gate (halts pipeline on any violation) ----------
  tar_target(validation,
             validate_data(metadata_all, joined_counts, long_counts, schema)),

  # ---- QC drop (depends on validation passing) ----------------------------
  tar_target(long_qc, {
    validation                       # force dependency: QC only after validate
    apply_qc(joined_counts, schema)
  }, format = "parquet"),

  # ---- drop report: which barcodes did NOT intersect (both directions) ----
  tar_target(drop_report,
             write_drop_report(long_counts, metadata_all, joined_counts, schema,
                               out_path = "results/dropped_barcodes.csv"),
             format = "file"),

  # ===== STAGE 3 — Step 1: taxonomic filtering (per marker) ================
  tar_target(tax_filtered, taxonomic_filter(long_qc)),

  # ===== STAGE 3 — Step 2: pre-QC assessment (STOP for thresholds) =========
  tar_target(preqc, preqc_assess(tax_filtered)),

  # file targets so the plots are tracked and land in results/
  tar_target(preqc_boxplot,    preqc$boxplot,     format = "file"),
  tar_target(preqc_rarefaction, preqc$rarefaction, format = "file"),

  # ===== STAGE 3 — Steps 3-5 (thresholds locked in config) =================
  tar_target(thresholds_file, "config/thresholds.yaml", format = "file"),
  tar_target(thresholds,      load_thresholds(thresholds_file)),

  # Step 3: apply sample-depth + OTU + prevalence filters (per marker)
  tar_target(thresholded, apply_thresholds(tax_filtered, thresholds)),

  # Step 3/4: master synchronized sample list
  tar_target(master_samples, master_sample_list(thresholded)),

  # Unified filtered-barcode audit: every barcode removed at ANY stage + reason
  tar_target(filtered_barcode_report,
             write_filtered_barcode_report(long_counts, metadata_all, joined_counts,
                                           tax_filtered, thresholded, schema, thresholds,
                                           out_path = "results/filtered_barcodes_all.csv"),
             format = "file"),

  # Step 4: post-QC assessment
  tar_target(postqc, postqc_summary(thresholded)),
  tar_target(postqc_boxplot, plot_postqc_boxplot(postqc$depth),
             format = "file"),
  tar_target(postqc_rarefaction,
             purrr::map_chr(unique(thresholded$marker),
                            ~ plot_postqc_rarefaction(thresholded, thresholds, .x)),
             format = "file"),

  # Step 5: synchronized rarefied + CLR tables (identical samples & OTUs)
  tar_target(norm_tables, build_rarefied_and_clr(thresholded, thresholds,
                                                 seed = 42)),

  # ===== STAGE 4 — alpha diversity + A/B/C/D stats =========================
  tar_target(comparisons_file, "config/comparisons.yaml", format = "file"),
  tar_target(comparisons, load_comparisons(comparisons_file)),

  # per-sample alpha metrics (Observed, Pielou, Shannon) from rarefied
  tar_target(alpha_div, compute_alpha(norm_tables, master_samples)),

  # alpha stats per marker (A/B/C/D unpaired tests)
  tar_target(alpha_stats_16s, run_alpha_stats(alpha_div, comparisons, "16S")),
  tar_target(alpha_stats_its, run_alpha_stats(alpha_div, comparisons, "ITS")),
  tar_target(alpha_stats_16s_csv, write_alpha_stats(alpha_stats_16s, "16S"),
             format = "file"),
  tar_target(alpha_stats_its_csv, write_alpha_stats(alpha_stats_its, "ITS"),
             format = "file"),

  # alpha plots routed into the depth-first per-universe tree (Results/)
  tar_target(report_tree, build_report_tree(alpha_div, comparisons,
                                            root = "Results"),
             format = "file"),

  # beta diversity (ordinations, dendrograms, PERMANOVA) into the same tree
  tar_target(beta_tree, build_beta_tree(alpha_div, norm_tables,
                                        root = "Results"),
             format = "file"),

  # relative abundance stacked bars (Goal B per field, Goal D pooled)
  tar_target(relabund_tree, build_relabund_tree(thresholded, norm_tables,
                                                master_samples, root = "Results",
                                                top_ns = c(10, 15)),
             format = "file"),

  # replicate-count tables per field (post-QC) -> Results/data_counts/<universe>/
  tar_target(count_tables, build_count_tables(master_samples, root = "Results"),
             format = "file"),

  # species matrices (separate files) -> species long restricted to master
  tar_target(species_files,
             list.files(c("data/raw/raw_mat_16s", "data/raw/raw_mat_its"),
                        pattern = "_species\\.tsv$", full.names = TRUE),
             format = "file"),
  tar_target(species_long,
             load_species_long(species_files, metadata_all, master_samples, schema)),

  # species stacked bars (Top-15 species within Top-15 genera) Goal B + D
  tar_target(species_tree,
             build_species_tree(thresholded, norm_tables, master_samples,
                                species_long, root = "Results", top_n = 15),
             format = "file"),

  # ===== STAGE 5 — Co-occurrence networks + Active Keystone identification =
  tar_target(network_tree,
             build_network_tree(norm_tables, master_samples, species_long,
                                root = "Results"),
             format = "file"),

  # ===== RAW EXPORTS — dated aggregated CSVs (no QC/filtering) =============
  tar_target(raw_exports,
             export_raw_data(schema,
                             mat_dirs  = c("data/raw/raw_mat_16s",
                                           "data/raw/raw_mat_its"),
                             meta_dirs = c("data/raw/metadata_16s",
                                           "data/raw/metadata_its"),
                             out_dir   = "Results/raw_exports"),
             format = "file")

  # ===== STAGE 5 COMPLETE ==================================================
)
