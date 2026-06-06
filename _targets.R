# =============================================================================
# _targets.R  —  Palm Oil Soil Metagenomics pipeline
# Stage 2 implemented (ingest + validation). Stage 3+ targets appended later.
# =============================================================================

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "yaml", "arrow", "vegan"),
  format   = "rds"     # default; long_qc overridden to parquet below
)

# load all function files
lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

list(

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
                        pattern = "\\.tsv$", full.names = TRUE),
             format = "file"),

  tar_target(mat_files_its,
             list.files("data/raw/raw_mat_its",
                        pattern = "\\.tsv$", full.names = TRUE),
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

  # QC dropped-sample audit (which samples removed at Step 3 and why)
  tar_target(qc_drop_report,
             write_qc_drop_report(tax_filtered, thresholded, thresholds,
                                  out_path = "results/dropped_samples_qc.csv"),
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
                                                 seed = 42))

  # ===== STAGE 4+ (alpha/beta diversity, stackbars, A/B/C/D) ===============
  # appended next, hanging off norm_tables (rarefied -> alpha, CLR -> beta).
)
