# Palm Oil Soil Metagenomics Pipeline

Reproducible `targets` pipeline for profiling soil microbiome under different
fertilizer treatments across multiple corporate fields (kebun), for two
amplicon markers (16S, ITS) and two plant stages (TM, Nursery).

The data is treated as **four isolated universes** — `16S_TM`, `16S_Nursery`,
`ITS_TM`, `ITS_Nursery` — which are never compared with one another.

## Setup

```bash
conda env create -f environment_r.yml      # R base + system libs
conda activate palmoil
R -e 'renv::restore()'                      # install exact package versions
```

## Run

```r
targets::tar_make()        # build everything (only stale targets rebuild)
targets::tar_visnetwork()  # view the dependency graph
```

New data = drop the metadata workbook and EPI2ME HTML reports into
`data/source/` and re-run `tar_make()`.
Nothing is merged by hand; the pipeline compiles by barcode string.

### Select analysis cohorts

Gold QC and recommendation files are always built for all cohorts. To limit
the downstream processing, dashboards, and analysis, edit
`config/analysis_cohorts.yaml`:

```yaml
selected: all                 # default
# selected: 16S_TM            # one cohort
# selected: [16S_TM, ITS_TM]  # multiple cohorts
```

Supported cohorts are `16S_TM`, `ITS_TM`, and `16S_Nursery`. Changing this
setting and running `targets::tar_make()` invalidates only the downstream
branch; the complete bronze, silver, and gold-QC layers remain unchanged. The
resolved selection is recorded in `Results/analysis/selected_cohorts.csv`.

## Project layout

```
config/        schema.yaml, thresholds.yaml, comparisons.yaml  (the data contract)
R/             pipeline functions (ingest, QC, alpha, beta, relabund, counts)
_targets.R     the DAG
docs/          DATA_CONTRACT.md  (rules for incoming data deliveries)
data/raw/      immutable input matrices + metadata (read-only, never edited)
Results/       generated outputs (see below)
```

## Pipeline stages

1. **Ingest** — per-batch TSVs compiled and joined to metadata on the
   globally-unique barcode string (batch number ignored). Strict validation.
2. **QC** — taxonomic filter (Bacteria/Archaea for 16S; Fungi for ITS;
   organelles removed), sample-depth + OTU + 5% prevalence thresholds.
3. **Normalize** — synchronized rarefied + CLR tables (identical sample/OTU set).
4. **Analysis**, per universe and per goal:
   - Alpha (Observed, Pielou, Shannon): unpaired Kruskal-Wallis / Wilcoxon,
     n>=2 gate, BH-adjusted.
   - Beta (Aitchison/Bray-Curtis/Jaccard): ordinations, Ward.D2 dendrograms,
     PERMANOVA. Generated for Goals B and D.
   - Relative abundance: Top-10 / Top-15 stacked bars (phylum + genus) +
     full CSVs.
   - Temporal differential abundance: ANCOM-BC2 at genus and species levels,
     one model per cohort and kebun with fertilizers pooled; exports only
     consecutive contrasts (T0 vs T1, T1 vs T2, ...) as CSVs and dumbbell plots.
   - FAPROTAX functional inference (16S only): normalized species profiles are
     mapped to putative prokaryotic functions, then consecutive timepoints are
     compared per cohort and kebun with fertilizers pooled. Wilcoxon tests use
     BH correction. A cohort-level effect heatmap summarizes direction,
     magnitude, and significance across kebun; selected-function trajectory
     panels show sample distributions and means through time; pairwise
     dumbbells retain detailed comparisons.
   - Replicate-count tables per field.

## The four analysis goals

| Goal | Audience | Question |
|------|----------|----------|
| A | Palm Oil Consultant | Fertilizers within one field at one timepoint |
| B | Palm Oil Consultant | T0->T1 change within one field |
| C | Fertilizer Consultant | Fertilizers across fields at one timepoint |
| D | Fertilizer Consultant | T0->T1 change across all fields |

## Output structure

```
Results/
├── data_counts/<universe>/<field>_counts.png   replicate inventory tables
├── dropped_barcodes.csv / filtered_barcodes_all.csv   QC audit trail
└── <universe>/
    ├── Goal_A_Intra_Snapshot/<field>/      alpha boxplots + sliced stats
    ├── Goal_B_Intra_Longitudinal/<field>/  alpha trajectories, beta, stacked bars
    ├── Goal_C_Cross_Snapshot/              pooled alpha snapshot
    ├── Goal_D_Cross_Longitudinal/          alpha trajectories, beta, stacked bars
    └── Differential_Abundance_Temporal/<field>/<T0_vs_T1>/<rank>/
                                             ANCOM-BC2 CSVs + dumbbell plot
    └── Functional_FAPROTAX/<field>/<T0_vs_T1>/
                                             functional CSVs + dumbbell plot
```

FAPROTAX annotation audit files are written to `data/gold/faprotax/`.
Cross-kebun run summaries and per-sample annotation coverage are written to
`Results/analysis/faprotax_temporal_summary.csv` and
`Results/analysis/faprotax_annotation_coverage.csv`. Each cohort's
`Functional_FAPROTAX/` root also contains `faprotax_effect_heatmap.png` and
`faprotax_trajectories.png`; the coverage QC plot is in `Results/analysis/`.
The same three visualizations are also generated inside every kebun directory,
using functions ranked from that kebun alone:
`faprotax_effect_heatmap.png`, `faprotax_trajectories.png`, and
`faprotax_annotation_coverage.png`.

## Important notes

- **Status: interim** — currently T0 and T1 only (study ongoing, up to T4).
- **Batch confound:** sequencing/extraction batch is confounded with timepoint.
  Within-timepoint comparisons (A, C) are clean; temporal comparisons (B, D)
  are reported descriptively with this caveat. ComBat is not applicable. ITS is
  more affected than 16S (fungal lysis sensitivity).
- **Replication:** 1-2 biological replicates per group; analyses are largely
  descriptive. Stats run only where a group has n>=2.
- **ITS is TM-only** — `ITS_Nursery` is skipped (no data).
- Stats use **raw p** for plot significance markers; CSVs carry raw + BH-adjusted.

## Reproducibility

Conda pins R + system libraries; `renv.lock` pins exact R package versions;
`targets` caches results and rebuilds only what changed. Rarefaction seed = 42.
