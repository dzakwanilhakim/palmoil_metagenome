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

New data = drop the batch files into `data/raw/` and re-run `tar_make()`.
Nothing is merged by hand; the pipeline compiles by barcode string.

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
    └── Goal_D_Cross_Longitudinal/          alpha trajectories, beta, stacked bars
```

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
