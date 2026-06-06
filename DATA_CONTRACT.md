# DATA CONTRACT — Palm Oil Soil Metagenomics

**Purpose.** Every data delivery (each new sequencing/extraction batch, each new
timepoint) must follow this structure so the analysis pipeline ingests it
automatically with no manual editing. If a delivery breaks these rules, the
pipeline halts at validation rather than producing wrong results.

Hand this document to whoever prepares metadata and exports matrices.

---

## 1. File organisation

Raw deliveries go under `data/raw/`, **named by batch, never edited after drop:**

```
data/raw/
├── metadata.xlsx                      # master (optional, human reference)
├── metadata_16s/metadata_16s_batch_<N>.tsv
├── metadata_its/metadata_its_batch_<N>.tsv
├── raw_mat_16s/wf_16s_batch_<N>_genus.tsv
├── raw_mat_16s/wf_16s_batch_<N>_phylum.tsv
├── raw_mat_its/wf_its_batch_<N>_genus.tsv
└── raw_mat_its/wf_its_batch_<N>_phylum.tsv
```

- `<N>` = sequencing batch number. **Filename is the only place batch number
  lives for matrices**, so it must be correct.
- A batch may contain **any mix** of timepoints and stages — that is fine.
- 16S and ITS need not have the same set of batches.

## 2. Matrix files (wide format, tab-separated)

- Rows = taxa. Columns = one per barcode, named exactly like `barcode01_97`.
- Must also contain taxonomy columns: `superkingdom, kingdom, phylum, class,
  order, family, genus` (genus file) and `tax`, plus `total`.
- **Barcodes reset every batch** (batch 1 and batch 2 both have barcode01).
  The barcode is only unique *within* its batch.

## 3. Metadata files (tab-separated, one row per sample per batch)

Required columns (exact Indonesian names — do not rename, do not reorder):

| Column | Meaning | Allowed values |
|---|---|---|
| `ID_Sampel` | biological sample id | e.g. TPP-OL0002 |
| `Jenis_Kebun` | plant stage | **`TM`** or **`Nursery`** only |
| `Kelompok_Kebun` | field / kebun | code (TPP, PSG, MJ, AMR, …) |
| `Kelompok_Pupuk` | fertilizer | **TM:** `PH1`–`PH6`  •  **Nursery:** `PHN01`–`PHN12` |
| `Fase_Treatment` | timepoint | `T0`, `T1`, `T2`, … |
| `Barcode_16S` | 16S barcode | must match a 16S matrix column header **for that batch**, or `NA` |
| `Barcode_ITS` | ITS barcode | must match an ITS matrix column header **for that batch**, or `NA` |
| `Batch_Ekstraksi` | extraction batch | integer |
| `Flag_16S` | 16S QC flag | `0` = keep, `1` = drop |
| `Flag_ITS` | ITS QC flag | `0` = keep, `1` = drop |

> **Critical:** fertilizer levels are different between stages. `TM` uses
> `PH1–PH6`; `Nursery` uses `PHN##`. Never put a Nursery sample under a `PH#`
> code or vice-versa.

## 4. The rules the pipeline enforces (validation will FAIL if broken)

1. Every metadata column above is present and spelled exactly as shown.
2. `Jenis_Kebun` ∈ {TM, Nursery}. No typos, no trailing spaces, no other values.
3. `Kelompok_Pupuk` matches its stage pattern (PH# for TM, PHN## for Nursery).
4. `Fase_Treatment` matches `T<number>`.
5. Each non-NA barcode in metadata matches a real column header in the matrix
   **of the same batch and marker**. A batch that matches **zero** columns
   halts the run (catches barcode-suffix mismatches early).
6. No duplicate `ID_Sampel` within the same marker + batch. A deliberate
   re-sequence must get a distinct id (e.g. `TPP-OL0002_rerun`).
7. `Batch_Ekstraksi` and `Flag` columns are populated (no blanks).

## 5. QC

A sample is dropped automatically if its marker `Flag == 1`. The
`Hasil_Sekuensing_*` text columns are kept for reference but not used for
filtering.

## 6. What you may add freely without code changes

- New batches (just drop the files).
- New timepoints (T2, T3, …) — anywhere, any batch.
- New fields / new fertilizer levels (PHN13, etc.) — recognised automatically.

## 7. What you must NOT do

- Do not edit or overwrite a previously delivered file.
- Do not rename columns or change `TM`/`Nursery`/`PH#` spellings.
- Do not merge batches into one file manually.
- Do not reuse an `ID_Sampel` for a different physical sample.
