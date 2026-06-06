library(dplyr)
library(readr)
library(readxl)
metadata_all <- read_excel("data/raw/metadata.xlsx",sheet="metadata")
dir.create("metadata_split", showWarnings = FALSE)

batches <- sort(unique(metadata_all$Batch_Ekstraksi))

# 16S
for (batch in batches) {

  df <- metadata_all %>%
    filter(`16S` == 1, Batch_Ekstraksi == batch) %>%
    select(
      -Barcode_ITS,
      -ITS,
      -Hasil_Sekuensing_ITS,
      -Flag_ITS,
      -`16S`
    )

  if (nrow(df) > 0) {
    write_tsv(
      df,
      file.path(
        "metadata_split",
        paste0("metadata_16s_batch_", batch, ".tsv")
      )
    )
  }
}

# ITS
for (batch in batches) {

  df <- metadata_all %>%
    filter(ITS == 1, Batch_Ekstraksi == batch) %>%
    select(
      -Barcode_16S,
      -`16S`,
      -Hasil_Sekuensing_16S,
      -Flag_16S,
      -ITS
    )

  if (nrow(df) > 0) {
    write_tsv(
      df,
      file.path(
        "metadata_split",
        paste0("metadata_its_batch_", batch, ".tsv")
      )
    )
  }
}