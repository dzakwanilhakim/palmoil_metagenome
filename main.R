library(targets)
targets::tar_make()

# the matched samples that survived
targets::tar_read(long_qc) |> dplyr::distinct(marker, id_sampel) |> dplyr::count(marker)

# the audit trail of everything dropped
readr::read_csv(targets::tar_read(drop_report))

net <- tar_visnetwork()
htmlwidgets::saveWidget(net, "pipeline.html")

# how many samples matched per marker?
targets::tar_read(long_qc) |> dplyr::distinct(marker, id_sampel) |> dplyr::count(marker)

# do the stage / fertilizer / timepoint distributions look right?
targets::tar_read(long_qc) |>
  dplyr::distinct(marker, id_sampel, stage, fertilizer, timepoint) |>
  dplyr::count(marker, stage, timepoint)

# what got dropped
readr::read_csv(targets::tar_read(drop_report)) |> dplyr::count(marker, side)

targets::tar_make()

source("R/functions_qc.R")
depth <- targets::tar_read(preqc)$depth_table

# ITS survivors at 10k, broken down by the groups your analysis will use
depth |>
  dplyr::filter(marker == "ITS", depth >= 10000) |>
  dplyr::count(stage)

targets::tar_read(tax_filtered) |>
  dplyr::filter(marker == "ITS") |>
  dplyr::group_by(id_sampel, stage, fertilizer) |>
  dplyr::summarise(depth = sum(count), .groups="drop") |>
  dplyr::filter(depth >= 10000) |>
  dplyr::count(stage, fertilizer) |> print(n = Inf)

targets::tar_read(postqc)$summary          # post-QC depths
norm <- targets::tar_read(norm_tables)
dim(norm$`16S`$rarefied);  dim(norm$`16S`$clr)   # must be identical
dim(norm$ITS$rarefied);    dim(norm$ITS$clr)     # must be identical
identical(rownames(norm$`16S`$rarefied), rownames(norm$`16S`$clr))  # TRUE