install.packages("targets")
install.packages("tarchetypes")
install.packages("tidyverse")
install.packages("yaml")
install.packages("rmarkdown")
install.packages("vegan")
install.packages("BH")
install.packages("BiocManager")
install.packages("Cairo")
install.packages("ggrastr")
install.packages("remotes")
remotes::install_version(
  "rbiom",
  version = "2.2.1",
  repos = "https://cran.r-project.org",
  upgrade = "never"
)
BiocManager::install("XVector")
BiocManager::install("Biostrings")
BiocManager::install("scater")
BiocManager::install("phyloseq")
BiocManager::install("mia")
renv::snapshot() 