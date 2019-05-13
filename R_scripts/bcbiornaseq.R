if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("remotes")
BiocManager::install("hbc/bcbioRNASeq")

library(bcbioRNASeq)