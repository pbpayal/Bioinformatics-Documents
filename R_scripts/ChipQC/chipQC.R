## ChIP-seq analysis in R
## October 2020

library(BiocManager)
install("ChIPQC")
## Load libraries
library(ChIPQC)

## Load sample data
samples <- read.csv('/Users/pbanerjee/Documents/Payal_Scripts/R/ChipQC/sample_metadata.csv')
View(samples)

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19") 
## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: Pilot Data Yihan Peng", reportFolder="ChIPQCreport")
