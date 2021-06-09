## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")

setwd("/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
library(clusterProfiler)
files <-  list(peak1 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample1/macs2_sample1__peaks.xls",
               peak3 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample3/macs2__peaks.xls",
               peak4 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample4/macs2__peaks.xls",
               peak6 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample6/macs2__peaks.xls",
               peak7 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample7/macs2__peaks.xls",
               peak8 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample8/macs2__peaks.xls",
               peak9 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample9/macs2__peaks.xls",
               peak11 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample11/macs2__peaks.xls",
               peak12 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample12/macs2__peaks.xls",
               peak13 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample13/macs2__peaks.xls",
               peak14 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample14/macs2__peaks.xls",
               peak16 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample16/macs2__peaks.xls",
               peak17 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample17/macs2__peaks.xls",
               peak18 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample18/macs2__peaks.xls",
               peak19 = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample19/macs2__peaks.xls")
files

peak4 = readPeakFile(files[[3]])
peak4

peak4 = readPeakFile("/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample4/macs2__peaks.xls")
peak4_cov_q_value_plot <- covplot(peak4, weightCol="X.log10.qvalue.", xlim=c(4.5e7, 5e7))
peak4_cov_foldenrichment_plot <- covplot(peak4, weightCol="fold_enrichment", xlim=c(4.5e7, 5e7))
# To look at specific chromosomes
# covplot(peak4, weightCol="X.log10.qvalue.", chrs=c("chr17", "chr18"), xlim=c(5.5e7, 6e7))


peak14 = readPeakFile(files[[11]])
peak14
covplot(peak14, weightCol="X.log10.qvalue.", xlim=c(4.5e7, 5e7))
covplot(peak14, weightCol="fold_enrichment", xlim=c(4.5e7, 5e7))

peakAnno_peak4 <- annotatePeak(peak4,
                         TxDb=txdb, annoDb="org.Hs.eg.db")



## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library("ChIPpeakAnno")

macs_peak4 <- system.file("extdata", "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample4/macs2__peaks.xls", package="ChIPpeakAnno")
macs_peak4 <- system.file("extdata", "macs2__peaks.xls", package="ChIPpeakAnno", mustWork = TRUE)
macs_peak4_Output <- toGRanges(peak4, format="MACS")

