library(DiffBind)
library(tidyverse)
library(cowplot)
samples <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/chip_seq_project/sample_metadata.csv")
dbObj <- dba(sampleSheet=samples)
dbObj
# 3 Samples, 321 sites in matrix (8581 total):
#   ID Replicate Caller Intervals
# 1 8-1         1 narrow       410
# 2 8-2         1 narrow       172
# 3 8-3         1 narrow      8489
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj
# FRIP score
# The FRiP score is defined as the fraction of reads that fall into a peak 
# and is often used as a measure of ChIP-seq quality.
# 3 Samples, 321 sites in matrix:
#   ID Replicate Caller Intervals FRiP
# 1 8-1         1 counts       321 0.02
# 2 8-2         1 counts       321 0.01
# 3 8-3         1 counts       321 0.04
dba.plotPCA(dbObj, label=DBA_ID)
plot(dbObj)
dbObj <- dba.contrast(dbObj, minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# Occupancy analysis
dba.show(dbObj, bContrast=T)

dba.plotVenn(dbObj$peaks, )
