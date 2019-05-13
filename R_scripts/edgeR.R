setwd("/Volumes/Drive\ 1_2/Payal/Kazue_Toru/G151")
load("G151.RData")
save.image("G151.RData")
library(limma)
library(edgeR)
#Merging HTSEQ output files for samples
files <- c("countschr_remsortedtrimmedKHT01_S1_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT02_S2_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT03_S3_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT04_S4_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT05_S5_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT06_S6_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT07_S7_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT08_S8_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT09_S9_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT10_S10_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT11_S11_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT12_S12_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT13_S13_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT14_S14_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT15_S15_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT16_S16_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT17_S17_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT18_S18_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT19_S19_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT20_S20_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT21_S21_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT22_S22_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT23_S23_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT24_S24_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT25_S25_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT26_S26_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT27_S27_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT29_S28_R1_001.fastq.gz.bam.txt","countschr_remsortedtrimmedKHT30_S29_R1_001.fastq.gz.bam.txt")

#readDGE read data from all files and merges it into one object-rawcounts
G151_counts <- readDGE(files, columns=c(1,2), header=FALSE)
G151_counts$samples
G151_counts_matrix <- as.matrix(G151_counts)

#TMM Normalization without cpm filtration
G151_normalized <- calcNormFactors(G151_counts, method="TMM")
G151_normalized$samples
G151_normalized_matrix <- as.matrix(G151_normalized)

#CPM Normalized Counts
G151_cpm_norm <- cpm(G151_normalized, normalized.lib.sizes=TRUE)

#Change the column names if required
colnames(G151_cpm_norm)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")
#colnames(G151_normCounts)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")
colnames(G151_counts_matrix)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")

#Write the required data matrix to an excel file
write.table(G151_cpm_norm, "G151_TMM_CPM.txt", sep="\t") #Normalized counts
write.table(G151_counts_matrix, "G151_raw_counts.txt", sep="\t") #Raw counts

#If we want the normalized pseudo-counts, useful for instance for cluster analysis, we can get them with the following commands.
scale <- G151_normalized$samples$lib.size*G151_normalized$samples$norm.factors
G151_normCounts <- round(t(t(G151_counts_matrix)/scale)*mean(scale))

#Boxplot
log2(G151_cpm_norm+1)
boxplot(log2(G151_cpm_norm+1), las =2, col= "royalblue2")