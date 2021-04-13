setwd("/Users/pbanerjee/Desktop")
load("edgeR_UDN_Hayk.RData")
save.image("edgeR_UDN_Hayk.RData")
library(limma)
library(edgeR)

#Father Vs Patient
files_father_vs_patient <- c("UDN287643_P_Blood_S1_counts.txt","UDN107024_UF_Blood_S2_counts.txt")
UDN_counts_FP <- readDGE(files_father_vs_patient, columns=c(1,2), header=FALSE)
UDN_counts_FP$samples
UDN_counts_FP_matrix <- as.matrix(UDN_counts_FP)
UDN_normalized_FP <- calcNormFactors(UDN_counts_FP, method="TMM")
UDN_normalized_FP$samples
UDN_normalized_FP_matrix <- as.matrix(UDN_normalized_FP)

group <- factor(c(1,2))
y <- DGEList(counts = UDN_counts_FP_matrix, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)

#Merging HTSEQ output files for samples
files <- c("UDN287643_P_Blood_S1_counts.txt","UDN107024_UF_Blood_S2_counts.txt","UDN642747_UM_Blood_S3_counts.txt")

#readDGE read data from all files and merges it into one object-rawcounts
UDN_counts <- readDGE(files, columns=c(1,2), header=FALSE)
UDN_counts$samples
UDN_counts_matrix <- as.matrix(UDN_counts)

#TMM Normalization without cpm filtration
UDN_normalized <- calcNormFactors(UDN_counts, method="TMM")
UDN_normalized$samples
UDN_normalized_matrix <- as.matrix(UDN_normalized)

#We can examine inter-sample relationships by producing a plot based on mutlidimensional scaling.
plotMDS(UDN_normalized)

#CPM Normalized Counts
UDN_cpm_norm <- cpm(UDN_normalized, normalized.lib.sizes=TRUE)

#Change the column names if required
colnames(UDN_cpm_norm)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")
#colnames(G151_normCounts)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")
colnames(UDN_counts_matrix)  <- c("KHT01_S1_R1_001","KHT02_S2_R1_001","KHT03_S3_R1_001","KHT04_S4_R1_001","KHT05_S5_R1_001","KHT06_S6_R1_001","KHT07_S7_R1_001","KHT08_S8_R1_001","KHT09_S9_R1_001","KHT10_S10_R1_001","KHT11_S11_R1_001","KHT12_S12_R1_001","KHT13_S13_R1_001","KHT14_S14_R1_001","KHT15_S15_R1_001","KHT16_S16_R1_001","KHT17_S17_R1_001","KHT18_S18_R1_001","KHT19_S19_R1_001","KHT20_S20_R1_001","KHT21_S21_R1_001","KHT22_S22_R1_001","KHT23_S23_R1_001","KHT24_S24_R1_001","KHT25_S25_R1_001","KHT26_S26_R1_001","KHT27_S27_R1_001","KHT29_S28_R1_001","KHT30_S29_R1_001")

#Write the required data matrix to an excel file
write.table(G151_cpm_norm, "G151_TMM_CPM.txt", sep="\t") #Normalized counts
write.table(G151_counts_matrix, "G151_raw_counts.txt", sep="\t") #Raw counts

#If we want the normalized pseudo-counts, useful for instance for cluster analysis, we can get them with the following commands.
scale <- G151_normalized$samples$lib.size*G151_normalized$samples$norm.factors
G151_normCounts <- round(t(t(G151_counts_matrix)/scale)*mean(scale))

#Boxplot
log2(G151_cpm_norm+1)