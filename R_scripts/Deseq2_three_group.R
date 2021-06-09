save.image("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/mouse_three_group_comparison.RData")
load("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/mouse_three_group_comparison.RData")

# If Deseq2 is not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
# Load DESeq2 library, if Deseq2 is already installed
library("DESeq2")

#Please check before installing if you have installed devtools, if not then only do the following otherwise just load the library
#install.packages("devtools")
library(devtools)
library(ggplot2)
library(gplots)
library(knitr)
library("biomaRt")
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
#listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")

setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/")

# ALL Matrix together
count_matrix_mouse <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Counts_matrix.csv', header =T, sep = ",")
rownames(count_matrix_mouse) <- count_matrix_mouse$GeneID
count_matrix_mouse <- subset(count_matrix_mouse, select = -c(GeneID))


sample_annotation_mouse <- read.table(file='/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/sample_annotation.csv', header =T, sep = ",")
Samplecondition_mouse <- c( "control","control","control","hpx_lps","hpx_lps","hpx_lps","hpx","hpx","hpx","lps","lps","lps")
condition = Samplecondition_mouse
sampleNames_mouse <- c("C1","C2","C3","DH1","DH2","DH3","H1M","H2M","H3M","L1","L2","L3")
sampleTable$condition <- factor(sampleTable$condition)
dds_mouse <- DESeqDataSetFromMatrix(countData = count_matrix_mouse,
                                  colData = sample_annotation_mouse, 
                                  design = ~ condition)
# This step is somehow essential because the order of the conditions remain as you direct 
# and is not sorted alphabetially
dds_mouse$condition <- factor( c( "control","control","control","hpx_lps","hpx_lps","hpx_lps","hpx","hpx","hpx","lps","lps","lps") , levels = c("control","hpx_lps","hpx","lps"))
# Pre-filtering, keeping only rows that have at least 10 reads total
keep <- rowSums(counts(dds_mouse)) >= 10
dds_mouse <- dds_mouse[keep,]
dds <- DESeq(dds_mouse)
resultsNames(dds)

normalized_dds <- as.data.frame(counts(dds,normalized =TRUE))
normalized_dds = as.matrix(normalized_dds)
normalized_dds = as.data.frame(normalized_dds)

# Control vs hpx_lps
res = results(dds, contrast=c("condition","control","hpx_lps"))
summary(res)
res <- res[order(res$padj),]
head(res)
kable(res[1:5,-(3:4)])

# Contrast 1
# Control vs. any hit (Control - Any hit)
# # control hpx_lps hpx lps
# control vs all other conditions with control as reference
res1 = results(dds, contrast=c(-1,1/3,1/3,1/3))
summary(res1)
res1 <- res1[order(res1$padj),]
head(res1)
res1_df <- as.data.frame(res1)
head(res1_df)
kable(res1[1:5,c(2,5:6)])
res1_df <- as.data.frame(res1)
res1_df$ensembl_gene_id <- rownames(res1_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res1_df$ensembl_gene_id, mart= ensembl)
m1 <- merge(res1_df, b, by="ensembl_gene_id")
filter_m1 <- m1$padj < 0.05
m1_subset_padj_0.05 <- m1[filter_m1,]
m1_subset_padj_0.05 <- as.matrix(m1_subset_padj_0.05)
m1_subset_padj_0.05_df <- as.data.frame(m1_subset_padj_0.05)
row.names(m1_subset_padj_0.05) <- m1_subset_padj_0.05_df$ensembl_gene_id
contrast1 <- merge(m1_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast1 <- subset(contrast1, select = -c(Row.names))
write.table(contrast1, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast1.txt", sep = "\t")



# Contrast 2
# Impact of hitting H (main effect, no hit - hit)
# control hpx_lps hpx lps
# The effect of hypoxia in hpx and hpx_lps over control and lps as reference
res2 = results(dds, contrast=c(-1/2,1/2,1/2,-1/2))
summary(res2)
res2 <- res2[order(res2$padj),]
head(res2)
kable(res2[1:5,c(2,5:6)])
# names(resdata6)[1] <- 'ensembl_gene_id'
res2_df <- as.data.frame(res2)
res2_df$ensembl_gene_id <- rownames(res2_df)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description" ), values=res2_df$ensembl_gene_id, mart= ensembl)
m <- merge(res2_df, b, by="ensembl_gene_id")
filter_m <- m$padj < 0.05
m_subset_padj_0.05 <- m[filter_m,]
m_subset_padj_0.05 <- as.matrix(m_subset_padj_0.05)
m_subset_padj_0.05_df <- as.data.frame(m_subset_padj_0.05)
row.names(m_subset_padj_0.05) <- m_subset_padj_0.05_df$ensembl_gene_id
contrast2 <- merge(m_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast2 <- subset(contrast2, select = -c(Row.names))
write.table(contrast2, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/contrast2.txt", sep = "\t")


# Contrast 3 
# Impact of hitting L (main effect, no hit - hit)
# control hpx_lps hpx lps
# The effect of lps in lps and hpx_lps over control and hpx as reference
res3 = results(dds, contrast=c(-1/2,1/2,-1/2,1/2))
summary(res3)
res3 <- res3[order(res3$padj),]
head(res3)
kable(res3[1:5,c(2,5:6)])
res3_df <- as.data.frame(res3)
res3_df$ensembl_gene_id <- rownames(res3_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res3_df$ensembl_gene_id, mart= ensembl)
m3 <- merge(res3_df, b, by="ensembl_gene_id")
filter_m3 <- m3$padj < 0.05
m3_subset_padj_0.05 <- m3[filter_m3,]
m3_subset_padj_0.05 <- as.matrix(m3_subset_padj_0.05)
m3_subset_padj_0.05_df <- as.data.frame(m3_subset_padj_0.05)
row.names(m3_subset_padj_0.05) <- m3_subset_padj_0.05_df$ensembl_gene_id
contrast3 <- merge(m3_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast3 <- subset(contrast3, select = -c(Row.names))
write.table(contrast3, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast3.txt", sep = "\t")


# Contrast 4
# Impact of hitting L vs hitting H (difference between main effects L and main effect of H)
# L vs H 
# control hpx_lps hpx lps
# Effect of LPS over hypoxia
res4 = results(dds, contrast=c(0,0,-1,1))
summary(res4)
res4 <- res4[order(res4$padj),]
head(res4)
kable(res4[1:5,c(2,5:6)])
res4_df <- as.data.frame(res4)
res4_df$ensembl_gene_id <- rownames(res4_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res4_df$ensembl_gene_id, mart= ensembl)
m4 <- merge(res4_df, b, by="ensembl_gene_id")
filter_m4 <- m4$padj < 0.05
m4_subset_padj_0.05 <- m4[filter_m4,]
m4_subset_padj_0.05 <- as.matrix(m4_subset_padj_0.05)
m4_subset_padj_0.05_df <- as.data.frame(m4_subset_padj_0.05)
row.names(m4_subset_padj_0.05) <- m4_subset_padj_0.05_df$ensembl_gene_id
contrast4 <- merge(m4_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast4 <- subset(contrast4, select = -c(Row.names,C1,C2,C3,DH1,DH2,DH3))
write.table(contrast4, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast4.txt", sep = "\t")

# Contrast 5
# Impact of hitting L when H was not hit
# control hpx_lps hpx lps
# Effect of LPS over control
res5 = results(dds, contrast=c(1,0,0,-1))
summary(res5)
res5 <- res5[order(res5$padj),]
head(res5)
kable(res5[1:5,c(2,5:6)])
res5_df <- as.data.frame(res5)
res5_df$ensembl_gene_id <- rownames(res5_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res5_df$ensembl_gene_id, mart= ensembl)
m5 <- merge(res5_df, b, by="ensembl_gene_id")
filter_m5 <- m5$padj < 0.05
m5_subset_padj_0.05 <- m5[filter_m5,]
m5_subset_padj_0.05 <- as.matrix(m5_subset_padj_0.05)
m5_subset_padj_0.05_df <- as.data.frame(m5_subset_padj_0.05)
row.names(m5_subset_padj_0.05) <- m5_subset_padj_0.05_df$ensembl_gene_id
contrast5 <- merge(m5_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast5 <- subset(contrast5, select = -c(Row.names,DH1,DH2,DH3,H1M,H2M,H3M))
write.table(contrast5, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast5.txt", sep = "\t")

# Contrast 6
# Impact of hitting L when H was hit
# control hpx_lps hpx lps
# Impact of L in hpx_lps over hpx
res6 = results(dds, contrast=c(0,1,-1,0))
summary(res6)
res6 <- res6[order(res6$padj),]
#head(res6)
kable(res6[1:5,c(2,5:6)])
res6_df <- as.data.frame(res6)
res6_df$ensembl_gene_id <- rownames(res6_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res6_df$ensembl_gene_id, mart= ensembl)
m6 <- merge(res6_df, b, by="ensembl_gene_id")
filter_m6 <- m6$padj < 0.05
m6_subset_padj_0.05 <- m6[filter_m6,]
m6_subset_padj_0.05 <- as.matrix(m6_subset_padj_0.05)
m6_subset_padj_0.05_df <- as.data.frame(m6_subset_padj_0.05)
row.names(m6_subset_padj_0.05) <- m6_subset_padj_0.05_df$ensembl_gene_id
contrast6 <- merge(m6_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast6 <- subset(contrast6, select = -c(Row.names,C1,C2,C3,L1,L2,L3))
write.table(contrast6, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast6.txt", sep = "\t")

# Contrast 7
# Impact of double hit compared to any single hit
# control hpx_lps hpx lps
# Effect of double hit hpx_lps over hypox and lps as reference
res7 = results(dds, contrast=c(0,1,-1/2,-1/2))
summary(res7)
res7 <- res7[order(res7$padj),]
res7_df <- as.data.frame(res7)
head(res7)
kable(res7[1:5,c(2,5:6)])
res7_df <- as.data.frame(res7)
res7_df$ensembl_gene_id <- rownames(res7_df)
b <- getBM(filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id","external_gene_name","description"), 
           values=res7_df$ensembl_gene_id, mart= ensembl)
m7 <- merge(res7_df, b, by="ensembl_gene_id")
filter_m7 <- m7$padj < 0.05
m7_subset_padj_0.05 <- m7[filter_m7,]
m7_subset_padj_0.05 <- as.matrix(m7_subset_padj_0.05)
m7_subset_padj_0.05_df <- as.data.frame(m7_subset_padj_0.05)
row.names(m7_subset_padj_0.05) <- m7_subset_padj_0.05_df$ensembl_gene_id
contrast7 <- merge(m7_subset_padj_0.05, normalized_dds, by = 'row.names', sort = FALSE)
contrast7 <- subset(contrast7, select = -c(Row.names,C1,C2,C3))
write.table(contrast7, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Three_way_contrasts/contrast7.txt", sep = "\t")




