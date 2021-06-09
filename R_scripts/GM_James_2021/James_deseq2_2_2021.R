source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# Load DESeq2 library
library("DESeq2")

#Please check before installing if you have installed devtools, if not then only do the following otherwise just load the library
#install.packages("devtools")
library(devtools)
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
library(pheatmap)
library(cowplot)


# Set directory
directory <- "/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021"
setwd(directory)
# ALL Matrix together
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021/counts_matrix_reverse.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021/sample_annotation.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$GeneID
count_matrix <- subset(count_matrix, select = -c(GeneID))
sampleCondition <- c("D2mdx-6wk","B10WT-6mo","D2mdx-6wk","B10WT-6mo",
                     "B10WT-6mo","D2WT-6mo","D2mdx-6wk","D2WT-6mo",
                     "B10mdx-6wk","B10mdx-6wk","B10mdx-6wk","B10mdx-6mo",
                     "B10mdx-6mo","D2mdx-6wk","B10WT-6mo","B10mdx-6mo",
                     "B10mdx-6mo","D2mdx-6wk","D2mdx-6mo","D2mdx-6mo",
                     "D2mdx-6mo","D2WT-6mo","B10mdx-6mo")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)

dds_vst <- vst(dds)
plotPCA(dds_vst)
pcaplot_vst<- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot_vst
ggsave(pcaplot_vst,filename = "PCA_reverse_plot_vst_with_labels.png")

dds2 <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds2, normalized=TRUE) + 1),
                           colData=colData(dds2))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
pcaplot2 <- plotPCA( DESeqTransform( se ) ) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot2
ggsave(pcaplot2, filename = "PCA_plot_reverse_with_normalization.png")


se2 <- SummarizedExperiment(log2(counts(dds2, normalized=FALSE) + 1),
                            colData=colData(dds2))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
pcaplot3 <- plotPCA( DESeqTransform( se2 ) ) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot3
ggsave(pcaplot3, filename = "PCA_plot_no_normalization.png")


# Differential Gene Expression
dds_deg <- DESeq(dds)
sizeFactors(dds_deg)

# Normalized Counts
normalized_counts <- counts(dds_deg, normalized=TRUE)
head(normalized_counts)
View(counts(dds_deg))
write.table(normalized_counts, file="/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021/normalized_reverse_counts.txt", sep="\t", quote=F, col.names=NA)

# Extrcat results
res <- results(dds_deg)
res
summary(res)

# Heirarachial clustering 
# https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/03_DGE_QC_analysis.md
dds_rlog <- rlog(dds, blind =TRUE)
### Extract the rlog matrix from the object
rld_mat <- assay(dds_rlog)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
# also possible to perform custom transformation:
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
### Load pheatmap package
library(pheatmap)
### Plot heatmap
heirarchial_heatmap <- pheatmap(rld_cor, cluster_rows = T, cluster_cols = T)
ggsave2(heirarchial_heatmap, filename = "heirarchial_heatmap.png")

## Blues color heatmap
heat.colors <- brewer.pal(6, "Blues")
heirarchial_heatmap_blue <- pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
                                     fontsize_row = 10, height=20)
ggsave2(heirarchial_heatmap_blue, filename = "heirarchial_heatmap_blue.png")


# PCA
# The calculation is done by a singular value decomposition of the (centered and possibly scaled) 
# data matrix, not by using eigen on the covariance matrix.
# Performs a principal components analysis on the given data matrix
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(sample_annotation, pca$x)
df <- data.frame(pca$x)
head(df)
pca_plot_prcomp <- ggplot(df) + geom_point(aes(x=PC4, y=PC5, color = condition, size =1)) + geom_text(aes(x=PC4, y=PC5, label=sample_annotation$sample), vjust=0.2)
pca_plot_prcomp
ggsave2(pca_plot_prcomp, filename = "pca_plot_PC4-5.png")

# PCA Plot to see sample distribution
# same results with and without DESeq
dds_vst <- vst(dds)
mat_vst_dds <- assay(dds_vst)
plotPCA(mat_vst_dds)

# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))
# mat_all = assay(dds_vst)
# mat_all = mat_all - rowMeans(mat_all)
# write.table(vst,"pig_vst.txt",sep = "\t")
pcaplot<- plotPCA(dds_rlog, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot
ggsave(pcaplot,filename = "PCA_plot_with_labels.png")

pcaplot_vst<- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot_vst
ggsave(pcaplot_vst,filename = "PCA_reverse_plot_vst_with_labels.png")

dds_deg_vst <- vst(dds_deg)
pcaplot_deg_vst<- plotPCA(dds_deg_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot_vst
ggsave(pcaplot_vst,filename = "PCA_reverse_plot_vst_with_labels.png")

##########################################
##########################################
# B10_mdx_1.5 vs D2_mdx_1.5
##########################################
##########################################
outputPrefix1 <- "B10_mdx_1.5_vs_D2_mdx_1.5"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles1 <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS1_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS3_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136250_TS8_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS20_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames1 <- c("TS10","TS11","TS12","TS1","TS3","TS8","TS20")
sampleCondition1 <- c("B10_mdx_1.5","B10_mdx_1.5","B10_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5")
sampleTable1 <- data.frame(sampleName = sampleNames1, fileName = sampleFiles1, condition = sampleCondition1)
treatments1 = c("B10_mdx_1.5","D2_mdx_1.5")
condition = sampleCondition1
design = ~ condition
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable1, directory = directory1, design = ~ condition)
colData(ddsHTSeq1)$condition <- factor(colData(ddsHTSeq1)$condition,levels = treatments1)

dds1 <- DESeq(ddsHTSeq1)
keep1 <- rowSums(counts(dds1)) >= 10
dds1 <- dds1[keep1,]
res1 <- results(dds1)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res1
summary(res1)
# out of 19839 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2876, 14%
# LFC < 0 (down)     : 2499, 13%
# outliers [1]       : 51, 0.26%
# low counts [2]     : 1154, 5.8%
# (mean count < 2)

# How many adjusted p-values were less than 0.1?
sum(res1$padj < 0.1, na.rm=TRUE)
# [1] 5375
sum(res1$padj < 0.05, na.rm=TRUE)
# [1] 3924
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata1)[1] <- 'ensembl_gene_id'
write.csv(resdata1, file = paste0(outputPrefix1, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_mdx_1.5_vs_D2_mdx_1.5-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_mdx_1.5_vs_D2_mdx_1.5-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix1, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot1 <- plotMA(dds1, ylim=c(-5,5), main = "B10-mdx 6 weeks vs D2-mdx 6 weeks")
ggsave(ma_plot1, file=paste0(outputPrefix1, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld1 <- rlogTransformation(dds1, blind=T)
vsd1 <- varianceStabilizingTransformation(dds1, blind=T)

# PCA Plot
condition <- treatments1
pcaplot1 <- plotPCA(vsd1, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot1
ggsave(pcaplot1,file=paste0(outputPrefix1, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat1 = assay(vsd1)[ head(order(res1$padj),100), ] # select the top 100 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df1 = as.data.frame(colData(vsd1)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df1) = "condition" # Rename the column header
rownames(df1) = colnames(mat1) # add rownames
# and plot the actual heatmap
pheatmap1 <- pheatmap(mat1, annotation_col=df1, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap1
ggsave(pheatmap1,file=paste0(outputPrefix1, "-Heatmap.png"))

##########################################
##########################################
# B10_mdx_1.5 vs D2_mdx_1.5 v2 without TS11 
##########################################
##########################################
outputPrefix5 <- "B10_mdx_1.5_vs_D2_mdx_1.5_v2"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles5 <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS1_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS3_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136250_TS8_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS20_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames5 <- c("TS10","TS12","TS1","TS3","TS8","TS20")
sampleCondition5 <- c("B10_mdx_1.5","B10_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5","D2_mdx_1.5")
sampleTable5 <- data.frame(sampleName = sampleNames5, fileName = sampleFiles5, condition = sampleCondition5)
treatments5 = c("B10_mdx_1.5","D2_mdx_1.5")
condition = sampleCondition5
design = ~ condition
ddsHTSeq5 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable5, directory = directory1, design = ~ condition)
colData(ddsHTSeq5)$condition <- factor(colData(ddsHTSeq5)$condition,levels = treatments5)
colData(ddsHTSeq5)
dds5 <- DESeq(ddsHTSeq5)
keep5 <- rowSums(counts(dds5)) >= 10
dds5 <- dds5[keep5,]
res5 <- results(dds5)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res5
summary(res5)
# out of 19439 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2013, 10%
# LFC < 0 (down)     : 1634, 8.4%
# outliers [1]       : 4, 0.021%
# low counts [2]     : 1508, 7.8%
# (mean count < 3)


# How many adjusted p-values were less than 0.1?
sum(res5$padj < 0.1, na.rm=TRUE)
# [1] 5375
sum(res5$padj < 0.05, na.rm=TRUE)
# [1] 3924
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata5)[1] <- 'ensembl_gene_id'
write.csv(resdata5, file = paste0(outputPrefix5, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_mdx_1.5_vs_D2_mdx_1.5_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_mdx_1.5_vs_D2_mdx_1.5_v2-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix5, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot5 <- plotMA(dds5, ylim=c(-5,5), main = "B10-mdx 6 weeks vs D2-mdx 6 weeks")
ggsave(ma_plot5, file=paste0(outputPrefix5, "-MAplot.png"))
save_plot(ma_plot5, file=paste0(outputPrefix5, "-MAplot.png"))
# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld5 <- rlogTransformation(dds5, blind=T)
vsd5 <- varianceStabilizingTransformation(dds5, blind=T)

# PCA Plot
condition <- treatments5
pcaplot5 <- plotPCA(vsd5, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot5
ggsave(pcaplot5,file=paste0(outputPrefix5, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat5 = assay(vsd5)[ head(order(res5$padj),100), ] # select the top 100 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df5 = as.data.frame(colData(vsd5)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df5) = "condition" # Rename the column header
rownames(df5) = colnames(mat5) # add rownames
# and plot the actual heatmap
pheatmap5 <- pheatmap(mat5, annotation_col=df5, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap5
ggsave(pheatmap5,file=paste0(outputPrefix5, "-Heatmap.png"))


##########################################
##########################################
# D2_wt_6 vs D2_mdx_6
##########################################
##########################################
outputPrefix2 <- "D2_wt_6_vs_D2_mdx_6"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles2 <- c("deseq2_counts_star_out_merged_TS6_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS9_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136246_TS25_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS21_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames2 <- c("TS6","TS9","TS25","TS21","TS22","TS23")
sampleCondition2 <- c("D2_wt_6","D2_wt_6","D2_wt_6","D2_mdx_6","D2_mdx_6","D2_mdx_6")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition = sampleCondition2)
treatments2 = c("D2_wt_6","D2_mdx_6")
condition = sampleCondition2
design = ~ condition
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory1, design = ~ condition)
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
colData(ddsHTSeq2)
dds2 <- DESeq(ddsHTSeq2)
keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep2,]
res2 <- results(dds2)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res2
summary(res2)
# out of 18966 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 709, 3.7%
# LFC < 0 (down)     : 1131, 6%
# outliers [1]       : 7, 0.037%
# low counts [2]     : 5516, 29%

# How many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
# [1] 1840
sum(res2$padj < 0.05, na.rm=TRUE)
# [1] 706
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, file = paste0(outputPrefix2, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "D2_wt_6_vs_D2_mdx_6-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "D2_wt_6_vs_D2_mdx_6-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix2, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot2 <- plotMA(dds2, ylim=c(-5,5), main = "D2-wt 6months vs D2-mdx 6months")
ggsave(ma_plot2, file=paste0(outputPrefix2, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot2
ggsave(pcaplot2,file=paste0(outputPrefix2, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat2 = assay(vsd2)[ head(order(res2$padj),100), ] # select the top 100 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df2 = as.data.frame(colData(vsd2)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df2) = "condition" # Rename the column header
rownames(df2) = colnames(mat2) # add rownames
# and plot the actual heatmap
pheatmap2 <- pheatmap(mat2, annotation_col=df2, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap2
ggsave(pheatmap2,file=paste0(outputPrefix2, "-Heatmap.png"))

##########################################
##########################################
# D2_wt_6 vs D2_mdx_6 v2 without TS21
##########################################
##########################################
outputPrefix7 <- "D2_wt_6_vs_D2_mdx_6_v2"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles7 <- c("deseq2_counts_star_out_merged_TS6_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS9_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136246_TS25_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames7 <- c("TS6","TS9","TS25","TS22","TS23")
sampleCondition7 <- c("D2_wt_6","D2_wt_6","D2_wt_6","D2_mdx_6","D2_mdx_6")
sampleTable7 <- data.frame(sampleName = sampleNames7, fileName = sampleFiles7, condition = sampleCondition7)
treatments7 = c("D2_wt_6","D2_mdx_6")
condition = sampleCondition7
design = ~ condition
ddsHTSeq7 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable7, directory = directory1, design = ~ condition)
colData(ddsHTSeq7)$condition <- factor(colData(ddsHTSeq7)$condition,levels = treatments7)
colData(ddsHTSeq7)
dds7 <- DESeq(ddsHTSeq7)
keep7 <- rowSums(counts(dds7)) >= 10
dds7 <- dds7[keep7,]
res7 <- results(dds7)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res7
summary(res7)
# out of 18669 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 49, 0.26%
# LFC < 0 (down)     : 88, 0.47%
# outliers [1]       : 0, 0%
# low counts [2]     : 4706, 25%

# How many adjusted p-values were less than 0.1?
sum(res7$padj < 0.1, na.rm=TRUE)
# [1] 137
sum(res7$padj < 0.05, na.rm=TRUE)
# [1] 58
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata7 <- merge(as.data.frame(res7), as.data.frame(counts(dds7,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata7)[1] <- 'ensembl_gene_id'
write.csv(resdata7, file = paste0(outputPrefix7, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "D2_wt_6_vs_D2_mdx_6_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "D2_wt_6_vs_D2_mdx_6_v2-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix7, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot7 <- plotMA(dds7, ylim=c(-5,5), main = "D2-wt 6months vs D2-mdx 6months")
ggsave(ma_plot7, file=paste0(outputPrefix7, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld7 <- rlogTransformation(dds7, blind=T)
vsd7 <- varianceStabilizingTransformation(dds7, blind=T)

# PCA Plot
condition <- treatments7
pcaplot7 <- plotPCA(vsd7, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot7
ggsave(pcaplot7,file=paste0(outputPrefix7, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat7 = assay(vsd7)[ head(order(res7$padj),50), ] # select the top 50 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df7 = as.data.frame(colData(vsd7)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df7) = "condition" # Rename the column header
rownames(df7) = colnames(mat7) # add rownames
# and plot the actual heatmap
pheatmap7 <- pheatmap(mat7, annotation_col=df7, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap7
ggsave(pheatmap7,file=paste0(outputPrefix7, "-Heatmap.png"))


##########################################
##########################################
# B10_mdx_6 vs D2_mdx_6
##########################################
##########################################
outputPrefix3 <- "B10_mdx_6_vs_D2_mdx_6"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles3 <- c("deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136247_TS27_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS21_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames3 <- c("TS13","TS14","TS18","TS27","TS21","TS22","TS23")
sampleCondition3 <- c("B10_6","B10_6","B10_6","B10_6","D2_6","D2_6","D2_6")
sampleTable3 <- data.frame(sampleName = sampleNames3, fileName = sampleFiles3, condition = sampleCondition3)
treatments3 = c("B10_6","D2_6")
condition = sampleCondition3
design = ~ condition
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable3, directory = directory1, design = ~ condition)
colData(ddsHTSeq3)$condition <- factor(colData(ddsHTSeq3)$condition,levels = treatments3)
#guts
dds3 <- DESeq(ddsHTSeq3)
keep3 <- rowSums(counts(dds3)) >= 10
dds3 <- dds3[keep3,]
res3 <- results(dds3)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res3
summary(res3)
# out of 19312 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5046, 26%
# LFC < 0 (down)     : 5213, 27%
# outliers [1]       : 9, 0.047%
# low counts [2]     : 0, 0%

# How many adjusted p-values were less than 0.1?
sum(res3$padj < 0.1, na.rm=TRUE)
# [1] 10259
sum(res3$padj < 0.05, na.rm=TRUE)
# [1] 9038
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata3)[1] <- 'ensembl_gene_id'
write.csv(resdata3, file = paste0(outputPrefix3, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_mdx_6_vs_D2_mdx_6-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_mdx_6_vs_D2_mdx_6-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix3, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_6 vs D2_6
ma_plot3 <- plotMA(dds3, ylim=c(-5,5), main = "B10-mdx 6months vs D2-mdx 6months")
ggsave(ma_plot3, file=paste0(outputPrefix3, "-MAplot.png"))


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld3 <- rlogTransformation(dds3, blind=T)
vsd3 <- varianceStabilizingTransformation(dds3, blind=T)

# PCA Plot
condition <- treatments3
pcaplot3 <- plotPCA(vsd3, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot3
ggsave(pcaplot3,file=paste0(outputPrefix3, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat3 = assay(vsd3)[ head(order(res3$padj),100), ] # select the top 50 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df3 = as.data.frame(colData(vsd3)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df3) = "condition" # Rename the column header
rownames(df3) = colnames(mat3) # add rownames
# and plot the actual heatmap
pheatmap3 <- pheatmap(mat3, annotation_col=df3, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap3
ggsave(pheatmap3,file=paste0(outputPrefix3, "-Heatmap.png"))


##########################################
##########################################
# B10_mdx_6 vs D2_mdx_6 without TS21
##########################################
##########################################
outputPrefix6 <- "B10_mdx_6_vs_D2_mdx_6_v2"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles6 <- c("deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136247_TS27_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames6 <- c("TS13","TS14","TS18","TS27","TS22","TS23")
sampleCondition6 <- c("B10_6","B10_6","B10_6","B10_6","D2_6","D2_6")
sampleTable6 <- data.frame(sampleName = sampleNames6, fileName = sampleFiles6, condition = sampleCondition6)
treatments6 = c("B10_6","D2_6")
condition = sampleCondition6
design = ~ condition
ddsHTSeq6 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable6, directory = directory1, design = ~ condition)
colData(ddsHTSeq6)$condition <- factor(colData(ddsHTSeq6)$condition,levels = treatments6)
colData(ddsHTSeq6)
dds6 <- DESeq(ddsHTSeq6)
keep6 <- rowSums(counts(dds6)) >= 10
dds6 <- dds6[keep6,]
res6 <- results(dds6)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res6
summary(res6)
# out of 19024 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4461, 23%
# LFC < 0 (down)     : 4789, 25%
# outliers [1]       : 8, 0.042%
# low counts [2]     : 0, 0%

# How many adjusted p-values were less than 0.1?
sum(res6$padj < 0.1, na.rm=TRUE)
# [1] 9250
sum(res6$padj < 0.05, na.rm=TRUE)
# [1] 8145
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata6 <- merge(as.data.frame(res6), as.data.frame(counts(dds6,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata6)[1] <- 'ensembl_gene_id'
write.csv(resdata6, file = paste0(outputPrefix6, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_mdx_6_vs_D2_mdx_6_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_mdx_6_vs_D2_mdx_6_v2-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix6, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot6 <- plotMA(dds6, ylim=c(-5,5), main = "B10-mdx 6months vs D2-mdx 6months")
ggsave(ma_plot6, file=paste0(outputPrefix6, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld6 <- rlogTransformation(dds6, blind=T)
vsd6 <- varianceStabilizingTransformation(dds6, blind=T)

# PCA Plot
condition <- treatments6
pcaplot6 <- plotPCA(vsd6, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot6
ggsave(pcaplot6,file=paste0(outputPrefix6, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat6 = assay(vsd6)[ head(order(res6$padj),100), ] # select the top 50 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df6 = as.data.frame(colData(vsd6)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df6) = "condition" # Rename the column header
rownames(df6) = colnames(mat6) # add rownames
# and plot the actual heatmap
pheatmap6 <- pheatmap(mat6, annotation_col=df6, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap6
ggsave(pheatmap6,file=paste0(outputPrefix6, "-Heatmap.png"))


##########################################
##########################################
# B10_wt_6 vs B10_mdx_6
##########################################
##########################################
outputPrefix4 <- "B10_wt_6_vs_B10_mdx_6"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles4 <- c("deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136241_TS2_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS4_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS16_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136247_TS27_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames4 <- c("TS2","TS4","TS16","TS13","TS14","TS18","TS27")
sampleCondition4 <- c("B10_wt_6","B10_wt_6","B10_wt_6","B10_mdx_6","B10_mdx_6","B10_mdx_6","B10_mdx_6")
sampleTable4 <- data.frame(sampleName = sampleNames4, fileName = sampleFiles4, condition = sampleCondition4)
treatments4 = c("B10_wt_6","B10_mdx_6")
condition = sampleCondition4
design = ~ condition
ddsHTSeq4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable4, directory = directory1, design = ~ condition)
colData(ddsHTSeq4)$condition <- factor(colData(ddsHTSeq4)$condition,levels = treatments4)
#guts
dds4 <- DESeq(ddsHTSeq4)
keep4 <- rowSums(counts(dds4)) >= 10
dds4 <- dds4[keep4,]
res4 <- results(dds4)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res4
summary(res4)
# out of 19112 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 583, 3.1%
# LFC < 0 (down)     : 247, 1.3%
# outliers [1]       : 35, 0.18%
# low counts [2]     : 3335, 17%

# How many adjusted p-values were less than 0.1?
sum(res4$padj < 0.1, na.rm=TRUE)
# [1] 830
sum(res4$padj < 0.05, na.rm=TRUE)
# [1] 538
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered3 <- res3[order(res3$pvalue),]
# save data results and normalized reads to csv
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata4)[1] <- 'ensembl_gene_id'
write.csv(resdata4, file = paste0(outputPrefix4, "-results-with-normalized.csv"))
# # Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_wt_6_vs_B10_mdx_6-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_wt_6_vs_B10_mdx_6-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = paste0(outputPrefix4, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

ma_plot4 <- plotMA(dds4, ylim=c(-5,5), main = "B10-wt 6months vs B10-mdx 6months")
ggsave(ma_plot4, file=paste0(outputPrefix4, "-MAplot.png"))


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld4 <- rlogTransformation(dds4, blind=T)
vsd4 <- varianceStabilizingTransformation(dds4, blind=T)

# PCA Plot
condition <- treatments4
pcaplot4 <- plotPCA(vsd4, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
pcaplot4
ggsave(pcaplot4,file=paste0(outputPrefix4, "-PCA_plot.png"))

# Heatmap
library("pheatmap")
mat4 = assay(vsd4)[ head(order(res4$padj),100), ] # select the top 100 genes with the lowest padj
# mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df4 = as.data.frame(colData(vsd4)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df4) = "condition" # Rename the column header
rownames(df4) = colnames(mat4) # add rownames
# and plot the actual heatmap
pheatmap4 <- pheatmap(mat4, annotation_col=df4, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
pheatmap4
ggsave(pheatmap4,file=paste0(outputPrefix4, "-Heatmap.png"))


