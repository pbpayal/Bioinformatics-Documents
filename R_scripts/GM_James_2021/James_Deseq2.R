save.image(file = "/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021/James_2021.RData")
load("/Users/pbanerjee/Documents/Payal_Scripts/R/GM_James_2021/James_2021.RData.RData")

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

# Nirmalized Counts
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

# Set the prefix for each output file name
##########################################
##########################################
# B10_1.5 vs D2_1.5
##########################################
##########################################
# v1
outputPrefix <- "B10_mdx_1.5_vs_D2_mdx_1.5_v1"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS3_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136250_TS8_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS20_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("TS10","TS11","TS12","TS3","TS8","TS20")
sampleCondition <- c("B10_1.5","B10_1.5","B10_1.5","D2_1.5","D2_1.5","D2_1.5")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("B10_1.5","D2_1.5")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory1, design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res
# # order results by padj value (most significant to least)
# res= subset(res, padj<0.05)
# res <- res[order(res$padj),]
# res
summary(res)
# out of 26871 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2740, 10%
# LFC < 0 (down)     : 2770, 10%
# outliers [1]       : 22, 0.082%
# low counts [2]     : 8395, 31%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
mcols(res, use.names=TRUE)
# DataFrame with 6 rows and 2 columns
# type                                         description
# <character>                                         <character>
#   baseMean       intermediate           mean of normalized counts for all samples
# log2FoldChange      results log2 fold change (MLE): condition D2 1.5 vs B10 1.5
# lfcSE               results         standard error: condition D2 1.5 vs B10 1.5
# stat                results         Wald statistic: condition D2 1.5 vs B10 1.5
# pvalue              results      Wald test p-value: condition D2 1.5 vs B10 1.5
# padj                results                                BH adjusted p-values


# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
# [1] 5510
sum(res$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'ensembl_gene_id'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig <- subset(resOrdered, padj < 0.1)
# resSig
# summary(resSig)
# write.csv(as.data.frame(resSig), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# resSig_05 <- subset(resOrdered, padj < 0.05)
# resSig_05
# summary(resSig_05)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.3.csv"))

# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_1.5_vs_D2_1.5_v1-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5_vs_D2_1.5_v1-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix, "normalized-Biomart-results.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_vs_D2_1.5
maplot <- plotMA(dds, ylim=c(-5,5),main = "B10_vs_D2_1.5")
ggsave(maplot,file=paste0(outputPrefix, "-MAplot.png"))



# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
#rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# PCA Plot
condition <- treatments
pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=0.7)
pcaplot
ggsave(pcaplot,file=paste0(outputPrefix, "-PCA.png"))

# # Heatmap
# library( "genefilter" )
# topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
# dev.copy(jpeg, paste0(outputPrefix, "-heatmap3.jpeg"))
# heatmap <- heatmap.2( assay(rld)[topVarGenes, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
# dev.off()


# v2 with TS10,TS11,TS12, TS1, TS3, TS8, TS20
outputPrefix <- "B10_1.5_vs_D2_1.5_v2"
directory1 <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse" 
sampleFiles <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_counts_star_out_merged_TS1_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS3_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136250_TS8_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS20_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("TS10","TS11","TS12","TS1","TS3","TS8","TS20")
sampleCondition <- c("B10_1.5","B10_1.5","B10_1.5","D2_1.5","D2_1.5","D2_1.5","D2_1.5")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("B10_1.5","D2_1.5")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory1, design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res

summary(res)

mcols(res, use.names=TRUE)


# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
# [1] 5510
sum(res$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# # We can order our results table by the smallest p value
# resOrdered <- res[order(res$pvalue),]

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'ensembl_gene_id'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# # followed by the write.csv function.
# resSig <- subset(resOrdered, padj < 0.1)
# resSig
# summary(resSig)
# write.csv(as.data.frame(resSig), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# resSig_05 <- subset(resOrdered, padj < 0.05)
# resSig_05
# summary(resSig_05)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.3.csv"))

# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_1.5_vs_D2_1.5_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5_vs_D2_1.5_v2-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
# listMarts()
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix, "normalized-Biomart-results.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_vs_D2_1.5
maplot <- plotMA(dds, ylim=c(-5,5),main = "B10_vs_D2_1.5")
ggsave(maplot,file=paste0(outputPrefix, "-MAplot.png"))



# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
#rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# PCA Plot
condition <- treatments
pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix, "-PCA.png"))

# Heatmap
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix, "-heatmap3.jpeg"))
heatmap <- heatmap.2( assay(rld)[topVarGenes, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()



##########################################
##########################################
# B10_1.5 vs B10_6
##########################################
##########################################

outputPrefix2 <- "B10_1.5 vs B10_6_v1"
sampleFiles2 <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames2 <- c("TS10","TS11","TS12","TS13","TS14","TS18")
sampleCondition2 <- c("B10_1.5","B10_1.5","B10_1.5","B10_6","B10_6","B10_6")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition = sampleCondition2)
treatments2 = c("B10_1.5","B10_6")
design = ~ condition
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory1, design = ~ condition)
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
#guts
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
# # order results by padj value (most significant to least)
# res2 = subset(res2, padj<0.05)
# res2 <- res2[order(res2$padj),]
res2
summary(res2)
# out of 26092 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1573, 6%
# LFC < 0 (down)     : 1810, 6.9%
# outliers [1]       : 26, 0.1%
# low counts [2]     : 9558, 37%
# How many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
# [1] 3383
sum(res2$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered2 <- res2[order(res2$pvalue),]
# save data results and normalized reads to csv
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, file = paste0(outputPrefix2, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig <- subset(resOrdered, padj < 0.1)
# resSig
# summary(resSig)
# write.csv(as.data.frame(resSig), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05 <- subset(resOrdered, padj < 0.05)
# resSig_05
# summary(resSig_05)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_1.5 vs B10_6_v1-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5 vs B10_6_v1-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix2, "normalized-Biomart-results.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

maplot <- plotMA(dds2, ylim=c(-5,5),main = "B10_1.5 vs B10_6")
ggsave(maplot,file=paste0(outputPrefix2, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
# rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name), vjust=0.7)
pcaplot2
ggsave(pcaplot2,file=paste0(outputPrefix2, "-PCAplot.png"))

# # Heatmap
# library( "genefilter" )
# topVarGenes2 <- head( order( rowVars( assay(rld2) ), decreasing=TRUE ), 35 )
# dev.copy(jpeg, paste0(outputPrefix2, "-heatmap3.jpeg"))
# heatmap <- heatmap.2( assay(rld2)[topVarGenes2, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
# dev.off()

outputPrefix2 <- "B10_1.5 vs B10_6_v2"
sampleFiles2 <- c("deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136247_TS27_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames2 <- c("TS10","TS11","TS12","TS13","TS14","TS18","TS27")
sampleCondition2 <- c("B10_1.5","B10_1.5","B10_1.5","B10_6","B10_6","B10_6","B10_6")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition = sampleCondition2)
treatments2 = c("B10_1.5","B10_6")
design = ~ condition
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory1, design = ~ condition)
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
#guts
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
# # order results by padj value (most significant to least)
# res2 = subset(res2, padj<0.05)
# res2 <- res2[order(res2$padj),]
res2
summary(res2)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2663, 10%
# LFC < 0 (down)     : 2834, 11%
# outliers [1]       : 51, 0.19%
# low counts [2]     : 8840, 33%
# How many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
sum(res2$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)
# # We can order our results table by the smallest p value
# resOrdered2 <- res2[order(res2$pvalue),]
# save data results and normalized reads to csv
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, file = paste0(outputPrefix2, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig <- subset(resOrdered, padj < 0.1)
# resSig
# summary(resSig)
# write.csv(as.data.frame(resSig), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05 <- subset(resOrdered, padj < 0.05)
# resSig_05
# summary(resSig_05)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_1.5 vs B10_6_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5 vs B10_6_v2-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix2, "normalized-Biomart-results.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

maplot <- plotMA(dds2, ylim=c(-5,5),main = "B10_1.5 vs B10_6")
# save_plot(maplot,file=paste0(outputPrefix2, "-MAplot.png"))
ggsave(maplot,file=paste0(outputPrefix2, "-MAplot.png"))

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
# rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name), vjust=0.7)
pcaplot2
ggsave(pcaplot2,file=paste0(outputPrefix2, "-PCAplot.png"))

# # Heatmap
# library( "genefilter" )
# topVarGenes2 <- head( order( rowVars( assay(rld2) ), decreasing=TRUE ), 35 )
# dev.copy(jpeg, paste0(outputPrefix2, "-heatmap3.jpeg"))
# heatmap <- heatmap.2( assay(rld2)[topVarGenes2, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
# dev.off()

##########################################
##########################################
# B10_6 vs D2_6
##########################################
##########################################
outputPrefix3 <- "B10_6 vs D2_6_v1"
sampleFiles3 <- c("deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS21_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames3 <- c("TS13","TS14","TS18","TS21","TS22","TS23")
sampleCondition3 <- c("B10_6","B10_6","B10_6","D2_6","D2_6","D2_6")
sampleTable3 <- data.frame(sampleName = sampleNames3, fileName = sampleFiles3, condition = sampleCondition3)
treatments3 = c("B10_6","D2_6")
condition = sampleCondition3
design = ~ condition
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable3, directory = directory1, design = ~ condition)
colData(ddsHTSeq3)$condition <- factor(colData(ddsHTSeq3)$condition,levels = treatments3)
#guts
dds3 <- DESeq(ddsHTSeq3)
res3 <- results(dds3)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res3
summary(res3)
# out of 26720 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4322, 16%
# LFC < 0 (down)     : 4518, 17%
# outliers [1]       : 5, 0.019%
# low counts [2]     : 7852, 29%

# How many adjusted p-values were less than 0.1?
sum(res3$padj < 0.1, na.rm=TRUE)
# [1] 8840
sum(res3$padj < 0.05, na.rm=TRUE)
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
results_csv <- "B10_6 vs D2_6_v1-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_6 vs D2_6_v1-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix3, "normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_6 vs D2_6
plotMA(dds3, ylim=c(-5,5),main = "B10_6 vs D2_6")
dev.copy(png, paste0(outputPrefix3, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
# rld3 <- rlogTransformation(dds3, blind=T)
vsd3 <- varianceStabilizingTransformation(dds3, blind=T)

# PCA Plot
condition <- treatments3
pcaplot3 <- plotPCA(vsd3, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
pcaplot3
ggsave(pcaplot3,file=paste0(outputPrefix3, "PCA_plot.pdf"))

# # Heatmap
# library( "genefilter" )
# topVarGenes3 <- head( order( rowVars( assay(rld3) ), decreasing=TRUE ), 35 )
# dev.copy(jpeg, paste0(outputPrefix3, "-heatmap.jpeg"))
# heatmap3 <- heatmap.2( assay(rld3)[topVarGenes3, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
# dev.off()

outputPrefix3 <- "B10_6 vs D2_6_v2"
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
sizeFactors(dds3)
res3 <- results(dds3)
# order results by padj value (most significant to least)
# res3 = subset(res3, padj<0.05)
# res3 <- res3[order(res3$padj),]
res3
summary(res3)
# out of 27286 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5068, 19%
# LFC < 0 (down)     : 5181, 19%
# outliers [1]       : 9, 0.033%
# low counts [2]     : 7032, 26%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# How many adjusted p-values were less than 0.1?
sum(res3$padj < 0.1, na.rm=TRUE)
# [1] 10249
sum(res3$padj < 0.05, na.rm=TRUE)
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
results_csv <- "B10_6 vs D2_6_v2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_6 vs D2_6_v2-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix3, "normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_6 vs D2_6
plotMA(dds3, ylim=c(-5,5),main = "B10_6 vs D2_6")
dev.copy(png, paste0(outputPrefix3, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
# rld3 <- rlogTransformation(dds3, blind=T)
vsd3 <- varianceStabilizingTransformation(dds3, blind=T)

# PCA Plot
condition <- treatments3
pcaplot3 <- plotPCA(vsd3, intgroup=c("condition")) + geom_text(aes(label=name), vjust=0.5)
pcaplot3
ggsave(pcaplot3,file=paste0(outputPrefix3, "PCA_plot.png"))

# # Heatmap
# library( "genefilter" )
# topVarGenes3 <- head( order( rowVars( assay(rld3) ), decreasing=TRUE ), 35 )
# dev.copy(jpeg, paste0(outputPrefix3, "-heatmap.jpeg"))
# heatmap3 <- heatmap.2( assay(rld3)[topVarGenes3, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
# dev.off()

##########################################
##########################################
# Wild Type B10_6 vs MDX B10_6
##########################################
##########################################


outputPrefix5 <- "Wd_B10_6_vs_Mdx_B10_6"
sampleFiles5 <- c("deseq2_counts_star_out_merged_TS5_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS4_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS16_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136242_TS13_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136243_TS14_L007_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136245_TS18_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames5 <- c("TS5","TS4","TS16","TS13","TS14","TS18")
sampleCondition5 <- c("Wd_B10_6","Wd_B10_6","Wd_B10_6","Mdx_B10_6","Mdx_B10_6","Mdx_B10_6")
sampleTable5 <- data.frame(sampleName = sampleNames5, fileName = sampleFiles5, condition = sampleCondition5)
treatments5 = c("Wd_B10_6","Mdx_B10_6")
condition = sampleCondition5
design = ~ condition
ddsHTSeq5 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable5, directory = directory, design = ~ condition)
colData(ddsHTSeq5)$condition <- factor(colData(ddsHTSeq5)$condition,levels = treatments5)
#guts
dds5 <- DESeq(ddsHTSeq5)
res5 <- results(dds5)
# order results by padj value (most significant to least)
res5 = subset(res5, padj<0.05)
res5 <- res5[order(res5$padj),]
res5
summary(res5)
# How many adjusted p-values were less than 0.1?
sum(res5$padj < 0.1, na.rm=TRUE)
sum(res5$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered5 <- res5[order(res5$pvalue),]

# save data results and normalized reads to csv
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata5)[1] <- 'ensembl_gene_id'
write.csv(resdata5, file = paste0(outputPrefix5, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# # Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "Wd_B10_6_vs_Mdx_B10_6-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Wd_B10_6_vs_Mdx_B10_6-results-with-normalized.txt"
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

# B10_6 vs D2_6
plotMA(dds5, ylim=c(-5,5),main = "Wd_B10_6 vs Mdx_B10_6")
dev.copy(png, paste0(outputPrefix5, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld5 <- rlogTransformation(dds5, blind=T)
vsd5 <- varianceStabilizingTransformation(dds5, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
# B10_vs_B10_6
condition <- treatments5
pcaplot <- plotPCA(vsd5, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix5, "-PCA_plot.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes5 <- head( order( rowVars( assay(rld5) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix5, "-heatmap.jpeg"))
heatmap <- heatmap.2( assay(rld5)[topVarGenes5, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()

##########################################
##########################################
# Wild Type D2_6 vs MDX D2_6
##########################################
##########################################


outputPrefix6 <- "Wd_D2_6_vs_Mdx_D2_6"
sampleFiles6 <- c("deseq2_counts_star_out_merged_TS6_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS9_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136246_TS25_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS21_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames6 <- c("TS6","TS9","TS25","TS21","TS22","TS23")
sampleCondition6 <- c("Wd_D2_6","Wd_D2_6","Wd_D2_6","Mdx_D2_6","Mdx_D2_6","Mdx_D2_6")
sampleTable6 <- data.frame(sampleName = sampleNames6, fileName = sampleFiles6, condition = sampleCondition6)
treatments6 = c("Wd_D2_6","Mdx_D2_6")
condition = sampleCondition6
design = ~ condition
ddsHTSeq6 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable6, directory = directory, design = ~ condition)
colData(ddsHTSeq6)$condition <- factor(colData(ddsHTSeq6)$condition,levels = treatments6)
#guts
dds6 <- DESeq(ddsHTSeq6)
res6 <- results(dds6)
# order results by padj value (most significant to least)
res6 = subset(res6, padj<0.05)
res6 <- res6[order(res6$padj),]
res6
summary(res6)
# How many adjusted p-values were less than 0.1?
sum(res6$padj < 0.1, na.rm=TRUE)
sum(res6$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered6 <- res6[order(res6$pvalue),]

# save data results and normalized reads to csv
resdata6 <- merge(as.data.frame(res6), as.data.frame(counts(dds6,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata6)[1] <- 'ensembl_gene_id'
write.csv(resdata6, file = paste0(outputPrefix6, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# # Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "Wd_D2_6_vs_Mdx_D2_6-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Wd_D2_6_vs_Mdx_D2_6-results-with-normalized.txt"
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
plotMA(dds6, ylim=c(-5,5),main = "Wd_D2_6 vs Mdx_D2_6")
dev.copy(png, paste0(outputPrefix6, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld6 <- rlogTransformation(dds6, blind=T)
vsd6 <- varianceStabilizingTransformation(dds6, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
# B10_vs_B10_6
condition <- treatments6
pcaplot <- plotPCA(vsd6, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix6, "-PCA_plot.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes6 <- head( order( rowVars( assay(rld6) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix6, "-heatmap.jpeg"))
heatmap <- heatmap.2( assay(rld6)[topVarGenes6, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()

##########################################
##########################################
# Wild Type B10_6 vs Mdx B10_1.5
##########################################
##########################################


outputPrefix7 <- "Wd_B10_6_vs_Mdx_B10_1.5"
sampleFiles7 <- c("deseq2_counts_star_out_merged_TS5_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS4_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS16_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS10_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS11_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS12_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames7 <- c("TS5","TS4","TS16","TS10","TS11","TS12")
sampleCondition7 <- c("Wd_B10_6","Wd_B10_6","Wd_B10_6","Mdx_B10_1.5","Mdx_B10_1.5","Mdx_B10_1.5")
sampleTable7 <- data.frame(sampleName = sampleNames7, fileName = sampleFiles7, condition = sampleCondition7)
treatments7 = c("Wd_B10_6","Mdx_B10_1.5")
condition = sampleCondition7
design = ~ condition
ddsHTSeq7 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable7, directory = directory, design = ~ condition)
colData(ddsHTSeq7)$condition <- factor(colData(ddsHTSeq7)$condition,levels = treatments7)
#guts
dds7 <- DESeq(ddsHTSeq7)
res7 <- results(dds7)
summary(res7)
# order results by padj value (most significant to least)
res7 = subset(res7, padj<0.05)
res7 <- res7[order(res7$padj),]
res7
summary(res7)
# How many adjusted p-values were less than 0.1?
sum(res7$padj < 0.1, na.rm=TRUE)
sum(res7$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered7 <- res7[order(res7$pvalue),]

# save data results and normalized reads to csv
resdata7 <- merge(as.data.frame(res7), as.data.frame(counts(dds7,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata7)[1] <- 'ensembl_gene_id'
write.csv(resdata7, file = paste0(outputPrefix7, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# # Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "Wd_B10_6_vs_Mdx_B10_1.5-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Wd_B10_6_vs_Mdx_B10_1.5-results-with-normalized.txt"
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
plotMA(dds7, ylim=c(-5,5),main = "Wd_B10_6_vs_Mdx_B10_1.5")
dev.copy(png, paste0(outputPrefix7, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld7 <- rlogTransformation(dds7, blind=T)
vsd7 <- varianceStabilizingTransformation(dds7, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
condition <- treatments7
pcaplot <- plotPCA(vsd7, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix7, "-PCA_plot.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes7 <- head( order( rowVars( assay(rld7) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix7, "-heatmap.jpeg"))
heatmap <- heatmap.2( assay(rld7)[topVarGenes7, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()


##########################################
##########################################
# Wild Type D2_6 vs Mdx D2_1.5
##########################################
##########################################


outputPrefix8 <- "Wd_D2_6_vs_Mdx_D2_1.5"
sampleFiles8 <- c("deseq2_counts_star_out_merged_TS6_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS9_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136246_TS25_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS3_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136250_TS8_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_counts_star_out_merged_TS20_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames8 <- c("TS6","TS9","TS25","TS3","TS8","TS20")
sampleCondition8 <- c("Wd_D2_6","Wd_D2_6","Wd_D2_6","Mdx_D2_1.5","Mdx_D2_1.5","Mdx_D2_1.5")
sampleTable8 <- data.frame(sampleName = sampleNames8, fileName = sampleFiles8, condition = sampleCondition8)
treatments8 = c("Wd_D2_6","Mdx_D2_1.5")
condition = sampleCondition8
design = ~ condition
ddsHTSeq8 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable8, directory = directory, design = ~ condition)
colData(ddsHTSeq8)$condition <- factor(colData(ddsHTSeq8)$condition,levels = treatments8)
#guts
dds8 <- DESeq(ddsHTSeq8)
res8 <- results(dds8)
summary(res8)
# order results by padj value (most significant to least)
res8 = subset(res8, padj<0.05)
res8 <- res8[order(res8$padj),]
res8
summary(res8)
# How many adjusted p-values were less than 0.1?
sum(res8$padj < 0.1, na.rm=TRUE)
sum(res8$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered8 <- res8[order(res8$pvalue),]

# save data results and normalized reads to csv
resdata8 <- merge(as.data.frame(res8), as.data.frame(counts(dds8,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata8)[1] <- 'ensembl_gene_id'
write.csv(resdata8, file = paste0(outputPrefix8, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
# resSig3 <- subset(resOrdered3, padj < 0.1)
# resSig3
# summary(resSig3)
# write.csv(as.data.frame(resSig3), file= paste0(outputPrefix,"-results-with-padj.1.csv"))
# 
# resSig_05_3 <- subset(resOrdered3, padj < 0.05)
# resSig_05_3
# summary(resSig_05_3)
# write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))
# # Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "Wd_D2_6_vs_Mdx_D2_1.5-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Wd_D2_6_vs_Mdx_D2_1.5-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix8, "-normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red
plotMA(dds8, ylim=c(-5,5),main = "Wd_D2_6_vs_Mdx_D2_1.5")
dev.copy(png, paste0(outputPrefix8, "-MAplot.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld8 <- rlogTransformation(dds8, blind=T)
vsd8 <- varianceStabilizingTransformation(dds8, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
# B10_vs_B10_6
condition <- treatments8
pcaplot <- plotPCA(vsd8, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix8, "-PCA_plot.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes8 <- head( order( rowVars( assay(rld8) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix8, "-heatmap.jpeg"))
heatmap <- heatmap.2( assay(rld8)[topVarGenes8, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()


