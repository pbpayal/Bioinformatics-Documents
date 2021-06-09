source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# Load DESeq2 library
library("DESeq2")

# Set directory
directory <- "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/James/Deseq2_counts/"
setwd(directory)

# Set the prefix for each output file name
##########################################
##########################################
# B10_1.5 vs D2_1.5
##########################################
##########################################

outputPrefix <- "B10_1.5_vs_D2_1.5_"
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
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)
# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
res
summary(res)
mcols(res, use.names=TRUE)
# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
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
# followed by the write.csv function.
resSig <- subset(resOrdered, padj < 0.1)
resSig
summary(resSig)
write.csv(as.data.frame(resSig), file= paste0(outputPrefix,"-results-with-padj.1.csv"))

resSig_05 <- subset(resOrdered, padj < 0.05)
resSig_05
summary(resSig_05)
write.csv(as.data.frame(resSig_05), file= paste0(outputPrefix,"-results-with-padj.05.3.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv <- "B10_1.5_vs_D2_1.5_-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5_vs_D2_1.5_-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix, "normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# B10_vs_D2_1.5
plotMA(dds, ylim=c(-5,5),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
condition <- treatments
pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))

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


outputPrefix2 <- "B10_1.5 vs B10_6_"
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
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory, design = ~ condition)
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
#guts
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
# order results by padj value (most significant to least)
res2 = subset(res2, padj<0.05)
res2 <- res2[order(res2$padj),]
res2
summary(res2)
# How many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
sum(res2$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered2 <- res2[order(res2$pvalue),]

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
results_csv <- "B10_1.5 vs B10_6_-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_1.5 vs B10_6_-results-with-normalized.txt"
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
write.csv(as.data.frame(m2),file = paste0(outputPrefix2, "normalized-Biomart-results_merged.csv"))


####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

plotMA(dds2, ylim=c(-5,5),main = "B10_1.5 vs B10_6")
dev.copy(png, paste0(outputPrefix2, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

# Log-Transformation
rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot2,file=paste0(outputPrefix2, "-ggplot2.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes2 <- head( order( rowVars( assay(rld2) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix2, "-heatmap3.jpeg"))
heatmap <- heatmap.2( assay(rld2)[topVarGenes2, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()


##########################################
##########################################
# B10_6 vs D2_6
##########################################
##########################################


outputPrefix3 <- "B10_6 vs D2_6_"
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
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable3, directory = directory, design = ~ condition)
colData(ddsHTSeq3)$condition <- factor(colData(ddsHTSeq3)$condition,levels = treatments3)
#guts
dds3 <- DESeq(ddsHTSeq3)
res3 <- results(dds3)
# order results by padj value (most significant to least)
res3 = subset(res3, padj<0.05)
res3 <- res3[order(res3$padj),]
res3
summary(res3)
# How many adjusted p-values were less than 0.1?
sum(res3$padj < 0.1, na.rm=TRUE)
sum(res3$padj < 0.05, na.rm=TRUE)
# # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
# res_05 <- results(dds, alpha=0.05)
# summary(res_05)

# We can order our results table by the smallest p value
resOrdered3 <- res3[order(res3$pvalue),]

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
results_csv <- "B10_6 vs D2_6_-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "B10_6 vs D2_6_-results-with-normalized.txt"
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
rld3 <- rlogTransformation(dds3, blind=T)
vsd3 <- varianceStabilizingTransformation(dds3, blind=T)

# PCA Plot
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
condition <- treatments3
pcaplot3 <- plotPCA(vsd3, intgroup=c("condition")) + geom_text(aes(label=name), vjust=5)
ggsave(pcaplot3,file=paste0(outputPrefix3, "PCA_plot.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes3 <- head( order( rowVars( assay(rld3) ), decreasing=TRUE ), 35 )
dev.copy(jpeg, paste0(outputPrefix3, "-heatmap.jpeg"))
heatmap3 <- heatmap.2( assay(rld3)[topVarGenes3, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()

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


