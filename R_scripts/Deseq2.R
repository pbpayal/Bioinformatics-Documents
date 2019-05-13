#If Deseq2 is not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
# Load DESeq2 library, if Deseq2 is already installed
library("DESeq2")

#Please check before installing if you have installed devtools, if not then only do the following otherwise just load the library
#install.packages("devtools")
library(devtools)

# Cell_Line1
directory <- "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Katrina_Cell2_redo"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Cell2_DESeq2"

sampleFiles<- c("OC1_S17.txt",
                "OC2_S18.txt",
                "OC3_S19.txt",
                "OC4_S20.txt",
                "OE1_S21.txt",
                "OE2_S22.txt",
                "OE3_S23.txt",
                "OE4_S24.txt")

sampleNames <- c("OC1","OC2","OC3","OC4","OE1","OE2","OE3","OE4")
sampleCondition <- c("control","control","control","control","exp","exp","exp","exp")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("control","exp")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
                                       
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)

#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)

# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
res
# Result for res
# log2 fold change (MLE): condition exp vs control 
# Wald test p-value: condition exp vs control 
# DataFrame with 53465 rows and 6 columns

summary(res)
# adjusted p-value < 0.1
# LFC > 0 (up)     : 52, 0.34% 
# LFC < 0 (down)   : 24, 0.16% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 10004, 66% 
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]
#write.csv(resOrdered, file = paste0(outputPrefix, "-results-ordered.csv"))

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)

# p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
res05 <- results(dds, alpha=0.05)
summary(res05)

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

resSig05 <- subset(resOrdered, padj < 0.05)
resSig05
summary(resSig05)
write.csv(as.data.frame(resSig05), file= paste0(outputPrefix,"-results-with-padj.05.csv"))


# send normalized counts to tab delimited file for GSEA, etc.
#write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
#mcols(res, use.names = T)
#write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

# ddsClean <- replaceOutliersWithTrimmedMean(dds)
# ddsClean <- DESeq(ddsClean)
# tab <- table(initial = results(dds)$padj < 0.05,
#              cleaned = results(ddsClean)$padj < 0.05)
# addmargins(tab)
# write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
# resClean <- results(ddsClean)
# resClean = subset(res, padj<0.05)
# resClean <- resClean[order(resClean$padj),]
# write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# Load biomaRt library
library("biomaRt")

# Convert final results .csv file into .txt file
results_csv <- "Cell1_DESeq2-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Cell1_DESeq2-results-with-normalized.txt"

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)

# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)

# listFilters()
# listAttributes()
b <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id","entrezgene", "description"), values=a$ensembl_gene_id, mart= ensembl)

#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"

m <- merge(a, b, by="ensembl_gene_id")

write.csv(as.data.frame(m),file = "CellLine1_DEseq2-Biomart-results_merged.csv")

####################################################################################
# Cell_Line2
directory <- "/Users/pbanerjee/Desktop/Katrina_counts/Cell_Line2/Deliverables/Raw_counts_cell_line2"
setwd(directory)
# Set the prefix for each output file name
outputPrefix2 <- "Cell2_DESeq2"
sampleFiles2 <- c("OC1_S17.txt",
                "OC2_S18.txt",
                "OC3_S19.txt",
                "OC4_S20.txt",
                "OE1_S21.txt",
                "OE2_S22.txt",
                "OE3_S23.txt",
                "OE4_S24.txt")
sampleNames2 <- c("OC1","OC2","OC3","OC4","OE1","OE2","OE3","OE4")
sampleCondition2 <- c("control","control","control","control","exp","exp","exp","exp")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition = sampleCondition2)
treatments2 = c("control","exp")
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory, design = ~ condition)
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
#guts
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
# order results by padj value (most significant to least)
res2= subset(res2, padj<0.05)
res2 <- res2[order(res2$padj),]
res2
summary(res2)
# How many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
sum(res2$padj < 0.05, na.rm=TRUE)
# p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
res2_05 <- results(dds2, alpha=0.05)
summary(res2_05)

# We can order our results table by the smallest p value
resOrdered2 <- res2[order(res2$pvalue),]

# save data results and normalized reads to csv
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, file = paste0(outputPrefix2, "-results-with-normalized.csv"))
# Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
# followed by the write.csv function.
resSig2 <- subset(resOrdered2, padj < 0.1)
resSig2
summary(resSig2)
write.csv(as.data.frame(resSig2), file= paste0(outputPrefix2,"-results-with-padj.1.csv"))

resSig2_05 <- subset(resOrdered2, padj < 0.05)
resSig2_05
summary(resSig2_05)
write.csv(as.data.frame(resSig2_05), file= paste0(outputPrefix2,"-results-with-padj.05.3.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
# Convert final results .csv file into .txt file
results_csv2 <- "Cell2_DESeq2-results-with-normalized.csv"
write.table(read.csv(results_csv2), gsub(".csv",".txt",results_csv2))
results_txt2 <- "Cell2_DESeq2-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt2, head=TRUE)
# listFilters()
# listAttributes()
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id","entrezgene", "description"), values=a2$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "Cell2_DEseq2-Biomart-results_merged.csv")


####################################################################################
####################################################################################
####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red
# Cell_Line1
plotMA(dds, ylim=c(-5,5),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()
# Cell_Line2
plotMA(dds2, ylim=c(-5,5),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix2, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
#Cell_Line1
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)
#Cell_Line2
rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
# Cell_Line1
condition <- treatments
pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=2)
ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))
# Cell_Line2
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name),vjust=2)
ggsave(pcaplot2,file=paste0(outputPrefix2, "-ggplot2.pdf"))

#If RColorBrewer not insatlled first install, if already installed just load library
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer")
#   library(RColorBrewer)
# }

# heatmap of data
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
# 1000 top expressed genes with heatmap.2
#Cell_Line1
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:100]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))
dev.off()


# Pheatmap
# install.packages("pheatmap", "RColorBrewer", "viridis")
library(pheatmap)
# Cell_line1
# select the top 20 genes with the lowest padj
topgenes <- head(rownames(resOrdered),20)
mat <- assay(rld)[topgenes,]
# Subtract the row means from each value
mat <- mat - rowMeans(mat)
# Optional, but to make the plot nicer:
#df <- as.data.frame(colData(dds))
df <- as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition" # Rename the column header
rownames(df) = colnames(mat) # add rownames
# and plot the actual heatmap
pheatmap(mat, annotation_col=df)

# Cell_line2
topgenes2 <- head(rownames(resOrdered2),3)
mat2 <- assay(rld)[topgenes2,]
mat2 <- mat2 - rowMeans(mat2)
#df2 <- as.data.frame(colData(dds2))
df2 <- as.data.frame(colData(rld2)[,c("condition")])
pheatmap(mat2, annotation_col=df2)




# save normalized values
# write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t'))
# write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))



# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
