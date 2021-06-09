# If Deseq2 is not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
# Load DESeq2 library, if Deseq2 is already installed
library("DESeq2")

#Please check before installing if you have installed devtools, if not then only do the following otherwise just load the library
#install.packages("devtools")
library(devtools)

# Pig
# Control(PN) Vs hypoxia and treated with EGFR inhibitor (PE)
# set diretory to where your counts files are
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN3_lane1_S21_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE1_lane1_S13_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames <- c("PN1","PN2","PN3","PE1","PE2","PE3")
sampleCondition <- c("control","control","control","exp","exp","exp")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("control","exp")

# Load the data in Deseq format
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)

#Perform the Differential Gene Expression using DESeq
dds <- DESeq(ddsHTSeq)
# Normalization
dds_norm <- estimateSizeFactors(dds)
sizeFactors(dds_norm)
normalized_counts <- counts(dds_norm, normalized=TRUE)
# save this normalized data matrix to file for later use:
write.table(normalized_counts, file="PN_vs_PE-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
# models the normalization inside the Generalized Linear Model (GLM). 
# These normalized counts will be useful for downstream visualization of results, 
# but cannot be used as input to DESeq2 or any other tools that peform 
# differential expression analysis which use the negative binomial model.

res <- results(dds)
res
summary(res)
# out of 18980 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 7, 0.037%
# LFC < 0 (down)     : 5, 0.026%

# order results by padj value (most significant to least)
res_05 = subset(res, padj<0.05)
res_05 <- res_05[order(res_05$padj),]
res_05
summary(res_05)
# out of 8 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4, 50%
# LFC < 0 (down)     : 4, 50%

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'ensembl_gene_id'
write.csv(resdata, "PN_vs_PE-results.csv")

# Install biomaRt package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# Load biomaRt library
library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PE-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PE-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PE-results_merged.csv")



####################################################################################
####################################################################################
# Visualisation
####################################################################################
####################################################################################
# 1) MA plot
# 2) PCA plot
# 3) Heatmap

# Packages that need to be installed
# If RColorBrewer not insatlled first install, if already installed just load library
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer")
#   library(RColorBrewer)
# }
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")

# MA plot
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# PN_vs_PE
plotMA(dds, ylim=c(-5,5),main = "Control VS Hypoxia with EGFR")
dev.copy(png, "PN_vs_PE-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition <- treatments
pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name), vjust=2)
pcaplot
ggsave(pcaplot,filename = "PN_vs_PE_PCA_plot.pdf")

library("pheatmap")
mat = assay(vsd)[ head(order(res$padj),30), ] # select the top 30 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df = as.data.frame(colData(vsd)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df) = "condition" # Rename the column header
rownames(df) = colnames(mat) # add rownames
# and plot the actual heatmap
pheatmap(mat, annotation_col=df, main = "Top 30 Gene Expression")

# Heatmap
# The assay function is used to extract the matrix of normalized values.
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 35 )
x <- assay(vsd)[topVarGenes,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd)[topVarGenes,],
          #cellnote = mat_data,  # same data set for cell labels
          main = "35 Top Expressed Genes", # heat map title
          #notecol="blue",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,16),     # widens margins around plot
          #          col=my_palette,       # use on color palette defined earlier
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA"            # turn off column clustering
)
dev.copy(png,"PN_vs_PE_HEATMAP.png")
dev.off()

# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000179094", intgroup="condition")
plotCounts(dds, gene="ENSG00000116584", intgroup="condition")






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
