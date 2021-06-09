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

setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Deseq2_files/")

# ALL Matrix together
count_matrix_pig <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/deseq2_data_matrix.csv', header =T, sep = ",")
sample_annotation_pig <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/sample_annotation_pig.csv", header =T, sep = ",")
# Remove sample PN3
sample_anno_subset <- sample_annotation_pig[-c(3,4), -c(3,4)]
row.names(count_matrix_pig) <- count_matrix_pig$GeneID
count_matrix_pig <- subset(count_matrix_pig, select = -c(GeneID,N3,E1))
sampleCondition_pig <- c("Normal","Normal",
                         "Hypoxia","Hypoxia","Hypoxia",
                     "Experiment","Experiment")
condition = sampleCondition_pig
# PN3 not included
dds_pig <- DESeqDataSetFromMatrix(countData = count_matrix_pig,
                              colData = sample_anno_subset, 
                              design = ~ condition)

# PCA Plot to see sample distribution
# same results with and without DESeq
dds_pig_rlog <- rlog(dds_pig)
dds_pig_vst <- vst(dds_pig)
mat_pig_all = assay(dds_pig_vst)
mat_pig_all = mat_pig_all - rowMeans(mat_pig_all)

write.table(vst,"pig_vst.txt",sep = "\t")
pcaplot_pig <- plotPCA(dds_pig_rlog, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=3)
pcaplot_pig
ggsave(pcaplot_pig,filename = "PCA_plot_pig_no_PN3.png")

# Differential Gene Expression
dds <- DESeq(dds_pig)
res <- results(dds)
res
summary(res)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 63, 0.33%
# LFC < 0 (down)     : 118, 0.62%

# Extract the normalized counts
dds_size_factors = estimateSizeFactors(dds_pig)
sizeFactors(dds_size_factors)
norm_counts = counts(dds, normalized = TRUE)
write.table(norm_counts, "Normalized_counts_all_pig.txt", sep="\t")

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DEGreport/inst/doc/DEGreport.html
# degQC(count_matrix_pig,  pvalue = res[["pvalue"]])

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'ensembl_gene_id'
write.csv(resdata, "All_pig_DEG-results.csv")
 


# Heatmap
# http://folk.uio.no/jonbra/MBV-INF4410_2017/exercises/2017-12-07_R_DESeq2_exercises_without_results.html
vsd <- varianceStabilizingTransformation(dds, blind=T)
library("pheatmap")
mat = assay(vsd)
#mat = assay(vsd)[ head(order(res$padj),900), ] # select the top 30 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df = as.data.frame(colData(vsd)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df) = "condition" # Rename the column header
rownames(df) = colnames(mat) # add rownames
# and plot the actual heatmap
pheatmap(mat, annotation_col=df, main = "Top 50 Gene Expression", 
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T)
write.table(mat, "heatmap_data2.txt", sep="\t")

library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")
a <- as.data.frame(mat_pig_all)
a$Ensembl_Id <- rownames(mat_pig_all)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$Ensembl_Id, mart= ensembl)
colnames(a)[colnames(a)=="Ensembl_Id"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "Heatmap_data_all_genes_vst.csv")



# Final heatmap
heatmap_data = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Heatmaps/heatmap_data_edited_Surajit.csv", sep = ",")
heatmap_data <- heatmap_data[c(1,4,7,11)]
row.names(heatmap_data) <- heatmap_data$X
heatmap_data <- subset(heatmap_data, select = -c(X))
heatmap_data <- as.matrix(heatmap_data)
# rename
colnames(heatmap_data)[1] <- "Normoxic"
colnames(heatmap_data)[2] <- "EGFRi"
colnames(heatmap_data)[3]<- "Hypoxic-Ischemic"
# reorder
heatmap_reorder <- heatmap_data[, c("Normoxic","Hypoxic-Ischemic","EGFRi")]
pheatmap(heatmap_reorder, main = "Top 800 Gene Expression", 
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-2, +2, length = 101))
# pheatmap(heatmap_data[1:600,], main = "Top 600 Gene Expression", 
#          show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
#          breaks = seq(-2, +2, length = 101))
# pheatmap(heatmap_data[1:100,], main = "Top 100 Gene Expression", 
#          show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T)
head(brewer.pal.info)
table(brewer.pal.info$category)
brewer.pal(11, "RdYlGn")

head(heatmap_data[, c("Normal","Hypoxia","Experiment")])


# Sub Heatmaps
synaptic_transmission_data = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Heatmaps/Heatmaps_2/synaptic_PNvsPH_heatmapper.txt", sep = "\t", header = TRUE)
syn_data_ensg <- synaptic_transmission_data[c(2,3,4)]
row.names(syn_data_ensg) <- syn_data_ensg$NAME
syn_data_ensg <- subset(syn_data_ensg, select = -c(NAME))
colnames(syn_data_ensg)[1] <- "Normoxic"
colnames(syn_data_ensg)[2] <- "Hypoxic-Ischemic"
syn_data_ensg <- as.matrix(syn_data_ensg)


pheatmap(syn_data_ensg, main = "Synaptic Transmission",
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-1, +1, length = 101))

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("yellow","white","blue"))
# my_palette <- colorRampPalette(c("yellow","blue"))


# # (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-2,-0.1,length=200),
               seq(0,0.1,length=200),
               seq(0.2,2,length=200))

# col_breaks = c(seq(-2,0.1,length=200),
#                seq(0.2,2,length=200))
#using a red and blue colour scheme without traces and scaling each row

library("RColorBrewer")
display.brewer.all()
brewer.pal.info
#col=brewer.pal(9,"YIGnBu")

heatmap.2(syn_data_ensg,
#          cellnote = mat_data,  # same data set for cell labels
          main = "Synaptic Transmission", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,9),     # widens margins around plot
#          col=my_palette,       # use on color palette defined earlier
          col=brewer.pal(9,"Blues"),
#          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          key = T,
#          symm=F,symkey=F,symbreaks=T, scale="none",
          cexRow = 0.8, 
          cexCol = 1, srtCol=0,adjCol = c(0.5,0.5),
          keysize=1, key.par = list(cex=0.5)) 

nervous_system_data = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Heatmaps/Heatmaps_2/nervous_system_heatmapper_PNvs_PH.txt", sep = "\t", header = TRUE)
nerv_data_ensg <- nervous_system_data[c(2,3,4)]
row.names(nerv_data_ensg) <- nerv_data_ensg$NAME
nerv_data_ensg <- subset(nerv_data_ensg, select = -c(NAME.Development))
colnames(nerv_data_ensg)[1] <- "Normoxic"
colnames(nerv_data_ensg)[2] <- "Hypoxic-Ischemic"
nerv_data_ensg <- as.matrix(nerv_data_ensg)


pheatmap(nerv_data_ensg, main = "Synaptic Transmission",
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-1, +1, length = 101))

# # creates a own color palette from red to green
# my_palette <- colorRampPalette(c("yellow","white","blue"))
# # my_palette <- colorRampPalette(c("yellow","blue"))


# # # (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(-2,-0.1,length=200),
#                seq(0,0.1,length=200),
#                seq(0.2,2,length=200))
# 
# # col_breaks = c(seq(-2,0.1,length=200),
# #                seq(0.2,2,length=200))
# #using a red and blue colour scheme without traces and scaling each row

library("RColorBrewer")
display.brewer.all()
brewer.pal.info
#col=brewer.pal(9,"YIGnBu")

heatmap.2(nerv_data_ensg,
          #          cellnote = mat_data,  # same data set for cell labels
          main = "Nervous System Development", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,9),     # widens margins around plot
          #          col=my_palette,       # use on color palette defined earlier
          col=brewer.pal(9,"YlOrRd"),
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          key = T,
          #          symm=F,symkey=F,symbreaks=T, scale="none",
          cexRow = 0.8, 
          cexCol = 1, srtCol=0, adjCol = c(0.5,0.5),
          keysize=1, key.par = list(cex=0.5))



syn_all_data = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Heatmaps/Heatmaps_2/synaptic_transmission_heatmap_genes_ALL.txt", sep = "\t", header = TRUE)
syn_all_data <- syn_all_data[c(1,3,4,5)]
row.names(syn_all_data) <- syn_all_data$NAME
syn_all_data <- subset(syn_all_data, select = -c(NAMES))
colnames(syn_all_data)[1] <- "Normoxic"
colnames(syn_all_data)[2] <- "Hypoxic-Ischemic"
colnames(syn_all_data)[3] <- "EGFRi"
syn_all_data <- as.matrix(syn_all_data)

heatmap.2(syn_all_data,
          #          cellnote = mat_data,  # same data set for cell labels
          main = "Synaptic Transmission", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,9),     # widens margins around plot
          #          col=my_palette,       # use on color palette defined earlier
          col=brewer.pal(11,"Blues"),
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          key = T,
          #          symm=F,symkey=F,symbreaks=T, scale="none",
          cexRow = 0.8, 
          cexCol = 1, srtCol=0, adjCol = c(0.5, 0.5),
          keysize=1, key.par = list(cex=0.5))

nerv_all_data = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Heatmaps/Heatmaps_2/nervous_transmission_heatmap_genes_ALL.txt", sep = "\t", header = TRUE)
nerv_all_data <- nerv_all_data[c(1,3,4,5)]
row.names(nerv_all_data) <- nerv_all_data$NAMES
nerv_all_data <- subset(nerv_all_data, select = -c(NAMES))
colnames(nerv_all_data)[1] <- "Normoxic"
colnames(nerv_all_data)[2] <- "Hypoxic-Ischemic"
colnames(nerv_all_data)[3] <- "EGFRi"
nerv_all_data <- as.matrix(nerv_all_data)

heatmap.2(nerv_all_data,
          #          cellnote = mat_data,  # same data set for cell labels
          main = "Nervous System Development", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,9),     # widens margins around plot
          #          col=my_palette,       # use on color palette defined earlier
          col=brewer.pal(9,"YlOrRd"),
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          key = T,
          #          symm=F,symkey=F,symbreaks=T, scale="none",
          cexRow = 0.8, 
          cexCol = 1, srtCol=0, adjCol = c(0.5, 0.5),
          keysize=1, key.par = list(cex=0.5))

# Pig
# Control(PN) Vs Hypoxia and treated with EGFR inhibitor (PE)
# set diretory to where your counts files are
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Deseq2_files"
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
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description","sscrofa_homolog_associated_gene_name" ), values=a$ensembl_gene_id, mart= ensembl)
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


################################################################################################
# PN_vs_PE_n2
sampleFiles2<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE1_lane1_S13_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames2 <- c("PN1","PN2","PE1","PE2")
sampleCondition2 <- c("control","control","exp","exp")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition2 = sampleCondition2)
treatments2 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory, design = ~ condition2)

colData(ddsHTSeq2)$condition2 <- factor(colData(ddsHTSeq2)$condition2,levels = treatments2)

#Perform the Differential Gene Expression using DESeq
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
res2
summary(res2)


# order results by padj value (most significant to least)
res2_05 = subset(res2, padj<0.05)
res2_05 <- res2_05[order(res2_05$padj),]
res2_05
summary(res_05)


# save data results and normalized reads to csv
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, "PN_vs_PE-n2-results.csv")

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PE-n2-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PE-n2-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PE-n2-results_merged.csv")

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

# PN_vs_PE_n2
plotMA(dds2, ylim=c(-5,5),main = "Control VS Acute Hypoxia N2")
dev.copy(png, "PN_vs_PE_n2-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
# PN_vs_PE_n2
condition2 <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition2")) + geom_text(aes(label=name), vjust=2)
pcaplot2
ggsave(pcaplot2,filename = "PN_vs_PE_PCA_n2_plot.pdf")

library("pheatmap")
mat2 = assay(vsd2)[ head(order(res2$padj),30), ] # select the top 30 genes with the lowest padj
mat2 = mat2 - rowMeans(mat2) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df2 = as.data.frame(colData(vsd2)[,c("condition2")]) # Create a dataframe with a column of the conditions
colnames(df2) = "condition" # Rename the column header
rownames(df2) = colnames(mat2) # add rownames
# and plot the actual heatmap
pheatmap(mat2, annotation_col=df2, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes2 <- head( order( rowVars( assay(vsd2) ), decreasing=TRUE ), 35 )
x2 <- assay(vsd2)[topVarGenes2,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd2)[topVarGenes2,],
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
dev.copy(png,"PN_vs_PE-n2_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res2, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res2, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res2, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000179094", intgroup="condition")
plotCounts(dds, gene="ENSG00000116584", intgroup="condition")




# Control(PN) Vs Acute hypoxia(PH)
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles3<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN3_lane1_S21_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames3 <- c("PN1","PN2","PN3","PH1","PH2","PH3")
sampleCondition3 <- c("control","control","control","exp","exp","exp")
sampleTable3 <- data.frame(sampleName3 = sampleNames3, fileName = sampleFiles3, condition3 = sampleCondition3)
treatments3 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable3, directory = directory, design = ~ condition3)

colData(ddsHTSeq3)$condition3 <- factor(colData(ddsHTSeq3)$condition3,levels = treatments3)

#Perform the Differential Gene Expression using DESeq
dds3 <- DESeq(ddsHTSeq3)
res3 <- results(dds3)
res3
summary(res3)
# out of 19050 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 9, 0.047%
# LFC < 0 (down)     : 9, 0.047%

# order results by padj value (most significant to least)
res3_05 = subset(res3, padj<0.05)
res3_05 <- res3_05[order(res3_05$padj),]
res3_05
summary(res3_05)
# out of 10 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4, 40%
# LFC < 0 (down)     : 6, 60%

# save data results and normalized reads to csv
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata3)[1] <- 'ensembl_gene_id'
write.csv(resdata3, "PN_vs_PH-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PH-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PH-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PH-results_merged.csv")


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
plotMA(dds3, ylim=c(-5,5),main = "Control VS Acute Hypoxia")
dev.copy(png, "PN_vs_PH-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd3 <- varianceStabilizingTransformation(dds3, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition3 <- treatments3
pcaplot3 <- plotPCA(vsd3, intgroup=c("condition3")) + geom_text(aes(label=name), vjust=2)
pcaplot3
ggsave(pcaplot3,filename = "PN_vs_PH_PCA_plot.pdf")

library("pheatmap")
mat3 = assay(vsd3)[ head(order(res3$padj),30), ] # select the top 30 genes with the lowest padj
mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df3 = as.data.frame(colData(vsd3)[,c("condition3")]) # Create a dataframe with a column of the conditions
colnames(df3) = "condition" # Rename the column header
rownames(df3) = colnames(mat3) # add rownames
# and plot the actual heatmap
pheatmap(mat3, annotation_col=df3, main = "Top 30 Gene Expression")



# Heatmap
# select 35 genes with the highest variance across samples
topVarGenes3 <- head( order( rowVars( assay(vsd3) ), decreasing=TRUE ), 35 )
x3 <- assay(vsd3)[topVarGenes3,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd3)[topVarGenes3,],
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
dev.copy(png,"PN_vs_PH_Heatmap.png")
dev.off()

# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res3, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res3, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res3, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# Control(PN) Vs Acute Hypoxia (PH)
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles4<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames4 <- c("PN1","PN2","PH1","PH2")
sampleCondition4 <- c("control","control","exp","exp")
sampleTable4 <- data.frame(sampleName = sampleNames4, fileName = sampleFiles4, condition4 = sampleCondition4)
treatments4 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable4, directory = directory, design = ~ condition4)

colData(ddsHTSeq4)$condition4 <- factor(colData(ddsHTSeq4)$condition4,levels = treatments4)

#Perform the Differential Gene Expression using DESeq
dds4 <- DESeq(ddsHTSeq4)
res4 <- results(dds4)
res4
summary(res4)
# out of 18643 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 233, 1.2%
# LFC < 0 (down)     : 187, 1%

# order results by padj value (most significant to least)
res4_05 = subset(res4, padj<0.05)
res4_05 <- res4_05[order(res4_05$padj),]
res4_05
summary(res4_05)
# out of 10 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4, 40%
# LFC < 0 (down)     : 6, 60%

# save data results and normalized reads to csv
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata4)[1] <- 'ensembl_gene_id'
write.csv(resdata4, "PN_vs_PH_n2-results.csv")
library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PH_n2-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PH_n2-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PH_n2-results_merged.csv")

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

# PN_vs_PH-n2
plotMA(dds4, ylim=c(-5,5),main = "Control VS Acute Hypoxia")
dev.copy(png, "PN_vs_PH-n2-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
#rld4 <- rlogTransformation(dds4, blind=T)
vsd4 <- varianceStabilizingTransformation(dds4, blind=T)

# PCA Plot
condition4 <- treatments4
pcaplot4 <- plotPCA(vsd4, intgroup=c("condition4")) + geom_text(aes(label=name), vjust=2)
pcaplot4
ggsave(pcaplot4,filename = "PN_vs_PH_n2_PCA_plot.pdf")

# Pheatmap
library("pheatmap")
mat = assay(rld)[ head(order(res$padj),30), ] # select the top 30 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df = as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df) = "condition" # Rename the column header
rownames(df) = colnames(mat) # add rownames
# and plot the actual heatmap
pheatmap(mat, annotation_col=df)

# Heatmap
# The assay function is used to extract the matrix of normalized values.
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes <- head( order( rowVars( assay(vsd4) ), decreasing=TRUE ), 35 )
x <- assay(vsd4)[topVarGenes,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd4)[topVarGenes,],
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
dev.copy(png,"PN_vs_PH-n2_HEATMAP.png")
dev.off()

# Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res4, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 



# Acute Hypoxia (PH) vs Hypoxia and treated with EGFR inhibitor (PE) 
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles5<- c("deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE1_lane1_S13_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames5 <- c("PH1","PH2","PH3","PE1","PE2","PE3")
sampleCondition5 <- c("control","control","control","exp","exp","exp")
sampleTable5 <- data.frame(sampleName = sampleNames5, fileName = sampleFiles5, condition5 = sampleCondition5)
treatments5 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq5 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable5, directory = directory, design = ~ condition5)

colData(ddsHTSeq5)$condition5 <- factor(colData(ddsHTSeq5)$condition5,levels = treatments5)

#Perform the Differential Gene Expression using DESeq
dds5 <- DESeq(ddsHTSeq5)
res5 <- results(dds5)
res5
summary(res5)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 9, 0.048%
# LFC < 0 (down)     : 4, 0.021%



# save data results and normalized reads to csv
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata5)[1] <- 'ensembl_gene_id'
write.csv(resdata5, "PH_vs_PE-results.csv")


library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PH_vs_PE-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PH_vs_PE-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PH_vs_PE-results_merged.csv")

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

plotMA(dds5, ylim=c(-5,5),main = "Control VS Acute Hypoxia")
dev.copy(png, "PH_vs_PE_MA_plot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
#rld4 <- rlogTransformation(dds4, blind=T)
vsd5 <- varianceStabilizingTransformation(dds5, blind=T)

# PCA Plot
condition5 <- treatments5
pcaplot5 <- plotPCA(vsd5, intgroup=c("condition5")) + geom_text(aes(label=name), vjust=2)
pcaplot5
ggsave(pcaplot5,filename = "PH_vs_PE_PCA_plot.pdf")

# Pheatmap
library("pheatmap")
mat5 = assay(vsd5)[ head(order(res5$padj),30), ] # select the top 30 genes with the lowest padj
mat5 = mat5 - rowMeans(mat5) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df5 = as.data.frame(colData(vsd5)[,c("condition5")]) # Create a dataframe with a column of the conditions
colnames(df5) = "condition" # Rename the column header
rownames(df5) = colnames(mat5) # add rownames
# and plot the actual heatmap
pheatmap(mat5, annotation_col=df5)

# Heatmap
# The assay function is used to extract the matrix of normalized values.
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes5 <- head( order( rowVars( assay(vsd5) ), decreasing=TRUE ), 35 )
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd5)[topVarGenes5,],
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
dev.copy(png,"PH_vs_PE_HEATMAP.png")
dev.off()

# Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res4, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

# PN_vs_PE_n2_1
sampleFiles6<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE1_lane1_S13_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames6 <- c("PN1","PN2","PE1","PE3")
sampleCondition6 <- c("control","control","exp","exp")
sampleTable6 <- data.frame(sampleName = sampleNames6, fileName = sampleFiles6, condition6 = sampleCondition6)
treatments6 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq6 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable6, directory = directory, design = ~ condition6)

colData(ddsHTSeq6)$condition6 <- factor(colData(ddsHTSeq6)$condition6,levels = treatments6)

#Perform the Differential Gene Expression using DESeq
dds6 <- DESeq(ddsHTSeq6)
res6 <- results(dds6)
res6
summary(res6)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 22, 0.12%
# LFC < 0 (down)     : 11, 0.059%


# save data results and normalized reads to csv
resdata6 <- merge(as.data.frame(res6), as.data.frame(counts(dds6,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata6)[1] <- 'ensembl_gene_id'
write.csv(resdata6, "PN_vs_PE-n2_1-results.csv")

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PE-n2_1-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PE-n2_1-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PE-n2_1-results_merged.csv")

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

# PN_vs_PE_n2_1
plotMA(dds6, ylim=c(-5,5),main = "Control VS Hypoxia with EGFR N2 (1)")
dev.copy(png, "PN_vs_PE_n2_1-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
# rld2 <- rlogTransformation(dds2, blind=T)
vsd6 <- varianceStabilizingTransformation(dds6, blind=T)

# PCA Plot
# PN_vs_PE_n2_1
condition6 <- treatments6
pcaplot6 <- plotPCA(vsd6, intgroup=c("condition6")) + geom_text(aes(label=name), vjust=2)
pcaplot6
ggsave(pcaplot6,filename = "PN_vs_PE_PCA_n2_1_plot.pdf")

library("pheatmap")
mat6 = assay(vsd6)[ head(order(res6$padj),30), ] # select the top 30 genes with the lowest padj
mat6 = mat6 - rowMeans(mat6) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df6 = as.data.frame(colData(vsd6)[,c("condition6")]) # Create a dataframe with a column of the conditions
colnames(df6) = "condition" # Rename the column header
rownames(df6) = colnames(mat6) # add rownames
# and plot the actual heatmap
pheatmap(mat6, annotation_col=df6, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes6 <- head( order( rowVars( assay(vsd6) ), decreasing=TRUE ), 35 )
#x2 <- assay(vsd6)[topVarGenes6,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd6)[topVarGenes6,],
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
dev.copy(png,"PN_vs_PE-n2_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res6, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res6, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000179094", intgroup="condition")
plotCounts(dds, gene="ENSG00000116584", intgroup="condition")



# PN_vs_PH_n2_1
sampleFiles7<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames7 <- c("PN1","PN2","PH1","PH3")
sampleCondition7 <- c("control","control","exp","exp")
sampleTable7 <- data.frame(sampleName = sampleNames7, fileName = sampleFiles7, condition7 = sampleCondition7)
treatments7 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq7 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable7, directory = directory, design = ~ condition7)

colData(ddsHTSeq7)$condition7 <- factor(colData(ddsHTSeq7)$condition7,levels = treatments7)

#Perform the Differential Gene Expression using DESeq
dds7 <- DESeq(ddsHTSeq7)
res7 <- results(dds7)
res7
summary(res7)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 157, 0.84%
# LFC < 0 (down)     : 112, 0.6%


# save data results and normalized reads to csv
resdata7 <- merge(as.data.frame(res7), as.data.frame(counts(dds7,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata7)[1] <- 'ensembl_gene_id'
write.csv(resdata7, "PN_vs_PH-n2_1-results.csv")

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PH-n2_1-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PH-n2_1-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PH-n2_1-results_merged.csv")

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

# PN_vs_PH_n2_1
plotMA(dds7, ylim=c(-5,5),main = "Control VS Acute Hypoxia N2 (1)")
dev.copy(png, "PN_vs_PH_n2_1-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
# rld2 <- rlogTransformation(dds2, blind=T)
vsd7 <- varianceStabilizingTransformation(dds7, blind=T)

# PCA Plot
# PN_vs_PE_n2_1
condition7 <- treatments7
pcaplot7 <- plotPCA(vsd7, intgroup=c("condition7")) + geom_text(aes(label=name), vjust=2)
pcaplot7
ggsave(pcaplot7,filename = "PN_vs_PH_PCA_n2_1_plot.pdf")

library("pheatmap")
mat7 = assay(vsd7)[ head(order(res7$padj),30), ] # select the top 30 genes with the lowest padj
mat7 = mat7 - rowMeans(mat7) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df7 = as.data.frame(colData(vsd7)[,c("condition7")]) # Create a dataframe with a column of the conditions
colnames(df7) = "condition" # Rename the column header
rownames(df7) = colnames(mat7) # add rownames
# and plot the actual heatmap
pheatmap(mat7, annotation_col=df7, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes7 <- head( order( rowVars( assay(vsd7) ), decreasing=TRUE ), 35 )
#x2 <- assay(vsd6)[topVarGenes6,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd7)[topVarGenes7,],
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
dev.copy(png,"PN_vs_PE-n2_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res7, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res7, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res7, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000179094", intgroup="condition")
plotCounts(dds, gene="ENSG00000116584", intgroup="condition")


# PH_vs_PE_n2_1

sampleFiles8<- c("deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames8 <- c("PH1","PH2","PE2","PE3")
sampleCondition8 <- c("control","control","exp","exp")
sampleTable8 <- data.frame(sampleName = sampleNames8, fileName = sampleFiles8, condition8 = sampleCondition8)
treatments8 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq8 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable8, directory = directory, design = ~ condition8)
colData(ddsHTSeq8)$condition8 <- factor(colData(ddsHTSeq8)$condition8,levels = treatments8)

#Perform the Differential Gene Expression using DESeq
dds8 <- DESeq(ddsHTSeq8)
# # Normalization
# dds_norm <- estimateSizeFactors(dds)
# sizeFactors(dds_norm)
# normalized_counts <- counts(dds_norm, normalized=TRUE)
# # save this normalized data matrix to file for later use:
# write.table(normalized_counts, file="PN_vs_PE-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# # NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
# # models the normalization inside the Generalized Linear Model (GLM). 
# # These normalized counts will be useful for downstream visualization of results, 
# # but cannot be used as input to DESeq2 or any other tools that peform 
# # differential expression analysis which use the negative binomial model.

res8 <- results(dds8)
res8
summary(res8)
# out of 18510 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 311, 1.7%
# LFC < 0 (down)     : 300, 1.6%

# # order results by padj value (most significant to least)
# res_05 = subset(res, padj<0.05)
# res_05 <- res_05[order(res_05$padj),]
# res_05
# summary(res_05)
# # out of 8 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 50%
# # LFC < 0 (down)     : 4, 50%

# save data results and normalized reads to csv
resdata8 <- merge(as.data.frame(res8), as.data.frame(counts(dds8,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata8)[1] <- 'ensembl_gene_id'
write.csv(resdata8, "PH_vs_PE_n2_1-results.csv")


library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PH_vs_PE_n2_1-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PH_vs_PE_n2_1-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PH_vs_PE-n2_1-results_merged.csv")

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

plotMA(dds8, ylim=c(-5,5),main = "Hypoxia VS Treatment N2 (1)")
dev.copy(png, "PH_vs_PE_n2_1-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
# rld2 <- rlogTransformation(dds2, blind=T)
vsd8 <- varianceStabilizingTransformation(dds8, blind=T)

# PCA Plot
# PN_vs_PE_n2_1
condition8 <- treatments8
pcaplot8 <- plotPCA(vsd8, intgroup=c("condition8")) + geom_text(aes(label=name), vjust=2)
pcaplot8
ggsave(pcaplot8,filename = "PH_vs_PE_PCA_n2_1_plot.pdf")

library("pheatmap")
mat8 = assay(vsd8)[ head(order(res8$padj),30), ] # select the top 30 genes with the lowest padj
mat8 = mat8 - rowMeans(mat8) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df8 = as.data.frame(colData(vsd8)[,c("condition8")]) # Create a dataframe with a column of the conditions
colnames(df8) = "condition" # Rename the column header
rownames(df8) = colnames(mat8) # add rownames
# and plot the actual heatmap
pheatmap(mat8, annotation_col=df8, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes8 <- head( order( rowVars( assay(vsd8) ), decreasing=TRUE ), 35 )
#x2 <- assay(vsd6)[topVarGenes6,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd8)[topVarGenes8,],
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
dev.copy(png,"PH_vs_PE-n2_1_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res8, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res8, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res8, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

# PH_vs_PE_n2_2

sampleFiles9<- c("deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames9 <- c("PH1","PH3","PE2","PE3")
sampleCondition9 <- c("control","control","exp","exp")
sampleTable9 <- data.frame(sampleName = sampleNames9, fileName = sampleFiles9, condition9 = sampleCondition9)
treatments9 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq9 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable9, directory = directory, design = ~ condition9)
colData(ddsHTSeq9)$condition9 <- factor(colData(ddsHTSeq9)$condition9,levels = treatments9)

#Perform the Differential Gene Expression using DESeq
dds9 <- DESeq(ddsHTSeq9)
# # Normalization
# dds_norm <- estimateSizeFactors(dds)
# sizeFactors(dds_norm)
# normalized_counts <- counts(dds_norm, normalized=TRUE)
# # save this normalized data matrix to file for later use:
# write.table(normalized_counts, file="PN_vs_PE-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# # NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
# # models the normalization inside the Generalized Linear Model (GLM). 
# # These normalized counts will be useful for downstream visualization of results, 
# # but cannot be used as input to DESeq2 or any other tools that peform 
# # differential expression analysis which use the negative binomial model.

res9 <- results(dds9)
res9
summary(res9)
# out of 18453 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 148, 0.8%
# LFC < 0 (down)     : 168, 0.91%

# # order results by padj value (most significant to least)
# res_05 = subset(res, padj<0.05)
# res_05 <- res_05[order(res_05$padj),]
# res_05
# summary(res_05)
# # out of 8 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 50%
# # LFC < 0 (down)     : 4, 50%

# save data results and normalized reads to csv
resdata9 <- merge(as.data.frame(res9), as.data.frame(counts(dds9,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata9)[1] <- 'ensembl_gene_id'
write.csv(resdata9, "PH_vs_PE_n2_2-results.csv")


library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PH_vs_PE_n2_2-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PH_vs_PE_n2_2-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PH_vs_PE_n2_2-results_merged.csv")

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

plotMA(dds9, ylim=c(-5,5),main = "Hypoxia VS Treatment N2 (2)")
dev.copy(png, "PH_vs_PE_n2_2-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
# rld2 <- rlogTransformation(dds2, blind=T)
vsd9 <- varianceStabilizingTransformation(dds9, blind=T)

# PCA Plot
# PN_vs_PE_n2_1
condition9 <- treatments9
pcaplot9 <- plotPCA(vsd9, intgroup=c("condition9")) + geom_text(aes(label=name), vjust=2)
pcaplot9
ggsave(pcaplot9,filename = "PH_vs_PE_PCA_n2_2_plot.pdf")

library("pheatmap")
mat9 = assay(vsd9)[ head(order(res9$padj),30), ] # select the top 30 genes with the lowest padj
mat9 = mat9 - rowMeans(mat9) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df9 = as.data.frame(colData(vsd9)[,c("condition9")]) # Create a dataframe with a column of the conditions
colnames(df9) = "condition" # Rename the column header
rownames(df9) = colnames(mat9) # add rownames
# and plot the actual heatmap
pheatmap(mat9, annotation_col=df9, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes9 <- head( order( rowVars( assay(vsd9) ), decreasing=TRUE ), 35 )
#x2 <- assay(vsd6)[topVarGenes6,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd9)[topVarGenes9,],
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
dev.copy(png,"PH_vs_PE-n2_1_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res9, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res9, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res9, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

# PH_vs_PE_n2_3

sampleFiles10<- c("deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames10 <- c("PH2","PH3","PE2","PE3")
sampleCondition10 <- c("control","control","exp","exp")
sampleTable10 <- data.frame(sampleName = sampleNames10, fileName = sampleFiles10, condition10 = sampleCondition10)
treatments10 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq10 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable10, directory = directory, design = ~ condition10)
colData(ddsHTSeq10)$condition10 <- factor(colData(ddsHTSeq10)$condition10,levels = treatments10)

#Perform the Differential Gene Expression using DESeq
dds10 <- DESeq(ddsHTSeq10)
# # Normalization
# dds_norm <- estimateSizeFactors(dds)
# sizeFactors(dds_norm)
# normalized_counts <- counts(dds_norm, normalized=TRUE)
# # save this normalized data matrix to file for later use:
# write.table(normalized_counts, file="PN_vs_PE-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# # NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
# # models the normalization inside the Generalized Linear Model (GLM). 
# # These normalized counts will be useful for downstream visualization of results, 
# # but cannot be used as input to DESeq2 or any other tools that peform 
# # differential expression analysis which use the negative binomial model.

res10 <- results(dds10)
res10
summary(res10)
# out of 18292 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 926, 5.1%
# LFC < 0 (down)     : 1104, 6%


# # order results by padj value (most significant to least)
# res_05 = subset(res, padj<0.05)
# res_05 <- res_05[order(res_05$padj),]
# res_05
# summary(res_05)
# # out of 8 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 50%
# # LFC < 0 (down)     : 4, 50%

# save data results and normalized reads to csv
resdata10 <- merge(as.data.frame(res10), as.data.frame(counts(dds10,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata10)[1] <- 'ensembl_gene_id'
write.csv(resdata10, "PH_vs_PE_n2_3-results.csv")


library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PH_vs_PE_n2_3-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PH_vs_PE_n2_3-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PH_vs_PE_n2_3-results_merged.csv")

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

plotMA(dds10, ylim=c(-5,5),main = "Hypoxia VS Treatment N2 (3)")
dev.copy(png, "PH_vs_PE_n2_3-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
# rld2 <- rlogTransformation(dds2, blind=T)
vsd10 <- varianceStabilizingTransformation(dds10, blind=T)

# PCA Plot
# PN_vs_PE_n2_1
condition10 <- treatments10
pcaplot10 <- plotPCA(vsd10, intgroup=c("condition10")) + geom_text(aes(label=name), vjust=2)
pcaplot10
ggsave(pcaplot10,filename = "PH_vs_PE_PCA_n2_3_plot.pdf")

library("pheatmap")
mat10 = assay(vsd10)[ head(order(res10$padj),30), ] # select the top 30 genes with the lowest padj
mat10 = mat10 - rowMeans(mat10) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df10 = as.data.frame(colData(vsd10)[,c("condition10")]) # Create a dataframe with a column of the conditions
colnames(df10) = "condition" # Rename the column header
rownames(df10) = colnames(mat10) # add rownames
# and plot the actual heatmap
pheatmap(mat10, annotation_col=df10, main = "Top 30 Gene Expression")

# Heatmap
# select 35 genes with the highest variance across samples
# a subset of most highly variable genes. Here, for demonstration, 
# select the 35 genes with the highest variance across samples
topVarGenes10 <- head( order( rowVars( assay(vsd10) ), decreasing=TRUE ), 35 )
#x2 <- assay(vsd6)[topVarGenes6,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd10)[topVarGenes10,],
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
dev.copy(png,"PH_vs_PE-n2_1_Heatmap.png")
dev.off()


# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res10, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res10, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res10, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# Pig
# Control(PN) Vs Hypoxia and treated with EGFR inhibitor (PE)
# set diretory to where your counts files are
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles11<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE1_lane1_S13_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames11 <- c("PN1","PN2","PE1","PE2","PE3")
sampleCondition11 <- c("control","control","exp","exp","exp")
sampleTable11 <- data.frame(sampleName = sampleNames11, fileName = sampleFiles11, condition11 = sampleCondition11)
treatments11 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq11 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable11, directory = directory, design = ~ condition11)
colData(ddsHTSeq11)$condition11 <- factor(colData(ddsHTSeq11)$condition11,levels = treatments11)

#Perform the Differential Gene Expression using DESeq
dds11 <- DESeq(ddsHTSeq11)
# # Normalization
# dds_norm11 <- estimateSizeFactors(dds11)
# sizeFactors(dds_norm)
# normalized_counts <- counts(dds_norm, normalized=TRUE)
# # save this normalized data matrix to file for later use:
# write.table(normalized_counts, file="PN_vs_PE-normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# # NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
# # models the normalization inside the Generalized Linear Model (GLM). 
# # These normalized counts will be useful for downstream visualization of results, 
# # but cannot be used as input to DESeq2 or any other tools that peform 
# # differential expression analysis which use the negative binomial model.

res11 <- results(dds11)
res11
summary(res11)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 94, 0.5%
# LFC < 0 (down)     : 44, 0.24%

# # order results by padj value (most significant to least)
# res_05 = subset(res, padj<0.05)
# res_05 <- res_05[order(res_05$padj),]
# res_05
# summary(res_05)
# # out of 8 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 50%
# # LFC < 0 (down)     : 4, 50%

# save data results and normalized reads to csv
resdata11 <- merge(as.data.frame(res11), as.data.frame(counts(dds11,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata11)[1] <- 'ensembl_gene_id'
write.csv(resdata11, "PN_vs_PE_n_2_3-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PE_n_2_3-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PE_n_2_3-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PE_n_2_3-results_merged.csv")


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
plotMA(dds11, ylim=c(-5,5),main = "Control VS Treatment")
dev.copy(png, "PN_vs_PE_n2_3-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd11 <- varianceStabilizingTransformation(dds11, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition11 <- treatments11
pcaplot11 <- plotPCA(vsd11, intgroup=c("condition11")) + geom_text(aes(label=name), vjust=2)
pcaplot11
ggsave(pcaplot11,filename = "PN_vs_PE_n2_3_PCA_plot.pdf")

library("pheatmap")
mat11 = assay(vsd11)[ head(order(res11$padj),30), ] # select the top 30 genes with the lowest padj
mat11 = mat11 - rowMeans(mat11) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df11 = as.data.frame(colData(vsd11)[,c("condition11")]) # Create a dataframe with a column of the conditions
colnames(df11) = "condition" # Rename the column header
rownames(df11) = colnames(mat11) # add rownames
# and plot the actual heatmap
pheatmap(mat11, annotation_col=df11, main = "Top 30 Gene Expression")



# Heatmap
# select 35 genes with the highest variance across samples
topVarGenes11 <- head( order( rowVars( assay(vsd11) ), decreasing=TRUE ), 35 )
x11 <- assay(vsd11)[topVarGenes11,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd11)[topVarGenes11,],
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
dev.copy(png,"PN_vs_PE_n2_3_Heatmap.png")
dev.off()

# Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res11, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res11, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res11, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# PN vs PH (n2 3)
sampleFiles12<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames12 <- c("PN1","PN2","PH1","PH2","PH3")
sampleCondition12 <- c("control","control","hypox","hypox","hypox")
sampleTable12 <- data.frame(sampleName = sampleNames12, fileName = sampleFiles12, condition12 = sampleCondition12)
treatments12 = c("control","hypox")

# Load the data in Deseq format
ddsHTSeq12 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable12, directory = directory, design = ~ condition12)

colData(ddsHTSeq12)$condition12 <- factor(colData(ddsHTSeq12)$condition12,levels = treatments12)

#Perform the Differential Gene Expression using DESeq
dds12 <- DESeq(ddsHTSeq12)
res12 <- results(dds12)
res12
summary(res12)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 629, 3.3%
# LFC < 0 (down)     : 418, 2.2%

# # order results by padj value (most significant to least)
# res3_05 = subset(res3, padj<0.05)
# res3_05 <- res3_05[order(res3_05$padj),]
# res3_05
# summary(res3_05)
# # out of 10 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 40%
# # LFC < 0 (down)     : 6, 60%

# save data results and normalized reads to csv
resdata12 <- merge(as.data.frame(res12), as.data.frame(counts(dds12,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata12)[1] <- 'ensembl_gene_id'
write.csv(resdata12, "PN_vs_PH_n2_3-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PH_n2_3-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PH_n2_3-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PH_n2_3-results_merged.csv")


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
plotMA(dds12, ylim=c(-5,5),main = "Control VS Acute Hypoxia")
dev.copy(png, "PN_vs_PH_n2_3-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd12 <- varianceStabilizingTransformation(dds12, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition12 <- treatments12
pcaplot12 <- plotPCA(vsd12, intgroup=c("condition12")) + geom_text(aes(label=name), vjust=2)
pcaplot12
ggsave(pcaplot12,filename = "PN_vs_PH_n2_3_PCA_plot.pdf")

library("pheatmap")
mat12 = assay(vsd12)
#mat12 = assay(vsd12)[ head(order(res12$padj),30), ] # select the top 30 genes with the lowest padj
mat12 = mat12 - rowMeans(mat12) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df12 = as.data.frame(colData(vsd12)[,c("condition12")]) # Create a dataframe with a column of the conditions
colnames(df12) = "condition" # Rename the column header
rownames(df12) = colnames(mat12) # add rownames
write.table(mat12, "heatmap_data2.txt", sep="\t")
# and plot the actual heatmap
pheatmap(mat12, annotation_col=df12, main = "Top 30 Gene Expression")



# Heatmap
# select 35 genes with the highest variance across samples
topVarGenes12 <- head( order( rowVars( assay(vsd12) ), decreasing=TRUE ), 35 )
x12 <- assay(vsd12)[topVarGenes12,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd12)[topVarGenes12,],
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
dev.copy(png,"PN_vs_PH_Heatmap.png")
dev.off()

# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res12, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res12, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res12, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 



# PH vs PE (n3 vs n3)
sampleFiles13<- c("deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PH3_lane1_S18_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames13 <- c("PH1","PH2","PH3","PE2","PE3")
sampleCondition13 <- c("hypox","hypox","hypox","exp","exp")
sampleTable13 <- data.frame(sampleName = sampleNames13, fileName = sampleFiles13, condition13 = sampleCondition13)
treatments13 = c("hypox","exp")

# Load the data in Deseq format
ddsHTSeq13 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable13, directory = directory, design = ~ condition13)

colData(ddsHTSeq13)$condition13 <- factor(colData(ddsHTSeq13)$condition13,levels = treatments13)

#Perform the Differential Gene Expression using DESeq
dds13 <- DESeq(ddsHTSeq13)
res13 <- results(dds13)
res13
summary(res13)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 502, 2.7%
# LFC < 0 (down)     : 663, 3.6%

# # order results by padj value (most significant to least)
# res3_05 = subset(res3, padj<0.05)
# res3_05 <- res3_05[order(res3_05$padj),]
# res3_05
# summary(res3_05)
# # out of 10 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 40%
# # LFC < 0 (down)     : 6, 60%

# save data results and normalized reads to csv
resdata13 <- merge(as.data.frame(res13), as.data.frame(counts(dds13,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata13)[1] <- 'ensembl_gene_id'
write.csv(resdata13, "PH_vs_PE_n3_vs_n2-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PH_vs_PE_n3_vs_n2-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PH_vs_PE_n3_vs_n2-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description","sscrofa_homolog_associated_gene_name" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PH_vs_PE_n3_vs_n2-results_merged.csv")


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
plotMA(dds13, ylim=c(-5,5),main = "Hypoxia VS Experiment")
dev.copy(png, "PH_vs_PE_n3_vs_n2-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd13 <- varianceStabilizingTransformation(dds13, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition13 <- treatments13
pcaplot13 <- plotPCA(vsd13, intgroup=c("condition13")) + geom_text(aes(label=name), vjust=2)
pcaplot13
ggsave(pcaplot13,filename = "PH_vs_PE_n3_vs_n2_PCA_plot.pdf")

library("pheatmap")
mat13 = assay(vsd13)
#mat13 = assay(vsd13)[ head(order(res13$padj),30), ] # select the top 30 genes with the lowest padj
mat13 = mat13 - rowMeans(mat13) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df13 = as.data.frame(colData(vsd13)[,c("condition13")]) # Create a dataframe with a column of the conditions
colnames(df13) = "condition" # Rename the column header
rownames(df13) = colnames(mat13) # add rownames
write.table(mat13, "heatmap_data_PH_vs_PE_vst.txt", sep="\t")
# and plot the actual heatmap
pheatmap(mat13, annotation_col=df13, main = "Top 30 Gene Expression")



# Heatmap
# select 35 genes with the highest variance across samples
topVarGenes13 <- head( order( rowVars( assay(vsd13) ), decreasing=TRUE ), 35 )
x13 <- assay(vsd13)[topVarGenes13,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd13)[topVarGenes13,],
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
dev.copy(png,"PN_vs_PH_Heatmap.png")
dev.off()

# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res13, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res13, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res13, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 





# # Batch Correction
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BatchQC")
# library(BatchQC)
# metadf <- read.tsv("./metadata.tsv")
# exprdata <- read.table(file = 'PN_vs_PE-normalized_counts.txt', sep = '\t', header = TRUE)
# exprdata_with_no_0_values <- exprdata[apply(df[,-1], 1, function(x) !all(x==0)),]
# metadf <- read.table(file = 'PN_vs_PE_sampleinfo.txt', sep = '\t', header = TRUE)
# batch = metadf$sample
# condition = metadf$condition
# 
# batchQC(dat=exprdata_with_no_0_values, batch=batch, condition=condition,
#         report_file="batchqc_report.html", report_dir=".",
#         report_option_binary="111111111",
#         view_report=TRUE, interactive=TRUE, batchqc_output=TRUE)


# PN vs PE (n2_2)
sampleFiles14<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PE2_lane1_S14_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_PE3_2_lane1_S15_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames14 <- c("PN1","PN2","PE2","PE3")
sampleCondition14 <- c("normox","normox","exp","exp")
sampleTable14 <- data.frame(sampleName = sampleNames14, fileName = sampleFiles14, condition14 = sampleCondition14)
treatments14 = c("normox","exp")

# Load the data in Deseq format
ddsHTSeq14 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable14, directory = directory, design = ~ condition14)

colData(ddsHTSeq14)$condition14 <- factor(colData(ddsHTSeq14)$condition14,levels = treatments14)

#Perform the Differential Gene Expression using DESeq
dds14 <- DESeq(ddsHTSeq14)
res14 <- results(dds14)
res14
summary(res14)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 176, 0.95%
# LFC < 0 (down)     : 96, 0.52%

# # order results by padj value (most significant to least)
# res3_05 = subset(res3, padj<0.05)
# res3_05 <- res3_05[order(res3_05$padj),]
# res3_05
# summary(res3_05)
# # out of 10 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 4, 40%
# # LFC < 0 (down)     : 6, 60%

# save data results and normalized reads to csv
resdata14 <- merge(as.data.frame(res14), as.data.frame(counts(dds14,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata14)[1] <- 'ensembl_gene_id'
write.csv(resdata14, "PN_vs_PE_n2_2-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="sscrofa_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "PN_vs_PE_n2_2-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PN_vs_PE_n2_2-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "PN_vs_PE_n2_2-results_merged.csv")


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
plotMA(dds14, ylim=c(-5,5),main = "Normox VS Experiment")
dev.copy(png, "PN_vs_PE_n2_2-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd14 <- varianceStabilizingTransformation(dds14, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition14 <- treatments14
pcaplot14 <- plotPCA(vsd14, intgroup=c("condition14")) + geom_text(aes(label=name), vjust=2)
pcaplot14
ggsave(pcaplot14,filename = "PN_vs_PE_n2_2_PCA_plot.png")

library("pheatmap")
mat14 = assay(vsd14)
#mat14 = assay(vsd14)[ head(order(res14$padj),30), ] # select the top 30 genes with the lowest padj
mat14 = mat14 - rowMeans(mat14) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df14 = as.data.frame(colData(vsd14)[,c("condition14")]) # Create a dataframe with a column of the conditions
colnames(df14) = "condition" # Rename the column header
rownames(df14) = colnames(mat14) # add rownames
# and plot the actual heatmap
pheatmap(mat14, annotation_col=df14, main = "Top 30 Gene Expression")



# Heatmap
# select 35 genes with the highest variance across samples
topVarGenes13 <- head( order( rowVars( assay(vsd13) ), decreasing=TRUE ), 35 )
x13 <- assay(vsd13)[topVarGenes13,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd13)[topVarGenes13,],
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
dev.copy(png,"PN_vs_PH_Heatmap.png")
dev.off()

# Volcano Plot
# PN_vs_PE
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res14, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res14, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res14, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 
