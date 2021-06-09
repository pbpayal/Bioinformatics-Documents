setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/")
save.image(file='Panos-mouse.RData')
load('Panos-mouse.RData')


# If Deseq2 is not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
# Load DESeq2 library, if Deseq2 is already installed
library("DESeq2")

#Please check before installing if you have installed devtools, if not then only do the following otherwise just load the library
#install.packages("devtools")
library(devtools)

setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/")
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Counts_matrix.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/sample_annotation.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$GeneID
count_matrix <- subset(count_matrix, select = -c(GeneID) )
sampleCondition <- c("control","control","control",
                     "hpx_lps","hpx_lps","hpx_lps",
                     "hyx","hyx","hpx",
                     "lps","lps","lps")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)


# Differential Gene Expression
dds_mouse <- DESeq(dds)
res <- results(dds_mouse)
res
# log2 fold change (MLE): condition lps vs control 
# Wald test p-value: condition lps vs control 
# DataFrame with 53463 rows and 6 columns
summary(res)
# out of 31117 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 277, 0.89%
# LFC < 0 (down)     : 97, 0.31%
# outliers [1]       : 71, 0.23%
# low counts [2]     : 13930, 45%

dds_mouse_vst <- vst(dds_mouse)
mat_mouse_all = assay(dds_mouse_vst)[ head(order(res$padj)), ]
mat_mouse_all = assay(dds_mouse_vst)

mat_mouse_all = mat_mouse_all - rowMeans(mat_mouse_all)
library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a <- as.data.frame(mat_mouse_all)
a$Ensembl_Id <- rownames(mat_mouse_all)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description" ), values=a$Ensembl_Id, mart= ensembl)
colnames(a)[colnames(a)=="Ensembl_Id"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "Heatmap_mouse_data_all_genes_vst.csv")

# Final heatmap
heatmap_data = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/Heatmap_mouse_data_all_genes_vst_1.csv", sep = ",", header = TRUE, stringsAsFactors=FALSE)
heatmap_data <- heatmap_data[c(1,6,10,14,18)]
# row.names(heatmap_data) <- heatmap_data$external_gene_name
# C730027H18Rik’, ‘Gm18645’ -2, ‘Gm35558’ - 2, ‘Gm38642’ -3, ‘Gm5089’ - 2 , ‘Jakmip1’ - 2, ‘Ndor1’ -3, ‘Zkscan7’-2
heatmap_data_no_duplicates <- heatmap_data[!duplicated(heatmap_data$external_gene_name),]
# Removed 9 genes
row.names(heatmap_data_no_duplicates) <- heatmap_data_no_duplicates$external_gene_name
# heatmap_data_subset has 9 duplicate genes less than heatmap_data 
heatmap_data_subset <- subset(heatmap_data_no_duplicates, select = -c(external_gene_name))
heatmap_data_subset <- as.matrix(heatmap_data_subset)

# rename
colnames(heatmap_data_subset)[1] <- "Normoxic"
colnames(heatmap_data_subset)[2]<- "Hypoxic-LPS"
colnames(heatmap_data_subset)[3]<- "Hypoxic"
colnames(heatmap_data_subset)[4]<- "LPS"
# reorder
heatmap_reorder <- heatmap_data_subset[, c("Normoxic","Hypoxic","LPS","Hypoxic-LPS")]
write.table(heatmap_reorder, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_reorder_mouse.txt", sep = "\t")
#heatmap
pheatmap(heatmap_reorder, main = "Gene Expression", 
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-0.5, +0.5, length = 101))

# Removed rows where the gene expression is 0 in all columns or samples in this case
#df2 <- heatmap_reorder[!rowSums(heatmap_reorder[,]) == 0, ] # its deleting rows where negative and positive values equates to 0 too. Won't work here.
# df3 <- heatmap_reorder[!rowSums(heatmap_reorder[,] == 0) == (ncol(heatmap_reorder)-1), ]
# non_zero_expression_genes <- heatmap_reorder[as.logical(rowSums(heatmap_reorder != 0)), ] # taking care of both negative and positive values
#write.table(non_zero_expression_genes, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_no_zero_mouse.txt", sep = "\t")
non_zero_expression_genes_apply_mehod <- heatmap_reorder[!apply(heatmap_reorder, 1, function(x) all(x == 0)), ]
#write.table(non_zero_expression_genes_apply_mehod, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_no_zero_apply_mouse.txt", sep = "\t")

mouse_all_genes_nopheatmap_non_zero_expression_genes <- pheatmap(non_zero_expression_genes_apply_mehod, main = "Gene Expression", 
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-0.5,+0.5, length = 101))
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/Mouse_heatmap_all_genes_with_no_zero_expression1.png", mouse_all_genes_nopheatmap_non_zero_expression_genes)
# mymat <- matrix()
# for(i in 1:dim(heatmap_reorder)[1]) {
#   for(j in 1:dim(heatmap_reorder)[2]) {
#     if(heatmap_reorder[i,j] == 0)
#     break()
#     else(mymat[i,j])
#   }
# }
# mat <- heatmap_reorder
# for(row in 1:nrow(mat)) {
#   for(col in 1:ncol(my_matrix)) {
#     print(my_matrix[row, col])
#   }
# }


# pheatmap(heatmap_reorder[2:600,], main = "Top 600 Gene Expression",
#          show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
#          breaks = seq(-2, +2, length = 101))
# pheatmap(heatmap_data[1:100,], main = "Top 100 Gene Expression", 
#          show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T)
head(brewer.pal.info)
table(brewer.pal.info$category)
brewer.pal(11, "RdYlGn")

head(heatmap_data[, c("Normal","Hypoxia","Experiment")])

# Sub Heatmaps DHvsH
sig_gene_data = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/DHvsH_sig_genes_heatmap_data.txt", sep = "\t", header = TRUE)
# sig_gene_data_sub <- sig_gene_data[c(2,3,4)]
row.names(sig_gene_data) <- sig_gene_data$DHvsH_significant_genes
sig_gene_data <- subset(sig_gene_data, select = -c(DHvsH_significant_genes))
colnames(sig_gene_data)[1] <- "Hypoxic-LPS"
colnames(sig_gene_data)[2] <- "Hypoxic"
sig_gene_data <- as.matrix(sig_gene_data)
heatmap_reorder <- sig_gene_data[, c("Hypoxic","Hypoxic-LPS")]
non_zero_expression_genes_apply_mehod_sub <- heatmap_reorder[!apply(heatmap_reorder, 1, function(x) all(x == 0)), ]
sub_heatmap_DHvsH <- pheatmap(non_zero_expression_genes_apply_mehod_sub, main = "Hypoxic vs Hypoxic-LPS",
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-0.7, +0.7, length = 101))
sub_heatmap_DHvsH
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/HvsDH_subheatmap.png",sub_heatmap_DHvsH)


# Sub Heatmaps CvsL
sig_gene_data_CvsL = read.table("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/CvsL_heatmap_data.txt", sep = "\t", header = TRUE)
# sig_gene_data_sub <- sig_gene_data[c(2,3,4)]
row.names(sig_gene_data_CvsL) <- sig_gene_data_CvsL$significant_gene_names
sig_gene_data_CvsL <- subset(sig_gene_data_CvsL, select = -c(significant_gene_names))
colnames(sig_gene_data_CvsL)[1] <- "Control"
colnames(sig_gene_data_CvsL)[2] <- "LPS"
sig_gene_data_CvsL <- as.matrix(sig_gene_data_CvsL)
#heatmap_reorder_CvsL <- sig_gene_data_CvsL[, c("Hypoxic","Hypoxix-LPS")]
non_zero_expression_genes_apply_mehod_sub_CvsL <- sig_gene_data_CvsL[!apply(sig_gene_data_CvsL, 1, function(x) all(x == 0)), ]
sub_heatmap_CvsL <- pheatmap(non_zero_expression_genes_apply_mehod_sub_CvsL, main = " Control vs LPS",
                              show_rownames=TRUE, cluster_rows=T, cluster_cols=F, show_colnames = T, fontsize_row=6,
                              angle_col = 0,
                              breaks = seq(-2, +2, length = 101))
sub_heatmap_CvsL
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Heatmaps/CvsL_subheatmap.png",sub_heatmap_CvsL)


dds_rlog <- rlog(dds_mouse)
pcaplot <- plotPCA(dds_rlog, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=3)
pcaplot
ggsave(pcaplot,filename = "PCA_plot.png")

# https://www.datacamp.com/community/tutorials/pca-analysis-r
counts.pca <- prcomp(count_matrix[,c(1:12)], center = TRUE,scale. = TRUE)
summary(counts.pca)
str(counts.pca)

install.packages('ggfortify')
library("ggfortify")
autoplot(counts.pca)


autoplot(counts.pca, colour = c("DH2","DN1"))

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)


# Mouse

# set diretory to where your counts files are
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/Deseq2_files/"
setwd(directory)

# C vs DH
sampleFiles1<- c("deseq2_fcounts_star_out_C1_lane1_S1_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_C2_lane1_S2_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_C3_lane1_S3_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_DH1_lane1_S4_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_DH2_lane1_S5_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_DH3_lane1_S6_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames1 <- c("C1","C2","C3","DH1","DH2","DH3")
sampleCondition1 <- c("control","control","control","DH","DH","DH")
sampleTable1 <- data.frame(sampleName = sampleNames1, fileName = sampleFiles1, condition1 = sampleCondition1)
treatments1 = c("control","DH")

# Load the data in Deseq format
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable1, directory = directory, design = ~ condition1)

colData(ddsHTSeq1)$condition1 <- factor(colData(ddsHTSeq1)$condition1,levels = treatments1)

#Perform the Differential Gene Expression using DESeq
dds1 <- DESeq(ddsHTSeq1)
res1 <- results(dds1)
res1
summary(res1)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 12
# LFC < 0 (down)     : 4

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
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata1)[1] <- 'ensembl_gene_id'
write.csv(resdata1, "C_vs_DH-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "C_vs_DH-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "C_vs_DH-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "C_vs_DH-results_merged.csv")


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
plotMA(dds1, ylim=c(-5,5),main = "C VS DH")
dev.copy(png, "C_vs_DH-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd1 <- varianceStabilizingTransformation(dds1, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition1 <- treatments1
pcaplot1 <- plotPCA(vsd1, intgroup=c("condition1")) + geom_text(aes(label=name), vjust=2)
pcaplot1
ggsave(pcaplot1,filename = "C_vs_DH_PCA_plot.pdf")

library("pheatmap")
mat1 = assay(vsd1)[ head(order(res1$padj),30), ] # select the top 30 genes with the lowest padj
mat1 = mat1 - rowMeans(mat1) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df1 = as.data.frame(colData(vsd1)[,c("condition1")]) # Create a dataframe with a column of the conditions
colnames(df1) = "condition" # Rename the column header
rownames(df1) = colnames(mat1) # add rownames
# and plot the actual heatmap
pheatmap(mat1, annotation_col=df1, main = "Top 30 Gene Expression")



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
with(res1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res1, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res1, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 

# C vs HM
sampleFiles2<- c("deseq2_fcounts_star_out_C1_lane1_S1_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_C2_lane1_S2_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_C3_lane1_S3_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H1M_lane1_S7_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H2M_lane1_S8_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H3M_lane1_S9_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames2 <- c("C1","C2","C3","H1M","H2M","H3M")
sampleCondition2 <- c("C","C","C","HM","HM","HM")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition2 = sampleCondition2)
treatments2 = c("C","HM")

# Load the data in Deseq format
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory, design = ~ condition2)

colData(ddsHTSeq2)$condition2 <- factor(colData(ddsHTSeq2)$condition2,levels = treatments2)

#Perform the Differential Gene Expression using DESeq
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
res2
summary(res2)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3
# LFC < 0 (down)     : 9

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
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'ensembl_gene_id'
write.csv(resdata2, "C_vs_HM-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "C_vs_HM-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "C_vs_HM-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "C_vs_HM-results_merged.csv")


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
plotMA(dds2, ylim=c(-5,5),main = "C VS HM")
dev.copy(png, "C_vs_HM-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition2 <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition2")) + geom_text(aes(label=name), vjust=2)
pcaplot2
ggsave(pcaplot2,filename = "C_vs_HM_PCA_plot.pdf")

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
with(res2, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res2, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res2, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# C vs L
sampleFiles3<- c("deseq2_fcounts_star_out_C1_lane1_S1_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_C2_lane1_S2_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_C3_lane1_S3_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_L1_lane1_S10_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_L2_lane1_S11_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_L3_lane1_S12_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames3 <- c("C1","C2","C3","L1","L2","L3")
sampleCondition3 <- c("C","C","C","L","L","L")
sampleTable3 <- data.frame(sampleName = sampleNames3, fileName = sampleFiles3, condition3 = sampleCondition3)
treatments3 = c("C","L")

# Load the data in Deseq format
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable3, directory = directory, design = ~ condition3)

colData(ddsHTSeq3)$condition3 <- factor(colData(ddsHTSeq3)$condition3,levels = treatments3)

#Perform the Differential Gene Expression using DESeq
dds3 <- DESeq(ddsHTSeq3)
res3 <- results(dds3)
res3
summary(res3)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 59
# LFC < 0 (down)     : 6

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
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata3)[1] <- 'ensembl_gene_id'
write.csv(resdata3, "C_vs_L-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "C_vs_L-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "C_vs_L-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "C_vs_L-results_merged.csv")


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
plotMA(dds3, ylim=c(-5,5),main = "C VS L")
dev.copy(png, "C_vs_L-MAplot.png")
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
ggsave(pcaplot3,filename = "C_vs_L_PCA_plot.pdf")

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
with(res3, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res3, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res3, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 



# DH vs HM
sampleFiles4 <- c("deseq2_fcounts_star_out_DH1_lane1_S4_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_DH2_lane1_S5_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_DH3_lane1_S6_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H1M_lane1_S7_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H2M_lane1_S8_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_H3M_lane1_S9_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames4 <- c("DH1","DH2","DH3","H1M","H2M","H3M")
sampleCondition4 <- c("Double-Hit","Double-Hit","Double-Hit",
                      "Hypoxia-only","Hypoxia-only","Hypoxia-only")
sampleTable4 <- data.frame(sampleName = sampleNames4, fileName = sampleFiles4, condition4 = sampleCondition4)
treatments4 = c("Double-Hit","Hypoxia-only")

# Load the data in Deseq format
ddsHTSeq4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable4, 
                                        directory = directory, 
                                        design = ~ condition4)

colData(ddsHTSeq4)$condition4 <- factor(colData(ddsHTSeq4)$condition4,
                                        levels = treatments4)

#Perform the Differential Gene Expression using DESeq
dds4 <- DESeq(ddsHTSeq4)
res4 <- results(dds4)
res4
summary(res4)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1259, 4.4%
# LFC < 0 (down)     : 1974, 6.8%

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
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata4)[1] <- 'ensembl_gene_id'
write.csv(resdata4, "DH_vs_HM-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "DH_vs_HM-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "DH_vs_HM-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description"), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "DH_vs_HM-DEG.csv")


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
plotMA(dds4, ylim=c(-5,5),main = "C VS L")
dev.copy(png, "DH_vs_HM-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd4 <- varianceStabilizingTransformation(dds4, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition4 <- treatments4
pcaplot4 <- plotPCA(vsd4, intgroup=c("condition4")) + geom_text(aes(label=name), vjust=2)
pcaplot4
ggsave(pcaplot4,filename = "DH_vs_HM_PCA_plot.png")

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
with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res4, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# DH vs L
sampleFiles5 <- c("deseq2_fcounts_star_out_DH1_lane1_S4_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_DH2_lane1_S5_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_DH3_lane1_S6_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L1_lane1_S10_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L2_lane1_S11_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L3_lane1_S12_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames5 <- c("DH1","DH2","DH3","L1","L2","L3")
sampleCondition5 <- c("Double_Hit","Double_Hit","Double_Hit",
                      "LPS_only","LPS_only","LPS_only")
sampleTable5 <- data.frame(sampleName = sampleNames5, 
                           fileName = sampleFiles5, 
                           condition5 = sampleCondition5)
treatments5 = c("Double_Hit","LPS_only")

# Load the data in Deseq format
ddsHTSeq5 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable5, 
                                        directory = directory, 
                                        design = ~ condition5)

colData(ddsHTSeq5)$condition5 <- factor(colData(ddsHTSeq5)$condition5,
                                        levels = treatments5)

#Perform the Differential Gene Expression using DESeq
dds5 <- DESeq(ddsHTSeq5)
res5 <- results(dds5)
res5
summary(res5)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1259, 4.4%
# LFC < 0 (down)     : 1974, 6.8%

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
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata5)[1] <- 'ensembl_gene_id'
write.csv(resdata5, "DH_vs_L-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "DH_vs_L-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "DH_vs_L-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "DH_vs_L-DEG.csv")


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
plotMA(dds5, ylim=c(-5,5),main = "C VS L")
dev.copy(png, "DH_vs_L-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd5 <- varianceStabilizingTransformation(dds5, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition5 <- treatments5
pcaplot5 <- plotPCA(vsd5, intgroup=c("condition5")) + geom_text(aes(label=name), vjust=2)
pcaplot5
ggsave(pcaplot5,filename = "DH_vs_L_PCA_plot.png")

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
with(res5, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res5, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res5, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


# HM vs L
sampleFiles6 <- c("deseq2_fcounts_star_out_H1M_lane1_S7_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_H2M_lane1_S8_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_H3M_lane1_S9_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L1_lane1_S10_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L2_lane1_S11_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                  "deseq2_fcounts_star_out_L3_lane1_S12_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames6 <- c("H1M","H2M","H3M","L1","L2","L3")
sampleCondition6 <- c("Hypoxia_only","Hypoxia_only","Hypoxia_only",
                      "LPS_only","LPS_only","LPS_only")
sampleTable6 <- data.frame(sampleName = sampleNames6, 
                           fileName = sampleFiles6, 
                           condition6 = sampleCondition6)
treatments6 = c("Hypoxia_only","LPS_only")

# Load the data in Deseq format
ddsHTSeq6 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable6, 
                                        directory = directory, 
                                        design = ~ condition6)

colData(ddsHTSeq6)$condition6 <- factor(colData(ddsHTSeq6)$condition6,
                                        levels = treatments6)

#Perform the Differential Gene Expression using DESeq
dds6 <- DESeq(ddsHTSeq6)
res6 <- results(dds6)
res6
summary(res6)
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2684, 9.3%
# LFC < 0 (down)     : 2029, 7%

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
resdata6 <- merge(as.data.frame(res6), as.data.frame(counts(dds6,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata6)[1] <- 'ensembl_gene_id'
write.csv(resdata6, "HM_vs_L-results.csv")

library("biomaRt")

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")


# Convert final results .csv file into .txt file
results_csv <- "HM_vs_L-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "HM_vs_L-results.txt"
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)
# listFilters()
# attributes <- listAttributes(ensembl)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description" ), values=a$ensembl_gene_id, mart= ensembl)
#colnames(a)[colnames(a)=="geneid"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "HM_vs_L-DEG.csv")


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
plotMA(dds6, ylim=c(-5,5),main = "C VS L")
dev.copy(png, "HM_vs_L-MAplot.png")
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
# PN_vs_PE
#rld3 <- rlogTransformation(dds3, blind=T)
vsd6 <- varianceStabilizingTransformation(dds6, blind=T)

# PCA Plot
# Shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the 
# overall effect of experimental covariates and batch effects.
condition6 <- treatments6
pcaplot6 <- plotPCA(vsd6, intgroup=c("condition6")) + geom_text(aes(label=name), vjust=2)
pcaplot6
ggsave(pcaplot6,filename = "HM_vs_L_PCA_plot.png")

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
with(res6, plot(log2FoldChange, -log10(pvalue), 
                pch=20, main="Volcano plot", 
                xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res6, padj<.01 ), 
     points(log2FoldChange, -log10(pvalue), 
            pch=20, col="blue"))
with(subset(res6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# Save the diagram directly by exporting it from R 


