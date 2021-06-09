source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# Load DESeq2 library
library("DESeq2")

# Cell_Line1 
directory <- "/Users/pbanerjee/Desktop/Katrina_G245_Cell_line1_Feature_count_DEseq2"
setwd(directory)
# Set the prefix for each output file name
outputPrefix2 <- "Cell1_DESeq2_Fcounts"
sampleFiles2 <- c("NC1.txt",
                  "NC2.txt",
                  "NC3.txt",
                  "NC4.txt",
                  "NE1.txt",
                  "NE2.txt",
                  "NE3.txt",
                  "NE4.txt")
sampleNames2 <- c("NC1","NC2","NC3","NC4","NE1","NE2","NE3","NE4")
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
x <- sum(res2$padj < 0.1, na.rm=TRUE)
y <- sum(res2$padj < 0.05, na.rm=TRUE)
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
write.csv(as.data.frame(resSig2_05), file= paste0(outputPrefix2,"-results-with-padj.05.csv"))
# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")
#If there is error with rlang version, install rlang and then install packages and then it will automatically restart.
#install.packages("rlang")
# Convert final results .csv file into .txt file
results_csv2 <- "Cell1_DESeq2_Fcounts-results-with-normalized.csv"
write.table(read.csv(results_csv2), gsub(".csv",".txt",results_csv2))
results_txt2 <- "Cell1_DESeq2_Fcounts-results-with-normalized.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
# In cases when server is slow, its better to point to a specific server, then use ensembl2
ensembl2=useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2009.archive.ensembl.org", path = "/biomart/martservice", archive = FALSE)
#ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a2 <- read.table(results_txt2, head=TRUE)
# filters = listFilters(ensembl)
# head(filters)
# listFilters()
# attributes = listAttributes(ensembl)
# head(attributes)
# listAttributes()
#b2 <- getBM(filters = "ensembl_gene_id", attributes = c("mgi_symbol","ensembl_gene_id", "description"), values = a2$ensembl_gene_id, mart = ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
b2 <- getBM(filters = "ensembl_gene_id", attributes = c("mgi_symbol","ensembl_gene_id", "description"), values = a2$ensembl_gene_id, mart = ensembl2)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "Cell1_DESeq2_Fcounts-Biomart-results_merged.csv")


# Convert final p-adjusted results .csv file into .txt file
results_csv3 <- "Cell1_DESeq2_Fcounts-results-with-padj.1.csv"
write.table(read.csv(results_csv3), gsub(".csv",".txt",results_csv3))
results_txt3 <- "Cell1_DESeq2_Fcounts-results-with-padj.1.txt"
# List available biomaRt databases
# Then select a database to use
listMarts()
# In cases when server is slow, its better to point to a specific server, then use ensembl2
ensembl2=useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2009.archive.ensembl.org", path = "/biomart/martservice", archive = FALSE)
#ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl2)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a3 <- read.table(results_txt3, head=TRUE)
# filters = listFilters(ensembl)
# head(filters)
# listFilters()
# attributes = listAttributes(ensembl)
# head(attributes)
# listAttributes()
#b2 <- getBM(filters = "ensembl_gene_id", attributes = c("mgi_symbol","ensembl_gene_id", "description"), values = a2$ensembl_gene_id, mart = ensembl)
#colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
b3 <- getBM(filters = "ensembl_gene_id", attributes = c("mgi_symbol","ensembl_gene_id", "description"), values = a3$ensembl_gene_id, mart = ensembl2)
m3 <- merge(a3, b3, by="ensembl_gene_id")
write.csv(as.data.frame(m3),file = "Cell1_DESeq2_Fcounts-padj-0.1-Biomart-results_merged.csv")


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
#If RColorBrewer not insatlled first install, if already installed just load library
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer")
#   library(RColorBrewer)
# }
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# The MA-plot provides a global view of the differential genes, 
# with the log2 fold change on the y-axis over the mean of normalized counts
# genes with padj < 0.1 are colored Red

# Cell_Line2
plotMA(dds2, ylim=c(-5,5),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix2, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
#Cell_Line2
rld2 <- rlogTransformation(dds2, blind=T)
vsd2 <- varianceStabilizingTransformation(dds2, blind=T)
ntd <- normTransform(dds2)

# PCA Plot
# Cell_Line2
condition <- treatments2
pcaplot2 <- plotPCA(vsd2, intgroup=c("condition")) + geom_text(aes(label=name),vjust=2)
ggsave(pcaplot2,file=paste0(outputPrefix2, "-ggplot2.pdf"))

# Heatmap
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld2) ), decreasing=TRUE ), 35 )
jpeg(file = "heatmap3.jpeg")
heatmap <- heatmap.2( assay(rld2)[topVarGenes, ], hclustfun = hclust,scale="none",trace="none", dendrogram="both", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(5,20))
dev.off()

# Pheatmap
# install.packages("pheatmap", "RColorBrewer", "viridis")
library("pheatmap")
mat = assay(rld2)[ head(order(res2$padj),30), ] # select the top 30 genes with the lowest padj
mat2 = mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df2 = as.data.frame(colData(rld2)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df2) = "condition" # Rename the column header
rownames(df2) = colnames(mat2) # add rownames
# and plot the actual heatmap
pheatmap(mat2, annotation_col=df)


