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
library(knitr)
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Load biomaRt library
library("biomaRt")


deseq2_unannotated <- function(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments){
  sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
  # condition = sampleCondition
  # design = ~ condition
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = sample_directory, design = ~ condition)
  colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
  colData(ddsHTSeq)
  dds <- DESeq(ddsHTSeq)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  res <- results(dds)
  # order results by padj value (most significant to least)
  # res3 = subset(res3, padj<0.05)
  # res3 <- res3[order(res3$padj),]
  print(res)
  print(summary(res))
  write.table(summary(res),file = paste0(outputPrefix, "-results-summary.txt"), sep = "\t")
  # out of 18669 with nonzero total read count
  # adjusted p-value < 0.1
  # LFC > 0 (up)       : 49, 0.26%
  # LFC < 0 (down)     : 88, 0.47%
  # outliers [1]       : 0, 0%
  # low counts [2]     : 4706, 25%
  # How many adjusted p-values were less than 0.1?
  print(sum(res$padj < 0.1, na.rm=TRUE))
  write.table(sum(res$padj < 0.1, na.rm=TRUE),file = paste0(output_directory,outputPrefix, "-results-summary-p-0.1.txt"), sep = "\t")
  # [1] 137
  print(sum(res$padj < 0.05, na.rm=TRUE))
  write.table(sum(res$padj < 0.05, na.rm=TRUE),file = paste0(output_directory,outputPrefix, "-results-summary-p-0.05.txt"), sep = "\t")
  # [1] 58
  # # p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
  # res_05 <- results(dds, alpha=0.05)
  # summary(res_05)
  # # We can order our results table by the smallest p value
  # resOrdered3 <- res3[order(res3$pvalue),]
  # save data results and normalized reads to csv
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
  names(resdata)[1] <- 'ensembl_gene_id'
  write.csv(resdata, file = paste0(output_directory, outputPrefix, "-results-with-normalized.csv"))
  # List available biomaRt databases
  # Then select a database to use
  # listMarts()
  ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
  # List available datasets in database
  # Then select dataset to use
  # listDatasets(ensembl)
  # Check the database for entries that match the IDs of the differentially expressed genes from the results file
  # a2 <- read.table(resdata, head=TRUE)
  # listFilters()
  # listAttributes()
  b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=resdata$ensembl_gene_id, mart= ensembl)
  #colnames(a)[colnames(a)=="gene"] <- "ensembl_gene_id"
  m2 <- merge(resdata, b2, by="ensembl_gene_id")
  write.csv(as.data.frame(m2),file = paste0(output_directory,outputPrefix, "-normalized-Biomart.csv"))
  # MA plot of RNAseq data for entire dataset
  # http://en.wikipedia.org/wiki/MA_plot
  # The MA-plot provides a global view of the differential genes, 
  # with the log2 fold change on the y-axis over the mean of normalized counts
  # genes with padj < 0.1 are colored Red
  save_plot(paste0(output_directory,outputPrefix, "-MAplot.png"),plotMA(dds, ylim=c(-5,5), main = outputPrefix))
  # ggsave(ma_plot, file=paste0(output_directory,outputPrefix, "-MAplot.png"))
  # transform raw counts into normalized values
  # DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
  # variance stabilization is very good for heatmaps, etc.
  # Log-Transformation
  rld <- rlogTransformation(dds, blind=T)
  vsd <- varianceStabilizingTransformation(dds, blind=T)
  # PCA Plot
  condition <- treatments
  pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.1)
  pcaplot
  save_plot(paste0(output_directory,outputPrefix, "-PCA_plot.png"),pcaplot)
  # ggsave(pcaplot,file=paste0(output_directory,outputPrefix1, "-PCA_plot.png"))
  # Heatmap
  library("pheatmap")
  mat = assay(vsd)[ head(order(res$padj),100), ] # select the top 100 genes with the lowest padj
  # mat3 = mat3 - rowMeans(mat3) # Subtract the row means from each value
  # Optional, but to make the plot nicer:
  df = as.data.frame(colData(vsd)[,c("condition")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  # and plot the actual heatmap
  pheatmap <- pheatmap(mat, annotation_col=df, main = "Top 100 Gene Expression", fontsize_col = 10, fontsize_row = 4)
  pheatmap
  save_plot(paste0(output_directory,outputPrefix, "-Heatmap.png"),pheatmap)
  # ggsave(pheatmap1,file=paste0(output_directory, outputPrefix1, "-Heatmap.png"))
}

# Inputs to function
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/deseq2_counts_reverse"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/James/D2_wt_6_vs_D2_mdx_6_test/"
outputPrefix <- "D2_wt_6_vs_D2_mdx_6_test"
sampleFiles <- c("deseq2_counts_star_out_merged_TS6_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS9_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_XSKCN_20190719_K00134_IL100136246_TS25_L008_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS22_R1.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_counts_star_out_merged_TS23_R1.fastq.gzAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("TS6","TS9","TS25","TS22","TS23")
sampleCondition <- c("D2_wt_6","D2_wt_6","D2_wt_6","D2_mdx_6","D2_mdx_6")
treatments = c("D2_wt_6","D2_mdx_6")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)






