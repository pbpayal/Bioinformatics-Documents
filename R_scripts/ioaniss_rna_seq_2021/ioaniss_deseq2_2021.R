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

# Set directory
directory <- "/Users/pbanerjee/Documents/Payal_Scripts/R/ioaniss_rna_seq_2021/"
setwd(directory)
# ALL Matrix together
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/Payal_Scripts/R/ioaniss_rna_seq_2021/counts_matrix.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/Payal_Scripts/R/ioaniss_rna_seq_2021/sample_annotation.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$ensemble_gene_id
count_matrix <- subset(count_matrix, select = -c(ensemble_gene_id))
sampleCondition <- c("sepsis_EVs","sepsis_EVs","sepsis_EVs",
                     "control_media","control_media","control_media",
                     "sepsis_media","sepsis_media", "sepsis_media",
                    "control_EV","control_EV","control_EV")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)
# VST
dds_vst <- vst(dds)
plotPCA(dds_vst)
pcaplot_vst<- plotPCA(dds_vst, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=0.1)
pcaplot_vst
ggsave(pcaplot_vst,filename = "PCA_plot_vst_with_labels.png")

# Differential Gene Expression
dds_deg <- DESeq(dds)
View(counts(dds_deg))
sizeFactors(dds_deg)

# Normalized Counts
normalized_counts <- counts(dds_deg, normalized=TRUE)
head(normalized_counts)
write.table(normalized_counts, file="/Users/pbanerjee/Documents/Payal_Scripts/R/ioaniss_rna_seq_2021/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

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

# Comparisons
deseq2_unannotated <- function(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments){
  sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
  # condition = sampleCondition
  # design = ~ condition
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = sample_directory, design = ~ condition)
  colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
  colData(ddsHTSeq)
  dds <- DESeq(ddsHTSeq)
  # keep <- rowSums(counts(dds)) >= 10
  # dds <- dds[keep,]
  res <- results(dds)
  # order results by padj value (most significant to least)
  # res3 = subset(res3, padj<0.05)
  # res3 <- res3[order(res3$padj),]
  print(res)
  print(summary(res))
  write.table(summary(res),file = paste0(outputPrefix, "-results-summary.txt"), sep = "\t")
  # How many adjusted p-values were less than 0.1?
  print(sum(res$padj < 0.1, na.rm=TRUE))
  write.table(sum(res$padj < 0.1, na.rm=TRUE),file = paste0(output_directory,outputPrefix, "-results-summary-p-0.1.txt"), sep = "\t")
  print(sum(res$padj < 0.05, na.rm=TRUE))
  write.table(sum(res$padj < 0.05, na.rm=TRUE),file = paste0(output_directory,outputPrefix, "-results-summary-p-0.05.txt"), sep = "\t")
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
  pcaplot <- plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name, vjust = 0.005), vjust=0.05)
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

# SMvsCE
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
# output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/SMvsCE/"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/"
outputPrefix <- "SMvsCE"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-12_S7_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-13_S8_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-14_S9_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                "deseq2_fcounts_star_out_final_merged_G738-15_S10_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-16_S11_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-17_S12_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("SM1","SM2","SM3","CE1","CE2","CE3")
sampleCondition <- c("sepsis_media","sepsis_media","sepsis_media","control_EV","control_EV","control_EV")
treatments = c("sepsis_media","control_EV")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)

# CMvsSM
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/CMvsSM/"
outputPrefix <- "CMvsSM"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-06_S4_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-07_S5_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-10_S6_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-12_S7_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-13_S8_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-14_S9_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("CM1","CM2","CM3","SM1","SM2","SM3")
sampleCondition <- c("control_media","control_media","control_media","sepsis_media","sepsis_media","sepsis_media")
treatments = c("control_media","sepsis_media")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)

# CMvsSE
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/CMvsSE/"
outputPrefix <- "CMvsSE"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-06_S4_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-07_S5_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-10_S6_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-01_S1_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-02_S2_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-03_S3_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("CM1","CM2","CM3","SE1","SE2","SE3")
sampleCondition <- c("control_media","control_media","control_media","sepsis_EVs","sepsis_EVs","sepsis_EVs")
treatments = c("control_media","sepsis_EVs")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)


# CEvsSE
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/CEvsSE/"
outputPrefix <- "CEvsSE"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-15_S10_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-16_S11_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-17_S12_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-01_S1_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-02_S2_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-03_S3_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("CE1","CE2","CE3","SE1","SE2","SE3")
sampleCondition <- c("control_EV","control_EV","control_EV","sepsis_EVs","sepsis_EVs","sepsis_EVs")
treatments = c("control_EV","sepsis_EVs")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)


# CMvsCE
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/CMvsCE/"
outputPrefix <- "CMvsCE"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-06_S4_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-07_S5_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-10_S6_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-15_S10_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-16_S11_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-17_S12_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("CM1","CM2","CM3","CE1","CE2","CE3")
sampleCondition <- c("control_media","control_media","control_media","control_EV","control_EV","control_EV")
treatments = c("control_media","control_EV")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)


# SMvsSE
sample_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/Deseq2_counts_files"
output_directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/2021/Genmed/Ioannis/SMvsSE/"
outputPrefix <- "SMvsSE"
sampleFiles <- c("deseq2_fcounts_star_out_final_merged_G738-12_S7_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-13_S8_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-14_S9_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-01_S1_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-02_S2_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_final_merged_G738-03_S3_L001_R1_001.fastqAligned.sortedByCoord.out.bam.txt")
sampleNames <- c("SM1","SM2","SM3","SE1","SE2","SE3")
sampleCondition <- c("sepsis_media","sepsis_media","sepsis_media","sepsis_EVs","sepsis_EVs","sepsis_EVs")
treatments = c("sepsis_media","sepsis_EVs")
# Run function
deseq2_unannotated(sample_directory, output_directory, outputPrefix, sampleFiles, sampleNames,sampleCondition, treatments)


