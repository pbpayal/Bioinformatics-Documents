setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/")
save.image(file='/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Panos-human2.RData')
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
library(biomaRt)
library(PCAtools)
library(cowplot)
library(pheatmap)

# All HC vs All HP
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/counts_matrix.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/sample_annotation.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$Gene_ID
count_matrix <- subset(count_matrix, select = -c(Gene_ID) )
sampleCondition <- c("control","control","control","control",
                     "control","control","control",
                     "exp", "exp", "exp","exp","exp","exp",
                     "exp","exp","exp","exp","exp")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)
dds$condition
dds$sample
# Differential Gene Expression
dds_human <- DESeq(dds)
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 1257 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
res <- results(dds_human)
res
# log2 fold change (MLE): condition exp vs control 
# Wald test p-value: condition exp vs control 
summary(res)
# out of 40279 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 0.0025%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

res_contrasts <- results(dds_human, contrast=c("condition","control","exp"))
res_contrasts
summary(res_contrasts)

# PCA Plot
dds_rlog <- rlog(dds_human)
pcaplot_human <- plotPCA(dds_rlog, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=3)
pcaplot_human
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/pcaplot_human.png", pcaplot_human)

# save data results and normalized reads to csv
res_with_normalized_val <- merge(as.data.frame(res), as.data.frame(counts(dds_human,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(res_with_normalized_val)[1] <- 'ensembl_gene_id'
write.csv(res_with_normalized_val, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/res_with_normalized_val.csv")
# Annotate the genes
# Convert final results .csv file into .txt file
results_csv <- "res_with_normalized_val.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "res_with_normalized_val.txt"
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a1 <- read.table(results_txt, head=TRUE)
b1 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description" ), values=a1$ensembl_gene_id, mart= ensembl)
# colnames(a1)[colnames(a1)=="Ensembl_Id"] <- "ensembl_gene_id"
m1 <- merge(a1, b1, by="ensembl_gene_id")
write.csv(as.data.frame(m1),file = "human_results.csv")

# Heatmap all genes
dds_human_vst <- vst(dds_human)
mat_human_all = assay(dds_human_vst)
mat_human_all = mat_human_all - rowMeans(mat_human_all)
# Annotate the genes
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
a <- as.data.frame(mat_human_all)
a$Ensembl_Id <- rownames(mat_human_all)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description" ), values=a$Ensembl_Id, mart= ensembl)
colnames(a)[colnames(a)=="Ensembl_Id"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "Heatmap_human_data_all_genes_vst.csv")
heatmap_data = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/Heatmap_human_data_all_genes_vst_1.csv", sep = ",", header = TRUE, stringsAsFactors=FALSE)
# heatmap_data_SUID <- heatmap_data_SUID[c(1,10,21)]
# Remove duplicate gene names
heatmap_data_no_duplicates <- heatmap_data[!duplicated(heatmap_data$external_gene_name),]
row.names(heatmap_data_no_duplicates) <- heatmap_data_no_duplicates$external_gene_name
# 52,063 genes from 52,986 genes 
heatmap_data_subset <- subset(heatmap_data_no_duplicates, select = -c(external_gene_name))
heatmap_data_subset <- as.matrix(heatmap_data_subset)
# rename
colnames(heatmap_data_subset)[1] <- "AVG-HC"
colnames(heatmap_data_subset)[2]<- "AVG-HP"
# # reorder
# heatmap_reorder <- heatmap_data_subset[, c("Normoxic","Hypoxic","LPS","Hypoxic-LPS")]
# write.table(heatmap_reorder, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_reorder_mouse.txt", sep = "\t")
# remove 0 expression from dataset, 39,705  from 52,063
non_zero_expression_genes_apply_mehod <- heatmap_data_subset[!apply(heatmap_data_subset, 1, function(x) all(x == 0)), ]
#heatmap
heatmap_all_genes_with_nozero_exp_human <- pheatmap(non_zero_expression_genes_apply_mehod, main = "HC vs HP", 
                                                         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
                                                         angle_col = 0,
                                                         breaks = seq(-0.1, +0.1, length = 101))
heatmap_all_genes_with_nozero_exp_human
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/heatmap_all_genes_with_nozero_exp_human.png",heatmap_all_genes_with_nozero_exp_human)


# non-SUID vs SUID
count_matrix_SUID <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/SUID/counts_matrix_SUID_1.csv', header =T, sep = ",")
sample_annotation_SUID <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/SUID/sample_annotation_SUID.csv", header =T, sep = ",")
row.names(count_matrix_SUID) <- count_matrix_SUID$Gene_ID
count_matrix_SUID <- subset(count_matrix_SUID, select = -c(Gene_ID) )
sampleCondition_SUID <- c("non_SUID","non_SUID","non_SUID","non_SUID","non_SUID","non_SUID",
                     "SUID", "SUID", "SUID","SUID","SUID","SUID",
                     "SUID","SUID","SUID","SUID","SUID","SUID")
condition = sampleCondition_SUID
dds_SUID <- DESeqDataSetFromMatrix(countData = count_matrix_SUID,
                              colData = sample_annotation_SUID, 
                              design = ~ condition)


dds_SUID$condition
dds_SUID$sample
# Differential Gene Expression
dds_human_SUID <- DESeq(dds_SUID)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 1189 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
res_SUID <- results(dds_human_SUID)
res_SUID
# log2 fold change (MLE): condition SUID vs non SUID 
# Wald test p-value: condition SUID vs non SUID 
summary(res_SUID)
# out of 39889 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 996, 2.5%
# LFC < 0 (down)     : 1305, 3.3%
# outliers [1]       : 71, 0.18%
# low counts [2]     : 19914, 50%
## When I include sample HP4 and HP5
# out of 40037 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 791, 2%
# LFC < 0 (down)     : 1024, 2.6%
# outliers [1]       : 85, 0.21%
# low counts [2]     : 20977, 52%


res_SUID_contrasts <- results(dds_human_SUID, contrast=c("condition","SUID","non_SUID"))
res_SUID_contrasts
# log2 fold change (MLE): condition SUID vs non_SUID 
# Wald test p-value: condition SUID vs non SUID 
summary(res_SUID_contrasts)
# out of 39889 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 996, 2.5%
# LFC < 0 (down)     : 1305, 3.3%

# save data results and normalized reads to csv
res_SUID_with_normalized_val <- merge(as.data.frame(res_SUID), as.data.frame(counts(dds_human_SUID,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(res_SUID_with_normalized_val)[1] <- 'ensembl_gene_id'
write.csv(res_SUID_with_normalized_val, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/res_SUID_with_normalized_val.csv")
# Annotate the genes
# Convert final results .csv file into .txt file
results_csv <- "res_SUID_with_normalized_val.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "res_SUID_with_normalized_val.txt"
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a1 <- read.table(results_txt, head=TRUE)
b1 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description" ), values=a1$ensembl_gene_id, mart= ensembl)
# colnames(a1)[colnames(a1)=="Ensembl_Id"] <- "ensembl_gene_id"
m1 <- merge(a1, b1, by="ensembl_gene_id")
write.csv(as.data.frame(m1),file = "human_SUID_results.csv")


# # order results by padj value (most significant to least)
res_SUID_05 = subset(res_SUID, padj<0.05)
res_SUID_05 <- res_SUID_05[order(res_SUID_05$padj),]
res_SUID_05
# log2 fold change (MLE): condition SUID vs non SUID 
# Wald test p-value: condition SUID vs non SUID 
# DataFrame with 1379 rows and 6 columns
summary(res_SUID_05)
# out of 1379 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 516, 37%
# LFC < 0 (down)     : 863, 63%



# PCA Plot
dds_rlog_SUID <- rlog(dds_human_SUID)
pcaplot_SUID <- plotPCA(dds_rlog_SUID, intgroup=c("condition")) + geom_text(aes(label=name), vjust=1) + geom_point(size=3)
pcaplot_SUID
save_plot()

# Heatmap all genes
dds_human_SUID_vst <- vst(dds_human_SUID)
mat_human_SUID_all = assay(dds_human_SUID_vst)
mat_human_SUID_all = mat_human_SUID_all - rowMeans(mat_human_SUID_all)
# Annotate the genes
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
a <- as.data.frame(mat_human_SUID_all)
a$Ensembl_Id <- rownames(mat_human_SUID_all)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description" ), values=a$Ensembl_Id, mart= ensembl)
colnames(a)[colnames(a)=="Ensembl_Id"] <- "ensembl_gene_id"
m <- merge(a, b, by="ensembl_gene_id")
write.csv(as.data.frame(m),file = "Heatmap_human_SUID_data_all_genes_vst.csv")
heatmap_data_SUID = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/Heatmap_human_SUID_data_all_genes_vst_1.csv", sep = ",", header = TRUE, stringsAsFactors=FALSE)
heatmap_data_SUID <- heatmap_data_SUID[c(1,10,21)]
# Remove duplicate gene names
heatmap_data_SUID_no_duplicates <- heatmap_data_SUID[!duplicated(heatmap_data_SUID$external_gene_name),]
row.names(heatmap_data_SUID_no_duplicates) <- heatmap_data_SUID_no_duplicates$external_gene_name
# 52,063 genes from 52,986 genes 
heatmap_data_SUID_subset <- subset(heatmap_data_SUID_no_duplicates, select = -c(external_gene_name))
heatmap_data_SUID_subset <- as.matrix(heatmap_data_SUID_subset)
# rename
colnames(heatmap_data_SUID_subset)[1] <- "Non-SUID"
colnames(heatmap_data_SUID_subset)[2]<- "SUID"
# # reorder
# heatmap_reorder <- heatmap_data_subset[, c("Normoxic","Hypoxic","LPS","Hypoxic-LPS")]
# write.table(heatmap_reorder, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_reorder_mouse.txt", sep = "\t")
# remove 0 expression from dataset, 39,334 from 52,063
non_zero_expression_genes_apply_mehod_SUID <- heatmap_data_SUID_subset[!apply(heatmap_data_SUID_subset, 1, function(x) all(x == 0)), ]
#heatmap
heatmap_all_genes_with_nozero_exp_human_SUID <- pheatmap(non_zero_expression_genes_apply_mehod_SUID, main = "Non-SUID vs SUID", 
         show_rownames=FALSE, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-0.2, +0.2, length = 101))
heatmap_all_genes_with_nozero_exp_human_SUID
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/heatmap_all_genes_with_nozero_exp_human_SUID.png",heatmap_all_genes_with_nozero_exp_human_SUID)

# Heatmap for genes with p-adj < 0.05
# filtered for gene names from the original Deseq2 results files wih padj < 0.05. 
# then mapped the genes for vst values from Heatmap_human_SUID_data_all_genes_vst_1.csv file
heatmap_data_SUID_padj = read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/Heatmap_human_SUID_padj_vst_1.csv", sep = ",", header = TRUE, stringsAsFactors=FALSE)
# heatmap_data_SUID_padj <- heatmap_data_SUID_padj[c(1,10,21)]
# Remove duplicate gene names - 2 genes , 1370 genes from 1372 genes
heatmap_data_SUID_padj_no_duplicates <- heatmap_data_SUID_padj[!duplicated(heatmap_data_SUID_padj$external_gene_name),]
row.names(heatmap_data_SUID_padj_no_duplicates) <- heatmap_data_SUID_padj_no_duplicates$external_gene_name
# 52,063 genes from 52,986 genes 
heatmap_data_SUID_padj_subset <- subset(heatmap_data_SUID_padj_no_duplicates, select = -c(external_gene_name))
heatmap_data_SUID_padj_subset <- as.matrix(heatmap_data_SUID_padj_subset)
# rename
colnames(heatmap_data_SUID_padj_subset)[1] <- "Non-SUID"
colnames(heatmap_data_SUID_padj_subset)[2]<- "SUID"
# # reorder
# heatmap_reorder <- heatmap_data_subset[, c("Normoxic","Hypoxic","LPS","Hypoxic-LPS")]
# write.table(heatmap_reorder, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Mouse/heatmap_reorder_mouse.txt", sep = "\t")
# # remove 0 expression from dataset, 39,334 from 52,063
# non_zero_expression_genes_apply_mehod_SUID <- heatmap_data_SUID_subset[!apply(heatmap_data_SUID_subset, 1, function(x) all(x == 0)), ]
#heatmap
heatmap_data_SUID_padj <- pheatmap(heatmap_data_SUID_padj_subset, main = "Non-SUID vs SUID", 
         show_rownames=F, cluster_rows=T, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         breaks = seq(-0.5, +0.5, length = 101))
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/heatmap_data_SUID_padj.png", heatmap_data_SUID_padj)


# All HC vs All HP
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/counts_matrix_1.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/sample_annotation_1.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$Gene_ID
count_matrix <- subset(count_matrix, select = -c(Gene_ID) )
sampleCondition <- c("control","control","control","control",
                     "control","control","control",
                     "exp", "exp", "exp","exp","exp",
                     "exp","exp","exp","exp","exp")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)
dds$condition
dds$sample
# Differential Gene Expression
dds_human <- DESeq(dds)
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 1257 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
res <- results(dds_human)
res
# log2 fold change (MLE): condition exp vs control 
# Wald test p-value: condition exp vs control 
summary(res)
# out of 40279 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 0.0025%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

res_contrasts <- results(dds_human, contrast=c("condition","control","exp"))
res_contrasts
summary(res_contrasts)

directory <- setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/")
# HP vs HC
sampleFiles1<- c("deseq2_fcounts_star_out_merged_HC1_1.fqAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC2_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC3_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC4_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC5_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC6_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HC7_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP1_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP2_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP3_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP4_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP5_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP6_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP7_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP8_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP9_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP10_1.fq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_HP11_1.fq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames1 <- c("HC1","HC2","HC3","HC4","HC5","HC6","HC7", "HP1", "HP2","HP3","HP4","HP5","HP6","HP7","HP8","HP9","HP10","HP11")
sampleCondition1 <- c("control","control","control","control","control","control","control", "pre", "pre","pre","pre","pre","pre","pre", "pre","pre","pre","pre")
sampleTable1 <- data.frame(sampleName = sampleNames1, fileName = sampleFiles1, condition1 = sampleCondition1)
treatments1 = c("control","pre")

# Load the data in Deseq format
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable1, directory = directory, design = ~ condition1)
colData(ddsHTSeq1)$condition1 <- factor(colData(ddsHTSeq1)$condition1,levels = treatments1)
#Perform the Differential Gene Expression using DESeq
dds1 <- DESeq(ddsHTSeq1)
res1 <- results(dds1)
res1
summary(res1)


# All HC vs All HP rerun with unstranded featurecounts
count_matrix <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/counts_matrix_rerun.csv', header =T, sep = ",")
sample_annotation <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/HCvsHP/sample_annotation.csv", header =T, sep = ",")
row.names(count_matrix) <- count_matrix$Gene_ID
count_matrix <- subset(count_matrix, select = -c(Gene_ID) )
sampleCondition <- c("control","control","control","control",
                     "control","control","control",
                     "exp", "exp", "exp","exp","exp","exp",
                     "exp","exp","exp","exp","exp")
condition = sampleCondition
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_annotation, 
                              design = ~ condition)
dds$condition
dds$sample
# Differential Gene Expression
# dds_human <- DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3087 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
res <- results(dds_human)
res
# log2 fold change (MLE): condition exp vs control 
# Wald test p-value: condition exp vs control 
summary(res)
# out of 41059 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 570, 1.4%

# non-SUID vs SUID rerun with unstranded featurecounts
count_matrix_SUID <- read.table(file = '/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/SUID/count_matrix_rerun.csv', header =T, sep = ",")
sample_annotation_SUID <- read.table(file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Human/Human2/Deseq2_files/SUID/sample_annotation_SUID.csv", header =T, sep = ",")
row.names(count_matrix_SUID) <- count_matrix_SUID$Gene_ID
count_matrix_SUID <- subset(count_matrix_SUID, select = -c(Gene_ID) )
sampleCondition_SUID <- c("non_SUID","non_SUID","non_SUID","non_SUID","non_SUID","non_SUID",
                          "SUID", "SUID", "SUID","SUID","SUID","SUID",
                          "SUID","SUID","SUID","SUID","SUID","SUID")
condition = sampleCondition_SUID
dds_SUID <- DESeqDataSetFromMatrix(countData = count_matrix_SUID,
                                   colData = sample_annotation_SUID, 
                                   design = ~ condition)

dds_SUID$condition
dds_SUID$sample
# Differential Gene Expression
dds_human_SUID <- DESeq(dds_SUID)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3383 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
res_SUID <- results(dds_human_SUID)
res_SUID
# log2 fold change (MLE): condition SUID vs non SUID 
# Wald test p-value: condition SUID vs non SUID 
summary(res_SUID)
# out of 40891 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 704, 1.7%
# LFC < 0 (down)     : 836, 2%
# outliers [1]       : 190, 0.46%
# low counts [2]     : 20119, 49%


