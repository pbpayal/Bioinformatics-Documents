library("DESeq2")

directory <- setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Pig/Deseq2_files")
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
# Factors are the data objects which are used to categorize the data and store it as levels.
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels = treatments)
#Perform the Differential Gene Expression using DESeq
dds <- DESeq(ddsHTSeq)
View(dds)

# estimation of size factors: estimateSizeFactors (median ratio method)
# estimation of dispersion: estimateDispersions
# Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
# final dispersion estimates
# fitting model and testing

# Normalization
dds_norm <- estimateSizeFactors(dds)
sizeFactors(dds_norm)
normalized_counts <- counts(dds_norm, normalized=TRUE)
normalized_counts <- as.data.frame(counts(dds,normalized =TRUE))

# Extract results from S4 object
res <- results(dds)
head(res)
res_df <- as.data.frame(res)
summary(res)



