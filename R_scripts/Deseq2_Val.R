source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# Load DESeq2 library
library("DESeq2")

directory <- "/Users/pbanerjee/Desktop/Val_counts/Part2_Analysis/6hrs"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Val_Results_6hrs"

sampleFiles<- c("m-6hrs-1.txt",
                "C-6hrs-2.txt",
                "C-6hrs-3.txt",
                "m-6hrs-1_part1.txt",
                "m-6hrs-2.txt",
                "m-6hrs-3_part1.txt")

sampleNames <- c("C1","C2","C3","M1","M2","M3")
sampleCondition <- c("control","control","control","exp","exp","exp")
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
summary(res)


# We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)
#write.csv(resOrdered, file = paste0(outputPrefix, "-results-ordered.csv"))

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)
sum(res$pvalue <0.05, na.rm=TRUE)

# p-adj factor alpha set to 0.05 instaed of 0.1 to be more stringent
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$pvalue <0.05, na.rm=TRUE)

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
results_csv <- "Val_Results-results-with-normalized.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Val_Results-results-with-normalized.txt"

# Convert padj= 0.05 results .csv file into .txt file
results_csv05 <- "Val_Results-results-with-padj.05.csv"
write.table(read.csv(results_csv05), gsub(".csv",".txt",results_csv05))
results_txt05 <- "Val_Results-results-with-padj.05.txt"

# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
# List available datasets in database
# Then select dataset to use
listDatasets(ensembl)

# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a05 <- read.table(results_txt05, head=TRUE)
colnames(a05)[colnames(a05)=="X"] <- "ensembl_gene_id"

# listFilters()
# listAttributes()
b05 <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol","ensembl_gene_id","entrezgene", "description"), values=a05$ensembl_gene_id, mart= ensembl)



m <- merge(a05, b05, by="ensembl_gene_id")

write.csv(as.data.frame(m),file = "Val_Results-Biomart-results_merged_padj05.csv")
