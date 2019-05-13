source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# Load DESeq2 library
library("DESeq2")

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