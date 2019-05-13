source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library("sva")

source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
library("DESeq2")

# Cell_Line2
directory <- "/Users/pbanerjee/Desktop/Katrina_counts/Cell_Line2/Deliverables/Raw_counts_cell_line2"
setwd(directory)
# Set the prefix for each output file name
outputPrefix2 <- "Cell2_DESeq2_123"
sampleFiles2 <- c("OC1_S17.txt",
                  "OC2_S18.txt",
                  "OC3_S19.txt",
                  "OE1_S21.txt",
                  "OE2_S22.txt",
                  "OE3_S23.txt")
sampleNames2 <- c("OC1","OC2","OC3","OE1","OE2","OE3")
sampleCondition2 <- c("control","control","control","exp","exp","exp")
sampleTable2 <- data.frame(sampleName = sampleNames2, fileName = sampleFiles2, condition = sampleCondition2)
treatments2 = c("control","exp")
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2, directory = directory, design = ~condition )
colData(ddsHTSeq2)$condition <- factor(colData(ddsHTSeq2)$condition,levels = treatments2)
#guts
dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)
summary(res2)
resOrdered2 <- res2[order(res2$pvalue),]
resSig2_05_2 <- subset(resOrdered2, padj < 0.05)
resSig2_05_2
summary(resSig2_05_2)
write.csv(as.data.frame(resOrdered2), file= paste0(outputPrefix2,"-results.csv"))
sum(res2$padj < 0.1, na.rm=TRUE)
resSig01 <- subset(resOrdered2, padj < 0.1)
resSig05 <- subset(resOrdered2, padj < 0.05)
sum(res2$padj < 0.05, na.rm=TRUE)

write.csv(as.data.frame(resSig01), file= paste0(outputPrefix2,"-results-with-padj.1.csv"))
write.csv(as.data.frame(resSig05), file= paste0(outputPrefix2,"-results-with-padj.05.csv"))

# # import raw data
# myData <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                      directory = pathToHTSeq,
#                                      design = ~Genotype)

# Make rawCounts variable
rawCounts <-data.frame( counts(dds) ) 
write.csv(rawCounts, file = "counts_merged_cellline2.csv")

#import needed library
library('sva')

# dat  <- counts(dds2, normalized = TRUE)
# idx  <- rowMeans(dat) > 1
# dat  <- dat[idx, ]
# mod  <- model.matrix(~ dex, colData(dds))
# mod0 <- model.matrix(~   1, colData(dds))
# svseq <- svaseq(dat, mod, mod0, n.sv = 2)


dds2 <- estimateSizeFactors(dds2)
dat <- counts(dds2, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
design = ~condition + batch
mod <- model.matrix(~condition, colData(dds2))
mod0 <- model.matrix(~ 1, colData(dds2))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
# 
ddssva <- dds2
dds2$SV1 <- svseq$sv[,1]
dds2$SV2 <- svseq$sv[,2]
design(dds2) <- ~ condition + SV1 + SV2

ddssva <- DESeq(ddssva)
#ressva <- results(ddssva,contrast=c("condition","control","experiment"))
res3 <- results(ddssva)
write.csv(as.data.frame(res3), file= paste0(outputPrefix2,"-results-with-batch.csv"))

resOrdered3 <- res3[order(res3$pvalue),]
resSig2_05_3 <- subset(resOrdered3, padj < 0.05)
resSig2_05_3
summary(resSig2_05_3)


vsd <- varianceStabilizingTransformation(ddssva, blind=T)
rld <- rlogTransformation(dds, blind=T)
condition <- treatments2
pcaplot3 <- plotPCA(ddssva, intgroup=c("condition")) + geom_text(aes(label=name),vjust=2)
ggsave(pcaplot2,file=paste0(outputPrefix2, "-ggplot2.pdf"))

# make a full model matrix
mod  <- model.matrix(~ condition, colData(dds2))

# make a null model to compare it to
mod0 <- model.matrix(~   1, colData(dds2))

# perform SVA without defining how many non-condition batch effects you think there are
svseq <- svaseq( as.matrix(rawCounts), mod, mod0, nSurr) 
print(svseq)


# # Normalization in Deseq2
# ddsHTSeq2_rlg <- assay(rlog(ddsHTSeq2))
# tst <- estimateSizeFactors(ddsHTSeq2)
# tst2<-counts(tst, normalized=TRUE)
# colData(tst2)$condition <- factor(colData(tst2)$condition,levels = treatments2)
# dds2 <- DESeq(tst2)
# 
# dds2$batch <- factor(rep(c("control","exp"),each=2))
# vsd <- vst(dds2)
# 
# plotPCA(tst)
# assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
# plotPCA(vsd, intgroup=c("condition")) + geom_text(aes(label=name))
# plotPCA(vsd, "batch")
# 
# batch <- sampleCondition2
# removeBatchEffect(dds2, batch)