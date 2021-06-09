library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(rtracklayer)
library(stringr)

#library(IlluminaHumanMethylationEPICmanifest)
#ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
testDir="X:/Suro/Methylation/Methylation_06232019"
illmn850.anno <- readRDS(file = "X:/Kocher_Limperopoulos_Methylation/Annotation/IlluminaEPIC_annotation_hg38.RDS")
targets <- read.metharray.sheet(testDir, pattern=".csv$")
library("ChAMP")
Sys.setenv(R_MAX_NUM_DLLS=150)
###Default Normalization

testDir="X:\\Suro\\Methylation\\Methylation_06232019"
myLoad<-champ.load(directory = testDir,
           method="minfi",
           methValue="M",
           autoimpute=TRUE,
           filterDetP=TRUE,
           ProbeCutoff=0,
           SampleCutoff=0.1,
           detPcut=0.01,
           filterBeads=TRUE,
           beadCutoff=0.05,
           filterNoCG=FALSE,
           filterSNPs=FALSE,
           population=NULL,
           filterMultiHit=FALSE,
           filterXY=TRUE,
           force=FALSE,
           arraytype="EPIC")	
champ.QC(beta=myLoad$beta,pheno=myLoad$pd$Sample_Group, resultsDir="./CHAMP_CTRL_CHD_P1/")
annotation(rgSet)["array"]="IlluminaHumanMethylationEPIC"
annotation(rgSet)["annotation"]="ilm10b4.hg19"
#getManifest(rgSet)
allAnn=getAnnotation(rgSet) 
annRanges= makeGRangesFromDataFrame(allAnn[,1:4],
                              keep.extra.columns=T,
                              ignore.strand=FALSE,
                              seqinfo=NULL,
                              seqnames.field=c("seqnames", "seqname",
                                               "chromosome", "chrom",
                                               "chr", "chromosome_name",
                                               "seqid"),
                              start.field="pos",
                              end.field="pos",
                              strand.field="strand",
                              starts.in.df.are.0based=FALSE)

chain=import.chain("X:/Kocher_Limperopoulos_Methylation/Annotation/hg19ToHg38.over.chain") # file downloaded from UCSC
hg38Locs=liftOver(annRanges,chain)
hg38LocsDF=data.frame(hg38Locs)
rownames(hg38LocsDF)=hg38LocsDF$group_name
pos38=start(unlist(hg38Locs))
allAnn=data.frame(allAnn,"pos.hg19"=allAnn$pos)
allAnn$pos=rep(NA,dim(allAnn)[1])
allAnn[hg38LocsDF$Name,"pos"]=hg38LocsDF[,"start"]
illmn850.anno <- readRDS(file = 'X:/Kocher_Limperopoulos_Methylation/Annotation/IlluminaEPIC_annotation_hg38.RDS')
targets <- read.metharray.sheet(testDir, pattern=".csv$")
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
detP <- detectionP(rgSet)
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport_07262019.pdf")
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[ ,keep]
mSetRaw <- preprocessRaw(rgSet)
#mSetSq <- preprocessSWAN(rgSet, mSetRaw) 

mSetSq <-preprocessFunnorm(rgSet, nPCs=2, bgCorr = TRUE,
                dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,
                verbose = TRUE)
mSetSq_GR <- mapToGenome(mSetSq, mergeManifest = TRUE) 
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
	   
detP <- detP[match(featureNames(mSetSq_GR),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq_GR) 
table(keep)
mSetSqFlt <- mSetSq_GR[keep,]
mSetSqFlt
keep <- !(featureNames(mSetSqFlt) %in% allAnn$Name[allAnn$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
#mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
 mVals <- getM(mSetSqFlt)
 bVals <- getBeta(mSetSqFlt)
cellType <- factor(targets$Sample_Group)
design <- model.matrix(~0+cellType, data=targets)
#design <- model.matrix(~0+Status:Type, data=targets)
fit <- lmFit(mVals, design)
###CHD_vs_CTL
contMatrix <- makeContrasts(cellTypeCHD-cellTypeCTL, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_up_CHD_CTL_P1 <- topTable(fit2, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_CHD_CTL_P1, file.path(outdir, "DMPs_CHD_CTL_P1.csv"), row.names = FALSE)
###PLA_vs_CTL
contMatrix_PLA_CTL <- makeContrasts(cellTypePLA-cellTypeCTL, levels = design)
fit3 <- contrasts.fit(fit, contMatrix_PLA_CTL)
fit3 <- eBayes(fit3)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
#allAnn=getAnnotation(mVals)
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_up_PLA_CTL_P1 <- topTable(fit3, num=Inf, coef=1, genelist = ann850kSub, p.value=0.1, lfc = 1)
write.csv(DMPs_up_PLA_CTL_P1, file.path(outdir, "DMPs_PLA_CTL_P1.csv"), row.names = FALSE)
###PRE_vs_CTL
 #contMatrix_PRE_CTL <- model.matrix(~Organ*Treatment,data=pData(eset))
#contMatrix_PRE_CTL <- makeContrasts(StatusMF2:TypeCTL-StatusMF2:TypeCHD,StatusP1:TypeCTL-StatusP1:TypeCHD, levels = design)
contMatrix_PRE_CTL <- makeContrasts(cellTypePRE-cellTypeCTL, levels = design)
fit4 <- contrasts.fit(fit, contMatrix_PRE_CTL)
fit4 <- eBayes(fit4)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_up_PRE_CTL_P1 <- topTable(fit4, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_PRE_CTL_P1, file.path(outdir, "DMPs_up_PRE_CTL_P1.csv"), row.names = FALSE)

###MF2 Calculation

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
testDir="X:/Suro/Methylation/Methylation_06232019"
targets <- read.metharray.sheet(testDir, pattern=".csv$")
rgSet <- read.metharray.exp(targets=targets)
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
detP <- detectionP(rgSet)
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport_07292019_MF2.pdf")
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[ ,keep]
mSetRaw <- preprocessRaw(rgSet)
#mSetSq <- preprocessSWAN(rgSet, mSetRaw) 

mSetSq <-preprocessFunnorm(rgSet, nPCs=2, bgCorr = TRUE,
                dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,
                verbose = TRUE)
mSetSq_GR <- mapToGenome(mSetSq) 
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
	   
detP <- detP[match(featureNames(mSetSq_GR),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.05) == ncol(mSetSq_GR) 
table(keep)
mSetSqFlt <- mSetSq_GR[keep,]
mSetSqFlt
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
 mVals <- getM(mSetSqFlt)
 bVals <- getBeta(mSetSqFlt)
cellType <- factor(targets$Sample_Group)
design <- model.matrix(~0+cellType, data=targets)
fit <- lmFit(mVals, design)
###CHD_vs_CTL
contMatrix <- makeContrasts(cellTypeCHD-cellTypeCTL, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]
DMPs_up_CHD_CTL_MF2 <- topTable(fit2, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_CHD_CTL_MF2	, file.path(outdir, "DMPs_up_CHD_CTL_MF2_withoutgDNA.csv"), row.names = FALSE)
###PLA_vs_CTL
contMatrix_PLA_CTL <- makeContrasts(cellTypePLA-cellTypeCTL, levels = design)
fit3 <- contrasts.fit(fit, contMatrix_PLA_CTL)
fit3 <- eBayes(fit3)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]
DMPs_up_PLA_CTL_MF2 <- topTable(fit3, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_PLA_CTL_MF2, file.path(outdir, "DMPs_PLA_CTL_MF2_withoutgDNA.csv"), row.names = FALSE)
###PRE_vs_CTL
 #contMatrix_PRE_CTL <- model.matrix(~Organ*Treatment,data=pData(eset))
'contMatrix_CHD_CTL <- makeContrasts(cellTypePRE-cellTypeCTL, levels = design)
fit4 <- contrasts.fit(fit, contMatrix_PRE_CTL)
fit4 <- eBayes(fit4)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]
DMPs_up_PRE_CTL_P1 <- topTable(fit4, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_PRE_CTL_P1, file.path(outdir, "DMPs_up_PRE_CTL_P1.csv"), row.names = FALSE)'
###TWO Way Anova

testDir="X:/Suro/Methylation/Methylation_06232019"
illmn850.anno <- readRDS(file = 'X:/Kocher_Limperopoulos_Methylation/Annotation/IlluminaEPIC_annotation_hg38.RDS')
targets <- read.metharray.sheet(testDir, pattern=".csv$")
rgSet <- read.metharray.exp(targets=targets)
annotation(rgSet)["array"]="IlluminaHumanMethylationEPIC"
annotation(rgSet)["annotation"]="ilm10b4.hg19"
getManifest(rgSet)
allAnn=getAnnotation(rgSet) #previously was: allAnn=getAnnotation(hgscNorm) but 'hgscNorm' should actually be RBset; I actually want the annotation of this object...
annRanges= makeGRangesFromDataFrame(allAnn[,1:4],
                              keep.extra.columns=T,
                              ignore.strand=FALSE,
                              seqinfo=NULL,
                              seqnames.field=c("seqnames", "seqname",
                                               "chromosome", "chrom",
                                               "chr", "chromosome_name",
                                               "seqid"),
                              start.field="pos",
                              end.field="pos",
                              strand.field="strand",
                              starts.in.df.are.0based=FALSE)

chain=import.chain("X:/Kocher_Limperopoulos_Methylation/Annotation/hg19ToHg38.over.chain") # file downloaded from UCSC
hg38Locs=liftOver(annRanges,chain)
hg38LocsDF=data.frame(hg38Locs)
rownames(hg38LocsDF)=hg38LocsDF$group_name
pos38=start(unlist(hg38Locs))
allAnn=data.frame(allAnn,"pos.hg19"=allAnn$pos)
allAnn$pos=rep(NA,dim(allAnn)[1])
allAnn[hg38LocsDF$Name,"pos"]=hg38LocsDF[,"start"]
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
detP <- detectionP(rgSet)
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport_07262019.pdf")
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[ ,keep]
mSetRaw <- preprocessRaw(rgSet)
#mSetSq <- preprocessSWAN(rgSet, mSetRaw) 

mSetSq <-preprocessFunnorm(rgSet, nPCs=2, bgCorr = TRUE,
                dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,
                verbose = TRUE)
mSetSq_GR <- mapToGenome(mSetSq, mergeManifest = TRUE) 
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)
	   
detP <- detP[match(featureNames(mSetSq_GR),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq_GR) 
table(keep)
mSetSqFlt <- mSetSq_GR[keep,]
mSetSqFlt
keep <- !(featureNames(mSetSqFlt) %in% allAnn$Name[allAnn$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))
 mVals <- getM(mSetSqFlt)
 bVals <- getBeta(mSetSqFlt)
#dupcor <- duplicateCorrelation(mVals, design=design, block=cellType,
                     #trim=0.15, weights=NULL)
design <- model.matrix(~0+Status:cellType , data = targets)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mVals,design)
design2 <- model.matrix(~0+cellType, data=targets)
fit6 <- lmFit(mVals, design2)
###CHD_vs_CTL
contMatrix_CHD_CTL_P1_MF2 <- makeContrasts(StatusMF2.cellTypeCTL-StatusMF2.cellTypeCHD,StatusP1.cellTypeCHD-StatusP1.cellTypeCTL, levels = make.names(colnames(design)))
fit2 <- contrasts.fit(fit, contMatrix_CHD_CTL_P1_MF2)
fit2 <- eBayes(fit2)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_CHD_CTL_P1_MF2 <- topTable(fit2, num=Inf, coef=c(1,2), genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_CHD_CTL_P1_MF2, file.path(outdir, "DMPs_CHD_CTL_P1_MF2.csv"), row.names = FALSE)
###PLA_vs_CTL
contMatrix_CHD_CTL <- makeContrasts(cellTypeCHD-cellTypeCTL, levels = design2)
fit3 <- contrasts.fit(fit6, contMatrix_CHD_CTL)
fit3 <- eBayes(fit3)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
#allAnn=getAnnotation(mVals)
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_up_CHD_CTL <- topTable(fit3, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_CHD_CTL, file.path(outdir, "contMatrix_CHD_CTL.csv"), row.names = FALSE)
'###PRE_vs_CTL
 #contMatrix_PRE_CTL <- model.matrix(~Organ*Treatment,data=pData(eset))

fit4 <- contrasts.fit(fit, contMatrix_PRE_CTL)
fit4 <- eBayes(fit4)
outdir <- "X:/Kocher_Limperopoulos_Methylation/July_2019_analysis/Suro/"
ann850kSub <- illmn850.anno[match(rownames(mVals),illmn850.anno$IlmnID),
                      c(1,6:11)]
DMPs_up_PRE_CTL_P1 <- topTable(fit4, num=Inf, coef=1, genelist = ann850kSub, p.value=0.05, lfc = 2)
write.csv(DMPs_up_PRE_CTL_P1, file.path(outdir, "DMPs_up_PRE_CTL_P1.csv"), row.names = FALSE)'