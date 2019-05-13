setwd("/Users/pbanerjee/Desktop/Bioinfo Unit_RNA seq/AB/")
load("Forbrain.RData")
save.image("Forbrain.RData")
library(EBSeq)

#write.csv(MultiPP_Forebrain_heatmap_datamat, "MultiPP_Forebrain_heatmap_datamat1.csv") ##incase to write the data out to a csv

#Forebrain - EBSEQ
Forbrain_data = read.table("all.genes.results.fb.txt", stringsAsFactors=F, row.names=1, header=T)
Forbrain_datamat = data.matrix(Forbrain_data)
str(Forbrain_datamat)

MultiSize <-  MedianNorm(Forbrain_datamat)
Conditions=c("untreated","untreated","untreated","stress","stress","stress","stress_drug","stress_drug","stress_drug","drug","drug","drug")
PosParti=GetPatterns(Conditions)
EBOut_Forbrain <- EBMultiTest(Forbrain_datamat,NgVector = NULL,Conditions = Conditions, AllParti = PosParti, sizeFactors = MultiSize, maxround = 1)
MultiFC_Forbrain <- GetMultiFC(EBOut_Forbrain, SmallNum = 0.01)
dfMultiFC_Forebrain <- data.frame(MultiFC_Forbrain)
MultiPP_Forbrain <- GetMultiPP(EBOut_Forbrain)
dfMultiPP_Forebrain <- data.frame(MultiPP_Forbrain)

#Fold change Heatmaps

##Stress_drug over Untreated 
Stress_drugoverUntreated <- dfMultiFC_Forebrain[,c(1,12)]
Stress_drugoverUntreated <- Stress_drugoverUntreated[-c(1)]
Stress_drugoverUntreated_heatmap_datamat <- data.matrix(Stress_drugoverUntreated)

#All fold change
Foldchange <- dfMultiFC_Forebrain[,c(7:12)]
Foldchange_heatmap_datamat <- data.matrix(Foldchange)
heatmap.2(Foldchange_heatmap_datamat, Rowv = FALSE, Colv = FALSE, hclustfun = hclust, scale = "column", trace = "none")

#Forebrain - HEATMAP
dfMultiPP_Forebrain_heatmap <- dfMultiPP_Forebrain[,c(0,2:15)]
MultiPP_Forebrain_heatmap_datamat <- data.matrix(dfMultiPP_Forebrain_heatmap)
na.omit(MultiPP_Forebrain_heatmap_datamat) #this is because while runing Heatmap, it gave errors. na.omit helps to find out if NA is there in data matrix
MultiPP_Forebrain_heatmap_datamat_NArmovd <- MultiPP_Forebrain_heatmap_datamat[-c(21863),]
heatmap(MultiPP_Forebrain_heatmap_datamat_NArmovd)

#Hindbrain
library(gplots)
hbdata=read.table("all.genes.results.hb.txt", stringsAsFactors=F, row.names=1, header=T)
hbdatamat = data.matrix(hbdata)
str(hbdatamat)

MultiSize=MedianNorm(hbdatamat)
Conditions=c("C1","C1","C1","C2","C2","C2","C3","C3","C3","C4","C4","C4")
PosParti=GetPatterns(Conditions)
EBOut_HB <- EBMultiTest(hbdatamat,NgVector = NULL,Conditions = Conditions, AllParti = PosParti, sizeFactors = MultiSize, maxround = 5)
HBMultiFC <- GetMultiFC(EBOut_HB, SmallNum = 0.01)
dfHBMultiFC <- data.frame(HBMultiFC)
HBMultiPP <- GetMultiPP(EBOut_HB)
dfHBMultiPP <- data.frame(HBMultiPP)

write.csv(dfFBMultiFC, "FBMultiFC.csv")



save.image(file='myEnvironment.RData')
load('myEnvironment.RData')