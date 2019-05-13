setwd("/Users/pbanerjee/Documents/Payal_Scripts/R")
load("/Users/pbanerjee/Documents/Payal_Scripts/R/Ebseq_UDN_Hayk_GRCh37.RData")
save.image("/Users/pbanerjee/Documents/Payal_Scripts/R/Ebseq_UDN_Hayk_GRCh37.RData")
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeq")
library(EBSeq)
packageVersion("EBSeq")


#Patient vs Father(GRCh37)
counts_pf = read.table("/Users/pbanerjee/Desktop/PatientvsFather_GRCh37.txt", stringsAsFactors=F, row.names=1, header=T)
counts_pf_datamat = data.matrix(counts_pf)
str(counts_pf_datamat)

#library size factor
Sizes_PF <-  MedianNorm(counts_pf_datamat)

#set condition matrix
Conditions_PF=as.factor(c("condition","normalF"))

#run Ebseq
EBout_PF <- EBTest(Data=counts_pf_datamat, NgVector = NULL,Conditions = Conditions_PF, sizeFactors = Sizes_PF, maxround = 2, Pool=TRUE )

#get Results
DE_PF <- GetDEResults(EBout_PF, FDR=0.05, Method="robust", FDRMethod="hard", Threshold_FC=0.7, Threshold_FCRatio=0.3, SmallNum=0.01)
res <- data.frame(DE_PF$PPMat)
#str(DE_PF)
#str(DE_PF$DEfound)
#head(DE_PF$PPMat)
#str(DE_PF$Status)

#calculate Fold Change
FC <- PostFC(EBout_PF)
#str(FC)

#head(FC$RealFC); head(FC$PostFC)

#prepare semi-final results, as some problems were encountered, which were solved below!!
res <- data.frame(DE_PF$PPMat)

#didn't work at first, so troble shooting below
res$PostFC <- FC$PostFC
res$RealFC <- FC$RealFC

#To obtain a matrix without these low expressers
out <- data.frame(DE_PF$PPMat[names(FC$PostFC),]) 
#head(out)
out$PostFC <- FC$PostFC 
out$RealFC <- FC$RealFC

#prepare Final Result after troublechooting 
res_all <- data.frame(out)

res_all_log2FC <- cbind(res_all, log2(res_all[,4]))

#Testing
head(res_all); dim(res_all); class(res_all)
length(FC$PostFC)



write.csv(res_all, "/Users/pbanerjee/Desktop/DE_PF_GRCh38.csv")
write.csv(res_all_log2FC, "/Users/pbanerjee/Desktop/DE_PF_log2FC_GRCh38.csv")

##################################################################################################################
#Patient vs Mother(GRCh37)
counts_pm = read.table("/Users/pbanerjee/Desktop/PatientvsMother.txt", stringsAsFactors=F, row.names=1, header=T)
counts_pm_datamat = data.matrix(counts_pm)
str(counts_pm_datamat)

#library size factor
Sizes_PM <-  MedianNorm(counts_pm_datamat)

#set condition matrix
Conditions_PM=as.factor(c("condition","normalM"))

#run Ebseq
EBout_PM <- EBTest(Data=counts_pm_datamat, NgVector = NULL,Conditions = Conditions_PM, sizeFactors = Sizes_PM, maxround = 2, Pool=TRUE )

#get Results
DE_PM <- GetDEResults(EBout_PM, FDR=0.05, Method="robust", FDRMethod="hard", Threshold_FC=0.7, Threshold_FCRatio=0.3, SmallNum=0.01)
res_m <- data.frame(DE_PM$PPMat)
#str(DE_PM)
#str(DE_PM$DEfound)
#head(DE_PM$PPMat)
#str(DE_PM$Status)

#calculate Fold Change
FC_M <- PostFC(EBout_PM)
#str(FC_M)

#head(FC_M$RealFC); head(FC_M$PostFC)

#prepare semi-final results, as some problems were encountered, which were solved below!!
res_m <- data.frame(DE_PM$PPMat)

#didn't work at first, so troble shooting below
res_m$PostFC <- FC_M$PostFC
res_m$RealFC <- FC_M$RealFC

#To obtain a matrix without these low expressers
out_m <- data.frame(DE_PM$PPMat[names(FC_M$PostFC),]) 
#head(out)
out_m$PostFC <- FC_M$PostFC 
out_m$RealFC <- FC_M$RealFC

#prepare Final Result after troublechooting 
res_all_m <- data.frame(out_m)

res_all_m_log2FC <- cbind(res_all_m, log2(res_all_m[,4]))

#Testing
head(res_all_m); dim(res_all_m); class(res_all_m)
length(FC_M$PostFC)

write.csv(res_all_m, "/Users/pbanerjee/Desktop/DE_PM_GRCh37.csv")
write.csv(res_all_m_log2FC, "/Users/pbanerjee/Desktop/DE_PM_log2FC_GRCh37.csv")

#Troubleshooting!!
#It turned out the PostFC doesn't output genes that have 0 counts for all samples. 
Diff <- setdiff(rownames(DE_PF$PPMat), names(FC[[1]]))
str(Diff)

summary(rowMeans(dat[Diff,])) #something is wrong with this command

sum(!is.na(DE_PF$PPMat[Diff,1])) 

table(DE_PF$Status[Diff]) 

listEnsembl()
listEnsembl(GRCh=37)
grch37 = useEnsembl(biomart="ensembl",GRCh=37)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
head(listFilters(ensembl))
head(listAttributes(ensembl))

