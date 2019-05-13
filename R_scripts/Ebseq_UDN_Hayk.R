setwd("/Users/pbanerjee/Documents/Payal_Scripts/R")
load("/Users/pbanerjee/Documents/Payal_Scripts/R/Ebseq_UDN_Hayk_GRCh38.RData")
save.image("/Users/pbanerjee/Documents/Payal_Scripts/R/Ebseq_UDN_Hayk_GRCh38.RData")
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeq")
library(EBSeq)
packageVersion("EBSeq")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

#write.csv(MultiPP_Forebrain_heatmap_datamat, "MultiPP_Forebrain_heatmap_datamat1.csv") ##incase to write the data out to a csv


#Patient vs Father(GRCh38)
counts_pf = read.table("/Users/pbanerjee/Desktop/PatientvsFather.txt", stringsAsFactors=F, row.names=1, header=T)
counts_pf_datamat = data.matrix(counts_pf)
str(counts_pf_datamat)

#library size factor
Sizes_PF <-  MedianNorm(counts_pf_datamat)

#set condition matrix
Conditions_PF=as.factor(c("condition","normalF"))

#run Ebseq
EBout_PF <- EBTest(Data=counts_pf_datamat, NgVector = NULL,Conditions = Conditions_PF, sizeFactors = Sizes_PF, maxround = 2, Pool=TRUE )

#get Results
DE_PF <- GetDEResults(EBout_PF, FDR=0.05)
res <- data.frame(DE_PF$PPMat)
str(DE_PF)
str(DE_PF$DEfound)
head(DE_PF$PPMat)
str(DE_PF$Status)

#calculate Fold Change
FC <- PostFC(EBout_PF)
str(FC)

head(FC$RealFC); head(FC$PostFC)

#prepare semi-final results, as some problems were encountered, which were solved below!!
res <- data.frame(DE_PF$PPMat)

#didn't work at first, so troble shooting below
res$PostFC <- FC$PostFC
res$RealFC <- FC$RealFC

#To obtain a matrix without these low expressers
out <- data.frame(DE_PF$PPMat[names(FC$PostFC),]) 
head(out)
out$PostFC <- FC$PostFC 
out$RealFC <- FC$RealFC

#prepare Final Result after troublechooting 
res_all <- data.frame(out)

#Testing
head(res); dim(res); class(res)
length(FC$PostFC)


write.csv(res_all, "/Users/pbanerjee/Desktop/DE_PF.csv")




#It turned out the PostFC doesn't output genes that have 0 counts for all samples. Troubleshooting!!
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

