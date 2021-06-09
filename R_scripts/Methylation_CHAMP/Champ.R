setwd("/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2")
# For a specific project
# load("~/Documents/Payal_Scripts/R/Methylation_Champ/.RData")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest","limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend","IlluminaHumanMethylationEPICmanifest","FEM","matrixStats","missMethyl","combinat"))

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")

library("ChAMP")

# This package is gometh specific for GSEA
# Intsall and load EPIC annoataion separately
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# testDir=system.file("extdata",package="ChAMPdata")
# mytsetLoad <- champ.load(testDir,arraytype="450K")

##############################################
####  LOAD  #### (Normalization - BMIQ, PBC)
##############################################
epic_load <- champ.load(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2_850k",
                         arraytype="EPIC")
CpG.GUI(CpG=rownames(set2$beta),arraytype="EPIC")


##############################################
####  IMPORT  #### (minfi) (Normalization - SWAN, Functional Normalization)
##############################################
set2_import <- champ.import(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2",
                            arraytype = "EPIC")
CpG.GUI(CpG=rownames(set2_import$beta),arraytype="EPIC")


##############################################
####  FILTER  #### (minfi)
##############################################
set2_filter <- champ.filter(set2_import$beta, pd=set2_import$pd, arraytype = "EPIC")


##############################################
####  QUALITY CHECK  ####
##############################################
set2_QC <- champ.QC(beta = set2_import$beta, 
                    pheno=set2_import$pd$Sample_Group,
                    mdsPlot=TRUE,
                    densityPlot=TRUE,
                    dendrogram=TRUE,
                    PDFplot=TRUE,
                    Rplot=TRUE,
                    Feature.sel="None",
                    resultsDir="/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results")
myLoad <- set2_import
QC.GUI(beta = myLoad$beta,arraytype="EPIC" )


##############################################
####  NORMALIZATION - BMIQ  ####
##############################################
#myLoad <- set2_import
myLoad <- set2_load
myNorm <- champ.norm(beta=myLoad$beta,
           rgSet=myLoad$rgSet,
           mset=myLoad$mset,
           resultsDir="/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results",
           method="BMIQ",
           plotBMIQ=TRUE,
           arraytype="EPIC",
           cores=3)
QC.GUI(beta = myNorm,arraytype="EPIC" )
write.csv2(myNorm,"/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/Normalized.csv" )

# champ.SVD(beta = myNorm,
#           rgSet=NULL,
#           pd=myLoad$pd,
#           RGEffect=FALSE,
#           PDFplot=TRUE,
#           Rplot=TRUE,
#           resultsDir="/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results")


##############################################
####  Differntial Methylation Probes  #### (limma)
##############################################
myDMP <- champ.DMP(beta = myNorm,
          pheno = myLoad$pd$Sample_Group,
          compare.group = NULL,
          adjPVal = 0.05,
          adjust.method = "BH",
          arraytype = "EPIC")
head(myDMP[[1]])
write.csv2(myDMP[[1]],"/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/DMP_850_grp1vsgrp2.csv" )
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)


##############################################
####  Differential Methylation Region  #### (minfi)
##############################################
myDMR <- champ.DMR(beta=myNorm,
                   pheno=myLoad$pd$Sample_Group,
                   method="Bumphunter")
head(myDMR[[1]])
write.csv2(myDMR[[1]], "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/DMR_850_grp1vsgrp2.csv")
DMR.GUI(DMR=myDMR)


##############################################
####  Gene Set Enrichment Analysis  ####
##############################################
myGSEA <- champ.GSEA(beta=myNorm,
                     DMP=myDMP[[1]], 
                     DMR=myDMR, 
                     arraytype="EPIC",
                     adjPval=0.05, 
                     method="gometh")
head(myGSEA$DMP)
head(myGSEA$DMP[1])
head(myGSEA$DMP[2])
write.csv2(myGSEA$DMP, "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/850K/GSEA_DMP_850K.csv")
write.csv2(myGSEA$DMP[6], "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/850K/GSEA_DMP_6_850K.csv")
head(myGSEA$DMR)
write.csv2(myGSEA$DMR, "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/GSEA_DMR_850K.csv")


##############################################
####  Ebays Gene Set Enrichment Analysis  ####
##############################################
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
head(myebayGSEA[[1]])
write.csv2(myebayGSEA[[1]], "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/GSEA_ebay_850K.csv")


########################################################
####  Differential Methylated Interaction Hotspots  ####
########################################################

myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group,
                         resultsDir="/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/850K/CHAMP_EpiMod",
                         PDFplot=TRUE,
                         arraytype="EPIC")



########################################################
####  Combine Arrays  ####
########################################################

fourfiftyK_load <- champ.load(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2_450K",
                        arraytype="450K",
                        method = "minfi")

epic_load <- champ.load(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2_850K",
                              arraytype="EPIC",
                              method = "minfi")

# epic_import <- champ.import(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2_850K",
#                         arraytype="EPIC")
# 
# fourfiftyK_import <- champ.import(directory = "/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2_450K",
#                                 arraytype="450K")

comb_array_load <- combineArrays(fourfiftyK_load$rgSet, epic_load$rgSet,
                            outType = c("IlluminaHumanMethylation450k",
                                        "IlluminaHumanMethylationEPIC"),
                            verbose = TRUE)

myLoad <- comb_array_load
##############################################
####  NORMALIZATION - ssNoob  ####
##############################################

noob_normalization <- preprocessNoob(comb_array_load)

##############################################
####  Diffrentially Methylated Probes - minfi(dmpfinder)  ####
##############################################

dmp <- dmpFinder(noob_normalization, pheno = noob_normalization$Sample_Group, qCutoff = 1, type = "continuous")
write.table(dmp,"/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/combined_dmp.txt", sep = "\t")

# dmrs <- bumphunter(noob_normalization, design = model.matrix(noob_normalization$Sample_Group), 
#                    B=0, type="Beta")

"IlluminaHumanMethylationEPICanno.ilm10b2.hg19"
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

myDMR <- champ.DMR(beta=noob_normalization,
                   pheno=noob_normalization$Sample_Group,
                   method="ProbeLasso")

myGSEA <- champ.GSEA(beta=noob_normalization,
                     DMP=dmp[[1]], 
                     DMR=myDMR, 
                     arraytype="EPIC",
                     adjPval=0.05, 
                     method="gometh")
