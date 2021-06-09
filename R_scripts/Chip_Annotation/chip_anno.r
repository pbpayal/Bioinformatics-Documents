## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library("ChIPpeakAnno")

macs_peak4 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Results/macs2_sample4/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_peak4_Output <- toGRanges(macs_peak4, format="MACS")

macs_peak4_Output[1:2]


source("https://bioconductor.org/biocLite.R")
biocLite("EnsDb.Hsapiens.v86")
library("EnsDb.Hsapiens.v86")
annoData <- toGRanges(EnsDb.Hsapiens.v86)
annoData[1:2]


## keep the seqnames in the same style
seqlevelsStyle(macs_peak4_Output) <- seqlevelsStyle(annoData)
## do annotation by nearest TSS
anno_peak4 <- annotatePeakInBatch(macs_peak4_Output, AnnotationData=annoData)
anno_peak4[1:2]


# A pie chart can be used to demonstrate the overlap features of the peaks.
pie1(table(anno_peak4$insideFeature))


source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

anno_peak4 <- addGeneIDs(anno_peak4, orgAnn="org.Hs.eg.db",
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
head(anno_peak4)
# 
# # Annotate the peaks with promoters provided by TxDb
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# 
# annoData_TXDB <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene)
# annoData_TXDB[1:2]


# Annotate the peaks in both sides with nearest transcription start sites within 5K bps.
anno_peak4 <- annotatePeakInBatch(macs_peak4_Output, AnnotationData=annoData, 
                            output="nearestBiDirectionalPromoters", 
                            bindingRegion=c(-5000, 500))
anno_peak4$symbol <- xget(anno_peak4$feature, org.Hs.egSYMBOL)
anno_peak4[anno_peak4$peak=="macs2__peak_902"]



macs_sample1 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample1_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample1_peak <- toGRanges(macs_sample1, format="MACS")
macs_sample1_peak[1:2]
seqlevelsStyle(macs_sample1_peak) <- seqlevelsStyle(annoData)
anno_sample1 <- annotatePeakInBatch(macs_sample1_peak, AnnotationData=annoData, 
                                  output="nearestBiDirectionalPromoters", 
                                  bindingRegion=c(-5000, 500))
anno_sample1 <- addGeneIDs(anno_sample1, orgAnn="org.Hs.eg.db",
                         feature_id_type="ensembl_gene_id",
                         IDs2Add=c("symbol"))
head(anno_sample1)
write.table(anno_sample1, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample1.txt")
anno_sample1$symbol <- xget(anno_sample1$feature, org.Hs.egSYMBOL)
anno_sample1[anno_sample1$peak=="macs2__peak_2"]
pie1(table(anno_sample1$insideFeature))

# Annotate with GO terms
anno_sample1_withGO <- getEnrichedGO(anno_sample1, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample1_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample1.txt")
write.table(anno_sample1_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample1.txt")
write.table(anno_sample1_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample1.txt")

# biological process
head(anno_sample1_withGO[["bp"]][, -3])
# molecular function
head(anno_sample1_withGO[["mf"]][, -3])
# cellular component
head(anno_sample1_withGO[["cc"]][, -3])

source("https://bioconductor.org/biocLite.R")
biocLite("KEGG.db")
library("KEGG.db")


source("https://bioconductor.org/biocLite.R")
biocLite("reactome.db")
library("reactome.db")


# Annotate for KEGG Pathways
anno_sample1_with_reactome <- getEnrichedPATH(anno_sample1, orgAnn="org.Hs.eg.db", 
                                       pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                       maxP=0.01, minPATHterm=10, multiAdjMethod="BH")

write.table(anno_sample1_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample1.txt")

macs_sample3 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample3.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample3_peak <- toGRanges(macs_sample3, format="MACS")
macs_sample3_peak[1:2]
seqlevelsStyle(macs_sample3_peak) <- seqlevelsStyle(annoData)
anno_sample3 <- annotatePeakInBatch(macs_sample3_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
anno_sample3 <- addGeneIDs(anno_sample3, orgAnn="org.Hs.eg.db",
                           feature_id_type="ensembl_gene_id",
                           IDs2Add=c("symbol"))
write.table(anno_sample3, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample3.txt")
anno_sample3$symbol <- xget(anno_sample3$feature, org.Hs.egSYMBOL)
anno_sample3[anno_sample3$peak=="macs2__peak_2"]
pie1(table(anno_sample3$insideFeature))

# Annotate with GO terms
anno_sample3_withGO <- getEnrichedGO(anno_sample3, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample3_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample3.txt")
write.table(anno_sample3_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample3.txt")
write.table(anno_sample3_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample3.txt")

# Annotate for KEGG Pathways
anno_sample3_with_reactome <- getEnrichedPATH(anno_sample3, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample3_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample3.txt")

macs_sample4 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample4_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample4_peak <- toGRanges(macs_sample4, format="MACS")
macs_sample4_peak[1:2]
seqlevelsStyle(macs_sample4_peak) <- seqlevelsStyle(annoData)
anno_sample4 <- annotatePeakInBatch(macs_sample4_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample4, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample4.txt")
anno_sample4$symbol <- xget(anno_sample4$feature, org.Hs.egSYMBOL)
anno_sample4[anno_sample4$peak=="macs2__peak_2"]
pie1(table(anno_sample4$insideFeature))

# Annotate with GO terms
anno_sample4_withGO <- getEnrichedGO(anno_sample4, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample4_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample4.txt")
write.table(anno_sample4_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample4.txt")
write.table(anno_sample4_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample4.txt")

# Annotate for KEGG Pathways
anno_sample4_with_reactome <- getEnrichedPATH(anno_sample4, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample4_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample4.txt")

macs_sample6 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample6_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample6_peak <- toGRanges(macs_sample6, format="MACS")
macs_sample6_peak[1:2]
seqlevelsStyle(macs_sample6_peak) <- seqlevelsStyle(annoData)
anno_sample6 <- annotatePeakInBatch(macs_sample6_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample6, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample6.txt")
anno_sample6$symbol <- xget(anno_sample6$feature, org.Hs.egSYMBOL)
anno_sample6[anno_sample6$peak=="macs2__peak_2"]
pie1(table(anno_sample6$insideFeature))

# Annotate with GO terms
anno_sample6_withGO <- getEnrichedGO(anno_sample6, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample6_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample6.txt")
write.table(anno_sample6_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample6.txt")
write.table(anno_sample6_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample6.txt")

# Annotate for KEGG Pathways
anno_sample6_with_reactome <- getEnrichedPATH(anno_sample6, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample6_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample6.txt")

macs_sample7 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample7_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample7_peak <- toGRanges(macs_sample7, format="MACS")
macs_sample7_peak[1:2]
seqlevelsStyle(macs_sample7_peak) <- seqlevelsStyle(annoData)
anno_sample7 <- annotatePeakInBatch(macs_sample7_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample7, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample7.txt")
anno_sample7$symbol <- xget(anno_sample7$feature, org.Hs.egSYMBOL)
anno_sample7[anno_sample7$peak=="macs2__peak_2"]
pie1(table(anno_sample7$insideFeature))

# Annotate with GO terms
anno_sample7_withGO <- getEnrichedGO(anno_sample7, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample7_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample7.txt")
write.table(anno_sample7_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample7.txt")
write.table(anno_sample7_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample7.txt")

# Annotate for KEGG Pathways
anno_sample7_with_reactome <- getEnrichedPATH(anno_sample7, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample7_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample7.txt")

macs_sample8 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample8_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample8_peak <- toGRanges(macs_sample8, format="MACS")
macs_sample8_peak[1:2]
seqlevelsStyle(macs_sample8_peak) <- seqlevelsStyle(annoData)
anno_sample8 <- annotatePeakInBatch(macs_sample8_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample8, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample8.txt")
anno_sample8$symbol <- xget(anno_sample8$feature, org.Hs.egSYMBOL)
anno_sample8[anno_sample8$peak=="macs2__peak_2"]
pie1(table(anno_sample8$insideFeature))

# Annotate with GO terms
anno_sample8_withGO <- getEnrichedGO(anno_sample8, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample8_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample8.txt")
write.table(anno_sample8_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample8.txt")
write.table(anno_sample8_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample8.txt")

# Annotate for KEGG Pathways
anno_sample8_with_reactome <- getEnrichedPATH(anno_sample8, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample8_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample8.txt")

macs_sample9 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample9_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample9_peak <- toGRanges(macs_sample9, format="MACS")
macs_sample9_peak[1:2]
seqlevelsStyle(macs_sample9_peak) <- seqlevelsStyle(annoData)
anno_sample9 <- annotatePeakInBatch(macs_sample9_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample9, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample9.txt")
anno_sample9$symbol <- xget(anno_sample9$feature, org.Hs.egSYMBOL)
anno_sample9[anno_sample9$peak=="macs2__peak_2"]
pie1(table(anno_sample9$insideFeature))

# Annotate with GO terms
anno_sample9_withGO <- getEnrichedGO(anno_sample9, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample9_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample9.txt")
write.table(anno_sample9_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample9.txt")
write.table(anno_sample9_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample9.txt")

# Annotate for KEGG Pathways
anno_sample9_with_reactome <- getEnrichedPATH(anno_sample9, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample9_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample9.txt")

macs_sample11 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample11_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample11_peak <- toGRanges(macs_sample11, format="MACS")
macs_sample11_peak[1:2]
seqlevelsStyle(macs_sample11_peak) <- seqlevelsStyle(annoData)
anno_sample11 <- annotatePeakInBatch(macs_sample11_peak, AnnotationData=annoData, 
                                    output="nearestBiDirectionalPromoters", 
                                    bindingRegion=c(-5000, 500))
write.table(anno_sample11, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample11.txt")
anno_sample11$symbol <- xget(anno_sample11$feature, org.Hs.egSYMBOL)
anno_sample11[anno_sample11$peak=="macs2__peak_2"]
pie1(table(anno_sample11$insideFeature))

# Annotate with GO terms
anno_sample11_withGO <- getEnrichedGO(anno_sample11, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample11_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample11.txt")
write.table(anno_sample11_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample11.txt")
write.table(anno_sample11_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample11.txt")

# Annotate for KEGG Pathways
anno_sample6_with_reactome <- getEnrichedPATH(anno_sample6, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample6_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample6.txt")


macs_sample12 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample12_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample12_peak <- toGRanges(macs_sample12, format="MACS")
macs_sample12_peak[1:2]
seqlevelsStyle(macs_sample12_peak) <- seqlevelsStyle(annoData)
anno_sample12 <- annotatePeakInBatch(macs_sample12_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample12, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample12.txt")
anno_sample12$symbol <- xget(anno_sample12$feature, org.Hs.egSYMBOL)
anno_sample12[anno_sample12$peak=="macs2__peak_27"]
pie1(table(anno_sample12$insideFeature))

# Annotate with GO terms
anno_sample12_withGO <- getEnrichedGO(anno_sample12, orgAnn="org.Hs.eg.db", 
                                     maxP=0.01, minGOterm=10, 
                                     multiAdjMethod="BH",
                                     condense=FALSE)
write.table(anno_sample12_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample12.txt")
write.table(anno_sample12_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample12.txt")
write.table(anno_sample12_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample12.txt")

# Annotate for KEGG Pathways
anno_sample12_with_reactome <- getEnrichedPATH(anno_sample12, orgAnn="org.Hs.eg.db", 
                                              pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                              maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample12_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample12.txt")


macs_sample13 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample13_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample13_peak <- toGRanges(macs_sample13, format="MACS")
macs_sample13_peak[1:2]
seqlevelsStyle(macs_sample13_peak) <- seqlevelsStyle(annoData)
anno_sample13 <- annotatePeakInBatch(macs_sample13_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample13, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample13.txt")
anno_sample13$symbol <- xget(anno_sample13$feature, org.Hs.egSYMBOL)
anno_sample13[anno_sample13$peak=="macs2__peak_27"]
pie1(table(anno_sample13$insideFeature))

# Annotate with GO terms
anno_sample13_withGO <- getEnrichedGO(anno_sample13, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample13_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample13.txt")
write.table(anno_sample13_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample13.txt")
write.table(anno_sample13_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample13.txt")

# Annotate for KEGG Pathways
anno_sample13_with_reactome <- getEnrichedPATH(anno_sample13, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample13_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample13.txt")

macs_sample14 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample14_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample14_peak <- toGRanges(macs_sample14, format="MACS")
macs_sample14_peak[1:2]
seqlevelsStyle(macs_sample14_peak) <- seqlevelsStyle(annoData)
anno_sample14 <- annotatePeakInBatch(macs_sample14_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample14, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample14.txt")
anno_sample14$symbol <- xget(anno_sample14$feature, org.Hs.egSYMBOL)
anno_sample14[anno_sample14$peak=="macs2__peak_27"]
pie1(table(anno_sample14$insideFeature))

# Annotate with GO terms
anno_sample14_withGO <- getEnrichedGO(anno_sample14, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample14_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample14.txt")
write.table(anno_sample14_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample14.txt")
write.table(anno_sample14_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample14.txt")

# Annotate for KEGG Pathways
anno_sample14_with_reactome <- getEnrichedPATH(anno_sample14, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample14_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample14.txt")

macs_sample16 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample16_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample16_peak <- toGRanges(macs_sample16, format="MACS")
macs_sample16_peak[1:2]
seqlevelsStyle(macs_sample16_peak) <- seqlevelsStyle(annoData)
anno_sample16 <- annotatePeakInBatch(macs_sample16_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample16, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample16.txt")
anno_sample16$symbol <- xget(anno_sample16$feature, org.Hs.egSYMBOL)
anno_sample16[anno_sample16$peak=="macs2__peak_27"]
pie1(table(anno_sample16$insideFeature))

# Annotate with GO terms
anno_sample16_withGO <- getEnrichedGO(anno_sample16, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample16_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample16.txt")
write.table(anno_sample16_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample16.txt")
write.table(anno_sample16_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample16.txt")

# Annotate for KEGG Pathways
anno_sample16_with_reactome <- getEnrichedPATH(anno_sample16, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample16_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample16.txt")


macs_sample17 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample17_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample17_peak <- toGRanges(macs_sample17, format="MACS")
macs_sample17_peak[1:2]
seqlevelsStyle(macs_sample17_peak) <- seqlevelsStyle(annoData)
anno_sample17 <- annotatePeakInBatch(macs_sample17_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample17, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample17.txt")
anno_sample17$symbol <- xget(anno_sample17$feature, org.Hs.egSYMBOL)
anno_sample17[anno_sample17$peak=="macs2__peak_2"]
pie1(table(anno_sample17$insideFeature))

# Annotate with GO terms
anno_sample17_withGO <- getEnrichedGO(anno_sample17, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample17_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample17.txt")
write.table(anno_sample17_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample17.txt")
write.table(anno_sample17_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample17.txt")

# Annotate for KEGG Pathways
anno_sample17_with_reactome <- getEnrichedPATH(anno_sample17, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample17_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample17.txt")


macs_sample18 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample18_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample18_peak <- toGRanges(macs_sample18, format="MACS")
macs_sample18_peak[1:2]
seqlevelsStyle(macs_sample18_peak) <- seqlevelsStyle(annoData)
anno_sample18 <- annotatePeakInBatch(macs_sample18_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
write.table(anno_sample18, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample18.txt")
anno_sample18$symbol <- xget(anno_sample18$feature, org.Hs.egSYMBOL)
anno_sample18[anno_sample18$peak=="macs2__peak_2"]
pie1(table(anno_sample18$insideFeature))

# Annotate with GO terms
anno_sample18_withGO <- getEnrichedGO(anno_sample18, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample18_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample18.txt")
write.table(anno_sample18_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample18.txt")
write.table(anno_sample18_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample18.txt")

# Annotate for KEGG Pathways
anno_sample18_with_reactome <- getEnrichedPATH(anno_sample18, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample18_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample18.txt")


macs_sample19 <- read.table(file="/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/Macs2_Broadpeaks_results/macs2_broad_sorted_aligned_merged_sample19_R1.fastq.gz.sam.bam/macs2__peaks.xls",  sep = "\t", header=TRUE)
macs_sample19_peak <- toGRanges(macs_sample19, format="MACS")
macs_sample19_peak[1:2]
seqlevelsStyle(macs_sample19_peak) <- seqlevelsStyle(annoData)
anno_sample19 <- annotatePeakInBatch(macs_sample19_peak, AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters", 
                                     bindingRegion=c(-5000, 500))
anno_sample19 <- addGeneIDs(macs_sample19_peak, orgAnn="org.Hs.eg.db",
                           feature_id_type="ensembl_gene_id",
                           IDs2Add=c("symbol"))
write.table(anno_sample19, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_sample19.txt")
anno_sample19$symbol <- xget(anno_sample19$feature, org.Hs.egSYMBOL)
anno_sample19[anno_sample19$peak=="macs2__peak_2"]
pie1(table(anno_sample19$insideFeature))

# Annotate with GO terms
anno_sample19_withGO <- getEnrichedGO(anno_sample19, orgAnn="org.Hs.eg.db", 
                                      maxP=0.01, minGOterm=10, 
                                      multiAdjMethod="BH",
                                      condense=FALSE)
write.table(anno_sample19_withGO$bp, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-bp_sample19.txt")
write.table(anno_sample19_withGO$mf, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-mf_sample19.txt")
write.table(anno_sample19_withGO$cc, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_GO-cc_sample19.txt")

# Annotate for KEGG Pathways
anno_sample19_with_reactome <- getEnrichedPATH(anno_sample19, orgAnn="org.Hs.eg.db", 
                                               pathAnn="reactome.db", feature_id_type="ensembl_gene_id",
                                               maxP=0.01, minPATHterm=10, multiAdjMethod="BH")
write.table(anno_sample19_with_reactome, "/Users/pbanerjee/Documents/CBU/CNMC_Projects/GenMed/Kate_Chip/annotated_pathway_sample19.txt")




