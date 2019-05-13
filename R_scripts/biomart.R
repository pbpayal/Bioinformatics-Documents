load("/Users/pbanerjee/Documents/Payal_Scripts/R/biomart.RData")
setwd("/Users/pbanerjee/Doesktop/")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(biomaRt)
listEnsembl()

data_list <- read.table('HI.2341.001.Index_6.586-5.isoforms.results.txt',  header=T)
data_df <- data.frame(data_list)
x <- data_df$gene_id
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene", "description"),values=x,mart= mart)



#grch37 = useEnsembl(biomart="ensembl",GRCh=37)
#listDatasets(grch37)[31:35,]
#ensembl_human_grch37 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
#listDatasets(ensembl_human_grch37)[31:35,]

save.image("/Users/pbanerjee/Documents/Payal_Scripts/R/biomart.RData")