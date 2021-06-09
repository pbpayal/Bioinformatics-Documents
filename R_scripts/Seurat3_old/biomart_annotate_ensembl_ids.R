# Install biomaRt package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
# Load biomaRt library
library("biomaRt")
# List available biomaRt databases
# Then select a database to use
listMarts()
ensembl=useMart("ensembl")
x <- listDatasets(ensembl)
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
attributes <- listAttributes(ensembl)


# Annotate old data with BiomartIds
deg_0_old_without_norm_integrated <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/rerun_2021/integrated_old_without_normalization/deg_0.csv", header = T)
colnames(deg_0_old_without_norm_integrated)[colnames(deg_0_old_without_norm_integrated)=="X"] <- "ensembl_gene_id"
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=deg_0_old_without_norm_integrated$ensembl_gene_id, mart= ensembl)
m <- merge(deg_0_old_without_norm_integrated, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_0_old_without_norm_integrated.csv")

