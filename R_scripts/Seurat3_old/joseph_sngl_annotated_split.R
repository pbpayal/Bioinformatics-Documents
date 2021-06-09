install.packages('Seurat')
library(Seurat)
install.packages("dplyr") 
# library(dplyr)
install.packages('Rtsne')
library(Rtsne)
install.packages("reticulate")
library(reticulate)
reticulate::py_install(packages = 'umap-learn')
use_condaenv(condaenv="Renv", conda="/home/druss/anaconda3/bin/conda")
library(umap)
library(cowplot)
library("biomaRt")
library("ggplot2")
library("magrittr")
combined_rpkm <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/joseph_combined.csv", 
                          stringsAsFactors = FALSE, header=TRUE, row.names = 1, check.names = FALSE)
str(combined_rpkm [,1:20])

# metadata <- readRDS(file = "/Users/pbanerjee/Downloads/pancreas_v3_files/pancreas_metadata.rds")
joseph_metadata <- read.csv(file = "/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/joseph_metatdata.csv", 
                            header = TRUE, row.names = 1, check.names = FALSE)
joseph_metadata_df <- data.frame(joseph_metadata)
comb <- CreateSeuratObject(counts = combined_rpkm, 
                           project = "Joseph_Comb", min.cells = 3, 
                           meta.data = joseph_metadata_df)


head(comb@meta.data)

# Visualize the number of cell counts per cell
no_of_cells <- joseph_metadata %>% 
  ggplot(aes(x=exp, fill=exp)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
save_plot("No_of_Cells_per_condition.png",no_of_cells )

# Add number of genes per UMI for each cell to metadata
comb$log10GenesPerUMI <- log10(comb$nFeature_RNA) / log10(comb$nCount_RNA)

# Add the colums to metadata
joseph_metadata$log10GenesPerUMI <- comb$log10GenesPerUMI
joseph_metadata$nCount_RNA <- comb$nCount_RNA
joseph_metadata$nFeature_RNA <- comb$nFeature_RNA

# Visualize the number UMIs/transcripts per cell
transcripts_per_cell <- joseph_metadata %>% 
  ggplot(aes(color=exp, x=nCount_RNA, fill= exp)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
save_plot("transcripts_per_cell.png", transcripts_per_cell)


# Visualize the distribution of genes detected per cell via histogram
distribution_of_gene_per_cell <- joseph_metadata %>% 
  ggplot(aes(color=exp, x=nFeature_RNA, fill= exp)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
save_plot("distribution_of_gene_per_cell.png",distribution_of_gene_per_cell)

# Visualize the distribution of genes detected per cell via boxplot
joseph_metadata %>% 
  ggplot(aes(x=exp, y=log10(nFeature_RNA), fill=exp)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
joseph_metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = exp, fill=exp)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


comb_split <- SplitObject(comb, split.by = "exp")


comb_split <- lapply(X = comb_split, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


comb.anchors <- FindIntegrationAnchors(object.list = comb_split, k.filter = 50)
comb.combined <- IntegrateData(anchorset = comb.anchors)

DefaultAssay(comb.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
comb.combined.scaled <- ScaleData(comb.combined, verbose = FALSE)
comb.combined.pca <- RunPCA(comb.combined.scaled, npcs = 30, verbose = FALSE)


# t-SNE and Clustering
immune.combined.UMAP <- RunUMAP(comb.combined.pca, reduction = "pca", dims = 1:20)
immune.combined.neighbors <- FindNeighbors(immune.combined.UMAP, reduction = "pca", dims = 1:20)
immune.combined.clusters <- FindClusters(immune.combined.neighbors, resolution = 1.2)


# Visualization
p1 <- DimPlot(immune.combined.clusters, reduction = "umap", group.by = "exp")
p1
save_plot("clusters_by_group.png",p1)
p2 <- DimPlot(immune.combined.clusters, reduction = "umap", label = TRUE)
p2
save_plot("clusters.png",p2)
save_plot("clusters_grid.png",plot_grid(p1, p2))
save_plot("control_vs_treatment.png",DimPlot(immune.combined.clusters, reduction = "umap", split.by = "exp"))

# Identify conserved cell type markers
DefaultAssay(immune.combined.clusters) <- "RNA"

# FindAllMarkers
comb_markers <- FindAllMarkers(immune.combined.clusters, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25)
colnames(comb_markers)[colnames(comb_markers)=="gene"] <- "ensembl_gene_id"

# Annotate the ENSEMBL GeneIDs
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=comb_markers$ensembl_gene_id, mart= ensembl)
m2 <- merge(comb_markers, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/all_markers_with_norm.csv")



# Annotation
levels(immune.combined.clusters)
# [1] "0" "1" "2" "3" "4"


# Make the dataet ready by adding custom column to accomodate celltype and experimnet information
immune.combined.clusters$celltype.exp <- paste(Idents(immune.combined.clusters), immune.combined.clusters$exp, sep = "_")
immune.combined.clusters$celltype <- Idents(immune.combined.clusters)
Idents(immune.combined.clusters) <- "celltype.exp"

# Oligodendroglial genes
deg_0 <- FindMarkers(immune.combined.clusters, ident.1 = "0_treatment", ident.2 = "0_control", verbose = FALSE)
head(deg_0, n = 15)

# GABA Receptor expressing Neurons
deg_1 <- FindMarkers(immune.combined.clusters, ident.1 = "1_treatment", ident.2 = "1_control", verbose = FALSE)
head(deg_1, n = 15)

# Migrating Neurons
deg_2 <- FindMarkers(immune.combined.clusters, ident.1 = "2_treatment", ident.2 = "2_control", verbose = FALSE)
head(deg_2, n = 15)

# Extracellular Matrix and Signaling components
deg_3 <- FindMarkers(immune.combined.clusters, ident.1 = "3_treatment", ident.2 = "3_control", verbose = FALSE)
head(deg_3, n = 15)

# Metabolic Transporters and channels
deg_4 <- FindMarkers(immune.combined.clusters, ident.1 = "4_treatment", ident.2 = "4_control", verbose = FALSE)
head(deg_4, n = 15)

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
# DEG cluster 0
deg_0$ensembl_gene_id <- rownames(deg_0)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_0), mart= ensembl)
m <- merge(deg_0, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_0.csv")
write.csv(as.data.frame(m),file = "deg_norm_0.csv")

# DEG cluster 1
deg_1$ensembl_gene_id <- rownames(deg_1)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_1), mart= ensembl)
m <- merge(deg_1, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_1.csv")
write.csv(as.data.frame(m),file = "deg_norm_1.csv")

# DEG cluster 2
deg_2$ensembl_gene_id <- rownames(deg_2)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_2), mart= ensembl)
m <- merge(deg_2, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_2.csv")
write.csv(as.data.frame(m),file = "deg_norm_2.csv")

# DEG cluster 3
deg_3$ensembl_gene_id <- rownames(deg_3)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_3), mart= ensembl)
m <- merge(deg_3, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_3.csv")
write.csv(as.data.frame(m),file = "deg_norm_3.csv")

# DEG cluster 4
deg_4$ensembl_gene_id <- rownames(deg_4)
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_4), mart= ensembl)
m <- merge(deg_4, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_4.csv")
write.csv(as.data.frame(m),file = "deg_norm_4.csv")


# Identifying conserved markers allows for identifying only those genes the are significantly 
# differentially expressed relative to the other clusters for all conditions. 
# This function performs differential gene expression testing for a single cluster against 
# all other clusters within each group and then combines the p-values 
# using meta-analysis methods from the MetaDE R package.

# Cluster 0
nk0.markers <- FindConservedMarkers(immune.combined.clusters, ident.1 = 0, grouping.var = "exp", verbose = FALSE)
head(nk0.markers)
write.csv(nk0.markers, "cluster0_markers.csv")
# Convert final results .csv file into .txt file
comb_markers_csv <- "cluster0_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "cluster0_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "cluster0_markers.csv")

# Cluster 1
nk1.markers <- FindConservedMarkers(immune.combined.clusters, ident.1 = 1, grouping.var = "exp", verbose = FALSE)
head(nk1.markers)
write.csv(nk1.markers, "cluster1_markers.csv")
# Convert final results .csv file into .txt file
comb_markers_csv <- "cluster1_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "cluster1_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "cluster1_markers.csv")

# Cluster 2
nk2.markers <- FindConservedMarkers(immune.combined.clusters, ident.1 = 2, grouping.var = "exp", verbose = FALSE)
head(nk2.markers)
write.csv(nk2.markers, "cluster2_markers.csv")
# Convert final results .csv file into .txt file
comb_markers_csv <- "cluster2_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "cluster2_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "cluster2_markers.csv")

# Cluster 3
nk3.markers <- FindConservedMarkers(immune.combined.clusters, ident.1 = 3, grouping.var = "exp", verbose = FALSE)
head(nk3.markers)
write.csv(nk3.markers, "cluster3_markers.csv")
# Convert final results .csv file into .txt file
comb_markers_csv <- "cluster3_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "cluster3_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "cluster3_markers.csv")



# Cluster 4
nk4.markers <- FindConservedMarkers(immune.combined.clusters, ident.1 = 4, grouping.var = "exp", verbose = FALSE)
head(nk4.markers)
write.csv(nk4.markers, "cluster4_markers.csv")
# Convert final results .csv file into .txt file
comb_markers_csv <- "cluster4_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "cluster4_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "cluster4_markers.csv")

# Calculatig the p-adj with Benjamini-Hocheberg

convert_pval_BH <- function(outputfilepath,inputfilepath, outputprefix){
  deg <- read.csv(inputfilepath)
  deg_BH <-  p.adjust(deg$p_val, method = "BH")
  write.csv(deg_BH, file = paste0(outputfilepath, outputprefix, ".csv"))
}

outputprefix = "deg_migrating_neurons_BH"
inputfilepath = "/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/rerun_2021/previous_annotated_clusters/deg_migrating_neurons.csv"
outputfilepath ="/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/rerun_2021/previous_annotated_clusters/"
convert_pval_BH(outputfilepath,inputfilepath, outputprefix)

