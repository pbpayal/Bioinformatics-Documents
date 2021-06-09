install.packages('Seurat')
library(Seurat)
install.packages("dplyr") 
library(dplyr)
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

combined_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_combined.csv", 
                          stringsAsFactors = FALSE, header=TRUE, row.names = 1, check.names = FALSE)
str(combined_rpkm [,1:20])

# metadata <- readRDS(file = "/Users/pbanerjee/Downloads/pancreas_v3_files/pancreas_metadata.rds")
joseph_metadata <- read.csv(file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_metatdata.csv", 
                            header = TRUE, row.names = 1, check.names = FALSE)
joseph_metadata_df <- data.frame(joseph_metadata)

comb <- CreateSeuratObject(counts = combined_rpkm, 
                           project = "Joseph_Comb", min.cells = 3, 
                           meta.data = joseph_metadata_df)


comb@meta.data

# Visualize the number of cell counts per cell
joseph_metadata %>% 
  ggplot(aes(x=exp, fill=exp)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Add number of genes per UMI for each cell to metadata
comb$log10GenesPerUMI <- log10(comb$nFeature_RNA) / log10(comb$nCount_RNA)


# Add the colums to metadata
joseph_metadata$log10GenesPerUMI <- comb$log10GenesPerUMI
joseph_metadata$nCount_RNA <- comb$nCount_RNA
joseph_metadata$nFeature_RNA <- comb$nFeature_RNA

# Visualize the number UMIs/transcripts per cell
joseph_metadata %>% 
  ggplot(aes(color=exp, x=nCount_RNA, fill= exp)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
joseph_metadata %>% 
  ggplot(aes(color=exp, x=nFeature_RNA, fill= exp)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

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
p2 <- DimPlot(immune.combined.clusters, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined.clusters, reduction = "umap", split.by = "exp")


# Identify conserved cell type markers
DefaultAssay(immune.combined.clusters) <- "RNA"

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

# FindAllMarkers
comb_markers <- FindAllMarkers(immune.combined.clusters, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25)

# Convert final results .csv file into .txt file
comb_markers_csv <- "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/Seurat_results/Split/all_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/Seurat_results/Split/all_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/Seurat_results/Split/all_markers.csv")



