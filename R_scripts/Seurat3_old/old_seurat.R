# install.packages('remotes')
# library("remotes")
# find.package('Seurat')
# remotes::install_version("Seurat", version = "3.0.0")
# install.packages("Seurat", lib.loc = "/Users/pbanerjee/Desktop/seurat-3.0.0")
# library()

# # Install biomaRt package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# # Load biomaRt library
# library("biomaRt")
# # List available biomaRt databases
# # Then select a database to use
# listMarts()
# ensembl=useMart("ensembl")
# x <- listDatasets(ensembl)
# ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
# attributes <- listAttributes(ensembl)

library("Seurat")
library("cowplot")
library("dplyr")
library("ggplot2")
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
comb_split <- SplitObject(comb, split.by = "exp")
comb_split <- lapply(X = comb_split, FUN = function(x) {
#  x <- NormalizeData(x)
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
# # Visualization
# p1 <- DimPlot(immune.combined.clusters, reduction = "umap", group.by = "exp")
# p1
# save_plot("clusters_by_group.png",p1)
# p2 <- DimPlot(immune.combined.clusters, reduction = "umap", label = TRUE)
# p2
# save_plot("clusters.png",p2)
# save_plot("clusters_grid.png",plot_grid(p1, p2))
# save_plot("control_vs_treatment.png",DimPlot(immune.combined.clusters, reduction = "umap", split.by = "exp"))

# Identify conserved cell type markers
DefaultAssay(immune.combined.clusters) <- "integrated"
saveRDS(immune.combined.clusters, file = "/Users/pbanerjee/Documents/Payal_Scripts/R/joseph-li-jin-2021/joseph-li-jin-2021.rds")

# Annotation
levels(immune.combined.clusters)
# [1] "0" "1" "2" "3" "4"
# After identification of cluster 1 and 2
# new.cluster.ids <- c("0","Neuronal progenitor cells", "Mitotic neurons", "3","4")
# names(new.cluster.ids) <- levels(immune.combined.clusters)
# immune.combined.clusters <- RenameIdents(immune.combined.clusters, new.cluster.ids)
# DimPlot(immune.combined.clusters, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

comb_markers <- FindAllMarkers(immune.combined.clusters, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25)
colnames(comb_markers)[colnames(comb_markers)=="gene"] <- "ensembl_gene_id"
write.csv(comb_markers,file = "integrated_old_all_markers_without_norm_integrated.csv")
top5_integrated <- comb_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.csv(top5_integrated, file = "top5_integrated_cluster_markers.csv")
DoHeatmap(immune.combined.clusters, features = top5_integrated$gene) + theme(axis.text.y = element_text(size = 6))
save_plot("top5_integrated_heatmap.png",DoHeatmap(immune.combined.clusters, features = top5_integrated$gene) + theme(axis.text.y = element_text(size = 6)))
DoHeatmap(immune.combined.clusters, features = "Grm3", "Sst") + theme(axis.text.y = element_text(size = 6))


# # RNA
# # Identify conserved cell type markers
# DefaultAssay(immune.combined.clusters) <- "RNA"
# # Annotation
# levels(immune.combined.clusters)
# # [1] "0" "1" "2" "3" "4"
# comb_markers <- FindAllMarkers(immune.combined.clusters, only.pos = TRUE, min.pct = 0.25,
#                                logfc.threshold = 0.25)
# colnames(comb_markers)[colnames(comb_markers)=="gene"] <- "ensembl_gene_id"
# write.csv(comb_markers,file = "integrated_old_all_markers_without_norm_RNA.csv")
# top10_rna <- comb_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(top10_rna, file = "top10_RNA_cluster_markers.csv")
# # almost all the markers have infinite value in top 10 markers, so no heatmap is showing up
# save_plot("top10_RNA_heatmap.png",DoHeatmap(immune.combined.clusters, features = top10_rna$ensembl_gene_id) + theme(axis.text.y = element_text(size = 4)))


# Change default assay to integrated again -
# This dataset gets singficant padj values only when switched to integarted
DefaultAssay(immune.combined.clusters) <- "integrated"

# Make the dataet ready by adding custom column to accomodate celltype and experimnet information
immune.combined.clusters$celltype.exp <- paste(Idents(immune.combined.clusters), immune.combined.clusters$exp, sep = "_")
immune.combined.clusters$celltype <- Idents(immune.combined.clusters)
Idents(immune.combined.clusters) <- "celltype.exp"
Idents(immune.combined.clusters)


# The cluster names are not annotated here. The annotation was done afterwards for the final heatmaps 
# 0
deg_0 <- FindMarkers(immune.combined.clusters, ident.1 = "0_treatment", ident.2 = "0_control", verbose = FALSE)
head(deg_0, n = 15)
write.csv(as.data.frame(deg_0),file = "deg_0.csv")

# 1
deg_1 <- FindMarkers(immune.combined.clusters, ident.1 = "1_treatment", ident.2 = "1_control", verbose = FALSE)
head(deg_1, n = 15)
write.csv(as.data.frame(deg_1),file = "deg_1.csv")

# 2
deg_2 <- FindMarkers(immune.combined.clusters, ident.1 = "2_treatment", ident.2 = "2_control", verbose = FALSE)
head(deg_2, n = 15)
write.csv(as.data.frame(deg_2),file = "deg_2.csv")


# 3
deg_3 <- FindMarkers(immune.combined.clusters, ident.1 = "3_treatment", ident.2 = "3_control", verbose = FALSE)
head(deg_3, n = 15)
write.csv(as.data.frame(deg_3),file = "deg_3.csv")


# 4
deg_4 <- FindMarkers(immune.combined.clusters, ident.1 = "4_treatment", ident.2 = "4_control", verbose = FALSE)
head(deg_4, n = 15)
write.csv(as.data.frame(deg_4),file = "deg_4.csv")

#########################################################
### Average expression of a gene in a cluster
#########################################################
# How can I calculate the average expression of all cells within a cluster?
cluster.averages <- AverageExpression(immune.combined.clusters)
head(cluster.averages[["RNA"]][, 1:5])
tail(cluster.averages[["RNA"]][, 1:5])
write.table(cluster.averages[["RNA"]], file = "average_expression_per_cluster_3.2.txt", quote = FALSE, sep = "\t", col.names = TRUE)
write.table(cluster.averages[["integrated"]], file = "average_expression_integrated_per_cluster_3.2.txt", quote = FALSE, sep = "\t", col.names = TRUE)

write.table(immune.combined.clusters@assays[["RNA"]]@counts, file='Gene_Count_per_Cell_3.2.tsv', quote=FALSE, sep='\t', col.names = TRUE)

#########################################################
### Pull number of cells in cluster from seurat object
#########################################################
install.packages("data.table")
library(data.table)
library(magrittr)
## extract meta data
metadata <- immune.combined.clusters@meta.data %>% as.data.table
# the resulting metadata object has one "row" per cell

## count the number of cells per unique combinations of "exp" and "seurat_clusters"
metadata[, .N, by = c("exp", "seurat_clusters")]

## with additional casting after the counting
cells_per_cluster_per_condition <- metadata[, .N, by = c("exp", "seurat_clusters")] %>% dcast(., exp ~ seurat_clusters, value.var = "N")
cells_per_cluster_per_condition
write.table(cells_per_cluster_per_condition, "cells_per_cluster_per_condition.txt",sep = "\t")

new.cluster.ids <- c("")
names(new.cluster.ids) <- levels(immune.combined.clusters)
immune.combined.clusters <- RenameIdents(immune.combined.clusters, new.cluster.ids)
DimPlot(immune.combined.clusters, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#########################################################
##### Subset cluster 1
#########################################################
Idents(immune.combined.clusters)
cluster_1_subset <- subset(immune.combined.clusters, idents = c("1_control", "1_treatment"), invert = FALSE)
cluster1_features = c("ENSMUSG00000003974","ENSMUSG00000004366","ENSMUSG00000006494","ENSMUSG00000006782","ENSMUSG00000008999","ENSMUSG00000015312","ENSMUSG00000015461","ENSMUSG00000021466","ENSMUSG00000021760","ENSMUSG00000023927","ENSMUSG00000024182","ENSMUSG00000026787","ENSMUSG00000027004","ENSMUSG00000029231","ENSMUSG00000030108","ENSMUSG00000030287","ENSMUSG00000031343","ENSMUSG00000031425","ENSMUSG00000041235","ENSMUSG00000043969","ENSMUSG00000044167","ENSMUSG00000074218")
DoHeatmap(cluster_1_subset, features = cluster1_features, size = 3, angle = 0,  draw.lines = T)

#########################################################
##### Subset cluster cluster 2
#########################################################
Idents(immune.combined.clusters)
cluster_2_subset <- subset(immune.combined.clusters, idents = c("2_control", "2_treatment"), invert = FALSE)
cluster2_features = c("ENSMUSG00000001755","ENSMUSG00000001827","ENSMUSG00000008540","ENSMUSG00000008999","ENSMUSG00000010066","ENSMUSG00000015839","ENSMUSG00000019891","ENSMUSG00000021760","ENSMUSG00000022136","ENSMUSG00000022175","ENSMUSG00000022817","ENSMUSG00000024182","ENSMUSG00000026628","ENSMUSG00000027004","ENSMUSG00000030108","ENSMUSG00000030793","ENSMUSG00000031714","ENSMUSG00000032348","ENSMUSG00000033016","ENSMUSG00000033538","ENSMUSG00000033581","ENSMUSG00000035107","ENSMUSG00000035505","ENSMUSG00000035547","ENSMUSG00000040697","ENSMUSG00000043969","ENSMUSG00000048450","ENSMUSG00000035783","ENSMUSG00000050953","ENSMUSG00000074218","ENSMUSG00000064215","ENSMUSG00000064373")
DoHeatmap(cluster_2_subset, features = cluster2_features, size = 3, angle = 0, draw.lines = T)

#########################################################
##### Subset cluster cluster 1 and 2 marker genes
#########################################################
cluster_1_2_subset <- subset(immune.combined.clusters, idents = c("1_control", "1_treatment","2_control", "2_treatment"), invert = FALSE)
features <- c("ENSMUSG00000008658","ENSMUSG00000020436","ENSMUSG00000021585","ENSMUSG00000022340","ENSMUSG00000027748","ENSMUSG00000043635","ENSMUSG00000049001","ENSMUSG00000049336","ENSMUSG00000000632","ENSMUSG00000005871","ENSMUSG00000015222","ENSMUSG00000016200","ENSMUSG00000020181","ENSMUSG00000022054","ENSMUSG00000023927","ENSMUSG00000024261","ENSMUSG00000026787","ENSMUSG00000027273","ENSMUSG00000030209","ENSMUSG00000030731","ENSMUSG00000057897")
DoHeatmap(cluster_1_2_subset, features = features, size = 3)

#########################################################
##### Feature plot cluster 1 
#########################################################
FeaturePlot(immune.combined.clusters, features = c("ENSMUSG00000008658","ENSMUSG00000020436","ENSMUSG00000021585","ENSMUSG00000022340","ENSMUSG00000027748","ENSMUSG00000043635","ENSMUSG00000049001","ENSMUSG00000049336"))
#########################################################
##### Feature plot cluster 2
#########################################################
FeaturePlot(immune.combined.clusters, features = c("ENSMUSG00000000632","ENSMUSG00000005871","ENSMUSG00000015222","ENSMUSG00000016200","ENSMUSG00000020181","ENSMUSG00000022054","ENSMUSG00000023927","ENSMUSG00000024261","ENSMUSG00000026787","ENSMUSG00000027273","ENSMUSG00000030209","ENSMUSG00000030731","ENSMUSG00000057897"))
