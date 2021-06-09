install.packages('Seurat')
library(Seurat)
install.packages("dplyr") 
library(dplyr)

# Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gzdata <- Read10X(data.dir = data_dir)seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
data_dir <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/1/outs/raw_feature_bc_matrix'
list.files(data_dir) 
expression_matrix1 <- Read10X(data.dir = data_dir)
seurat_object1 = CreateSeuratObject(counts = expression_matrix1)

#What does data in a count matrix look like?
expression_matrix1[c("PSMB1", "PHF10", "THBS2"), 1:30]

#seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-")

head(seurat_object1@meta.data)

#violin_plot_sample1 <- VlnPlot(seurat_object1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot2 <- FeatureScatter(seurat_object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Filter or subset data
seurat_object1_subset <- subset(seurat_object1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seurat_object1_subset[["RNA"]]

#Normalize
seurat_object1_subset_normalized <- NormalizeData(seurat_object1_subset, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Features
seurat_object1_subset_normalized_variable_features <- FindVariableFeatures(seurat_object1_subset_normalized, selection.method = "vst", nfeatures = 2000)
seurat_object1_subset_normalized_variable_features[["RNA"]]

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object1_subset_normalized_variable_features), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(seurat_object1_subset_normalized_variable_features)
plot3 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0.0, ynudge = 0.0)
CombinePlots(plots = list(plot1, plot2))

#Scaling
all.genes <- rownames(seurat_object1_subset_normalized_variable_features)
scaled_seurat_object1_subset_normalized_variable_features <- ScaleData(seurat_object1_subset_normalized_variable_features, features = all.genes)


scaled <- ScaleData(seurat_object1_subset_normalized_variable_features)

#Linear Dimension Reduction
dim_scaled_norm <- RunPCA(scaled, features = VariableFeatures(object = scaled))

# Examine and visualize PCA results a few different ways
print(dim_scaled_norm[["pca"]], dims = 1:10, nfeatures = 10)


VizDimLoadings(dim_scaled_norm, dims = 1:2, reduction = "pca")

DimPlot(dim_scaled_norm, reduction = "pca")

DimHeatmap(dim_scaled_norm, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(dim_scaled_norm, dims = 1:15, cells = 500, balanced = TRUE)

# Clustering
neighbors <- FindNeighbors(dim_scaled_norm, dims = 1:10)
clusters <- FindClusters(neighbors, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(clusters), 5)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(clusters, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(clusters, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(clusters, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster.markers <- FindAllMarkers(clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)
cluster_markers <- cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Violin Plot (shows expression probability distributions across clusters)
VlnPlot(clusters, features = c("IL10RB", "ATAD2"))
VlnPlot(clusters, features = c("ATAD2"))

# you can plot raw counts as well
VlnPlot(clusters, features = c("IL10RB", "ATAD2"), slot = "counts", log = TRUE)

# FeaturePlot (visualizes feature expression on a tSNE or PCA plot)
FeaturePlot(clusters, features = c("IL10RB", "ATAD2","FRMD4A", "XYLT1"))

# DoHeatmap generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(clusters, features = top10$gene) + NoLegend()


saveRDS(clusters, file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/1.rds")
readRDS(file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/1.rds" )
