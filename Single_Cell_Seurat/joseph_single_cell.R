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

####### Seurat ######
#####################

####### Control ######
#####################

# Read RPKM values
control_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_control.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(control_rpkm [,1:20]) # restrict to the first 20 columns (cells)

# Create Seurat Object
ctrl <- CreateSeuratObject(counts = control_rpkm, project = "Joseph_CTRL", min.cells = 3)

#Finding Variable Features
# a subset of features that exhibit high cell-to-cell variation in the dataset 
# (i.e, they are highly expressed in some cells, and lowly expressed in others). 
ctrl_variable_features <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
write.csv(ctrl_variable_feature_id, "ctrl_variable_features")

# Identify the 10 most highly variable genes
ctrl_variable_feature_id <- VariableFeatures(ctrl_variable_features)
top10_ctrl_variable_feature_id <- head(VariableFeatures(ctrl_variable_features), 10)
write.csv(top10_ctrl_variable_feature_id, "top10_ctrl_variable_features")

# plot variable features with and without labels
plot_ctrl_variable_features <- VariableFeaturePlot(ctrl_variable_features)
plot_ctrl_variable_features_top10 <- LabelPoints(plot = plot_ctrl_variable_features, points = top10_ctrl_variable_feature_id, repel = TRUE)


#Scale
all.genes <- rownames(ctrl_variable_features)
scaled_ctrl_variable_features <- ScaleData(ctrl_variable_features, features = all.genes)
#PCA
pca_ctrl <- RunPCA(scaled_ctrl_variable_features, features = VariableFeatures(object = scaled_ctrl_variable_features))
# Examine and visualize PCA results a few different ways
print(pca_ctrl[["pca"]], dims = 1:5, nfeatures = 5)

#Plots and Visualization
vizdimload_ctrl_plot <- VizDimLoadings(pca_ctrl, dims = 1:2, reduction = "pca")
dimplot_ctrl_plot <- DimPlot(pca_ctrl, reduction = "pca")

#FindNeighbors
ctrl_neighbors <- FindNeighbors(pca_ctrl, dims = 1:10)
ctrl_clusters <- FindClusters(ctrl_neighbors, resolution = 0.5)
ctrl_clusters <- FindClusters(object = ctrl_neighbors, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(ctrl_clusters), 25)

ctrl_tsne <- RunTSNE(object = ctrl_clusters, dims.use = 1:5, do.fast = TRUE, perplexity = 1)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = ctrl_tsne)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = ctrl_tsne, ident.1 = 0, min.pct = 0.25)

# find markers for every cluster compared to all remaining cells, report only the positive ones
ctrl_markers <- FindAllMarkers(ctrl_tsne, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

ctrl_cluster1.markers <- FindMarkers(ctrl_clusters, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

####### Treatment ######
#####################


# treatment <- read.table(file ="/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_treatment.csv", sep = ",")

treatment_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_treatment.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(treatment_rpkm [,1:20]) # restrict to the first 20 columns (cells)

# Set up control
trt <- CreateSeuratObject(counts = treatment_rpkm, project = "Joseph_TRT", min.cells = 3)

VlnPlot(trt, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Pearson correlation between the two features
correlation_plot <- FeatureScatter(comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
correlation_plot
# ctrl$treat <- "CTRL"
# ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
# ctrl <- NormalizeData(ctrl, verbose = FALSE) # since this is already normaliozed
trt_variable_features <- FindVariableFeatures(trt, selection.method = "vst", nfeatures = 2000)
trt_variable_feature_id <- VariableFeatures(trt_variable_features)
write.csv(trt_variable_feature_id, "trt_variable_features")
top10_trt_variable_feature_id <- head(VariableFeatures(trt_variable_features), 10)
write.csv(top10_trt_variable_feature_id, "top10_trt_variable_features")
# plot variable features with and without labels
plot_trt_variable_features <- VariableFeaturePlot(trt_variable_features)
plot_trt_variable_features_top10 <- LabelPoints(plot = plot_trt_variable_features, points = top10_trt_variable_feature_id, repel = TRUE)
#CombinePlots(plots = list(plot_ctrl_variable_features, plot_ctrl_variable_features_top10))

#Scale
all.genes.trt <- rownames(trt_variable_features)
scaled_trt_variable_features <- ScaleData(trt_variable_features, features = all.genes.trt)
#PCA
pca_trt <- RunPCA(scaled_trt_variable_features, features = VariableFeatures(object = scaled_trt_variable_features))
# Examine and visualize PCA results a few different ways
print(pca_trt[["pca"]], dims = 1:5, nfeatures = 5)

#Plots and Visualization
vizdimload_trt_plot <- VizDimLoadings(pca_trt, dims = 1:2, reduction = "pca")
dimplot_trt_plot <- DimPlot(pca_trt, reduction = "pca")

#FindNeighbors
trt_neighbors <- FindNeighbors(pca_trt, dims = 1:10)
trt_clusters <- FindClusters(trt_neighbors, resolution = 0.5)
trt_clusters <- FindClusters(object = trt_neighbors, reduction.type = "pca", dims.use = 1:10, 
                              resolution = 0.6, print.output = 0, save.SNN = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(trt_clusters), 25)

trt_tsne <- RunTSNE(object = trt_clusters, dims.use = 1:5, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = trt_tsne)

# find all markers of cluster 1
trt.cluster1.markers <- FindMarkers(object = trt_tsne, ident.1 = 1, min.pct = 0.25)
write.csv(trt.cluster1.markers, "treatment_cluster1_markers")
# find markers for every cluster compared to all remaining cells, report only the positive ones
trt_markers <- FindAllMarkers(trt_tsne, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(trt_markers, "treatment_markers")


################ Combined ###########################
#####################################################


combined_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_combined.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(combined_rpkm [,1:20]) # restrict to the first 20 columns (cells)

# Create Seurat Object
comb <- CreateSeuratObject(counts = combined_rpkm, project = "Joseph_Comb", min.cells = 10, min.features = 100)

# Visualize QC metrics as a violin plot
# nCount is the number of UMI counts taken across all cells from the RNA@data matrix  
# nFeatures counts any gene with at least 1 UMI count.

VlnPlot(comb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Pearson correlation between the two features
correlation_plot <- FeatureScatter(comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
correlation_plot

# comb_subset <- subset(comb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
# VlnPlot(comb_subset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# correlation_plot_subset <- FeatureScatter(comb_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# correlation_plot_subset


# No Normalization, since this is RPKM data (Tophat > Cufflinks pipeline)
# comb_normalized <- NormalizeData(comb, normalization.method = "LogNormalize", scale.factor = 10000)
# VlnPlot(comb_normalized, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Find Variable Features
# measure of a singlecell dispersion after controlling for mean expression, 
# selected the 2,000 genes with the highest standardized variance as "highly variable"
# exhibit high cell-to-cell variation in the dataset 
comb_variable_features <- FindVariableFeatures(comb, selection.method = "vst", num.bin = 20)

# Get and set variable feature information 
comb_variable_feature_id <- VariableFeatures(comb_variable_features)
write.csv(comb_variable_feature_id, "comb_variable_features")
top10_comb_variable_feature_id <- head(VariableFeatures(comb_variable_features), 10)
write.csv(top10_comb_variable_feature_id, "top10_trt_variable_features")
# plot variable features with and without labels
plot_comb_variable_features <- VariableFeaturePlot(comb_variable_features)
plot_comb_variable_features_top10 <- LabelPoints(plot = plot_comb_variable_features, points = top10_comb_variable_feature_id, repel = TRUE)

#Scale
# linear transformation 
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1 
# This step gives equal weight in downstream analyses, 
# so that highly-expressed genes do not dominate

all.genes.comb <- rownames(comb_variable_features)
scaled_comb_variable_features <- ScaleData(comb_variable_features, features = all.genes.comb, model.use = "negbinom", use.umi= TRUE )


# PCA
# linear dimensional reduction
pca_comb <- RunPCA(scaled_comb_variable_features, features = VariableFeatures(object = scaled_comb_variable_features))
# Examine and visualize PCA results a few different ways
print(pca_comb[["pca"]], dims = 1:5, nfeatures = 5)

#Plots and Visualization
vizdimload_comb_plot <- VizDimLoadings(pca_comb, dims = 1:2, reduction = "pca")
dimplot_comb_plot <- DimPlot(pca_comb, reduction = "pca")
DimHeatmap(pca_comb, nfeatures = 20, dims = 3, cells = 500, balanced = TRUE)
DimHeatmap(pca_comb, dims = 1:15, cells = 500, balanced = TRUE)

pca_comb_jackstraw <- JackStraw(pca_comb)
pca_comb_jackstraw <- ScoreJackStraw(pca_comb_jackstraw)
JackStrawPlot(pca_comb_jackstraw)

# Clustering
# graph-based clustering approach
# distance metric which drives the clustering analysis, (based on previously identified PCs) 
# K-nearest neighbor (KNN) graph, based on euclidean distance in PCA space
# edges drawn between cells with similar feature expression patterns
# refine the edge weights between any two cells, 
# based on the shared overlap in their local neighborhoods (Jaccard similarity)
# cluster the cells, we next apply modularity optimization techniques,
# such as the Louvain algorithm (default) or SLM 
# partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities

#FindNeighbors
comb_neighbors <- FindNeighbors(pca_comb, dims = 1:5)
#comb_clusters <- FindClusters(comb_neighbors, resolution = 0.5)
comb_clusters <- FindClusters(object = comb_neighbors, reduction.type = "pca", dims.use = 1:5, 
                             resolution = 1.2, print.output = 0, save.SNN = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(comb_clusters), 25)

comb_tsne <- RunTSNE(object = comb_clusters, dims.use = 1:5, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = comb_tsne)


# Cluster biomarkers
# Finds markers (differentially expressed genes) for each cluster
# identifes positive and negative markers of a single cluster (specified in ident.1), 
# compared to all other cellscompared to all other cells
# default - Wilcoxon Rank Sum test

# find all markers of cluster 1
comb.cluster1.markers <- FindMarkers(object = comb_tsne, ident.1 = 1, min.pct = 0.25)
write.csv(comb.cluster1.markers, "comb_cluster1_markers")
# find markers for every cluster compared to all remaining cells, report only the positive ones
comb_markers <- FindAllMarkers(comb_tsne, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5, return.thresh=0.05)
write.csv(comb_markers, "comb_markers.csv")

library("biomaRt")
# Convert final results .csv file into .txt file
comb_markers_csv <- "comb_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "comb_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "comb_annotated.csv")



library(SeuratData)

