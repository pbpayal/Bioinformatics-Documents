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

# Read expression files 
# In joseph's case these are FPKM normalized results from Tophat - Cufflinks pipeline

# control <- read.table(file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_control.csv", sep = ",")
# str(control [,1:20])
# treatment <- read.table(file ="/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_treatment.csv", sep = ",")

# 
# control_rpkm_log <- log(control_rpkm + 1,2) # put on the log scale
# 
# hist(as.vector(as.matrix(control_rpkm)))
# dim(control_rpkm)
# head(rownames(control_rpkm))
# head(colnames(control_rpkm))
# extract top 1000 variable genes
# gene.var <- apply(control_rpkm, 1, function(x) var(x[x>0]))
# control_rpkm.top1000 <- control_rpkm[which(rank(-gene.var)<=1000),]
# 
# rpkm.pca <- prcomp(control_rpkm.top1000,
#                    center = TRUE,
#                    scale. = TRUE) 
# summary(rpkm.pca)$importance[,1:5]
# 
# plot(rpkm.pca, type="l", main="Top 10 PCs")
# 
# #it's often advantageous to run a quick association of the top components of variation with your known variables
# #Association with PC1
# PC1 = rpkm.pca$rotation[,1]
# model.pc1 <- anova(lm(PC1 ~. , df[,-c(2)]))
# #Association with PC2
# PC2 = rpkm.pca$rotation[,2]
# model.pc2 <- anova(lm(PC2 ~. , df[,-c(1)]))
# summary(model.pc2)


####### Seurat ######
#####################

####### Control ######
#####################

control_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_control.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(control_rpkm [,1:20]) # restrict to the first 20 columns (cells)

# Set up control
# ctrl_assay <- CreateAssayObject(counts = control, min.cells = 3)
ctrl <- CreateSeuratObject(counts = control_rpkm, project = "Joseph_CTRL", min.cells = 3)
# ctrl$treat <- "CTRL"
# ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
# ctrl <- NormalizeData(ctrl, verbose = FALSE) # since this is already normaliozed
ctrl_variable_features <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
write.csv(ctrl_variable_feature_id, "ctrl_variable_features")
# Identify the 10 most highly variable genes
ctrl_variable_feature_id <- VariableFeatures(ctrl_variable_features)
top10_ctrl_variable_feature_id <- head(VariableFeatures(ctrl_variable_features), 10)
write.csv(top10_ctrl_variable_feature_id, "top10_ctrl_variable_features")
# plot variable features with and without labels
plot_ctrl_variable_features <- VariableFeaturePlot(ctrl_variable_features)
plot_ctrl_variable_features_top10 <- LabelPoints(plot = plot_ctrl_variable_features, points = top10_ctrl_variable_feature_id, repel = TRUE)
#CombinePlots(plots = list(plot_ctrl_variable_features, plot_ctrl_variable_features_top10))

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

# Set up control
comb <- CreateSeuratObject(counts = combined_rpkm, project = "Joseph_Comb", min.cells = 3, meta.data = joseph_metadata_df)
# ctrl$treat <- "CTRL"
# ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
# ctrl <- NormalizeData(ctrl, verbose = FALSE) # since this is already normaliozed
comb_variable_features <- FindVariableFeatures(comb, selection.method = "vst", nfeatures = 2000)
comb_variable_feature_id <- VariableFeatures(comb_variable_features)
# write.csv(comb_variable_feature_id, "comb_variable_features")
top10_comb_variable_feature_id <- head(VariableFeatures(comb_variable_features), 10)
# write.csv(top10_comb_variable_feature_id, "top10_trt_variable_features")
# plot variable features with and without labels
plot_comb_variable_features <- VariableFeaturePlot(comb_variable_features)
plot_comb_variable_features_top10 <- LabelPoints(plot = plot_comb_variable_features, points = top10_comb_variable_feature_id, repel = TRUE)

#Scale
all.genes.comb <- rownames(comb_variable_features)
scaled_comb_variable_features <- ScaleData(comb_variable_features, features = all.genes.comb)
#PCA
pca_comb <- RunPCA(scaled_comb_variable_features, features = VariableFeatures(object = scaled_comb_variable_features))
# Examine and visualize PCA results a few different ways
print(pca_comb[["pca"]], dims = 1:5, nfeatures = 5)

#Plots and Visualization
vizdimload_comb_plot <- VizDimLoadings(pca_comb, dims = 1:2, reduction = "pca")
vizdimload_comb_plot
dimplot_comb_plot <- DimPlot(pca_comb, reduction = "pca")
dimplot_comb_plot
DimHeatmap(pca_comb, nfeatures = 20, dims = 3, cells = 500, balanced = TRUE)
DimHeatmap(pca_comb, dims = 1:15, cells = 500, balanced = TRUE)

#FindNeighbors
comb_neighbors <- FindNeighbors(pca_comb, dims = 1:10)
#comb_clusters <- FindClusters(comb_neighbors, resolution = 0.5)
comb_clusters <- FindClusters(object = comb_neighbors, reduction.type = "pca", dims.use = 1:10, 
                             resolution = 0.6, print.output = 0, save.SNN = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(comb_clusters), 25)

comb_tsne <- RunTSNE(object = comb_clusters, dims.use = 1:5, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = comb_tsne)

# find all markers of cluster 1
comb.cluster1.markers <- FindMarkers(object = comb_tsne, ident.1 = 1, min.pct = 0.25)
#write.csv(comb.cluster1.markers, "comb_cluster1_markers")
# find markers for every cluster compared to all remaining cells, report only the positive ones
comb_markers <- FindAllMarkers(comb_tsne, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(comb_markers, "comb_markers")



library("biomaRt")
# Convert final results .csv file into .txt file
comb_markers_csv <- "comb_markers.csv"
write.table(read.csv(comb_markers_csv), gsub(".csv",".txt",comb_markers_csv))
comb_markers_txt <- "comb_markers.txt"
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
a2 <- read.table(comb_markers_txt, head=TRUE)
colnames(a2)[colnames(a2)=="X"] <- "ensembl_gene_id"
b2 <- getBM(filters= "ensembl_gene_id", attributes= c("mgi_symbol","ensembl_gene_id", "description"), values=a2$ensembl_gene_id, mart= ensembl)
m2 <- merge(a2, b2, by="ensembl_gene_id")
write.csv(as.data.frame(m2),file = "comb_annotated.csv")



##############################
##############################

metadata <- readRDS(file = "/Users/pbanerjee/Downloads/pancreas_v3_files/pancreas_metadata.rds")
joseph_metadata <- read.csv(file = "/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Joseph/joseph_metatdata.csv")

joseph_metadata_df <- data.frame(joseph_metadata)



comb <- CreateSeuratObject(counts = combined_rpkm, project = "Joseph_Comb", 
                           min.cells = 3, meta.data = joseph_metadata_df)
comb[["RNA"]]@data

comb_split <- SplitObject(comb, split.by = "exp")


comb_variable_features <- FindVariableFeatures(comb, selection.method = "vst", nfeatures = 2000)
comb_variable_feature_id <- VariableFeatures(comb_variable_features)
scaled_comb <- ScaleData(comb_variable_features, features = all.genes.comb)

pca_comb <- RunPCA(scaled_comb_variable_features, features = comb_variable_feature_id)
comb_tsne <- RunTSNE(object = pca_comb, dims.use = 1:5, do.fast = TRUE)
comb_neighbors <- FindNeighbors(comb_tsne, dims = 1:5)
comb_clusters <- FindClusters(object = comb_neighbors, reduction.type = "pca", dims.use = 1:10, 
                              resolution = 0.6, print.output = 0, save.SNN = TRUE)

TSNEPlot(object = comb_clusters)

# Visualization
p1 <- DimPlot(comb_clusters, reduction = "tsne", group.by = "exp", cols.highlight = "#DE2D26")
p2 <- DimPlot(comb_clusters, reduction = "tsne", label = TRUE, cols.highlight = "#DE2D26")
plot_grid(p1, p2)

DimPlot(comb_clusters, reduction = "tsne", split.by = "exp")
