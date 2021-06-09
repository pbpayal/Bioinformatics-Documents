install.packages('Seurat')
library(Seurat)
install.packages("dplyr") 
library(dplyr)

pbmc <- readRDS(file = "/Users/pbanerjee/Desktop/pbmc3k_final.rds")
set.seed(42)
pbmc$replicate <- sample(c("replicate1_1", "replicate_2"), size = ncol(pbmc), replace = TRUE)

# Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
data_dir1 <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/1/raw_feature_bc_matrix'
list.files(data_dir1) 
rep1_expression_matrix <- Read10X(data.dir = data_dir1)
rep_1 = CreateSeuratObject(counts = rep1_expression_matrix)

# data_dir2 <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/2/raw_feature_bc_matrix'
# list.files(data_dir2) 
# rep2_expression_matrix <- Read10X(data.dir = data_dir2)
# rep_2 = CreateSeuratObject(counts = rep2_expression_matrix)

data_dir3 <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/3/raw_feature_bc_matrix'
list.files(data_dir3) 
rep3_expression_matrix <- Read10X(data.dir = data_dir3)
rep_3 = CreateSeuratObject(counts = rep3_expression_matrix)

data_dir5 <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/5/raw_feature_bc_matrix'
list.files(data_dir5) 
rep5_expression_matrix <- Read10X(data.dir = data_dir5)
rep_5 = CreateSeuratObject(counts = rep5_expression_matrix)

data_dir7 <-'/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/Results/7/raw_feature_bc_matrix'
list.files(data_dir7) 
rep7_expression_matrix <- Read10X(data.dir = data_dir7)
rep_7 = CreateSeuratObject(counts = rep7_expression_matrix)

normoxia <- merge(x = rep_1, y = c(rep_3, rep_5, rep_7), add.cell.ids = c("Nx1","Nx3","Nx5","Nx7"), project = "Normoxia")

head(colnames(normoxia))
tail(colnames(normoxia))

unique(sapply(X = strsplit(colnames(normoxia), split = "_"), FUN = "[", 1))
table(normoxia$orig.ident)

head(normoxia@meta.data)

#Filter or subset data
normoxia_feature_subset <- subset(normoxia, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
normoxia_feature_subset[["RNA"]]

# normoxia_count_subset <- subset(normoxia, subset = nCount_RNA > 200 & nCount_RNA < 2500)
# normoxia_count_subset[["RNA"]]


#Normalize
normoxia_feature_subset_norm <- NormalizeData(normoxia_feature_subset, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Features
normoxia_variable_features <- FindVariableFeatures(normoxia_feature_subset_norm, selection.method = "vst", nfeatures = 2000)
normoxia_variable_features[["RNA"]]

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(normoxia_variable_features), 10)
head(TopCells(object = normoxia))

#Plots
feature_scatter_plot <- FeatureScatter(normoxia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

all.markers_normoxia <- FindAllMarkers(object = normoxia_variable_features)
head(x = all.markers_normoxia)

FindClusters(object = normoxia)

FindMarkers(normoxia_variable_features)

