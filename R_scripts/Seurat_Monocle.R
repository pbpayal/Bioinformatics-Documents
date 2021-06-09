setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/")
save.image(file='Camille_seurat_monocle.RData')
load('Camille_seurat_monocle.RData')

################################################
################################################
##########      INSTALLATION         ###########
################################################
################################################

# To install Bioconductor, open R and run:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
library(BiocManager)
# Next, install a few Bioconductor dependencies that aren't automatically installed:
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))

#Now, install monocle3 through the cole-trapnell-lab GitHub, execute:
install.packages("devtools")
library(devtools)
devtools::install_github('cole-trapnell-lab/leidenbase')
# install.packages("sf")
# library("sf")
devtools::install_url('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz')
# BiocManager::install("batchelor")
devtools::install_github('cole-trapnell-lab/monocle3')
library("monocle3")
library(dplyr)
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

# https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html

################################################
################################################
##########      Data Load  Normox      #########
################################################
################################################

normox <- Read10X("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/normoxia_1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix")
normox[1:10, 1:3]
object.size(normox) # size in bytes
object.size(as.matrix(normox)) # size in bytes
counts_per_cell <- Matrix::colSums(normox)
head(counts_per_cell)
counts_per_gene <- Matrix::rowSums(normox)
head(counts_per_gene)
genes_per_cell <- Matrix::colSums(normox>0) # count gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(normox>0) # only count cells where the gene is expressed
# We are also transforming the values to log10 before plotting, this is done with the log10 method. 
# When logging count data, the + 1 is used to avoid log10(0) which is not defined.
hist(log10(counts_per_cell+1), main='counts per gene', col = 'wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
#histogram of counts per gene in log10 scale?
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

################################################
##########      Seurat Object           ########
################################################

normox_seurat <- CreateSeuratObject(counts = normox, 
                            min.cells = 3, 
                            min.features = 200, 
                            project = "normox")

# report number of genes (rows) and number of cells (columns)
dim(normox_seurat) 
# 18349 (genes) 82503 (cells)
#  We can view the first 10 rows (genes) and the first 10 columns (cells).
normox_seurat@assays$RNA@counts[1:10,1:10]

# compute the proportion of transcripts that are of mitochondrial origin for every cell (percent.mito)
normox_seurat[["percent.mt"]] <- PercentageFeatureSet(normox_seurat, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(normox_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Add number of genes per UMI for each cell to metadata
normox_seurat$log10GenesPerUMI <- log10(normox_seurat$nFeature_RNA) / log10(normox_seurat$nCount_RNA)

################################################
################################################
##########      Data Load  Hypox       #########
################################################
################################################

hypox <- Read10X("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/hypoxia_2_6_8_merge/outs/raw_feature_bc_matrix")
hypox[1:10, 1:3]
object.size(hypox) # size in bytes
# 3069180328 bytes

counts_per_cell <- Matrix::colSums(hypox)
head(counts_per_cell)
counts_per_gene <- Matrix::rowSums(hypox)
head(counts_per_gene)
genes_per_cell <- Matrix::colSums(hypox>0) # count gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(hypox>0) # only count cells where the gene is expressed
# We are also transforming the values to log10 before plotting, this is done with the log10 method. 
# When logging count data, the + 1 is used to avoid log10(0) which is not defined.
hist(log10(counts_per_cell+1), main='counts per gene', col = 'wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
#histogram of counts per gene in log10 scale
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

################################################
###########      Seurat Object      ############
################################################
hypox_seurat <- CreateSeuratObject(counts = hypox, 
                                    min.cells = 3, 
                                    min.features = 200, 
                                    project = "hypox")

# report number of genes (rows) and number of cells (columns)
dim(hypox) 
# 25880 (genes) 20384640 (cells) 
#  We can view the first 10 rows (genes) and the first 10 columns (cells).
hypox_seurat@assays$RNA@counts[1:10,1:10]

# compute the proportion of transcripts that are of mitochondrial origin for every cell (percent.mito)
hypox_seurat[["percent.mt"]] <- PercentageFeatureSet(hypox_seurat, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(hypox_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Add number of genes per UMI for each cell to metadata
hypox_seurat$log10GenesPerUMI <- log10(hypox_seurat$nFeature_RNA) / log10(hypox_seurat$nCount_RNA)

################################################
################################################
##########   Data Merge Hypox Normox   #########
################################################
################################################

normox_hypox <- merge(normox_seurat, y = hypox_seurat,add.cell.ids = c("nromox", "hypox"), project = "normox_hypox")
head(normox_hypox)
dim(normox_hypox)
# 18738 246355
head(colnames(normox_hypox))
head(features@counts)
head(normox_hypox@meta.data)
#  We can view the first 10 rows (genes) and the first 10 columns (cells).
normox_hypox@assays$RNA@counts[1:10,1:10]

head(x = rownames(x = normox_hypox))
head(x = colnames(x = normox_hypox))
tail(x = colnames(x = normox_hypox))

# Gene Names in matrix
gene_names = rownames(x = normox_hypox)
write.table(gene_names, "normox_hypox_gene_names.txt")
library(stringr)
str_detect(gene_names, "^MT")
match(1, str_detect(gene_names, "^MT"))

# Cell Ids
cell_ids = colnames(x = normox_hypox)
write.table(cell_ids, "normox_hypox_cell_ids.txt")

# Normox_Hypox_metadata
normox_hypox_metadata <- read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/Seurat_Normox_Hypox/normox_hypox_metadata.csv")
normox_hypox_metadata_df <- data.frame(normox_hypox_metadata)

table(normox_hypox$orig.ident)
# hypox normox 
# 163852  82503

# Visualize the number of cell counts per cell
normox_hypox_metadata_df %>% 
  ggplot(aes(x=exp, fill=exp)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Split by condition
normox_hypox_split <- SplitObject(normox_hypox, split.by = "orig.ident")
# log-normalization and finding variable features
normox_hypox_split <- lapply(X = normox_hypox_split, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})


normox_hypox_anchors <- FindIntegrationAnchors(object.list = normox_hypox_split, k.filter = 50, dims=1:20)
normox_hypox_combined <- IntegrateData(anchorset = normox_hypox_anchors)

DefaultAssay(comb.combined) <- "integrated"

# compute the proportion of transcripts that are of mitochondrial origin for every cell (percent.mito)
normox_hypox[["percent.mt"]] <- PercentageFeatureSet(normox_hypox, pattern = "^MT-")
head(normox_hypox@meta.data, 5)
tail(normox_hypox@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(normox_hypox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(normox_hypox, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(normox_hypox, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
CombinePlots(plots = list(plot1, plot2))
