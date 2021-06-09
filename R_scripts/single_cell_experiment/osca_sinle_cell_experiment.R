##################################
###### OSCA Tutorial #############
##################################
# http://bioconductor.org/books/release/OSCA/data-infrastructure.html

BiocManager::install('SingleCellExperiment')
library(SingleCellExperiment)
BiocManager::install(c('scater', 'scran', 'uwot'))

counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix) # must be a matrix object!
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
sce
# class: SingleCellExperiment 
# dim: 10 3 
# metadata(0):
# assays(1): counts
# rownames(10): gene_1 gene_2 ... gene_9 gene_10
# rowData names(0):
# colnames(3): cell_1 cell_2 cell_3
# colData names(0):
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(0):

# The most general method, where we can supply the name of the assay as the second argument.
assay(sce, "counts")
# This is a short-cut for the above, but only works for assays with the special name "counts"
counts(sce)


#################################
######## scRNAseq ###############
#################################
# Here, we use the a droplet-based retina dataset from Macosko et al. (2015), provided 
# in the scRNAseq package. This starts from a count matrix and finishes with clusters 
# (Figure 5.2) in preparation for biological interpretation. Similar workflows are available 
# in abbreviated form in the Workflows.,

# install.packages("scRNAseq")
# library(scRNAseq)
# # Warning in install.packages :
# package ‘scRNAseq’ is not available (for R version 3.6.2)
# sce <- MacoskoRetinaData()
library(Seurat)
sce <- Read10X("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/normoxia_1_5_merge/1_5_merge_monocle3/outs/filtered_feature_bc_matrix")

# Quality control.
library(scater)
# This is pig data "^MT" won't work here, but just running it to test the pipeline
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)
# Error in (function (classes, fdef, mtable)  : 
#             unable to find an inherited method for function ‘logNormCounts’ for signature ‘"dgCMatrix"’
          
### stopped after this ######


# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

# Visualization.
plotUMAP(sce, colour_by="label")