install.packages("Seurat")
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


setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/Joseph_final_results/Jospeh_DEG")
combined_rpkm <- read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/joseph_combined.csv", 
                          stringsAsFactors = FALSE, header=TRUE, row.names = 1, check.names = FALSE)
str(combined_rpkm [,1:20])

# metadata <- readRDS(file = "/Users/pbanerjee/Downloads/pancreas_v3_files/pancreas_metadata.rds")
joseph_metadata <- read.csv(file = "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/joseph_metatdata.csv", 
                            header = TRUE, row.names = 1, check.names = FALSE)
joseph_metadata_df <- data.frame(joseph_metadata)

comb <- CreateSeuratObject(counts = combined_rpkm, 
                           project = "Joseph_Comb", min.cells = 3, 
                           meta.data = joseph_metadata_df)


head(comb@meta.data)

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

# Annotation
levels(immune.combined.clusters)
# [1] "0" "1" "2" "3" "4"
# Assigning cell type identity to clusters
new.cluster.ids <- c("Oligodendroglial genes", "GABA Receptor expressing Neurons", 
                     "Migrating Neurons", "Extracellular Matrix and Signaling components", 
                     "Metabolic Transporters and channels")
names(new.cluster.ids) <- levels(immune.combined.clusters)
immune.combined.clusters.anno <- RenameIdents(immune.combined.clusters, new.cluster.ids)
DimPlot(immune.combined.clusters.anno, reduction = "umap", split.by = "exp")

#DimPlot(immune.combined.clusters.anno, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# No of cells per cluster, most info is collectected in meta.data
#table(immune.combined.clusters.anno@active.ident, immune.combined.clusters.anno@meta.data$orig.ident)

#Identify differential expressed genes across conditions
oligos <- subset(immune.combined.clusters.anno, idents = "Oligodendroglial genes")
Idents(oligos) <- "exp"
avg.oligos <- log1p(AverageExpression(oligos, verbose = FALSE)$RNA)
head(avg.oligos)
avg.oligos$gene <- rownames(avg.oligos)
p1 <- ggplot(avg.oligos, aes(control, treatment)) + geom_point() + ggtitle("Oligodendroglial genes")
# SOX10 is ENSMUSG00000033006
# need to annotate the ENSEMBL IDs
genes.to.label = c("ENSMUSG00000033006")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
plot(p1)

# Make the dataet ready by adding custom column to accomodate celltype and experimnet information
immune.combined.clusters.anno$celltype.exp <- paste(Idents(immune.combined.clusters.anno), immune.combined.clusters.anno$exp, sep = "_")
immune.combined.clusters.anno$celltype <- Idents(immune.combined.clusters.anno)
Idents(immune.combined.clusters.anno) <- "celltype.exp"

# Oligodendroglial genes
deg_oligos <- FindMarkers(immune.combined.clusters.anno, ident.1 = "Oligodendroglial genes_treatment", ident.2 = "Oligodendroglial genes_control", verbose = FALSE)
head(deg_oligos, n = 15)

# GABA Receptor expressing Neurons
deg_gaba <- FindMarkers(immune.combined.clusters.anno, ident.1 = "GABA Receptor expressing Neurons_treatment", ident.2 = "GABA Receptor expressing Neurons_control", verbose = FALSE)
head(deg_gaba, n = 15)

# Migrating Neurons
deg_migrating_neurons <- FindMarkers(immune.combined.clusters.anno, ident.1 = "Migrating Neurons_treatment", ident.2 = "Migrating Neurons_control", verbose = FALSE)
head(deg_migrating_neurons, n = 15)

# Extracellular Matrix and Signaling components
deg_em <- FindMarkers(immune.combined.clusters.anno, ident.1 = "Extracellular Matrix and Signaling components_treatment", ident.2 = "Extracellular Matrix and Signaling components_control", verbose = FALSE)
head(deg_em, n = 15)

# Metabolic Transporters and channels
deg_metabolic <- FindMarkers(immune.combined.clusters.anno, ident.1 = "Metabolic Transporters and channels_treatment", ident.2 = "Metabolic Transporters and channels_control", verbose = FALSE)
head(deg_metabolic, n = 15)

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
b <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol", "description" ),
           values=rownames(deg_metabolic), mart= ensembl)
deg_metabolic$ensembl_gene_id <- rownames(deg_metabolic)
m <- merge(deg_metabolic, b, by="ensembl_gene_id")
names(m)[names(m)=="pct.1"] <- "pct.1(treatment)"
names(m)[names(m)=="pct.2"] <- "pct.2(control)"
write.csv(as.data.frame(m),file = "deg_metabolic.csv")

# # Since we will do differential expression and gene symbols are more human readable 
# # than Ensembl gene IDs, we will get the corresponding gene symbols from Ensembl.
# # library(BUSpaRse)
# gns <- tr2g_ensembl(species = "Mus musculus", use_gene_name = TRUE, 
#                     ensembl_version = 97)[,c("gene", "gene_name")] %>% 
#   distinct()

FeaturePlot(immune.combined.clusters.anno, features = c("SOX10"), split.by = "exp", max.cutoff = 3, 
            cols = c("grey", "red"))


saveRDS(immune.combined.clusters.anno, file = "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/Joseph_final_results/joseph_final.rds")
# sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.4
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_3.1.5
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.1          tidyr_1.0.2         jsonlite_1.6.1      viridisLite_0.3.0  
# [5] splines_3.6.2       lsei_1.2-0          leiden_0.3.3        gtools_3.8.2       
# [9] assertthat_0.2.1    yaml_2.2.1          ggrepel_0.8.2       sessioninfo_1.1.1  
# [13] globals_0.12.5      pillar_1.4.3        lattice_0.20-41     glue_1.4.0         
# [17] reticulate_1.15     digest_0.6.25       RColorBrewer_1.1-2  colorspace_1.4-1   
# [21] cowplot_1.0.0       htmltools_0.4.0     Matrix_1.2-18       plyr_1.8.6         
# [25] pkgconfig_2.0.3     tsne_0.1-3          listenv_0.8.0       purrr_0.3.4        
# [29] patchwork_1.0.0     scales_1.1.0        RANN_2.6.1          gdata_2.18.0       
# [33] Rtsne_0.15          tibble_3.0.1        ggplot2_3.3.0       ellipsis_0.3.0     
# [37] withr_2.2.0         ROCR_1.0-7          pbapply_1.4-2       lazyeval_0.2.2     
# [41] cli_2.0.2           survival_3.1-12     magrittr_1.5        crayon_1.3.4       
# [45] fansi_0.4.1         future_1.17.0       nlme_3.1-147        MASS_7.3-51.6      
# [49] gplots_3.0.3        ica_1.0-2           tools_3.6.2         fitdistrplus_1.0-14
# [53] data.table_1.12.8   lifecycle_0.2.0     stringr_1.4.0       plotly_4.9.2.1     
# [57] munsell_0.5.0       cluster_2.1.0       irlba_2.3.3         compiler_3.6.2     
# [61] rsvd_1.0.3          caTools_1.18.0      rlang_0.4.5         grid_3.6.2         
# [65] ggridges_0.5.2      rstudioapi_0.11     RcppAnnoy_0.0.16    htmlwidgets_1.5.1  
# [69] igraph_1.2.5        bitops_1.0-6        npsurv_0.4-0        gtable_0.3.0       
# [73] codetools_0.2-16    reshape2_1.4.4      R6_2.4.1            gridExtra_2.3      
# [77] zoo_1.8-7           dplyr_0.8.5         uwot_0.1.8          future.apply_1.5.0 
# [81] KernSmooth_2.23-17  ape_5.3             stringi_1.4.6       parallel_3.6.2     
# [85] Rcpp_1.0.4.6        vctrs_0.2.4         sctransform_0.2.1   png_0.1-7          
# [89] tidyselect_1.0.0    lmtest_0.9-37      
