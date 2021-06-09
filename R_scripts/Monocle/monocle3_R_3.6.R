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

# 
# library(Matrix)
# matrix_dir = " /Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix"
# barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
# features.path <- paste0(matrix_dir, "/features.tsv.gz")
# matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
# mat <- readMM(file = matrix.path)
# feature.names = read.delim("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix/features.tsv.gz", 
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# # Assign gene_short_name to the feature.names datframe 
# # as thats required for monocle3 gene_metadata
# names(feature.names)[2] <- "gene_short_name"
# # names(feature.names)[1] <- "ensemble_id"
# # Assign the first ensembl id column as the rowname for the dataframe
# rownames(feature.names) <- feature.names[,1]
# feature.names[,1] <- NULL
# barcode.names = read.delim("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix/barcodes.tsv.gz", 
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# 
# # Provide the path to the Cell Ranger output. 
# # The directory should have the outs folder from cellranger 
# cds_10x <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3")
# # check out the expression data using exprs()
# # 25880 genes across 17877 single cells for Camille 1_5_merge dataset
# # dim(cds)
# # dim(exprs(cds))
# # checking out first five genes in the first five cells
# # exprs(cds)[1:5, 1:5]
# # Extract cell metadata from the expression matrix
# # cell_metadata <- pData(cds)
# # gene_metadata <- fData(cds)
# # Make the new cell dataset with the cell_metadata and the gene_metadata
# 
# # cds <- new_cell_data_set(cds,cell_metadata = pData(cds),
# #                                gene_metadata = fData(cds))
# 
# 
# # cds <- new_cell_data_set(mat,cell_metadata = barcode.names,
# #                          gene_metadata = feature.names)
# 
# # Your CDS object contains cells with zero reads. 
# # This causes size factor calculation to fail. 
# # Please remove the zero read cells using 
# # cds <- cds[,Matrix::rowSums(exprs(cds)) != 0] and then run 
# # cds <- estimate_size_factors(cds)
# 
# cds_test <- cds[,Matrix::colSums(exprs(cds)) != 0]
# cds_ult <- estimate_size_factors(cds_test)
# 
# 
# # Step 1: Normalize and pre-process the data
# cds_final <- preprocess_cds(cds_ult, num_dim = 30)
# 
# # Step 2: Reduce the dimensionality of the data 
# # Then, apply a further round of (nonlinear) dimensionality reduction using UMAP:
# cds_dim <- reduce_dimension(cds_final, reduction_method = 'UMAP')
# 
# # Step 3: Cluster the cells
# cds_cluster <- cluster_cells(cds_dim)
# 
# # Step 4: Learn a graph
# cds <- learn_graph(cds)
# 
# # Step 5: Order cells
# cds <- order_cells(cds)
# 
# plot_cells(cds)

# test_cell_metadata <- cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
# test_gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))


# Provide the path to the Cell Ranger output. 
# The directory should have the outs folder from cellranger 
cds_nmx_monocle3 <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3")

cds_m_preproc <- preprocess_cds(cds_monocle3, num_dim = 100)
plot_pc_variance_explained(cds_m_preproc)
cds_red_dim <- reduce_dimension(cds_m_preproc, reduction_method = 'UMAP', preprocess_method = 'PCA') #preprocess_method = 'PCA'
plot_cells(cds_red_dim)
# cluster_cells uses a technique called community detection to group cells. 
# This approach was introduced by Levine et al as part of the phenoGraph algorithm.
# cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
# using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm.
cds_cluster <- cluster_cells(cds_red_dim)
cds_graph <- learn_graph(cds_cluster)
plot_cells(cds_graph, show_trajectory_graph = FALSE)
# Chose cluster 7 as the root, since they are the 
# Neural progenitor cells (NPCs), according to Camille's annotation
cds_ordered <- order_cells(cds_graph)
plot_cells(cds_graph,show_trajectory_graph = TRUE)
plot_cells(cds_graph,show_trajectory_graph = FALSE)

# Once cells have been clustered, 
# what genes makes them different from one another.
# Find top markers by cluster
top_markers_cluster <- top_markers(cds_graph, group_cells_by = "cluster",
            genes_to_test_per_group = 25, reduction_method = "UMAP",
            marker_sig_test = TRUE, reference_cells = NULL, cores = 6,
            verbose = FALSE)
write.table(top_markers_cluster, file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/top_marker_cluster_dim100.txt", sep = ",")

# Find top markers by partition
top_markers_partition <- top_markers(cds_graph, group_cells_by = "partition",
                                     genes_to_test_per_group = 25, reduction_method = "UMAP",
                                     marker_sig_test = TRUE, reference_cells = NULL, cores = 6,
                                     verbose = FALSE)
write.table(top_markers_partition, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/top_marker_partition_dim100.txt", sep = ",")

# library dyplr for %%
library(dplyr)

# Rank markers according to pseudo_R2, groupby cluster
top_specific_markers <- top_markers_cluster %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
write.table(top_specific_markers, file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/top_specific_marker_cluster_dim100.txt", sep = "\t")
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# Rank markers according to pseudo_R2, groupby partition
top_specific_markers_partition <- top_markers_partition %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
write.table(top_specific_markers_partition, "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/top_specific_marker_partition_dim100.txt", sep = "\t")
top_specific_marker_ids_partition <- unique(top_specific_markers_partition %>% pull(gene_id))

# Plot the expression and fraction of cells that express each marker in each group 
# color of the dot to represent the percentage and size of the dot the mean expression
plot_genes_by_group(cds_graph,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3 )

# color of the dot to represent the percentage and size of the dot the mean expression
plot_genes_by_group(cds_graph,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3 )

# color of the dot to represent the percentage and size of the dot the mean expression
plot_genes_by_group(cds_graph,
                    top_specific_marker_ids_partition,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3 )

# color of the dot to represent the percentage and size of the dot the mean expression
plot_genes_by_group(cds_graph,
                    top_specific_marker_ids_partition,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3 )



## More than one marker
top_specific_markers_more_than_one_marker <- top_specific_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids__than_one_marker <- unique(top_specific_markers_more_than_one_marker %>% pull(gene_id))

plot_genes_by_group(cds_graph,
                    top_specific_marker_ids__than_one_marker,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)



#### PseudoTime #####
#####################

cds_10x <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CNMC_Projects/
                                Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3")

# Preprocess
# Method - PCA
cds_preprocess <- preprocess_cds(cds_10x, num_dim = 50)
plot_pc_variance_explained(cds_preprocess)
# Align
# Method - removes batch effects using mutual nearest neighbor alignment, 
# a technique introduced by John Marioni's lab. 
# Monocle 3 does so by calling Aaron Lun's excellent package batchelor.
cds_align <- align_cds(cds_preprocess)
# Reduce Dimension
# Method - UMAP (default)
cds_reduce_dim<- reduce_dimension(cds_align)
# Cluster the cells
# Method - community detection
# introduced by Levine et al as part of the phenoGraph algorithm.
cds_cluster_again <- cluster_cells(cds_reduce_dim)
plot_cells(cds_cluster_again, color_cells_by = "partition")
# Learn a graph
cds_learn <- learn_graph(cds_cluster_again)
plot_cells(cds_learn,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
# Order cells 
# Have to choose root
cds_ordered <- order_cells(cds_learn, reduction_method = "UMAP", root_pr_nodes = NULL,
                           root_cells = NULL, verbose = FALSE)
# Plot Cells
plot_cells(cds_learn)
plot_cells(cds_learn, label_groups_by_cluster=FALSE)
plot_cells(cds_learn,show_trajectory_graph = FALSE)

# Annoatate the cell clusters
# Add a new column to colData(cds_ordered) and initialize it with the values of clusters(cds)
colData(cds_ordered)$assigned_cell_type <- as.character(partitions(cds_ordered))
# Now, we can use the dplyrpackage's recode() function 
# to remap each cluster to a different cell type
colData(cds_ordered)$assigned_cell_type = dplyr::recode(colData(cds_ordered)$assigned_cell_type,
                                                "1"="OL1",
                                                "2"="Med Spiny Neurons",
                                                "3"="Med Spiny Neurons",
                                                "4"="OL2",
                                                "5"="Immune Cells",
                                                "6"="Astrocytes",
                                                "7"="NSPCs",
                                                "8"="post-NSPCs",
                                                "9"="post-NSPCs",
                                                "10"="OPCs",
                                                "11"="Unclassified neurons",
                                                "12"="OPCs",
                                                "13"="OL1",
                                                "14"="post-NSPCs",
                                                "15"="Immune Cells",
                                                "16"="Ependymal Cells",
                                                "17"="Unclassified neurons",
                                                "18"="Neurons (IN)",
                                                "19"="Unclassified neurons",
                                                "20"="post-NSPCs",
                                                "21"="Med Spiny Neurons",
                                                "22"="Vessels",
                                                "23"="Vessels",
                                                "24"="Astrocytes",
                                                "25"="Immune Cells",
                                                "26"="Neurons (IN)",
                                                "27"="Immune Cells",
                                                "28"="Immune Cells",
                                                "29"="Unclassified neurons",
                                                "30"="Unclassified neurons",
                                                "31"="Neurons (IN)",
                                                "32"="OL1",
                                                "33"="OL1")
# Let's see how the new annotations look:

plot_cells(cds_ordered, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = FALSE)

plot_cells(cds_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
