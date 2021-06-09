################################################
################################################
### SAVE AND LOAD WORKSPACE 
################################################
################################################

setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/")
save.image(file='Camille_pseudotime.RData')
load('Camille_pseudotime.RData')

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


##########################################################################################################
##########################################################################################################
##########################################################################################################
# Test Embryo Data
##########################################################################################################
##########################################################################################################
##########################################################################################################
# expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
# cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
# gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))
# 
# test_embryo_cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)
# test_embryo_cds_preproc <- preprocess_cds(test_embryo_cds, num_dim = 50, preprocess_method = "PCA")
# test_embryo_cds_align <- align_cds(test_embryo_cds_preproc, alignment_group = "batch", 
#                  residual_model_formula_str = "~ bg.300.loading + bg.400.loading + 
#                  bg.500.1.loading + bg.500.2.loading + bg.r17.loading + 
#                  bg.b01.loading + bg.b02.loading")
# test_embryo_cds_redim <- reduce_dimension(test_embryo_cds_align,  reduction_method = "UMAP", preprocess_method = "PCA")
# plot_cells(test_embryo_cds_redim, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# plot_cells(test_embryo_cds_redim,
#            color_cells_by = "cell.type",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)

##########################################################################################################
##########################################################################################################
##########################################################################################################
# Normox
##########################################################################################################
##########################################################################################################
##########################################################################################################
# Provide the path to the Cell Ranger output. 
# The directory should have the outs folder from cellranger 
cds_nmx_monocle3 <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/normoxia_1_5_merge/1_5_merge_monocle3")

cds_nmx_m_preproc <- preprocess_cds(cds_nmx_monocle3, num_dim = 50)
plot_pc_variance_explained(cds_nmx_m_preproc)
cds_nmx_red_dim <- reduce_dimension(cds_nmx_m_preproc, 
                                    reduction_method = 'UMAP', preprocess_method = 'PCA') 
plot_cells(cds_nmx_red_dim)
# cluster_cells uses a technique called community detection to group cells. 
# This approach was introduced by Levine et al as part of the phenoGraph algorithm.
# cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
# using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm.
cds_nmx_cluster <- cluster_cells(cds_nmx_red_dim)
plot_cells(cds_nmx_cluster)
#### Annoate the cell clusters
# Add a new column to colData(cds_ordered) and initialize it with the values of clusters(cds)
colData(cds_nmx_cluster)$assigned_cell_type <- as.character(partitions(cds_nmx_cluster))
# # Dim 30 
# colData(cds_nmx_cluster)$assigned_cell_type = dplyr::recode(colData(cds_nmx_cluster)$assigned_cell_type,
#                                                         "1"="OL1",
#                                                         "2"="Med Spiny Neurons",
#                                                         "3"="Med Spiny Neurons",
#                                                         "4"="OL2",
#                                                         "5"="Immune Cells",
#                                                         "6"="Astrocytes",
#                                                         "7"="NSPCs",
#                                                         "8"="post-NSPCs",
#                                                         "9"="post-NSPCs",
#                                                         "10"="OPCs",
#                                                         "11"="Unclassified neurons",
#                                                         "12"="OPCs",
#                                                         "13"="OL1",
#                                                         "14"="post-NSPCs",
#                                                         "15"="Immune Cells",
#                                                         "16"="Ependymal Cells",
#                                                         "17"="Unclassified neurons",
#                                                         "18"="Neurons (IN)",
#                                                         "19"="Unclassified neurons",
#                                                         "20"="post-NSPCs",
#                                                         "21"="Med Spiny Neurons",
#                                                         "22"="Vessels",
#                                                         "23"="Vessels",
#                                                         "24"="Astrocytes",
#                                                         "25"="Immune Cells",
#                                                         "26"="Neurons (IN)",
#                                                         "27"="Immune Cells",
#                                                         "28"="Immune Cells",
#                                                         "29"="Unclassified neurons",
#                                                         "30"="Unclassified neurons",
#                                                         "31"="Neurons (IN)",
#                                                         "32"="OL1",
#                                                         "33"="OL1")
## Dim 50
#Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
colData(cds_nmx_cluster)$assigned_cell_type = dplyr::recode(colData(cds_nmx_cluster)$assigned_cell_type,
                                                            "1"="OL1",
                                                            "2"="Med Spiny Neurons",
                                                            "3"="Med Spiny Neurons",
                                                            "4"="OL2",
                                                            "5"="Immune Cells",
                                                            "6"="post-NSPCs",
                                                            "7"="NSPCs",
                                                            "8"="Astrocytes/NPCs",
                                                            "9"="OPCs",
                                                            "10"="post-NSPCs",
                                                            "11"="OPCs",
                                                            "12"="OL1",
                                                            "13"="Immune Cells",
                                                            "14"="NPCs",
                                                            "15"="OL1",
                                                            "16"="Ependymal Cells",
                                                            "17"="Med Spiny Neurons",
                                                            "18"="Unclassified neurons",
                                                            "19"="post-NSPCs",
                                                            "20"="Neurons (IN)",
                                                            "21"="Med Spiny Neurons",
                                                            "22"="Med Spiny Neurons",
                                                            "23"="Vessels",
                                                            "24"="Vessels",
                                                            "25"="Neurons (IN)",
                                                            "26"="Vessels",
                                                            "27"="Neurons (IN)",
                                                            "28"="Immune Cells",
                                                            "29"="OL1",
                                                            "30"="Unclassified neurons")

plot_cells(cds_nmx_cluster, group_cells_by="partition", 
                         color_cells_by="assigned_cell_type",           
                         show_trajectory_graph = FALSE)

# find all possible partitions
all_partitions <- unique(cds_nmx_cluster@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]

# set all partitions to 1
cds_nmx_cluster@clusters$UMAP$partitions[cds_nmx_cluster@clusters$UMAP$partitions %in% all_partitions] <- "1"

# # with partition
# cds_nmx_graph <- learn_graph(cds_nmx_cluster, use_partition = TRUE, close_loop = TRUE)
# plot_cells(cds_nmx_graph,show_trajectory_graph = TRUE)
# plot_cells(cds_nmx_graph,show_trajectory_graph = FALSE)
# plot_cells(cds_nmx_graph, group_cells_by="partition", 
#            color_cells_by="assigned_cell_type",           
#            show_trajectory_graph = TRUE)
# without partition
# root_group = colnames(cds_nmx_cluster)[clusters(cds_nmx_cluster) == 1]
cds_nmx_graph_without_partition <- learn_graph(cds_nmx_cluster, use_partition = FALSE, close_loop = TRUE)
plot_cells(cds_nmx_graph_without_partition,show_trajectory_graph = TRUE)
plot_cells(cds_nmx_graph_without_partition,show_trajectory_graph = FALSE)
plot_cells(cds_nmx_graph_without_partition, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = FALSE)
# # Since, with partition, I select root nodes for each cluster
# #root_group = colnames(cds_nmx_cluster)[clusters(cds_nmx_cluster) == 1]
# #cds_ordered = order_cells(cds_nmx_graph, root_cells = root_group)
# cds_nmx_ordered <- order_cells(cds_nmx_graph, root_pr_nodes= )
# plot_cells(cds_nmx_ordered,show_trajectory_graph = TRUE)
# plot_cells(cds_nmx_ordered,show_trajectory_graph = FALSE)
# plot_cells(cds_nmx_ordered,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5)
# Since without partition, Chose cluster 7 as the root, since they are the 
# Neural STEM progenitor cells (NSPCs), according to Camille's annotation
cds_nmx_without_partition_ordered <- order_cells(cds_nmx_graph_without_partition)
plot_cells(cds_nmx_without_partition_ordered,show_trajectory_graph = TRUE)
plot_cells(cds_nmx_without_partition_ordered,show_trajectory_graph = FALSE)
plot_cells(cds_nmx_without_partition_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


interesting_genes <- c("DDX6",
                       "DDX5",
                       "PDE10A",
                       "ADCY5")
plot_cells(cds_nmx_cluster,
           genes=interesting_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


# 3D trajectories
cds_nmx_red_dim_3d <- reduce_dimension(cds_nmx_m_preproc, max_components = 3,preprocess_method = 'PCA',reduction_method = 'UMAP')
cds_nmx_cluster_3d <- cluster_cells(cds_nmx_red_dim_3d)
cds_nmx_graph_3d <- learn_graph(cds_nmx_cluster_3d)
plot_cells_3d(cds_nmx_graph_3d, color_cells_by="partition", show_trajectory_graph = FALSE, reduction_method = "UMAP")
cds_nmx_without_partition_3d <- learn_graph(cds_nmx_graph_3d, use_partition = FALSE, close_loop = TRUE)
plot_cells_3d(cds_nmx_without_partition_3d, color_cells_by="partition", show_trajectory_graph =FALSE)
# TO DO - have to ascertain root nodes
cds_nmx_ordered_3d <- order_cells(cds_nmx_graph_3d , root_pr_nodes=get_earliest_principal_node(cds_nmx_cluster_3d ))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")


########################################
#Choose subset
########################################
cds_nmx_cluster <- cluster_cells(cds_nmx_red_dim)
colData(cds_nmx_cluster)$assigned_cell_type <- as.character(partitions(cds_nmx_cluster))
#Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
colData(cds_nmx_cluster)$assigned_cell_type = dplyr::recode(colData(cds_nmx_cluster)$assigned_cell_type,
                                                            "1"="OL1",
                                                            "2"="Med Spiny Neurons",
                                                            "3"="Med Spiny Neurons",
                                                            "4"="OL2",
                                                            "5"="Immune Cells",
                                                            "6"="post-NSPCs",
                                                            "7"="NSPCs",
                                                            "8"="Astrocytes/NPCs",
                                                            "9"="OPCs",
                                                            "10"="post-NSPCs",
                                                            "11"="OPCs",
                                                            "12"="OL1",
                                                            "13"="Immune Cells",
                                                            "14"="NPCs",
                                                            "15"="OL1",
                                                            "16"="Ependymal Cells",
                                                            "17"="Med Spiny Neurons",
                                                            "18"="Unclassified neurons",
                                                            "19"="post-NSPCs",
                                                            "20"="Neurons (IN)",
                                                            "21"="Med Spiny Neurons",
                                                            "22"="Med Spiny Neurons",
                                                            "23"="Vessels",
                                                            "24"="Vessels",
                                                            "25"="Neurons (IN)",
                                                            "26"="Vessels",
                                                            "27"="Neurons (IN)",
                                                            "28"="Immune Cells",
                                                            "29"="OL1",
                                                            "30"="Unclassified neurons")
plot_cells(cds_nmx_cluster, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = FALSE)
cds_normox_subset <- choose_cells(cds_nmx_cluster)
# pr_graph_test_res <- graph_test(cds_normox_subset, neighbor_graph="knn", cores=8)
# pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
# cds_subset_cluster <- cluster_cells(cds_normox_subset, resolution=1e-2)
# plot_cells(cds_subset_cluster, color_cells_by="cluster")
# find all possible partitions
all_partitions_subset <- unique(cds_normox_subset@clusters$UMAP$partitions)
all_partitions_subset <- all_partitions_subset[all_partitions_subset != "1"]
# set all partitions to 1
cds_normox_subset@clusters$UMAP$partitions[cds_normox_subset@clusters$UMAP$partitions %in% all_partitions_subset] <- "1"
# cds_subset_cluster <- cluster_cells(cds_normox_subset, resolution=1e-2)
# plot_cells(cds_subset_cluster, color_cells_by="cluster")
cds_nmx_subcluster_graph <- learn_graph(cds_normox_subset, use_partition = FALSE, close_loop = TRUE)
plot_cells(cds_nmx_subcluster_graph,show_trajectory_graph = TRUE)
cds_nmx_subcluster_graph_ordered <- order_cells(cds_nmx_subcluster_graph)
plot_cells(cds_nmx_subcluster_graph_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



##########################################################################################################
setwd("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/hypoxia_2_6_8_merge/")
##########################################################################################################
##########################################################################################################
##########################################################################################################
# Hypox
##########################################################################################################
##########################################################################################################
##########################################################################################################
# Provide the path to the Cell Ranger output. 
# The directory should have the outs folder from cellranger 
cds_hpx_monocle3 <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/hypoxia_2_6_8_merge")
cds_hpx_preproc <- preprocess_cds(cds_hpx_monocle3, num_dim = 50)
plot_pc_variance_explained(cds_hpx_preproc)
cds_hpx_red_dim <- reduce_dimension(cds_hpx_preproc, 
                                    reduction_method = 'UMAP', preprocess_method = 'PCA') 
plot_cells(cds_nmx_red_dim)
# cluster_cells uses a technique called community detection to group cells. 
# This approach was introduced by Levine et al as part of the phenoGraph algorithm.
# cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
# using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm.
cds_hpx_cluster <- cluster_cells(cds_hpx_red_dim)
# with partition
cds_hpx_graph <- learn_graph(cds_hpx_cluster, use_partition = TRUE, close_loop = TRUE)
plot_cells(cds_hpx_graph,show_trajectory_graph = TRUE)
plot_cells(cds_hpx_graph,show_trajectory_graph = FALSE)
# After annotation
plot_cells(cds_hpx_graph, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = TRUE)
# without partition
cds_hpx_graph_without_partition <- learn_graph(cds_hpx_cluster, use_partition = FALSE, close_loop = TRUE)
plot_cells(cds_hpx_graph_without_partition,show_trajectory_graph = TRUE)
plot_cells(cds_hpx_graph_without_partition,show_trajectory_graph = FALSE)
plot_cells(cds_hpx_graph_without_partition, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = TRUE)
# Since, with partition, I select root nodes for each cluster
cds_hpx_ordered <- order_cells(cds_hpx_graph)
plot_cells(cds_hpx_ordered,show_trajectory_graph = TRUE)
plot_cells(cds_hpx_ordered,show_trajectory_graph = FALSE)
plot_cells(cds_hpx_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# Since without partition, Chose cluster 7 as the root, since they are the 
# Neural progenitor cells (NPCs), according to Camille's annotation
cds_nmx_without_partition_ordered <- order_cells(cds_nmx_graph_without_partition)
plot_cells(cds_nmx_without_partition_ordered,show_trajectory_graph = TRUE)
plot_cells(cds_nmx_without_partition_ordered,show_trajectory_graph = FALSE)
plot_cells(cds_nmx_without_partition_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


top_hpx_markers_cluster <- top_markers(cds_hpx_graph, group_cells_by = "cluster",
                                   genes_to_test_per_group = 25, reduction_method = "UMAP",
                                   marker_sig_test = TRUE, reference_cells = NULL, cores = 6,
                                   verbose = FALSE)
write.table(top_hpx_markers_cluster, file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/hypoxia_2_6_8_merge/top_hpx_marker_cluster.txt", sep = "\t")

# Rank markers according to pseudo_R2, groupby cluster
top_specific_hpx_markers_cluster <- top_hpx_markers_cluster %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
write.table(top_specific_hpx_markers_cluster, file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/hypoxia_2_6_8_merge/top_hpx_specific_marker_cluster.txt", sep = "\t")
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


#### Annoate the cell clusters
# Add a new column to colData(cds_ordered) and initialize it with the values of clusters(cds)
colData(cds_nmx_cluster)$assigned_cell_type <- as.character(partitions(cds_nmx_cluster))
## Dim 50
#Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
colData(cds_nmx_cluster)$assigned_cell_type = dplyr::recode(colData(cds_nmx_cluster)$assigned_cell_type,
                                                            "1"="OL1",
                                                            "2"="Med Spiny Neurons",
                                                            "3"="Med Spiny Neurons",
                                                            "4"="OL2",
                                                            "5"="Immune Cells",
                                                            "6"="post-NSPCs",
                                                            "7"="NSPCs",
                                                            "8"="Astrocytes/NPCs",
                                                            "9"="OPCs",
                                                            "10"="post-NSPCs",
                                                            "11"="OPCs",
                                                            "12"="OL1",
                                                            "13"="Immune Cells",
                                                            "14"="NPCs",
                                                            "15"="OL1",
                                                            "16"="Ependymal Cells",
                                                            "17"="Med Spiny Neurons",
                                                            "18"="Unclassified neurons",
                                                            "19"="post-NSPCs",
                                                            "20"="Neurons (IN)",
                                                            "21"="Med Spiny Neurons",
                                                            "22"="Med Spiny Neurons",
                                                            "23"="Vessels",
                                                            "24"="Vessels",
                                                            "25"="Neurons (IN)",
                                                            "26"="Vessels",
                                                            "27"="Neurons (IN)",
                                                            "28"="Immune Cells",
                                                            "29"="OL1",
                                                            "30"="Unclassified neurons")                           


interesting_genes <- c("DDX6",
                       "DDX5",
                       "PDE10A",
                       "ADCY5")
plot_cells(cds_hpx_cluster,
           genes=interesting_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# 3D trajectories
cds_hpx_red_dim_3d <- reduce_dimension(cds_hpx_preproc, max_components = 3,preprocess_method = 'PCA',reduction_method = 'UMAP')
cds_hpx_cluster_3d <- cluster_cells(cds_hpx_red_dim_3d)
cds_hpx_graph_3d <- learn_graph(cds_hpx_cluster_3d)
plot_cells_3d(cds_hpx_graph_3d, color_cells_by="partition", reduction_method = "UMAP")
cds_hpx_without_partition_3d <- learn_graph(cds_hpx_graph_3d, use_partition = FALSE, close_loop = TRUE)
plot_cells_3d(cds_hpx_without_partition_3d, color_cells_by="partition", show_trajectory_graph =FALSE)
# TO DO - have to ascertain root nodes
cds_nmx_ordered_3d <- order_cells(cds_nmx_graph_3d , root_pr_nodes=get_earliest_principal_node(cds_nmx_cluster_3d ))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

################
# Dimension 10 clusters
top_markers_cluster_dim10 <- top_markers(cds_nmx_graph_without_partition, group_cells_by = "cluster",
                                   genes_to_test_per_group = 25, reduction_method = "UMAP",
                                   marker_sig_test = TRUE, reference_cells = NULL, cores = 6,
                                   verbose = FALSE)
write.table(top_markers_cluster_dim10, file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/normoxia_1_5_merge/top_marker_cluster_dim10.txt", sep = ",")

#### Annoate the cell clusters
# Add a new column to colData(cds_ordered) and initialize it with the values of clusters(cds)
colData(cds_nmx_cluster)$assigned_cell_type <- as.character(partitions(cds_nmx_cluster))
## Dim 10
#Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
colData(cds_nmx_cluster)$assigned_cell_type = dplyr::recode(colData(cds_nmx_cluster)$assigned_cell_type,
                                                            "2"="OL1",
                                                            "3"="post-NSPCs",
                                                            "4"="Ependymal Cells",
                                                            "5"="Differentiating Oligodendrocytes/Immune Cells",
                                                            "6"="OL2",
                                                            "8"="Astrocytes/NPCs",
                                                            "10"="OL1",
                                                            "11"="OPCs",
                                                            "12"="Neurons",
                                                            "14"="NSPCs",
                                                            "19"="OL1",
                                                            "20"="Mature Oligos")


cds_nmx_graph_without_partition <- learn_graph(cds_nmx_cluster, use_partition = FALSE, close_loop = TRUE)
plot_cells(cds_nmx_graph_without_partition, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = FALSE)
plot_cells(cds_nmx_graph_without_partition, show_trajectory_graph = TRUE)
cds_nmx_without_partition_ordered <- order_cells(cds_nmx_graph_without_partition)
plot_cells(cds_nmx_without_partition_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# find all possible partitions
all_partitions_subset <- unique(cds_subset_cluster@clusters$UMAP$partitions)
all_partitions_subset <- all_partitions_subset[all_partitions_subset != "1"]
# set all partitions to 1
cds_subset_cluster@clusters$UMAP$partitions[cds_subset_cluster@clusters$UMAP$partitions %in% all_partitions_subset] <- "1"

#Choose subset
########################################
cds_normox_subset <- choose_cells(cds_nmx_cluster)
# pr_graph_test_res <- graph_test(cds_normox_subset, neighbor_graph="knn", cores=8)
# pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
cds_subset_cluster <- cluster_cells(cds_normox_subset, resolution=1e-2)
plot_cells(cds_subset_cluster, color_cells_by="cluster" )
# find all possible partitions
all_partitions_subset <- unique(cds_subset_cluster@clusters$UMAP$partitions)
all_partitions_subset <- all_partitions_subset[all_partitions_subset != "1"]
# set all partitions to 1
cds_subset_cluster@clusters$UMAP$partitions[cds_subset_cluster@clusters$UMAP$partitions %in% all_partitions_subset] <- "1"
cds_subset_cluster <- cluster_cells(cds_normox_subset, resolution=1e-2)
plot_cells(cds_subset_cluster, color_cells_by="cluster")
cds_nmx_subcluster_graph <- learn_graph(cds_subset_cluster, use_partition = FALSE, close_loop = TRUE)
plot_cells(cds_nmx_subcluster_graph,show_trajectory_graph = TRUE)
cds_nmx_subcluster_graph_ordered <- order_cells(cds_nmx_subcluster_graph)
plot_cells(cds_nmx_subcluster_graph_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


