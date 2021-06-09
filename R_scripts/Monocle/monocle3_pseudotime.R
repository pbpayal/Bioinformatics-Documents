#### Load the following librararies
library("monocle3")
library(devtools)
library(dplyr)


#### Load the data
cds_monocle3 <- load_cellranger_data("../../outs")
#### Pre-process the data
# Normalize the data using Principal Components Analysis since this is RNA seq analysis
cds_m_preproc <- preprocess_cds(cds_monocle3, num_dim = 50)
#### Reduction of Dimensionality
cds_red_dim <- reduce_dimension(cds_m_preproc, reduction_method = 'UMAP', preprocess_method = 'PCA')
#### Cluster the Cells
# cluster_cells uses a technique called community detection to group cells. 
# This approach was introduced by Levine et al as part of the phenoGraph algorithm.
# cluster_cells() also divides the cells into larger, more well separated groups 
# called partitions, using a statistical test from Alex Wolf et al, introduced as part of 
# their PAGA algorithm.
cds_cluster <- cluster_cells(cds_red_dim)
# cds_partition <- partitions(cds_cluster)
#### Learn the trajectory graph
cds_graph <- learn_graph(cds_cluster, use_partition = FALSE, close_loop = TRUE)
cds_graph_with_partition <- learn_graph(cds_cluster, use_partition = TRUE, close_loop = TRUE)
plot_cells(cds_graph)
# Plot the cells from the learn graph to see the clusters, with and without the trajectory graph to detect the root cluster/node
plot_cells(cds_graph, show_trajectory_graph = FALSE)
# In this case,I have preliminary annotations according to cluster number and thus selected ```show_trajectory_graph = FALSE```
cds_ordered <- order_cells(cds_graph)
cds_ordered_with_partition <- order_cells(cds_graph_with_partition)

plot_cells(cds_ordered,show_trajectory_graph = TRUE)
#### Annoate the cell clusters
# Add a new column to colData(cds_ordered) and initialize it with the values of clusters(cds)
colData(cds_cluster)$assigned_cell_type <- as.character(partitions(cds_cluster))
## Dim 50
#Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
colData(cds_cluster)$assigned_cell_type = dplyr::recode(colData(cds_cluster)$assigned_cell_type,
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


# ## Dim 100
# #Now, we can use the dplyrpackage's recode() function to remap each cluster to a different cell type
# colData(cds_cluster)$assigned_cell_type = dplyr::recode(colData(cds_cluster)$assigned_cell_type,
#                                                         "1"="Unclassified neurons",
#                                                         "2"="Unclassified neurons",
#                                                         "3"="Unclassified neurons",
#                                                         "4"="Unclassified neurons",
#                                                         "5"="Unclassified neurons",
#                                                         "6"="post-NSPCs",
#                                                         "7"="post-NSPCs",
#                                                         "8"="Unclassified neurons",
#                                                         "9"="OPCs",
#                                                         "10"="OPCs",
#                                                         "11"="OPCs",
#                                                         "12"="OL1",
#                                                         "13"="Immune Cells",
#                                                         "14"="OL1",
#                                                         "15"="Ependymal Cells",
#                                                         "16"="Astrocytes",
#                                                         "17"="Neurons (IN)",
#                                                         "18"="Med Spiny Neurons",
#                                                         "19"="Med Spiny Neurons",
#                                                         "20"="Vessels",
#                                                         "21"="Vessels",
#                                                         "22"="OL2",
#                                                         "23"="Vessels",
#                                                         "24"="Immune cells/Astrocytes",
#                                                         "25"="Neurons (IN)",
#                                                         "26"="Unclassified neurons",
#                                                         "27"="OL1",
#                                                         "28"="Immune Cells",
#                                                         "29"="OPCs",
#                                                         "30"="Unclassified neurons")

# Let's see how the new annotations look:
plot_cells(cds_cluster, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = FALSE)

plot_cells(cds_graph, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = TRUE)

plot_cells(cds_graph_with_partition, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = TRUE)

#### Learn the trajectory graph
cds_graph_after_annotation <- learn_graph(cds_cluster, use_partition = FALSE)
cds_ordered_after_annotation <- order_cells(cds_graph_after_annotation)
## I chose multiple roots for ech cluster and multiple roots for subclusters for testing
plot_cells(cds_ordered_after_annotation, group_cells_by="partition", 
           color_cells_by="assigned_cell_type",           
           show_trajectory_graph = TRUE)

#### Plotting the cells and coloring them by pseudotime shows how they were ordered
plot_cells(cds_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds_ordered_with_partition,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

differentially_expr_genes_nmx <- graph_test(cds_ordered, neighbor_graph="knn", cores=8)
nmx_deg_ids <- row.names(subset(differentially_expr_genes_nmx, q_value < 0.05))

deg_nmx_pseudotime <- graph_test(cds_ordered,neighbor_graph="principal_graph",cores=4 )
deg_nmx_pseudotime_ids <- row.names(subset(deg_nmx_pseudotime, q_value < 0.05))
# As before, we can collect the trajectory-variable genes into modules:
gene_module_df_nmx <- find_gene_modules(cds_ordered[deg_nmx_pseudotime_ids,], resolution=c(0,10^seq(-6,-1)))
plot_cells(cds_ordered,
           genes=gene_module_df_nmx %>% filter(module %in% c(27, 10, 7, 30)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# set of genes that vary in some interesting way across the clusters
gene_module_df_nmx <- find_gene_modules(cds_ordered[nmx_deg_ids,], resolution=1e-2)
plot_cells(cds_ordered, genes=c("PSMB1", "ENSSSCG0000003737", "FAM120B", "ENSSSCG00000027274","DLL1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
## Choose a subset of cells

####################################################################################
cds_NMX_vs_HPX <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/normoxia_vs_hypoxia/")
cds_NMX_vs_HPX_preproc <- preprocess_cds(cds_NMX_vs_HPX, num_dim = 50)
cds_NMX_vs_HPX_red_dim <- reduce_dimension(cds_NMX_vs_HPX_preproc, reduction_method = 'UMAP', preprocess_method = 'PCA')
cds_NMX_vs_HPX_cluster <- cluster_cells(cds_NMX_vs_HPX_red_dim)
plot_cells(cds_NMX_vs_HPX_cluster,show_trajectory_graph = FALSE)

cds_NMX_vs_HPX_graph <- learn_graph(cds_NMX_vs_HPX_cluster, use_partition = FALSE, close_loop = TRUE)
cds_NMX_vs_HPX_graph_with_partition <- learn_graph(cds_NMX_vs_HPX_cluster, use_partition = TRUE, close_loop = TRUE)
plot_cells(cds_NMX_vs_HPX_graph)
plot_cells(cds_NMX_vs_HPX_graph_with_partition)

cds_NMX_vs_HPX_ordered <- order_cells(cds_NMX_vs_HPX_graph)
cds_NMX_vs_HPX_ordered_with_partition <- order_cells(cds_NMX_vs_HPX_graph_with_partition)

plot_cells(cds_NMX_vs_HPX_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds_NMX_vs_HPX_ordered_with_partition,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# Working with 3D trajectories
cds_NMX_vs_HPX_red_dim_3d <- reduce_dimension(cds_NMX_vs_HPX_preproc, max_components = 3,preprocess_method = 'PCA',reduction_method = 'UMAP')
cds_NMX_vs_HPX_cluster_3d <- cluster_cells(cds_NMX_vs_HPX_3d)
cds_NMX_vs_HPX_graph_3d <- learn_graph(cds_NMX_vs_HPX_cluster_3d)
plot_cells_3d(cds_NMX_vs_HPX_graph_3d, color_cells_by="partition")
cds_NMX_vs_HPX_graph_without_partition_3d <- learn_graph(cds_NMX_vs_HPX_graph_3d, use_partition = FALSE, close_loop = TRUE)
plot_cells_3d(cds_NMX_vs_HPX_graph_without_partition_3d, color_cells_by="partition")

cds_NMX_vs_HPX_ordered_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")