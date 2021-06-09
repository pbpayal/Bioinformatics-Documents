library(Seurat)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
clutsers_filtered <- readRDS(file ="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/PART2/Pipeline_Final_Results/Sample_1_5_6_8/DEG/filtered_clusters_RNA_1_5_6_8.rds")
# Check if Default Assay id RNA
DefaultAssay(clutsers_filtered)
# If DefaultAssay is not set to "RNA" set it to RNA 
# DefaultAssay(clutsers_filtered) <- "RNA"
levels(clutsers_filtered)
# Cluster 0
deg_cluster0 <- FindConservedMarkers(clutsers_filtered, ident.1 = 0, grouping.var = "sample_id" , verbose = FALSE)
head(deg_cluster0, n = 15)

# Function for Conserved Markers and write out corresponding files
for(i in 0:28) {
  print(i)
  #deg_cluster <- FindMarkers(filtered_clusters, ident.1 = normox_1_5, ident.2 = hypox_6_8, verbose = FALSE)
  markers <- FindConservedMarkers(clutsers_filtered, ident.1 = i, grouping.var = "sample_id" , verbose = FALSE)
  write.table(markers, paste("conservedmarkers",i,".txt",sep = ""), sep = "\t", row.names = TRUE, col.names = TRUE)
}
