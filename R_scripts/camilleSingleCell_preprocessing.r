
###Reading hypox data
hypox.data <- Read10X(data.dir = "V:/Payal_Projects/Camille/Camille_data_Payal/hypox_2_6_8_merge/outs/filtered_feature_bc_matrix")
###Creating Seurat object
hypox.seurat <- CreateSeuratObject(counts = hypox.data, project = "hypox")
###Reading a text file containing the gene list
r1 <- read.table("X:/Suro/SingleCellRNAseq/Camille/gene_result.txt", header = TRUE)
###getting the assay information
counts1 <- GetAssayData(hypox.seurat, assay = "RNA")
###matching with the gene list.
m1 <- match(gene, rownames(counts1))
###Deleting matches from the counts
counts2 <- counts1[-m1,]
###Putting it back into the seurat object
hypox.seurat <- subset(hypox.seurat, features = rownames(counts2))
###writing the seurat object to be functional for cloupe usage
write10xCounts(
     path = "X:/Suro/SingleCellRNAseq/Camille/Results/DE/",
     x = counts2,
     barcodes = colnames(counts2),
     gene.id = rownames(counts),
     gene.type = "Gene Expression",
     overwrite = FALSE,
     type = c("sparse"),
     genome = "SusScrofa",
     version = c("3"),
     chemistry = "Single Cell 3' v3",
     original.gem.groups = 1L,
     library.ids = "custom"
)
