# Read 10X data
sample1 <- Read10X("../..")
# Seurat object
sample1_seurat <- CreateSeuratObject(counts = sample1, 
                                     min.cells = 3, 
                                     min.features = 200, 
                                     project = "sample1")
# Read 10X data
sample2 <- Read10X("../..")
# Seurat object
sample2_seurat <- CreateSeuratObject(counts = sample2, 
                                     min.cells = 3, 
                                     min.features = 200, 
                                     project = "sample2")
# Merge
merge_1_2 <- merge(sample1_seurat, y = sample2_seurat, add.cell.ids = c("sample1", "sample2"), project = "sample_1_2")
# Plot and save violin plot
unfiltered_violin_plot <- VlnPlot(merge_1_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
unfiltered_feature_scatter_plot <- FeatureScatter(merge_1_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Gene Level filetring - Filtering specific Mitochondrial genes for Camille project 
# Get the counts from assay
counts <- GetAssayData(merge_1_2, assay = "RNA")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_cell_gene <- CreateSeuratObject(filtered_counts, meta.data = merge_1_2@meta.data)
# Again extract counts
counts1 <- GetAssayData(filtered_cell_gene, assay = "RNA")
# Remove the list of genes from the counts
no_mito <- counts1[-(which(rownames(counts1) %in% c('ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5',
                                                      'ND6', 'CYTB','ENSSSCG00000018060','ENSSSCG00000018061',
                                                      'ENSSSCG00000018062','ENSSSCG00000018063',
                                                      'ENSSSCG00000018064','ENSSSCG00000018065',
                                                      'ENSSSCG00000018066','ENSSSCG00000018067',
                                                      'ENSSSCG00000018068','ENSSSCG00000018069',
                                                      'ENSSSCG00000018070','ENSSSCG00000018071',
                                                      'ENSSSCG00000018072','ENSSSCG00000018073',
                                                      'ENSSSCG00000018074','ENSSSCG00000018075',
                                                      'ENSSSCG00000018076','ENSSSCG00000018077',
                                                      'ENSSSCG00000018078','ENSSSCG00000018079',
                                                      'ENSSSCG00000018080','ENSSSCG00000018081',
                                                      'ENSSSCG00000018082','ENSSSCG00000018083',
                                                      'ENSSSCG00000018084','ENSSSCG00000018085',
                                                      'ENSSSCG00000018086','ENSSSCG00000018087',
                                                      'ENSSSCG00000018088','ENSSSCG00000018089',
                                                      'ENSSSCG00000018090','ENSSSCG00000018091',
                                                      'ENSSSCG00000018092','ENSSSCG00000018093',
                                                      'ENSSSCG00000018094','ENSSSCG00000018095',
                                                      'ENSSSCG00000018096'))),]
# Putting it back into the seurat object
filtered_mito_object <- subset(filtered_cell_gene, features = rownames(no_mito))
# number of genes per UMI for each cell 
filtered_mito_object$log10GenesPerUMI <- log10(filtered_mito_object$nFeature_RNA) / log10(filtered_mito_object$nCount_RNA)
# Cell Level Filtering
filtered_mito.cell_object <- subset(x = filtered_mito_object, 
                                      subset= (nCount_RNA >= 500) & 
                                        (nFeature_RNA >= 250) & 
                                        (log10GenesPerUMI > 0.80))
filtered_feature_scatter_plot <- FeatureScatter(filtered_mito.cell_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
filtered_violin_plot <- VlnPlot(filtered_mito.cell_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


  