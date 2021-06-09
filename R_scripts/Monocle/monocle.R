source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("monocle")
library(monocle)
# version used for this post
packageVersion('monocle')
data_dir <- '/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/raw_feature_bc_matrix'

library(Matrix)
matrix_dir = " /Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix"
barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix/features.tsv.gz", 
                           header = FALSE,
                           stringsAsFactors = FALSE)
names(feature.names)[2] <- "gene_short_name"
names(feature.names)[1] <- "ensemble_id"
barcode.names = read.delim("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3/outs/raw_feature_bc_matrix/barcodes.tsv.gz", 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


final_mat <- cbind(mat,feature.names$V2)



#rows
head(mat[,1])
#columns9
head(mat[1,])

mat[3000:3003,3000:3003]

# final_datamat <- read.csv2("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/raw_feature_bc_matrix/final_matrix.csv", sep = ",", stringsAsFactors = FALSE)
# HSMM <- newCellDataSet(as.matrix(final_datamat, "sparseMatrix"),
#                        lowerDetectionLimit = 0.5,
#                        expressionFamily = negbinomial.size())

# HSMM_1 <- newCellDataSet(as.matrix(mat, "sparseMatrix"),
#                        lowerDetectionLimit = 0.5,
#                        expressionFamily = negbinomial.size())
# HSMM <- detectGenes(HSMM, min_expr = 0.1)
# print(head(fData(HSMM)))
# HSMM <- detectGenes(HSMM, min_expr = 0.1)
# print(head(fData(HSMM)))
# expressed_genes <- row.names(subset(fData(HSMM),
#                                     X1 >= 10))
# print(head(pData(HSMM)))










## external 
cell_meta <- data.frame(cell = colnames(final_datamat),
                        cell_type = annot$labels[match(colnames(mat_filtered), rownames(annot))],
                        stringsAsFactors = FALSE)
rownames(cell_meta) <- colnames(mat_filtered)
gene_meta <- gns %>% 
  filter(gene %in% rownames(mat_filtered), !is.na(gene)) %>% 
  rename(id = gene, gene_short_name = gene_name)
rownames(gene_meta) <- gene_meta$id
gene_meta <- gene_meta[rownames(mat_filtered),]
cds <- newCellDataSet(mat_filtered, cell_metadata = cell_meta, gene_metadata = gene_meta)

