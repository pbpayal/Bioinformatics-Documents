# First install Bioconductor and Monocle 3
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install()

# Next install a few more dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))

# Install Monocle3 if not already installed
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')devtools::install_github('cole-trapnell-lab/monocle3')

# Install a few Garnett dependencies:
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db','org.Ss.eg.db'))

# Install Garnett
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
library(garnett)
library(org.Ss.eg.db)
set.seed(260)


garnett_markers <- top_markers_cluster %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)

generate_garnett_marker_file(garnett_markers, file="/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/marker_file.txt")


ss_brain_classifier <- train_cell_classifier(cds = cds_monocle3,
                                         marker_file = "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Camille_Single_cell/marker_file.txt",
                                         db=org.Ss.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         classifier_gene_id_type = "ENSEMBL",
                                         marker_file_gene_id_type = "ENSEMBL")

cds_monocle3 <- load_cellranger_data("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/Camille_Single_cell/1_5_merge/1_5_merge_monocle3")

classified_cds_monocle3 <- classify_cells(cds_monocle3, pbmc_classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")
