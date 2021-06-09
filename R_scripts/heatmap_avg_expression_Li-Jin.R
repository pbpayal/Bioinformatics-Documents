#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
library(cowplot)
library(pheatmap)

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/Joseph_final_results/Li-Jin_diagrams/heatmap_avg_expression.csv", comment.char="#")
rownames(data) <- data$Gene_names
# 'Coasy’, ‘Coq6’, ‘Mtor’ present in more than 1 cluster
data <- read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/Joseph_final_results/Li-Jin_diagrams/avg_expression_log_transformed_final.csv", comment.char="#")
rownames(data) <- data$Gene_names

# ‘Coq6’ in ECM and Migrant Neurons downregulated in ECM and upregulated in Migrant Neurons
# removed 'Coq6' from migrant neurons and kept the one in ECM n because the p-adj value is better there

# 'Coasy’ in ECM and GABA upregulated in both
# removed 'Coasy’ from ECM and kept the one in GABA because the p-adj value is better there

# 'Mtor' in ECM and Metabolic upregulated in both
# removed 'Mtor' from metabolic and kept the one in ECM because the p-adj value is better there
data_subset <- subset(data, select = -c(Gene_names, Cluster))

# # (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=70),
               seq(0.1,1,length=70),
               seq(1.1,2,length=70))

heatmap1 <- pheatmap(data_subset, main = "Average Gene Expression P8 GAD65GFP", 
         show_rownames=T, cluster_rows=F, cluster_cols=F, show_colnames = T,
         angle_col = 0,
         fontsize_row = 6,
         breaks = seq(-1.5,2.5, length = 95))
         #breaks= col_breaks)
heatmap1
save_plot("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Joseph/Joseph_final_results/Li-Jin_diagrams/all_gene_avg_expression_1.png",heatmap1)
#########################################################
### C) Customizing and plotting the heat map
#########################################################

#display all colour schemes
display.brewer.all()

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green","red"))
   


# # Red-Green Palette
# heatmap.2(mat_data,
#           #cellnote = mat_data,  # same data set for cell labels
#           main = "Correlation", # heat map title
#           #notecol="blue",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(12,9),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="row",     # only draw a row dendrogram
#           Colv="NA",            # turn off column clustering
#           key=T,
#           cexRow = 0.75,
#           cexCol = 0.75,
#           trace="none",
#           lhei=c(2,4), lwid=c(2,3.5), 
#           keysize=0.5, key.par = list(cex=0.5))            

data_subset_matrix <- as.matrix(data_subset)
#My preferred color Palette
heatmap.2(data_subset_matrix,
          #cellnote = mat_data,  # same data set for cell labels
          main = "Average Gene Expression", # heat map title
          notecol="blue",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=brewer.pal(11,"RdBu"), # use on color palette defined earlier
          scale="none",
#          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          Colv="NA",           # turn off column clustering
          cexRow = 0.7,
          cexCol = 1)
