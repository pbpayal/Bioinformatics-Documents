#########################################################
### link - http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
#########################################################

#########################################################
### A) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/Maria_dendogram_andheatmap/dendogram_data_maria.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names

png("/Users/pbanerjee/Desktop/maria_dendogram1.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 12*300,
    res = 400,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

d <- dist((mat_data), method = "euclidean") # distance matrix
hc <- hclust(d, method="ward.D2")
# display dendogram
plot(hc, labels = NULL, hang = 0.1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = NULL, ylab = "Average_Expression", cex = 0.6)

# Optional
# Draw dendogram with red borders around the 5 clusters
groups <- cutree(hc, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(hc, k=5, border="red")
plot(hc, labels = NULL, hang = 0.1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = NULL, ylab = "Average_Expression", cex = 0.6)

#Dendrograms with ggdendro
# remember to install the package: install.packages("ggdendro")
library(ggplot2)
#install.packages("ggdendro")
library(ggdendro)

# basic option
plot(ggdendrogram(hc, theme_dendro = FALSE),cex = 0.6 )

#another option
plot(ggdendrogram(hc, rotate = TRUE, size = 4, theme_dendro = FALSE), hang = -1, cex = 0.6)

# # tweeking some parameters
# op = par(bg = "#DDE3CA")
# plot(fit, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
#      col.axis = "#F38630", lwd = 3, lty = 3, sub = "", hang = -1, axes = FALSE)
# # add axis
# axis(side = 2, at = seq(0, 400, 100), col = "#F38630", labels = FALSE, 
#      lwd = 2)
# # add text in margin
# mtext(seq(0, 400, 100), side = 2, at = seq(0, 400, 100), line = 1, 
#       col = "#A38630", las = 2)
# par(op)
# 
# 
# # Alternative dendrograms
# # 
# # An alternative way to produce dendrograms is to specifically convert hclust objects into dendrograms objects.
# 
# # using dendrogram objects
# hcd = as.dendrogram(fit)
# # alternative way to get a dendrogram
# plot(hcd)
# # using dendrogram objects
# plot(hcd, type = "triangle")
# 
# # Zooming-in on dendrograms
# # 
# # Another very useful option is the ability to inspect selected parts of a given tree. For instance, if we wanted to examine the top partitions of the dendrogram, we could cut it at a height of 75
# 
# # plot dendrogram with some cuts
# op = par(mfrow = c(2, 1))
# plot(cut(hcd, h = 50)$upper, main = "Upper tree of cut at h=50")
# plot(cut(hcd, h = 50)$lower[[2]], main = "Second branch of lower tree with cut at h=50")
# par(op)

# 
# Phylogenetic trees
# 
# A very nice tool for displaying more appealing trees is provided by the R package ape. In this case, what we need is to convert the hclust objects into phylo pbjects with the funtions as.phylo

# load package ape; remember to install it: install.packages('ape')
library(ape)
# plot basic tree
plot(as.phylo(hc), cex = 0.9, label.offset = 1)

# Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)
# Default plot
plot(hcd, type = "rectangle", ylab = "Average_Expression")
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")

# Horizontal plot
plot(hcd,  xlab = "Average_Expression",
     nodePar = nodePar, horiz = TRUE)

# Unrooted
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)

# Fan
plot(as.phylo(hc), type = "fan")


# Radial
plot(as.phylo(hc), type = "radial")

# Cut the dendrogram into 3 clusters
colors = c("red", "blue", "green")
clus3 = cutree(hc, 3)
plot(as.phylo(hc), type = "fan", tip.color = colors[clus3],
     label.offset = 1, cex = 0.7)