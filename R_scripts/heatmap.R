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


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/josh_new_diagram/Maria_heatmap2.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names


#########################################################
### C) Customizing and plotting the heat map
#########################################################

#display all colour schemes
display.brewer.all()

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green","red"))

# # (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(4,6,length=200),
               seq(7,10,length=200))

# creates a 5 x 5 inch image
png("/Users/pbanerjee/Documents/CBU/CNMC_Projects/Neuroscience/josh_new_diagram/maria_red_green_heatmap_testing.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 12*300,
    res = 400,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

# Red-Green Palette
heatmap.2(mat_data,
          #cellnote = mat_data,  # same data set for cell labels
          main = "Correlation", # heat map title
          #notecol="blue",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",            # turn off column clustering
          key=T,
          cexRow = 0.75,
          cexCol = 0.75,
          trace="none",
          lhei=c(2,4), lwid=c(2,3.5), 
          keysize=0.5, key.par = list(cex=0.5))            

# #My preferred color Palette
# heatmap.2(mat_data,
#           #cellnote = mat_data,  # same data set for cell labels
#           main = "Periodic Gene Expression", # heat map title
#           notecol="blue",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           col=brewer.pal(11,"RdBu"), # use on color palette defined earlier
#           scale="column",
# #          breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="row",     # only draw a row dendrogram
#           Colv="NA",           # turn off column clustering
#           cexRow = 0.7, 
#           cexCol = 1)

dev.off()               # close the PNG device