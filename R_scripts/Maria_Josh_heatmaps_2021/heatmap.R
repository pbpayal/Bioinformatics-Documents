#########################################################
### Load previous image
#########################################################
load("/Users/pbanerjee/Documents/Payal_Scripts/R/Maria_Josh_heatmaps_2021/data.RData") 

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
library(ggplot2)

# Female
#########################################################
### B) Reading in data and transform it into matrix format
#########################################################
data <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/Maria_Josh_heatmaps_2021/Comp1_Females.csv", comment.char="#")
rownames(data) <- data$SYMBOL
data <- subset(data, select = -c(SYMBOL))
data_mat <- as.matrix(data)
# Females
# non-unique value when setting 'row.names': ‘PSD3’
# SYMBOL mtF1.wtF1 mtF2.wtF2
# PSD3   -0.082    -0.323
# PSD3  -0.050    -0.335
# removed second PSD3
#########################################################
### C) Customizing and plotting the heat map
#########################################################
# #display all colour schemes
# display.brewer.all()
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green","red"))
# # # (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(4,6,length=200),
#                seq(7,10,length=200))
# Red-Green Palette
# png("female_heatmap.png",
#     width = 3*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)
# female_heatmap <- heatmap.2(data_mat,
#           #cellnote = mat_data,  # same data set for cell labels
#           main = "Female", # heat map title
#           #notecol="blue",      # change font color of cell labels to black
#           # density.info="none",  # turns off density plot inside color legend
#           # trace="none",         # turns off trace lines inside the heat map
#           margins =c(2,5),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           # breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="row",     # only draw a row dendrogram
#           Colv="NA",         # turn off column clustering
#           key=T,
#           cexRow = 0.39,
#           cexCol = 0.75,
#           trace="none",
#           lhei=c(2,8), lwid=c(2,3.5), 
#           keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
# dev.off()
############
## Tiff
############
tiff("female_heatmap.tiff",
     width = 4*300,        # 5 x 300 pixels
     height = 7*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)
female_heatmap <- heatmap.2(data_mat,
                            #cellnote = mat_data,  # same data set for cell labels
                            main = "Female", # heat map title
                            #notecol="blue",      # change font color of cell labels to black
                            # density.info="none",  # turns off density plot inside color legend
                            # trace="none",         # turns off trace lines inside the heat map
                            margins =c(2,5),     # widens margins around plot
                            col=my_palette,       # use on color palette defined earlier
                            # breaks=col_breaks,    # enable color transition at specified limits
                            dendrogram="none",     # only draw a row dendrogram
                            Colv="NA",         # turn off column clustering
                            key=T,
                            cexRow = 0.39,
                            cexCol = 0.75,
                            trace="none",
                            lhei=c(2,8), lwid=c(2,3.5),
                            keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
dev.off()
############
### PDF 
############
# ExportPlot <- function(gplot, filename, width=2, height=1.5) {
#   # Export plot in PDF and EPS.
#   # Notice that A4: width=11.69, height=8.27
#   ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
#   postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
#   print(gplot)
#   dev.off()
# }
pdf("female_heatmap.pdf",
    width = 3.5,        
    height = 6,
    # res = 300,            # 300 pixels per inch
    pointsize = 8)
female_heatmap <- heatmap.2(data_mat,
                            #cellnote = mat_data,  # same data set for cell labels
                            main = "Female", # heat map title
                            #notecol="blue",      # change font color of cell labels to black
                            density.info="none",  # turns off density plot inside color legend
                            # trace="none",         # turns off trace lines inside the heat map
                            margins =c(2,5),     # widens margins around plot
                            col=my_palette,       # use on color palette defined earlier
                            # breaks=col_breaks,    # enable color transition at specified limits
                            # Rowv = FALSE, # turn off column clustering
                            dendrogram="row",     # only draw a row dendrogram
                            Colv="NA",
                            key=T,
                            cexRow = 0.39,
                            cexCol = 0.75,
                            trace="none",
                            lhei=c(2,8), lwid=c(2,3.5), 
                            keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
# ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height, dpi=dpi)
dev.off()

# #############
# ## EPS
# #############
# # gplot=female_heatmap
# # filename="female_heatmap"
# # width=2
# # height=6
# # dpi=300
# # postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
# setEPS()
# postscript("female_heatmap.eps", width = 4, height = 11)
# female_heatmap <- heatmap.2(data_mat,
#                             #cellnote = mat_data,  # same data set for cell labels
#                             main = "Female", # heat map title
#                             #notecol="blue",      # change font color of cell labels to black
#                             # density.info="none",  # turns off density plot inside color legend
#                             # trace="none",         # turns off trace lines inside the heat map
#                             margins =c(2,5),     # widens margins around plot
#                             col=my_palette,       # use on color palette defined earlier
#                             # breaks=col_breaks,    # enable color transition at specified limits
#                             dendrogram="row",     # only draw a row dendrogram
#                             Colv="NA",         # turn off column clustering
#                             key=T,
#                             cexRow = 0.39,
#                             cexCol = 0.75,
#                             trace="none",
#                             lhei=c(2,8), lwid=c(2,3.5), 
#                             keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
# dev.off()

# Male
data2 <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/Maria_Josh_heatmaps_2021/Comp2_Male.csv", comment.char="#")
rownames(data2) <- data2$SYMBOL
data2 <- subset(data2, select = -c(SYMBOL))
data_mat2 <- as.matrix(data2)
# Females
# non-unique value when setting 'row.names': ‘PSD3’
# SYMBOL mtF1.wtF1 mtF2.wtF2
# PSD3  0.15	0.38
# PSD3  0.13	0.18

# removed second PSD3
#########################################################
### C) Customizing and plotting the heat map
#########################################################
# #display all colour schemes
# display.brewer.all()
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green","red"))
# # # (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(4,6,length=200),
#                seq(7,10,length=200))
# Red-Green Palette
# png("male_heatmap.png",
#     width = 3*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)
# male_heatmap <- heatmap.2(data_mat2,
#                             #cellnote = mat_data,  # same data set for cell labels
#                             main = "Male", # heat map title
#                             #notecol="blue",      # change font color of cell labels to black
#                             # density.info="none",  # turns off density plot inside color legend
#                             # trace="none",         # turns off trace lines inside the heat map
#                             margins =c(2,5),     # widens margins around plot
#                             col=my_palette,       # use on color palette defined earlier
#                             # breaks=col_breaks,    # enable color transition at specified limits
#                             dendrogram="row",     # only draw a row dendrogram
#                             Colv="NA",         # turn off column clustering
#                             key=T,
#                             cexRow = 0.39,
#                             cexCol = 0.75,
#                             trace="none",
#                             lhei=c(2,7), lwid=c(2,3.5), 
#                             keysize=0.5, key.par = list(cex=0.7), srtCol = 360)    
# dev.off()
############
## Tiff
############
tiff("male_heatmap.tiff",
     width = 4*300,        # 5 x 300 pixels
     height = 7*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)
male_heatmap <- heatmap.2(data_mat2,
                          #cellnote = mat_data,  # same data set for cell labels
                          main = "Male", # heat map title
                          #notecol="blue",      # change font color of cell labels to black
                          # density.info="none",  # turns off density plot inside color legend
                          # trace="none",         # turns off trace lines inside the heat map
                          margins =c(2,5),     # widens margins around plot
                          col=my_palette,       # use on color palette defined earlier
                          # breaks=col_breaks,    # enable color transition at specified limits
                          dendrogram="row",     # only draw a row dendrogram
                          Colv="NA",         # turn off column clustering
                          key=T,
                          cexRow = 0.39,
                          cexCol = 0.75,
                          trace="none",
                          lhei=c(2,7), lwid=c(2,3.5),
                          keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
dev.off()
############
### PDF 
############
pdf("male_heatmap.pdf",
    width = 3.5,        
    height = 6,
    # res = 300,            # 300 pixels per inch
    pointsize = 8)
male_heatmap <- heatmap.2(data_mat2,
                          #cellnote = mat_data,  # same data set for cell labels
                          main = "Male", # heat map title
                          #notecol="blue",      # change font color of cell labels to black
                          # density.info="none",  # turns off density plot inside color legend
                          # trace="none",         # turns off trace lines inside the heat map
                          margins =c(2,5),     # widens margins around plot
                          col=my_palette,       # use on color palette defined earlier
                          # breaks=col_breaks,    # enable color transition at specified limits
                          dendrogram="row",     # only draw a row dendrogram
                          Colv="NA",         # turn off column clustering
                          key=T,
                          cexRow = 0.39,
                          cexCol = 0.75,
                          trace="none",
                          lhei=c(2,7), lwid=c(2,3.5),
                          keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
dev.off()
# Female vs Male
data3 <- read.csv("/Users/pbanerjee/Documents/Payal_Scripts/R/Maria_Josh_heatmaps_2021/Comp3_Females_vs_Males.csv", comment.char="#")
rownames(data3) <- data3$SYMBOL
# non-unique value when setting 'row.names': ‘MBP’ 
# SYMBOL2	mtF1/mtM1	mtF2/mtM2
# MBP	-0.01	-0.3
# MBP	-0.14	-0.23
# removed MBP	-0.01	-0.3
data3 <- subset(data3, select = -c(SYMBOL))
data_mat3 <- as.matrix(data3)
# Females
# non-unique value when setting 'row.names': ‘PSD3’
# SYMBOL mtF1.wtF1 mtF2.wtF2
# PSD3  0.15	0.38
# PSD3  0.13	0.18

# removed second PSD3
#########################################################
### C) Customizing and plotting the heat map
#########################################################
# #display all colour schemes
# display.brewer.all()
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green","red"))
# # # (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(4,6,length=200),
#                seq(7,10,length=200))
# Red-Green Palette
# png("female_vs_male_heatmap.png",
#     width = 3*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)
# female_vs_male_heatmap <- heatmap.2(data_mat3,
#                             #cellnote = mat_data,  # same data set for cell labels
#                             main = "Female vs Male", # heat map title
#                             #notecol="blue",      # change font color of cell labels to black
#                             # density.info="none",  # turns off density plot inside color legend
#                             # trace="none",         # turns off trace lines inside the heat map
#                             margins =c(2,5),     # widens margins around plot
#                             col=my_palette,       # use on color palette defined earlier
#                             # breaks=col_breaks,    # enable color transition at specified limits
#                             dendrogram="row",     # only draw a row dendrogram
#                             Colv="NA",         # turn off column clustering
#                             key=T,
#                             cexRow = 0.39,
#                             cexCol = 0.75,
#                             trace="none",
#                             lhei=c(2,7), lwid=c(2,3.5), 
#                             keysize=0.5, key.par = list(cex=0.7), srtCol = 360)    
# dev.off()
############
## Tiff
############
tiff("female_vs_male_heatmap.tiff",
     width = 4*300,        # 5 x 300 pixels
     height = 5*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)
female_vs_male_heatmap <- heatmap.2(data_mat3,
                                    #cellnote = mat_data,  # same data set for cell labels
                                    main = "Female vs Male", # heat map title
                                    #notecol="blue",      # change font color of cell labels to black
                                    # density.info="none",  # turns off density plot inside color legend
                                    # trace="none",         # turns off trace lines inside the heat map
                                    margins =c(2,5),     # widens margins around plot
                                    col=my_palette,       # use on color palette defined earlier
                                    # breaks=col_breaks,    # enable color transition at specified limits
                                    dendrogram="row",     # only draw a row dendrogram
                                    Colv="NA",         # turn off column clustering
                                    key=T,
                                    cexRow = 0.39,
                                    cexCol = 0.75,
                                    trace="none",
                                    lhei=c(2,7), lwid=c(2,3.5),
                                    keysize=0.5, key.par = list(cex=0.7), srtCol = 360)
dev.off()
############
### PDF 
############
pdf("female_vs_male_heatmap.pdf",
    width = 3.5,        
    height = 6,
    # res = 300,            # 300 pixels per inch
    pointsize = 8)
female_vs_male_heatmap <- heatmap.2(data_mat3,
                                    #cellnote = mat_data,  # same data set for cell labels
                                    main = "Female vs Male", # heat map title
                                    #notecol="blue",      # change font color of cell labels to black
                                    # density.info="none",  # turns off density plot inside color legend
                                    # trace="none",         # turns off trace lines inside the heat map
                                    margins =c(2,5),     # widens margins around plot
                                    col=my_palette,       # use on color palette defined earlier
                                    # breaks=col_breaks,    # enable color transition at specified limits
                                    dendrogram="row",     # only draw a row dendrogram
                                    Colv="NA",         # turn off column clustering
                                    key=T,
                                    cexRow = 0.39,
                                    cexCol = 0.75,
                                    trace="none",
                                    lhei=c(2,7), lwid=c(2,3.5), 
                                    keysize=0.5, key.par = list(cex=0.7), srtCol = 360)    
dev.off()
#########################################################
### Save Image workspace
save.image()
#########################################################

