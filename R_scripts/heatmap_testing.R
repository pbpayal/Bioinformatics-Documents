library("DESeq2")
library(devtools)
library("RColorBrewer")
# If installed, just load the library, if not installed first install gplots
#install.packages("gplots")
library("gplots")
library("ggplot2")
directory <- "/Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Panos/Deseq2_files"
setwd(directory)

sampleFiles4<- c("deseq2_fcounts_star_out_PN1_lane1_S19_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PN2_lane1_S20_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH1_lane1_S16_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt",
                 "deseq2_fcounts_star_out_PH2_lane1_S17_L001_R1_001.fastq.gzAligned.sortedByCoord.out.bam.txt")

sampleNames4 <- c("PN1","PN2","PH1","PH2")
sampleCondition4 <- c("control","control","exp","exp")
sampleTable4 <- data.frame(sampleName = sampleNames4, fileName = sampleFiles4, condition4 = sampleCondition4)
treatments4 = c("control","exp")

# Load the data in Deseq format
ddsHTSeq4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable4, directory = directory, design = ~ condition4)

colData(ddsHTSeq4)$condition4 <- factor(colData(ddsHTSeq4)$condition4,levels = treatments4)

#Perform the Differential Gene Expression using DESeq
dds4 <- DESeq(ddsHTSeq4)
res4 <- results(dds4)
res4
summary(res4)
vsd4 <- varianceStabilizingTransformation(dds4, blind=T)
# select 35 genes with the highest variance across samples
topVarGenes <- head( order( rowVars( assay(vsd4) ), decreasing=TRUE ), 35 )
x <- assay(vsd4)[topVarGenes,]
# The heatmap becomes more interesting if we do not 
# look at absolute expression strength but rather at the 
# amount by which each gene deviates in a specific sample 
# from the gene’s average across all samples. 
# Hence, we center and scale each genes’ values across samples, 
# and plot a heatmap.
heatmap.2(assay(vsd4)[topVarGenes,],
          #cellnote = mat_data,  # same data set for cell labels
          main = "35 Top Expressed Genes", # heat map title
          #notecol="blue",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,16),     # widens margins around plot
#          col=my_palette,       # use on color palette defined earlier
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA"            # turn off column clustering
          )

library("pheatmap")
mat4 = assay(vsd4)[ head(order(res4$padj),30), ] # select the top 30 genes with the lowest padj
mat4 = mat4 - rowMeans(mat4) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df4 = as.data.frame(colData(vsd4)[,c("condition4")]) # Create a dataframe with a column of the conditions
colnames(df4) = "condition" # Rename the column header
rownames(df4) = colnames(mat4) # add rownames
# and plot the actual heatmap
pheatmap(mat4, annotation_col=df4, main = "Top 30 Gene Expression")