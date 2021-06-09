# TSCAN
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/TSCAN")
library(TSCAN)
TSCANui()

data(lpsdata)

# Preprocess gene expression data
procdata <- preprocess(lpsdata)
# Constructing Pseudotime
lpsmclust <- exprmclust(procdata)
plotmclust(lpsmclust)
# Use the TSCANorder function to obtain the TSCAN ordering.
lpsorder <- TSCANorder(lpsmclust)
lpsorder
# Testing Differentially Expressed Genes
diffval <- difftest(procdata,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
head(row.names(diffval)[diffval$qval < 0.05])
# Use singlegeneplot function to plot the expression value of 
# a single gene againsta given pseudotime.  
# Notice that here orderonly should be set to FALSE.
STAT2expr <- log2(lpsdata["STAT2",]+1)
singlegeneplot(STAT2expr, TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))