if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

library(fgsea)
library(dplyr)

# Read your .rnk file
rnk_file <- read.table("WM_24hr_JMI_fgsea.rnk", sep = "\t", header=TRUE, colClasses = c("character", "numeric"))
str(rnk_file)
# Set vlaues of the names acoording to this format
ranks_final <- setNames(rnk_file$padj,rnk_file$gene)
str(ranks_final)
View(ranks_final)

# Load the pathways
pathways <- gmtPathways("msigdb.v7.0.symbols.gmt")
str(head(pathways))

# Run fgsea
fgseaRes <- fgsea(pathways, ranks_final, minSize = 15, maxSize = 5000, nperm = 1000)

