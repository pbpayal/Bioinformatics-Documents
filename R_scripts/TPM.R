# First, import the GTF-file that you have also used as input for htseq-count
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/Users/pbanerjee/Documents/Genomes/Homo_sapiens.GRCh37.75.gtf",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
write.table(exonic.gene.sizes, "/Users/pbanerjee/Desktop/genecounts_GCRh37", sep="\t")

x <- unlist(exonic.gene.sizes)
matrix <- matrix(x)
# RPKM versus TPM
# 
# RPKM and TPM are both normalized for library size and gene length.

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

genes <- data.frame(
  Gene = c("A","B","C","D","E"),
  Length = c(100, 50, 25, 5, 1)
)

counts <- data.frame(
  S1 = c(80, 10,  6,  3,   1),
  S2 = c(20, 20, 10, 50, 400)
)


tpms <- apply(counts, 2, function(x) tpm(x, genes$Length))
