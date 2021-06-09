#Session 4: Introduction of scRNA-seq Data Analysis (Hands-on Lab)

#load required packages
library(SingleR)
library(Seurat)
library(cluster)
library(umap)
library(pheatmap)
library(viper)
library(MAST)
library(ggplot2)
library(ggrepel)
library(infercnv)
library(MASS)
library(SeuratData)
library(Hmisc)

#Seurat walkthrough-- data loading and cleaning
AvailableData()

InstallData("ifnb")
LoadData("ifnb")
table(ifnb$stim)

#Downsample dataset for ease of analysis
set.seed(1234)
ifnb=ifnb[,sample(1:ncol(ifnb),2000)]
table(ifnb$stim)

#visualize data quality distributions
VlnPlot(ifnb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size=0,group.by="stim")
ncol(ifnb)
##filter out low-quality cells and potential doublets
ifnb <- subset(ifnb, subset = nCount_RNA > 1000 & nCount_RNA < 25000)
VlnPlot(ifnb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size=0,group.by="stim")
ncol(ifnb)

#split data by experiment batch
ifnb.list <- SplitObject(ifnb, split.by = "stim")

#Normalize Data
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)

#Integrate to Perform Batch Correction
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", 
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
##NOTE: To speed up this process for large datasets and maintain consistency as you collect new data, can set one high-quality sample as a reference

#PCA Dimensionality Reduction
immune.combined.sct <- RunPCA(immune.combined.sct, features=VariableFeatures(object=immune.combined.sct))
ElbowPlot(immune.combined.sct)
immune.combined.sct <- RunLDA(immune.combined.sct, features=VariableFeatures(object=immune.combined.sct),labels = immune.combined.sct$seurat_annotations)

#Visualization
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- RunTSNE(immune.combined.sct, reduction = "pca", dims = 1:30)

DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="seurat_annotations") + NoLegend()
DimPlot(immune.combined.sct, reduction = "tsne",label = TRUE,repel=T,group.by="seurat_annotations") + NoLegend()
DimPlot(immune.combined.sct, reduction = "pca",label = TRUE,repel=T,group.by="seurat_annotations") + NoLegend()
DimPlot(immune.combined.sct, reduction = "lda",label = TRUE,repel=T,group.by="seurat_annotations") + NoLegend()

DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="seurat_annotations",split.by="stim") + NoLegend()
DimPlot(immune.combined.sct, reduction = "tsne",label = TRUE,repel=T,group.by="seurat_annotations",split.by="stim") + NoLegend()
DimPlot(immune.combined.sct, reduction = "pca",label = TRUE,repel=T,group.by="seurat_annotations",split.by="stim") + NoLegend()
DimPlot(immune.combined.sct, reduction = "lda",label = TRUE,repel=T,group.by="seurat_annotations",split.by="stim") + NoLegend()

##density plots
ggplot(data.frame(UMAP_1=immune.combined.sct@reductions$umap@cell.embeddings[,1],UMAP_2=immune.combined.sct@reductions$umap@cell.embeddings[,2],cluster=immune.combined.sct$seurat_annotations),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+xlim(-15,15)+ylim(-15,20)

##Louvain Clustering
immune.combined.sct <- FindNeighbors(immune.combined.sct, dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution=seq(0.1,1,by=0.1), verbose = FALSE,algorithm=1) 
###compare high-resolution vs low-resolution Louvain Clustering
DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="integrated_snn_res.0.1") + NoLegend()
DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="integrated_snn_res.1") + NoLegend()

#PAM Clustering
nclust=5
clust5=pam(immune.combined.sct@reductions$pca@cell.embeddings,nclust)
immune.combined.sct$pam5=clust5$clustering
DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="pam5") + NoLegend()
nclust=20
clust20=pam(immune.combined.sct@reductions$pca@cell.embeddings,nclust)
immune.combined.sct$pam20=clust20$clustering
DimPlot(immune.combined.sct, reduction = "umap",label = TRUE,repel=T,group.by="pam20") + NoLegend()

#Silhouette Score Evaluation
plot(silhouette(clust5,dist(immune.combined.sct@reductions$pca@cell.embeddings)),col=1:5,border=NA)
plot(silhouette(clust20,dist(immune.combined.sct@reductions$pca@cell.embeddings)),col=1:20,border=NA)
summary(silhouette(clust5,dist(immune.combined.sct@reductions$pca@cell.embeddings)))$avg.width
summary(silhouette(clust20,dist(immune.combined.sct@reductions$pca@cell.embeddings)))$avg.width

silhouette_scores=sapply(5:10,function(x){
  clust=pam(immune.combined.sct@reductions$pca@cell.embeddings,x)
  summary(silhouette(clust,dist(immune.combined.sct@reductions$pca@cell.embeddings)))$avg.width
})
plot(5:10,silhouette_scores,pch=18)
lines(5:10,silhouette_scores)

##Hybrid Louvain Approach With Silhouette Scoring

#' Function to select optimal Louvain clustering of single-cell matrix from 10 
#' alternative resolution values. Sub-samples 1000 cells 100 times at each resolution
#' value to compute mean and standard deviation of silhouette score.
#' @param mat: matrix with rows as principal component vectors and columns as samples
#' @param clust: matrix with rows as samples as each column as a clustering vector for a given resolution
#' outputs list of mean silhouette scores and standard deviations of silhouette scores for each clustering. 
sil_subsample=function(mat,clust){
  out=as.data.frame(matrix(rep(NA,100*ncol(clust)),nrow=100))
  for(x in 1:100){
    i=sample(1:ncol(mat),min(1000,ncol(mat)))
    d=as.dist(1 - cor(mat[,i], method = "pearson"))
    for(j in 1:ncol(clust)){
      if(length(table(clust[i,j]))==1){out[x,j]=0}
      if(length(table(clust[i,j]))>1){
        sil=silhouette(as.numeric(clust[i,j]),d)
        out[x,j]=mean(sil[, "sil_width"])}}
  }
  means=apply(out,2,mean)
  sd=apply(out,2,sd)
  return(list(means,sd))
}
clust=immune.combined.sct@meta.data[,which(grepl("integrated_snn_res.",colnames(immune.combined.sct@meta.data)))]
mat=as.data.frame(t(immune.combined.sct$pca@cell.embeddings))
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.1,1,by=0.1)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
immune.combined.sct$seurat_clusters=immune.combined.sct@meta.data[,which(colnames(immune.combined.sct@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(immune.combined.sct) <- "seurat_clusters"
plot(DimPlot(immune.combined.sct, reduction = "umap",label = TRUE) + NoLegend())

##Differential Expression
#MAST, wilcox, roc

#'Function to plot heatmap of custom gene list grouped by cluster
#'truncates color scale to >5th percentile and <95th percentile, shades by quantile. 
#' @param dat: a matrix input with genes as rows and samples as columns
#' @param clust: a vector of cluster labels
#' @param genes: a vector of genes to plot
#' @param genes_by_cluster: a boolean indicator of whether genes are grouped by cluster identity (e.g. top5 genes in cluster 1 followed by top5 genes in cluster2, etc.)
#' @param n_top_genes_per_cluster: a number of top genes per cluster being plotted if genes_by_cluster=T. Otherwise ignored. 
#' @param color_palette: A manual color palette for clusters. Defaults to hue_pal(). If provided, must be vector of colors same length as the number of unique clusters. 
#' @param scaled: a boolean indicator of whether data are already scaled. Defaults to F, in which case row-wise z-score scaling is applied.
#' The function outputs a custom heatmap of manually specified genes, grouped by cluster, default
#' behavior is to apply row-wise z-score scaling and plot the top 5 genes per cluster. 
geneHeatmap_plot=function(dat,clust,genes,genes_by_cluster=T,n_top_genes_per_cluster=5,color_palette=NA,scaled=F){
  identities <- levels(clust)
  if(is.na(color_palette)){my_color_palette <- hue_pal()(length(identities))}
  else{my_color_palette=color_palette}
  features=genes
  i=sample(1:ncol(dat),min(10000,ncol(dat)),replace = F)
  x=dat[features,i]
  df <- data.frame(clust[i],clust[i])
  rownames(df)=colnames(x)
  colnames(df)=c("cluster","cluster2")
  anno_colors <- list(cluster = my_color_palette)
  names(anno_colors$cluster) <- levels(df$cluster)
  o=order(df$cluster)
  x=x[,o]
  df=df[o,]
  df=df[1]
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  if(scaled==F){t=as.matrix(apply(x,1,function(x){(x-mean(x))/sd(x)}))}
  if(scaled==T){t=as.matrix(x)}
  mat_breaks <- c(quantile_breaks(t[which(t<0)], n = 10),0,quantile_breaks(t[which(t>0)], n = 10))
  mat_breaks=mat_breaks[2:(length(mat_breaks)-1)] #restrict range of data to quantiles 5%-95%, extreme values excluded
  if(genes_by_cluster){
    anno_colors$group=anno_colors$cluster
    anno_row=data.frame(group=unlist(lapply(unique(df$cluster),function(x){rep(x,n_top_genes_per_cluster)})))
    gene_names=rownames(x)
    rownames(x)=1:nrow(x)
    if(!scaled){return(pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 10,show_colnames = F,annotation_colors = anno_colors,scale="row",gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F))}
    if(scaled){return(pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 10,show_colnames = F,annotation_colors = anno_colors,gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F))}
  }
  else{
    if(!scaled){return(pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 8,show_colnames = F,annotation_colors = anno_colors,scale="row"))}
    if(scaled){return(pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 8,show_colnames = F,annotation_colors = anno_colors))}
  }
}
library(dplyr)
library(scales)

markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "wilcox")
head(markers)
markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "roc")
head(markers)
markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
head(markers)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneHeatmap_plot(immune.combined.sct@assays$SCT@counts,as.factor(immune.combined.sct$seurat_clusters),top10$gene,n_top_genes_per_cluster = 5)


#Run SingleR cell type inference
immune.combined.sct.singler=CreateSinglerObject(immune.combined.sct[["SCT"]]@counts, annot = NULL, 
                                                project.name = "ifnb", min.genes = 0, 
                                                technology = "10X", species = "Human", citation = "", 
                                                do.signatures = F, clusters = NULL, numCores = numCores,
                                                fine.tune=F,variable.genes = "de",reduce.file.size = T,do.main.types = T)
immune.combined.sct$hpca_labels=immune.combined.sct.singler$singler[[1]][[1]][[2]]
immune.combined.sct$hpca_main_labels=immune.combined.sct.singler$singler[[1]][[4]][[2]]
immune.combined.sct$blueprint_labels=immune.combined.sct.singler$singler[[2]][[1]][[2]]
immune.combined.sct$blueprint_main_labels=immune.combined.sct.singler$singler[[2]][[4]][[2]]
immune.combined.sct$hpca_pvals=immune.combined.sct.singler$singler[[1]][[1]][[3]]
immune.combined.sct$hpca_main_pvals=immune.combined.sct.singler$singler[[1]][[4]][[3]]
immune.combined.sct$blueprint_pvals=immune.combined.sct.singler$singler[[2]][[1]][[3]]
immune.combined.sct$blueprint_main_pvals=immune.combined.sct.singler$singler[[2]][[4]][[3]]

#Filter SingleR labels to labels with p<0.1 and number of cells per label >50

l=immune.combined.sct$blueprint_main_labels
l[which(immune.combined.sct$blueprint_main_pvals>0.05)]=NA
table(l)
l[which(l %in% names(which(table(l)<50)))]=NA
immune.combined.sct$celltype=l
plot(DimPlot(immune.combined.sct, reduction = "umap",label=TRUE,group.by="celltype",repel=T))
#saveRDS(immune.combined.sct,"immune.combined.sct.rds")

table(immune.combined.sct$celltype,immune.combined.sct$seurat_clusters)

