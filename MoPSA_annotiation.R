library(Seurat)
library(flexmix)
library(ggplot2)
library(patchwork)
library(SeuratObject)
options(future.globals.maxSize = 1e13)
options(Seurat.object.assay.version = "v5")
setwd("/media/simon/SC1")
data=readRDS("data_umap.rds")
# calculate cluster marker genes
# previous found too many clusters, now decrease resolution
data=FindClusters(data,resolution =0.2)
gc()
#only run findmarkers once to save time
#markers=FindAllMarkers(data, only.pos = F, min.pct =0.7, logfc.threshold = 1)
#write.csv(markers,"cluster_markers.csv")
#
#data[["RNA"]]$scale.data=NULL
for (i in 0:length(table(data$seurat_clusters))-1) {
  print(paste0("Cluster_",i))
  print(length(data$histo[which(data$seurat_clusters==i)]))
  print(table(data$histo[which(data$seurat_clusters==i)]))
}
# now load annotiation file to annotiate clusters
anno=read.csv("cluster_annotiation.csv")
Idents(data)=data$seurat_clusters
ids=as.character(anno$name)
names(ids)=levels(data)
data=RenameIdents(data,ids)
data$cell_type=Idents(data)
Idents(data)=data$cell_type
data=subset(x=data,idents=c("Low_quality","SV"),invert=T)
data=subset(data, subset =nCount_RNA>1000)
data=RunUMAP(data,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
data=FindNeighbors(data, reduction = "umap",dims = 1:2)
data=FindClusters(data, resolution = 0.3, cluster.name = "clusters")
saveRDS(data,"data_annotiated.rds")
#markers=FindAllMarkers(data, only.pos = TRUE, min.pct = 0.7, logfc.threshold = 1)
#write.csv(markers,"markers_new_lcuster.csv")
# now labe cells with further annotiation
anno2=read.csv("annotiation_for_new_clusters.csv")
ids=as.character(anno2$name)
names(ids)=levels(data)
data=RenameIdents(data,ids)
data$cell_type2=Idents(data)
# now only subset cancer associated clusters to annotiate the populations
adpc=subset(data,idents=c("Spink1","Krt7","AdPCa","NEPCa","Proliferating","Pou2f3"),invert=F)
adpc=RunUMAP(adpc,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
adpc=FindNeighbors(adpc,reduction = "umap",dims = 1:2)
adpc=FindClusters(adpc,resolution = 0.3)
markers_adpc=FindAllMarkers(adpc, only.pos = TRUE, min.pct = 0.7, logfc.threshold = 1)
write.csv(markers_adpc,"AdPCa_markers.csv")
for (i in 0:length(table(adpc$seurat_clusters))-1) {
  print(paste0("Cluster_",i))
  print(length(adpc$histo[which(adpc$seurat_clusters==i)]))
  print(table(adpc$histo[which(adpc$seurat_clusters==i)]))
}
# now annotiate AdPC populations
Idents(adpc)=adpc$seurat_clusters
anno3=read.csv("annotiation_for_adpc.csv")
ids=as.character(anno3$name)
names(ids)=levels(adpc)
adpc=RenameIdents(adpc,ids)
adpc$cell_type3=Idents(adpc)
# now transfer the label to whole populations
data$cell_type3=as.character(data$cell_type2)
data$cell_type3[rownames(adpc@meta.data)]=as.character(adpc$cell_type3)
Idents(data)=data$cell_type3
gc()
saveRDS(data,"data_annotiated_refined.rds")
xtabs(~ cell_type3 + histo,data = data@meta.data)
