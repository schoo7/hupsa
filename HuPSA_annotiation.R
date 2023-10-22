library(Seurat)
library(flexmix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(SeuratObject)
library(SeuratData)
library(SeuratWrappers)
library(sctransform)
library(DoubletFinder)
library(fishpond)
library(SummarizedExperiment)
library(Azimuth)
library(UCell)
library(GSEABase)
library(escape)
library(SCpubr)
options(future.globals.maxSize = 1e13)
options(Seurat.object.assay.version = "v5")
data=readRDS("data_umap.rds")
# perform cell cycle analysis
s.genes=cc.genes$s.genes
g2m.genes=cc.genes$g2m.genes
data=CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# calculate cluster marker genes
# previous found too many clusters, now decrease resolution
data=FindClusters(data,resolution =1)
Idents(data)=data$clusters
gc()
#only run findmarkers once to save time
#markers=FindAllMarkers(data, only.pos = F, min.pct =0.7, logfc.threshold = 0.7)
#write.csv(markers,"cluster_markers.csv")
#
#data[["RNA"]]$scale.data=NULL
for (i in 0:length(table(data$clusters))-1) {
  print(paste0("Cluster_",i))
  print(length(data$histo[which(data$clusters==i)]))
  print(table(data$histo[which(data$clusters==i)]))
}
# now load annotiation file to annotiate clusters
anno=read.csv("cluster_annotiation.csv")
Idents(data)=data$clusters
ids=as.character(anno$name)
names(ids)=levels(data)
data=RenameIdents(data,ids)
data$cell_type=Idents(data)
Idents(data)=data$cell_type
data=subset(x=data,idents="Low_quality",invert=T)
data=RunUMAP(data,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
saveRDS(data,"data_annotiated.rds")
# now use sc-type to annotiated populations
library(EasyCellType)
data=FindNeighbors(data,reduction = "umap",dims = 1:2)
data=FindClusters(data,resolution = 0.5)
#markers=FindAllMarkers(data, only.pos = TRUE, min.pct = 0.7, logfc.threshold = 1)
#write.csv(markers,"markers_new_lcuster.csv")
markers=read.csv("markers_new_lcuster.csv")
library("org.Hs.eg.db")
library(AnnotationDbi)
markers$entrezid <- mapIds(org.Hs.eg.db,
                           keys=markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
markers <- na.omit(markers)
library(dplyr)
markers_sort <- data.frame(gene=markers$entrezid, cluster=markers$cluster, 
                           score=markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(markers_sort[, 1:3])
annot.GSEA <- easyct(input.d, db="cellmarker", species="Human", 
                     tissue=c("Blood", "Peripheral blood","Prostate"), p_cut=0.3,
                     test="GSEA")
# using markers to define Tcell populations

list=list()
list$Dendritic_cells=c("CORO1A","MS4A6A","ITGB2","GPR183","HLA-DRB1","HLA-DPB1","HLA-DPA1","HLA-DQB1","HLA-DQA1","HLA-DMA")
list$Monocytes=c("FCER1G","S100A9","S100A12","VCAN","COTL1","S100A8","CD14","LST1","FCN1","AIF1")
list$NK=c("GZMA","CD7","CCL4","CST7","NKG7","GNLY","CTSW","CCL5","GZMB","PRF1")
list$CD8=c("CD8A","CD3E","CCL4","CD2","CXCR4","GZMA","NKG7","IL32","CD3D","CCL5")
list$CD4=c("CORO1A","KLRB1","CD3E","LTB","CXCR4","IL7R","TRAC","IL32","CD2","CD3D")
list$Macrophage=c("C1QB","C1QA","HLA-DRB1","AIF1","LYZ","CTSZ","CTSL","FCER1G","C1QC","LAPTM5")
data=AddModuleScore_UCell(data, features = list,maxRank = 3000)
for (i in 0:length(table(data$seurat_clusters))-1) {
  print(paste0("Cluster_",i))
  print(length(data$histo[which(data$seurat_clusters==i)]))
  print(table(data$histo[which(data$seurat_clusters==i)]))
}
# now labe cells with further annotiation
anno2=read.csv("annotiation_for_new_clusters.csv")
ids=as.character(anno2$name)
names(ids)=levels(data)
data=RenameIdents(data,ids)
data$cell_type2=Idents(data)
# now only subset AdPCa to annotiate the populations
adpc=subset(data,idents="AdPCa",invert=F)
adpc=RunUMAP(adpc,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
adpc=FindNeighbors(adpc,reduction = "umap",dims = 1:2)
adpc=FindClusters(adpc,resolution = 0.3)
gc()
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
# now run different signaling markers
ne=list(NEPCa_up = c("NCAM1","SOX2","CHGA","FOXA2","INSM1","SRRM4","ASCL1","SYP","EZH2"))
# AR signaling score https://doi.org/10.3389/fonc.2022.955166
ar=list(AR_signaling=c("KLK3", "KLK2", "STEAP4", "TMPRSS2", "FKBP5", "ALDH1A3", "NKX3-1", "PPAP2A", "PMEPA1","PART1"))
GS.hallmark=getGeneSets(library = "H")
hmark=geneIds(GS.hallmark)
hmark=hmark[c(3,4,7,12,13,14,19,20,22,23,24,25,28,29,25,37,38,44,49)]
sig=c(notch,ne,ar)
gc()
data=AddModuleScore_UCell(data, features = sig,maxRank = 3000)
gc()
data=AddModuleScore_UCell(data, features = hmark,maxRank = 3000)
gc()
Idents(data)=data$cell_type3
gc()
cell_marker=FindAllMarkers(data, only.pos = TRUE, min.pct = 0.7, logfc.threshold = 1.5)
write.csv(cell_marker,"Cell_type_markers.csv")
saveRDS(data,"data_annotiated_refined_pathwayscored.rds")