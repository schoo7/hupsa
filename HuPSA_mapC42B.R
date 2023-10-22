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
library(UCell)
library(GSEABase)
library(escape)
library(SCpubr)
library(decoupleR)
library(ggsankey)
library(extrafont)
library(extrafontdb)
options(Seurat.object.assay.version = "v5")
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration")
hupsa=readRDS("data_annotiated_refined_pathwayscored.rds")
hupsa=RunUMAP(hupsa,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap",return.model = T)
data=Read10X_h5("C42B.h5")
data=CreateSeuratObject(counts = data)
data$group="C42B_Xeno"
gc()
data=subset(data, nCount_RNA > 200 & nFeature_RNA < 5000)
data=NormalizeData(data)
data=FindVariableFeatures(data)
data=ScaleData(data)
data=RunPCA(data)
data=FindNeighbors(data, dims = 1:13)
data=FindClusters(data, resolution = 0.5)
data=RunUMAP(data, dims = 1:30)
SCpubr::do_DimPlot(sample = data,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,label.color = "black")
ggsave("C42B_unprojected_withlabel.tiff",width = 30,height = 25,units = "cm",dpi = 1200,compression = "lzw+p")



# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
gc()
anchor <- FindTransferAnchors(
  reference = hupsa,
  query = data,
  reference.reduction = "integrated.rpca",
  normalization.method = "LogNormalize",
  dims = 1:50
)
gc()
data=MapQuery(
  anchorset = anchor,
  query = data,
  reference = hupsa,
  reference.reduction = "integrated.rpca",
  refdata = list(
    cell_type3 = "cell_type3"
  ),
  reduction.model = "umap"
)
notch=list(Notch=c("CD44","DTX1","EPHB3","HES1","HES4","HES5","HES7","HEY1","HEY2","HEYL","MYC","NFKB2","NOX1","NRARP","PBX1","PIN1","PLXND1","SOX9"))
data=AddModuleScore_UCell(data, features = notch,maxRank = 3000)
saveRDS(data,"Projected_C42B.rds")
# start to visulize
Idents(data)=data$predicted.cell_type3
tiff("Projected_C42B_dimplot_withoutlabel.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_DimPlot(sample = data,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,label.color = "black")
dev.off()
tiff("Projected_C42B_dimplot_withlabel.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_DimPlot(sample = data,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,label.color = "black")
dev.off()
tiff("Projected_C42B_INSM1_NebulosaPlot.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_NebulosaPlot(data,"INSM1")
dev.off()
tiff("Projected_C42B_DLL3_NebulosaPlot.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_NebulosaPlot(data,"DLL3")
dev.off()
tiff("Projected_C42B_CHGA_NebulosaPlot.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_NebulosaPlot(data,"CHGA")
dev.off()
data_small=subset(x=data,idents=c("AdPCa_AR+","NEPCa"))
tiff("Projected_C42B_small_dimplot_withoutlabel.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_DimPlot(sample = data_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,label.color = "black")
dev.off()
tiff("Projected_C42B_small_dimplot_withlabel.tif",width = 30,height = 25,units = "cm",res = 1200,compression = "lzw+p")
SCpubr::do_DimPlot(sample = data_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,label.color = "black")
dev.off()




SCpubr::do_DimPlot(sample = hupsa,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,label.color = "black")



