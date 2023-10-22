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
library(decoupleR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(stringr)
options(future.globals.maxSize = 1e13)
options(Seurat.object.assay.version = "v5")
# now first rename clusters on HuPSA and MoPSA
hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored.rds")

new.cluster.ids <- c("Plasma","AdPCa_Proliferating","Dendritic","Monocyte","AdPCa_AR+_1","Macrophage","Fibroblast","Pericyte","AdPCa_ARlo","T_proliferating",
                     "SM","T_CD8","Club","Endothelial","Normal","T_CD4","B","NK","AdPCa_ARhi","Mast","KRT7","AdPCa_AR+_2","Basal","AdPCa_AR+_3",
                     "Progenitor_like","NEPCa")
names(new.cluster.ids) <- levels(hupsa)
hupsa <- RenameIdents(hupsa, new.cluster.ids)
Idents(hupsa)=factor(Idents(hupsa),levels = c("Normal","AdPCa_AR+_1","AdPCa_AR+_2","AdPCa_AR+_3","AdPCa_ARhi","AdPCa_ARlo","AdPCa_Proliferating","Basal","Club","KRT7","Progenitor_like",
                                              "NEPCa","Fibroblast","SM","Pericyte","Endothelial","T_CD4","T_CD8","T_proliferating","NK","B","Plasma","Mast","Monocyte","Macrophage","Dendritic"))
saveRDS(hupsa,"/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE.rds")
new.cluster.ids <- c("LE/AdPCa_1","SM","Basal","AdPCa_2","Schwann","Fibroblast","Endothelial","Pericyte","Spink1","AdPCa_Proliferating","Krt7",
                     "Macrophage","Monocyte","B","Erythroblasts","T","Dendritic","Neutrophils","NEPCa","AdPCa_3","Pou2f3")
names(new.cluster.ids) <- levels(mopsa)
mopsa <- RenameIdents(mopsa, new.cluster.ids)
Idents(mopsa)=factor(Idents(mopsa),levels = c("LE/AdPCa_1","AdPCa_2","AdPCa_3","Spink1","AdPCa_Proliferating","Basal","Krt7","NEPCa",
                                              "Pou2f3","Fibroblast","SM","Pericyte","Endothelial","T","Neutrophils",
                                              "B","Monocyte","Macrophage","Dendritic","Erythroblasts","Schwann"))
saveRDS(mopsa,"/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")

# Do dot plot
hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")

hmarker=list("Luminal"=c("FOXA1","HOXB13","NKX3-1","AR","KLK3"),
             "Proliferation"=c("TOP2A","MKI67"),
             "Basal"=c("KRT5","KRT14","TP63"),
             "Club"=c("MMP7","WFDC2"),
             "KRT7"=c("KRT7","S100A9"),
             "Progenitor"=c("SOX2","FOXA2"),
             "NEPCa"=c("NCAM1","INSM1"),
             "Fibroblast"=c("LUM","COL1A2","PDGFRA"),
             "SM"=c("TAGLN","ACTA2"),
             "Pericyte"=c("RGS5","NDUFA4L2"),
             "Endothelial"=c("FLT1","CALCRL"),
             "T_CD4"=c("IL7R"),
             "T_CD8"=c("CCL5","CCL4"),
             "NK"=c("GNLY","NKG7"),
             "B"=c("CD83","CD37"),
             "Plasma"=c("IGLL5","IGJ"),
             "Monocyte"=c("EREG","AQP9"),
             "Macrophage"=c("C1QB","C1QC"),
             "Dendritic"=c("CSF2RA","CD86"))
p=SCpubr::do_DotPlot(sample = hupsa, 
                     features = hmarker, 
                     cluster.idents = F, 
                     plot.title = "HuPSA")
ggsave("hupsa_marker_dot_plot.tiff",plot = p,width = 55,height = 25,units = "cm",dpi = 1200,compression = "lzw")
gc()
mmarker=list("Luminal"=c("Foxa1","Hoxb13","Nkx3-1","Ar","Pbsn"),
             "Proliferation"=c("Top2a","Mki67"),
             "Basal"=c("Krt5","Krt14","Trp63"),
             "KRT7"=c("Krt7","Clu"),
             "NEPCa"=c("Ncam1","Insm1"),
             "Pou2f3"=c("Pou2f3","Rgs13"),
             "Fibroblast"=c("Mgp","Col1a2","Pdgfra"),
             "SM"=c("Tagln","Acta2"),
             "Pericyte"=c("Rgs5","Ndufa4l2"),
             "Endothelial"=c("Flt1","Calcrl"),
             "T"=c("Cd3g","Trbc2"),
             "B"=c("Igkc","Ighm"),
             "Monocyte"=c("Lyz2","Ccl9"),
             "Macrophage"=c("C1qb","C1qc"),
             "Dendritic"=c("Cd74","Cd86"),
             "Erythroblasts"=c("Hba-a1","Hbb-b1"),
             "Schwann"=c("Csmd1","Plp1"))
p=SCpubr::do_DotPlot(sample = mopsa, 
                     features = mmarker, 
                     cluster.idents = F, 
                     plot.title = "MoPSA")
ggsave("mopsa_marker_dot_plot.tiff",plot = p,width = 60,height = 20,units = "cm",dpi = 1200,compression = "lzw")
# analysis AdPCa populations
hupsa_ad=subset(x=hupsa,idents=c("Normal","AdPCa_AR+_1","AdPCa_AR+_2","AdPCa_AR+_3","AdPCa_ARhi","AdPCa_ARlo","AdPCa_Proliferating"))
mopsa_ad=subset(x=mopsa,idents=c("LE/AdPCa_1","AdPCa_2","AdPCa_3","Spink1","AdPCa_Proliferating"))
rm(hupsa)
rm(mopsa)
hupsa_ad=RunUMAP(hupsa_ad,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
mopsa_ad=RunUMAP(mopsa_ad,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")

p1=SCpubr::do_DimPlot(sample = hupsa_ad,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p2=SCpubr::do_DimPlot(sample = mopsa_ad,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p=p1 | p2
p
ggsave("AdPCa_clusters_umap.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200,compression = "lzw")
genes=c("KLK3","AR","FOXA1","NKX3-1","HOXB13")
p1=SCpubr::do_DotPlot(sample = hupsa_ad, 
                      features = genes, 
                      cluster.idents = T, 
                      legend.position = "left")
genes_m=c("Pbsn","Ar","Foxa1","Nkx3-1","Hoxb13")
p2=SCpubr::do_DotPlot(sample = mopsa_ad, 
                      features = genes_m, 
                      cluster.idents = T, 
                      legend.position = "left")

p <- p1 | p2
ggsave("Genexpression_dot.tiff",plot = p,width = 40,height = 20,units = "cm",dpi = 1200,compression = "lzw")
hupsa_ad$cell_type3=Idents(hupsa_ad)
p1=SCpubr::do_BarPlot(hupsa_ad, 
                   group.by = "histo", 
                   split.by = "cell_type3",
                   plot.title = "HuPSA",
                   position = "fill")
mopsa_ad$cell_type3=Idents(mopsa_ad)
p2=SCpubr::do_BarPlot(mopsa_ad, 
                      group.by = "histo", 
                      split.by = "cell_type3",
                      plot.title = "MoPSA",
                      position = "fill")
p=p1|p2
ggsave("AdPCa_clusters_histo_distribution.tiff",plot = p,width = 22,height = 30,units = "cm",dpi = 1200,compression = "lzw")

# calculate markers for AdPCa subclusters
ad_h_markers=FindAllMarkers(hupsa_ad,min.pct = 0.7,logfc.threshold = 1)
write.csv(ad_h_markers,"AdPCa_clusters_markers_HuPSA.csv")
gc()
ad_m_markers=FindAllMarkers(mopsa_ad,min.pct = 0.7,logfc.threshold = 1)
write.csv(ad_m_markers,"AdPCa_clusters_markers_MoPSA.csv")


p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                      features = "KLK3",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p2=SCpubr::do_ViolinPlot(sample = mopsa_ad,
                         features = "Pbsn",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p=p1|p2
ggsave("KLK3_in_AdPCa.tiff",plot = p,width = 22,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "AR",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p2=SCpubr::do_ViolinPlot(sample = mopsa_ad,
                         features = "Ar",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p=p1|p2
ggsave("AR_in_AdPCa.tiff",plot = p,width = 22,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "FOXA1",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p2=SCpubr::do_ViolinPlot(sample = mopsa_ad,
                         features = "Foxa1",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p=p1|p2
ggsave("FOXA1_in_AdPCa.tiff",plot = p,width = 22,height = 8,units = "cm",dpi = 800)
gc()

p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "NKX3-1",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p2=SCpubr::do_ViolinPlot(sample = mopsa_ad,
                         features = "Nkx3-1",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p=p1|p2
ggsave("NKX3-1_in_AdPCa.tiff",plot = p,width = 22,height = 8,units = "cm",dpi = 800)
gc()

p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "HOXB13",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p2=SCpubr::do_ViolinPlot(sample = mopsa_ad,
                         features = "Hoxb13",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
p=p1|p2
ggsave("HOXB13_in_AdPCa.tiff",plot = p,width = 22,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "NPY",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
ggsave("NPY_in_AdPCa.tiff",plot = p1,width = 11,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "COL1A1",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
ggsave("COL1A1_in_AdPCa.tiff",plot = p1,width = 11,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "TAGLN",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
ggsave("TAGLN_in_AdPCa.tiff",plot = p1,width = 11,height = 8,units = "cm",dpi = 800)
gc()
p1=SCpubr::do_ViolinPlot(sample = hupsa_ad,
                         features = "ACTA2",legend.position = "none",pt.size = 0.1,plot_boxplot =F)
ggsave("ACTA2_in_AdPCa.tiff",plot = p1,width = 11,height = 8,units = "cm",dpi = 800)
gc()

# now visulize NEPCa clusters
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
hupsa_ne=subset(x=hupsa,idents="NEPCa")
mopsa_ne=subset(x=mopsa,idents="NEPCa")
rm(hupsa)
rm(mopsa)
hupsa_ne=RunUMAP(hupsa_ne,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
mopsa_ne=RunUMAP(mopsa_ne,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
p1=SCpubr::do_DimPlot(sample = hupsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "histo")
p2=SCpubr::do_DimPlot(sample = mopsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "histo")
p=p1 | p2
p
ggsave("NEPCa_clusrter_histo.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)
# now human NEPCa has 2 significant clusters
# identify clusters and annotiate
hupsa_ne=FindNeighbors(hupsa_ne,reduction = "umap",dims = 1:2)
hupsa_ne=FindClusters(hupsa_ne,resolution = 0.1)
DimPlot(hupsa_ne)
hupsa_ne_markers <- FindAllMarkers(hupsa_ne, only.pos = T, min.pct = 0.7, logfc.threshold = 1)
new.cluster.ids <- c("CHGA","CHGA","NEUROD1","CHGA","CHGA","TOX","STAT3","TOX","G2M")
names(new.cluster.ids) <- levels(hupsa_ne)
hupsa_ne <- RenameIdents(hupsa_ne, new.cluster.ids)
hupsa_ne_markers <- FindAllMarkers(hupsa_ne, only.pos = T, min.pct = 0.7, logfc.threshold = 1)
write.csv(hupsa_ne_markers,"HuPSA_NE_subclusters_markers.csv")

# predict pathway activity
network <- decoupleR::get_progeny(organism = "human")
activities=decoupleR::run_wmean(mat = as.matrix(hupsa_ne[["RNA"]]$data),
                                network = network,
                                .source = "source",
                                .targe = "target",
                                .mor = "weight",
                                times = 100,
                                minsize = 5)
out=SCpubr::do_PathwayActivityPlot(sample = hupsa_ne,
                                   activities = activities)
tiff("HuPSA_NEPCa_subclusters_pathway_activity.tiff",width = 20,height = 10,units = "cm",res = 1200)
out$heatmaps$average_scores
dev.off()
network=decoupleR::get_dorothea(organism = "human",
                                levels = c("A", "B"))
activities=decoupleR::run_wmean(mat =as.matrix(hupsa_ne[["RNA"]]$data),
                                network = network,
                                .source = "source",
                                .targe = "target",
                                .mor = "mor",
                                times = 100,
                                minsize = 5)
out=SCpubr::do_TFActivityPlot(sample = hupsa_ne,
                              activities = activities,n_tfs = 40)
tiff("HuPSA_NEPCa_subclusters_TFactivity.tiff",width = 60,height = 10,units = "cm",res = 1200)
out$heatmaps$average_scores
dev.off()

# redo dimplot
p1=SCpubr::do_DimPlot(sample = hupsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p2=SCpubr::do_DimPlot(sample = mopsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p=p1 | p2
p
ggsave("NEPCa_clusters_umap.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)

# do feature plots for both hupsa and mopsa
p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "CHGA",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Chga",legend.position ="top")
p=p1 |p2
ggsave("CHGA_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "NEUROD1",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Neurod1",legend.position ="top")
p=p1 |p2
ggsave("NEUROD1_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "STAT3",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Stat3",legend.position ="top")
p=p1 |p2
ggsave("STAT3_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "TOX",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Tox",legend.position ="top")
p=p1 |p2
ggsave("TOX_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "TOP2A",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Top2a",legend.position ="top")
p=p1 |p2
ggsave("TOP2A_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "MKI67",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Mki67",legend.position ="top")
p=p1 |p2
ggsave("MKI67_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=do_FeaturePlot(hupsa_ne, features = "ASCL1",legend.position ="top")
p2=do_FeaturePlot(hupsa_ne, features = "NEUROD1",legend.position ="top")
p= p1 | p2

p1=SCpubr::do_FeaturePlot(hupsa_ne, features = "ASCL1",legend.position ="top")
p2=SCpubr::do_FeaturePlot(mopsa_ne, features = "Ascl1",legend.position ="top")
p=p1 |p2
ggsave("ASCL1_NEPCa_subclusters.tiff",plot = p,width = 36,height = 20,units = "cm",dpi = 1200)

p1=SCpubr::do_DimPlot(sample = hupsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "study")
p2=SCpubr::do_DimPlot(sample = mopsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "study")
p=p1 | p2
p
ggsave("NEPCa_clusters_umap_study.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)


p1=SCpubr::do_DimPlot(sample = hupsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "sample")
p2=SCpubr::do_DimPlot(sample = mopsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "sample")
p=p1 | p2
p
ggsave("NEPCa_clusters_umap_sample.tiff",plot = p,width = 35,height = 25,units = "cm",dpi = 1200)

p1=SCpubr::do_DimPlot(sample = hupsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "Phase")
p2=SCpubr::do_DimPlot(sample = mopsa_ne,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "Phase")
p=p1 | p2
p
ggsave("NEPCa_clusters_phase.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)



# now analysis KRT7 clusters
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
gc()
# now visulize KRT7 expression in clusters
p1=do_GeyserPlot(hupsa,features = "KRT7",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Krt7",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("KRT7_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "AR",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Ar",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("AR_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "KLK3",scale_type = "categorical",order_by_mean =F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Pbsn",scale_type = "categorical",order_by_mean =F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("KLK3_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "FOXA1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Foxa1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("FOXA1_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")

gc()
p1=do_GeyserPlot(hupsa,features = "KRT8",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Krt8",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("KRT8_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "KRT18",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Krt18",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("KRT18_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "KRT5",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Krt5",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("KRT5_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "AR_signaling_UCell",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "AR_signaling_UCell",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("AR_signaling_UCell_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "HALLMARK_ANDROGEN_RESPONSE_UCell",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "HALLMARK_ANDROGEN_RESPONSE_UCell",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("HALLMARK_ANDROGEN_RESPONSE_UCell_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "HOXB13",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Hoxb13",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("HOXB13_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "NKX3-1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Nkx3-1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("NKX3-1_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "TP63",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Trp63",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("TP63_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "MMP7",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Mmp7",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("MMP7_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "FOXA2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Foxa2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("FOXA2_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "SOX2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Sox2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("SOX2_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "INSM1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Insm1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("INSM1_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "CHGA",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Chga",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("CHGA_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
p1=do_GeyserPlot(hupsa,features = "NCAM1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Ncam1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("NCAM1_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
# river Sankeyplotp
hupsa$cell_type3=Idents(hupsa)
hupsa$cell_type3=factor(hupsa$cell_type3,levels = levels(Idents(hupsa)))
mopsa$cell_type3=Idents(mopsa)
p=SCpubr::do_SankeyPlot(sample = hupsa,
                         first_group = "cell_type3",
                         middle_groups = "histo",
                         last_group = "study",
                         type = "sankey")
ggsave("sankey_hupsa.tiff",plot = p,width = 22,height = 22,units = "cm",dpi = 1200,compression = "lzw")






hupsa_krt7=subset(x=hupsa,idents="AdPCa_KRT7")
mopsa_krt7=subset(x=mopsa,idents="AdPCa_Krt7")
rm(hupsa)
rm(mopsa)
hupsa_krt7=RunUMAP(hupsa_krt7,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")
mopsa_krt7=RunUMAP(mopsa_krt7,reduction = "integrated.rpca",dims = 1:30,reduction.name = "umap")

p1=SCpubr::do_DimPlot(sample = hupsa_krt7,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p2=SCpubr::do_DimPlot(sample = mopsa_krt7,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom")
p=p1 | p2
p
ggsave("KRT7_clusters_phase.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)

p1=SCpubr::do_DimPlot(sample = hupsa_krt7,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "histo")
p2=SCpubr::do_DimPlot(sample = mopsa_krt7,pt.size = 0.1,raster = F,border.size = 20,legend.position = "bottom",group.by = "histo")
p=p1 | p2
p
ggsave("KRT7_clusters_phase_histo.tiff",plot = p,width = 18,height = 10,units = "cm",dpi = 1200)

# bar plot to show the composition of KRT7 clusters
do_BarPlot(hupsa_krt7,group.by = "histo",position = "stack")
do_BarPlot(mopsa_krt7,group.by = "histo",position = "stack")

# Now analysis MMP7
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE.rds")
hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored.rds")
p1=do_GeyserPlot(hupsa,features = "WFDC2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Wfdc2",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("WFDC2_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()
hupsa_ne=subset(x=hupsa,idents="NEPCa")
mopsa_ne=subset(x=mopsa,idents="NEPCa")
rm(hupsa)
rm(mopsa)


# load TF activity files and visulize for some populations
tf_h=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/TFactivity.csv")
tf_m=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/TFactivity.csv")

rownames(tf_h)=tf_h$X
tf_h=tf_h[,c("SREBF1","TCF7L1","SREBF2","GRHL2","ASCL1","ONECUT1")]
p=pheatmap(t(tf_h),cluster_cols = F,breaks = seq(0,2,0.05),colorRampPalette(c("green", "black", "red"))(41),border_color = NA,cluster_rows = F,angle_col = 45,cellwidth = 15,cellheight = 15,fontsize = 14)
ggsave("Upstream_TF_for_progenotor_like.tiff",plot = p,width = 36,height = 5,dpi = 800,compression = "lzw")



# bubble plot to visulize marker genes expression
#hmarker=c("IGLL5","IGJ","KIAA0101","BIRC5","CSF2RA","CD86","EREG","AQP9","KLK3","KRT8","AR","C1QB","C1QC","LUM","COL1A2","RGS5","NDUFA4L2","TAGLN","ACTA2","CCL5","CCL4","MMP7","WFDC2","FLT1","CALCRL","CD83","GNLY","TPSAB1","KRT7","S100A9","EEF1A2","NPY","KRT5","TACSTD2","FOXA2","SOX2","CHGA","INSM1")

hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
hmarker=list("Plasma"=c("IGLL5","IGJ"),
             "Dendritic"=c("CSF2RA","CD86"),
             "Monocyte"=c("EREG","AQP9"),
             "AdPCa"=c("AR","KLK3","HOXB13","FOXA1","NKX3-1"),
             "Macrophage"=c("C1QB","C1QC"),
             "Fibroblast"=c("LUM","COL1A2"),
             "Pericyte"=c("RGS5","NDUFA4L2"),
             "SM"=c("TAGLN","ACTA2"),
             "T_CD8"=c("CCL5","CCL4"),
             "T_CD4"=c("IL7R"),
             "Club"=c("MMP7","WFDC2"),
             "Endothelial"=c("FLT1","CALCRL"),
             "B"=c("CD83","CD37"),
             "NK"=c("GNLY","NKG7"),
             "Mast"=c("TPSAB1","TPSB2"),
             "KRT7"=c("KRT7","S100A9"),
             "Basal"=c("KRT5","KRT14"),
             "Progenitor"=c("SOX2","FOXA2"),
             "NEPCa"=c("CHGA","INSM1"))
p=SCpubr::do_DotPlot(sample = hupsa, 
                   features = hmarker, 
                   cluster.idents = TRUE, 
                   plot.title = "Clustered")
ggsave("hupsa_marker_dot_plot.tiff",plot = p,width = 55,height = 20,units = "cm",dpi = 1200,compression = "lzw")

rm(hupsa)
mmarker=list("Dendritic"=c("Cd74","Cd86"),
             "Monocyte"=c("Lyz2","Ccl9"),
             "AdPCa"=c("Ar","Pbsn","Hoxb13","Foxa1","Nkx3-1"),
             "Macrophage"=c("C1qb","C1qc"),
             "Fibroblast"=c("Mgp","Col1a2"),
             "Pericyte"=c("Rgs5","Ndufa4l2"),
             "SM"=c("Tagln","Acta2"),
             "T"=c("Cd3g","Trbc2"),
             "B"=c("Igkc","Ighm"),
             "Neutrophils"=c("S100a8","S100a9"),
             "Krt7"=c("Krt7","Clu"),
             "Basal"=c("Krt5","Krt14"),
             "Schwann"=c("Csmd1","Plp1"),
             "Endothelial"=c("Flt1","Calcrl"),
             "Spink1"=c("Spink1","Sbp"),
             "Erythroblasts"=c("Hba-a1","Hbb-b1"),
             "Pou2f3"=c("Pou2f3","Rgs13"),
             "NEPCa"=c("Chga","Insm1"))
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
p=SCpubr::do_DotPlot(sample = mopsa, 
                     features = mmarker, 
                     cluster.idents = TRUE, 
                     plot.title = "Clustered")
ggsave("mopsa_marker_dot_plot.tiff",plot = p,width = 60,height = 20,units = "cm",dpi = 1200,compression = "lzw")



# castration and regression analysis
# GSE146811
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
mopsa$cell_type4=Idents(mopsa)
Idents(mopsa)=mopsa$study
mopsa=subset(x = mopsa, idents = "GSE146811")
Idents(mopsa)=mopsa$cell_type4
mopsa_small=subset(x = mopsa, idents = c("LE/AdPCa_1","Krt7"))
mopsa_small=RunUMAP(mopsa_small,reduction = "integrated.rpca",dims = 1:50)
mopsa_small$annotiation=gsub("Intact_EPCAM_negtive_sorted","Intact",mopsa_small$annotiation)
mopsa_small$annotiation=gsub("Intact_EPCAM_positive_sorted","Intact",mopsa_small$annotiation)
mopsa_small$annotiation=factor(mopsa_small$annotiation,levels = c("Intact","CX1D","CX7D","CX14D","CX28D","R1D","R2D","R3D","R7D","R14D","R28D"))
p1=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "right")
p2=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "right",group.by = "annotiation")
p=p1/p2
p
ggsave("KRT7_regression_umap.tiff",plot = p,width = 18,height = 22,units = "cm",dpi = 1200)
p=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,split.by = "annotiation")
ggsave("KRT7_regression_split_umap.tiff",plot = p,width = 36,height = 26,units = "cm",dpi = 1200)

# Enlarge scopre to non-basal epithelium cells
mopsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration/data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
mopsa$cell_type4=Idents(mopsa)
Idents(mopsa)=mopsa$study
mopsa=subset(x = mopsa, idents = "GSE146811")
Idents(mopsa)=mopsa$cell_type4
mopsa_small=subset(x = mopsa, idents = c("LE/AdPCa_1","Krt7","AdPCa_2","Spink1","AdPCa_Proliferating"))
mopsa_small=RunUMAP(mopsa_small,reduction = "integrated.rpca",dims = 1:50)
mopsa_small$annotiation=gsub("Intact_EPCAM_negtive_sorted","Intact",mopsa_small$annotiation)
mopsa_small$annotiation=gsub("Intact_EPCAM_positive_sorted","Intact",mopsa_small$annotiation)
mopsa_small$annotiation=factor(mopsa_small$annotiation,levels = c("Intact","CX1D","CX7D","CX14D","CX28D","R1D","R2D","R3D","R7D","R14D","R28D"))
mopsa_small$annotiation2=mopsa_small$annotiation
p1=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "right")
p2=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "right",group.by = "annotiation2")
p=p1/p2
p
ggsave("KRT7_regression_umap_new.tiff",plot = p,width = 18,height = 22,units = "cm",dpi = 1200,compression = "lzw")
p=SCpubr::do_DimPlot(sample = mopsa_small,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,split.by = "annotiation")
ggsave("KRT7_regression_split_umap_new.tiff",plot = p,width = 36,height = 26,units = "cm",dpi = 1200,compression = "lzw")



p3=SCpubr::do_BarPlot(mopsa_small,
                      group.by = "annotiation",
                      split.by = "cell_type4",
                      plot.title = "Number of cells per cluster in each sample",
                      position = "stack")
p3




# tem
p1=do_GeyserPlot(hupsa,features = "ABCC1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Abcc1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("ABCC1_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

p1=do_GeyserPlot(hupsa,features = "CASR",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p2=do_GeyserPlot(mopsa,features = "Casr",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
p=p1 | p2
ggsave("CASR_geyserplot.tiff",plot = p,width = 36,height = 10,units = "cm",dpi = 1200,compression = "lzw")
gc()

##load mapped C42B

data=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/Projected_C42B.rds")
data_small=subset(x = data, idents = c("AdPCa_AR+","AdPCa_AR-","AdPCa_AR+_2","AdPCa_Proliferating","NEPCa"))
DimPlot(data_small)
do_GeyserPlot(data_small,features = "AR",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
ggsave("AR_C42B.tiff",width = 15,height = 8,units = "cm",dpi = 1200,compression = "lzw")
do_GeyserPlot(data_small,features = "CHGA",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
ggsave("CHGA_C42B.tiff",width = 15,height = 8,units = "cm",dpi = 1200,compression = "lzw")
do_GeyserPlot(data_small,features = "KLK3",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
ggsave("KLK3_C42B.tiff",width = 15,height = 8,units = "cm",dpi = 1200,compression = "lzw")
do_GeyserPlot(data_small,features = "TOP2A",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
ggsave("TOP2A_C42B.tiff",width = 15,height = 8,units = "cm",dpi = 1200,compression = "lzw")
do_GeyserPlot(data_small,features = "INSM1",scale_type = "categorical",order_by_mean = F,legend.position = "none",jitter = 0.45,pt.size = 0.1)
ggsave("INSM1_C42B.tiff",width = 15,height = 8,units = "cm",dpi = 1200,compression = "lzw")









