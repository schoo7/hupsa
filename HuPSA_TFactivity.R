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
options(future.globals.maxSize = 1e13)
options(Seurat.object.assay.version = "v5")
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration")
data=readRDS("data_annotiated_refined_pathwayscored.rds")
network=decoupleR::get_progeny(organism = "human")
data_small=subset(x = data, downsample = 1500)
rm(data)
gc()
activities=decoupleR::run_wmean(mat = as.matrix(data_small[["RNA"]]$data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)
out=SCpubr::do_PathwayActivityPlot(sample = data_small,
                                      activities = activities)
tiff("pathway_analysis.tiff",width = 20,height = 40,units = "cm",res = 1200,compression = "lzw+p")
out$heatmaps$average_scores
dev.off()
network=decoupleR::get_dorothea(organism = "human",
                                   levels = c("A", "B"))
activities=decoupleR::run_wmean(mat =as.matrix(data_small[["RNA"]]$data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 100,
                                   minsize = 5)
out=SCpubr::do_TFActivityPlot(sample = data_small,
                                 activities = activities,n_tfs = 200)
tiff("TF_activity.tiff",width = 400,height = 40,units = "cm",res = 600,compression = "lzw+p")
out$heatmaps$average_scores
dev.off()
network=decoupleR::get_dorothea(organism = "human",
                                levels = c("A", "B","C"))
activities=decoupleR::run_wmean(mat =as.matrix(data_small[["RNA"]]$data),
                                network = network,
                                .source = "source",
                                .targe = "target",
                                .mor = "mor",
                                times = 100,
                                minsize = 5)
out=SCpubr::do_TFActivityPlot(sample = data_small,
                              activities = activities,n_tfs = 400,max.cutoff = 3)
write.csv(out[["heatmaps"]][["average_scores"]]@ht_list[["TF activity"]]@matrix,"TFactivity.csv")

