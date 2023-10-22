library(Seurat)
library(SeuratObject)
library(flexmix)
library(ggplot2)
library(patchwork)
library(sctransform)
library(SummarizedExperiment)
library(GSEABase)
library(escape)
library(SCpubr)
library(decoupleR)
library(infercnv)
library(dplyr)
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration")
data=readRDS("data_annotiated_refined_pathwayscored_combinedLE_renamed.rds")
data_small=subset(x = data, downsample = 2000)
data_small$cell_type3=Idents(data_small)
data_smallv3=as(object = data_small[["RNA"]], Class = "Assay")
rm(data)
gc()
meta=data_small[[]]
meta$cell_type_cnv=data_small$cell_type3
meta$cell_type_cnv <- as.character(meta$cell_type_cnv)
meta$cell_type_cnv[meta$cell_type_cnv %in% c("Basal", "Fibroblast", "SM", "Pericyte", "Endothelial",
                                             "T", "Neutrophils", "B", "Monocyte", "Macrophage", "Dendritic",
                                             "Erythroblasts", "Schwann")] <- "normal"
meta$cell_type_cnv=factor(meta$cell_type_cnv,levels = c("normal","Spink1","LE/AdPCa_1","AdPCa_2","AdPCa_3","AdPCa_Proliferating","Krt7","NEPCa","Pou2f3"))
meta=as.matrix(meta)
mat=as.matrix(data_smallv3@counts)
write.table(meta[,"cell_type_cnv"],"meta_for_cnv.txt",sep = "\t",row.names = TRUE,col.names = F)
write.table(mat,"matrix_for_cnv.txt",sep = "\t",row.names = TRUE,col.names = T)
# generate a new meta_for_CNV to separate normal and cancer cells
gc()
infercnv_obj=CreateInfercnvObject(raw_counts_matrix = "C:/Users/13831/OneDrive - LSU Health Shreveport/MoPSA/integration/matrix_for_cnv.txt",
                                  annotations_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/MoPSA/integration/meta_for_cnv.txt",
                                  gene_order_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/MoPSA/integration/GRCm38_gene_order.txt",
                                  delim = "\t",
                                  ref_group_names=c("normal"),
                                  max_cells_per_group = 500,
                                  chr_exclude =c("chrM","chrY"))
options(scipen = 100)
gc()
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="C:/Users/13831/OneDrive - LSU Health Shreveport/MoPSA/integration/inferCNV_new",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             cluster_references = F,
                             denoise=T,
                             num_threads=4,
                             output_format = "pdf",
                             useRaster = T,
                             HMM=T)






