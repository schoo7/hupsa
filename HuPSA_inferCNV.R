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
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration")
data=readRDS("data_annotiated_refined_pathwayscored_renamed.rds")
data_small=subset(x = data, downsample = 2000)
data_smallv3=as(object = data_small[["RNA"]], Class = "Assay")
rm(data)
gc()
data_small$cell_type3=Idents(data_small)
meta=data_small[[]]
meta$cell_type_cnv=meta$cell_type3  
meta$cell_type_cnv=gsub("Normal","Normal_LE",meta$cell_type_cnv)
meta[, "cell_type_cnv"]=ifelse(meta[, "cell_type_cnv"] %in% c("B", "SM", "NK", "Endothelial",
                                                              "Fibroblast", "Plasma", "Dendritic",
                                                              "Monocyte", "Macrophage", "Pericyte",
                                                              "T_proliferating", "T_CD8","T_CD4", 
                                                              "Mast", "Basal"),
                               "normal", meta[, "cell_type_cnv"])
meta=as.matrix(meta)
mat=as.matrix(data_smallv3@counts)
tem=GetAssayData(data_small,slot = "data")
tem2=as(tem, "sparseMatrix")
write.table(meta[,"cell_type_cnv"],"meta_for_cnv.txt",sep = "\t",row.names = TRUE,col.names = F)
write.table(mat,"matrix_for_cnv.txt",sep = "\t",row.names = TRUE,col.names = T)
infercnv_obj=CreateInfercnvObject(raw_counts_matrix = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/matrix_for_cnv.txt",
                                     annotations_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/meta_for_cnv.txt",
                                     gene_order_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/hg38_gencode_v27.txt",
                                     delim = "\t",
                                     ref_group_names=c("normal"),
                                     max_cells_per_group = 500,
                                  chr_exclude =c("chrM","chrY"))
options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/inferCNV_new",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             cluster_references = F,
                             denoise=T,
                             num_threads=4,
                             output_format = "pdf",
                             useRaster = T,
                             HMM=T)




out <- SCpubr::do_CopyNumberVariantPlot(sample = data,
                                        infercnv_object = infercnv_object,
                                        using_metacells = TRUE,
                                        metacell_mapping = sample$metacell_mapping,
                                        chromosome_locations = human_chr_locations,
                                        chromosome_focus = NULL)

# only access CNV on AdPCa sub clusters
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration")
data=readRDS("data_annotiated_refined_pathwayscored_renamed.rds")
data=subset(x = data, idents = c("B", "SM", "NK","Endothelial","Fibroblast","Plasma","Dendritic","Monocyte","Macrophage","Pericyte",
                                                                               "T_proliferating", "T_CD8",
                                                                               "T_CD4", "Mast", "Basal","Normal","AdPCa_AR+_1","AdPCa_AR+_2","AdPCa_AR+_3","AdPCa_ARhi","AdPCa_ARlo","AdPCa_Proliferating"), invert = F)


data_small=subset(x = data, downsample = 1000)
data_smallv3=as(object = data_small[["RNA"]], Class = "Assay")
rm(data)
gc()
meta=data_small[[]]
meta$cell_type_cnv=meta$cell_type3  
meta[, "cell_type_cnv"]=ifelse(meta[, "cell_type_cnv"] %in% c("B", "SM", "NK", "Endothelial",
                                                              "Fibroblast", "Plasma", "Dendritic",
                                                              "Monocyte", "Macrophage", "Pericyte",
                                                              "T_proliferating", "T_CD8",
                                                              "T_CD4", "Mast", "Basal"),
                               "normal", meta[, "cell_type_cnv"])

meta=as.matrix(meta)
mat=as.matrix(data_smallv3@counts)
write.table(meta[,"cell_type_cnv"],"meta_for_cnv_AdPCa.txt",sep = "\t",row.names = TRUE,col.names = F)
write.table(mat,"matrix_for_cnv_AdPCa.txt",sep = "\t",row.names = TRUE,col.names = T)
infercnv_obj=CreateInfercnvObject(raw_counts_matrix = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/matrix_for_cnv_AdPCa.txt",
                                  annotations_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/meta_for_cnv_AdPCa.txt",
                                  gene_order_file = "C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/hg38_gencode_v27.txt",
                                  delim = "\t",
                                  ref_group_names=c("normal"),
                                  max_cells_per_group = 1000)
options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="C:/Users/13831/OneDrive - LSU Health Shreveport/HuPSA/HuPSA_integration/inferCNV_AdPCa_new",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             cluster_references = F,
                             denoise=F,
                             num_threads=8,
                             output_format = "pdf",
                             useRaster = T,
                             HMM=F)






