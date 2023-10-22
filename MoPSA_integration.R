# data integration for HuPSA, using Seurat V5
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
library(org.Mm.eg.db)
options(future.globals.maxSize = 1e15)
options(Seurat.object.assay.version = "v5")
#
sample_annotation <- read.csv("sample_annotiation.csv")
index=duplicated(sample_annotation$sample)
summary(index)
sample_annotation=sample_annotation[!index,]
# load samples
seurat_list <- list()
for (i in 1:length(sample_annotation$ID)) {
  folder_path <- paste0("sample", i, "_alevin_results/af_quant")
  print(paste0("Working on Sample",i))
  # Load the sample
  sample <- loadFry(folder_path, outputFormat = "snRNA")
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = assay(sample, "counts"), min.cells = 3, min.features = 200)
  # Assign metadata
  seurat_obj$sample <- sample_annotation$ID[i]
  seurat_obj$histo <- sample_annotation$Histo[i]
  seurat_obj$chemistry <- sample_annotation$Chemistry[i]
  seurat_obj$study <- sample_annotation$Study[i]
  seurat_obj$annotiation <- sample_annotation$Annotiation[i]
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &nCount_RNA>1000)
  print(ncol(seurat_obj))
  # Add the Seurat object to the list
  seurat_list[[i]] <- seurat_obj
  gc()
}
rm(list = setdiff(ls(), "seurat_list"))
gc()
saveRDS(seurat_list,"seuratlist.rds")
cell_counts <- sapply(seurat_list, function(obj) ncol(obj@assays$RNA))
indices_to_remove <- which(cell_counts < 500)
seurat_list <- seurat_list[-indices_to_remove]

# Merge all Seurat objects into a single object called "data"
data=merge(seurat_list[[1]], y =c(seurat_list[2:138]), project = "MoPSA")
rm(list = setdiff(ls(), "data"))
gc()
saveRDS(data,"merged.rds")
gc()
data=JoinLayers(data)
saveRDS(data,"Join.rds")
gc()
count=GetAssayData(data[["RNA"]],layer="counts")
gene_symbols <- mapIds(org.Mm.eg.db, rownames(count), "SYMBOL", "ENSEMBL")
filtered_rows <- !is.na(gene_symbols) & !duplicated(gene_symbols)
count <- count[filtered_rows, ]
filtered_gene_symbols <- gene_symbols[filtered_rows]
rownames(count) <- filtered_gene_symbols
MoPSA=CreateSeuratObject(counts = count)
MoPSA[[]]=data[[]]
rm(list = setdiff(ls(), "MoPSA"))
gc()
MoPSA[["RNA"]]=split(MoPSA[["RNA"]], f = MoPSA$study)
gc()
MoPSA=NormalizeData(MoPSA)
gc()
MoPSA=FindVariableFeatures(MoPSA)
gc()
MoPSA=ScaleData(MoPSA)
gc()
#saveRDS(MoPSA,"data_scaled.rds")
gc()
MoPSA=RunPCA(MoPSA)
gc()
#saveRDS(MoPSA,"data_afterPCA.rds")
gc()
# use RPCA integration for big dataset
MoPSA=IntegrateLayers(
  object = MoPSA, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = T
)
gc()
#saveRDS(MoPSA,"data_integrated.rds")
gc()
MoPSA=RunUMAP(MoPSA, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap")
gc()
MoPSA=FindNeighbors(MoPSA, reduction = "umap",dims = 1:2)
gc()
MoPSA=FindClusters(MoPSA, resolution = 2, cluster.name = "clusters")
gc()
MoPSA=JoinLayers(MoPSA)
gc()
saveRDS(MoPSA,"data_umap.rds")





