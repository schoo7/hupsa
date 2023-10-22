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
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1e13)
options(Seurat.object.assay.version = "v5")
data=readRDS("data_annotiated_refined_pathwayscored_combinedLE.rds")
# Now LR analysis
labels <- Idents(data)
data.input <- GetAssayData(data, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 50)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
tiff("Cell_interations_overview.tif",width = 20,height = 23,units = "cm",res = 600,compression = "lzw+p")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
mat <- cellchat@net$weight
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
tiff("Cell_interations_overview_heatmap.tif",width = 20,height = 40,units = "cm",res = 600,compression = "lzw+p")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 18,width = 7)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 18,width = 7)
ht1 + ht2
dev.off()
# Find all significant pathways
cellchat@netP$pathways
