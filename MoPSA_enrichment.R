# add signaling activities to MoPSA
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
library(nichenetr)
library(UCell)
library(GSEABase)
library(escape)
library(SCpubr)
options(future.globals.maxSize = 1e15)
options(Seurat.object.assay.version = "v5")
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/MoPSA/integration")
data=readRDS("data_annotiated_refined.rds")
# calculate cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes=na.omit(convert_human_to_mouse_symbols(s.genes))
g2m.genes=na.omit(convert_human_to_mouse_symbols(g2m.genes))
data=CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
ne=list(NEPCa_up = c("NCAM1","SOX2","CHGA","FOXA2","INSM1","SRRM4","ASCL1","SYP","EZH2"))
# AR signaling score https://doi.org/10.3389/fonc.2022.955166
ar=list(AR_signaling=c("KLK3", "KLK2", "STEAP4", "TMPRSS2", "FKBP5", "ALDH1A3", "NKX3-1", "PPAP2A", "PMEPA1","PART1"))
GS.hallmark=getGeneSets(library = "H")
hmark=geneIds(GS.hallmark)
hmark=hmark[c(3,4,7,12,13,14,19,20,22,23,24,25,28,29,25,37,38,44,49)]
sig=c(ne,ar,hmark)
converted_sig <- vector("list", length(sig))
for (i in 1:length(sig)) {
  converted_sig[[i]]=na.omit(convert_human_to_mouse_symbols(sig[[i]]))
}
names(converted_sig)=names(sig)
gc()
data=AddModuleScore_UCell(data, features = converted_sig,maxRank = 2500)
gc()
Idents(data)=data$cell_type3
gc()
saveRDS(data,"data_annotiated_refined_pathwayscored.rds")
# The LE population is combined into AdPCa, their transcriptome is not significantlly different to AdPCa,
#because the AdPCa populations also contains a huge amount of LE
data$cell_type3=gsub("LE","AdPCa",data$cell_type3)
Idents(data)=data$cell_type3
saveRDS(data,"data_annotiated_refined_pathwayscored_combinedLE.rds")
