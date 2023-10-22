# establish human prostate bulk RNAseq atlas
library(recount3)
library(CSCDRNA)
library(stringr)
library(SeuratObject)
library(reshape2)
library(ggplot2)
library(singscore)
library(preprocessCore)
library(pheatmap)
library(fst)
options(scipen = 999)
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration")
tcga=recount3::create_rse_manual(
  project = "PRAD",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
meta=colData(tcga)
meta=meta@listData
assay(tcga, "counts")=transform_counts(tcga)
count=tcga@assays@data$counts
tpm=recount::getTPM(tcga)
gene_anno=tcga@rowRanges@elementMetadata@listData
identical(rownames(count),gene_anno$gene_id)
count=as.data.frame(count)
count$symbol=gene_anno$gene_name
index=duplicated(count$symbol)
count=count[!index,]
tpm=as.data.frame(tpm)
tpm=tpm[rownames(count),]
rownames(count)=count$symbol
rownames(tpm)=rownames(count)
# now load previous re-aligned SRA data
count_pre=read.csv("PCa_patients_atlas_gene_counts.csv")
rownames(count_pre)=count_pre$X
count_pre=count_pre[,-1]
rowname=intersect(rownames(count_pre),rownames(count))
count_pre=count_pre[rowname,]
count=count[rowname,]
# now change column names
tcga_anno=as.data.frame(matrix(nrow = 558))
tcga_anno$ID_ori=as.character(meta[["external_id"]])
tcga_anno$group=as.character(meta[["tcga.cgc_sample_sample_type"]])

tcga_anno[, "group"]=ifelse(tcga_anno[, "group"] %in% c("Solid Tissue Normal"),
                               "normal", tcga_anno[, "group"])
tcga_anno[, "group"]=ifelse(tcga_anno[, "group"] %in% c("Primary Tumor","Metastatic"),
                            "AdPCa", tcga_anno[, "group"])
table(tcga_anno$group)
tcga_anno=tcga_anno[,-1]
start=713
n=558
tcga_anno$ID=paste0("Sample", start:(start+n-1))
count=count[,tcga_anno$ID_ori]
colnames(count)=tcga_anno$ID
tpm=tpm[,tcga_anno$ID_ori]
colnames(tpm)=tcga_anno$ID
# Now combine counts
identical(rownames(count_pre),rownames(count))
hPro_count=cbind(count_pre,count)
tpm_pre=read.csv("PCa_patients_atlas_gene_tpm.csv")
rownames(tpm_pre)=tpm_pre$X
tpm_pre=tpm_pre[,-1]
tpm_pre=tpm_pre[rowname,]
tpm=tpm[rowname,]
identical(rownames(tpm_pre),rownames(tpm))
hPro_tpm=cbind(tpm_pre,tpm)
saveRDS(hPro_count,"hProAtlas_withTCGA_count.rds")
saveRDS(hPro_tpm,"hProAtlas_withTCGA_tpm.rds")
write.csv(tcga_anno,"TCGA_sample_anno.csv")
# now start to decomposition
sc=readRDS("data_for_inferCNV.rds") # this is a subset of HuPSA
table(Idents(sc))
hPro_count_matrix=as.matrix(hPro_count)
bulk.eset=Biobase::ExpressionSet(assayData =hPro_count_matrix)
sc.counts.matrix=as.matrix(LayerData(sc,layer = "counts"))
sample.ids=colnames(sc.counts.matrix)
individual.labels=sc@meta.data$sample
cell.type.labels=sc@meta.data$cell_type3
cell.type.labels=ifelse(cell.type.labels %in% c("AdPCa_AR-","AdPCa_AR+","AdPCa_AR+_2","AdPCa_AR+_3","AdPCa_ARHigh","AdPCa_Proliferating","Normal_AR-"),
                               "AdPCa", cell.type.labels)
cell.type.labels=ifelse(cell.type.labels %in% c("T_CD4","T_CD8","T_proliferating"),
                        "T", cell.type.labels)
cell.type.labels=ifelse(cell.type.labels %in% c("Dendritic","Monocyte","Macrophage"),
                        "Myeloid", cell.type.labels)
cell.type.labels=ifelse(cell.type.labels %in% c("Fibroblast","Pericyte","SM"),
                        "Mesenchymal", cell.type.labels)

sc.pheno=data.frame(check.names=FALSE, check.rows=FALSE,
                       stringsAsFactors=FALSE,row.names=sample.ids,
                       SubjectName=individual.labels,cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),
                      row.names=c("SubjectName","cellType"))
sc.pdata=new("AnnotatedDataFrame",data=sc.pheno, varMetadata=sc.meta)
sc.eset=Biobase::ExpressionSet(assayData=sc.counts.matrix,phenoData=sc.pdata)
options(Seurat.object.assay.version = "v3")
#######
# need to remove.package() to all seurat packages
#reinstall the seurat V3
analysis=CSCD(bulk.eset=bulk.eset,sc.eset= sc.eset,
                 min.p=0.3,markers=NULL,cell.types="cellType",
                 subj.names="SubjectName",verbose=TRUE)
prob=as.data.frame(t(analysis$bulk.props))
# now start to assign groups based on decomposition result and markers enrichments
#read in samples histology data
histo=read.csv("hProAtlas_raw_anno.csv")
rownames(histo)=histo$X
histo=histo[,-1]
identical(rownames(histo),rownames(prob))
prob=cbind(prob,histo)
identical(rownames(prob),colnames(hPro_tpm))
rankData=rankGenes(hPro_tpm)
ne_score=simpleScore(rankData = rankData,upSet = c("INSM1","ASCL1","CHGA","SOX2","NCAM1","EZH2","SYP","DLL3"))
identical(rownames(prob),rownames(ne_score))
prob$NE_score=ne_score$TotalScore
# use MsigDB androgen response signature
AR_msigdb=c("ABCC4","ABHD2","ACSL3","ACTN1","ADAMTS1","ADRM1","AKAP12","AKT1","ALDH1A3","ANKH","APPBP2","ARID5B","AZGP1","B2M","B4GALT1","BMPR1B","CAMKK2","CCND1","CCND3","CDC14B",
            "CDK6","CENPN","DBI","DHCR24","DNAJB9","ELK4","ELL2","ELOVL5","FADS1","FKBP5","GNAI3","GPD1L","GSR","GUCY1A1","H1-0","HERC3","HMGCR","HMGCS1","HOMER2","HPGD","HSD17B14",
            "IDI1","INPP4B","INSIG1","IQGAP2","ITGAV","KLK2","KLK3","KRT19","KRT8","LIFR","LMAN1","MAF","MAK","MAP7","MERTK","MYL12A","NCOA4","NDRG1","NGLY1","NKX3-1","PA2G4","PDLIM5",
            "PGM3","PIAS1","PLPP1","PMEPA1","PTK2B","PTPN21","RAB4A","RPS6KA3","RRP12","SAT1","SCD","SEC24D","SELENOP","SGK1","SLC26A2","SLC38A2","SMS","SORD","SPCS3","SPDEF","SRF","SRP19","STEAP4",
            "STK39","TARP","TMEM50A","TMPRSS2","TNFAIP8","TPD52","TSC22D1","UAP1","UBE2I","UBE2J1","VAPA","XRCC5","XRCC6","ZBTB10","ZMIZ1")
AR_msigdb_score=simpleScore(rankData = rankData,upSet =AR_msigdb)
identical(rownames(prob),rownames(AR_msigdb_score))
prob$AR_score_msigdb=AR_msigdb_score$TotalScore
# Also need to add AR expression TPM 
ar_expression=as.data.frame(t(hPro_tpm["AR",rownames(prob)]))
identical(rownames(ar_expression),rownames(prob))
prob$AR_expression=ar_expression$AR
# Now calculate prostatic HOX code correlation
hox=read.csv("refined HOX code based on HOX pattern.csv")
hox=na.omit(hox)
rownames(hox)=hox$ID
hox_tpm=hPro_tpm[rownames(hox),]
identical(rownames(hox_tpm),rownames(hox))
hox_cor=vector(length = 877)
for (n in 1:877) {
  tem1=as.numeric(hox_tpm[,n])
  tem2=as.numeric(hox$prostate)
  tem3=cor.test(as.numeric(tem1),as.numeric(tem2))
  hox_cor[n]=tem3$estimate
}

hox_cor=as.data.frame(hox_cor)
rownames(hox_cor)=colnames(hox_tpm)
hox_cor=hox_cor[rownames(prob),]
prob$Hox_code=hox_cor
# Now assign sample groups
# First build a new column 
prob$class="unclassified"
prob[which(prob$group=="normal"),"class"]="Normal"
prob[which(prob$group=="BPH"),"class"]="Normal"
prob[which((prob$AdPCa>0.1)&(prob$class=="unclassified")),"class"]="AdPCa"
prob[which((prob$AR_expression>30)&(prob$class=="AdPCa")),"class"]="AdPCa_ARHigh"
prob[which((prob$AR_expression<=30)&(prob$class=="AdPCa")),"class"]="AdPCa_ARlow"
# Now identify NEPCa samples
prob[which((prob$NEPCa>0.1)&(prob$NE_score>0.2)),"class"]="NEPCa"
prob[which((prob$class=="NEPCa")&(prob$AR_expression>10)),"class"]="NEPCa_AdPCa_mix"
# now identify Progenitor_lie
prob[which(prob$Progenitor_like>0.15),"class"]="Progenitor_like"
prob[which((prob$AdPCa>0.1)&(prob$class=="Progenitor_like")),"class"]="Progenitor_like_AdPCa_mix"
# now identify KRT7
prob[which((prob$AdPCa_KRT7>0.15)&(prob$class=="unclassified")),"class"]="KRT7"
prob[which((prob$AdPCa>0.1)&(prob$class=="KRT7")),"class"]="KRT7_AdPCa_mix"
prob[which((prob$AdPCa_KRT7>0.1)&(prob$class=="AdPCa_ARHigh")),"class"]="KRT7_AdPCa_mix"
prob[which((prob$AdPCa_KRT7>0.1)&(prob$class=="AdPCa_ARlow")),"class"]="KRT7_AdPCa_mix"
# Now identify non cancer contamination

prob[which(prob$Basal>0.2),"class"]="contaminate"
prob[which(prob$Club>0.2),"class"]="contaminate"
prob[which(prob$B>0.2),"class"]="contaminate"
prob[which(prob$Endothelial>0.3),"class"]="contaminate"
prob[which(prob$Mesenchymal>0.3),"class"]="contaminate"
prob[which(prob$Myeloid>0.3),"class"]="contaminate"
prob[which(prob$Mast>0.3),"class"]="contaminate"
prob[which(prob$NK>0.3),"class"]="contaminate"
prob[which(prob$Plasma>0.3),"class"]="contaminate"
prob[which(prob$T>0.3),"class"]="contaminate"

# Also these are many mix samples assigned into AdPCa
write.csv(prob,"Assigned_sample_groups.csv")
prob$class=factor(prob$class,levels = c("Normal","AdPCa_ARHigh","AdPCa_ARlow","KRT7","Progenitor_like","NEPCa","KRT7_AdPCa_mix","Progenitor_like_AdPCa_mix","NEPCa_AdPCa_mix","unclassified","contaminate"))
hm=prob[order(prob$class,decreasing = F),]
anno=hm[,c("NE_score","AR_score_msigdb","Hox_code","group","class")]
hm=as.data.frame(t(hm[,c(1:14)]))
hm=hm[c(1,2,4,5,13,10,8,6,14,11,3,12,7,9),]
tiff("hPCa_atlas_composition_heatmap.tiff",width = 120,height = 40,units = "cm",res = 800,compression = "lzw+p")
pheatmap(hm,cluster_cols = F,annotation_col = anno,gaps_col = c(67,360,570,579,586,629,630,631,632,799),
         cluster_rows = F,show_colnames = F,
         cellheight = 8,cellwidth = 1,breaks = seq(0,0.8,0.01),color = colorRampPalette(c("navy", "white", "firebrick3"))(25),
         na_col = "firebrick3",border_color = NA)
dev.off()
# Now normalize TPM value and then generate data for ShinyAPP
gene_tpm_norm=normalize.quantiles(as.matrix(hPro_tpm),keep.names = T)
gene_tpm_norm2=as.data.frame(t(gene_tpm_norm))
gene_tpm_norm2$ID=rownames(gene_tpm_norm2)
prob$ID=rownames(prob)
result=merge(gene_tpm_norm2,prob[,c("group","class","ID")],by="ID")
result=melt(result,id=c("group","class","ID"))
write_fst(result,"Input_human_gene.fst")



# generate non-normalized hPro_tpm
# I tried, it is better to not perform normalization on this dataset
result_raw=as.data.frame(t(hPro_tpm))
result_raw$ID=rownames(result_raw)
prob$ID=rownames(prob)
result_raw=merge(result_raw,prob[,c("group","class","ID")],by="ID")
result_raw=melt(result_raw,id=c("group","class","ID"))
write_fst(result_raw,"Input_human_gene_raw.fst")
save.image("Enviroment_for_deconvolution.RData")
