library(plotly)
library(htmlwidgets)
library("rmarkdown")
library(Seurat)
library(monocle3)
library(scales)
library(igraph)
library(ggraph)
library(graphlayouts)
library(qgraph)
library(AUCell)
library(MetaNeighbor)
library(Matrix)
library(dendextend)

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

subplots <- function(sobj,colName,filename,raster=FALSE){
classes = unique(sobj@meta.data[colName])[,1]
pdf(filename)
for(i in classes){
sobj@meta.data[i]="Other"
sobj@meta.data[which(sobj@meta.data[colName]==i),i]=i
print(DimPlot(sobj, reduction = "umap", group.by = i,order=c(i,'other'),raster=raster) + ggplot2::labs(title=i))
}
dev.off()
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

MM2HSref = read.csv('/local/projects-t3/idea/bherb/annotation/MM2HS_EG100.csv')
colnames(MM2HSref) = c('MM_ID','MM_Gene','HS_ID','HS_Gene')
HS2MMref = read.csv('/local/projects-t3/idea/bherb/annotation/HS2MM_EG100.csv')
colnames(HS2MMref) = c('HS_ID','HS_Gene','MM_ID','MM_Gene')

HStoMM <- function(x){
    tmpind = match(as.character(x),HS2MMref$HS_Gene)
return(as.character(HS2MMref$MM_Gene[tmpind]))
}

MMtoHS <- function(x){
    tmpind = match(as.character(x),MM2HSref$MM_Gene)
return(as.character(MM2HSref$HS_Gene[tmpind]))
}

HSTF = read.csv('/local/projects-t3/idea/bherb/annotation/Hsapiens/Human_TF.csv')
NPlist = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20
NPlist = MMtoHS(as.character(NPlist[,1]))
TF_NP = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/TF_NP_list_from_Moffitt.csv') #direct from Moffitt suppmental materials 

for(i in 1:ncol(TF_NP)){
    TF_NP[,i] = MMtoHS(TF_NP[,i])
} ## human convention 

NPlist2 = unique(c(na.omit(as.character(TF_NP$Neuropeptides)),NPlist))

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

getTips = function(tree=NA,node=NA){
tmp = getDescendants(tree,node)
tmp2 = tmp[which(tmp<=length(tree$tip.label))]
return(tmp2)
}

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

### new integration of mouse and human July 2023

MMHSneu_int = readRDS('./SeuratObj/Hypothalamus_HumanMouseInt.rds')
treeTrim = readRDS('./Analysis/tree80_L15_treeTrim.rds')
treeTrim = treeTrim[intersect(colnames(MMHSneu_int),rownames(treeTrim)),]

## res 5 is K116 - actually 108 groups

MMHSneu_int@meta.data$HSatlas_K6 = NA
MMHSneu_int@meta.data[match(rownames(treeTrim),colnames(MMHSneu_int)),'HSatlas_K6'] = treeTrim[,'K6']

MMHSneu_int@meta.data$HSatlas_K16 = NA
MMHSneu_int@meta.data[match(rownames(treeTrim),colnames(MMHSneu_int)),'HSatlas_K16'] = treeTrim[,'K16']

MMHSneu_int@meta.data$HSatlas_K55 = NA
MMHSneu_int@meta.data[match(rownames(treeTrim),colnames(MMHSneu_int)),'HSatlas_K55'] = treeTrim[,'K55']

MMHSneu_int@meta.data$HSatlas_K116 = NA
MMHSneu_int@meta.data[match(rownames(treeTrim),colnames(MMHSneu_int)),'HSatlas_K116'] = treeTrim[,'K116']

MMHSneu_int@meta.data$Species = 'MouseAdult'

MMHSneu_int@meta.data$Species[grep('Dev',MMHSneu_int@meta.data$Dataset)] = 'MouseDev'
#MMHSneu_int@meta.data$Species[grep('P1',MMHSneu_int@meta.data$sample)] = 'MouseDev'
#MMHSneu_int@meta.data$Species[grep('P4',MMHSneu_int@meta.data$sample)] = 'MouseDev'
#MMHSneu_int@meta.data$Species[grep('P8',MMHSneu_int@meta.data$sample)] = 'MouseDev'
MMHSneu_int@meta.data$Species[grep('Siletti',MMHSneu_int@meta.data$Dataset)] = 'HumanAdult'
MMHSneu_int@meta.data$Species[grep('Zhou',MMHSneu_int@meta.data$Dataset)] = 'HumanDev'
MMHSneu_int@meta.data$Species[grep('Herb',MMHSneu_int@meta.data$Dataset)] = 'HumanDev'

## just checking integration 

PrintColCat = c('Species','Dataset','C66','HSatlas_K55') 

pdf(file=paste("./TestPlots/MMHSneu_July2023_Check_Clustering.pdf",sep=''),width=12,height=8)
for(k in PrintColCat){
print(DimPlot(MMHSneu_int, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
dev.off()

subplots(MMHSneu_int,colName='Species',filename="./TestPlots/MMHSneu_July2023_IndvSpecies.pdf",raster=TRUE)

pdf(file=paste("./TestPlots/MMHSneu_July2023_Check_POMC.pdf",sep=''),width=12,height=8)
print(FeaturePlot(MMHSneu_int, features = 'POMC', pt.size = 0.2))
dev.off()

## metaneighbor  - running into memory issues above C185/K116

mouse_clust = 'C465'
human_clust = 'K444'

MMHSneu_int@meta.data[,paste0('HSatlas_',human_clust)] = NA
cellOv = intersect(rownames(treeTrim),colnames(MMHSneu_int))
#MMHSneu_int@meta.data[match(rownames(treeTrim),colnames(MMHSneu_int)),'HSatlas_K55'] = treeTrim[,'K55']
MMHSneu_int@meta.data[cellOv,paste0('HSatlas_',human_clust)] = treeTrim[cellOv,human_clust]
HSatlasInd = which(!is.na(MMHSneu_int@meta.data[,paste0('HSatlas_',human_clust)]))
HSatlasMat = MMHSneu_int@assays$integrated@scale.data[,HSatlasInd]
HSatlasPheno = data.frame(Sample_ID=colnames(MMHSneu_int)[HSatlasInd],Celltype = paste0(human_clust,'_',MMHSneu_int@meta.data[HSatlasInd,paste0('HSatlas_',human_clust)]),Study_ID = "HSatlas")
#HSatlasPheno$Celltype = paste0('Human_',HSatlasPheno$Celltype)
 #81471
#Celltype = MMHSneu_int@meta.data[HSatlasInd,paste0('HSatlas_',human_clust)]
## get nuclei and neuron type? - predicted? 

hypoMapInd = which(!is.na(MMHSneu_int@meta.data[,paste0(mouse_clust,'_named')])) 
hypoMapMat = MMHSneu_int@assays$integrated@scale.data[,hypoMapInd]
hypoMapPheno = data.frame(Sample_ID=colnames(MMHSneu_int)[hypoMapInd],Celltype = MMHSneu_int@meta.data[hypoMapInd,paste0(mouse_clust,'_named')],Study_ID = "hypoMap") ## 54458 cells 

## for simplicity, set min to 20 cells per type
hypoMapOK = hypoMapPheno$Celltype%in%names(table(hypoMapPheno$Celltype))[which(table(hypoMapPheno$Celltype)>=20)] ## kicks out 43 - and now down to 45 groups
hypoMapPheno2 = hypoMapPheno[hypoMapOK,]
hypoMapMat2 = hypoMapMat[,hypoMapOK]

HSatlasOK = HSatlasPheno$Celltype%in%names(table(HSatlasPheno$Celltype))[which(table(HSatlasPheno$Celltype)>=20)] ## kicks out 0
HSatlasPheno2 = HSatlasPheno[HSatlasOK,]
HSatlasMat2 = HSatlasMat[,HSatlasOK]

rm(HSatlasPheno,hypoMapPheno)
gc()
rm(HSatlasMat,hypoMapMat)
gc()

varGene =  MMHSneu_int@assays$integrated@var.features

gc()
HSatlasMat2 = as.matrix(HSatlasMat2[varGene,])
hypoMapMat2 = as.matrix(hypoMapMat2[varGene,])

#compMat2(xmat=HSatlasMat,ymat=hypoMapMat,xpheno=HSatlasPheno,ypheno=hypoMapPheno,cexRow=0.1,cexCol=0.1)

totPheno = rbind(HSatlasPheno2,hypoMapPheno2)
rownames(totPheno)= totPheno$Sample_ID
totPheno  = totPheno[,-1]
colnames(totPheno) = c('cell_type','study_id')
gc()

totMat = SingleCellExperiment(Matrix(cbind(HSatlasMat2,hypoMapMat2),sparse=TRUE),colData=totPheno)

rm(HSatlasMat2,hypoMapMat2)
gc()
rm(HSatlasPheno2,hypoMapPheno2)
gc()

celltype_NV = MetaNeighborUS(var_genes=varGene,dat=totMat,study_id=totMat$study_id,cell_type=totMat$cell_type,fast_version=TRUE)

cols=rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))

breaks = seq(0,1,length=101)
pdf(paste0('./TestPlots/HSatlas',human_clust,'_VS_hypoMap',mouse_clust,'_MetaNeighbor_20min_July_2023_IntDat.pdf'),width=110,height=110)
print(gplots::heatmap.2(celltype_NV,margins=c(40,40),keysize=1,key.xlab="AUROC",key.title="NULL",trace="none",density.info="none",col=cols,breaks=breaks,offsetRow=0.1,offsetCol=0.1,cexRow=2,cexCol=2))
dev.off()

top_hits_IntDat=topHits(cell_NV = celltype_NV,dat=totMat,study_id=totMat$study_id,cell_type=totMat$cell_type,threshold=0.8)

#topHitsByStudy_IntDat=topHitsByStudy(auroc = celltype_NV,threshold=0.8)

mn_output_IntDat = print(gplots::heatmap.2(celltype_NV,margins=c(25,25),keysize=1,key.xlab="AUROC",key.title="NULL",trace="none",density.info="none",col=cols,breaks=breaks,offsetRow=0.1,offsetCol=0.1,cexRow=0.5,cexCol=0.5))

pdf(paste0('./TestPlots/HSatlas',human_clust,'_VS_hypoMap',mouse_clust,'_MetaNeighbor_Dendogram_test_July_2023_IntDat.pdf'))
plot(mn_output_IntDat[["rowDendrogram"]])
dev.off()

celltype_NV = celltype_NV

top_hits_IntDat = top_hits_IntDat
saveRDS(celltype_NV,file=paste0('./Analysis/MN_July_2023_celltype_NV_human',human_clust,'_mouse',mouse_clust,'.rds'))
saveRDS(top_hits_IntDat,file=paste0('./Analysis/MN_July_2023_top_hits_NV_human',human_clust,'_mouse',mouse_clust,'.rds'))

