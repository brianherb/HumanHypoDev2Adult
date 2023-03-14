## to replicate conda environment used in paper, install conda, and then call the following in terminal: conda env create -f environment_droplet.yml - yml file here: https://github.com/brianherb/HumanHypothalamusDev/blob/main/environment_droplet.yml

library(Seurat)
library(limma)
library(plotly)
library(htmlwidgets)
library("rmarkdown")
library(scales)
library(monocle3)
library(projectR)
library(plotrix) 
library(magrittr)
library(volcano3D)
library(ggrepel)
library(RColorBrewer)
library(igraph)

devtools::source_url('https://github.com/brianherb/HumanHypothalamusDev/blob/main/HypoPub_functions_reference.R?raw=TRUE')

## Paths in these scripts mirror directories used in https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/ We suggest that you establish a single working directory and add the following folders: Analysis, Counts, Networks, Reference_Files, SeuratObj, and TestPlots - a directory can be created in R like so: mkdir('./Counts')


## Kriegstein data 
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_PostSCT.rda')) #HypoNormDat

## GW18 Cortex
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18CtxSamples_PostSCT.rda')) #GW18CtxNormDat

## GW18 GE
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18GESamples_PostSCT.rda')) #GW18GENormDat

GW18AllGE4KCombDat = c(HypoNormDat[['GW18']],GW18CtxNormDat,GW18GENormDat)
names(GW18AllGE4KCombDat)[1] = 'GW18_Hypo'

GW18AllGE4KSamples = names(GW18AllGE4KCombDat)

for(i in GW18AllGE4KSamples){
GW18AllGE4KCombDat[[i]] = subset(GW18AllGE4KCombDat[[i]],cells = colnames(GW18AllGE4KCombDat[[i]])[sample(c(1:ncol(GW18AllGE4KCombDat[[i]])),4000)])
}

for (i in GW18AllGE4KSamples) {
    GW18AllGE4KCombDat[[i]] <- NormalizeData(GW18AllGE4KCombDat[[i]], verbose = FALSE)
    GW18AllGE4KCombDat[[i]] <- FindVariableFeatures(GW18AllGE4KCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

GW18AllGE4K.anchors <- FindIntegrationAnchors(object.list = GW18AllGE4KCombDat, dims = 1:30)
GW18AllGE4K.integrated <- IntegrateData(anchorset = GW18AllGE4K.anchors, dims = 1:30)
DefaultAssay(GW18AllGE4K.integrated) <- "integrated"
GW18AllGE4K.integrated <- ScaleData(GW18AllGE4K.integrated, verbose = FALSE)
GW18AllGE4K.integrated <- RunPCA(GW18AllGE4K.integrated, npcs = 30, verbose = FALSE)
GW18AllGE4K.integrated <- RunUMAP(GW18AllGE4K.integrated, reduction = "pca", dims = 1:30)

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18AllGE4K.integrated.rda'

## GW19
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19_NonHySamples_PostSCT.rda')) # GW19_NormDat

GW19AllGE4KCombDat = c(HypoNormDat[['GW19_hypo']],GW19_NormDat)
names(GW19AllGE4KCombDat)[1] = 'GW19_Hypo'

GW19AllGE4KSamples = names(GW19AllGE4KCombDat)

## random sampling 

for(i in GW19AllGE4KSamples){
if(length(colnames(GW19AllGE4KCombDat[[i]])) > 4000){
GW19AllGE4KCombDat[[i]] = subset(GW19AllGE4KCombDat[[i]],cells = colnames(GW19AllGE4KCombDat[[i]])[sample(c(1:ncol(GW19AllGE4KCombDat[[i]])),4000)])
}
}

for (i in GW19AllGE4KSamples) {
    GW19AllGE4KCombDat[[i]] <- NormalizeData(GW19AllGE4KCombDat[[i]], verbose = FALSE)
    GW19AllGE4KCombDat[[i]] <- FindVariableFeatures(GW19AllGE4KCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

GW19AllGE4K.anchors <- FindIntegrationAnchors(object.list = GW19AllGE4KCombDat, dims = 1:30)
GW19AllGE4K.integrated <- IntegrateData(anchorset = GW19AllGE4K.anchors, dims = 1:30)
DefaultAssay(GW19AllGE4K.integrated) <- "integrated"
GW19AllGE4K.integrated <- ScaleData(GW19AllGE4K.integrated, verbose = FALSE)
GW19AllGE4K.integrated <- RunPCA(GW19AllGE4K.integrated, npcs = 30, verbose = FALSE)
GW19AllGE4K.integrated <- RunUMAP(GW19AllGE4K.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/GW19AllGE4K_Merged.pdf',width=12,height=8)
print(DimPlot(GW19AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("sample")){
print(DimPlot(GW19AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19AllGE4K_2D_integrated.rds'

## GW20_34 - 

load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34_NonHySamples_PostSCT.rda')) # GW20_34_NormDat

GW20_34AllGE4KCombDat = c(HypoNormDat[['GW20_34_hypo']],GW20_34_NormDat)
names(GW20_34AllGE4KCombDat)[1] = 'GW20_34_Hypo'

GW20_34AllGE4KSamples = names(GW20_34AllGE4KCombDat)

## random sampling 

for(i in GW20_34AllGE4KSamples){
if(length(colnames(GW20_34AllGE4KCombDat[[i]])) > 4000){
GW20_34AllGE4KCombDat[[i]] = subset(GW20_34AllGE4KCombDat[[i]],cells = colnames(GW20_34AllGE4KCombDat[[i]])[sample(c(1:ncol(GW20_34AllGE4KCombDat[[i]])),4000)])
}
}

for (i in GW20_34AllGE4KSamples) {
    GW20_34AllGE4KCombDat[[i]] <- NormalizeData(GW20_34AllGE4KCombDat[[i]], verbose = FALSE)
    GW20_34AllGE4KCombDat[[i]] <- FindVariableFeatures(GW20_34AllGE4KCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

GW20_34AllGE4K.anchors <- FindIntegrationAnchors(object.list = GW20_34AllGE4KCombDat, dims = 1:30)
GW20_34AllGE4K.integrated <- IntegrateData(anchorset = GW20_34AllGE4K.anchors, dims = 1:30)
DefaultAssay(GW20_34AllGE4K.integrated) <- "integrated"
GW20_34AllGE4K.integrated <- ScaleData(GW20_34AllGE4K.integrated, verbose = FALSE)
GW20_34AllGE4K.integrated <- RunPCA(GW20_34AllGE4K.integrated, npcs = 30, verbose = FALSE)
GW20_34AllGE4K.integrated <- RunUMAP(GW20_34AllGE4K.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/GW20_34AllGE4K_Merged.pdf',width=12,height=8)
print(DimPlot(GW20_34AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("sample")){
print(DimPlot(GW20_34AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34AllGE4K_2D_integrated.rds'



## Integrate GW18, GW19 and GW20

options(future.globals.maxSize = 4000 * 1024^2)

AllCells = c(colnames(GW18AllGE4K.integrated),colnames(GW19AllGE4K.integrated),colnames(GW20_34AllGE4K.integrated))

AllGE4KCombDat = c(GW18AllGE4KCombDat,GW19AllGE4KCombDat,GW20_34AllGE4KCombDat)

AllGE4KSamples = names(AllGE4KCombDat)


for(i in AllGE4KSamples){
AllGE4KCombDat[[i]] = subset(AllGE4KCombDat[[i]],cells = intersect(AllCells,colnames(AllGE4KCombDat[[i]])))
}

for (i in AllGE4KSamples) {

    DefaultAssay(AllGE4KCombDat[[i]]) = 'RNA'
    AllGE4KCombDat[[i]] <- subset(AllGE4KCombDat[[i]], percent.mt<=10 & nFeature_RNA>=500)
    AllGE4KCombDat[[i]] <- NormalizeData(AllGE4KCombDat[[i]], verbose = FALSE)
    AllGE4KCombDat[[i]] <- FindVariableFeatures(AllGE4KCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}


features <- SelectIntegrationFeatures(object.list = AllGE4KCombDat)

for (i in 1:length(AllGE4KCombDat)) {
 AllGE4KCombDat[[i]]<- ScaleData(AllGE4KCombDat[[i]], features = features, verbose = FALSE)
    AllGE4KCombDat[[i]] <- RunPCA(AllGE4KCombDat[[i]], features = features, verbose = FALSE)
}


AllGE4K.anchors <- FindIntegrationAnchors(object.list = AllGE4KCombDat,  anchor.features = features, reduction = "rpca")
AllGE4K.integrated <- IntegrateData(anchorset = AllGE4K.anchors, dims = 1:30)
DefaultAssay(AllGE4K.integrated) <- "integrated"
AllGE4K.integrated <- ScaleData(AllGE4K.integrated, verbose = FALSE)
AllGE4K.integrated <- RunPCA(AllGE4K.integrated, npcs = 30, verbose = FALSE)
AllGE4K.integrated <- RunUMAP(AllGE4K.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/AllGE4K_Merged_mt10n500.pdf',width=12,height=8)
print(DimPlot(AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("sample")){
print(DimPlot(AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18_GW20_34AllGE4K_2D_integrated_mt10n500.rds'

DeConCG = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/DeCon_ALL.cells_CoGAPS_24samples_24552genes_1sets_5patterns_20000iterations_5fP.rds'))

DeConCGpat = DeConCG@featureLoadings

huGene=MMtoHS(rownames(DeConCGpat))

DeConCGpat = DeConCGpat[-which(is.na(huGene)),]

rownames(DeConCGpat)=MMtoHS(rownames(DeConCGpat))

AllGE4K.integrated.data = as.matrix(AllGE4K.integrated@assays$SCT[1:21810,1:95107])

tmpMat = DeConCGpat
uniqGenes = names(table(rownames(tmpMat)))[which(table(rownames(tmpMat))==1)]
tmpMat = tmpMat[match(uniqGenes,rownames(tmpMat)),]

proRes=projectR(data=AllGE4K.integrated.data, loadings = tmpMat, full = FALSE) ##
AllGE4K.integrated@meta.data = cbind(AllGE4K.integrated@meta.data,t(proRes))

clrs00=c("black","darkred","red","orange","yellow")# "originColor" - 1st to last correspond to lo to hi values

pdf('./TestPlots/All_CtxHyGE_mt10n500_projection_DeCon_5pattern.pdf',width=12,height=8)
for(j in c("sample","Cell_Type","Maj_Clust_From_Ref")){
print(DimPlot(AllGE4K.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE,raster=TRUE,shuffle=TRUE))
}
print(FeaturePlot(AllGE4K.integrated,features='nCount_RNA', pt.size = 0.4,raster=TRUE))
for(k in rownames(proRes)){
clr=color.scale(proRes[k,],extremes=clrs00)
print(FeaturePlot(AllGE4K.integrated,features=k, pt.size = 0.4,raster=TRUE))
}
dev.off()

## plot saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/TestPlots/All_CtxHyGE_mt10n500_projection_DeCon_5pattern.pdf'

## proRes saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/All_CtxHyGE_mt10n500projectionResults_DeCon_5pattern.rda'


### plot 3D UMAP

AllGE4K_3D.integrated <- RunUMAP(AllGE4K.integrated, reduction = "pca", dims = 1:30,n.components = 3L)

tmpcells = gsub('GW18_','GW18_hypo_',colnames(AllGE4K_3D.integrated))

tmpcells[grep('GW18',tmpcells)] = paste(tmpcells[grep('GW18',tmpcells)] ,'-1',sep='')

AllGE4K_3D.integrated@meta.data$HypoCS13_GW25_SC = HypoCS13_GW25_mt10@meta.data[tmpcells,'seurat_clusters']

AllGE4K_3D.integrated@meta.data$Stage = 'GW18'

AllGE4K_3D.integrated@meta.data$Stage[grep('GW19',AllGE4K_3D.integrated@meta.data$sample)] = 'GW19'

AllGE4K_3D.integrated@meta.data$Stage[grep('GW20_34',AllGE4K_3D.integrated@meta.data$sample)] = 'GW20_34'

AllGE4K_3D.integrated@meta.data$BrainRegion = NA

AllGE4K_3D.integrated@meta.data$BrainRegion[grep('GW18_',AllGE4K_3D.integrated@meta.data$sample)] =  'Cortex'

AllGE4K_3D.integrated@meta.data$BrainRegion[grep('GW19_',AllGE4K_3D.integrated@meta.data$sample)] =  'Cortex'

AllGE4K_3D.integrated@meta.data$BrainRegion[grep('GW20_34_',AllGE4K_3D.integrated@meta.data$sample)] =  'Cortex'

AllGE4K_3D.integrated@meta.data$BrainRegion[which(AllGE4K_3D.integrated@meta.data$sample == 'GW18')] = 'Hypo'

AllGE4K_3D.integrated@meta.data$BrainRegion[which(AllGE4K_3D.integrated@meta.data$sample == 'GW19_hypo')] = 'Hypo'

AllGE4K_3D.integrated@meta.data$BrainRegion[which(AllGE4K_3D.integrated@meta.data$sample == 'GW20_34_hypo')] = 'Hypo'

AllGE4K_3D.integrated@meta.data$BrainRegion[grep('GE',AllGE4K_3D.integrated@meta.data$sample)] =  'GE'


### call seurat clusters in 3D

AllGE4K_3D.integrated@meta.data$old_seurat_clusters = AllGE4K_3D.integrated@meta.data$seurat_clusters

AllGE4K_3D.integrated <- FindNeighbors(AllGE4K_3D.integrated, dims = 1:30, verbose = FALSE)
AllGE4K_3D.integrated <- FindClusters(AllGE4K_3D.integrated, verbose = FALSE, resolution=0.5)

### Calculate differentially expressed genes between seurat clusters (to find cluster marker genes, regardless of brain region) 
## Sup table 19 

DefaultAssay(AllGE4K_3D.integrated)='RNA'

AllGE4K_3D.integrated@active.ident=as.factor(AllGE4K_3D.integrated$seurat_clusters)

AllGE4K_3D_deg = FindAllMarkers(AllGE4K_3D.integrated,features=rownames(AllGE4K_3D.integrated))

save(AllGE4K_3D_deg,file='./Analysis/AllGE4K_3D_deg_mt10n500.rda',compress=TRUE)

write.csv(AllGE4K_3D_deg,file='./Analysis/AllGE4K_deg_bw_cluster_mt10n500.csv')

write.csv(AllGE4K_3D_deg[which(AllGE4K_3D_deg$p_val_adj<=0.05),],file='./Analysis/AllGE4K_sig_deg_bw_cluster_mt10n500.csv')

ClustDifGenes = unique(AllGE4K_3D_deg$gene[which(AllGE4K_3D_deg$p_val_adj<=0.05)]) #7590


### Calculate differentially expressed genes between brain regions within each seurat cluster
## Sup table 20

RegionDifBySC = vector(mode='list',length=length(unique(AllGE4K_3D.integrated$seurat_clusters)))

names(RegionDifBySC) = unique(AllGE4K_3D.integrated$seurat_clusters)

DefaultAssay(AllGE4K_3D.integrated)='RNA'

AllGE4K_3D.integrated@active.ident=as.factor(AllGE4K_3D.integrated$BrainRegion)

for(i in names(RegionDifBySC)){

tmpSC = subset(AllGE4K_3D.integrated,cells=colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$seurat_clusters==i)])

RegionDifBySC[[i]] = FindAllMarkers(tmpSC,logfc.threshold=0,min.pct = 0,min.cells.feature=0, min.cells.group=0, return.thresh = 1)

cat(i)
}

save(RegionDifBySC,file='./Analysis/All_CtxHyGE_deg_bySC_mt10n500.rda',compress=TRUE)


for(i in names(RegionDifBySC)){
if(nrow(RegionDifBySC[[i]])>0){
RegionDifBySC[[i]]$seurat_cluster = i   
}    
}

RegionDifBySCmat = do.call(rbind,RegionDifBySC)

RegionDifBySCmatSig = RegionDifBySCmat[which(RegionDifBySCmat$p_val_adj<=0.05),]

write.csv(RegionDifBySCmatSig,file='./Analysis/AllGE4K_deg_BrainRegion_within_cluster_mt10n500.csv')

for(i in 1:length(RegionDifBySC)){
if(i==1){
RegDifGenes = unique(RegionDifBySC[[i]]$gene[which(RegionDifBySC[[i]]$p_val_adj<=0.05)])
} else {
RegDifGenes = unique(c(RegDifGenes,RegionDifBySC[[i]]$gene[which(RegionDifBySC[[i]]$p_val_adj<=0.05)]))
}
}

DifGenes = unique(c(RegDifGenes,ClustDifGenes))
## 8359 - these were used for GENIE3 network and Gene module creation

save(DifGenes,file='./Analysis/All_CtxHyGE_Clustering_mt10n500_DifGenes.rda')

write.csv(data.frame(Name=DifGenes),file='./Analysis/All_CtxHyGE_Clustering_mt10n500_DifGenes.csv')

write.csv(data.frame(Name=intersect(TFs,DifGenes)),file='./Analysis/All_CtxHyGE_Clustering_mt10n500_DifTFs.csv',row.names=FALSE,quote=FALSE)

## DifGenes object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/All_CtxHyGE_Clustering_mt10n500_DifGenes.rda'

### Use polar coordinate code to help find commonalities in gene expression across brain regions. 

ComDifDEG = vector(mode='list',length=length(unique(AllGE4K_3D_deg$cluster)))

names(ComDifDEG) = c(0:max(as.numeric(as.character(AllGE4K_3D_deg$cluster))))

for(j in names(ComDifDEG)){

TMPind = which(AllGE4K_3D.integrated@meta.data$seurat_clusters==j)

TMPmarkers = AllGE4K_3D_deg[which(AllGE4K_3D_deg$cluster==j & AllGE4K_3D_deg$avg_log2FC>0 & AllGE4K_3D_deg$p_val_adj<=0.05),]

rownames(TMPmarkers) = TMPmarkers$gene

TMPmarkers$CtxPerEx = 0
TMPmarkers$GePerEx = 0
TMPmarkers$HyPerEx = 0

for(i in 1:nrow(TMPmarkers)){

dummy = c(0,0,0)
names(dummy) = c('Cortex','GE','Hypo')

dummy2 = c(0,0,0)
names(dummy2) = c('Cortex','GE','Hypo')

sampCounts=table(AllGE4K_3D.integrated$BrainRegion[intersect(TMPind,which(AllGE4K_3D.integrated@assays$RNA@counts[TMPmarkers$gene[i],]>0))])

dummy[names(sampCounts)] = sampCounts

totCells = table(AllGE4K_3D.integrated$BrainRegion[TMPind])

dummy2[names(totCells)] = totCells

tmp = round((dummy/dummy2)*100,2)

TMPmarkers$CtxPerEx[i] = tmp['Cortex']
TMPmarkers$GePerEx[i] = tmp['GE']
TMPmarkers$HyPerEx[i] = tmp['Hypo']

cat(i)
}
## values from polar plot
tmpSC = subset(AllGE4K_3D.integrated,cells=colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$seurat_clusters==j)])

if(length(unique(tmpSC$BrainRegion))<3) {
ComDifDEG[[j]] = TMPmarkers
next
}

if(1%in%table(AllGE4K_3D.integrated$BrainRegion[TMPind])) {
ComDifDEG[[j]] = TMPmarkers
next
}

out5 = tmpSC$BrainRegion

dat5 = tmpSC@assays$RNA@counts[intersect(rownames(tmpSC),TMPmarkers$gene),]

test=polar_coords(outcome=out5,data=t(dat5))

tmpLab = test@df[[1]]

TMPmarkers$RadialCall = tmpLab[rownames(TMPmarkers),'lab']

TMPmarkers$RadialSig = tmpLab[rownames(TMPmarkers),'pvalue']

ComDifDEG[[j]] = TMPmarkers
}

save(ComDifDEG,file='./Analysis/ComDifDEG.rda')

## ComDifDEG object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/ComDifDEG.rda'


### Investigate bulk cell type expression to find brain region specific gene expression per cell type. 

## Assign seurat clusters to cell type

AllGE4K_3D.integrated$PseudoBulkGroups = 'Other'

#21
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(21)),'PseudoBulkGroups'] = 'Tanycyte_Hypo'

#18
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(18)),'PseudoBulkGroups'] = 'Ependymal_Hypo' ## there are 2 cortex cells in here - other

## 13
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(13)),'PseudoBulkGroups'] = 'HOPXpos-Progenitor_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(13)),'PseudoBulkGroups'] = 'HOPXpos-Progenitor_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(13)),'PseudoBulkGroups'] = 'HOPXpos-Progenitor_Hypo'

## 16
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(16)),'PseudoBulkGroups'] = 'Oligodendrocytes_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(16)),'PseudoBulkGroups'] = 'Oligodendrocytes_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(16)),'PseudoBulkGroups'] = 'Oligodendrocytes_Hypo'

## 17
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(17)),'PseudoBulkGroups'] = 'Astrocytes_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(17)),'PseudoBulkGroups'] = 'Astrocytes_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(17)),'PseudoBulkGroups'] = 'Astrocytes_Hypo'

#19
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$seurat_clusters%in%c(19)),'PseudoBulkGroups'] = 'Microglia'

#20
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$seurat_clusters%in%c(20)),'PseudoBulkGroups'] = 'Epithelial'


#23
AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$seurat_clusters%in%c(23)),'PseudoBulkGroups'] = 'Blood'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(10,9,6)),'PseudoBulkGroups'] = 'DivProgenitor_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(10,9,6)),'PseudoBulkGroups'] = 'DivProgenitor_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(10,9,6)),'PseudoBulkGroups'] = 'DivProgenitor_Hypo'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(4)),'PseudoBulkGroups'] = 'RadialGlia_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(4)),'PseudoBulkGroups'] = 'RadialGlia_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(4)),'PseudoBulkGroups'] = 'RadialGlia_Hypo'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(11)),'PseudoBulkGroups'] = 'EarlyNeuPro_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(11)),'PseudoBulkGroups'] = 'EarlyNeuPro_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(11)),'PseudoBulkGroups'] = 'EarlyNeuPro_Hypo'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(12)),'PseudoBulkGroups'] = 'LateExcNeuPro_Ctx'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(1)),'PseudoBulkGroups'] = 'LateInhibNeuPro_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(1)),'PseudoBulkGroups'] = 'LateInhibNeuPro_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(1)),'PseudoBulkGroups'] = 'LateInhibNeuPro_Hypo'


AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(5)),'PseudoBulkGroups'] = 'LHX6pos-InhibNeu_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(5)),'PseudoBulkGroups'] = 'LHX6pos-InhibNeu_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(5)),'PseudoBulkGroups'] = 'LHX6pos-InhibNeu_Hypo'


AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(7)),'PseudoBulkGroups'] = 'FOXP1pos-InhibNeu_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(7)),'PseudoBulkGroups'] = 'FOXP1pos-InhibNeu_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(7)),'PseudoBulkGroups'] = 'FOXP1pos-InhibNeu_Hypo'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(8)),'PseudoBulkGroups'] = 'CALB2pos-InhibNeu_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(8)),'PseudoBulkGroups'] = 'CALB2pos-InhibNeu_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(8)),'PseudoBulkGroups'] = 'CALB2pos-InhibNeu_Hypo'


AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(3,0,2)),'PseudoBulkGroups'] = 'ExcItNeu_Cortex'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='GE' & AllGE4K_3D.integrated$seurat_clusters%in%c(3,0,2)),'PseudoBulkGroups'] = 'ExcItNeu_GE'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Hypo' & AllGE4K_3D.integrated$seurat_clusters%in%c(3,0,2)),'PseudoBulkGroups'] = 'ExcItNeu_Hypo'

AllGE4K_3D.integrated@meta.data[which(AllGE4K_3D.integrated$BrainRegion=='Cortex' & AllGE4K_3D.integrated$seurat_clusters%in%c(14,15)),'PseudoBulkGroups'] = 'ExcEtNeu_Cortex'


tmpGrps = gsub('_Hypo','',AllGE4K_3D.integrated@meta.data$PseudoBulkGroups)

tmpGrps = gsub('_Cortex','',tmpGrps)

tmpGrps = gsub('_GE','',tmpGrps)

AllGE4K_3D.integrated@meta.data$PseudoBulkGroupsComb = tmpGrps

## write out metadata for Sup table 18

write.csv(AllGE4K_3D.integrated@meta.data[,c(4,16,15,2,3,5,10,21,24,25)],file='./Analysis/All_CtxHyGE_mt10n500_CellLevel_metadata.csv',quote=FALSE)

## Metadata object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/All_CtxHyGE_mt10n500_CellLevel_metadata.csv'

###################
### Figure 4a #####
###################

plot.data <- FetchData(object = AllGE4K_3D.integrated, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "PseudoBulkGroups","PseudoBulkGroupsComb","BrainRegion"))

p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~PseudoBulkGroupsComb,colors = hue_pal()(length(unique(plot.data$PseudoBulkGroupsComb))),type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, "./TestPlots/All_CtxHyGE_Clustering_mt10n500_3D_PseudoBulkGroupsComb.html")


p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~PseudoBulkGroups,colors = hue_pal()(length(unique(plot.data$PseudoBulkGroups))),type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, "./TestPlots/All_CtxHyGE_Clustering_mt10n500_3D_PseudoBulkGroups.html")


###################
### Figure 4b #####
###################

p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~BrainRegion,colors = hue_pal()(length(unique(plot.data$BrainRegion))),type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, "./TestPlots/All_CtxHyGE_Clustering_mt10n500_3D_BrainRegion.html")

## also plot in 2D - 

AllGE4K.integrated$PseudoBulkGroups = AllGE4K_3D.integrated@meta.data[colnames(AllGE4K.integrated),'PseudoBulkGroups']


AllGE4K.integrated$PseudoBulkGroupsComb = AllGE4K_3D.integrated@meta.data[colnames(AllGE4K.integrated),'PseudoBulkGroupsComb']

AllGE4K.integrated$BrainRegion = AllGE4K_3D.integrated@meta.data[colnames(AllGE4K.integrated),'BrainRegion']

AllGE4K.integrated$Phase = AllGE4K_3D.integrated@meta.data[colnames(AllGE4K.integrated),'Phase']

AllGE4K.integrated$sc3D = AllGE4K_3D.integrated@meta.data[colnames(AllGE4K.integrated),'seurat_clusters']


pdf(paste('./TestPlots/All_CtxHyGE_mt10n500_BrainRegions_Groups_2D.pdf',sep=''),width=12,height=8)
for(j in c("BrainRegion","Phase","sc3D","PseudoBulkGroups","PseudoBulkGroupsComb")){
print(DimPlot(AllGE4K.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE,raster=TRUE,shuffle=TRUE))
}
dev.off()


pdf(paste('./TestPlots/All_CtxHyGE_mt10n500_BrainRegions_NoLabel_2D.pdf',sep=''),width=12,height=8)
for(j in c("BrainRegion")){
print(DimPlot(AllGE4K.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE,raster=TRUE,shuffle=TRUE))
}
dev.off()


###################
### Figure 4c #####
###################

AllGE4K_3D.integrated@active.ident=as.factor(AllGE4K_3D.integrated$PseudoBulkGroups)

Bulk = AverageExpression(AllGE4K_3D.integrated)

### BrainRegions per group

data2= table(AllGE4K_3D.integrated@meta.data$BrainRegion,AllGE4K_3D.integrated$PseudoBulkGroupsComb)

data3 = data.table::melt(data2)

# Stacked

pdf(file=paste("./TestPlots/All_CtxHyGE_mt10n500_PseudoBulkGroupsComb_Composition.pdf",sep=''),width=16,height=8)

ggplot(data3, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="fill", stat="identity") + theme(axis.text.x=element_text(angle=45, hjust=1,size=12))
dev.off()


###################
### Figure 4d #####
###################


DefaultAssay(AllGE4K_3D.integrated)='RNA'

AllGE4K_3D.integrated@active.ident=as.factor(AllGE4K_3D.integrated$BrainRegion)


TFs = intersect(rownames(AllGE4K_3D.integrated),unique(c(HSTF$Name,na.omit(TF_NP$Transcription_factors))))

for(i in c(0:23)){

tmpSC = subset(AllGE4K_3D.integrated,cells=colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$seurat_clusters==i)])

if(length(unique(tmpSC$BrainRegion))<3) next

out5 = tmpSC$BrainRegion

dat5 = tmpSC@assays$RNA@counts[intersect(rownames(tmpSC),TFs),]

test=polar_coords(outcome=out5,data=t(dat5))

tmpLab = test@df[[1]]
tmpLab = tmpLab[which(tmpLab$r>0.5),]
if(nrow(tmpLab)==0)  next
tmpLab$lab=rownames(tmpLab)

pdf(paste('./TestPlots/CxtGeHy_mt10n500_Test_polarCoords_Clust',i,'.pdf',sep=''))

print(radial_ggplot(test,colours=c("grey60", "green3", "cyan","blue", "purple", "red", "gold2")) + geom_text_repel(
    data = tmpLab,label=tmpLab$lab))
dev.off()

cat(paste(i,', ',sep=''))
}

###################
### Figure 4e #####
###################

AllGE4K_3D.integrated.hm = AllGE4K_3D.integrated

AllGE4K_3D.integrated.hm@active.ident=factor(AllGE4K_3D.integrated.hm$PseudoBulkGroupsComb,levels=unique(AllGE4K_3D.integrated.hm$PseudoBulkGroupsComb)[c(6,11,4,15,19,8,18,1,7,12,3,14,2,10,9,13,5,17,16)])

DefaultAssay(AllGE4K_3D.integrated.hm) <- 'RNA'
AllGE4K_3D.integrated.hm <- NormalizeData(AllGE4K_3D.integrated.hm)

grps = unique(AllGE4K_3D.integrated.hm$PseudoBulkGroups)

plotGroups = grps[c(6,36,22,4,31,19,15,37,20)]

AllGE4K_3D.integrated.hm2 = subset(AllGE4K_3D.integrated,cells=colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated@meta.data$PseudoBulkGroups%in%plotGroups)])

AllGE4K_3D.integrated.hm2@active.ident=factor(AllGE4K_3D.integrated.hm2$PseudoBulkGroupsComb)

test=FindMarkers(AllGE4K_3D.integrated.hm2,ident.1 ='DivProgenitor',ident.2='EarlyNeuPro')

for( i in c('4','6','9','10','11','1','12')){

gn = ComDifDEG[[i]]$gene[which(ComDifDEG[[i]]$gene%in%c(TFs))][1:50]

pdf(paste('./TestPlots/All_CtxHyGE_mt10n500_earlyNeuron_TopMarkers_Cluster',i,'_BulkHeatmap.pdf',sep=''),width=16,height=8)
print(plotBulkHeatmap(AllGE4K_3D.integrated.hm2,gn))
dev.off()
cat(i)
}
}

### check non-dividing precursors 

DefaultAssay(AllGE4K_3D.integrated)='RNA'

test2 = FindMarkers(AllGE4K_3D.integrated,ident.1 = colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$seurat_clusters=='11' & AllGE4K_3D.integrated$Phase=='G1')],ident.2 = colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$seurat_clusters=='4' & AllGE4K_3D.integrated$Phase=='G1')])

AllGE4K_3D.integrated.g1 = subset(AllGE4K_3D.integrated,cells=colnames(AllGE4K_3D.integrated)[which(AllGE4K_3D.integrated$Phase=='G1' & AllGE4K_3D.integrated$PseudoBulkGroupsComb%in%c('RadialGlia','EarlyNeuPro','LateExcNeuPro_Ctx','LateInhibNeuPro','ExcEtNeu','ExcItNeu','CALB2pos-InhibNeu','FOXP1pos-InhibNeu','LHX6pos-InhibNeu'))])

DefaultAssay(AllGE4K_3D.integrated.g1)='RNA'

Idents(AllGE4K_3D.integrated.g1) = 'PseudoBulkGroupsComb'

testTF = FindAllMarkers(AllGE4K_3D.integrated.g1,features=intersect(TFs,rownames(AllGE4K_3D.integrated.g1)))

plotGr = unique(AllGE4K_3D.integrated.g1$PseudoBulkGroups)[c(3,17,9,7,20,10,13,5,8,12,1,18,15,4,19,14,6,21,11,2,22,16,23)]

AllGE4K_3D.integrated.g1@active.ident=factor(AllGE4K_3D.integrated.g1$PseudoBulkGroups,levels=plotGr)

AllGE4K_3D.integrated.g1 <- NormalizeData(AllGE4K_3D.integrated.g1)

testTFsig = testTF[which(testTF$p_val_adj<=0.001 & testTF$avg_log2FC>1 & testTF$pct.2<0.15),]

for(i in c('RadialGlia','EarlyNeuPro','LateExcNeuPro_Ctx','ExcItNeu','ExcEtNeu','LateInhibNeuPro','CALB2pos-InhibNeu','LHX6pos-InhibNeu','FOXP1pos-InhibNeu')){

if(i=='RadialGlia'){
gn = na.omit(testTFsig$gene[which(testTFsig$cluster==i)][1:5])
} else {
gn = c(gn,na.omit(testTFsig$gene[which(testTFsig$cluster==i)][1:5]))
}
}

gn=unique(gn)

## add SOX2, PAX6, ASCL1

gn2 = c(gn[1:5],'SOX2','PAX6','ASCL1',gn[7:29])

pdf('./TestPlots/All_CtxHyGE_mt10n500__NeuronLin_TopMarkers_BulkHeatmap_minZero.pdf',width=8,height=16)
plotBulkHeatmap(AllGE4K_3D.integrated.g1,gn2)
dev.off()

## final heatmap for figure 4E
pdf('./TestPlots/All_CtxHyGE_mt10n500__NeuronLin_TopMarkers_BulkHeatmap_minZero_noEGR1.pdf',width=8,height=16)
plotBulkHeatmap(AllGE4K_3D.integrated.g1,gn2[-2])
dev.off()


#############
## Networks - construct GENIE3 networks and calculate gene expression modules

#############
### Genie3 

## Prepare data matricies 

samples = unique(AllGE4K_3D.integrated$sample)

for(i in samples){

cellind = which(AllGE4K_3D.integrated@meta.data$sample==i)

dat8359 = as.matrix(AllGE4K_3D.integrated@assays$RNA@counts[DifGenes,cellind])

write.csv(t(dat8359),paste('./Networks/',i,'_DifGene8359_Counts_ExpressionMatrix.csv',sep=''))

dist8359 = as.matrix(dist(AllGE4K_3D.integrated@reductions[['pca']]@cell.embeddings[cellind,1:30],diag=TRUE))

dat8359knn5 = smoother_aggregate_nearest_nb(mat = dat8359, D = dist8359, k = 5)

write.csv(t(dat8359knn5),paste('./Networks/',i,'_DifGene8359_KNN5_ExpressionMatrix.csv',sep=''))

cat(paste('\n\n',i,'\n\n',sep=''))

}

## Expression matricies saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Networks/*_DifGene8359_KNN5_ExpressionMatrix.csv'


### example Genie3 script - submitted separately - 

library(GENIE3)
set.seed(123) # For reproducibility of results

Sample = 'GW18_MGE_DifGene8359_KNN5'

tmpdir = paste('/local/scratch/bherb/',Sample,'_TMP',sep='')
if(!dir.exists(tmpdir)){
dir.create(tmpdir)
}
setwd(tmpdir)

hstf =read.csv('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/Analysis/All_CtxHyGE_Clustering_mt10n500_DifTFs.csv',header=TRUE)

exprMatr=read.csv(paste('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/Networks/',Sample,'_ExpressionMatrix.csv',sep=''),row.names=1)

exprMatr = t(as.matrix(exprMatr))

regulators = unique(intersect(hstf$Name,rownames(exprMatr)))

weightMat <- GENIE3(exprMatr,regulators=regulators, nCores=12, verbose=TRUE)

saveRDS(weightMat,file=paste('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/Networks/',Sample,'_WeightMatrix.rds',sep=''))

linkList <- getLinkList(weightMat)

saveRDS(linkList,file=paste('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/Networks/',Sample,'_linkList.rds',sep=''))


## Link lists saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Networks/*_DifGene8359_KNN5_linkList.rds'


###########
## read in genie result 

WtMats=list.files(path='./Networks/',pattern = 'DifGene8359_KNN5_WeightMatrix.rds',full.names=TRUE)

LkLists=list.files(path='./Networks/',pattern = 'DifGene8359_KNN5_linkList.rds',full.names=TRUE)

Links = vector(mode='list',length=length(LkLists))

names(Links) = gsub('_DifGene8359_KNN5_linkList.rds','',list.files(path='./Networks/',pattern = 'DifGene8359_KNN5_linkList.rds'))

for(i in 1:length(LkLists)){
Links[[i]] = readRDS(LkLists[i])
}

## evaluate similarities 

for(i in 1:length(Links)){

Links[[i]]$ID = paste(Links[[i]]$regulatoryGene,Links[[i]]$targetGene,sep='_')

}


## GE 

GEsamples = c("GW18_CGE","GW19_CGE","GW20_34_CGE","GW18_LGE","GW19_LGE","GW20_34_LGE","GW18_MGE","GW19_MGE","GW20_34_MGE")

GEind = matrix(0,nrow=nrow(Links[[GEsamples[1]]]),ncol=length(GEsamples))

colnames(GEind) = GEsamples

for(i in GEsamples){
GEind[,i] = match(Links[[GEsamples[1]]]$ID,Links[[i]]$ID)
}

GEtmp = Links[GEsamples]

for(j in GEsamples){
GEtmp[[j]] = GEtmp[[j]][GEind[,j],]
}

GEtmpWeight = GEtmp[[1]]$weight

for(k in 2:length(GEtmp)){
GEtmpWeight = cbind(GEtmpWeight,GEtmp[[k]]$weight)
}

GEnet = Links[[GEsamples[1]]]

GEnet$weight = rowMeans(GEtmpWeight)

GEnet = GEnet[order(GEnet$weight,decreasing=TRUE),]


## Cortex 

CTXsamples = c("GW18_motor","GW18_parietal","GW18_PFC","GW18_somatosensory","GW18_V1","GW19_motor","GW19_parietal","GW19_PFC","GW19_somato","GW20_34_somato","GW20_34_PFC","GW20_34_V1")

CTXind = matrix(0,nrow=nrow(Links[[CTXsamples[1]]]),ncol=length(CTXsamples))

colnames(CTXind) = CTXsamples

for(i in CTXsamples){
CTXind[,i] = match(Links[[CTXsamples[1]]]$ID,Links[[i]]$ID)
}

CTXtmp = Links[CTXsamples]

for(j in CTXsamples){
CTXtmp[[j]] = CTXtmp[[j]][CTXind[,j],]
}

CTXtmpWeight = CTXtmp[[1]]$weight

for(k in 2:length(CTXtmp)){
CTXtmpWeight = cbind(CTXtmpWeight,CTXtmp[[k]]$weight)
}

CTXnet = Links[[CTXsamples[1]]]

CTXnet$weight = rowMeans(CTXtmpWeight)

CTXnet = CTXnet[order(CTXnet$weight,decreasing=TRUE),]

## hypo

ind = match(Links[["GW18"]]$ID,Links[["GW19_hypo"]]$ID)

ind2 = match(Links[["GW18"]]$ID,Links[["GW20_34_hypo"]]$ID)

cor(Links[["GW18"]]$weight[1:1000],Links[["GW19_hypo"]]$weight[ind[1:1000]])

cor(Links[["GW18"]]$weight[1:1000],Links[["GW20_34_hypo"]]$weight[ind2[1:1000]])

hyponet = Links[["GW18"]]

hyponet$weight = rowMeans(cbind(Links[["GW18"]]$weight,Links[["GW19_hypo"]]$weight[ind],Links[["GW20_34_hypo"]]$weight[ind2]))

hyponet = hyponet[order(hyponet$weight,decreasing=TRUE),]


save(GEnet,CTXnet,hyponet,file='./Networks/All_CtxHyGE_BrainRegion_Networks.rda',compress=TRUE)


## Sup table 21
write.csv(hyponet[1:500000,1:3], file = "./Analysis/Hypothalamus_Genie3_network.csv",quote=FALSE, row.names=FALSE)

## Sup table 22
write.csv(GEnet[1:500000,1:3], file = "./Analysis/GE_Genie3_network.csv",quote=FALSE, row.names=FALSE)

## Sup table 23
write.csv(CTXnet[1:500000,1:3], file = "./Analysis/Cortex_Genie3_network.csv",quote=FALSE, row.names=FALSE)


## plotting HMGB2 and a few target genes to test overlap
pdf(paste('./TestPlots/All_CtxHyGE_mt10n500_HMGB2_targets.pdf',sep=''),width=12,height=8)
for(j in c("sample","Cell_Type","Maj_Clust_From_Ref")){
print(DimPlot(AllGE4K.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE,raster=TRUE,shuffle=TRUE))
}
for(k in c('HMGB2','NUSAP1','NUF2','TUBA1B','TOP2A')){
print(FeaturePlot(AllGE4K.integrated,features=k, pt.size = 0.4,raster=TRUE))
}
dev.off()


CellType = 'EarlyNeuPro'

HypoCellType = paste(CellType,'_Hypo',sep='')

GeCellType = paste(CellType,'_GE',sep='')

CtxCellType = paste(CellType,'_Cortex',sep='')

HypoGenes = rownames(Bulk[[1]])[which(Bulk[[1]][,HypoCellType]>quantile(Bulk[[1]][,HypoCellType],0.95))]

tmphyponet = hyponet[sort(intersect(which(hyponet[,1]%in%HypoGenes),which(hyponet[,2]%in%HypoGenes))),]

tmphyponet$Tissue = 'Hypo'

GeGenes = rownames(Bulk[[1]])[which(Bulk[[1]][,GeCellType]>quantile(Bulk[[1]][,GeCellType],0.95))]

tmpGEnet = GEnet[sort(intersect(which(GEnet[,1]%in%GeGenes),which(GEnet[,2]%in%GeGenes))),]

tmpGEnet$Tissue = 'GE'

CtxGenes = rownames(Bulk[[1]])[which(Bulk[[1]][,CtxCellType]>quantile(Bulk[[1]][,CtxCellType],0.95))]

tmpCTXnet = CTXnet[sort(intersect(which(CTXnet[,1]%in%CtxGenes),which(CTXnet[,2]%in%CtxGenes))),]

tmpCTXnet$Tissue = 'CTX'

tmpNet = rbind(tmphyponet,tmpGEnet,tmpCTXnet)

tmpNet = tmpNet[order(tmpNet$weight,decreasing=TRUE),]

tmpNet = tmpNet[1:10000,]

colorInd = data.frame(ID=unique(tmpNet$ID),Tissues = NA,color=NA)

for(i in 1:nrow(colorInd)){
 colorInd[i,'Tissues'] = paste(sort(tmpNet$Tissue[tmpNet$ID%in%colorInd[i,'ID']]),collapse='_')
if(i%%1000==0) cat(paste(i,', ',sep=''))
}

colorInd[which(colorInd$Tissues=='Hypo'),'color'] = "blue"
colorInd[which(colorInd$Tissues=='GE'),'color'] = "green3"
colorInd[which(colorInd$Tissues=='CTX'),'color'] = "red"
colorInd[which(colorInd$Tissues=='CTX_GE'),'color'] = "gold2"
colorInd[which(colorInd$Tissues=='CTX_Hypo'),'color'] = "purple"
colorInd[which(colorInd$Tissues=='CTX_GE_Hypo'),'color'] = "grey60"
colorInd[which(colorInd$Tissues=='GE_Hypo'),'color'] = "cyan"
tmpNet$color = colorInd$color[match(tmpNet$ID,colorInd$ID)]

assign(paste(CellType,'_Net',sep=''),tmpNet)
tmpNet2 = unique(tmpNet[,c('regulatoryGene','targetGene','color')])

TopTF = names(sort(table(tmpNet2[,1]),decreasing=TRUE)[1:20])

## plot 
nodes=data.frame(name=unique(c(tmpNet2[,1],tmpNet2[,2])))

nodes$type = 'target'
nodes[nodes$name%in%tmpNet2[,1],'type']='TF'
nodes$size=0.01
nodes[nodes$name%in%tmpNet2[,1],'size']=1
nodes[nodes$name%in%TopTF,'size']=5
relations = data.frame(from=tmpNet2[,1],to=tmpNet2[,2],color=tmpNet2$color)
net.tmp = graph_from_data_frame(d=relations, vertices=nodes, directed=F)

V(net.tmp)$label <- ""
E(net.tmp)$arrow.mode <- 0
#V(net.tmp)$size <- 8

nodelabels = rep("",nrow(nodes))

nodelabels[which(nodes[,1]%in%names(sort(table(tmpNet2[,1]),decreasing=TRUE)[1:20]))] = as.character(nodes[which(nodes[,1]%in%names(sort(table(tmpNet2[,1]),decreasing=TRUE)[1:20])),1])

l <- layout_with_fr(net.tmp)

pdf(paste('./TestPlots/All_CtxHyGE_mt10n500_testNetwork_',CellType,'_with_fr.pdf',sep=''),width=12,height=8)
plot(net.tmp,edge.color=relations$color,vertex.label=nodelabels,layout=l,edge.curved=0.2,edge.width=0.5,vertex.label.dist=1)
dev.off()

##############################
### Calculate Gene modules ###
##############################

##############################
### Figure Sup 6 ###
##############################

load('./Analysis/All_CtxHyGE_Clustering_mt10n500_DifGenes.rda')

obj2 = subset(AllGE4K.integrated,features=DifGenes)

DefaultAssay(obj2) = 'RNA'
cellinfo = obj2@meta.data[,c(4,10,6)]
colnames(cellinfo)[2]='seurat_clusters'
id='CxtGeHy_mt10n500'
k=100

metaplot=c('sample','seurat_clusters','Maj_Clust_From_Ref') ## metadata to plot from the obj2

geneinfo = data.frame( genes = rownames(obj2) )
rownames(geneinfo) = rownames(obj2)

counts = obj2@assays$RNA@counts

logcpm = cpm( counts , log = T )
datExpr=sweep(logcpm,1,apply(logcpm,1,mean),"-")
indx.sd=(apply(logcpm,1,sd))==0 # these will produce NAs
datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
datExpr[indx.sd,]=0
if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
#kmeans

cl=kmeans(x=datExpr,centers=k,iter.max=10000,nstart=10,alg="Lloyd")
saveRDS(cl , file = paste('./Analysis/',id,'_8359DifGenes_k',k,'CoreCellGene.rds',sep='') )

clorig=cl

minModSize = 10
kAEtoStay = 0.3
cutHeight = 0.25

#clusters = readRDS('k25.5perc.rds')
clusters = cl

#datExpr = readRDS('cicero_agg_counts.rds')
datExpr = counts
datExpr = as.matrix( datExpr )
datExpr = datExpr[ names(clusters$cluster) , ]

kAE = cor( t(clusters$centers) , t(datExpr) )

colors = clusters$cluster

for( i in 1:length(colors) ) {
  r.i = kAE[ colors[i] , i ]
  if( r.i < kAEtoStay ) colors[i] = 0
}

size = table( colors )
too.small = as.numeric(names(which(size<minModSize)))
colors[ colors %in% too.small ] = 0


centers = sapply( sort(unique(colors)) , function(i)
  colMeans(datExpr[ colors == i , ]) )

colnames(centers) = paste( 'AE' , sort(unique(colors)) , sep = '' )

r = cor( centers )

d = as.dist( 1 - r )

cl = cutree( hclust( d , method = 'average' ) , h = cutHeight )

mergeColors = rep(NA,length(colors))
for( i in 1:max(cl) ) {
  idx = as.numeric( gsub( 'AE','', names(cl)[ cl == i ] ))
  mergeColors[ colors %in% idx ] = i
}
mergeColors = mergeColors - 1
names(mergeColors) = names(colors)

saveRDS( mergeColors , file = paste('./Analysis/',id,'_mergeColors._8359DifGenes_k',k,'.cutHeight0.25_CoreCellGene.rds',sep='') )

AE = sapply( sort(unique(mergeColors)) , function(i)
  colMeans(datExpr[ mergeColors == i , ]) )
colnames(AE) = paste( 'AE' , 0:max(mergeColors) , sep = '' )

kAE = cor( t(datExpr) , AE )
sort( kAE[ mergeColors == 1 , 'AE1' ] , decreasing = T )[1:50]

saveRDS( kAE , file = paste('./Analysis/',id,'_mergeColors_8359DifGenes_k',k,'.cutHeight0.25.kAE_CoreCellGene.rds',sep='') )

mtx = as.matrix(GetAssayData(obj2))
mtx = mtx[ names(colors) , ] ##genes intersecting with genie3 (essentially is genie3 genes )
datExpr2 = t(mtx)

MEs = moduleEigengenes( datExpr2 , colors , softPower = 1 , excludeGrey = T )

saveRDS( MEs , file = paste('./Analysis/',id,'_mergeColors_8359DifGenes.k',k,'.MEs.rds',sep='') )


for(n in colnames(MEs[['eigengenes']])){

obj2@meta.data[,n] <- MEs[['eigengenes']][colnames(obj2),n]
cat(paste(n,', ',sep=''))
} 

pdf(paste('./TestPlots/',id,'_Kmean_8359DifGenes',ncol(MEs[['eigengenes']]),'groupsEigengenes.pdf',sep=''),width=12,height=8)
#print(DimPlot(obj2, label = TRUE) + NoLegend())
for(x in metaplot){
print(DimPlot(obj2, reduction = "umap", group.by = x, label = TRUE, repel = TRUE,raster=TRUE))
}
for(i in colnames(MEs[['eigengenes']])){
print(FeaturePlot(obj2, features = i, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red'),raster=TRUE))
}
dev.off()













