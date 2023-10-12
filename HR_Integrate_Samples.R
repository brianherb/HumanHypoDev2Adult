## to replicate conda environment used in paper, install conda, and then call the following in terminal: conda env create -f environment_droplet.yml - yml file here: https://github.com/brianherb/HumanHypothalamusDev/blob/main/environment_droplet.yml

library(Seurat)
library(limma)

devtools::source_url('https://github.com/brianherb/HumanHypoDev2Adult/blob/main/HR_functions_reference.R?raw=TRUE')

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

######################
### load datasets ####
######################

## Lein hypo samples 
#load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoSamples_PostSCT.rda')) #EdHypoNormDat

load('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/SeuratObj/HSatlasSamples_PostSCT.rda') # HSatlasNormDat

## Kriegstein data 
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_PostSCT.rda')) #HypoNormDat

load('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/SeuratObj/HypoSamples_PostSCT.rda') 

## Non-hypothalamus Human Embryonic samples 

## GW18 Cortex
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18CtxSamples_PostSCT.rda')) #GW18CtxNormDat

## GW18 GE
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18GESamples_PostSCT.rda')) #GW18GENormDat

## GW19
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19_NonHySamples_PostSCT.rda')) # GW19_NormDat

## GW20
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34_NonHySamples_PostSCT.rda')) # GW20_34_NormDat

## Zhou data 

load('/local/projects-t3/idea/bherb/Hypothalamus/Zhou/SeuratObj/ZhouSamples_PostSCT.rda') # ZhouNormDat

## Mouse atlas 

load('./SeuratObj/Ref2Samples_PostSCT.rda') #Ref2NormDat

load('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/SeuratObj/RefSamples_PostSCT.rda') #RefNormDat

############################################
### Integrate CS13 - GW25, all cells    ####
############################################

combDat = HypoNormDat[c("CS13","CS14","CS15","CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")]  
IntName = 'HypoCS13_GW25_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','NEUROD6','SLC17A6')
#40927 cells in integrated object

HypoCS13_GW25_mt10 = readRDS('./SeuratObj/HypoCS13_GW25_mt10_integrated.rds')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoCS13_GW25_mt10_integrated.rds'

######################################################
### Integrate All Lein adult samples, all cells   ####
######################################################
 
combDat = HSatlasNormDat	  
IntName = 'HSatlas_mt10' 
PrintColCat = c('SampleID','Donor','Chemistry','Age','orig.ident','Roi','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoAll_mt10_integrated.rds' - used for determining Neuronal vs Non-neuronal cells

######################################################
### Integrate all Lein adult samples, neurons     ####
######################################################


HSatlas.integrated = readRDS('./SeuratObj/HSatlas_mt10_integrated.rds')

DefaultAssay(HSatlas.integrated) = 'RNA'

pdf('./TestPlots/HSatlas_CellTypeMarkers_SC_VlnPlots.pdf',width=12,height=8)
for(i in c('GAD1','SLC17A6','SLC32A1','PTPRC','CLDN5','LUM','SOX10','AQP','FOXJ1')){
print(VlnPlot(HSatlas.integrated, features = i, group.by='seurat_clusters'))
}
dev.off()

combDat = HSatlasNormDat  

exInd = unique(which(is.na(match(HSatlas.integrated@meta.data$seurat_clusters,c(2,4,12,13,38,18,19,41,36))))) 

HSatlasNeuronBarcodes_mt10 = colnames(HSatlas.integrated)[exInd]

save(HSatlasNeuronBarcodes_mt10,file='./SeuratObj/HSatlasNeuronBarcodes_mt10.rda',compress=TRUE) 


for(i in names(combDat)){
combDat[[i]] = subset(combDat[[i]],cells=intersect(colnames(combDat[[i]]),HSatlasNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

combDat = combDat[-c(13:14)] #10X362_5 , 10X362_6 -low number of cells are neurons  - skip 

IntName = 'HSatlasNeuro_mt10'  
PrintColCat = c('SampleID','Roi','Age') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')


###########################################################
### Integrate all Kriegstein neurons, starting at CS22 ####
###########################################################

## this is done after creating 'HypoCS13_GW25_mt10' - using seurat groups that express GAD1/2 and SLC17A6 (and are more mature - ignoring progenitors )

combDat = HypoNormDat[c("CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")] 

#KaNeuronBarcodes_mt10 = colnames(HypoCS13_GW25_mt10)[which(!is.na(match(HypoCS13_GW25_mt10@meta.data$seurat_clusters,c(0,4,15,19,20,25))))]

#save(KaNeuronBarcodes_mt10,file='./SeuratObj/KaNeuronBarcodes_mt10.rda',compress=TRUE)


HannahCells=read.csv('./Analysis/CleanedClusters_Figure1_19DEC22.csv',row.names=1)

HannahNeuronBarcodes_mt10 = HannahCells[which(HannahCells$Pop=='Neuronal'),'Row.names']



## KaNeuronBarcodes_mt10 saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/KaNeuronBarcodes_mt10.rda'

for(i in names(combDat)){
combDat[[i]] = subset(combDat[[i]],cells=intersect(colnames(combDat[[i]]),KaNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

IntName = 'KaHypoNeuro_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/KaHypoNeuro_mt10_integrated.rds'

##############################################
### Integrate Zhou, Kriegstein and Lein neurons ####
##############################################

HannahCells=read.csv('./Analysis/CleanedClusters_Figure1_19DEC22.csv',row.names=1)

HannahNeuronBarcodes_mt10 = HannahCells[which(HannahCells$Pop=='Neuronal'),'Row.names']

combDatKA = HypoNormDat[c("CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")] 

for(i in names(combDatKA)){
combDatKA[[i]] = subset(combDatKA[[i]],cells=intersect(colnames(combDatKA[[i]]),HannahNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

combDatEd = HSatlasNormDat[-c(13:14)] #10X362_5 , 10X362_6 -low number of cells are neurons  - skip  

load('./SeuratObj/HSatlasNeuronBarcodes_mt10.rda') #HSatlasNeuronBarcodes_mt10

for(i in names(combDatEd)){
combDatEd[[i]] = subset(combDatEd[[i]],cells=intersect(colnames(combDatEd[[i]]),HSatlasNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}


for(i in names(combDatEd)){
combDatEd[[i]]$sample = combDatEd[[i]]$SampleID
cat(paste(i,', ',sep=''))
}


## Zhou neurons - GW10 and up? GW7 is about 25% neurons, but GW8 is 80%

combDatZhou = ZhouNormDat[-c(1:4)] # Skip GW8 GW9

for(i in names(combDatZhou)){
combDatZhou[[i]] = subset(combDatZhou[[i]],cells=intersect(colnames(combDatZhou[[i]]),HannahNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
} ## smallest sample - 78 cells 


combDat = c(combDatKA,combDatZhou,combDatEd)

for(i in names(combDat)){
if(any(colnames(combDat[[i]])%in%NAcells$X)){
combDat[[i]] = subset(combDat[[i]],cells=setdiff(colnames(combDat[[i]]),NAcells$X))
}
cat(paste(i,', ',sep=''))
}

IntName = 'EdKaZhouHypoNeurons_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdKaHypoNeurons_mt10_3D_integrated.rds'

##############################################
### Integrate Kriegstein and Zhou all cells ##
##############################################

for(i in names(HypoNormDat)){
	HypoNormDat[[i]]@meta.data$Study = 'Kriegstein'
}

for(i in names(ZhouNormDat)){
	ZhouNormDat[[i]]@meta.data$Study = 'Zhou'
}

combDat = c(HypoNormDat,ZhouNormDat)

IntName = 'KaZhouAll_mt10'  
PrintColCat = c('sample','Study','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')


KaZhouAll_mt10 = readRDS('./SeuratObj/KaZhouAll_mt10_integrated.rds')


##############################################
### Mouse pub data 				  		  ####
##############################################

names(RefNormDat) = gsub('_','-',names(RefNormDat))
names(RefNormDat) = gsub('.','-',names(RefNormDat),fixed=TRUE)


for(i in names(RefNormDat)){
	RefNormDat[[i]]@meta.data$Study = flexsplit(i,'-')[1]
}

for(i in names(Ref2NormDat)){
	Ref2NormDat[[i]]@meta.data$Study = flexsplit(i,'-')[1]
}

combDat = c(RefNormDat,Ref2NormDat)

IntName = 'MouseRefAll_mt10'  
PrintColCat = c('sample','Study','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## neurons marked by GAD2 and SLC17A6 are in clusters - 1,2,5,6,7,11,12,13,15,16,21,23,25,28

MouseRefAll.integrated =readRDS('./SeuratObj/MouseRefAll_mt10_integrated.rds')

exInd = unique(which(is.na(match(MouseRefAll.integrated@meta.data$seurat_clusters,c(0,3,4,8,9,10,14,17,18,19,20,22,24,26,27,29,30,31))))) # #  137483 neurons in c(1,2,5,6,7,11,12,13,15,16,21,23,25,28)

MouseRefAllNeuronBarcodes_mt10 = colnames(MouseRefAll.integrated)[exInd] 

save(MouseRefAllNeuronBarcodes_mt10,file='./SeuratObj/MouseRefAllNeuronBarcodes_mt10.rda',compress=TRUE) 

##############################################
### Mouse and Human neurons 	  		  ####
##############################################

## Mouse References 
combDatMMsc = c(RefNormDat,Ref2NormDat)

for(i in names(combDatMMsc)){
	tmpcells = intersect(colnames(combDatMMsc[[i]]),MouseRefAllNeuronBarcodes_mt10)
if(length(tmpcells)>0){
combDatMMsc[[i]] = subset(combDatMMsc[[i]],cells=tmpcells)
} else {
combDatMMsc[[i]] = NA
}
cat(paste(i,', ',sep=''))
}


naInd=which(is.na(combDatMMsc))

combDatMMsc = combDatMMsc[-naInd] # Black-E15L

for(i in names(combDatKA)){
combDatKA[[i]]$study='DevHuKa'
}

for(i in names(combDatZhou)){
combDatZhou[[i]]$study='DevHuZhou'
}

for(i in names(combDatEd)){
combDatEd[[i]]$study='AdultHu'
}

## remove NA cells from Lein samples 

## cells Hannah deamed as Nucleus accumbens 

NAcells = read.csv("./Analysis/NACBarcs.csv")


combDat = c( combDatMMsc, combDatKA,combDatZhou,combDatEd) #combDatLow,

for(i in names(combDat)){
if(any(colnames(combDat[[i]])%in%NAcells$X)){
combDat[[i]] = subset(combDat[[i]],cells=setdiff(colnames(combDat[[i]]),NAcells$X))
}
cat(paste(i,', ',sep=''))
}

for(i in c(1:length(combDat))){
	if(i==1){
totsum=ncol(combDat[[i]]) 
} else {
totsum=totsum + ncol(combDat[[i]]) 
	}
}



IntName = 'MusRefEdKaZhouHypoNeurons_mt10'  
PrintColCat = c('study','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6','POMC')


## calc per mt for 

for(i in 1:length(combDat)){

combDat[[i]] <- PercentageFeatureSet(combDat[[i]], pattern = "^MT-", col.name = "percent.mt")

}

##############################################
### Mouse and Human POMC neurons 		  ####
##############################################

#####  based on Hannah's assignments from cluster 13 


HG_POMCneurons = read.csv('./Analysis/Hannah-POMCAssignments_2_7_23.csv',row.names=1)

combDat = c( combDatMMsc, combDatKA,combDatZhou,combDatEd)

for(i in names(combDat)){
	tmpcells = intersect(colnames(combDat[[i]]),HG_POMCneurons[,1])
if(length(tmpcells)>0){
combDat[[i]] = subset(combDat[[i]],cells=tmpcells)
} else {
combDat[[i]] = NA
}
cat(paste(i,', ',sep=''))
}


naInd=which(is.na(combDat))

combDat = combDat[-naInd] # Black-E15L

#combDat = c(combDat,combDatLow)


IntName = 'Hannah_MMHSpomcNeurons_mt10'  
PrintColCat = c('study','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6','POMC')


##############################################
### Integrate Ctx GE Hypo 				  ####
##############################################

## Done in Fig4_HumanHypoCtxGE.R 

############################################
### Common code for integration 	    ####
############################################

## count cells 

for(i in 1:length(combDat)){
if(i==1){
	cellCounts = ncol(combDat[[i]])
} else {
	cellCounts = c(cellCounts,ncol(combDat[[i]]))
}
}


if(length(which(cellCounts<=20))>0){
	combDat = combDat[-which(cellCounts<=20)]
} ## set to 30 for POMC integration , set to 20 for Hannah POMC


for (i in 1:length(combDat)) {
	DefaultAssay(combDat[[i]]) = 'RNA'
	combDat[[i]] <- subset(combDat[[i]], percent.mt<=10)
	if('FOXG1' %in% rownames(combDat[[i]]@assays$RNA) ){
		combDat[[i]] <- subset(combDat[[i]], cells=colnames(combDat[[i]])[which(combDat[[i]]@assays$RNA@counts['FOXG1',]==0)])
	}
	if('NEUROD6' %in% rownames(combDat[[i]]@assays$RNA) ){
		combDat[[i]] <- subset(combDat[[i]], cells=colnames(combDat[[i]])[which(combDat[[i]]@assays$RNA@counts['NEUROD6',]==0)])
	}
    combDat[[i]] <- NormalizeData(combDat[[i]], verbose = FALSE)
    combDat[[i]] <- FindVariableFeatures(combDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    cat(paste(i,', ',sep=''))
}



features <- SelectIntegrationFeatures(object.list = combDat)

for (i in 1:length(combDat)) {
 combDat[[i]]<- ScaleData(combDat[[i]], features = features, verbose = FALSE)
    combDat[[i]] <- RunPCA(combDat[[i]], features = features, verbose = FALSE,npcs=15 ) #,npcs=15 for POMC integration
}

All.anchors <- FindIntegrationAnchors(object.list = combDat, anchor.features = features, reduction = "rpca",dims=1:15,k.score=15,k.anchor=15,k.filter=15) # for Mouse (HypoMap)+Human ,reference=c(133,154,99), ### for all mouse ref (no low) + human c(60,33,2,68,91,93) ## for all mm, ,reference=c(2,33,60) ## for Hannah pomc neurons - reduction = "rpca",dims=1:15,,k.score=15,k.anchor=15,k.filter=15
All.integrated <- IntegrateData(anchorset = All.anchors, dims = 1:15, k.weight =15) # ## dims = 1:30, k.weight =50 for integrating Ka, Ed, Zhou human neurons ,  for all MM + HS, with  "Romanov-P23" "Black-P45"   "Chen"        "GW16_hypo"   "10X192_7"    "10X193_3"  as reference  k. weight = 30, dims = 1:15  / for all MM, Romanov P23, Blackshaw/Kim P45, Chen , k. weight = 30, dims = 1:15 , for pomc integration, k.weight=15
DefaultAssay(All.integrated) <- "integrated"
All.integrated <- ScaleData(All.integrated, verbose = FALSE)
All.integrated <- RunPCA(All.integrated, npcs = 30, verbose = FALSE)
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:30)

All.integrated <- FindNeighbors(All.integrated, dims = 1:30, verbose = FALSE)
All.integrated <- FindClusters(All.integrated, verbose = FALSE, resolution=0.5)
saveRDS(All.integrated,file=paste('./SeuratObj/',IntName,'_integrated.rds',sep=''),compress=TRUE) 

DefaultAssay(All.integrated) = 'RNA'

pdf(file=paste("./TestPlots/",IntName,"_Check_Clustering.pdf",sep=''),width=12,height=8)
for(k in PrintColCat){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
for(m in PrintColCont){
print(FeaturePlot(All.integrated, features = m, pt.size = 0.2))
}

dev.off()
