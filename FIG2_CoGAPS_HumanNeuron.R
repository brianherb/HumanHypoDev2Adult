
# qlogin -P "owhite-nemo" -l mem_free=150G -q interactive.q -pe thread 4

## conda activate /autofs/burnsfs/projects-t3/idea/amentlab_software/conda_envs/r_4.1.1



library(CoGAPS)
library(Seurat)



flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

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



setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

humanNeu = readRDS('./SeuratObj/EdKaZhouHypoNeurons_mt10_integrated.rds')





#CoGAPS(GIST.matrix, nIterations=10000, outputFrequency=5000, nThreads=1, seed=5)

#params <- setDistributedParams(params, nSets=3)

#CoGAPS(GISTCsvPath, params, distributed="single-cell", messages=FALSE, transposeData=TRUE, nIterations=1000)

## create single cell experiment 


DefaultAssay(humanNeu) = 'RNA'

#TMPsce = Seurat::as.SingleCellExperiment(humanNeu)

#@assays$RNA[rownames(humanNeu@assays$integrated)]

TMPdat = as.matrix(humanNeu@assays$RNA[rownames(humanNeu@assays$integrated),])

### from Carlo 

rm(humanNeu,TMPsce)
gc()

#library("CoGAPS")

nThrd=16 # was 1
cells <- colnames(TMPdat)
genes <- rownames(TMPdat)
params <- new("CogapsParams")
params <- setParam(params, "sparseOptimization", FALSE)
params <- setParam(params, "nIterations", 10000) #was 40000
params <- setParam(params, "nPatterns", 15) # was 30
params <- setDistributedParams(params, nSets=1)

xxCoGAPS <- CoGAPS(data=TMPdat,params=params,geneNames=genes,sampleNames=cells,messages=TRUE,transposeData=FALSE,nPatterns=15,nIterations=10000,nThreads=16)



save(xxCoGAPS,file='./Analysis/HSatlas_CoGAPS_15Pat_10KIter_Interactive.rda')


### plotting 

load('./Analysis/HumanNeu_CoGAPS_15Pat_10KIter_Interactive.rda') # xxCoGAPS


xxCoGAPS2= readRDS('./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter.rds')

pat1=xxCoGAPS@sampleFactors

colnames(pat1) = paste0('CoGap1_',colnames(pat1))

pat2=xxCoGAPS2@sampleFactors

colnames(pat2) = paste0('CoGap2_',colnames(pat2))



humanNeu@meta.data= cbind(humanNeu@meta.data,pat1[colnames(humanNeu),])

humanNeu@meta.data= cbind(humanNeu@meta.data,pat2[colnames(humanNeu),])


DefaultAssay(humanNeu)='RNA'

pdf(file="./TestPlots/EdKaZhouHypoNeurons_CoGaps1_15pat.pdf",width=12,height=8)
print(DimPlot(humanNeu, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
#print(DimPlot(KaHypoNeuro_mt10, reduction = "umap", group.by = 'TopRegion', label = TRUE, repel = TRUE))

for(i in c('GAD2','SLC17A6',colnames(pat1))){
print(FeaturePlot(humanNeu, features = i, pt.size = 0.2))

}

dev.off()




pdf(file="./TestPlots/EdKaZhouHypoNeurons_CoGaps2_15pat.pdf",width=12,height=8)
print(DimPlot(humanNeu, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
#print(DimPlot(KaHypoNeuro_mt10, reduction = "umap", group.by = 'TopRegion', label = TRUE, repel = TRUE))

for(i in c('GAD2','SLC17A6',colnames(pat2))){
print(FeaturePlot(humanNeu, features = i, pt.size = 0.2))

}

dev.off()

 write.csv(humanNeu@meta.data[,grep('CoGap1_Pattern',colnames(humanNeu@meta.data))],file='./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter_sampleFactors.csv')

 write.csv(xxCoGAPS@featureLoadings,file='./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter_featureLoadings.csv')




 write.csv(humanNeu@meta.data[,grep('CoGap2_Pattern',colnames(humanNeu@meta.data))],file='./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter_JustPat.csv')


### try 3D plots and lineage trees 


### reconstructing 

## EdKaZa_3D = read.csv('./SeuratObj/EdKaZhouNeu_3D_coordinates.csv',row.names=1)

## colnames(EdKaZa_3D) = c('UMAP_1','UMAP_2','UMAP_3')

## EdKaZa_3D = as.matrix(EdKaZa_3D)

## humanNeu_3D = humanNeu

## humanNeu_3D@reductions$umap@cell.embeddings = EdKaZa_3D

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = humanNeu_3D, vars = c("UMAP_1", "UMAP_2", "UMAP_3", paste0('CoGap1_Pattern_',seq(from=1,to=15,by=1))))


for(i in paste0('CoGap1_Pattern_',seq(from=1,to=15,by=1))){

p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = paste0('~',i),type = "scatter3d", mode= "markers",  marker = list(size = 1, width=1),colors='Blues')

htmlwidgets::saveWidget(p, paste0("./TestPlots/EdKaZhouNeu_mt10_3D_",i,".html"))

}




p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Zhou_SC,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, "./TestPlots/KaZhouAll_mt10_3D_Zhou_SC.html")














### can split into another script - but try on all embryonic samples, with more cores 


# qlogin -P "owhite-nemo" -l mem_free=200G -q interactive.q -pe thread 16


KaZhouAll_mt10 = readRDS('./SeuratObj/KaZhouAll_mt10_integrated.rds')


DefaultAssay(KaZhouAll_mt10) = 'RNA'

#TMPsce = Seurat::as.SingleCellExperiment(KaZhouAll_mt10)

#@assays$RNA[rownames(KaZhouAll_mt10@assays$integrated)]

TMPdat = as.matrix(KaZhouAll_mt10@assays$RNA[rownames(KaZhouAll_mt10@assays$integrated),])

### from Carlo 

rm(KaZhouAll_mt10)
gc()

#library("CoGAPS")

nThrd=1 # was 1
cells <- colnames(TMPdat)
genes <- rownames(TMPdat)
params <- new("CogapsParams")
params <- setParam(params, "sparseOptimization", FALSE)
params <- setParam(params, "nIterations", 10000) #was 40000
params <- setParam(params, "nPatterns", 15) # was 30
params <- setDistributedParams(params, nSets=1)

xxCoGAPS <- CoGAPS(data=TMPdat,params=params,geneNames=genes,sampleNames=cells,messages=TRUE,transposeData=FALSE,nPatterns=15,nIterations=10000,nThreads=1)

save(xxCoGAPS,file='./Analysis/HSembryo_CoGAPS_15Pat_10KIter_Interactive.rda')


### plotting  - from submitted job 
xxCoGAPS = readRDS('./Analysis/EmbryonicSamp_CoGAPS_15pat_10Kiter.rds')

pat=xxCoGAPS@sampleFactors

colnames(pat) = paste0('CoGap_',colnames(pat))

KaZhouAll_mt10@meta.data= cbind(KaZhouAll_mt10@meta.data,pat[colnames(KaZhouAll_mt10),])

DefaultAssay(KaZhouAll_mt10)='RNA'

pdf(file="./TestPlots/EdKaZhouHypoNeurons_CoGaps_15pat.pdf",width=12,height=8)
print(DimPlot(KaZhouAll_mt10, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
#print(DimPlot(KaHypoNeuro_mt10, reduction = "umap", group.by = 'TopRegion', label = TRUE, repel = TRUE))

for(i in c('GAD2','SLC17A6',colnames(pat))){
print(FeaturePlot(KaZhouAll_mt10, features = i, pt.size = 0.2))

}

dev.off()


pdf(file="./TestPlots/KaZhouHypo_ASCL1_NEUROG2.pdf",width=12,height=8)
print(DimPlot(KaZhouAll_mt10, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
#print(DimPlot(KaHypoNeuro_mt10, reduction = "umap", group.by = 'TopRegion', label = TRUE, repel = TRUE))

for(i in c('ASCL1','NEUROG2','GAD2','SLC17A6')){
print(FeaturePlot(KaZhouAll_mt10, features = i, pt.size = 0.2))

}

dev.off()


