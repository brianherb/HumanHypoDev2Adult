#  suggested comp env - mem_free=150G thread 4
## conda environment can be created from https://github.com/brianherb/HumanHypoDev2Adult/blob/main/Hypo_R_Env.yml

library(CoGAPS)
library(Seurat)

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

humanNeu = readRDS('./SeuratObj/EdKaZhouHypoNeurons.rds')

## create single cell experiment 
DefaultAssay(humanNeu) = 'RNA'
TMPdat = as.matrix(humanNeu@assays$RNA[rownames(humanNeu@assays$integrated),])

rm(humanNeu,TMPsce)
gc()

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

pat1=xxCoGAPS@sampleFactors

colnames(pat1) = paste0('CoGap1_',colnames(pat1))

humanNeu@meta.data= cbind(humanNeu@meta.data,pat1[colnames(humanNeu),])

DefaultAssay(humanNeu)='RNA'

pdf(file="./TestPlots/EdKaZhouHypoNeurons_CoGaps1_15pat.pdf",width=12,height=8)
print(DimPlot(humanNeu, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))

for(i in c('GAD2','SLC17A6',colnames(pat1))){
print(FeaturePlot(humanNeu, features = i, pt.size = 0.2))
}
dev.off()

write.csv(humanNeu@meta.data[,grep('CoGap1_Pattern',colnames(humanNeu@meta.data))],file='./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter_sampleFactors.csv')
write.csv(xxCoGAPS@featureLoadings,file='./Analysis/EdKaZhouHypoNeurons_CoGAPS_15pat_10Kiter_featureLoadings.csv')











