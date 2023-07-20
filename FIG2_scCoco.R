# conda activate /autofs/burnsfs/projects-t3/idea/amentlab_software/conda_envs/r_4.1.1

library(scCoco)
library(Seurat)


## average expression per cluster? 

## HSatlas 

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

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



HSatlasNeuron_mt10 = readRDS('./SeuratObj/HSatlasNeuro_mt10_integrated.rds')



## functions:

#  aba_ids

testGenes = unique(na.omit(HStoMM(VariableFeatures(HSatlasNeuron_mt10))))

testABA = aba_ids(testGenes)

## back again to get gene info? 

## need clusters 


humanABA = unique(MMtoHS(names(testABA)))

#HSneu_sc_res = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/HSneu_sc_res.rds'))

treeTrim = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/tree80_L15_treeTrim.rds'))


HSatlasNeuron_mt10@meta.data$K55 = treeTrim[colnames(HSatlasNeuron_mt10),'K55']


humanExp = AverageExpression(
  HSatlasNeuron_mt10,
  assays = 'RNA',
  features = humanABA,
  return.seurat = FALSE,
  group.by = "K55")

rownames(humanExp[['RNA']]) = HStoMM(rownames(humanExp[['RNA']]))

K55list=vector(mode='list',length=ncol(humanExp[['RNA']]))

names(K55list) = paste0('clust_',colnames(humanExp[['RNA']]))

for(k in 1:length(K55list)){

testVec = as.vector(humanExp[['RNA']][,k])

names(testVec) = rownames(humanExp[['RNA']])

K55list[[k]] = names(testVec)[order(testVec,decreasing=T)[1:20]]

}


#test=coco_expression(testABA[1]) ## 67 x 2378

#test2=region_annotation(t(humanExp))

test3 = findRegions_genesets(gene_set=K55list)

test3[["scores_per_leaf_region_all"]][order(test3[["scores_per_leaf_region_all"]]$clust_1,decreasing=FALSE),][1:10,1:10]


test4=summariseRegions_genesets(test3)


## Run (or pull) DEG information 

DefaultAssay(HSatlasNeuron_mt10) = 'RNA'

HSatlasNeuron_mt10@active.ident=as.factor(HSatlasNeuron_mt10$K55)

testDEG = FindAllMarkers(HSatlasNeuron_mt10,features=humanABA,verbose=FALSE)

saveRDS(testDEG,'./Analysis/HSneu_K55_Cluster_DEG.rds')


testDEG = readRDS('./Analysis/HSneu_K55_Cluster_DEG.rds')


testGenes = unique(na.omit(HStoMM(testDEG$gene))) #618

testABA = aba_ids(testGenes)


humanABA = unique(MMtoHS(names(testABA)))


## for each cluster, grab DEG's (pos sig)

K55clusts=paste0('clust_',unique(testDEG$cluster))

K55deg=vector(mode='list',length=length(K55clusts))

names(K55deg) = K55clusts



for(i in K55clusts){

tmpDEG = testDEG[which(testDEG$cluster==as.numeric(gsub('clust_','',i)) & testDEG$avg_log2FC>0 & testDEG$p_val_adj<=0.05),]


tmpGene = unique(HStoMM(tmpDEG$gene))

K55deg[[i]] = tmpGene


cat(paste0(i,', '))

}


test1 = findRegions_genesets(K55deg,target_structure_id = "1097",exclude_ids = c("338"),target_level = "8")

saveRDS(test1,file='./Analysis/scCoco_K55.rds')


#,target_structure_id = "997",exclude_ids = c("338"),target_level = "8"

test1[["scores_per_leaf_region_all"]][order(test1[["scores_per_leaf_region_all"]]$clust_1,decreasing=FALSE),][1:10,1:10]


test2=summariseRegions_genesets(test1)

## check against Hannah's 

nuclei = read.csv('./Analysis/AllHumanAnnots_ByCell_12JUL23.csv')


K55= unique(nuclei[,c('K55','K55_Class','K55Nuclei')])

rownames(K55) = K55$K55

K55 = K55[order(K55$K55,decreasing=FALSE),]

rownames(test2) = gsub('clust_','',rownames(test2))

test2$Hannah_nuclei = K55$K55Nuclei

test2 = test2[,c(1,7,2:6)]

test2$K55_cluster = rownames(K55)

saveRDS(test2,file='./Analysis/scCoco_K55_wHannahNuclei.rds')

write.csv(test2,file='./Analysis/scCoco_K55_wHannahNuclei.csv',quote=FALSE)


saveRDS(K55deg,file='./Analysis/scCoco_K55_DEGs.rds')

saveRDS(humanABA,file='./Analysis/scCoco_K55_humanABAgenes.rds')





