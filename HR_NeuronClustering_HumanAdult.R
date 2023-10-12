library(Seurat)
library(clustree)
library(dplyr)
library(plyr)
library(ggtree)
library(MetaNeighbor)
library(Matrix)
library(dendextend)
library(SingleCellExperiment)
library(scrattch.hicat)
library(phylogram)
library(mrtree)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

HSatlasNeuron_mt10 = readRDS('./SeuratObj/HSatlasNeuro_mt10_integrated.rds')

## what is the range? 
for(i in seq(from=6, to=20, by=1)){
HSatlasNeuron_mt10 <- FindClusters(HSatlasNeuron_mt10, verbose = FALSE, resolution=i)
cat(paste0(i,', '))
}

pdf('./TestPlots/HSatlasNeuron_Check_cluster_resolution_6_20.pdf',width=12,height=8)
for (k in seq(from=6,to=20,by=1)){
print(DimPlot(HSatlasNeuron_mt10, reduction = "umap", group.by = paste("integrated_snn_res.",k,sep=""), label = TRUE, repel = TRUE))
}
dev.off()

HSneu_sc_res=HSatlasNeuron_mt10@meta.data[,grep('integrated_snn_res',colnames(HSatlasNeuron_mt10@meta.data))] 
HSneu_sc_res=HSneu_sc_res[,order(as.numeric(gsub('integrated_snn_res.','',colnames(HSneu_sc_res))))]

saveRDS(HSneu_sc_res,file='./Analysis/HSneu_sc_res.rds')

HSneu_sc_res = readRDS('./Analysis/HSneu_sc_res.rds')

HSTF = read.csv('/local/projects-t3/idea/bherb/annotation/Hsapiens/Human_TF.csv')
## Hannah's assignments for old HiCat

HG_HiCat = read.csv('./Analysis/Hicat_Assignments_4FEB23.csv')


## try multithreading - 24 cores

tree80_L15=mrtree(HSneu_sc_res[,c('integrated_snn_res.0.01','integrated_snn_res.0.1','integrated_snn_res.0.2','integrated_snn_res.0.5','integrated_snn_res.1','integrated_snn_res.5','integrated_snn_res.10','integrated_snn_res.15','integrated_snn_res.20','integrated_snn_res.30','integrated_snn_res.40','integrated_snn_res.50','integrated_snn_res.60','integrated_snn_res.70','integrated_snn_res.80')],n.cores=24)

saveRDS(tree80_L15,file='./Analysis/HSatlasNeu_MRTree_15_level_res.rds')

### Identifying "true" nodes
## are at least two samples present?
## is there a minimum number of DMG's? 

library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 180000 * 1024^2)

DefaultAssay(HSatlasNeuron_mt10)='RNA'
dir.create('./Analysis/tree80_L15_DEGs')

## expressed genes and no mt, RPS/RPL
skipgenes = unique(rownames(HSatlasNeuron_mt10)[c(which(rowSums(HSatlasNeuron_mt10@assays$RNA@counts)==0),grep('^MT-',rownames(HSatlasNeuron_mt10)),grep('^MT-',rownames(HSatlasNeuron_mt10)),grep('^MT-',rownames(HSatlasNeuron_mt10)))])

testgenes = setdiff(rownames(HSatlasNeuron_mt10),skipgenes)
treePre=tree80_L15$labelmat.mrtree
#cellbc = rownames(treeTest)
treeTrim=treePre
treeComment = treePre
treeComment[,] = NA
treeDonorStat = treePre
treeDonorStat[,] = NA

## DB = Donor Bias (>80%)
## MB = Lacking gene expression differences - "merge branch"
## SC = Single child 
## GE = true branch - sufficent DEG

treeDonor = HSatlasNeuron_mt10@meta.data[rownames(treeComment),'Donor']

for( i in rev(1:(ncol(treePre)-1))){
HSatlasNeuron_mt10@meta.data$tmpNodes = treePre[colnames(HSatlasNeuron_mt10),i+1]
HSatlasNeuron_mt10@active.ident=as.factor(HSatlasNeuron_mt10$tmpNodes)
## which are 1:1 
childCount = tapply(X=treePre[,i+1],IND=treePre[,i],FUN=function(x){length(unique(x))})
singleChild = which(childCount[as.character(treePre[,i])]==1)
treeComment[singleChild,i] = "SC"

## donor bias 
tmpDonor = tapply(X=treeDonor,IND=treePre[,i+1],FUN=function(x){any(prop.table(table(x))>0.80)})
treeDonorStat[which(tmpDonor[as.character(treePre[,i+1])]),i+1] = "DB"
geneTest = c(1:nrow(treePre))[-which(treeComment[,i]=="SC")]

checkNodes = unique(treePre[geneTest,i])

for( j in checkNodes){
tmpInd = which(treePre[,i]==j)

numChild = childCount[as.character(j)]

if(numChild==2){

tmpid = paste(treePre[tmpInd[1],1:i],collapse='_')
tmpChild = unique(treePre[tmpInd,i+1])
tmpInd1 = which(treePre[tmpInd,i+1]==tmpChild[1])
tmpInd2 = which(treePre[tmpInd,i+1]==tmpChild[2])

if(min(c(length(tmpInd1),length(tmpInd2)))<5){
treeComment[tmpInd,i+1] = "MB"
next
} ## need at least 5 cells in 1:1 comp

testDEG = FindMarkers(HSatlasNeuron_mt10,ident.1=names(tmpInd1),ident.2=names(tmpInd2),features=testgenes)

saveRDS(testDEG,file=paste0('./Analysis/tree80_L15_DEGs/DEG_',tmpid,'.rds'))
lenDEG = length(which(abs(testDEG$avg_log2FC)>=0.5 & testDEG$p_val_adj<=0.05))

if(lenDEG>=10){
treeComment[tmpInd,i] = "GE"
} else {
treeComment[tmpInd,i] = "MB"
}
} else {
## 3 or more
tmpid = paste(treePre[tmpInd[1],1:i],collapse='_')
tmpChild = unique(treePre[tmpInd,i+1])

tmpObj = subset(HSatlasNeuron_mt10,cells=names(tmpInd))
testDEG = FindAllMarkers(tmpObj,features=testgenes,verbose=FALSE)

saveRDS(testDEG,file=paste0('./Analysis/tree80_L15_DEGs/DEG_',tmpid,'.rds'))
lenDEG = length(which(abs(testDEG$avg_log2FC)>=0.5 & testDEG$p_val_adj<=0.05))

if(lenDEG>=(10*length(tmpChild))){
treeComment[tmpInd,i] = "GE"
} else {
treeComment[tmpInd,i] = "MB"
}
} ## end single node work
} #j - nodes within resolution
cat(paste0('\n\n','resolution ',i,' done','\n\n'))
} #i - each resolution

saveRDS(treeComment,'./Analysis/tree80_L15_treeComment.rds')

treeComment = readRDS('./Analysis/tree80_L15_treeComment.rds')
## only thing not considered is for nodes with more than two children I did not count the children 
### eval merging - merge MB, 

treeTrim = treePre

for( i in rev(1:(ncol(treePre)-1))){
childCount = tapply(X=treePre[,i+1],IND=treePre[,i],FUN=function(x){length(unique(x))})
tmpBranches = names(childCount)[which(childCount>1)]
DBind = which(treeDonorStat[,i+1]=='DB')
DBparent = unique(treePre[DBind,i])
dbBranches = intersect(tmpBranches,DBparent)

for(k in dbBranches){

indPar = which(treePre[,i]==k)
for(h in (i+1):ncol(treePre)){
tmpChild = unique(treePre[indPar,h])
treeTrim[which(treePre[,h]%in%tmpChild),h] = min(tmpChild)
}
}

MBind = which(treeComment[,i]=='MB')
mbBranches = unique(treePre[MBind,i])

for(k in mbBranches){
indPar = which(treePre[,i]==k)

for(h in (i+1):ncol(treePre)){
tmpChild = unique(treePre[indPar,h])
treeTrim[which(treePre[,h]%in%tmpChild),h] = min(tmpChild)
}
}

cat(paste0('\n\n','resolution ',i,' done','\n\n'))
} #i - each resolution


saveRDS(treeTrim,'./Analysis/tree80_L15_treeTrim.rds')
#treeTrim = readRDS('./Analysis/tree80_L15_treeTrim.rds')

## trimmed tree

labelmat=treeTrim

n = nrow(labelmat)
p = ncol(labelmat)
labelmat = matrix(paste(matrix(rep(colnames(labelmat), each = n),nrow = n), labelmat, sep = "-"), nrow = n)
df = as.data.frame(unique(labelmat), stringsAsFactors = F)
df$pathString = apply(df, 1, function(x) paste(c("all", x),collapse = "/"))
tree.datatree = data.tree::as.Node(df)
treeTrim.phylo = data.tree::as.phylo.Node(tree.datatree)

circ <- ggtree(treeTrim.phylo, layout = "circular") 

pdf("./TestPlots/HSatlas_SC80_L15_TreeTrim_check.pdf",width=100,height=100)
print(circ)
dev.off()


## Hannah's assignments as of 2/22/23

K560_asn = read.csv('./Analysis/tree80_L15_560_assignments.csv')


rownames(K560_asn) <- paste0('K560-',K560_asn$K560)
K560_asn = K560_asn[treeTrim.phylo$tip.label,]


circ <- ggtree(treeTrim.phylo, layout = "circular")
p = gheatmap(circ, K560_asn[19], offset=0, width=.3,colnames_angle=90, colnames_offset_y = .25,font.size =25) 

pdf("./TestPlots/HSatlas_Tree80_L15_mrtreeTrim_CircTest_Hannah_Class_only.pdf",width=100,height=100)
print(p)
dev.off()

## Hannah's colors 


LightPink =   "#ffc2d4"
MedPink =     "#ff9ebb"
DarkPink =    "#ff7aa2"
DarkRed =     "#9e1b44"
LightRed =    "#d64050"
DarkOrange =  "#f36c44"
LightOrange = "#faad60"
Mustard =     "#fecf29"
LightYellow = "#fee08b"
GreenLighest ="#aad576"
GreenMedLight="#78a02d"
GreenMedDark ="#538d22"
GreenDark =   "#245501"
BlueDark =    "#133c55"
BlueMedDark = "#386fa4"
BlueMedLight ="#59a5d8"
BlueLighest = "#84d2f6"
LightPurple = "#b185db"
DarkPurple =  "#7251b5"
Grey =        "#666666"
 
NucleiColors = c("TM"= LightPink, "ARC"= DarkPink, "PVH"= DarkRed, "VMH"= DarkOrange, "DMH"= LightOrange, "LH"= Mustard, "SCN"= GreenLighest, "PO"= GreenMedDark, "AH"= BlueMedDark, "SMN" = BlueLighest, "MN"= LightPurple, "ZI" = DarkPurple,  "Unassigned" = Grey, "Ukn" = Grey,  "Ambiguous" = Grey)
ClassColors = c("Glutaminergic" = LightRed, "GABAergic" = BlueMedDark, "Histaminergic" = GreenLighest, "Unknown" = Grey, "Ambiguous" = Grey)



circ <- ggtree(treeTrim.phylo, layout = "fan",size=2,open.angle=15)
p1 = gheatmap(circ, K560_asn[19], offset=-4, width=.1,colnames_angle=90, colnames_offset_y = 0,font.size =25) +
    scale_fill_manual(breaks=c(names(ClassColors)), 
        values=ClassColors, name="Class")

p2 = gheatmap(p1, K560_asn[17], offset=7, width=.1,colnames_angle=90, colnames_offset_y = 0,font.size =25) +
    scale_fill_manual(breaks=c(names(ClassColors),names(NucleiColors)), 
        values=c(ClassColors,NucleiColors), name="Nuclei")

pdf("./TestPlots/HSatlas_Tree80_L15_mrtreeTrim_CircTest_Hannah_ClassNuclei_colors.pdf",width=100,height=100)
print(p2)
dev.off()



flattree <- ggtree(treeTrim.phylo, size=2)
p1 = gheatmap(flattree, K560_asn[19], offset=-4, width=.1,colnames_angle=90, colnames_offset_y = 0,font.size =25) +
    scale_fill_manual(breaks=c(names(ClassColors)), 
        values=ClassColors, name="Class")

p2 = gheatmap(p1, K560_asn[17], offset=7, width=.1,colnames_angle=90, colnames_offset_y = 0,font.size =25) +
    scale_fill_manual(breaks=c(names(ClassColors),names(NucleiColors)), 
        values=c(ClassColors,NucleiColors), name="Nuclei")

pdf("./TestPlots/HSatlas_Tree80_L15_mrtreeTrim_FlatTest_Hannah_ClassNuclei_colors.pdf",width=100,height=100)
print(p2)
dev.off()


## check agreement among disections 


### depreciated - tried to do all possible combinations 
for(i in 1:length(libraries)){

lib1 = readRDS(paste0('./Analysis/HSatlasNeuro_mt10_Merge_clusters_250_',libraries[1],'.rds'))

Lib1Ind = which(colnames(HSatlasNeuron_mt10)%in%names(lib1$cl))

Lib1Mat = HSatlasNeuron_mt10@assays$integrated@scale.data[,Lib1Ind]

Lib1Mat = Lib1Mat[,names(lib1$cl)]

Lib1Pheno = data.frame(Sample_ID=names(lib1$cl),Celltype = lib1$cl,Study_ID = "Lib1") #745

for(j in 1:length(libraries)){

lib2 = readRDS(paste0('./Analysis/HSatlasNeuro_mt10_Merge_clusters_250_',libraries[2],'.rds'))

Lib2Ind = which(colnames(HSatlasNeuron_mt10)%in%names(lib2$cl))

Lib2Mat = HSatlasNeuron_mt10@assays$integrated@scale.data[,Lib2Ind]

Lib2Mat = Lib2Mat[,names(lib2$cl)]

Lib2Pheno = data.frame(Sample_ID=names(lib2$cl),Celltype = lib2$cl,Study_ID = "Lib2") #745

totPheno = rbind(Lib1Pheno,Lib2Pheno)
rownames(totPheno)= totPheno$Sample_ID
totPheno  = totPheno[,-1]
colnames(totPheno) = c('cell_type','study_id')

totMat = SingleCellExperiment(Matrix(cbind(Lib1Mat,Lib2Mat),sparse=TRUE),colData=totPheno)

celltype_NV_Lib = MetaNeighborUS(var_genes=rownames(Lib1Mat),dat=totMat,study_id=totMat$study_id,cell_type=totMat$cell_type,fast_version=TRUE)



cols=rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))

breaks = seq(0,1,length=101)
pdf("./TestPlots/HiCat_VS_HiCat_MetaNeighbor.pdf",width=100,height=100)
print(gplots::heatmap.2(celltype_NV_Lib,margins=c(25,25),keysize=1,key.xlab="AUROC",key.title="NULL",trace="none",density.info="none",col=cols,breaks=breaks,offsetRow=0.1,offsetCol=0.1,cexRow=2,cexCol=2))
dev.off()

top_hits_IntDat_13pomc=topHits(cell_NV = celltype_NV_13pomc,dat=totMat_13pomc,study_id=totMat_13pomc$study_id,cell_type=totMat_13pomc$cell_type,threshold=0)







## enrichments 

tapply(X=tmp3$Freq,INDEX=tmp3,FUN=sum)


V <- crossprod(table(cbind(tmp$cl,tmp2[,1])))


diag(V) <- 0

## enrichment test for two sample's hicat in seurat clusters? 

## enrichment score for Hannah's marker genes? Can I use funneling of gene enrichment in clusters to inverslely provide a probablity of clusters to region / neuron type? 

## further, 





#test = Seurat:::NNHelper(HSatlasNeuron_mt10)


## 





### can use hicat compare_annotate function to compare suerat clusters to hicat clusters 

