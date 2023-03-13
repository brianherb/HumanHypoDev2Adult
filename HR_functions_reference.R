library(Seurat)

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

MM2HSref = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/MM2HS_EG100.csv')
colnames(MM2HSref) = c('MM_ID','MM_Gene','HS_ID','HS_Gene')
HS2MMref = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/HS2MM_EG100.csv')
colnames(HS2MMref) = c('HS_ID','HS_Gene','MM_ID','MM_Gene')


HStoMM <- function(x){
    tmpind = match(as.character(x),HS2MMref$HS_Gene)
return(as.character(HS2MMref$MM_Gene[tmpind]))
}

MMtoHS <- function(x){
    tmpind = match(as.character(x),MM2HSref$MM_Gene)
return(as.character(MM2HSref$HS_Gene[tmpind]))
}

HSTF = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/Human_TF.csv')

NPlist = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20

NPlist = MMtoHS(as.character(NPlist[,1]))


TF_NP = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/TF_NP_list_from_Moffitt.csv') #direct from Moffitt suppmental materials 

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

plotBulkHeatmap <- function(sobj = sobj,genes = genes){
  combined.tmp <- ScaleData(sobj, features = genes)
  x <- intersect(genes, rownames(GetAssayData(combined.tmp, slot = 'data')))
  mat <- AverageExpression(combined.tmp, features = x, slot = 'data')
  mat1 <- mat$RNA
  re <- pheatmap::pheatmap(mat1,color=colorRampPalette(brewer.pal(n = 7, name =
  "YlOrRd"))(100),angle_col = 45,  border = NA,cluster_col=FALSE,cluster_row=FALSE)
  return(re)
}

## slot was scale.data

plotBulkHeatmap2 <- function(sobj = sobj,genes = genes){
  x <- intersect(genes, rownames(GetAssayData(sobj, slot = 'data')))
  mat <- AverageExpression(sobj, features = x, slot = 'data')
  mat1 <- t(scale(t(mat$RNA)))
  re <- pheatmap::pheatmap(mat1, angle_col = 45, border = NA)
  return(re)
}

## Hannah

DEGs = function(Input){
Input.markers <- FindAllMarkers(Input, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
return(Input.markers)}

CheckUMAP = function(SeuFilez){
CheckMeta = as.data.frame(c(rep("RoI", length(CheckInput))))
colnames(CheckMeta) = "Pop"
row.names(CheckMeta) = CheckInput
SeuFilez = AddMetaData(SeuFilez, CheckMeta, "CheckMeta")
Idents(SeuFilez) = "CheckMeta"
DimPlot(SeuFilez, reduction="umap")
}

GenerateMetaData = function(ListMeta){
MetaOutput = as.data.frame(matrix(ncol = 1, nrow =0))  
for(x in seq(1,length(ListMeta),1)){
if(class(ListMeta[[x]]) == "data.frame"){
Temp = as.data.frame(rep(names(ListMeta)[[x]], length(row.names(ListMeta[[x]]))))  
colnames(Temp) = "Pop"
row.names(Temp) = row.names(ListMeta[[x]])
}else{
Temp = as.data.frame(rep(names(ListMeta)[[x]], length(colnames(ListMeta[[x]]))))  
colnames(Temp) = "Pop"
row.names(Temp) = colnames(ListMeta[[x]]) 
}
MetaOutput = rbind(MetaOutput, Temp)  
}
return(MetaOutput)
}

## GeneLists used by Hannah

load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/GeneLists.rda'))

ClusterFunc_All_RNA = function(SeuFile){
Filename = as.character(substitute(SeuFile))

for(y in set.kparam){
for(z in set.dim){
for(v in set.res){
DefaultAssay(SeuFile) = "integrated"
SeuFile <- FindNeighbors(SeuFile, k.param=y, dims=1:z)
SeuFile <- FindClusters(SeuFile, resolution = v)
DefaultAssay(SeuFile) = "RNA"

pdf(paste("ALLFETAL_", Filename, "_res", v, "_k", y, "_dim", z, "_umapSML.pdf", sep=""), width=12, height=10)
dimplot = DimPlot(SeuFile, reduction="umap", label=T)
print(dimplot)
dev.off()


pdf(paste("ALLFETAL_", Filename, "_res", v, "_k", y, "_dim", z, "_umapLGE.pdf", sep=""), width=25, height=25)
dimplot = DimPlot(SeuFile, reduction="umap", label=T)
print(dimplot)
dev.off()

###VLNS
AllVlns = list()

for(d in names(GeneLists)){
Genes = GeneLists[[d]]
StorePlots = list()
  for(x in Genes[1]){
            plotA <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
            plotA <- plotA + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.y = element_blank(), legend.position = "none", plot.title = element_text(size=12))+ labs(title = d, subtitle = Genes[1])
            StorePlots[[x]] = plotA 
            }
  for(x in Genes[2:length(Genes)]){
            plotB <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
            plotB <- plotB + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), 
                      axis.ticks.y = element_blank(), legend.position = "none", axis.text.y = element_blank(), plot.title = element_text(size=12))
           StorePlots[[x]] = plotB
           }
AllVlns[[d]] <- ggarrange(plotlist = StorePlots, widths=c(1.4, rep(1, length(Genes)-1)), ncol = 40,  nrow = 1)  
}
pdf(paste("ALLFETAL_", Filename, "_res", v, "_k", y, "_dim", z, "_AllMultiVlns.pdf", sep=""), width=40, height=length(unique(SeuFile@active.ident)))
print(AllVlns)
dev.off()

#Cell No
CellNo = as.data.frame(table(SeuFile@meta.data$seurat_clusters))
write.csv(CellNo, paste("ALLFETAL_", Filename, "_counts_k", y, "_dim", z, "_res", v, ".csv", sep=""), row.names = F)
}}

#Feature Plots  
FPList = list()  

for(d in names(GeneLists)){
Genes = GeneLists[[d]]

FPSinglePage = list()
FPSinglePage[[1]] = FeaturePlot(SeuFile, Genes[1], reduction="umap") + labs(title=paste(d, Genes[1]))
for(p in seq(2, length(Genes), 1)){
FPSinglePage[[p]] = FeaturePlot(SeuFile, Genes[p], reduction="umap") 
}
FPList[[d]] = ggarrange(plotlist = FPSinglePage, ncol=5, nrow=8)
}

pdf(paste("ALLFETAL_", Filename, "_dim", z, "_AllFPs.pdf", sep=""), width = 25, height = 45)
print(FPList)
dev.off()
}}


PlotFunc = function(SeuFile){
Filename = as.character(substitute(SeuFile))

pdf(paste("ALLFETAL_", Filename, "_umapSML.pdf", sep=""), width=12, height=10)
dimplot = DimPlot(SeuFile, reduction="umap", label=T)
print(dimplot)
dev.off()


pdf(paste("ALLFETAL_", Filename,  "_umapLGE.pdf", sep=""), width=25, height=25)
dimplot = DimPlot(SeuFile, reduction="umap", label=T)
print(dimplot)
dev.off()

###VLNS
AllVlns = list()

for(d in names(GeneLists)){
Genes = GeneLists[[d]]
StorePlots = list()
  for(x in Genes[1]){
            plotA <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
            plotA <- plotA + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.y = element_blank(), legend.position = "none", plot.title = element_text(size=12))+ labs(title = d, subtitle = Genes[1])
            StorePlots[[x]] = plotA 
            }
  for(x in Genes[2:length(Genes)]){
            plotB <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
            plotB <- plotB + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), 
                      axis.ticks.y = element_blank(), legend.position = "none", axis.text.y = element_blank(), plot.title = element_text(size=12))
           StorePlots[[x]] = plotB
           }
AllVlns[[d]] <- ggarrange(plotlist = StorePlots, widths=c(1.4, rep(1, length(Genes)-1)), ncol = 40,  nrow = 1)  
}
pdf(paste("ALLFETAL_", Filename, "_AllMultiVlns.pdf", sep=""), width=40, height=length(unique(SeuFile@active.ident)))
print(AllVlns)
dev.off()

#Cell No
CellNo = as.data.frame(table(SeuFile@meta.data$seurat_clusters))
write.csv(CellNo, paste("ALLFETAL_", Filename, "_counts.csv", sep=""), row.names = F)

#Feature Plots  
FPList = list()  

for(d in names(GeneLists)){
Genes = GeneLists[[d]]

FPSinglePage = list()
FPSinglePage[[1]] = FeaturePlot(SeuFile, Genes[1], reduction="umap") + labs(title=paste(d, Genes[1]))
for(p in seq(2, length(Genes), 1)){
FPSinglePage[[p]] = FeaturePlot(SeuFile, Genes[p], reduction="umap") 
}
FPList[[d]] = ggarrange(plotlist = FPSinglePage, ncol=5, nrow=8)
}

pdf(paste("ALLFETAL_", Filename, "_AllFPs.pdf", sep=""), width = 25, height = 45)
print(FPList)
dev.off()
}
