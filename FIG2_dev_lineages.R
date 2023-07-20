library(Seurat)
library(BiocGenerics)
library(S4Vectors)

library(plotly)
library(htmlwidgets)
library("rmarkdown")
library(Seurat)
library(monocle3)
library(scales)
library(igraph)
library(ggraph)
#library(networkdata) #
library(graphlayouts)
library(qgraph)
library(AUCell)
#library(clustree)
#library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')




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

## Load monocle lineages and nodes - 

load('./Analysis/Fig2_Trajectories_GLOTMP_17JAN_20_20_3_k100.RData') ## loaded on a 200Gb machine - needed? Got error, but seemed to work

## object "Pull" has monocle3 node assignments 

## human neurons

EdKaZhouHypoNeurons = readRDS("./SeuratObj/EdKaZhouHypoNeurons_mt10_integrated.rds")


EdKaZhouHypoNeurons@meta.data$Vertex = Pull[colnames(EdKaZhouHypoNeurons),'Vertex']

## mouse, human neurons 
MMHSneu.integrated = readRDS('./SeuratObj/MusRefEdKaZhouHypoNeurons_mt10_integrated.rds')


lin=monocle3:::get_principal_path(M3Seu,reduction_method = "UMAP",starting_cell = 'Y_376',end_cells = 'Y_3')$nodes



slot(M3Seu, 'preprocess_aux', check=FALSE) <- SimpleList(c())  ##fix 


## get vertex labels for reference

pdf('./TestPlots/Hannah_July2023_neuron_lineage.pdf',width=12,height=8)
plot(net.NODE,layout=lrt.NODE,vertex.label.cex=0.5)
dev.off()



ind = which(EdKaZhouHypoNeurons@meta.data$Vertex%in%nodesNODE$id)

tmpGene = EdKaZhouHypoNeurons@assays$RNA@counts['ASCL1',ind]

tmpGeneVertex = tapply(X=tmpGene,INDEX=EdKaZhouHypoNeurons@meta.data$Vertex[ind],FUN=mean)


tmpGeneVertex2 = round(scales::rescale(tmpGeneVertex, to = c(1, 100)),0)

fun_color_range <- colorRampPalette(c("#1b98e0", "red"))
my_colors <- fun_color_range(100)

tmpGeneColors = my_colors[tmpGeneVertex2]

names(tmpGeneColors) = names(tmpGeneVertex2)

V(net.NODE)$color <- tmpGeneColors[as.vector(net.NODE$id)]


pdf('./TestPlots/EdKaZhouHypoNeurons_lineage_ASCL1exp.pdf',width=12,height=8)
plot(net.NODE,layout=lrt.NODE,vertex.label.cex=0.5)
dev.off()


DefaultAssay(EdKaZhouHypoNeurons) = 'RNA'

pdf(file="./TestPlots/EdKaZhouHypoNeurons_ASCL1.pdf",width=12,height=8)
print(DimPlot(EdKaZhouHypoNeurons, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
#print(DimPlot(KaHypoNeuro_mt10, reduction = "umap", group.by = 'TopRegion', label = TRUE, repel = TRUE))

for(i in c('GAD2','SLC17A6','ASCL1','HDC')){
print(FeaturePlot(EdKaZhouHypoNeurons, features = i, pt.size = 0.2))

}

dev.off()

## other plots for sup 

nuclei = read.csv('./Analysis/AllHumanAnnots_ByCell_12JUL23.csv')

treeTrim=readRDS('./Analysis/tree80_L15_treeTrim.rds')


EdKaZhouHypoNeurons@meta.data$K116 = NA

EdKaZhouHypoNeurons@meta.data[rownames(treeTrim),'K116']=treeTrim[,'K116']


matInd = match(EdKaZhouHypoNeurons@meta.data$K116,nuclei$K116)

EdKaZhouHypoNeurons@meta.data$K116_nuclei = nuclei[matInd,'K_116Nuclei']

EdKaZhouHypoNeurons@meta.data$K116_class = nuclei[matInd,'K116_Class']


EdKaZhouHypoNeurons$Stage = "Trimester1"

EdKaZhouHypoNeurons@meta.data[EdKaZhouHypoNeurons$sample%in%c('GW15_A','GW15_M','GW15_P','GW16_hypo','GW18_A','GW18_hypo','GW18_Lane1','GW18_Lane2','GW18_Lane3','GW18_M','GW18_P','GW19_hypo','GW20_34_hypo','GW20_A','GW20_M','GW20_P','GW22T_hypo1','GW25_3V_hypo'),'Stage'] = 'Trimester2'

EdKaZhouHypoNeurons@meta.data[grep('10X',EdKaZhouHypoNeurons$sample),'Stage'] = 'Adult'


pdf(file="./TestPlots/EdKaZhouHypoNeurons_K116_nuclei_class.pdf",width=12,height=8)
print(DimPlot(EdKaZhouHypoNeurons, reduction = "umap", group.by = 'K116_nuclei', label = TRUE, repel = TRUE))
print(DimPlot(EdKaZhouHypoNeurons, reduction = "umap", group.by = 'K116_class', label = TRUE, repel = TRUE))
print(DimPlot(EdKaZhouHypoNeurons, reduction = "umap", group.by = 'Stage', label = TRUE, repel = TRUE))
dev.off()
### Hannah did prepare rtti287m, which is actually rooted in Y_376

## POMC / KISS1 branch is Y_534 - Y_359


for( i in 1:nrow(rtti287m)){

FROM=rtti287m$root[i]
TO=rtti287m$leaf[i]

gc()

#if(file.exists(paste("./Analysis/EdKaHypoNeuron_mt10_Monocle3_3D_LineageGenes_From_",FROM,"_To_",TO,".csv",sep=''))) next

#lins = all_simple_paths(M3Seu287@principal_graph[['UMAP']],FROM,TO)
#lengs = unlist(lapply(lins,FUN=length))
#shortInd = which(lengs==min(lengs))

## Hijack for POMC/KISS1
# FROM = 'Y_534'
# TO = 'Y_359'


## new for revision July 2023 - 

## ARC: NPs - POMC, KISS1 and AGRP. Vertexes: Y_131,  Y_472,  Y_134,  Y_453,  Y_389,  Y_377,  Y_380,  Y_468,  Y_101,  Y_170,  Y_8,  Y_337,  Y_167,  Y_336,  Y_378,  Y_300,  Y_356,  Y_358,  Y_379, Y_357, Y_463,  Y_301, Y_462,  Y_359
# FROM = 'Y_131'
# TO = 'Y_359'

NPs = c('POMC', 'KISS1')
region='ARC'


#PVH: NPs – OXT, CRH, AVP. Vertexes: Y_82, Y_499, Y_76, Y_518, Y_2, Y_495, Y_399, Y_424, Y_398, Y_89, Y_185, Y_409, Y_48, Y_533, Y_136, Y_124, Y_9, Y_316, Y_448, Y_407, Y_540, Y_199, Y_251, Y_220, Y_12, Y_291, Y_387, Y_425, Y_465, Y_427, Y_75, Y_541, Y_239, Y_237, Y_259, Y_64, Y_228

## termini - Y_64(use in paper), Y_228, Y_251
## branched Y_82 and Y_89

# FROM = 'Y_82'
# TO = 'Y_64'

NPs = c('OXT', 'CRH', 'AVP')
region='PVH'

#SCN: NPs - VIP. Vertexes: Y_87, Y_544, Y_559, Y_548, Y_253, Y_216, Y_527, Y_56, Y_54, Y_530, Y_397, Y_413, Y_434, Y_20, Y_116, Y_126, Y_125, Y_189, Y_528, Y_52, Y_439, Y_137, Y_417, Y_286, Y_31, Y_526, Y_467, Y_531, Y_4, Y_543, Y_451, Y_410, Y_469, Y_546, Y_29, Y_265, Y_209, Y_53, Y_11, Y_3

## termini - Y_3 (use in paper), Y_125 (single node branch) , Y_286

# FROM = 'Y_87'
# TO = 'Y_3'

NPs = c('VIP')
region='SCN'


lin=monocle3:::get_principal_path(M3Seu,reduction_method = "UMAP",starting_cell = FROM,end_cells = TO)$nodes
lin = as.numeric(gsub('Y_','',lin))


TMPcds = M3Seu[,which(M3Seu@principal_graph_aux[['UMAP']]['pr_graph_cell_proj_closest_vertex'][[1]][,1]%in%lin)]

TMPcds <- order_cells(TMPcds, root_pr_nodes=FROM)

## more memory
p=plot_cells_3d(TMPcds,color_cells_by = "pseudotime")

htmlwidgets::saveWidget(p, paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_UMAP_pseudotime_From_",FROM,"_To_",TO,".html",sep=''))

## DEG
gene_fits <- fit_models(TMPcds, model_formula_str = "~pseudotime")
fit_coefs <- as.data.frame(coefficient_table(gene_fits))

write.csv(fit_coefs[,-c(4,5)],file=paste("./Analysis/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_LineageGenes_From_",FROM,"_To_",TO,".csv",sep=''))

fit_coefs_pseu = fit_coefs[which(fit_coefs$term=='pseudotime'),]
fit_coefs_pseu = fit_coefs_pseu[order(fit_coefs_pseu$q_value,decreasing=FALSE),]
Top_genes_up <- fit_coefs_pseu$gene_short_name[which(fit_coefs_pseu$estimate>0)]
MTindUp = grep('MT-',Top_genes_up)

if(length(MTindUp)>0)  Top_genes_up=Top_genes_up[-MTindUp]
Top_genes_up=Top_genes_up[1:25]
Top_lineage_cds_up <- TMPcds[rowData(TMPcds)$gene_short_name %in% Top_genes_up,]
Top_genes_dn <- fit_coefs_pseu$gene_short_name[which(fit_coefs_pseu$estimate<0)]
MTindDn = grep('MT-',Top_genes_dn)

if(length(MTindDn)>0)  Top_genes_dn=Top_genes_dn[-MTindDn]
Top_genes_dn=Top_genes_dn[1:25]
Top_lineage_cds_dn <- TMPcds[rowData(TMPcds)$gene_short_name %in% Top_genes_dn,]

pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_Top25Genes_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(Top_lineage_cds_up, color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(Top_lineage_cds_dn, color_cells_by="sample", min_expr=0.5))
dev.off()



pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_",region,"_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% NPs,], color_cells_by="sample", min_expr=0.5))
dev.off()

## POMC - KISS1 transition



fit_coefs_pseu_TF = fit_coefs_pseu[which(fit_coefs_pseu$gene_short_name%in%unique(c(HSTF$Name,TF_NP$Transcription_factors))),]


fit_coefs_pseu_TF = fit_coefs_pseu_TF[which(fit_coefs_pseu_TF$q_value<=0.05),] #309

fit_coefs_pseu_TF = fit_coefs_pseu_TF[which(fit_coefs_pseu_TF$num_cells_expressed>500),]



pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_TopTFs_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[1:5],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[6:10],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[11:15],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[16:20],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[21:25],], color_cells_by="sample", min_expr=0.5))
dev.off()

pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_TopTFsNext25_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[26:30],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[31:35],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[36:40],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[41:45],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_TF$gene_id[46:50],], color_cells_by="sample", min_expr=0.5))
dev.off()


### just check for other NP's



fit_coefs_pseu_NP = fit_coefs_pseu[which(fit_coefs_pseu$gene_short_name%in%NPlist2),]


fit_coefs_pseu_NP = fit_coefs_pseu_NP[which(fit_coefs_pseu_NP$q_value<=0.05),] #309

fit_coefs_pseu_NP = fit_coefs_pseu_NP[which(fit_coefs_pseu_NP$num_cells_expressed>100),]



pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_TopNPs_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[1:5],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[6:10],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[11:15],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[16:20],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[21:25],], color_cells_by="sample", min_expr=0.5))
dev.off()

pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_TopNPsNext25_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=12)
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[26:30],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[31:35],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[36:40],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[41:45],], color_cells_by="sample", min_expr=0.5))
print(plot_genes_in_pseudotime(TMPcds[rowData(TMPcds)$gene_short_name %in% fit_coefs_pseu_NP$gene_id[46:50],], color_cells_by="sample", min_expr=0.5))
dev.off()







## try heatmap and order by exp, try to demonstrate ones aligned with POMC (first row) or KISS1 (last row )
#genes     <- read.delim("data/TF_heatmap.txt",header=T,sep="\t",stringsAsFactors=F)[,2]


## genes=c('POMC','TSHZ2','ZFPM2','ISL1','RORB','BCL11A','TCF7L2','INSM1','SOX3','DACH1','SIX3','GAS7','HDAC9','PLAGL1','ELL2','EEA1','TCF4','TCF12','BACH2','NR5A2','L3MBTL4','ESR1','AR','KISS1') ##ARC

# genes=c('ZIC1','ZIC4','POU2F2','ZIC2','ZFHX3','ZFHX4','FOXP2','AVP','CREB3L2','ZNF704','FOSL2','POU3F2','CRH','HIPK2','ZNF385B','OXT','ELL2','ZNF804A','ZBTB20','POU2F1','SH3D19','STAT5B','NPAS2','SIM1') ##PVH 82-64

# genes=c('HDAC2','ZIC1','ZIC2','AVP','PTMA','SOX4','TEAD1','YBX1','CRH','NACA','BTF3','OXT','ANK1','CUX2','PPARGC1A','PPP1R12B','TSHZ2') ##PVH 82-228

# genes=c('ANK1','NRF1','PRMT3','BACH2','ZNF536','GLIS3','MYT1L','GLIS1','THRB','CUX2','SIX3','MEIS2','POU6F2','SIX6','LMO4','VIP','RORB','NFIB','LHX1','NFIA','NFIX','ESR1','SATB2','INSM1') ## SCN 87-3



pt.matrix <- as.matrix(assays(TMPcds)$counts[match(genes,rowData(TMPcds)[,1]),order(pseudotime(TMPcds))])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

## order heatmap by max z score placement along pseudotime

maxZ = apply(pt.matrix,1,FUN=function(x){mean(which(x==max(x)))})
pt.matrix=pt.matrix[order(maxZ),]

## close, but going manual 

## pt.matrix=pt.matrix[c('SOX3','INSM1','BCL11A','GAS7','TSHZ2','TCF7L2','ISL1','ZFPM2','DACH1','SIX3','HDAC9','RORB','POMC','NR5A2','EEA1','BACH2','ESR1','KISS1','PLAGL1','ELL2','TCF4','TCF12','L3MBTL4','AR'),] ## original POMC + KISS1 - 

## pt.matrix=pt.matrix[c('SOX3','INSM1','GAS7','TSHZ2','ISL1','DACH1','BCL11A','TCF7L2','ZFPM2','SIX3','HDAC9','RORB','POMC','NR5A2','EEA1','BACH2','ESR1','PLAGL1','ELL2','TCF4','TCF12','L3MBTL4','AR','KISS1'),] ##new ARC without AGRP 

## pt.matrix=pt.matrix[c('ZIC1','ZIC2','ZIC4','POU2F2','ZFHX3','ZFHX4','CRH','AVP','POU3F2','ZBTB20','FOXP2','CREB3L2','ZNF704','FOSL2','HIPK2','ZNF385B','ELL2','ZNF804A','POU2F1','SH3D19','STAT5B','NPAS2','SIM1'),] ##PVH 82-64





## pt.matrix=pt.matrix[c('MEIS2','NRF1','ANK1','PRMT3','BACH2','MYT1L','GLIS1','THRB','GLIS3','CUX2','ZNF536','POU6F2','VIP','RORB','SIX3','SIX6','LHX1','INSM1','NFIX','NFIB','NFIA','LMO4','ESR1','SATB2'),] ## SCN 87-3


ht <- Heatmap( pt.matrix,name = "z-score",col= colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),show_row_names = TRUE,show_column_names=FALSE,row_names_gp= gpar(fontsize = 10),km = 1, row_title_rot = 0,cluster_rows = FALSE,cluster_row_slices= FALSE,cluster_columns = FALSE) #

pdf(file=paste("./TestPlots/EdKaZhouHypoNeuron_mt10_Monocle3_3D_Root376_",region,"_Heatmap_From_",FROM,"_To_",TO,".pdf",sep=''),width=12,height=8)
print(ht)
dev.off()











TBX3 and SIX3, SIX6, PRDM12 and HMGN3  ISL1, DRAP1, PIN1, NHLH2, CXXC4, ZFHX4, TBPL1 and TEAD1



## when do certain NP's first be expressed? 

## load KaZhouAll_mt10_integrated

KaZhouAll.integrated = readRDS('./SeuratObj/KaZhouAll_mt10_integrated.rds')

table(KaZhouAll.integrated$sample[which(KaZhouAll.integrated@assays$RNA@counts['GHRH',]>0)])

POMC, GHRH, GAL,  - CS22 (GW7)
KISS1 - more in GW10? 
AGRP - GW15
NPY - Zhou shows some expression in GW7, but almost nothing in GW8 and KA CS22
TRH - GW7 (strong early exp)
OXT - some expression early, most tri2
AVP - strong exp in GW8
CRH - truley late GW12

 - CS22 (GW7)


