#r4.1

library(dplyr)
library(chooseR)
library(matrixStats)
library(tictoc)
library(boot)
library(ggplot2)

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}


boot_median <- function(x, interval = 0.95, R = 25000, type = "bca") {
  # Define median to take data and indices for use with boot::
  med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }

  # Calculate intervals
  boot_data <- boot::boot(data = x, statistic = med, R = R)
  boot_ci <- boot::boot.ci(boot_data, conf = interval, type = type)

  # Extract desired statistics
  ci <- list(
    low_med = boot_ci$bca[4],
    med = boot_ci$t0,
    high_med = boot_ci$bca[5]
  )
  return(ci)
}


setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/')

#HSatlasNeuron_mt10 = readRDS('./SeuratObj/HSatlasNeuro_mt10_integrated.rds')

#totCells = colnames(HSatlasNeuron_mt10)

#saveRDS(totCells,'./Analysis/chooseR/HSatlasCells.rds')

totCells = readRDS('./Analysis/chooseR/HSatlasCells.rds')

#rndNames = gsub('_Cells.rds','',list.files(path='./Analysis/chooseR/',pattern='_Cells.rds'))


## original tree trim 

treeTrim = readRDS('./Analysis/tree80_L15_treeTrim.rds')

#rndNames = gsub('_HSneuMRTree15trim.rds','',list.files(path='./Analysis/chooseR/',pattern='_HSneuMRTree15trim.rds'))

## first one was OS32RV, might have got messed up

#rndNames = rndNames[-grep('OS32RV',rndNames)]

#rndNames = rndNames[1:100] ## did a few more iterations, just ignore

#saveRDS(rndNames,'./Analysis/chooseR/rndNames_used.rds')

rndNames = readRDS('./Analysis/chooseR/rndNames_used.rds')


Trees = vector(mode='list',length=length(rndNames))

names(Trees) = rndNames

for(i in rndNames){

Trees[[i]] = readRDS(paste0('./Analysis/chooseR/',i,'_HSneuMRTree15trim.rds'))

}

#ResScore = vector(mode='list',length=ncol(Trees[[1]]))

#names(ResScore) = colnames(treeTrim)




for(k in 1:nrow(treeTrim)){ ## per clust res

origClust = sort(unique(treeTrim[totCells,k]))

origClustInd = treeTrim[totCells,k]

## construct matrix of clusters from iterations

tmpClust = matrix(NA,nrow=length(totCells),ncol=length(names(Trees)))

rownames(tmpClust) = totCells
colnames(tmpClust) = names(Trees)

for(j in names(Trees)){

tmpClust[rownames(Trees[[j]]),j] = Trees[[j]][,k]

}

## based on original clusters - records average distance from cell to all other cells within each cluster
tmpClustMat = matrix(nrow=length(totCells),ncol=length(origClust),0)

rownames(tmpClustMat) = totCells
colnames(tmpClustMat) = origClust

tic()

for(m in 1:length(totCells)){

tmpIdent = t(matrix(nrow=length(names(Trees)),ncol=length(totCells),tmpClust[m,]))

tmpBoo = tmpIdent==tmpClust

tmpNAcount =rowSums(is.na(tmpBoo))

tmpEvid = ncol(tmpBoo)-tmpNAcount

tmpBi = matrix(nrow=length(totCells),ncol=length(names(Trees)),0)

tmpBi[tmpBoo] = 1

tmpSums = rowSums(tmpBi)

tmpFreq = tmpSums/tmpEvid

## average per original cluster
tmpDist = 1-tmpFreq

tmpDistClust= tapply(X=tmpDist,IND=origClustInd, FUN=mean) 

tmpClustMat[m,] = tmpDistClust

if(m%%100==0) cat(paste0(m,', '))

}

toc()

## reapply origClustInd to tmpClustMat to calculate Silhouette scores - one score per pair - presumably compare against closest cluster (keep lowest score)

## can use 

tmpClustMatInd =matrix(nrow=length(totCells),ncol=length(origClust),FALSE)
colnames(tmpClustMatInd) = origClust
rownames(tmpClustMatInd) = totCells


for(n in 1:length(totCells)){

tmpClustMatInd[n,as.character(origClustInd[n])]=TRUE
if(n%%100==0) cat(paste0(n,', '))

}

WithinVal = rep(0,length(totCells))

for(f in 1:length(totCells)){
WithinVal[f] = tmpClustMat[f,][tmpClustMatInd[f,]]
if(f%%100==0) cat(paste0(f,', '))

}

ClosestVal  = rep(0,length(totCells))

for(f in 1:length(totCells)){
ClosestVal[f] = min(tmpClustMat[f,][-which(tmpClustMatInd[f,])])
if(f%%100==0) cat(paste0(f,', '))

}

Silhouettes = (ClosestVal-WithinVal)/rowMaxs(cbind(ClosestVal,WithinVal))

saveRDS(Silhouettes,file=paste0('./Analysis/chooseR/Indv_SilhouetteScores_',colnames(treeTrim)[k],'.rds'))

SilMean = tapply(X=Silhouettes,IND=origClustInd,FUN=mean)

saveRDS(SilMean,file=paste0('./Analysis/chooseR/Avg_SilhouetteScores_',colnames(treeTrim)[k],'.rds'))

saveRDS(tmpClustMat,file=paste0('./Analysis/chooseR/Avg_Distances_',colnames(treeTrim)[k],'.rds'))


} ## end per resolution


### plotting 


AvgSilPerClust = vector(mode='list',length=ncol(treeTrim))

names(AvgSilPerClust) = colnames(treeTrim)

for(k in names(AvgSilPerClust)){

AvgSilPerClust[[k]] = readRDS(paste0('./Analysis/chooseR/Avg_SilhouetteScores_',k,'.rds'))

}

scores = do.call(c,AvgSilPerClust)

scoreDF = data.frame(SilhouetteScore = scores,Resolution=flexsplit(names(scores),'.')[,1])

scoreDF$Resolution <- factor(scoreDF$Resolution, levels = unique(scoreDF$Resolution))

Resolution = rownames(bootTot)




## based on plotting code from chooseR

threshold <- max(bootTot$low_med)
choice <- as.character(
  bootTot %>%
  dplyr::filter(med >= threshold) %>%
  dplyr::arrange(n_clusters) %>%
  tail(n = 1) %>%
  dplyr::pull(Resolution)
)


choice = 'K116'


ggplot(bootTot, aes(factor(Resolution, levels=unique(Resolution)), med)) +
  geom_crossbar(aes(ymin = low_med, ymax = high_med),fill = "grey", size = 0.25) +  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +   geom_jitter(
    data = scoreDF,
    aes(factor(Resolution), SilhouetteScore),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "SilhouetteScore",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

  ggsave(
  filename = paste0("./TestPlots/silhouette_distribution_plot.png"),
  dpi = 300,
  height = 3.5,
  width = 5,
  units = "in",
  bg='white'
)



## 








### K116 (k=6) and K169 (k=7) don't look like they finished - check avg dist 

K116dist = readRDS('./Analysis/chooseR/Avg_Distances_K116.rds') # zero at tail 

K169dist = readRDS('./Analysis/chooseR/Avg_Distances_K169.rds') # zero at tail 













### Junque 

resolutions = colnames(treeTrim)

library(purrr)

results_path = "/autofs/burnsfs/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite/Analysis/chooseR/"

scores <- purrr::map(
  paste0(results_path, "Avg_SilhouetteScores_", resolutions, ".rds"),
  readRDS
)

scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()
meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)









## testing

## just need to construct matrix per resolution

tree1 = matrix(nrow=500,ncol=100,sample(1:12,50000,replace=TRUE))

test1 = t(matrix(nrow=100,ncol=500,tree1[1,]))

boo1 = tree1==test1

res1 = matrix(nrow=500,ncol=100,0)

res1[boo1] = 1

freq1 = rowSums(res1)/100





samp100 = sample(1:20,100,replace=TRUE)


tmp = tapply(X=samp100,IND=rep(1:10,each=10),FUN=mean)


list_sum <- Reduce("+", my_list)  


