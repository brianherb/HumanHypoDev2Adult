source("/local/projects-t3/idea/bherb/software/chooseR/pipeline.R")
library(Seurat)
library(ggplot2)
`%>%` <- magrittr::`%>%`

# Load Seurat object
# Make sure this has already been normalised and that PCA has been performed
# The data included here is the Ding, et al., Nature Biotechnology, 2020
# human PBMC Smart-Seq data from our manuscript
#obj <- readRDS("data/ding_smartseq_pbmc_preprocessed.rds")

#obj = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HSatlasNeuro_mt10_integrated.rds'))

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite')

obj = readRDS('./SeuratObj/HSatlasNeuro_mt10_integrated.rds')

# Define the number of PCs to use, and which assay and reduction to use.
# We recommend testing a broad range of resolutions
# For more on picking the correct number of PCs, see:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
npcs <- 30
resolutions <- c(0.01,0.1,0.2,0.5,1, 5, 10, 15, 20, 30, 40, 50,60,70,80)
assay <- "integrated"
reduction <- "pca"
results_path <- "./Analysis"

#DefaultAssay(obj) <- "integrated"
#obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# Run pipeline
for (res in resolutions) {
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")

  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  obj <- find_clusters(obj,reduction = reduction,assay = assay,resolution = res)
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]

  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(obj, n = 100, size = 0.8, npcs = npcs, res = res,reduction = reduction, assay = assay)

  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  # This is the more time efficient vectorisation
  # However, it exhausts vector memory for (nearly) all datasets
  # matches <- purrr::map(columns, find_matches, df = results)
  # matches <- purrr::reduce(matches, `+`)
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }

  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )

  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "...")) ## here - vector too large

### break - testing 
treeTrim = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/tree80_L15_treeTrim.rds'))

## write out for python

write.csv(Embeddings(object = obj[[reduction]])[, 1:npcs],file='./Analysis/HSatlasNeuro_Top30pcs.csv',row.names = FALSE)

write.csv(treeTrim,file='./Analysis/tree80_L15_treeTrim.csv',row.names = FALSE)



test = as.matrix(mtchs)

library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = obj[[reduction]])[, 1:npcs])
#clusters <- dataset$celltype
clusters <- treeTrim[,1]
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)


### end break 


  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))

  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
}

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))

# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
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

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
  dplyr::filter(med >= threshold) %>%
  dplyr::arrange(n_clusters) %>%
  tail(n = 1) %>%
  dplyr::pull(res)
)

# And plot!
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
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
  filename = paste0(results_path, "silhouette_distribution_plot.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Finally, a dot plot of silhouette scores to help identify less robust clusters
# The initial pipe is to order the clusters by silhouette score
scores %>%
  dplyr::filter(res == choice) %>%
  dplyr::arrange(dplyr::desc(avg_sil)) %>%
  dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
  ggplot(aes(factor(cluster), avg_sil)) +
    geom_point() +
    scale_x_discrete("Cluster") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_grid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )

ggsave(
  filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)



### other cluster eval - bootstrapping 

library(scran)

sce <- as.SingleCellExperiment(obj)
SingleCellExperiment::reducedDim(sce, 'PCA_sub') <- SingleCellExperiment::reducedDim(sce, 'PCA')[,1:30, drop = FALSE]

sce@colData = cbind(sce@colData,treeTrim)


ass_prob <- bluster::bootstrapStability(sce, FUN = function(x) {
    g <- buildSNNGraph(x, use.dimred = 'PCA_sub')
    igraph::cluster_walktrap(g)$membership
  },clusters = sce$K6)


saveRDS(ass_prob,file='./Analysis/HSatlas_neu_bootstrap_prob.rds')

ass_prob = readRDS('./Analysis/HSatlas_neu_bootstrap_prob.rds')

conda activate /autofs/burnsfs/projects-t3/idea/amentlab_software/conda_envs/r_4.2.3

library(tidyverse) #not working in r4.1
library(patchwork)
library(Seurat)
library(scran) #not avaliable for 4.2.3
library(patchwork)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)




p <- ass_prob %>%
  as_tibble() %>%
  rownames_to_column(var = 'cluster_1') %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = 'cluster_2',
    values_to = 'probability'
  ) %>%
  mutate(
    cluster_1 = as.character(as.numeric(cluster_1) - 1),
    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
  ) %>%
  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
  scale_x_discrete(name = 'Cluster', position = 'top') +
  scale_y_discrete(name = 'Cluster') +
  scale_fill_gradient(
    name = 'Probability', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
    limits = c(0,1),
    guide = guide_colorbar(
      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
      title.theme = element_text(hjust = 1, angle = 90),
      barwidth = 0.75, barheight = 10
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = 'right',
    panel.grid.major = element_blank()
  )

ggsave('./TestPlots/cluster_stability.png', p, height = 6, width = 7)




