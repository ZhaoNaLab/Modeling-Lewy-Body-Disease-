# ********
# 2c-reint-markers.r
# Code for computing marker genes for initial integration results.
source('Analysis/util.r')
library(Seurat)
library(presto)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(MAST)

# Compute marker genes for the initial integration (SCT.seurat_snn_res.0.1 clusters)
#x <- readRDS("Data/Processed/alphasyn-triplication-initial-integrated.rds")
#re <- readRDS('Results/Integration-RoundOne/roundone-integrated-list.rds')
DefaultAssay(x) <- "SCT"

x <- PrepSCTFindMarkers(x)
Idents(x) <- 'SCT.seurat_snn_res.0.1'
all.markers <- FindAllMarkers(
  x,
  slot = 'scale.data',
  only.pos = TRUE,
  min.pct = 0.2
)
write.table(all.markers, 'Results/Initial/DF round one results/Na-YJ-SNCA-roundone-cluster-markers.txt', sep='\t', row.names=T, col.names=T, quote=F)
top_genes <- all.markers %>% 
  group_by(cluster) %>%
  arrange(cluster, p_val_adj, desc(pct.1), pct.2) %>% 
  filter(pct.1 > 0.5) %>% 
  slice_head(n=10) %>%
  ungroup() %>% 
  select(gene) %>%
  unlist

# In addition, compute marker genes for the subclusters obtained within each SCT.seurat_snn_res.0.1 cluster
markers <- list()
for (i in 0:10) {
  re[[paste0("g", i)]] <- PrepSCTFindMarkers(re[[paste0('g', i)]])
  markers[[paste0("g", i)]] <- FindMarkers(
    re[[paste0("g", i)]],
    slot = 'scale.data', 
    min.pct = 0.2,
    only.pos = TRUE
  )
  write.table(markers[[paste0("g", i)]], paste0('Results/Initial/DF round one results/subcluster marker genes/', n, ' subcluster markers.txt'), sep='\t', col.names=T, row.names=FALSE, quote=FALSE)
}

reg = loadRegistry('ExperimentRegistries/reg-integration')
p <- DoHeatmap(x, features=top_genes, group.by='seurat_clusters')
cowplot::save_plot(plot = p,
                   filename='Results/Initial/Initial integration/YJ cluster heatmap.png',
                   base_asp=1.8,
                   base_height=20)

# ---------- Assign initial clusters and celltypes ------------
celltype = list(
  `0` = 'Intermediate cell',
  `1` = 'CN', # Cortical neuron
  `2` = 'CN', # Neural-epithelial cell
  `3` = 'NEC', # non-neuronal
  `4` = 'GPC',
  `5` = 'Neuron (immature)',
  `6` = 'CN',
  `7` = 'AS', # no doubt; 7.4 contaminating
  `8` = 'Inhibitory neuron', # 0 and 3 Inh; 4 is contaminating
  `9` = 'PGC', # no doubt; # maybe non-neuronal or neuronal
  `10` = 'GPC', # glia progenitor cell
)

x@meta.data %<>% mutate(
  celltype = case_when(
    seurat_clusters == 0 ~ 'Intermediate neuron',
    seurat_clusters == 1 ~ 'Cortical neuron', # Cortical neuron
    seurat_clusters == 2 ~ 'Cortical neuron', # Neural-epithelial cell
    seurat_clusters == 3 ~ 'Neural epithelial cells', # non-neuronal
    seurat_clusters == 4 ~ 'Intermediate neuronal',
    seurat_clusters == 5 ~ 'Intermediate neuronal',
    seurat_clusters == 6 ~ 'Eye neuron',
    seurat_clusters == 7 ~ 'Astrocyte', # Divide into AS and CBC; remove 4 and possibly remove 2
    seurat_clusters == 8 ~ 'Inhibitory neuron', # 0 and 3 Inh; 4 is contaminating
    seurat_clusters == 9 ~ 'Proteoglycan cluster', # no doubt; # maybe non-neuronal or neuronal
    seurat_clusters == 10 ~ 'GPC', # glia progenitor cell 
  ),
  neuronal = case_when(
    seurat_clusters == 0 ~ 1,
    seurat_clusters == 1 ~ 1,
    seurat_clusters == 2 ~ 0.5,
    seurat_clusters == 3 ~ 0,
    seurat_clusters == 4 ~ 0.5,
    seurat_clusters == 5 ~ 0.5,
    seurat_clusters == 6 ~ 1,
    seurat_clusters == 7 ~ 0,
    seurat_clusters == 8 ~ 1,
    seurat_clusters == 9 ~ 0.5,
    seurat_clusters == 10 ~ 0,
  )
)
p <- DimPlot(x, group.by='celltype')
cowplot::save_plot(plot=p, filename='Results/RoundOne-Integration/celltype-plot.png', base_asp=1.2, base_height=5)

p <- DimPlot(x, group.by='neuronal')
cowplot::save_plot(plot=p, filename='Results/RoundOne-Integration/neuronal-plot.png', base_asp=1.2, base_height=5)
