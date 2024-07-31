source('Analysis/util.r')
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)

x <- readRDS("Data/Processed/alphasyn-triplication-initial-integrated.rds")
meta <- data.table::fread("Results/Initial/Initial integration/doublet-finder-metadata.csv") %>% tibble::column_to_rownames("V1")
meta <- meta[meta$DF.classification.upper == 'Singlet',]
x <- subset(x, cells=rownames(meta)) # keep only singlet
DefaultAssay(x) <- 'RNA'
re <- list() # for saving into result list

# Apply dimensionality reduction and clustering within each cluster to find additional contaminating groups of cells
for (cluster in 0:10) { # the levels of SCT.seurat_snn_res.0.1
  if (cluster == 8) {
    dims = 1:10 # due to lack of cells
  } else {
    dims = 1:30
  }
  xs <- subset(x, SCT.seurat_snn_res.0.1 == cluster)
  
  xs@meta.data[['line.id']] <- factor(xs=xs@meta.data[['line.id']])
  xs <- SCTransform(object=xs, assay="RNA", vars.to.regress=c('percent.mt', 'percent.rp', 'nFeature_RNA'), method='glmGamPoi', vst.flavor='v2') # This should run SCTransform on each of the 8 SCT layers we prepared earlier -- if not, then split by line.id and do separately.
  xs <- RunPCA(xs) %>%
    FindNeighbors(dims=dims) %>%
    FindClusters(resolution=0.1, graph.name = c('SCT_nn', 'SCT_snn')) %>%
    RunUMAP(dims=dims, reduction='pca')
  
  saveRDS(xs, paste0("Data/Processed/Integrated-RoundOne/", cluster, ".rds"))
  re[[paste0("g", cluster)]] <- xs
}
saveRDS(re, 'Results/Integration-RoundOne/roundone-integrated-list.rds')

#### Load Results ####
f <- function(x) {
  DefaultAssay(x) <- 'SCT'
  p.list = FeaturePlot(x, cols=c('gray96', 'red'), features=c('STMN2', 'GAP43', 'SNCA', 'DCX', 'VIM', 'HES1', 'SOX2', 'TBR1', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'EOMES', 'TOP2A', 'MKI67', 'APOE', 'S100B', 'SLC1A3', 'MNS1', 'NPHP1', 'BMP4', 'MSX1', 'MYL1', 'MYH3', 'BGN', 'DCN', 'DDIT3'), combine = FALSE, max.cutoff = 'q95')
  return(p.list)
}
g <- function(x) {
  DefaultAssay(x) <- 'SCT'
  #x@active.ident = x$SCT.seurat_snn_res.0.1
  p = VlnPlot(x, features=c('STMN2', 'GAP43', 'SNCA', 'DCX', 'VIM', 'HES1', 'SOX2', 'TBR1', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'EOMES', 'TOP2A', 'MKI67', 'APOE', 'S100B', 'SLC1A3', 'MNS1', 'NPHP1', 'BMP4', 'MSX1', 'MYL1', 'MYH3', 'BGN', 'DCN', 'DDIT3'), ncol = 5)
  return(p)
}
for (i in 1:length(re)) {
  p.list <- f(re[[i]])
  p <- g(re[[i]])
  cowplot::save_plot(plot = wrap_plots(p.list, ncol=5), filename = paste0('Results/Initial/DF round one results/cluster feature plots/cluster  ', i-1, '-feature-plots.png'), base_height=10, base_asp=1.2)
  cowplot::save_plot(plot = p, filename = paste0('Results/Initial/DF round one results/cluster violin plots/cluster ', i-1, '-violin-plots.png'), base_height=10, base_asp=1.2)
}
