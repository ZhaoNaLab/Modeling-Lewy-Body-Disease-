# ********
# 1a-integration.r
# Code for preliminary integration and clustering of the data by subclone.

source("Analysis/util.r")
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)

x = readRDS('Data/Processed/alphasyn-triplication-initial-filtered-seurat-wo13.rds')
x$percent.mt <- PercentageFeatureSet(x, '^MT-')
x$percent.rp <- PercentageFeatureSet(x, '^RP(S|L)')
x$orig.ident <- factor(x$orig.ident, paste0('BU-SNCA-', 1:24))
x$Line <- stringr::str_split_fixed(x$line.id, ' ', 2)[,1]

library(future.apply)
x@meta.data[['line.id']] <- factor(x=x@meta.data[['line.id']])
x <- SplitObject(object=x, split.by='line.id')
x <- lapply(x, function(y)  SCTransform(object=y, assay="RNA", vars.to.regress=c('percent.mt', 'percent.rp'), method='glmGamPoi', vst.flavor='v2'))
for (i in 1:length(x)) {
  levels(x[[i]][['SCT']]) <- names(x)[[i]]
}
features <- SelectIntegrationFeatures(object.list=x, nfeatures=3000)
x <- PrepSCTIntegration(object.list=x, anchor.features=features)
anchors <- FindIntegrationAnchors(object.list=x, anchor.features=features, normalization.method="SCT", reference=1, dims=1:50, scale=FALSE)
x <- IntegrateData(anchorset=anchors, normalization.method='SCT', dims=1:50)
x <- RunPCA(x) %>%
  FindNeighbors(dims=1:50) %>%
  FindClusters(resolution=0.1, graph.name = c('SCT.seurat_nn', 'SCT.seurat_snn')) %>%
  RunUMAP(dims=1:50, reduction='pca')

saveRDS(x, "Data/Processed/alphasyn-triplication-initial-integrated.rds")