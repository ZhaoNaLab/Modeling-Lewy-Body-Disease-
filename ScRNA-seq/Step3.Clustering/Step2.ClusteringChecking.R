library("Seurat")
library("tidyverse")

setwd("./Projects/SNCSynuclein")

# Intergrated datasets at the individual level
Seurat.SNC.ind <- readRDS("Step3.Clustering/Seurat.SNC.Ind.rds")

# Dimplot by individual
DimPlot(Seurat.SNC.ind, label = T, split.by = "orig.ident") &
  NoLegend()

# Featureplot
FeaturePlot(Seurat.SNC.ind, features = c("nFeature_RNA", "decontX_contamination", "percent.mt"))


# Intergrated datasets at the batch level
Seurat.SNC.batch <- readRDS("Step3.Clustering/Seurat.SNC.batch.rds")

# Dimplot by individual
DimPlot(Seurat.SNC.batch, label = T, split.by = "orig.ident") &
  NoLegend()

# Featureplot
FeaturePlot(Seurat.SNC.batch, features = c("nFeature_RNA", "decontX_contamination", "percent.mt"))

## Note: based on the result, we can use the batch integration for downstream analysis








