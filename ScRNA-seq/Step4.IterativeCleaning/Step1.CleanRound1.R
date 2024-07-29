library("Seurat")
library("tidyverse")

setwd("./Projects/SNCSynuclein")

# Intergrated datasets at the batch level
Seurat.SNC.Ind <- readRDS("Step3.Clustering/Seurat.SNC.Ind.rds")

## Recluster with a higher resolution
Seurat.SNC.Ind <- Seurat.SNC.Ind %>% 
  FindClusters(resolution = 0.9) %>% 
  identity()

# visualize the clusters
DimPlot(Seurat.SNC.Ind, label = T)
table(Seurat.SNC.Ind$RNA_snn_res.0.9)

## Remove clusters expressing marker genes of different cell type
DefaultAssay(Seurat.SNC.Ind) <- "RNA"
Seurat.SNC.Ind <- NormalizeData(Seurat.SNC.Ind)

# marker gene expression
DotPlot(Seurat.SNC.Ind, features = c("RBFOX3","GABRB2", "SATB2", "SLC17A7", "GAD1", "GAD2", "GFAP", "AQP4",
                                      "PLP1", "MOBP", "MBP", "VCAN", "CSF1R", "CD74", "C3", "CLDN5",
                                      "FLT1", "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

## Define the major clusters based on markers
Seurat.SNC.Ind$MajorCluster <- NA
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(0, 1, 4, 10, 17)] <- "Olig"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(2, 5, 6, 9, 11, 14, 18, 19, 21, 24:28, 30)] <- "Ex"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(3, 16)] <- "Ast"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(8, 23)] <- "Mic"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(12, 13, 15, 20, 22, 31)] <- "In"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(7)] <- "OPC"
Seurat.SNC.Ind$MajorCluster[Seurat.SNC.Ind$seurat_clusters %in% c(29, 32)] <- "Vascular"

## save the data
saveRDS(Seurat.SNC.Ind, "Step3.Clustering/Seurat.SNC.Ind.rds")

## end of the code ##



