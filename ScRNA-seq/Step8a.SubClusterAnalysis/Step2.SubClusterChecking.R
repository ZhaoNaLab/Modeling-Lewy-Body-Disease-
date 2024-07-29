library("Seurat")

setwd("./Step8a.SubClusterAnalysis")

## Excitatory neurons
Ex.SNC <- readRDS("Ex/Ex.Harmony.Ind.rds")
Idents(Ex.SNC) <- Ex.SNC$RNA_snn_res.0.2
DimPlot(Ex.SNC)
DimPlot(Ex.SNC, split.by = "orig.ident")
Ex.SNC$seurat_clusters <- Ex.SNC$RNA_snn_res.0.2
Ex.SNC$Subcluster <- paste0(Ex.SNC$MajorCluster, ".", Ex.SNC$seurat_clusters)
saveRDS(Ex.SNC, "Ex/Ex.Harmony.Ind.rds")

## Inhibitory neurons
In.SNC <- readRDS("In/In.Harmony.Ind.rds")
Idents(In.SNC) <- In.SNC$RNA_snn_res.0.2
DimPlot(In.SNC)
DimPlot(In.SNC, split.by = "orig.ident")
In.SNC$seurat_clusters <- In.SNC$RNA_snn_res.0.2
In.SNC$Subcluster <- paste0(In.SNC$MajorCluster, ".", In.SNC$seurat_clusters)
saveRDS(In.SNC, "In/In.Harmony.Ind.rds")

## Oligodendrocytes
Olig.SNC <- readRDS("Olig/Olig.Harmony.Ind.rds")
Idents(Olig.SNC) <- Olig.SNC$RNA_snn_res.0.2
DimPlot(Olig.SNC)
DimPlot(Olig.SNC, split.by = "orig.ident")

Olig.SNC <- Olig.SNC %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.1, 0.15)) %>% 
  identity()




Olig.SNC$seurat_clusters <- Olig.SNC$RNA_snn_res.0.1
Olig.SNC$Subcluster <- paste0(Olig.SNC$MajorCluster, ".", Olig.SNC$seurat_clusters)
saveRDS(Olig.SNC, "Olig/Olig.Harmony.Ind.rds")

## OPC
OPC.SNC <- readRDS("OPC/OPC.Harmony.Ind.rds")
Idents(OPC.SNC) <- OPC.SNC$RNA_snn_res.0.2
DimPlot(OPC.SNC, split.by = "orig.ident")

OPC.SNC$seurat_clusters <- OPC.SNC$RNA_snn_res.0.2
OPC.SNC$Subcluster <- paste0(OPC.SNC$MajorCluster, ".", OPC.SNC$seurat_clusters)
saveRDS(OPC.SNC, "OPC/OPC.Harmony.Ind.rds")


## Astrocytes
Ast.SNC <- readRDS("Ast/Ast.Harmony.Ind.rds")
Idents(Ast.SNC) <- Ast.SNC$RNA_snn_res.0.2
DimPlot(Ast.SNC, split.by = "orig.ident")
DimPlot(Ast.SNC, split.by = "TYPE")

FeaturePlot(Ast.SNC, features = c("APOE"), split.by = "TYPE")
Ast.SNC$seurat_clusters <- Ast.SNC$RNA_snn_res.0.2
Idents(Ast.SNC) <- Ast.SNC$seurat_clusters
saveRDS(Ast.SNC, "Ast/Ast.Harmony.Ind.rds")


## Microglia
Mic.SNC <- readRDS("Mic/Mic.Harmony.Ind.rds")
Idents(Mic.SNC) <- Mic.SNC$RNA_snn_res.0.2
DimPlot(Mic.SNC, split.by = "orig.ident")
DimPlot(Mic.SNC, split.by = "TYPE")
Mic.SNC$seurat_clusters <- Mic.SNC$RNA_snn_res.0.2
Idents(Mic.SNC) <- Mic.SNC$seurat_clusters
saveRDS(Mic.SNC, "Mic/Mic.Harmony.Ind.rds")

# ## vasculature
# Vascular.SNC <- readRDS("Vascular/Vascular.Harmony.Ind.rds")
# Idents(Vascular.SNC) <- Vascular.SNC$RNA_snn_res.0.2
# DimPlot(Vascular.SNC, split.by = "orig.ident")
# Vascular.SNC$seurat_clusters <- Vascular.SNC$RNA_snn_res.0.2
# Idents(Vascular.SNC) <- Vascular.SNC$seurat_clusters
# saveRDS(Vascular.SNC, "Vascular/Vascular.Harmony.Ind.rds")


























