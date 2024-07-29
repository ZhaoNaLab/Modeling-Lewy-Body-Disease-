library("Seurat")
library("tidyverse")

# load the data
setwd("./Step4.IterativeCleaning")

# Ex
Ex <- readRDS("Recluster.R1/Ex.rds")

# inspect the data
DimPlot(Ex, label = T)

# marker gene expression
NormalizeData(Ex, assay = "RNA")

DotPlot(Ex, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                         "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                         "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                         "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                         "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

Ex <- subset(Ex, seurat_clusters %in% c(0:11, 13, 15:17, 19, 23, 26:28))
saveRDS(Ex, "CleanedClusters.R1/Ex.rds")


# In
In <- readRDS("Recluster.R1/In.rds")

# inspect the data
DimPlot(In, label = T)

# marker gene expression
NormalizeData(In, assay = "RNA")
DotPlot(In, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                         "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                         "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                         "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                         "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

In <- subset(In, seurat_clusters %in% c(0:3, 5:7, 9:20, 22, 26))
saveRDS(In, "CleanedClusters.R1/In.rds")


# OPC
OPC <- readRDS("Recluster.R1/OPC.rds")

# inspect the data
DimPlot(OPC, label = T)

# marker gene expression
NormalizeData(OPC, assay = "RNA")
DotPlot(OPC, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                          "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                          "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                          "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                          "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

OPC <- subset(OPC, seurat_clusters %in% c(0:2, 4, 7))
saveRDS(OPC, "CleanedClusters.R1/OPC.rds")

# Ast
Ast <- readRDS("Recluster.R1/Ast.rds")

# inspect the data
DimPlot(Ast, label = T)

# marker gene expression
NormalizeData(Ast, assay = "RNA")
DotPlot(Ast, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                          "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                          "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                          "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                          "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

Ast <- subset(Ast, seurat_clusters %in% c(0:2, 4, 6, 10:12, 15))
saveRDS(Ast, "CleanedClusters.R1/Ast.rds")


# Mic
Mic <- readRDS("Recluster.R1/Mic.rds")

# inspect the data
DimPlot(Mic, label = T)

# marker gene expression
NormalizeData(Mic, assay = "RNA")
DotPlot(Mic, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                          "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                          "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                          "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                          "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

Mic <- subset(Mic, seurat_clusters %in% c(0, 2:6, 8:9, 13, 15:17))
saveRDS(Mic, "CleanedClusters.R1/Mic.rds")

# Olig
Olig <- readRDS("Recluster.R1/Olig.rds")

# inspect the data
DimPlot(Olig, label = T)

# marker gene expression
NormalizeData(Olig, assay = "RNA")
DotPlot(Olig, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                           "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                           "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                           "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                           "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

Olig <- subset(Olig, seurat_clusters %in% c(0:13))
saveRDS(Olig, "CleanedClusters.R1/Olig.rds")


# Vascular.rds
Vascular <- readRDS("Recluster.R1/Vascular.rds")

# inspect the data
DimPlot(Vascular, label = T)

# marker gene expression
NormalizeData(Vascular, assay = "RNA")
DotPlot(Vascular, features = c("RBFOX3","GABRB2", "SATB2", "SLC17A7", "GAD1", "GAD2", "GFAP", "AQP4",
                           "PLP1", "MOBP", "MBP", "VCAN", "CSF1R", "CD74", "C3", "CLDN5",
                           "FLT1", "PDGFRB"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

Vascular <- subset(Vascular, seurat_clusters %in% c(0:4, 6, 8:9))
saveRDS(Vascular, "CleanedClusters.R1/Vascular.rds")


