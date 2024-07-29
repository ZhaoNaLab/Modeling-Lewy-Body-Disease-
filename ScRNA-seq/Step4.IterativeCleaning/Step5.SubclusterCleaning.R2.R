library("Seurat")
library("tidyverse")

# load the data
setwd("./Step4.IterativeCleaning")

# Ex
Ex <- readRDS("Recluster.R2/Ex.rds")

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

# no need for further cleaning
FeaturePlot(Ex, features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))
saveRDS(Ex, "CleanedClusters.R2/Ex.rds")


# In
In <- readRDS("Recluster.R2/In.rds")

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

FeaturePlot(In, features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))


In <- subset(In, seurat_clusters %in% c(0:20))
saveRDS(In, "CleanedClusters.R2/In.rds")


# OPC
OPC <- readRDS("Recluster.R2/OPC.rds")

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

OPC <- subset(OPC, seurat_clusters %in% c(0:5))
saveRDS(OPC, "CleanedClusters.R2/OPC.rds")

# Ast
Ast <- readRDS("Recluster.R2/Ast.rds")

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

Ast <- subset(Ast, seurat_clusters %in% c(0:9))
saveRDS(Ast, "CleanedClusters.R2/Ast.rds")


# Mic
Mic <- readRDS("Recluster.R2/Mic.rds")

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

Mic <- subset(Mic, seurat_clusters %in% c(0:6))
saveRDS(Mic, "CleanedClusters.R2/Mic.rds")

# Olig
Olig <- readRDS("Recluster.R2/Olig.rds")

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

# no need for further cleaning
FeaturePlot(Olig, features = c("nFeature_RNA", "decontX_contamination", "percent.mt",
                               "SLC17A7"))

Olig <- subset(Olig, seurat_clusters %in% c(0:9, 12, 15))
saveRDS(Olig, "CleanedClusters.R2/Olig.rds")


# Vascular.rds
Vascular <- readRDS("Recluster.R2/Vascular.rds")

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

Vascular <- subset(Vascular, seurat_clusters %in% c(0:5, 7:8))
saveRDS(Vascular, "CleanedClusters.R2/Vascular.rds")




