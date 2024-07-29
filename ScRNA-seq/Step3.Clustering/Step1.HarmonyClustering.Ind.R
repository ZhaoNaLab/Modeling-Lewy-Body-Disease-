library("harmony")
library("Seurat")
library("tidyverse")

setwd("./Projects/SNCSynuclein")

# load the Seurat object
Seurat.SNC <- readRDS("Step2.QC/Seurat.SNC.QC.R1.rds")

# run harmony for data integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
Seurat.SNC <- FindVariableFeatures(Seurat.SNC, 
                                  selection.method = "vst", 
                                  nfeatures = 3000,
                                  verbose = F)
Seurat.SNC <- ScaleData(Seurat.SNC, 
                       vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                       verbose = FALSE)

Seurat.SNC <- RunPCA(Seurat.SNC, 
                    npcs = 50, 
                    pc.genes = Seurat.SNC@var.genes,
                    verbose = TRUE)


## integration based on individual
Seurat.SNC <- Seurat.SNC %>% RunHarmony("orig.ident", plot_convergence = TRUE)


## Identify clusters
Seurat.SNC <- Seurat.SNC %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

saveRDS(Seurat.SNC, "Step3.Clustering/Seurat.SNC.Ind.rds")

## end of the code
