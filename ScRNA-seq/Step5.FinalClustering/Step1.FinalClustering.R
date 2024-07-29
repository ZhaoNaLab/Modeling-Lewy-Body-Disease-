library("harmony")
library("Seurat")
library("tidyverse")

setwd("./Step4.IterativeCleaning/CleanedClusters.R2")
SeuratObjs <- list.files(pattern = ".rds")

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.SNC <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:7, function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# further remove of contaminated cells
Seurat.SNC <- subset(Seurat.SNC, decontX_contamination < 0.4 & percent.mt <5)
Seurat.SNC <- DietSeurat(Seurat.SNC, assays = c("RNA", "decontX"), counts = TRUE)

# run harmony for data integration
options(repr.plot.height = 2.5, repr.plot.width = 6)

DefaultAssay(Seurat.SNC) <- "RNA"
Seurat.SNC <- FindVariableFeatures(Seurat.SNC, 
                                  selection.method = "vst", 
                                  nfeatures = 3000,
                                  verbose = F)
Seurat.SNC <- ScaleData(Seurat.SNC, 
                       vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                       verbose = FALSE)

Seurat.SNC <- RunPCA(Seurat.SNC, 
                    npcs = 30, 
                    pc.genes = Seurat.SNC@var.genes,
                    verbose = TRUE)


## integration based on individual
Seurat.SNC <- Seurat.SNC %>% RunHarmony("orig.ident", plot_convergence = TRUE)


## Identify clusters
Seurat.SNC <- Seurat.SNC %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()


setwd("./Projects/SNCSynuclein")
saveRDS(Seurat.SNC, "Step5.FinalClustering/Seurat.SNC.final.rds")

## end of the code
