library("Seurat")
library("tidyverse")

setwd("./Step1.DataIntegration/SeuratObj.Ind")
SeuratObjs <- list.files()

# load the files
Seurat.list <- lapply(SeuratObjs, readRDS)

# Merge into one single Seurat object
Seurat.SNC <- merge(x = Seurat.list[[1]], y = unlist(lapply(2:8, function(x){Seurat.list[[x]]})))
rm(Seurat.list)

# add the individual information to the metadata
Meta.cell <- Seurat.SNC@meta.data %>% 
  rownames_to_column("Cells") 

## load the metadata
Metadata.SNC <- read.csv("./Step1.DataIntegration/Metadata.SNC.csv")
Metadata.SNC$zAge <- scale(Metadata.SNC$Age)
Metadata.SNC$orig.ident <- paste0("S", Metadata.SNC$ID)

# add sample information to the metadata
Meta.cell <- Meta.cell %>% left_join(Metadata.SNC, by = "orig.ident") %>%
  column_to_rownames("Cells")


## Add the metadata to the Seurat object
Seurat.SNC <- AddMetaData(Seurat.SNC, metadata = Meta.cell)

# remove some varibles
Seurat.SNC$nFeature_decontX <- NULL
Seurat.SNC$nCount_decontX <- NULL

## make the new varible cngeneson (scaled cell detection rate)
Seurat.SNC$cngeneson <- scale(Seurat.SNC$nFeature_RNA)

## reorder the varibles
Seurat.SNC@meta.data <- Seurat.SNC@meta.data %>% dplyr::select(-ID)

setwd("./Step1.DataIntegration")
saveRDS(Seurat.SNC, "Seurat.SNC.Single.Double.rds")

## The end of the code







