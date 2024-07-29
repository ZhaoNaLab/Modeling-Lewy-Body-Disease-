library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("./Projects/SNCSynuclein")

# load the seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# set the levels for factors: Pathological Dx
Seurat.SNC.final$Dx <- NA
Seurat.SNC.final$Dx[grepl("LBD", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "LBD"
Seurat.SNC.final$Dx[grepl("Normal", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "Ctrl"
Seurat.SNC.final$Dx <- factor(Seurat.SNC.final$Dx, levels = c("Ctrl", "LBD"))


# load the gene filters
GeneSelect <- readRDS("Step6.DEG/FeatureSelect.rds")

# loop the transformation 
MajorClusters <- unique(Seurat.SNC.final$MajorCluster)
scaObjs <- list()

for (i in 1:length(MajorClusters)){
  SeuratObj <- subset(Seurat.SNC.final, MajorCluster == MajorClusters[i])
  
  # filter the seurat object
  DefaultAssay(SeuratObj) <- 'RNA'
  SeuratObj <- subset(SeuratObj, features = GeneSelect[[MajorClusters[i]]])
  
  ################################# construct singlecellassay object ####################################
  latent.vars <- SeuratObj@meta.data
  latent.vars$wellKey <- rownames(x = latent.vars)
  
  # make sure the expression matrix and the metadata have the same order
  latent.vars <- latent.vars[match(colnames(SeuratObj), latent.vars$wellKey), ]
  
  # prepare the fdat
  fdat <- data.frame(rownames(x = SeuratObj))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  
  # construct the SingleCellAssay object
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(SeuratObj[['RNA']]@counts),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  
  # double check the factor levels
  scaObjs[[i]] <- sca
  
}

names(scaObjs) <- MajorClusters
saveRDS(scaObjs, "Step6.DEG/scaObjs.rds")



