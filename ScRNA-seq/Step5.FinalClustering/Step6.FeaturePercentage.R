library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("./Projects/SNCSynuclein")

############ Convert the Seurat object to SCA obejct for Marker calculation #########
# load the seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")
DefaultAssay(Seurat.SNC.final) <- 'RNA'

# Calcualte the percentage of cells expression the features
ExpMat <- GetAssayData(Seurat.SNC.final, assay = "RNA", slot = "counts")      
PercMatOverall <- as.matrix(rowMeans(ExpMat > 0))*100
colnames(PercMatOverall) <- "Percentage.Overall"

# Calcualte the percentage of Major clusters expression the features
DefaultAssay(Seurat.SNC.final) <- "RNA"
MajorCellType <- unique(Seurat.SNC.final$MajorCluster)
FeaturePerc <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(Seurat.SNC.final, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- paste0(MajorCellType[i], ".Percentage")
  FeaturePerc[[i]] <- PercMat
}

FeaturePerc.Comb <- do.call(cbind, FeaturePerc)

# all.equal(rownames(PercMatOverall), rownames(FeaturePerc.Comb))
FeatureAllinOne <- cbind(PercMatOverall, FeaturePerc.Comb)

## save the data
saveRDS(FeatureAllinOne, "Step5.FinalClustering/Feature.Cell.Percent.rds")


