library('RColorBrewer')
library('MAST')
library('data.table')
library('Seurat')
library('biomaRt')
library('lme4')

setwd("./Projects/SNCSynuclein")

# load the seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# set the levels for factors:
Seurat.SNC.final$TYPE <- factor(Seurat.SNC.final$TYPE, levels = c("Ctrl", "LBD"))

# filter the seurat object
DefaultAssay(Seurat.SNC.final) <- 'RNA'

################################# construct singlecellassay object ####################################
latent.vars <- Seurat.SNC.final@meta.data
latent.vars$wellKey <- rownames(x = latent.vars)

# make sure the expression matrix and the metadata have the same order
latent.vars <- latent.vars[match(colnames(Seurat.SNC.final), latent.vars$wellKey), ]

# prepare the fdat
fdat <- data.frame(rownames(x = Seurat.SNC.final))
colnames(x = fdat)[1] <- "primerid"
rownames(x = fdat) <- fdat[, 1]

# construct the SingleCellAssay object
sca <- MAST::FromMatrix(
  exprsArray = as.matrix(Seurat.SNC.final[['RNA']]@counts),
  check_sanity = FALSE,
  cData = latent.vars,
  fData = fdat
)

saveRDS(sca, "Step5.FinalClustering/scaObjs.rds")



