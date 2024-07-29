library("batchtools")


# prepare for the list of directory for data import
DirList <- list()

for (i in 1:8){
  DirList[[i]] <- paste0("./mrnaseq/221006-T3-NZ/", i, "/")
}

DirList <- unlist(DirList)

## use the bath tools for RNA contamination detection and doublet removal
setwd("./Step1.DataIntegration")
reg <- makeExperimentRegistry("SnRNAseq.Registry",
                              packages = c("tidyverse", "Seurat", "data.table",
                                           "singleCellTK", "purrr", "DoubletFinder",
                                           "celda", "ggpubr", "biomaRt"), seed = 20221005)

reg <- loadRegistry("SnRNAseq.Registry", writeable = TRUE)

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  sampleID <- i
  sce <- importCellRanger(sampleDirs = data[i])
  sce.raw <- importCellRanger(sampleDirs = data[i], dataType = "raw")
  sce <- decontX(sce, background = sce.raw)
  rownames(sce) <- rowData(sce)$feature_name
  Seurat.obj <- as.Seurat(sce, counts = "counts", data = NULL)
  
  # add another assay to the Seurat object
  Seurat.obj[["decontX"]] <- CreateAssayObject(counts = assay(sce, "decontXcounts"))
  
  # rename the assays
  Seurat.obj <- RenameAssays(object = Seurat.obj, originalexp = 'RNA')
  
  # rename orig.ident
  Seurat.obj$orig.ident <- paste0("S", Seurat.obj$orig.ident)
  Seurat.obj@meta.data <- Seurat.obj@meta.data[, c(1, 12, 11, 7)]
  list(data = data, Seurat.obj = Seurat.obj, sampleID = sampleID)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "snRNAseqDat", 
           data = DirList, 
           fun = subsample) 


# Define tasks
snRNAseqClean <- function(data, job, instance) { # instance is inherited from the last step
  data = instance$Seurat.obj
  sampleID = instance$sampleID
  
  ####### remove cells with features <200 and features expressed by less than 5 cells #######
  SelectedCells <- WhichCells(data, expression = nFeature_RNA > 200)
  SelectedFeatures <- rownames(data)[Matrix::rowSums(data) > 5]
  data <- subset(data, features = SelectedFeatures, cells = SelectedCells)
  
  # remove cell with a mitochondria percentage >= 15%
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data[["percent.rb"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
  SelectedMito <- WhichCells(data, expression = percent.mt < 15)
  data <- subset(data, cells = SelectedMito)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  data <- ScaleData(data, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = FALSE)
  data <- RunPCA(data, verbose = F, npcs = 20)
  data <- RunUMAP(data, dims = 1:20, verbose = F)

  # run parameter optimization with paramSweep
  sweep.res <- paramSweep_v3(data) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # define the expected number of doublet cellscells.
  nExp <- round(ncol(data) * 0.04)  # expect 4% doublets
  data <- doubletFinder_v3(data, pN = 0.25, pK = mpK, nExp = nExp, PCs = 1:15)
  colnames(data@meta.data)[grepl("DF.classification", colnames(data@meta.data))] <- "DFClass"
  colnames(data@meta.data)[grepl("pANN", colnames(data@meta.data))] <- "pANN"
  
  # Save the individual Seurat file
  SaveDir <- "./Step1.DataIntegration/SeuratObj.Ind"
  saveRDS(data, paste0(SaveDir, "/Seurat.obj.S", sampleID, ".rds"))
}


addAlgorithm(name = "snRNAseqClean", fun = snRNAseqClean)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:8)

# final problem design
pdes <- list(snRNAseqDat = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(snRNAseqClean = data.table())

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()

# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs() # ids = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56)
getStatus()
getJobTable()
getErrorMessages()
findExpired() 
