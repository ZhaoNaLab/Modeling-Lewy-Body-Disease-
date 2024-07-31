# Doublet Finder registry
source("Analysis/util.r")
library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(data.table)
x <- readRDS("Data/Processed/alphasyn-triplication-initial-integrated.rds")
x$orig.ident <- factor(x$orig.ident, paste0('BU-SNCA-', 1:24))

x_list <- SplitObject(x, split.by='orig.ident')
meta_list <- list()

num.cores <- 1
PCs <- 1:30
sct <- TRUE
GT <- FALSE
pN <- 0.25
label <- 'SCT.seurat_snn_res.0.1' # Alternatively, "seurat_clusters"
stopifnot(labels %in% colnames(x@meta.data))

for (orig_ident in paste0("BU-SNCA-", 1:24)) {
  xs <- x_list[[orig_ident]]
  DefaultAssay(xs) <- 'SCT'
  
  # Detect doublets with DoubletFinder
  ## pK identification (no ground-truth)
  
  sweep.res.list <- paramSweep_v3(x, PCs=PCs, sct=sct, num.cores=num.cores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=GT)
  bcmvn <- find.pK(sweep.stats)
  #print(bcmvn)
  pK = bcmvn[which.max(bcmvn$BCmetric),'pK']
  #print(pK)
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- x@meta.data[, labels]
  homotypic.prop <- modelHomotypic(annotations)

  # estimated multiplet rate for v3.1 10X chemistry
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry

  dat <- data.frame(
    multiplet_rate=c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6) * 0.01,
    number_cells_loaded=c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    number_cells_recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )

  # create a linear model from the above table, assuming the number of cells corresponds to number of cells recovered.
  # should work alright even if we extrapolate
  fit <- lm(multiplet_rate ~ number_cells_recovered, data=dat)
  a = as.numeric(fit$coefficients[1])
  b = as.numeric(fit$coefficients[2])

  print(paste('ncells:', ncells))
  tenX_doublet_prop = a + b*ncells
  print(paste('tenX_doublet_prop:', tenX_doublet_prop))
  nExp_poi <- round(tenX_doublet_prop*nrow(x@meta.data)) ## Assuming the above table of doublet proportions for v3.1
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies
  x <- DoubletFinder:::doubletFinder_v3(xs, PCs=PCs, pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=sct)
  pANN = colnames(xs@meta.data)[stringr::str_detect(colnames(xs@meta.data), 'pANN')]
  x <- DoubletFinder:::doubletFinder_v3(xs, PCs=PCs, pN=pN, pK=pK, nExp=nExp_poi.adj, reuse.pANN=pANN, sct=sct)
 
  # Clean up results and save metadata
  xs <- clean_doubletfinder(xs)
  meta_list[[orig_ident]] <- xs@meta.data
}
meta.data <- rbind(meta_list)
meta.data <- meta.data[colnames(x), ]
data.table::fwrite(meta.data, "Results/Initial/Initial integration/doublet-finder-metadata.csv", row.names=TRUE)

# Clean up doublet Finder results.
clean_doubletfinder = function(x, ...) {
  library(stringr)
  pANN = colnames(x@meta.data)[str_detect(colnames(x@meta.data), 'pANN')]
  x$DF.pANN <- x@meta.data[, pANN]
  DF.classifications = colnames(x@meta.data)[str_detect(colnames(x@meta.data), 'DF.classification')] 
  if (length(DF.classifications) == 2) {
    x$DF.classification.upper <- x@meta.data[, DF.classifications[1]]
    x$DF.classification.lower <- x@meta.data[, DF.classifications[2]]
  } else {
    x$DF.classification.upper <- x@meta.data[, DF.classifications]
  }
  x@meta.data <- x@meta.data[, !colnames(x@meta.data) %in% c(pANN, DF.classifications)]
  return(x)
}  