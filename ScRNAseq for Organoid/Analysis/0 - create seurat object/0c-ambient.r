# *********
# 0c-ambient.r
# Script for estimating ambient RNA contamination in the dataset with DecontX.
# Expected within the folder Data/secondary is are folders BU-SNCA-1, BU-SNCA-2, ..., BU-SNCA-24,
# within which we expect the folder structure of:
# outs/
#   - filtered_feature_bc_matrix
#   - raw_feature_bc_matrix
# from 10X scRNAseq cellranger
library(celda)
library(data.table)
library(ggplot2)
library(magrittr)
library(dplyr)

orig.idents <- list.files('Data/secondary/') # paste0("BU-SNCA-", 1:24)
for (orig.ident in orig.idents) {
  filt_counts <- Read10X(data.dir=paste0('Data/secondary/', orig.ident, '/outs/filtered_feature_bc_matrix'))
  x <- CreateSeuratObject(counts=filt_counts, project=orig.ident)
  raw_counts <- Read10X(data.dir=paste0('Data/secondary/', orig.ident, '/outs/raw_feature_bc_matrix'))
  raw <- CreateSeuratObject(counts=raw_counts, project=orig.ident)
  raw <- as.SingleCellExperiment(raw)
  sce <- as.SingleCellExperiment(x)
  
  sce <- decontX(sce, background=raw)
  x <- as.Seurat(x=sce)
  x[['RNA_decontX']] <- CreateAssayObject(counts=SummarizedExperiment::assay(sce, 'decontXcounts'))
  saveRDS(x, paste0("Data/Initial/decontX/", orig.ident, '.rds'))
}