# *******
# 0a-import.r
# Code for initial construction of the seurat object holding UMI
# Note: Expected within Data/secondary are individual secondary 10X outputs for each organoid scRNAseq sample.
library(tidyverse)
library(magrittr)
library(Seurat)
source("Analysis/util/import_10x.r")
meta <- readxl::read_xlsx('Data/Meta/Sample information-SNCA triplicaiton project -- updated.xlsx', range = 'B2:H26')
meta <- meta %>%
  rename(orig.ident=ID,
         collection.date=`Date of collection`,
         line.id=`iPSC Line and subclone ID`,
         flow.cell=`Flow Cell`)
x <- create_merged_seurat_object(secondary.path='Data/secondary/', meta=meta)
x$orig.ident <- factor(x$orig.ident, paste0('BU-SNCA-',1:24))
x$collection.date <- as.Date(x$collection.date / (3600*24), '1970-01-01')
x$collection.date <- factor(x$collection.date, as.Date(sort(unique(x$collection.date)), '1970-01-01'))
saveRDS(x, 'Data/Initial/alphasyn-triplication-initial-seurat.rds')