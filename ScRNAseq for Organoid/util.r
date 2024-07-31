# util.r
library(dplyr)
library(magrittr)
library(Seurat)
library(ggplot2)
library(forcats)
library(stringr)
library(patchwork)
library(RColorBrewer)
library(grDevices)

celltype.levels = c('Neuron 0', 'Neuron 5', 'CN 1', 'CN 2', 'CN 3/Photo 6', 'Inhibitory neuron 8', 'proRG 3', 'RG 4', 'RG 9', 'Astro 7')
cell.broadtype.levels = c( 'Neuron (ExN immature) 0', 'Neuron (InN immature) 5', 'ExN 1', 'ExN 2', 'ExN 3/Photo 6', 'InN 8', 'proRG 3', 'RG 4', 'RG 9', 'AS 0', 'Astro 7')
Genotype.levels = c('SNCA trip', 'Ctrl')

load_sampleinfo <- function() {
  readxl::read_xlsx('Data/Meta/Sample information-SNCA triplicaiton project.xlsx', range = 'B2:H26') %>%
  rename(
    orig.ident=ID,
    collection.date=`Date of collection`,
    line.id=`iPSC Line and subclone ID`,
    flow.cell=`Flow Cell`
    ) %>% 
  mutate(
    orig.ident = factor(orig.ident, paste0('BU-SNCA-',1:24)),
    collection.date = lubridate::as_date(collection.date)
    #collection.date = as.Date(collection.date / (3600*24), '1970-01-01'),
    #collection.date = factor(collection.date, as.Date(sort(unique(collection.date)), '1970-01-01')) 
    ) %>% mutate(line.id = ifelse(orig.ident %in% paste0('BU-SNCA-', 4:6), 'MC0117 #7', line.id))
}

load_integrated_metadata <- function(meta.data = NULL) {
  if (is.null(meta.data)) {
    meta.data <- data.table::fread('Results/Processed/cell_metadata.csv') %>%
      tibble::column_to_rownames('V1')
  }
  meta.data <- meta.data %>%
    mutate(cell.broadtype = paste0(celltype, ' ', SCT.seurat_snn_res.0.1),
                               cell.broadtype = case_when(
                                 cell.broadtype == 'CN 6' ~ 'CN 3/Photo 6',
                                 cell.subtype == 'AS 0' ~ 'AS 0',
                                 #stringr::str_detect(cell.broadtype, 'AS') ~ cell.subtype,
                                 TRUE ~ cell.broadtype),
           cell.broadtype = case_when(
             cell.broadtype == 'Neuron (immature) 5' ~ 'Neuron 5',
             cell.broadtype == 'Intermediate cell 0' ~ 'Neuron 0',
             cell.broadtype == 'NEC 3' ~ 'DV 3',
             TRUE ~ cell.broadtype
             ),
           # keep these levels -- seems to be specific to above naming scheme
           cell.broadtype = factor(cell.broadtype, levels=c( 'Neuron 0', 'Neuron 5', 'CN 1', 'CN 2', 'CN 3/Photo 6', 'Inhibitory neuron 8', 'GPC 4', 'AS 0', 'AS 7', 'DV 3', 'PGC 9')),
           orig.ident = factor(orig.ident, paste0("BU-SNCA-", c(1:24))),
           orig.ident = droplevels(orig.ident)) 
  return(meta.data)
}

load_integrated <- function() {
  x <- readRDS("Data/Roundtwo-Integration/alphasyn-triplication-final-roundtwo-integrated.rds")
  x@meta.data <- load_integrated_metadata(re@meta.data) 
  x <- subset(x, cell.broadtype != 'AS 0')
  x$cell.broadtype <- droplevels(x$cell.broadtype)
  x$cell.broadtype <- stringr::str_replace(x$cell.broadtype, ' ', '.')
  return(x)
}
