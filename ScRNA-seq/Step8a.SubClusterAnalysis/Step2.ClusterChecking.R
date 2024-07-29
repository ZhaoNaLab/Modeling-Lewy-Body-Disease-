library("Seurat")
library("tidyverse")

setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis")

############################## Ex ###############################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/Ex")
Ex.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(Ex.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(Ex.Seurat, DimPlot)

## After a comprehensive comparision Ex.Seurat Ex.Harmony.2.30 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/Ex")
Ex.Seurat <- readRDS("Ex.Harmony.rds")
DimPlot(Ex.Seurat)

IDs <- paste0("S", 1:length(unique(Ex.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(Ex.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}

############################## In ###############################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/In")
In.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(In.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(In.Seurat, DimPlot)

## After a comprehensive comparision In.Seurat In.Harmony.2.30 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/In")
In.Seurat <- readRDS("In.Harmony.rds")
DimPlot(In.Seurat)

IDs <- paste0("S", 1:length(unique(In.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(In.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}

############################## Ast ###############################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/Ast")
Ast.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(Ast.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(Ast.Seurat, DimPlot)

## take Ast.Harmony.1.50 as the final 
## After a comprehensive comparision Ast.Seurat Ast.Harmony.1.50 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/Ast")
Ast.Seurat <- readRDS("Ast.Harmony.rds")
DimPlot(Ast.Seurat)

IDs <- paste0("S", 1:length(unique(Ast.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(Ast.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}

############################## Microglia #############################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/Mic")
Mic.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(Mic.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(Mic.Seurat, DimPlot)

## take Mic.Harmony.1.50 as the final
Mic.Seurat <- Mic.Seurat$Mic.Harmony.2.30

## After a comprehensive comparision Mic.Seurat Mic.Harmony.2.30 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/Mic")
Mic.Seurat <- readRDS("Mic.Harmony.rds")
DimPlot(Mic.Seurat)
FeaturePlot(Mic.Seurat, 
            features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))


IDs <- paste0("S", 1:length(unique(Mic.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(Mic.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID.test/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}


################################################### Olig ####################################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/Olig")
Olig.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(Olig.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(Olig.Seurat, DimPlot)
Olig.Seurat <- Olig.Seurat$Olig.Harmony.2.20

## After a comprehensive comparision Olig.Seurat Olig.Harmony.2.20 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/Olig")
Olig.Seurat <- readRDS("Olig.Harmony.rds")
DimPlot(Olig.Seurat)
FeaturePlot(Olig.Seurat, 
            features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))

IDs <- paste0("S", 1:length(unique(Olig.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(Olig.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}

############################## OPC ###############################
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA.SNC/OPC")
OPC.Seurat <- lapply(list.files(pattern = "*.rds"), readRDS)
names(OPC.Seurat) <- gsub(".rds", "", list.files(pattern = "*.rds"))

P <- lapply(OPC.Seurat, DimPlot)
OPC.Seurat <- OPC.Seurat$OPC.Harmony.2.50

## After a comprehensive comparision OPC.Seurat OPC.Harmony.2.50 is the optimal one
setwd("./Projects/snRNAseq/Step12.SubclusterAnalysis/OPC")
OPC.Seurat <- readRDS("OPC.Harmony.rds")

DimPlot(OPC.Seurat)
FeaturePlot(OPC.Seurat, 
            features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))

IDs <- paste0("S", 1:length(unique(OPC.Seurat$orig.ident)))
IDchuncks <- split(IDs, ceiling(seq_along(1:length(IDs))/4))

# make the figure look better
IDchuncks$`15` <- c("S55", "S56", "S57", "S58")
for (i in 1:length(IDchuncks)){
  dat <- subset(OPC.Seurat, orig.ident %in% unlist(IDchuncks[i]))
  Idents(dat) <- dat$seurat_clusters
  
  jpeg(paste0("QC/SubclusterByID/DimPlo.Chunck", i,".jpeg"), width = 16, height = 5, unit = 'in', res = 600)
  
  p <- DimPlot(dat, split.by = 'orig.ident', label = T, ncol = 4,
               # label.size = 8.0,
               raster = FALSE) +
    xlab("UMAP-1") +
    ylab("UMAP-2") +
    theme(axis.title = element_text(size = 20, color = 'black'),
          axis.text= element_text(size = 20, color = 'black'),
          axis.line = element_line(linewidth = 1.0),
          axis.ticks = element_line(linewidth = 1.0),
          axis.ticks.length = unit(0.25, "cm")) &
    NoLegend()
  
  print(p)
  dev.off()
}
