library("Seurat")
library("tidyverse")

setwd("./Step5.FinalClustering")
Seurat.SNC.final <- readRDS("Seurat.SNC.final.rds")

# Step1: Update the marker genes for microglia
# marker gene expression
NormalizeData(Seurat.SNC.final, assay = "RNA")
Seurat.SNC.final$MajorCluster <- factor(Seurat.SNC.final$MajorCluster, 
                                       levels = c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Vascular"))

DefaultAssay(Seurat.SNC.final) <- "RNA"
Idents(Seurat.SNC.final) <- Seurat.SNC.final$MajorCluster

pdf("MarkerPlot1.pdf", width = 10.0, height = 5.0)
DotPlot(Seurat.SNC.final, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7",
                                      "GAD1", "GAD2", "GFAP", "AQP4", "GJA1",
                                      "ALDH1L1", "MOBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSPG4", "BCAN", "CSF1R", "CX3CR1", 
                                      "TMEM119","CLDN5", "FLT1", "PDGFRB"),
        cols = c("lightgrey", "#c04b57"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

dev.off()

## Export the metadata for SNC
Seurat.SNC.final$Dx <- NA
Seurat.SNC.final$Dx[grepl("LBD", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "LBD"
Seurat.SNC.final$Dx[grepl("Normal", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "Ctrl"
Seurat.SNC.final$Dx <- factor(Seurat.SNC.final$Dx, levels = c("Ctrl", "LBD"))

Metadata.Cell <- Seurat.SNC.final@meta.data
write.csv(Metadata.Cell, "Metadata.Cell.csv", row.names = TRUE)

## Visualize SCN and other genes
# a. by Dx
GeneofInterest <- c("SNCA", "LRRK2", "VPS35", "EIF4G1", "DNAJC13", "CHCHD2", "PARKIN", 
                    "PINK1", "PARK7", "TMEM175", "GBA", "BCL7C", "STX1B", "GABRB3",
                    "HERC5", "HERC3", "NAP1L5", "FAM13A", "TIGD2", "GPRIN3", "SNCA", 
                    "MMRN1", "PIGY", "FAM13AOS", "LOC644248", "PYURI")

# select genes that are in the rownames
GeneofInterest <- intersect(rownames(Seurat.SNC.final), GeneofInterest)

for (i in 1:length(GeneofInterest)){
  p <- VlnPlot(Seurat.SNC.final,
               features = GeneofInterest[i],
               group.by = "MajorCluster",
               split.by = "Dx",
               split.plot = TRUE,
               pt.size = 0,
               ncol = 1)
  
  p <- p + geom_jitter(data = p$data, 
                       aes(color = split),
                       position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
                       size = 0.1, alpha = 0.15) +
    guides(color = "none") +
    scale_color_manual(values = c("black", 'black'))
  
  ggsave(paste0("Figures/", GeneofInterest[i], ".Vlnplot.by.Dx.pdf"), width = 6.5, height = 4.0, plot = p)
}


## Feature plot
for (i in 1:length(GeneofInterest)){
  p <- FeaturePlot(Seurat.SNC.final,
                   features = GeneofInterest[i],
                   split.by = "Dx",
                   cols = c("gray", "#ff0000"),
                   pt.size = 0.01, 
                   order = TRUE,
                   ncol = 1) &
    xlab("UMAP-1") &
    ylab("UMAP-2")
  
  ggsave("Figures/", paste0(GeneofInterest[i], ".UMAP.by.Dx.pdf"), 
         width = 8.0, height = 4.0, plot = p)
}


# b. by PathDx
Seurat.SNC.final$PathDx <- factor(Seurat.SNC.final$PathDx, levels = c("Normal", "DLBD",
                                                                    'TLBD', "Dupl_DLBD",
                                                                    "Tri_DLBD"))
for (i in 1:length(GeneofInterest)){
  p <- VlnPlot(Seurat.SNC.final,
                features = GeneofInterest[i],
                group.by = "MajorCluster",
                split.by = "PathDx",
                # split.plot = TRUE,
                pt.size = 0,
                ncol = 1) +
    ggsci::scale_fill_futurama()
  
  p <- p + geom_jitter(data = p$data, 
                         aes(color = split),
                         position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
                         size = 0.05, alpha = 0.1) +
    guides(color = "none") +
    scale_color_manual(values = c("black", 'black', 'black', 'black', 'black'))
  
  
  ggsave(paste0("Figures/", GeneofInterest[i], ".Vlnplot.by.PathoDx.pdf"), 
         width = 8.0, height = 4.0, plot = p)
  
}


## Feature plot
for (i in 1:length(GeneofInterest)){
  p <- FeaturePlot(Seurat.SNC.final,
                      features = GeneofInterest[i],
                      split.by = "PathDx",
                      cols = c("gray", "#ff0000"),
                      pt.size = 0.01, 
                      order = TRUE,
                      ncol = 3) &
    xlab("UMAP-1") &
    ylab("UMAP-2")
  
  ggsave(paste0("Figures/", GeneofInterest[i], ".UMAP.by.PathoDx.pdf"), 
         width = 15, height = 3.0, plot = p)
}

table(Seurat.SNC.final$MajorCluster, Seurat.SNC.final$orig.ident)

## Combine DLBD and TLBD
Metadata.Cell$PathDx2 <- Metadata.Cell$PathDx
Metadata.Cell$PathDx2[Metadata.Cell$PathDx %in% c("DLBD", "TLBD")] <- "LBD"

## update the metadata
Seurat.SNC.final <- AddMetaData(Seurat.SNC.final, metadata = Metadata.Cell)
Seurat.SNC.final$PathDx2 <- factor(Seurat.SNC.final$PathDx2, levels = c("Normal", "LBD",
                                                                      "Dupl_DLBD", "Tri_DLBD"))
## 3. by PathDx2
for (i in 1:length(GeneofInterest)){
  p <- VlnPlot(Seurat.SNC.final,
               features = GeneofInterest[i],
               group.by = "MajorCluster",
               split.by = "PathDx2",
               # split.plot = TRUE,
               pt.size = 0,
               ncol = 1) +
    ggsci::scale_fill_futurama()
  
  p <- p + geom_jitter(data = p$data, 
                       aes(color = split),
                       position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
                       size = 0.05, alpha = 0.1) +
    guides(color = "none") +
    scale_color_manual(values = c("black", 'black', 'black', 'black', 'black'))
  
  
  ggsave(paste0("Figures/", GeneofInterest[i], ".Vlnplot.by.PathoDx2.pdf"), 
         width = 7.5, height = 4.0, plot = p)
  
}


## Feature plot
for (i in 1:length(GeneofInterest)){
  p <- FeaturePlot(Seurat.SNC.final,
                   features = GeneofInterest[i],
                   split.by = "PathDx2",
                   cols = c("gray", "#ff0000"),
                   pt.size = 0.01, 
                   order = TRUE,
                   ncol = 3) &
    xlab("UMAP-1") &
    ylab("UMAP-2")
  
  ggsave(paste0("Figures/", GeneofInterest[i], ".UMAP.by.PathoDx2.pdf"), 
         width = 14, height = 3.0, plot = p)
}
