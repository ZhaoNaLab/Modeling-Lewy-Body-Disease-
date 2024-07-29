library("Seurat")
library("tidyverse")

setwd("./Step5.FinalClustering")
Seurat.SNC.final <- readRDS("Seurat.SNC.final.rds")

## DimPlot of the data
Idents(Seurat.SNC.final) <- Seurat.SNC.final$seurat_clusters

jpeg("Clusters.jpeg", width = 6.0, height = 6.0, unit = "in", res = 600)
DimPlot(Seurat.SNC.final, label = T, label.size = 6.0) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(color = "black", size = 20),
        axis.text =  element_text(color = "black", size = 20)) & 
  NoLegend() 
dev.off()


jpeg("Clusters.By.Orig.ident.jpeg", width = 12.0, height = 6.0, unit = "in", res = 600)
DimPlot(Seurat.SNC.final, label = T, label.size = 3.0,
        split.by = "orig.ident",
        ncol = 4) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(color = "black", size = 10),
        axis.text =  element_text(color = "black", size = 10)) & 
  NoLegend() 
dev.off()



DimPlot(Seurat.SNC.final, label = T, split.by = 'orig.ident')
FeaturePlot(Seurat.SNC.final, features = c("nFeature_RNA", "percent.mt", "decontX_contamination"))

# marker gene expression
NormalizeData(Seurat.SNC.final, assay = "RNA")

DotPlot(Seurat.SNC.final, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7", "SYT1", 
                                      "SYP", "STX1A", "GAD1", "GAD2", "GFAP", "AQP4",
                                      "PLP1", "MOBP", "MBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSF1R", "CD74", "C3", "CLDN5", "FLT1", 
                                      "PDGFRB"), assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

## Rename the SubClusters
Seurat.SNC.final$SubCluster <- NA
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 0] <- "Olig1"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 1] <- "Olig2"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 2] <- "Ex1"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 5] <- "Ex2"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 8] <- "Ex3"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 10] <- "Ex4"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 12] <- "Ex5"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 13] <- "Ex6"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 14] <- "Ex7"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 16] <- "Ex8"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 17] <- "Ex9"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 18] <- "Ex10"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 3] <- "Ast"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 4] <- "In1"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 9] <- "In2"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 11] <- "In3"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 15] <- "In4"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 6] <- "Mic"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 7] <- "OPC"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 19] <- "Endo"
Seurat.SNC.final$SubCluster[Seurat.SNC.final$seurat_clusters == 20] <- "Peri"

## Define Major clusters
Seurat.SNC.final$MajorCluster[grepl("Ex", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Ex"
Seurat.SNC.final$MajorCluster[grepl("In", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "In"
Seurat.SNC.final$MajorCluster[grepl("Olig", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Olig"
Seurat.SNC.final$MajorCluster[grepl("Mic", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Mic"
Seurat.SNC.final$MajorCluster[grepl("Ast", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Ast"
Seurat.SNC.final$MajorCluster[grepl("OPC", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "OPC"
Seurat.SNC.final$MajorCluster[grepl("Endo", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Endo"
Seurat.SNC.final$MajorCluster[grepl("Peri", Seurat.SNC.final$SubCluster, fixed = TRUE)] <- "Peri"

## Dot plot to show marker gene expression
# marker gene expression
Seurat.SNC.final$SubCluster <- factor(Seurat.SNC.final$SubCluster, 
                                     levels = c("Ex1", "Ex2", "Ex3", "Ex4", "Ex5", 
                                                "Ex6", "Ex7", "Ex8", "Ex9", "Ex10",
                                                "In1", "In2", "In3","In4", "Olig1",
                                                "Olig2", "OPC", "Ast", "Mic", 
                                                "Endo", "Peri"))


Idents(Seurat.SNC.final) <- Seurat.SNC.final$SubCluster
NormalizeData(Seurat.SNC.final, assay = "RNA")

pdf("MarkerPlot.pdf", width = 10.0, height = 5.0)
DotPlot(Seurat.SNC.final, features = c("RBFOX3","GABRB2","SATB2", "SLC17A7",
                                      "GAD1", "GAD2", "GFAP", "AQP4", "GJA1",
                                      "ALDH1L1", "MOBP", "MAG", "MOG", "VCAN", 
                                      "PDGFRA", "CSPG4", "BCAN", "CSF1R", "CD74", 
                                      "C3","CLDN5", "FLT1", "PDGFRB"),
        cols = c("lightgrey", "#c04b57"),
        assay = 'RNA') +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") 

dev.off()


# set the subcluser as the Ident
Idents(Seurat.SNC.final) <- Seurat.SNC.final$MajorCluster
jpeg("Clusters.Major.jpeg", width = 7.0, height = 6.0, unit = "in", res = 600)
DimPlot(Seurat.SNC.final, label = T, label.size = 6.0) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(color = "black", size = 20),
        axis.text =  element_text(color = "black", size = 20)) & 
  NoLegend() 
dev.off()

# By Major cell type
jpeg("Clusters.Major.By.Orig.ident.jpeg", width = 12.0, height = 6.0, unit = "in", res = 600)
DimPlot(Seurat.SNC.final, label = T, label.size = 3.0,
        split.by = "orig.ident",
        ncol = 4) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  theme(axis.title = element_text(color = "black", size = 10),
        axis.text =  element_text(color = "black", size = 10)) & 
  NoLegend() 
dev.off()

# FeaturePlot
jpeg("QC.Features.jpeg", width = 12.0, height = 8.0, unit = "in", res = 600)
FeaturePlot(Seurat.SNC.final, 
            features = c("nFeature_RNA", "decontX_contamination", "percent.mt", "percent.rb"),
            cols = c("lightgrey", "#ef4b86")) &
  xlab("UMAP-1") &
  ylab("UMAP-2")

dev.off()


#### Additional dot plot
NormalizeData(Seurat.SNC.final, assay = "RNA")

## set the Major Cluster as the idents
Seurat.SNC.final$MajorCluster <- factor(Seurat.SNC.final$MajorCluster,
                                       levels = c("Ex", "In", "Ast",
                                                  "Olig", "OPC", "Mic",
                                                  "Endo", "Peri"))

Idents(Seurat.SNC.final) <- Seurat.SNC.final$MajorCluster

pdf("AdditonalDotPlot.pdf", width = 4.0, height = 3.0)

DotPlot(Seurat.SNC.final, 
        features = c("BCL11B","STAB2","TBR1", "SOX2", "PAX6", "SLC17A7"), 
        assay = 'RNA', 
        cols = c("lightgrey", "#be3428"),) +
  coord_flip() +
  guides(fill = "none") +
  xlab("Cluster ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_line(linetype = "dashed"),
        axis.text = element_text(color = "black", size = 10))


dev.off()



## Cell Composition by Ind
CellNumbr <- Seurat.SNC.final@meta.data %>%
  dplyr::select(orig.ident, SubCluster) %>%
  group_by(orig.ident, SubCluster) %>%
  mutate(Number.SubCluster = n()) %>%
  unique() %>%
  group_by(SubCluster) %>%
  mutate(Numbr.Toal.Cell = sum(Number.SubCluster)) %>%
  mutate(PercOfInd = Number.SubCluster/Numbr.Toal.Cell *100)
  
## bar chart
p1 <- ggplot(aes(x = SubCluster, y = PercOfInd, fill = orig.ident), data = CellNumbr) +
  geom_bar(stat = "identity") +
  ggsci::scale_fill_aaas()


## Cell Composition by Ind
CellNumbr <- Seurat.SNC.final@meta.data %>%
  dplyr::select(orig.ident, MajorCluster) %>%
  group_by(orig.ident, MajorCluster) %>%
  mutate(Number.MajorCluster = n()) %>%
  unique() %>%
  group_by(MajorCluster) %>%
  mutate(Numbr.Toal.Cell = sum(Number.MajorCluster)) %>%
  mutate(PercOfInd = Number.MajorCluster/Numbr.Toal.Cell *100)

## bar chart
p2 <- ggplot(aes(x = MajorCluster, y = PercOfInd, fill = orig.ident), 
             data = CellNumbr) +
  geom_bar(stat = "identity") +
  ggsci::scale_fill_aaas() +
  ylab("Contribution from each sample (%)") +
  theme_bw() +
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

ggsave("CellComp.pdf", width = 10.0, height = 4.0, plot = p2)


## Cell Composition by Ind
CellNumbr <- Seurat.SNC.final@meta.data %>%
  dplyr::select(orig.ident, MajorCluster) %>%
  group_by(orig.ident, MajorCluster) %>%
  mutate(Number.MajorCluster = n()) %>%
  unique() %>%
  group_by(orig.ident) %>%
  mutate(Numbr.Toal.Cell = sum(Number.MajorCluster)) %>%
  mutate(PercOfCluster = Number.MajorCluster/Numbr.Toal.Cell *100)

## bar chart
p3 <- ggplot(aes(x = orig.ident, y = PercOfCluster, fill = MajorCluster), 
             data = CellNumbr) +
  geom_bar(stat = "identity") +
  ggsci::scale_fill_aaas() +
  ylab("Contribution from each sample (%)") +
  theme_bw() +
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

ggsave("CellComp.byID.pdf", width = 10.0, height = 4.0, plot = p2)

# ## Update the Seurate objects
# saveRDS(Seurat.SNC.final, "Seurat.SNC.final.rds")

## Cell compostion boxplot
## Now calculate the percentage of the cells
Metadata <- Seurat.SNC.final@meta.data
CellNumber <- Metadata %>%
  group_by(orig.ident, MajorCluster, PathDx, TYPE) %>%
  summarise(CellNumberByCellType = n()) %>%
  group_by(orig.ident) %>%
  mutate(TotalByID = sum(CellNumberByCellType)) %>%
  mutate(Percent = CellNumberByCellType/TotalByID * 100)


CellNumber$MajorCluster <- factor(CellNumber$MajorCluster,
                                  levels = c("Ex", "In", "Ast", "Olig",
                                             "OPC", "Mic", "Endo", "Peri"))
CellNumber <- CellNumber %>%
  mutate(group = case_when(PathDx %in% c("DLBD", "TLBD") ~ "LBD",
                           PathDx == "Normal" ~ "Ctrl",
                           PathDx == "Dupl_DLBD" ~ "SNCA Dup",
                           PathDx == "Tri_DLBD" ~ "SNCA Tri"))


CellNumber$group <- factor(CellNumber$group, levels = c("Ctrl", "LBD", "SNCA Dup", "SNCA Tri"))

## Let's visualize the data
p <- ggplot(aes(x = MajorCluster, y = Percent),
            data = CellNumber) +
  geom_boxplot(aes(fill = group), outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.5, aes(color = group), dodge.width = 0.75) +
  scale_color_manual(values = c("black", "black", "black", "black", "black")) +
  scale_x_discrete(labels = c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Endo", "Peri")) +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("#f3756e", "#7bb042", "#1cbdc2", "#a681ba")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 10))

ggsave("CellPercentage.Human.Update.pdf", width = 5.5, height = 3.0, plot = p)

