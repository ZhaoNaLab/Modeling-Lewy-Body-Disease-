setwd("./Projects/SNCSynuclein")
library("Seurat")
library("tidyverse")
library("DESeq2")
library("vsn")
library("factoextra")
library("ggpubr")

# load the Seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# normalize the data before plotting
DefaultAssay(Seurat.SNC.final) <- "RNA"
Seurat.SNC.final <- NormalizeData(Seurat.SNC.final, assay = 'RNA')

## get the metadata
MetaData.ind <- Seurat.SNC.final@meta.data %>% 
  dplyr::select(orig.ident, NPID, TYPE, PathDx, Sex, Age, APOE, RIN) %>% 
  unique() %>%
  remove_rownames() %>%
  column_to_rownames("orig.ident")


## Save the data
# write.csv(MetaData.ind, "Step6.DEG/MetaData.ind.csv", row.names = TRUE)
MetaData.ind <- read.csv("Step6.DEG/MetaData.ind.csv", row.names = 1)

## Step1: produce the individual-wise Pseudo-bulk
Seurat.SNC.PseudoBulk <- AggregateExpression( ## Calculate the sum, can be used for DEG analysis
  Seurat.SNC.final,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = c("orig.ident", "MajorCluster"),
  add.ident = NULL,
  data = "data",
  verbose = TRUE)


## Split the data by cell type
CellType <- c("Ex", "In", "Olig", "OPC", "Ast", "Mic", "Peri", "Endo")
PseudoBulkMat <- list()
PseudoBulkMat2 <- list()

for (i in 1:length(CellType)){
  Dat <- Seurat.SNC.PseudoBulk$RNA %>% as.data.frame() %>%
    dplyr::select(ends_with(CellType[i]))
  
  names(Dat) <- gsub(paste0("_", CellType[i]), "", names(Dat))
  PseudoBulkMat[[i]] <- Dat
  names(PseudoBulkMat)[i] <- CellType[i]
  
  Dat2 <- t(Dat) %>% as.data.frame() %>% rownames_to_column("orig.ident")
  Dat2$MajorCelltype <- CellType[i]
  PseudoBulkMat2[[i]] <- Dat2
}


# saveRDS(PseudoBulkMat, "Step6.DEG/PseudoBulkMatForPCA.rds")
PseudoBulkMat <- readRDS("Step6.DEG/PseudoBulkMatForPCA.rds")


## For the following code, I will follow the DEGseq2 pipline for PCA analysis
plot_list <- list()
for (i in 1:length(PseudoBulkMat)){
  dds <- DESeqDataSetFromMatrix(countData = PseudoBulkMat[[i]],
                                colData = MetaData.ind,
                                design = ~ TYPE)
  # Normalization
  dds <- estimateSizeFactors(dds)

    # Variance stabilization
  vsd <- vst(dds, blind=TRUE)
  meanSdPlot(as.matrix(PseudoBulkMat[[i]])) 
  meanSdPlot(assay(vsd)) 
  
  # PCA plot
  pcaData <- plotPCA(vsd, intgroup=c("PathDx"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pcaData$group <- NA
  pcaData$group[pcaData$PathDx == "Normal"] <- "Normal"
  pcaData$group[pcaData$PathDx %in% c("DLBD", "TLBD")] <- "LBD"
  pcaData$group[pcaData$PathDx == "Dupl_DLBD"] <- "SNCA_Dup"
  pcaData$group[pcaData$PathDx == "Tri_DLBD"] <- "SNCA_Tri"
  
  pcaData$group <- factor(pcaData$group, 
                          levels = c("Normal", "LBD", "SNCA_Dup", "SNCA_Tri"))
  
  plot_list[[names(PseudoBulkMat)[i]]] <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw(base_size = 12.0) +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(color = "black"))
  
  # ggsave(paste0("Step6.DEG/", names(PseudoBulkMat)[i], "_pcaPlot.pdf"), width = 6.0, height = 6.0)
}

pComb <- do.call(ggarrange, c(plot_list, 
                              list(ncol = 2, 
                                   nrow = length(plot_list)/2,
                                   common.legend = TRUE,
                                   align = "hv")))


ggsave("Step6.DEG/CombinedPCA.pdf", width = , height = 1.0, plot = pComb)

## Let's plot PCA plot with circlues
setwd("./Step12.Revision/HumanBrainPCA/")

for (i in 1:length(PseudoBulkMat)){
  dds <- DESeqDataSetFromMatrix(countData = PseudoBulkMat[[i]],
                                colData = MetaData.ind,
                                design = ~ TYPE)
  # Normalization
  dds <- estimateSizeFactors(dds)
  
  # Variance stabilization
  vsd <- vst(dds, blind=TRUE)
  
  dataMatrix <- assay(vsd)
  RowVariance <- matrixStats::rowVars(dataMatrix)
  
  # Select the top N most variable proteins for PCA analysis
  Protein.Select <- order(RowVariance, decreasing=TRUE)[seq_len(min(500, length(RowVariance)))]
  
  
  # Perform a PCA on the data for the selected proteins
  DatMat <- t(dataMatrix[Protein.Select, ])
  
  km.res <- kmeans(DatMat, 3, nstart = 10)
  result <- fviz_cluster(km.res, 
                         DatMat, 
                         geom = c("point"),
                         ellipse.type = "convex",
                         fill = "green",
                         pointsize = 2.0,
                         ellipse.level = 0.95,
                         show.clust.cent = FALSE,
                         stand = TRUE) +
    theme_bw()
  
  
  plot.data <- result$data %>%
    left_join(MetaData.ind %>% rownames_to_column("name"), by = "name")
  
  plot.data$group <- NA
  plot.data$group[plot.data$PathDx == "Normal"] <- "Normal"
  plot.data$group[plot.data$PathDx %in% c("DLBD", "TLBD")] <- "LBD"
  plot.data$group[plot.data$PathDx == "Dupl_DLBD"] <- "SNCA_Dup"
  plot.data$group[plot.data$PathDx == "Tri_DLBD"] <- "SNCA_Tri"
  
  plot.data$group <- factor(plot.data$group, 
                          levels = c("Normal", "LBD", "SNCA_Dup", "SNCA_Tri"))
  
  p <- ggpubr::ggscatter(plot.data, "x", "y", 
                         color = "cluster",
                         shape = "group",
                         ellipse = TRUE, 
                         ellipse.type = "norm") +
    xlab(str_replace(result$labels$x, "Dim1", "PC1")) +
    ylab(str_replace(result$labels$y, "Dim2", "PC2")) +
    guides(color = "none", fill = "none") +
    scale_shape_manual(values = c(0, 1, 2, 3, 5, 7, 8, 9)) +
    ggsci::scale_color_futurama() +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 10))
  
  ggsave(paste0(names(PseudoBulkMat)[i], "_pcaPlot_withCircle.pdf"), width = 4.0, height = 3.0, plot = p)
}


