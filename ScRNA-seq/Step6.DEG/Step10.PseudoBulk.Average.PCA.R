setwd("./Projects/SNCSynuclein")
library("Seurat")
library("tidyverse")
library("DESeq2")
library("vsn")

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
write.csv(MetaData.ind, "Step6.DEG/MetaData.ind.csv", row.names = FALSE)

## Step1: produce the individual-wise Pseudo-bulk
Seurat.SNC.PseudoBulk <- AverageExpression( ## Calculate the sum, can be used for DEG analysis
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


# saveRDS(PseudoBulkMat, "Step6.DEG/PseudoBulk.Average.MatForPCA.rds")
PseudoBulkMat <- readRDS("Step6.DEG/PseudoBulk.Average.MatForPCA.rds")


## robust PCA analysis
# Calculate Row Variance
for (i in 1:length(PseudoBulkMat)){
  dataMatrix <- PseudoBulkMat[[i]]
  
  RowVariance <- matrixStats::rowVars(as.matrix(dataMatrix))
  
  # Select the top 1000 most variable proteins for PCA analysis
  Gene.Select <- order(RowVariance, decreasing=TRUE)[seq_len(min(500, length(RowVariance)))]
  
  # Perform a PCA on the data for the selected proteins
  Norm.PCA <- prcomp(t(dataMatrix[Gene.Select, ]))
  
  # Contribution to the total variance for each component
  percentVar <- Norm.PCA$sdev^2 / sum(Norm.PCA$sdev^2)
  
  # 2D-PCA plot to show the samples
  DatVis <- data.frame(PC1 = Norm.PCA$x[, 1], PC2 = Norm.PCA$x[, 2])
  DatVis <- cbind(DatVis, MetaData.ind)
  
  # Generate PCA plots
  p.PCA <- ggplot(aes(x = PC1, y = PC2, color = PathDx), data = DatVis) +
    geom_point() +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    # scale_color_manual(values = c("#e54c35", "#2b4e98", "#03a087")) +
    theme(plot.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black', size = 12),
          axis.title = element_text(color = 'black', size = 12),
          legend.title = element_blank())
  p.PCA
}


