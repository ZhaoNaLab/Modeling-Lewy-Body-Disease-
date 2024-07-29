library("Seurat")
library("biomaRt")

setwd("./Projects/SNCSynuclein")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")

ProteinEncodingGenes <- biomaRt::getBM(attributes = c("hgnc_symbol",
                                                      "transcript_biotype"), 
                                       filters = "transcript_biotype",
                                       values = "protein_coding", 
                                       mart = mart)

saveRDS(ProteinEncodingGenes, "Step6.DEG/ProteinEncodingGenes.rds")

# load the data
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# Calculate the percentage of cells expressing the genes
DefaultAssay(Seurat.SNC.final) <- "RNA"

MajorCellType <- unique(Seurat.SNC.final$MajorCluster)

FeatureSelect <- list()

for (i in 1:length(MajorCellType)){
  SeuratObj <- subset(Seurat.SNC.final, MajorCluster == MajorCellType[i])
  ExpMat <- GetAssayData(SeuratObj, assay = "RNA", slot = "counts")      
  PercMat <- as.matrix(rowMeans(ExpMat > 0))*100
  colnames(PercMat) <- "Percentage"

  ## include genes expressed by expressed by at least 15% of the cells
  FeaturesIncluded <- rownames(subset(as.data.frame(PercMat), Percentage >= 10))
  
  ## Use only protein-encoding genes
  FeaturesIncluded <- FeaturesIncluded[FeaturesIncluded %in% ProteinEncodingGenes$hgnc_symbol]
  FeatureSelect[[i]] <- FeaturesIncluded
  
}

names(FeatureSelect) <- MajorCellType

# save the data
saveRDS(FeatureSelect, "Step6.DEG/FeatureSelect.rds")











