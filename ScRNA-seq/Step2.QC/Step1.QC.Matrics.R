library("Seurat")
library("tidyverse")
library("ggpubr")
library("purrr")
library("biomaRt")


############# step1: load the Single cell experiment objects #############
setwd("./")
Seurat.SNC.Single.Double <- readRDS("Step1.DataIntegration/Seurat.SNC.Single.Double.rds")

setwd("./Step2.QC")
## QC for doublet identification
FeatureToPlot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.SNC.Single.Double,
          features = FeatureToPlot[i],
          split.by = "DFClass", 
          group.by = "orig.ident",
          pt.size = 0,
          ncol = 1) 
  p <- p + geom_jitter(data = p$data, 
                       position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9),
                       size = 0.1, alpha = 0.1)
   ggsave(paste0("DoubletQC.", FeatureToPlot[i], ".pdf"), width = 8.0, height = 4.0, plot = p)
}


## subset the singlet and remove cell with decontX_contamination >= 0.6 
Seurat.SNC <- subset(Seurat.SNC.Single.Double, DFClass == "Singlet" & decontX_contamination < 0.6)

FeatureToPlot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination")
for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.SNC,
               features = FeatureToPlot[i],
               group.by = "orig.ident",
               pt.size = 0,
               ncol = 1) 
  p <- p + geom_jitter(data = p$data, 
                       position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9),
                       size = 0.1, alpha = 0.1)
  ggsave(paste0("QC.", FeatureToPlot[i], ".pdf"), width = 8.0, height = 4.0, plot = p)
}

# feature scatter plots
pdf("QC.FeatureScatter.pdf", width = 6, height = 6)
FeatureScatter(Seurat.SNC, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()

######################### step 2b: Filtering: remove cells with less than 200 features ######################## 
# Check the most abundant features (plot the median percentage of the features)
AbundantMat <- Seurat.SNC@assays$RNA@counts
AbundantMat <- Matrix::t(Matrix::t(AbundantMat)/Matrix::colSums(AbundantMat)) * 100
HigestExpressed <- order(apply(AbundantMat, 1, median), decreasing = T)[20:1]

# boxplot show the median
par(mar = c(6, 10, 2, 1))
tiff("QC.FeatureAbundance.tif", width = 6.0, height = 4.0, res = 600)
boxplot(t(as.matrix(AbundantMat[HigestExpressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

## MALAT1 is the most abundant gene, check the expression of MALAT1 in each 
VlnPlot(Seurat.SNC,
        features = "MALAT1",
        group.by = "orig.ident",
        pt.size = 0.01,
        ncol = 1) 


# optional: filter genes
# Filter MALAT1: be noted that Seurat autmoatically udpate n_Feature and n_count 
# after the removal of the gene
Seurat.SNC <- Seurat.SNC[!grepl("MALAT1", rownames(Seurat.SNC)), ]

# Filter Mitocondrial
Seurat.SNC <- Seurat.SNC[!grepl("^MT-", rownames(Seurat.SNC)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) Seurat.SNC
# Seurat.SNC <- Seurat.SNC[ ! grepl('^RP[SL]', rownames(Seurat.SNC)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
Seurat.SNC <- Seurat.SNC[!grepl("^HB[^(P)]", rownames(Seurat.SNC)), ]

## update the gene detection rate
Seurat.SNC$cngeneson <- scale(Seurat.SNC$nFeature_RNA)


######################## step 2d: Sample sex information confirmation ########################
# # Get the chorosome information of each gene
# mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# 
# # fetch chromosome info plus some other annotations
# genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
#                                                  "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
#                                   useCache = F))

genes.table <- read.csv("genes.table.csv")
chrY.gene <- genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
SelectY <- intersect(rownames(Seurat.SNC), chrY.gene)

# calculate the percentage of Y chrosome gene
Seurat.SNC$pct_chrY <- colSums(Seurat.SNC@assays$RNA@counts[SelectY, ])/colSums(Seurat.SNC@assays$RNA@counts)
FeatureScatter(Seurat.SNC, feature1 = "XIST", feature2 = "pct_chrY")

pdf("SexConfirmation.pdf", width = 12, height = 6.0)
VlnPlot(Seurat.SNC, features = c("XIST", "pct_chrY"), group.by = "orig.ident",
        ncol = 1)
dev.off()


######################## step 2e: Calculate cell-cycle scores ########################
# Before running CellCycleScoring the data need to be normalized and
# logtransformed.
Seurat.SNC <- NormalizeData(Seurat.SNC)
Seurat.SNC <- CellCycleScoring(object = Seurat.SNC, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

pdf("CellCycleScore.pdf", width = 12, height = 6.0)
VlnPlot(Seurat.SNC, features = c("S.Score", "G2M.Score"), 
        group.by = "orig.ident", ncol = 1, pt.size = 0.1)
dev.off()

## save the Seurat object for clustering analysis
saveRDS(Seurat.SNC, "Seurat.SNC.QC.R1.rds")

## end of the code ##



