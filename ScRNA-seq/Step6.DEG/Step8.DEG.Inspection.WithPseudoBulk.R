setwd("./Projects/SNCSynuclein")
library("Seurat")
library("tidyverse")

# load the Seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# normalize the data before plotting
DefaultAssay(Seurat.SNC.final) <- "RNA"
Seurat.SNC.final <- NormalizeData(Seurat.SNC.final, assay = 'RNA')

## get the metadata
MetaData.ind <- Seurat.SNC.final@meta.data %>% 
  dplyr::select(orig.ident, NPID, TYPE, PathDx, Sex, Age, APOE, RIN) %>% 
  unique() %>%
  remove_rownames()


## Save the data
write.csv(MetaData.ind, "Step6.DEG/MetaData.ind.csv", row.names = FALSE)
MetaData.ind <- read.csv("Step6.DEG/MetaData.ind.csv")

## Step1: produce the individual-wise Pseudo-bulk
Seurat.SNC.PseudoBulk <- AverageExpression( ## Calculate the average for data visualization
  Seurat.SNC.final,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = c("orig.ident", "MajorCluster"),
  add.ident = NULL,
  layer = "data",
  verbose = TRUE)


# Seurat.SNC.PseudoBulk2 <- AggregateExpression( ## Calculate the sum, can be used for DEG analysis
#   Seurat.SNC.final,
#   assays = "RNA",
#   return.seurat = FALSE,
#   group.by = c("orig.ident", "MajorCluster"),
#   add.ident = NULL,
#   data = "data",
#   verbose = TRUE)


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


saveRDS(PseudoBulkMat, "Step6.DEG/PseudoBulkMatForVis.rds")
PseudoBulkMat <- readRDS("Step6.DEG/PseudoBulkMatForVis.rds")
# Combine all PseudoBulk into one dataframe and add APOE genotye, 
# Dx, and Sex information into PseudoBulkMat.df

PseudoBulkMat.df <- do.call(rbind, PseudoBulkMat2) %>%
  left_join(MetaData.ind, by = "orig.ident")

# saveRDS(PseudoBulkMat.df, "Step6.DEG/PseudoBulkMat.df.rds")
PseudoBulkMat.df <- readRDS("Step6.DEG/PseudoBulkMat.df.rds")

## Regroup 
PseudoBulkMat.df$Group <- PseudoBulkMat.df$PathDx
PseudoBulkMat.df$Group[PseudoBulkMat.df$PathDx %in% c("DLBD", "TLBD")] <- "LBD"
PseudoBulkMat.df$Group[PseudoBulkMat.df$PathDx %in% c("Dupl_DLBD")] <- "SNCA Dup"
PseudoBulkMat.df$Group[PseudoBulkMat.df$PathDx %in% c("Tri_DLBD")] <- "SNCA Tri"

# reorder Group
PseudoBulkMat.df$Group <- factor(PseudoBulkMat.df$Group, levels = c("Normal", "LBD", "SNCA Dup", "SNCA Tri"))



### Visual inspection for the gene of interest
FeatureToPlot <- unlist(str_split("SLC11A1/CASP4/SLC15A2/ANKRD17/CCK/FCGR3A/EXOC1/HDAC4/LGALS9/CCR1/CD74/MNDA/C1QA/PPM1B/HRH2/CALCOCO2/FCGR2A/APP/SOD1/CD4/WIPI2/WNT5A/HLA-C/CSF3R/CLOCK/G3BP2/FOXP1/C1QC/TLR5/ATM/MX2/RAB14/CTSC/PRKAR1B/ADAR/XRCC6/HLA-DRA/CLEC7A/PDCD4/EIF2AK4/SLA/CAMK2D/ZC3HAV1/DHX36/TMIGD3/IL1RAP/PRKACA/ZBTB1/FAU/TMF1/PRKACB/RAB27A/RHBDF2/HIF1A/HCK/IL4R/TRIM56/DDX1/HLA-DRB1/STXBP3/NCF1/MX1/RIPK1/PTPRC/HLA-E/PRKAR1A/NDFIP1/RPS19/CD226/SNCA/HSP90AB1/HTRA1/NLRC5/LY96/IFI16/ERAP1/TYROBP/EXT1/SP100/PRDX1/SRPK2/IFI44L/CST3/NAIP/AKNA/LAPTM5/LILRB4/RAB20/USP15/TNFRSF1B/IFNAR1/CREB1/STAT3/HMGB1/SEC14L1/STAT5B/ARID5A/CD81/HLA-A/RBPJ/CYBA/HLA-B/PCBP2/STAT6/IRAK4/IRF2/CDK19/C3/BST2/TNFRSF1A/LRRK2/TRIM14/PPARD/RPL13A/XAF1/LRSAM1/PLSCR1/TLR2/GAPDH/CYBC1/PARP9/TLR1/B2M/ALOX5/TRIM38/ARRB2/TRIM25/PARP14/VPS35/PIK3AP1/TRIM22/STAT2/APPL1/MAPK14/CUL1/PARP4/STAT1/AOAH/IRAK3/CXCL16/ALPK1/LGALS8/HLA-DPB1/NOTCH2/CIITA/CYLD/IFNGR2/BNIP3L",
                                  pattern = "\\/", n = Inf))

FeatureToPlot <- c("SYT7", "ETFA", "FADS1", "PRKN", "ALDH5A1", "ABCB10", "GSR", "MTCH2")
for (i in 1:length(FeatureToPlot)){
  p1 <- ggplot(aes(x = Group, y = !!sym(FeatureToPlot[i]), fill = Group), 
              data = subset(PseudoBulkMat.df, MajorCelltype == "Ex")) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.75, size = 0.5) +
    scale_y_continuous(limits = c(0, NA)) +
    # facet_wrap(~ apoe, nrow =3) +
    guides(fill = "none") +
    scale_fill_manual(values = c("#364a9a", "#e36844", "#ec3e39", "#803133")) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 10))
  
  # p2 <- ggplot(aes(x = MajorCelltype, y = !!sym(FeatureToPlot[i]), fill = apoe), 
  #              data = subset(PseudoBulkMat.df, MajorCelltype == "Olig")) +
  #   geom_boxplot() +
  #   ggbeeswarm::geom_quasirandom(dodge.width = 0.75, size = 0.5) +
  #   facet_wrap(~ bg, nrow =3) +
  #   theme_classic()
  
  # p <- cowplot::plot_grid(p1, p2)
  # p
  
  ggsave(paste0("Step6.DEG/DEGinspection.PseudoBulk/", FeatureToPlot[i], ".pdf"),
         width = 3.5, height = 3.0, plot = p1)
}



ggplot(aes(x = Group, y = GOSR1, fill = Group), 
       data = subset(PseudoBulkMat.df, MajorCelltype == "Ex")) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(dodge.width = 0.75, size = 0.5) +
  theme_classic() +
  theme(legend.title = element_blank())






