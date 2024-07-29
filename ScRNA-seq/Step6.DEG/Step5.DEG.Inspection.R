setwd("./Projects/SNCSynuclein")
library("Seurat")
library("tidyverse")

# load the DEG files
DEGFileCombind <- readRDS("Step6.DEG/DEGs/DEGFileCombind.rds")

Clusters <- names(DEGFileCombind)
for (i in 1:length(DEGFileCombind)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame()
  write.csv(dat, paste0("Step6.DEG/DEGs/", Clusters[i], ".Random.DEG.csv"))
}

# load the Seurat object
Seurat.SNC.final <- readRDS("Step5.FinalClustering/Seurat.SNC.final.rds")

# Recode Dx
Seurat.SNC.final$Dx <- NA
Seurat.SNC.final$Dx[grepl("LBD", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "LBD"
Seurat.SNC.final$Dx[grepl("Normal", Seurat.SNC.final$PathDx, fixed = TRUE)] <- "Ctrl"

## Visual inspection for the top DEGs
VlnPlot(Seurat.SNC.final, features = ,
        group.by = "MajorCluster",
        split.by = "Dx",
        split.plot = TRUE)


## Plot 2a: Violin Plot
Seurat.SNC.final$MajorCluster <- factor(Seurat.SNC.final$MajorCluster,
                                       levels = c("Ex", "In", "Ast", "Olig",
                                                  "OPC", "Mic", "Vascular"))

FeatureToPlot <- str_split("SLC11A1/CASP4/SLC15A2/ANKRD17/CCK/FCGR3A/EXOC1/HDAC4/LGALS9/CCR1/CD74/MNDA/C1QA/PPM1B/HRH2/CALCOCO2/FCGR2A/APP/SOD1/CD4/WIPI2/WNT5A/HLA-C/CSF3R/CLOCK/G3BP2/FOXP1/C1QC/TLR5/ATM/MX2/RAB14/CTSC/PRKAR1B/ADAR/XRCC6/HLA-DRA/CLEC7A/PDCD4/EIF2AK4/SLA/CAMK2D/ZC3HAV1/DHX36/TMIGD3/IL1RAP/PRKACA/ZBTB1/FAU/TMF1/PRKACB/RAB27A/RHBDF2/HIF1A/HCK/IL4R/TRIM56/DDX1/HLA-DRB1/STXBP3/NCF1/MX1/RIPK1/PTPRC/HLA-E/PRKAR1A/NDFIP1/RPS19/CD226/SNCA/HSP90AB1/HTRA1/NLRC5/LY96/IFI16/ERAP1/TYROBP/EXT1/SP100/PRDX1/SRPK2/IFI44L/CST3/NAIP/AKNA/LAPTM5/LILRB4/RAB20/USP15/TNFRSF1B/IFNAR1/CREB1/STAT3/HMGB1/SEC14L1/STAT5B/ARID5A/CD81/HLA-A/RBPJ/CYBA/HLA-B/PCBP2/STAT6/IRAK4/IRF2/CDK19/C3/BST2/TNFRSF1A/LRRK2/TRIM14/PPARD/RPL13A/XAF1/LRSAM1/PLSCR1/TLR2/GAPDH/CYBC1/PARP9/TLR1/B2M/ALOX5/TRIM38/ARRB2/TRIM25/PARP14/VPS35/PIK3AP1/TRIM22/STAT2/APPL1/MAPK14/CUL1/PARP4/STAT1/AOAH/IRAK3/CXCL16/ALPK1/LGALS8/HLA-DPB1/NOTCH2/CIITA/CYLD/IFNGR2/BNIP3L",
                           pattern = "\\/", n = Inf)

# normalize the data before plotting
DefaultAssay(Seurat.SNC.final) <- "RNA"
Seurat.SNC.final <- NormalizeData(Seurat.SNC.final, assay = 'RNA')

for (i in 1:length(FeatureToPlot)){
  p <- VlnPlot(Seurat.SNC.final,
               features = FeatureToPlot[i],
               group.by = "MajorCluster",
               split.by = "Dx",
               split.plot = TRUE,
               pt.size = 0,
               ncol = 1) +
    NoLegend()
  
  p <- p + geom_jitter(data = p$data, 
                       aes(color = split),
                       position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
                       size = 0.1, alpha = 0.15) +
    scale_color_manual(values = c("black", 'black'))
  
  ggsave(paste0("Step8.DEG/DEGinspection/", FeatureToPlot[i], ".jpeg"),
         width = 14, height = 4.0, unit = "in", dpi = 600, plot = p)
}





