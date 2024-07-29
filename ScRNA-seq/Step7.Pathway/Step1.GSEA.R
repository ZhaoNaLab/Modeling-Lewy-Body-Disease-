library("tidyverse")
library("clusterProfiler")
library('data.table')
library("BiocParallel")


# prepare from GMT files download from gprofile
ReadGMT_GSEA <- function (gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  fetchName <- function(x){
    unlist(x)[2]
  }
  names(pathways) <- sapply(pathwayLines, fetchName)
  
  temp <- list()
  for (i in 1:length(pathways)){
    temp[[i]] <- as.data.frame(cbind(GO_ID = names(pathways)[i],
                                     gene = pathways [[i]]))
  }
  gmt <- do.call(rbind, temp)
  return(gmt)
}

# load the GMT files
# test with GSEA function 
setwd("/research/labs/moleneurosci/bug/data/ZH_processing_data/SnRNA/GMT/gProfiler")
BPGMT <- ReadGMT_GSEA('hsapiens.GO_BP.name.gmt')
KeggGMT <- read.csv('KEGGMig.csv', row.names = 1)


# load the DEG datasets
setwd("./Projects/SNCSynuclein")

# load the DEG files
DEGFileCombind <- readRDS("Step6.DEG/DEGs/DEGFileCombind.rds")
Cluser <- names(DEGFileCombind)

for (i in 1:length(Cluser)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame()
  dat <- setorder(dat, -z)
  dat$index <- 1:nrow(dat)
  DEGList <- with(dat[, c('primerid', 'z')], setNames(z, primerid))
  DEGList <- DEGList[!is.na(DEGList)]
  dir.create(paste0("./Step7.Pathway/GSEA.", Cluser[i]), recursive = T, showWarnings = F)
  Dir <- paste0("./Step7.Pathway/GSEA.", Cluser[i]) 
    
  # perform GSEA analysis
  setwd(Dir)
  set.seed(2022110)
  GSEA.BP <- GSEA(DEGList, TERM2GENE =BPGMT, verbose = FALSE,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  BPPARAM=MulticoreParam(workers = 4))
  GSEA.BP.Table <- GSEA.BP@result
  
  saveRDS(GSEA.BP, paste0("GSEA.BP.", Cluser[i], ".rds"))
  write.csv(GSEA.BP.Table, paste0("GSEA.BP.Table.", Cluser[i], ".csv"))
  
  set.seed(2022110)
  GSEA.KEGG <- GSEA(DEGList, TERM2GENE = KeggGMT, verbose = FALSE,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    BPPARAM=MulticoreParam(workers = 4))
  GSEA.KEGG.Table <- GSEA.KEGG@result
  
  saveRDS(GSEA.KEGG, paste0("GSEA.KEGG.", Cluser[i], ".rds"))
  write.csv(GSEA.KEGG.Table, paste0("GSEA.KEGG.Table.", Cluser[i], ".csv"))
}


## Let's also save the complete results
for (i in 1:length(Cluser)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame()
  dat <- setorder(dat, -z)
  dat$index <- 1:nrow(dat)
  DEGList <- with(dat[, c('primerid', 'z')], setNames(z, primerid))
  DEGList <- DEGList[!is.na(DEGList)]
  dir.create(paste0("./Step7.Pathway/GSEA.", Cluser[i]), recursive = T, showWarnings = F)
  Dir <- paste0("./Step7.Pathway/GSEA.", Cluser[i]) 
  
  # perform GSEA analysis
  setwd(Dir)
  set.seed(2022110)
  GSEA.BP <- GSEA(DEGList, TERM2GENE =BPGMT, verbose = FALSE,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1.0,
                  BPPARAM=MulticoreParam(workers = 4))
  GSEA.BP.Table <- GSEA.BP@result
  
  # saveRDS(GSEA.BP, paste0("GSEA.BP.", Cluser[i], ".rds"))
  write.csv(GSEA.BP.Table, paste0("GSEA.BP.Table.Complete.", Cluser[i], ".csv"))
  
  set.seed(2022110)
  GSEA.KEGG <- GSEA(DEGList, TERM2GENE = KeggGMT, verbose = FALSE,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1.0,
                    BPPARAM=MulticoreParam(workers = 4))
  GSEA.KEGG.Table <- GSEA.KEGG@result
  
  # saveRDS(GSEA.KEGG, paste0("GSEA.KEGG.", Cluser[i], ".rds"))
  write.csv(GSEA.KEGG.Table, paste0("GSEA.KEGG.Table.Complete.", Cluser[i], ".csv"))
}



