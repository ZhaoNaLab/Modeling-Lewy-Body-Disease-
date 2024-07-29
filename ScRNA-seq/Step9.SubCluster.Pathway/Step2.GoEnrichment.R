library("tidyverse")
library("gprofiler2")
library("data.table")
library("BiocParallel")


# load the DEG datasets
setwd("./Projects/SNCSynuclein")

# load the DEG files
DEGFileCombind <- readRDS("Step8.SubCluster.DEG/DEG.random/DEGFileCombind.rds")
Cluser <- names(DEGFileCombind)

for (i in 1:length(Cluser))
  tryCatch({dat <- DEGFileCombind[[i]] %>% as.data.frame()
  Up.DEG <- subset(dat, DEG == "up")$primerid
  Down.DEG <- subset(dat, DEG == "down")$primerid
  
  # For upregulated DEGs
  gostres.up <- gost(query = Up.DEG, 
                     organism = "hsapiens", ordered_query = TRUE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "g_SCS", 
                     domain_scope = "annotated", custom_bg = NULL, 
                     numeric_ns = "", sources = NULL, as_short_link = FALSE)
  Result.up <- apply(gostres.up$result, 2, as.character) %>% as.data.frame()
  
  # For download DEGs
  gostres.down <- gost(query = Down.DEG, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = TRUE, 
                       user_threshold = 0.05, correction_method = "g_SCS", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE)
  Result.down <- apply(gostres.down$result, 2, as.character) %>% as.data.frame()
  
  write.csv(Result.up, paste0("Step9.SubCluster.Pathway/GOenrcih.result/", Cluser[i], ".up.csv"))
  write.csv(Result.down, paste0("Step9.SubCluster.Pathway/GOenrcih.result/", Cluser[i], ".down.csv"))}, 
           error = function(e){})
  
  

